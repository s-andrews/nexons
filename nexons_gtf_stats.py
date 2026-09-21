#!/usr/bin/env python3
"""Summarise retained GTF exons using nexons-compatible TSL filtering.

Counts describe genes/transcripts with at least one retained exon. Coordinates
are inclusive. Pairs are unordered and restricted to a gene and strand; repeated
boundaries in different transcripts contribute separately. Identical terminal
boundaries contribute a zero distance; identical internal ends do not.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import html
from itertools import combinations
import json
from pathlib import Path
import re
import shlex
import sys


GOOD_TAGS = ('MANE_Select', 'Ensembl_Canonical', 'gencode_primary', 'gencode_basic')


def attributes(text):
    result = defaultdict(list)
    for item in text.split(';'):
        parts = item.strip().split(None, 1)
        if len(parts) == 2:
            result[parts[0]].append(parts[1].strip().strip('"'))
    return result


def first(attrs, *keys):
    return next((attrs[k][0] for k in keys if attrs.get(k)), None)


def passes_tsl(attrs, raw, maximum):
    """Match nexons.read_gtf's exon-level TSL and preferred-tag rules."""
    tsl = None
    for value in attrs.get('transcript_support_level', []):
        words = value.split()
        if words and words[0].isdigit():
            tsl = int(words[0])
    if any(tag in raw for tag in GOOD_TAGS):
        tsl = 1
    return maximum in (None, 0) or (tsl is not None and tsl <= maximum)


def read_gtf(path, max_tsl=2, chromosome=None):
    genes = {}
    gene_types, transcript_types = {}, {}
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    seen_chromosome = False
    with opener(path, 'rt', encoding='utf8') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if chromosome is not None:
                if fields[0] != chromosome:
                    if seen_chromosome:
                        break
                    continue
                seen_chromosome = True
            if len(fields) != 9:
                print(f'Warning: skipping malformed GTF line {number}', file=sys.stderr)
                continue
            chrom, _, feature, start, end, _, strand, _, raw = fields
            if feature not in ('gene', 'transcript', 'exon'):
                continue
            attrs = attributes(raw)
            gid = first(attrs, 'gene_id', 'gene_name')
            tid = first(attrs, 'transcript_id', 'transcript_name')
            if gid is None:
                continue
            gkey = (chrom, gid)
            tkey = (chrom, gid, tid)
            gb = first(attrs, 'gene_biotype', 'gene_type')
            tb = first(attrs, 'transcript_biotype', 'transcript_type')
            if gb:
                gene_types[gkey] = gb
            if tb and tid is not None:
                transcript_types[tkey] = tb
            if feature != 'exon' or not passes_tsl(attrs, raw, max_tsl):
                continue
            if tid is None or strand not in ('+', '-'):
                print(f'Warning: skipping exon without transcript ID or strand at line {number}', file=sys.stderr)
                continue
            start, end = int(start), int(end)
            if start < 1 or end < start:
                raise ValueError(f'Invalid exon coordinates at line {number}')
            transcripts = genes.setdefault(gkey, {})
            transcript = transcripts.setdefault(tid, {'strand': strand, 'exons': []})
            if transcript['strand'] != strand:
                raise ValueError(f'Inconsistent transcript strand at line {number}')
            transcript['exons'].append((start, end))
    for gkey, transcripts in genes.items():
        for tid, transcript in transcripts.items():
            transcript['exons'].sort(reverse=transcript['strand'] == '-')
            transcript['biotype'] = transcript_types.get((*gkey, tid), gene_types.get(gkey, 'unknown'))
    return genes, gene_types


def pair_distances(groups, counter, cap, width=1, exclude_equal=False,
                   distinct_positions=False):
    zero_count = 0
    seen_pairs = set()
    for group, entries in groups.items():
        if distinct_positions:
            # Collapse terminal positions within each compatible boundary group
            # before pairing, avoiding quadratic work in transcript duplicates.
            entries = list(enumerate(sorted({position for _, position in entries})))
        for (tid_a, end_a), (tid_b, end_b) in combinations(entries, 2):
            if tid_a == tid_b:
                continue
            if distinct_positions:
                # A position pair can qualify through multiple shared boundaries;
                # count it only once per gene and strand, not once per boundary.
                pair = (group[0], end_a, end_b)
                if pair in seen_pairs:
                    continue
                seen_pairs.add(pair)
            distance = abs(end_a - end_b)
            if distance == 0:
                zero_count += 1
            if exclude_equal and distance == 0:
                continue
            if distance >= cap:
                key = cap
            elif distance == 0:
                key = 0
            else:
                key = (distance - 1) // width * width + 1
            counter[key] += 1
    return zero_count


def collect_metrics(genes, gene_types, zero_counts=None, distinct_terminal_positions=False):
    if zero_counts is None:
        zero_counts = Counter()
    metrics = {name: Counter() for name in (
        'Gene Count', 'Transcript Count', 'Transcripts Per Gene',
        'Exons Per Transcript', 'Mature Transcript Length',
        'Alternate Splice Length', 'Transcript Start Length', 'Transcript End Length')}
    metrics['Gene Count']['All Genes'] = len(genes)
    metrics['Transcript Count']['All Transcripts'] = sum(map(len, genes.values()))
    for gkey, transcripts in genes.items():
        metrics['Gene Count']['Biotype: ' + gene_types.get(gkey, 'unknown')] += 1
        metrics['Transcripts Per Gene'][len(transcripts)] += 1
        internal, starts, ends = defaultdict(list), defaultdict(list), defaultdict(list)
        for tid, transcript in transcripts.items():
            exons, strand = transcript['exons'], transcript['strand']
            metrics['Transcript Count']['Biotype: ' + transcript['biotype']] += 1
            metrics['Exons Per Transcript'][len(exons)] += 1
            length = sum(end - start + 1 for start, end in exons)
            metrics['Mature Transcript Length'][length // 100 * 100] += 1
            # Express each exon as (5-prime boundary, 3-prime boundary).
            directed = [(a, b) if strand == '+' else (b, a) for a, b in exons]
            for start, end in directed[1:-1]:
                internal[(strand, start)].append((tid, end))
            if len(directed) > 1:
                start, end = directed[0]
                starts[(strand, end)].append((tid, start))
                start, end = directed[-1]
                ends[(strand, start)].append((tid, end))
        pair_distances(internal, metrics['Alternate Splice Length'], 100, exclude_equal=True)
        for metric, groups, cap, width in (
            ('Transcript Start Length', starts, 1000, 5),
            ('Transcript End Length', ends, 10000, 50),
        ):
            zero_counts[metric] += pair_distances(
                groups, metrics[metric], cap, width,
                distinct_positions=distinct_terminal_positions)
    return metrics


def distance_bins(cap, width, include_zero=True):
    """Positive bins start at 1; the cap is a separate overflow class."""
    return ([0] if include_zero else []) + list(range(1, cap, width)) + [cap]


def write_metrics(metrics, handle):
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(('Metric', 'Value', 'Count'))
    caps = {'Alternate Splice Length': (100, 1),
            'Transcript Start Length': (1000, 5), 'Transcript End Length': (10000, 50)}
    for metric, counts in metrics.items():
        if metric in caps:
            cap, width = caps[metric]
            keys = distance_bins(cap, width, metric != 'Alternate Splice Length')
        else:
            keys = sorted(counts)
        for key in keys:
            label = str(key)
            if metric in caps:
                cap, width = caps[metric]
                label = f'{key}+bp' if key == cap else (
                    f'{key}bp' if width == 1 or key == 0 else
                    f'{key}-{min(key + width - 1, cap - 1)}bp')
            elif metric == 'Mature Transcript Length':
                label = f'{key}-{key + 99}bp'
            writer.writerow((metric, label, counts[key]))


def write_html(metrics, args, command, zero_counts=None):
    """Embed report data safely in the Bootstrap/Chart.js report template."""
    doughnuts = []
    for metric in ('Gene Count', 'Transcript Count'):
        counts = metrics[metric]
        labels = sorted((key for key in counts if key.startswith('Biotype: ')),
                        key=lambda key: (-counts[key], key))
        doughnuts.append({'title': metric.replace('Count', 'Biotypes'),
                         'labels': [key.removeprefix('Biotype: ') for key in labels],
                         'counts': [counts[key] for key in labels]})
    lines = []
    terminal_counting = ('Each distinct pair of positions is counted once per gene.'
                         if args.distinct else 'Each qualifying transcript pair is counted once.')
    descriptions = {
        'Transcripts Per Gene': 'How many transcripts does each gene have.',
        'Exons Per Transcript': 'How many exons are in each transcript.',
        'Mature Transcript Length': 'The overall length of mature transcripts after splicing (100bp resolution).',
        'Alternate Splice Length': 'Absolute distance between splice donor sites for exons of different isoforms which start at the same position.',
        'Transcript Start Length': 'Absolute distance between transcript start sites for first exons with identical splice donor site positions. ' + terminal_counting,
        'Transcript End Length': 'Absolute distance between terminator positions for terminal exons with identical splice acceptor site positions. ' + terminal_counting,
    }
    for metric, xlabel, ylabel, width, cap in (
        ('Transcripts Per Gene', 'Transcripts per gene', 'Genes', 1, None),
        ('Exons Per Transcript', 'Exons per transcript', 'Transcripts', 1, None),
        ('Mature Transcript Length', 'Mature transcript length (bp)', 'Transcripts', 100, None),
        ('Alternate Splice Length', 'Distance to alternate donor (bp)', 'Exon pairs', 1, 100),
        ('Transcript Start Length', 'Distance to alternate TSS (bp)', 'Transcript pairs', 5, 1000),
        ('Transcript End Length', 'Distance to alternate Terminator (bp)', 'Transcript pairs', 50, 10000),
    ):
        counts = metrics[metric]
        keys = sorted(counts)
        if cap is not None:
            keys = distance_bins(cap, width, metric != 'Alternate Splice Length')
        else:
            # Keep large sparse distributions compact without drawing across empty bins.
            padded = set(keys)
            for left, right in zip(keys, keys[1:]):
                if right - left > width:
                    padded.update((left + width, right - width))
            keys = sorted(padded)
        lines.append({'title': metric, 'description': descriptions[metric],
                      'xlabel': xlabel, 'ylabel': ylabel,
                      'width': width, 'cap': cap,
                      'zeroCount': (zero_counts or {}).get(metric, 0),
                      'points': [{'x': key, 'y': counts[key]} for key in keys]})
        if args.distinct and metric in ('Transcript Start Length', 'Transcript End Length'):
            lines[-1]['ylabel'] = 'Distinct terminal position pairs'
    summary = [('File analysed', args.gtf), ('Command line', command),
               ('Maximum TSL', args.maxtsl if args.maxtsl else 'Disabled (0)'),
               ('Chromosome', args.chromosome or 'All chromosomes'),
               ('Output prefix', args.outbase),
               ('Terminal counting', 'Distinct terminal position pairs per gene'
                if args.distinct else 'All qualifying transcript pairs'),
               ('Total genes', metrics['Gene Count']['All Genes']),
               ('Total transcripts', metrics['Transcript Count']['All Transcripts'])]
    rows = '\n'.join(f'<tr><th scope="row">{html.escape(label)}</th>'
                     f'<td>{html.escape(str(value))}</td></tr>' for label, value in summary)
    data = json.dumps({'doughnuts': doughnuts, 'lines': lines}).replace('<', '\\u003c')
    template = Path(__file__).parent / 'templates/nexons_gtf_stats_template.html'
    report = template.read_text(encoding='utf8')
    replacements = {'%%SUMMARY%%': rows, '%%DATA%%': data}
    report = re.sub(r'%%(?:SUMMARY|DATA)%%', lambda match: replacements[match[0]], report)
    Path(str(args.outbase) + '.html').write_text(report, encoding='utf8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gtf', help='Input GTF or gzipped GTF')
    parser.add_argument('--maxtsl', type=int, default=2,
                        help='Maximum transcript support level (default 2; 0 disables filtering)')
    parser.add_argument('--chromosome', '--chrom', help='Process one chromosome (exact name; input must be grouped by chromosome)')
    parser.add_argument('--distinct', action='store_true',
                        help='Count each distinct terminal position pair once per gene for start/end distances; shared-boundary rules still apply; no zero-distance pairs')
    parser.add_argument('-o', '--outbase', default='nexons_gtf_stats',
                        help='Output prefix; appends .txt and .html (default nexons_gtf_stats)')
    args = parser.parse_args()
    genes, gene_types = read_gtf(args.gtf, args.maxtsl, args.chromosome)
    zero_counts = Counter()
    metrics = collect_metrics(genes, gene_types, zero_counts, args.distinct)
    with open(args.outbase + '.txt', 'w', encoding='utf8', newline='') as handle:
        write_metrics(metrics, handle)
    write_html(metrics, args, shlex.join(sys.argv), zero_counts)


if __name__ == '__main__':
    main()
