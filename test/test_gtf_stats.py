"""Regression coverage for strand-aware GTF statistics."""
import gzip
import io
import json
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import nexons_gtf_stats as stats


class GtfStatsTests(unittest.TestCase):
    def test_distinct_terminals_by_gene_and_strand(self):
        for strand in ('+', '-'):
            transcripts = {}
            for tid, start, boundary, end in (
                ('a', 10, 100, 500), ('duplicate', 10, 100, 500),
                ('b', 20, 100, 550), ('c', 30, 100, 600),
                ('d', 10, 110, 500), ('e', 20, 110, 550),
                ('incompatible', 40, 120, 700),
            ):
                exons = [(start, boundary), (boundary + 200, end)]
                if strand == '-':
                    exons = [(1000 - b, 1000 - a) for a, b in exons]
                transcripts[tid] = {'strand': strand, 'biotype': 'unknown', 'exons': exons}
            genes = {('chr1', 'g'): transcripts}
            original = stats.collect_metrics(genes, {})
            distinct = stats.collect_metrics(genes, {}, distinct_terminal_positions=True)
            self.assertEqual(distinct['Transcript Start Length'], {6: 2, 16: 1})
            self.assertEqual(distinct['Transcript End Length'], {1: 2, 51: 1})
            self.assertGreater(sum(original['Transcript Start Length'].values()), 3)
            for metric in original:
                if metric not in ('Transcript Start Length', 'Transcript End Length'):
                    self.assertEqual(original[metric], distinct[metric])
            genes[('chr1', 'other_gene')] = transcripts
            twice = stats.collect_metrics(genes, {}, distinct_terminal_positions=True)
            self.assertEqual(twice['Transcript Start Length'], {6: 4, 16: 2})

    def test_terminal_bin_boundaries_and_labels(self):
        for metric, width, cap in (('Transcript Start Length', 5, 1000),
                                   ('Transcript End Length', 50, 10000)):
            counts = stats.Counter()
            distances = [0, 1, width, width + 1, cap - 1, cap, cap + 1]
            groups = {i: [('a', 0), ('b', distance)] for i, distance in enumerate(distances)}
            stats.pair_distances(groups, counts, cap, width)
            self.assertEqual(counts, {0: 1, 1: 2, width + 1: 1, cap - width + 1: 1, cap: 2})
            output = io.StringIO()
            stats.write_metrics({metric: counts}, output)
            self.assertIn(f'{metric}\t0bp\t1\n', output.getvalue())
            self.assertIn(f'{metric}\t1-{width}bp\t2\n', output.getvalue())
            self.assertIn(f'{metric}\t{cap - width + 1}-{cap - 1}bp\t1\n', output.getvalue())
            self.assertIn(f'{metric}\t{cap}+bp\t2\n', output.getvalue())

    def test_exact_zero_counts_preserve_small_nonzero_distances(self):
        genes = {('chr1', 'g'): {
            tid: {'strand': '+', 'biotype': 'unknown',
                  'exons': [(10, 20), (40, end)]}
            for tid, end in [('a', 50), ('b', 50), ('c', 55)]}}
        zeros = stats.Counter()
        metrics = stats.collect_metrics(genes, {}, zeros)
        self.assertEqual(metrics['Transcript End Length'], {0: 1, 1: 2})
        self.assertEqual(zeros['Transcript End Length'], 1)
        self.assertEqual(zeros['Transcript Start Length'], 3)
        self.assertEqual(metrics['Transcript End Length'][1], 2)

    def test_cli_reports_and_embedded_data(self):
        script = Path(stats.__file__).resolve()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input<&.gtf'
            path.write_text(self.fixture().replace('processed_transcript', '</script><unsafe>'))
            for options, prefix in (([], 'nexons_gtf_stats'),
                                    (['--outbase', 'custom', '--chromosome', 'chr1'], 'custom'),
                                    (['--outbase', 'distinct', '--distinct'], 'distinct')):
                subprocess.run([sys.executable, str(script), str(path), *options],
                               cwd=directory, check=True, capture_output=True, text=True)
                text = (Path(directory) / (prefix + '.txt')).read_text()
                expected = io.StringIO()
                stats.write_metrics(stats.collect_metrics(*stats.read_gtf(path)), expected)
                self.assertEqual(text, expected.getvalue())
                report = (Path(directory) / (prefix + '.html')).read_text()
                if prefix == 'distinct':
                    self.assertIn('Distinct terminal position pairs per gene', report)
                self.assertIn('input&lt;&amp;.gtf', report)
                self.assertNotIn('</script><unsafe>', report)
                payload = report.split('<script id="stats-data" type="application/json">')[1].split('</script>')[0]
                data = json.loads(payload)
                self.assertEqual(data['doughnuts'][1]['counts'], [2])
                self.assertEqual(data['doughnuts'][1]['labels'], ['</script><unsafe>'])
                self.assertEqual(len(data['lines']), 6)
                lengths = next(s for s in data['lines'] if s['title'] == 'Mature Transcript Length')
                self.assertEqual(lengths['points'], [{'x': 200, 'y': 1}, {'x': 300, 'y': 1}])
                ends = data['lines'][-1]
                self.assertEqual(ends['width'], 50)
                self.assertEqual(ends['points'][:3], [{'x': 0, 'y': 0}, {'x': 1, 'y': 1}, {'x': 51, 'y': 0}])
                starts = data['lines'][-2]
                self.assertEqual(starts['width'], 5)
                self.assertEqual(next(p for p in starts['points'] if p['x'] == 16)['y'], 1)
                self.assertEqual(ends['points'][-1], {'x': 10000, 'y': 0})
                self.assertIn('chart.resetZoom()', report)

    def fixture(self, strand='+'):
        rows = ['chr1\tx\tgene\t1\t20000\t.\t' + strand + '\t.\tgene_id "g"; gene_type "protein_coding";\n']
        for tid, exons in [('a', [(100, 200), (300, 400), (500, 600)]),
                           ('b', [(120, 200), (300, 380), (500, 635)])]:
            for start, end in exons:
                if strand == '-':
                    start, end = 20000 - end, 20000 - start
                rows.append(f'chr1\tx\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id "g"; transcript_id "{tid}"; transcript_support_level "2"; transcript_type "processed_transcript";\n')
        return ''.join(rows)

    def test_strands_and_compression(self):
        for strand in ('+', '-'):
            for zipped in (False, True):
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / ('input.gtf.gz' if zipped else 'input.gtf')
                    opener = gzip.open if zipped else open
                    with opener(path, 'wt') as handle:
                        handle.write(self.fixture(strand))
                    metrics = stats.collect_metrics(*stats.read_gtf(path))
                self.assertEqual(metrics['Gene Count']['All Genes'], 1)
                self.assertEqual(metrics['Gene Count']['Biotype: protein_coding'], 1)
                self.assertEqual(metrics['Transcript Count']['Biotype: processed_transcript'], 2)
                self.assertEqual(metrics['Transcripts Per Gene'], {2: 1})
                self.assertEqual(metrics['Exons Per Transcript'], {3: 2})
                self.assertEqual(metrics['Mature Transcript Length'], {300: 1, 200: 1})
                self.assertEqual(metrics['Alternate Splice Length'], {20: 1})
                self.assertEqual(metrics['Transcript Start Length'], {16: 1})
                self.assertEqual(metrics['Transcript End Length'], {1: 1})

    def test_tsl_rules(self):
        for raw in ('', 'transcript_support_level "NA";', 'transcript_support_level "3 (assigned)";'):
            self.assertFalse(stats.passes_tsl(stats.attributes(raw), raw, 2))
            self.assertTrue(stats.passes_tsl(stats.attributes(raw), raw, 0))
            for tag in stats.GOOD_TAGS:
                tagged = raw + f' tag "{tag}";'
                self.assertTrue(stats.passes_tsl(stats.attributes(tagged), tagged, 1))
        raw = 'transcript_support_level "2"; transcript_support_level "NA";'
        self.assertTrue(stats.passes_tsl(stats.attributes(raw), raw, 2))

    def test_chromosome_early_stop(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input.gtf'
            path.write_text(self.fixture().replace('chr1', 'chr0') + self.fixture() +
                            'chr2\tx\texon\tbad\tbad\t.\t+\t.\tgene_id "g";\n' + self.fixture())
            genes, _ = stats.read_gtf(path, chromosome='chr1')
            self.assertEqual(list(genes), [('chr1', 'g')])
            self.assertEqual(len(genes[('chr1', 'g')]['a']['exons']), 3)
            self.assertEqual(stats.read_gtf(path, chromosome='missing')[0], {})

    def test_caps_pair_multiplicity_and_zero(self):
        groups = {'boundary': [('a', 0), ('b', 0), ('c', 20000)]}
        for cap, width in ((100, 1), (1000, 5), (10000, 50)):
            counts = stats.Counter()
            stats.pair_distances(groups, counts, cap, width)
            self.assertEqual(counts, {0: 1, cap: 2})
        counts = stats.Counter()
        stats.pair_distances(groups, counts, 100, exclude_equal=True)
        self.assertEqual(counts, {100: 2})
        output = io.StringIO()
        stats.write_metrics(stats.collect_metrics({}, {}), output)
        self.assertIn('Transcript End Length\t10000+bp\t0\n', output.getvalue())
        self.assertIn('Gene Count\tAll Genes\t0\n', output.getvalue())

    def test_terminal_exclusion_and_single_exons(self):
        genes = {('chr1', 'g'): {
            tid: {'strand': '+', 'biotype': 'unknown', 'exons': exons}
            for tid, exons in [('a', [(10, 20)]), ('b', [(10, 30)]),
                               ('c', [(10, 20), (40, 50)])]}}
        metrics = stats.collect_metrics(genes, {})
        for metric in ('Alternate Splice Length', 'Transcript Start Length', 'Transcript End Length'):
            self.assertEqual(metrics[metric], {})


if __name__ == '__main__':
    unittest.main()
