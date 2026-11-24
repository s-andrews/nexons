"""CLI integration coverage for concurrent BAM files and ordered aggregation."""
from pathlib import Path
import json
import re
import subprocess
import sys
import tempfile
import unittest

import pysam

SCRIPT = Path(__file__).resolve().parents[1] / 'nexons.py'


class ParallelFilesTests(unittest.TestCase):
    def test_serial_and_parallel_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gtf = root / 'genes.gtf'
            gtf.write_text('chr1\tx\texon\t101\t200\t.\t+\t.\tgene_id "g"; transcript_id "t"; transcript_support_level "1";\n')
            inputs = []
            for name, count in [('z', 7), ('a', 3)]:
                path = root / f'{name}.bam'
                inputs.append(str(path))
                with pysam.AlignmentFile(path, 'wb', header={
                    'HD': {'VN': '1.6', 'SO': 'coordinate'},
                    'SQ': [{'SN': 'chr1', 'LN': 1000}]}) as bam:
                    for i in range(count):
                        read = pysam.AlignedSegment(bam.header)
                        read.query_name = f'{name}_{i}'
                        read.reference_id = 0
                        read.reference_start = 100
                        read.cigarstring = '100M'
                        read.query_sequence = 'A' * 100
                        bam.write(read)
            for jobs in (1, 2):
                result = subprocess.run([sys.executable, str(SCRIPT), str(gtf),
                    *inputs, '--startflex', '20', '--endflex', '50', '--parallel', str(jobs), '--outbase', str(root / f'out{jobs}'),
                    *(['--allqc'] if jobs == 2 else [])],
                    capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
            for jobs in (1, 2):
                report = (root / f'out{jobs}_combined_qc.html').read_text()
                self.assertNotIn('%%', report)
                data = json.loads(re.search(r'<script id="qc-data" type="application/json">(.*?)</script>', report, re.S).group(1))
                self.assertEqual(data['names'], ['z', 'a'])
                for sample, count in zip(data['samples'], (7, 3)):
                    self.assertEqual(set(map(int, sample['start_flex'])), set(range(-20, 21)))
                    self.assertEqual(set(map(int, sample['end_flex'])), set(range(-50, 51)))
                    self.assertEqual(sample['start_flex']['0'], count)
                    self.assertEqual(sample['end_flex']['0'], count)

                self.assertEqual([sample['outcomes']['Total_Reads'] for sample in data['samples']], [7, 3])
                self.assertEqual(report.count('<canvas '), 14)
            for category in ('unique', 'partial', 'gene'):
                serial = (root / f'out1_{category}.txt').read_text()
                self.assertEqual(serial, (root / f'out2_{category}.txt').read_text())
                self.assertEqual(serial.splitlines()[0].split('\t')[-2:], inputs)
                self.assertEqual(serial.splitlines()[1].split('\t')[-2:], ['7', '3'])
            for name in ('z', 'a'):
                self.assertEqual((root / f'out1_{name}_nexons_stats.txt').read_text(),
                                 (root / f'out2_{name}_nexons_stats.txt').read_text())
                self.assertTrue((root / f'out2_{name}_qc.html').exists())
                individual = (root / f'out2_{name}_qc.html').read_text()
                self.assertNotIn('%%', individual)
                self.assertIn('id="startflexchart"', individual)
                stats = json.loads((root / f'out2_{name}_nexons_stats.txt').read_text())
                self.assertEqual(set(map(int, stats['start_flex'])), set(range(-20, 21)))
                self.assertIn('labels: ' + str(list(range(-20, 21))), individual)
                self.assertIn('labels: ' + str(list(range(-50, 51))), individual)

                self.assertFalse((root / f'out1_{name}_qc.html').exists())
                with pysam.AlignmentFile(root / f'out1_{name}.bam', 'rb') as a, \
                     pysam.AlignmentFile(root / f'out2_{name}.bam', 'rb') as b:
                    self.assertEqual([r.to_string() for r in a], [r.to_string() for r in b])

    def test_soft_clipping_uses_original_read_orientation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gtf = root / 'genes.gtf'
            gtf.write_text('chr1\tx\texon\t101\t200\t.\t+\t.\tgene_id "g"; transcript_id "t"; transcript_support_level "1";\n')
            path = root / 'clips.bam'
            reads = [
                (0, '500H9S100M20S600H', 129),
                (16, '500H100M10S', 110),
                (0, '100M', 100),
                (0, '999S100M1000S', 2099),
                (0, '1000S100M1001S', 2101),
                (0, '1001S100M999S', 2100),
                (256, '400S100M400S', 900),
                (2048, '400S100M400S', 900),
                (4, None, 100),
            ]
            with pysam.AlignmentFile(path, 'wb', header={
                'SQ': [{'SN': 'chr1', 'LN': 1000}]
            }) as bam:
                for index, (flag, cigar, length) in enumerate(reads):
                    read = pysam.AlignedSegment(bam.header)
                    read.query_name = str(index)
                    read.flag = flag
                    if cigar:
                        read.reference_id = 0
                        read.reference_start = 100
                        read.cigarstring = cigar
                    read.query_sequence = 'A' * length
                    bam.write(read)

            subprocess.run([sys.executable, str(SCRIPT), str(gtf), str(path),
                            '--allqc', '--outbase', str(root / 'out')],
                           capture_output=True, text=True, check=True)
            stats = json.loads((root / 'out_clips_nexons_stats.txt').read_text())
            combined = (root / 'out_combined_qc.html').read_text()
            individual = (root / 'out_clips_qc.html').read_text()
            payload = re.search(r'<script id="qc-data" type="application/json">(.*?)</script>',
                                combined, re.S).group(1)
            sample = json.loads(payload)['samples'][0]
            expected = {
                'start_soft_clip_lengths': {0: 2, 10: 1, 990: 1, 1000: 2},
                'end_soft_clip_lengths': {0: 2, 20: 1, 990: 1, 1000: 2},
            }
            for metric, counts in expected.items():
                bins = [[length, counts.get(length, 0)] for length in range(0, 1001, 10)]
                self.assertEqual(stats[metric], bins)
                self.assertEqual(sample[metric], bins)
            for report in (individual, combined):
                self.assertLess(report.index('<h2>Read Lengths and Coverage</h2>'),
                                report.index('<h2>Soft Clipping</h2>'))
                self.assertIn('id="startsoftclipchart"', report)
                self.assertIn('id="endsoftclipchart"', report)
                self.assertIn('10 bp bins, final bin 1,000+', report)
                self.assertNotIn('%%', report)
            for metric in expected:
                labels, counts = map(list, zip(*stats[metric]))
                self.assertIn(f'{labels}, {counts});', individual)

    def test_invalid_parallel(self):
        for args in (['--parallel', '0', 'x.gtf', 'a.bam'],
                     ['--parallel', '2', 'x.gtf', 'one/a.bam', 'two/a.bam']):
            result = subprocess.run([sys.executable, str(SCRIPT), *args],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn('--parallel', result.stderr)


if __name__ == '__main__':
    unittest.main()
