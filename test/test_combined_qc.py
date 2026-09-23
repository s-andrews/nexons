import argparse
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import nexons


class CombinedQcTests(unittest.TestCase):
    def test_empty_and_single_sample(self):
        sample = {'file': 'sample<&".bam', 'outcomes': {'Total_Reads': 0},
                  'read_lengths': [], 'start_flex': {-1: 2, 0: 3, 1: 1}, 'end_flex': {-2: 0, 0: 0},
                  'inner_flex': {0: 0}, 'coverage': [0] * 101}
        with tempfile.TemporaryDirectory() as directory:
            prefix = str(Path(directory) / 'run')
            for samples in ([], [sample]):
                nexons.write_combined_qc_report(samples,
                    argparse.Namespace(gtf='</script><script>alert(1)</script>'), prefix)
                report = Path(prefix + '_combined_qc.html').read_text()
                self.assertNotIn('%%', report)
                self.assertNotIn('<script>alert(1)</script>', report)
                payload = re.search(r'<script id="qc-data" type="application/json">(.*?)</script>', report, re.S).group(1)
                data = json.loads(payload)
                self.assertEqual(len(data['samples']), len(samples))
                if samples:
                    self.assertEqual(data['names'], ['sample<&"'])
                    self.assertEqual(data['samples'][0]['start_flex'], {'-1': 2, '0': 3, '1': 1})
                    self.assertIn('id="startflexchart"', report)
                    self.assertIn('sample&lt;&amp;&quot;', report)
                    self.assertNotIn('<', payload)


if __name__ == '__main__':
    unittest.main()
