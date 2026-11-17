import tempfile
from pathlib import Path
import unittest

import h5py
import numpy as np

from pipeline.rebuild import column
from pipeline.inputs import DEFINITIONS
from pipeline.io import read_json


class RebuildTests(unittest.TestCase):
    def test_categorical_missing_and_order(self):
        with tempfile.TemporaryDirectory() as temp:
            with h5py.File(Path(temp) / 'test.h5', 'w') as f:
                group = f.create_group('label')
                group.create_dataset('categories', data=np.array([b'A', b'B']))
                group.create_dataset('codes', data=np.array([1, -1, 0]))
                self.assertEqual(column(f, b'label').tolist(), ['B', '', 'A'])

    def test_invalid_categorical_code_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            with h5py.File(Path(temp) / 'test.h5', 'w') as f:
                group = f.create_group('label')
                group.create_dataset('categories', data=np.array([b'A']))
                group.create_dataset('codes', data=np.array([2]))
                with self.assertRaisesRegex(ValueError, 'categorical'):
                    column(f, 'label')

    def test_source_objects_pinned(self):
        manifest = read_json(DEFINITIONS / 'raw_reference_sources.json')
        self.assertEqual(len(manifest['files']), 11)
        self.assertEqual(len({f['name'] for f in manifest['files']}), 11)
        for source in manifest['files']:
            self.assertEqual(Path(source['name']).name, source['name'])
            self.assertRegex(source['sha256'], r'^[0-9a-f]{64}$')
            self.assertTrue(source['url'].startswith('https://'))


if __name__ == '__main__':
    unittest.main()
