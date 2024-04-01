import unittest
import os
from pathlib import Path
import tempfile
import shutil


#from bcm_spectra.utils import get_files
from bcm_spectra.utils import get_files


def test_get_files():
    pass


class TestGetFiles(unittest.TestCase):

    def setUp(self):
        # Setup before each test method
        self.test_dir = Path(__file__).parent / 'test_files'
        if not self.test_dir.exists():
            self.test_dir.mkdir()
        # Create test files
        self.test_mzml = self.test_dir / 'file1.mzML'
        self.test_pepxml = self.test_dir / 'file2.pepXML'
        self.test_mzml.touch(exist_ok=True)
        self.test_pepxml.touch(exist_ok=True)

    def tearDown(self):
        # Cleanup after each test method
        self.test_mzml.unlink(missing_ok=True)
        self.test_pepxml.unlink(missing_ok=True)
        if self.test_dir.exists():
            shutil.rmtree(self.test_dir)

    def test_get_files_single_base_filename(self):
        # Test with a single base filename
        results = get_files(['file1'], str(self.test_dir))
        self.assertIn('file1', results)
        self.assertTrue(len(results['file1']['mzml_files']) > 0)
        self.assertTrue(len(results['file1']['pepxml_files']) == 0)


    def test_get_files_with_matches(self):
        # Test that get_files finds the correct files
        base_filenames = ['file1', 'file2']
        result = get_files(base_filenames, basepath=self.test_dir)

        self.assertIn('file1', result)
        self.assertIn('file2', result)
        self.assertTrue(any(self.test_dir / "file1.mzML" == file for file in result['file1']['mzml_files']))
        self.assertTrue(any(self.test_dir / "file2.pepXML" == file for file in result['file2']['pepxml_files']))

    def test_get_files_no_matches(self):
        # Test that get_files returns an empty list for no matches
        base_filenames = ['nonexistent']
        result = get_files(base_filenames, basepath=self.test_dir)

        self.assertIn('nonexistent', result)
        self.assertEqual(len(result['nonexistent']['mzml_files']), 0)
        self.assertEqual(len(result['nonexistent']['pepxml_files']), 0)

    def test_get_files_empty_base_filenames(self):
        # Test that get_files behaves correctly with an empty list
        base_filenames = []
        result = get_files(base_filenames, basepath=self.test_dir)

        self.assertEqual(result, {})



if __name__ == '__main__':
    unittest.main()
