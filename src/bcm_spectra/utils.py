# utils.py

import os
from collections import defaultdict
import argparse
import logging
#import glob
from pathlib import Path

def get_files(base_filenames, basepath="."):
    results = dict()
    for base_filename in base_filenames:
        mzml_files = list(Path(basepath).rglob(f'{base_filename}*.mzML'))
        pepxml_files = list(Path(basepath).rglob(f'{base_filename}*.pepXML'))
        results[base_filename] = dict(mzml_files = mzml_files,
                                      pepxml_files = pepxml_files)
    return results

