# utils.py

import os
from collections import defaultdict
import argparse
import logging
#import glob
from pathlib import Path


from . import parser

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_files(base_filenames, basepath=".") -> dict:
    results = dict()
    for base_filename in base_filenames:
        mzml_files = list(Path(basepath).rglob(f'{base_filename}*.mzML'))
        pepxml_files = list(Path(basepath).rglob(f'{base_filename}*.pepXML'))
        results[base_filename] = dict(mzml_files = mzml_files,
                                      pepxml_files = pepxml_files)
    return results

# def get_files(base_filenames: List[str], searchdir=".") -> Dict[str, Dict[str, List[Path]]]:
#     results: Dict[str, Dict[str, List[Path]]] = {}
#     for base_filename in base_filenames:
#         mzml_files = list(Path(searchdir).rglob(f'{base_filename}*.mzML'))
#         pepxml_files = list(Path(searchdir).rglob(f'{base_filename}*.pepXML'))
#         results[base_filename] = {'mzml_files': mzml_files, 'pepxml_files': pepxml_files}
#     return results


def get_filescans(files, df) -> dict:
    filescans = defaultdict()
    # mzml_info = parse_mzml_files(files['53640_1_EXP_MCF7_EGFRa_LF_phos']['mzml_files'])
    # pepxml_info = parse_mzml_files(files['53640_1_EXP_MCF7_EGFRa_LF_phos']['pepxml_files'])

    for file in files.keys():
        scaninfo = df[ df.SpectrumFile == file ]
        mzml_files = files[file]['mzml_files']
        pepxml_files = files[file]['pepxml_files']
        # TODO fix for if there is more than one file # or maybe not necessary if basefilename is unique by fraction (believe it should be)
        if len(mzml_files) == 0:
            logger.error(f"no mzml file found for {file}")
            continue
        if len(pepxml_files) == 0:
            logger.error(f"no pepxml file found for {file}")
            continue

        # I believe this will not happen as basefilename should be unique by fraction & run

        if len(mzml_files) > 1:
            logger.warning(f"more than one mzml file found for {file}")
        if len(pepxml_files) > 1:
            logger.warning(f"more than one pepxml file found for {file}")

        targetscans = scaninfo.FragScanNumber.tolist()

        mzml_file = mzml_files[0]
        pepxml_file = pepxml_files[0]

        logger.info(f"processing {mzml_file} and {pepxml_file}")
        # note here we can calculate the UID for each scan query, check in db for presence of data,
        # if not we use parser.get_scans_from_files to get the data
        # then insert into db
        scans = parser.get_scans_from_files(mzml_file, pepxml_file, targetscans=targetscans)

        filescans[file] = scans

    return filescans

# def get_all_filescans(files, df) -> dict:
#     for file in files.keys