# utils.py

import os
from collections import defaultdict
import argparse
import logging

# import glob
from pathlib import Path

import collections.abc # for isinstance check of iterable


import sqlalchemy
from sqlalchemy.orm import Session

from . import parser
from . import db
from . import models
from . import crud

from spectrum_utils import fragment_annotation as fa, proforma #, utils

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_files(base_filenames, basepath=".") -> dict:
    results = dict()
    for base_filename in base_filenames:
        mzml_files = list(Path(basepath).rglob(f"{base_filename}*.mzML"))
        pepxml_files = list(Path(basepath).rglob(f"{base_filename}*.pepXML"))
        results[base_filename] = dict(mzml_files=mzml_files, pepxml_files=pepxml_files)
    return results


# def get_files(base_filenames: List[str], searchdir=".") -> Dict[str, Dict[str, List[Path]]]:
#     results: Dict[str, Dict[str, List[Path]]] = {}
#     for base_filename in base_filenames:
#         mzml_files = list(Path(searchdir).rglob(f'{base_filename}*.mzML'))
#         pepxml_files = list(Path(searchdir).rglob(f'{base_filename}*.pepXML'))
#         results[base_filename] = {'mzml_files': mzml_files, 'pepxml_files': pepxml_files}
#     return results


def get_filescans(mzml_file, pepxml_file, session=None, runobj=None, targetscans=None):
    """
    here is logic for fetching scan data from db if present, else fetch from files and save to db
    this is for ms2 scans
    note this gets scan info AND SearchEngine info
    multimodal creation
    """

    logger.info(f"processing {mzml_file} and {pepxml_file}")

    # note here we can calculate the UID for each scan query, check in db for presence of data,
    # if not we use parser.get_scans_from_files to get the data
    # then insert into db

    if session is None:
        raise ValueError("session not provided")
        # session = db.get_session()

    if runobj is None:
        logger.warning("runobj not present")

    scans = dict()

    for targetscan in targetscans:
        # scan = crud.get_scan_from_run_by_scan_number(runobj, targetscan)
        scan = crud.get_scan_from_run_by_scan_number(
            session=session, run=runobj, scan_number=targetscan
        )
        logger.info(f"looking for {targetscan} in db")
        if scan is not None:
            logger.info(f"fetching scan {targetscan} from db")
            scans[targetscan] = scan

    missing_scans = set(targetscans) - set(scans.keys())

    missing_scans_res = dict()
    if len(missing_scans) > 0:
        missing_scans_res = parser.get_scans_from_files(
            mzml_file, pepxml_file, targetscans=missing_scans
        )
        # this takes a while
        # no scans are saved to db until after all scans are processed
        # this could be changed / modified

    new_obj_collection = list()

    for scan in missing_scans_res.values():
        # actually we make new Scan objects here
        new_objs = parser.prepare_ms2_objects(scan, runobj=runobj)
        new_obj_collection.append(new_objs)

    # scans.update(missing_scans_res)

    # add new scans to db
    import ipdb; ipdb.set_trace()
    crud.commit_all_objects(new_obj_collection, session)

    # for new_objs in new_obj_collection:
    #     for objname, obj in new_objs.items():
    #         if isinstance(obj, models.Base):
    #             session.add(obj)
    #         if objname == "scan":  # we want to keep these and return them
    #             scans[obj.scan_number] = obj # these guys have links to everything we might need
    #         elif isinstance(obj, collections.abc.Iterable) and not isinstance(obj, (str, bytes)):  # Handling any iterable object but excluding strings/bytes
    #             for item in obj:
    #                 if isinstance(item, models.Base):
    #                     session.add(item)
    #                 else:
    #                     raise TypeError(f"Non-model object in iterable: {type(item).__name__}")
    #         else:
    #             raise TypeError(f"Unsupported type in new_obj_collection: {type(obj).__name__}")
    #         else:
    #             for o in obj:  # this is for search hits.
    #                 session.add(o)

    session.commit()

    return scans


def get_all_filescans(files, df, session: Session) -> dict:
    filescans = defaultdict()
    # mzml_info = parse_mzml_files(files['53640_1_EXP_MCF7_EGFRa_LF_phos']['mzml_files'])
    # pepxml_info = parse_mzml_files(files['53640_1_EXP_MCF7_EGFRa_LF_phos']['pepxml_files'])
    for file in files.keys():

        runobj = crud.get_or_create_run(filename=file, session=session)

        mzml_files = files[file]["mzml_files"]
        pepxml_files = files[file]["pepxml_files"]
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
        mzml_file = mzml_files[0]
        pepxml_file = pepxml_files[0]

        scaninfo = df[df.SpectrumFile == file]
        targetscans = scaninfo.FragScanNumber.tolist()

        scans_for_file = get_filescans(
            mzml_file,
            pepxml_file,
            session=session,
            runobj=runobj,
            targetscans=targetscans,
        )

        filescans[file] = scans_for_file

    return filescans


# def get_all_filescans(files, df) -> dict:
#     for file in files.keys
