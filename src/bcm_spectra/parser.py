# parser.py

import os
import re
from collections import defaultdict
from typing import List, Dict
from pathlib import Path


import argparse
import logging
#import glob
from pathlib import Path
import pandas as pd

from pyteomics import mzml, pepxml
from pyteomics import mass
import spectrum_utils
from matplotlib import gridspec

import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup
import spectrum_utils.iplot as supi
from tqdm import tqdm
import altair as alt

import matplotlib.pyplot as plt
import matplotlib as mpl


from . import models

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def check_cols(df):
    assert "SurveyScanNumber" in df.columns
    assert "FragScanNumber" in df.columns
    assert "SpectrumFile" in df.columns
    return df

# def get_files(base_filenames):
#     results = dict()
#     for base_filename in base_filenames:
#         mzml_files = list(Path('.').rglob(f'{base_filename}*.mzML'))
#         pepxml_files = list(Path('.').rglob(f'{base_filename}*.pepXML'))
#         results[base_filename] = dict(mzml_files = mzml_files,
#                                       pepxml_files = pepxml_files)
#     return results


annotation_settings = {
    "fragment_tol_mass": 0.05,
    "fragment_tol_mode": "Da",
    "ion_types": "by",
    "max_ion_charge": 2,
    "neutral_losses": {"NH3": -17.026549, "H2O": -18.010565, "H3PO4": -97.976896},
}

def parse_mzml_files(mzml_files):
    # not being used
    # This function will return a dictionary with scan number as key
    # and the parsed mzML information as values
    mzml_info = {}
    for mzml_file in mzml_files:
        with mzml.MzML(mzml_file.__str__()) as reader:
            for spectrum in reader:
                scan_number = spectrum['id']
                mzml_info[scan_number] = spectrum
    return mzml_info

def parse_pepxml_files(pepxml_files):
    # not being used
    # Similar to parse_mzml_files, return a dictionary with scan number as key
    pepxml_info = {}
    for pepxml_file in pepxml_files:
        with pepxml.read(pepxml_file.__str__()) as reader:
            for psm in reader:
                scan_number = psm['start_scan']
                pepxml_info[scan_number] = psm
    return pepxml_info


def extract_ion_annotation(msms: spectrum_utils.spectrum.MsmsSpectrum):
    if msms.annotation is None:
        return None
    annotation = msms.annotation
    ixs = [ i for i, x in enumerate(annotation) if str(x) != '?' ]
    annotation_txt = [ str(x) for x in annotation if str(x) != '?' ]
    annotation_txt = [ x.replace(",", '\n') for x in annotation_txt]
    mzs = [ msms.mz[i] for i in ixs ]
    intensities = [ msms.intensity[i] for i in ixs ]
    mzs = map(lambda x: f"{x:.4f}", mzs)
    intensities = map(lambda x: f"{x:.4f}", intensities)
    df = pd.DataFrame.from_dict(dict(mz=mzs,
                                     annotation=annotation_txt,
                                     intensity=intensities,
                                        ))
    return df





def prepare_ms2_objects(scan: dict, runobj=None):
    """
    this is for ms2 scans mzml and pepxml
    scan is a dictionary with keys 'mzml' and 'pepxml'
    """

    mzmlobj = scan.get('mzml')
    pepxmlobj = scan.get('pepxml')
    # ==

    if mzmlobj is None:
        logger.error("mzmlobj not present")
        raise ValueError("mzmlobj not present")

    if pepxmlobj is None:
        logger.error("pepxmlobj not present")
        raise ValueError("pepxmlobj not present")

    if runobj is None:
        logger.error("runobj not present")
        raise ValueError("runobj not present")


    if not mzmlobj['ms level'] == 2:
        print(f"scan not msms")
        return


    mzarray = mzmlobj['m/z array']
    intensityarray = mzmlobj['intensity array']
    _name = mzmlobj['id'] # concat other info to make more comprehensive
    scanno = mzmlobj['index'] + 1 # concat other info to make more comprehensive
    start_scan = pepxmlobj['start_scan']


    assert scanno == start_scan
    #


    assert mzmlobj['precursorList']['count'] == 1
    precursor = mzmlobj['precursorList']['precursor'][0]

    specref = precursor['spectrumRef']
    _precursor_id_regex = re.match(r'.*scan=(\d+)', specref)
    if _precursor_id_regex is None:
        raise ValueError(f"precursor id not found in {specref}")
    precursor_id = _precursor_id_regex.group(1)


    selected_ion = precursor['selectedIonList']['selectedIon'][0]
    precursor_mz = selected_ion['selected ion m/z']
    precursor_charge = selected_ion['charge state']
    precursor_intensity = selected_ion['peak intensity']


    scanobj = models.Scan(
        scan_number = scanno,
        ms_level = 2,
        mz_array = mzarray,
        intensity_array = intensityarray,
        run = runobj
    )


    fragmentobj = models.Fragment(
        id = scanno,
        scan = scanobj,
    )



    search_hits = pepxmlobj['search_hit']

    search_hit_objects = list()
    for search_hit in search_hits:
        search_hit_obj = process_search_hit(search_hit)
        search_hit_obj.scan = scanobj
        search_hit_obj.fragment = fragmentobj
        search_hit_objects.append(search_hit_obj)


    return {'scan' : scanobj, 'fragment' : fragmentobj, 'search_hits' : search_hit_objects}

def process_search_hit(search_hit):

    def get_massdiff():
        out = dict()
        for ix in modifications_aa:
            aa = modifications_aa[ix]
            if aa not in mass.std_aa_mass:
                ValueError(f"aa {aa} not in mass.std_aa_mass")
            massdiff = mass_shift_dict[ix] - mass.std_aa_mass.get(aa, 0)
            out[ix] = massdiff
        return out

    hit_rank = search_hit['hit_rank']
    hit_massdiff = search_hit['massdiff']

    peptide = search_hit['peptide']
    mass_shift_dict = {mod['position']: mod['mass'] for mod in search_hit['modifications']}
    modifications_aa = {k: peptide[k-1]  for k in mass_shift_dict.keys()}  #1 based index
    massdiff_dict = get_massdiff()
    rank = search_hit['hit_rank']
    score = search_hit['search_score'] # a dictionary {'hyperscore': float, 'next_score': float, 'expect': float}

    proforma_sequence = ''
    for i, aa in enumerate(peptide, start=1):
        proforma_sequence += aa
        if i in massdiff_dict:
            # Add the modification in square brackets
            mass_diff = massdiff_dict[i]
            proforma_sequence += f'[+{mass_diff:.4f}]'



    search_hit_obj = models.SearchResult(
        peptide_sequence = peptide,
        proforma_sequence = proforma_sequence,
        rank = rank,
        score = score,
        mass_shifts = mass_shift_dict,
        mass_diffs = massdiff_dict,
    )


    return search_hit_obj



def other():

    search_hits = pepxmlobj['search_hit']

    for search_hit in search_hits:
        # this has many many little ways to get tricked up. so many things to map

        #dict_keys(['peptide', 'num_missed_cleavages', 'num_tot_proteins', 'tot_num_ions', 'hit_rank', 'num_matched_ions', 'search_score', 'modified_peptide', 'massdiff', 'calc_neutral_pep_mass', 'is_rejected', 'proteins', 'modifications'])
        hit_rank = search_hit['hit_rank']
        hit_massdiff = search_hit['massdiff']

        peptide = search_hit['peptide']
        modifications = {mod['position']: mod['mass'] for mod in search_hit['modifications']}
        modifications_aa = {k: peptide[k-1]  for k in modifications.keys()}  #1 based index
        massdiff_dict = get_massdiff()

        proforma_sequence = ''
        for i, aa in enumerate(peptide, start=1):
            proforma_sequence += aa
            if i in massdiff_dict:
                # Add the modification in square brackets
                mass_diff = massdiff_dict[i]
                proforma_sequence += f'[+{mass_diff:.4f}]'

        msms = sus.MsmsSpectrum(identifier = _name,
                        precursor_mz = float(precursor_mz),
                        precursor_charge = precursor_charge,
                        mz = mzarray,
                        intensity = intensityarray,
                        retention_time = 1.,
                        #peptide = peptide,
                        #modifications = modifications
                        )

        msms.annotate_proforma(proforma_sequence,
                               **annotation_settings
                               )

        # fig, ax = plt.subplots()
        # sup.spectrum(msms, ax=ax, grid=False)
        msms_annot_table = extract_ion_annotation(msms)


        outpath = Path("./spectra")
        outpath.mkdir(exist_ok=True)

def check_cols(df):
    assert "SurveyScanNumber" in df.columns
    assert "FragScanNumber" in df.columns
    assert "SpectrumFile" in df.columns
    return df



# def main(targetfile: str, searchdir: str):
#     """
#     :param targetfile: str: path to the target csv file that contains the survey scan number, fragment scan number, and spectrum file
#     :param searchdir: str: path to the directory where the mzML and pepXML files are located
#     """
#     df = pd.read_csv(targetfile)
#     df = check_cols(df)
#     base_filenames = df.SpectrumFile.unique()
#     files = get_files(base_filenames, searchdir=searchdir)
#
#     for file in files:
#         fileinfo = files[file]
#         mzml_files = fileinfo.get('mzml_files')
#         pepxml_files = fileinfo.get('pepxml_files')
#         mzml_file = mzml_files[0]
#         pepxml_file = pepxml_files[0]
#
#
#
#         scaninfo = df[ df.SpectrumFile == file ]
#         targetscans = scaninfo.FragScanNumber.tolist()
#
#         scans = get_scans_from_files(mzml_file, pepxml_file, targetscans=targetscans)
#
#


def get_scans_from_files(mzml_file, pepxml_file, targetscans=None, _all=False):
    """
    :param fileinfo: dict: dictionary containing the mzml_files and pepxml_files
    :param targetscans: list: list of scan numbers to extract. MS2 only right now
    """
    logger.info(f"processing {mzml_file} and {pepxml_file}")

    scans = defaultdict(dict)

    mzmlinfo = dict()
    pepxmlinfo = dict()

    with mzml.MzML(mzml_file.__str__()) as reader:

        logger.info(f"processing {mzml_file}")

        for spectrum in tqdm(reader):
            scan_number = spectrum['index'] + 1
            if scan_number in targetscans:
                print(f"got scan {scan_number}")
                #scans[scan_number]['mzml'] = spectrum
                mzmlinfo[scan_number] = spectrum

            if all([x in mzmlinfo for x in targetscans]) and not _all:
                break

    with pepxml.PepXML(pepxml_file.__str__()) as reader:

        logger.info(f"processing {pepxml_file}")

        for spectrum in tqdm(reader):
            scan_number = spectrum['start_scan']
            if scan_number in targetscans:
                print(f"got scan {scan_number}")
                #scans[scan_number]['pepxml'] = spectrum
                pepxmlinfo[scan_number] = spectrum

            if all([x in pepxmlinfo for x in targetscans]) and not _all:
                break

    for key in mzmlinfo:
        scans[key]['mzml'] = mzmlinfo[key]
    for key in pepxmlinfo:
        scans[key]['pepxml'] = pepxmlinfo[key]

    return scans