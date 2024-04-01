# main.py

import os
from collections import defaultdict
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

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def check_cols(df):
    assert "SurveyScanNumber" in df.columns
    assert "FragScanNumber" in df.columns
    assert "SpectrumFile" in df.columns
    return df

def get_files(base_filenames):
    results = dict()
    for base_filename in base_filenames:
        mzml_files = list(Path('.').rglob(f'{base_filename}*.mzML'))
        pepxml_files = list(Path('.').rglob(f'{base_filename}*.pepXML'))
        results[base_filename] = dict(mzml_files = mzml_files,
                                      pepxml_files = pepxml_files)
    return results

def get_scan_info(df, files):
    df.SpectrumFile.unique()

# def parse_mzml_file(mzml_files):
#     # This function will return a dictionary with scan number as key
#     # and the parsed mzML information as values
#     mzml_info = {}
#     for mzml_file in mzml_files:
#         with mzml.MzML(mzml_file) as reader:
#             for spectrum in reader:
#                 scan_number = spectrum['id']
#                 mzml_info[scan_number] = spectrum
#     return mzml_info
#
# def parse_pepxml_files(pepxml_files):
#     # Similar to parse_mzml_files, return a dictionary with scan number as key
#     pepxml_info = {}
#     for pepxml_file in pepxml_files:
#         with pepxml.read(pepxml_file) as reader:
#             for psm in reader:
#                 scan_number = psm['start_scan']
#                 pepxml_info[scan_number] = psm
#     return pepxml_info

annotation_settings = {
    "fragment_tol_mass": 0.05,
    "fragment_tol_mode": "Da",
    "ion_types": "by",
    "max_ion_charge": 2,
    "neutral_losses": {"NH3": -17.026549, "H2O": -18.010565, "H3PO4": -97.976896},
}



#def get_scan_info(df, scan_numbers)

def handle_scan(scan: dict):

    mzmlobj = scan.get('mzml')
    pepxmlobj = scan.get('pepxml')

    if not mzmlobj['ms level'] == 2:
        print(f"scan['id'] not msms")
        return


    mzarray = mzmlobj['m/z array']
    intensityarray = mzmlobj['intensity array']
    _name = mzmlobj['id'] # concat other info to make more comprehensive
    scanno = mzmlobj['index'] + 1 # concat other info to make more comprehensive
    start_scan = pepxmlobj['start_scan']

    assert scanno == start_scan
    assert mzmlobj['precursorList']['count'] == 1

    precursor = mzmlobj['precursorList']['precursor'][0]
    precursor_id = precursor['spectrumRef']
    selected_ion = precursor['selectedIonList']['selectedIon'][0]
    precursor_mz = selected_ion['selected ion m/z']
    precursor_charge = selected_ion['charge state']


    def get_massdiff():
        out = dict()
        for ix in modifications_aa:
            aa = modifications_aa[ix]
            if aa not in mass.std_aa_mass:
                ValueError(f"aa {aa} not in mass.std_aa_mass")
            massdiff = modifications[ix] - mass.std_aa_mass.get(aa, 0)
            out[ix] = massdiff
        return out


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

        # fig, axs = plt.subplots(1, 1, figsize=(10.5, 7))

        fig = plt.figure(figsize=(14, 6))
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[3, 1])
        ax1 = fig.add_subplot(gs[0, 0]) # First row, first column
        ax2 = fig.add_subplot(gs[1, 0]) # Second row, first column
        ax3 = fig.add_subplot(gs[:, 1]) # Both rows, second column
        ax3.axis('off')
        sup.spectrum(msms, ax=ax1, grid=False)
        sup.mass_errors(msms, ax=ax2, plot_unknown=False)
        table = pd.plotting.table(ax3, msms_annot_table, loc='center')
        table.auto_set_font_size(True)
        search_score = str(search_hit['search_score'])
        title = (f"scan {scanno} {proforma_sequence}\n"
                    f"hit rank: {hit_rank} massdiff: {hit_massdiff:4f}\n"
                    f"search score: {search_score}"
                    )
        fig.suptitle(proforma_sequence, fontsize=16, weight='bold', color='black')
        fig.axes[0].title.set_text(title)
        table.set_fontsize(12)


        #fig = sup.facet(
        #    spec_top=msms,

        #    spec_mass_errors=msms,
        #    mass_errors_kws={"plot_unknown": False},
        #    height=7,
        #    width=10.5,
        #)


        plt.tight_layout()
        outname = f"scan_{scanno}_{peptide}_hit_{hit_rank}"
        fulloutname = outpath.joinpath(outname+".png")
        if not os.path.exists(fulloutname):
            plt.savefig(fulloutname, dpi=300, bbox_inches='tight')
        plt.close(fig)


        # ============ altair chart

        achart = supi.spectrum(msms)
        #achart.properties(width=640, height=400)
        achart.width=1200
        achart.height=400
        achart.title = proforma_sequence
        # Define the footnote text
        footnote_text = f"Neutral losses:" + str(annotation_settings['neutral_losses'])

        # Create a text chart for the footnote
        afootnote = alt.Chart(
            {"values": [{"text": footnote_text}]}
        ).mark_text(align='left', baseline='top', fontSize=10).encode(
            text='text:N'
        )

        final_chart = alt.vconcat(achart, afootnote, spacing=5)
        #altair_table = make_table_altair(msms_annot_table)


        outname = outname
        fulloutname = outpath.joinpath(outname+".html")
        if not os.path.exists(fulloutname):
            final_chart.interactive().save(fulloutname)

def make_table_mpl(ax, table_data):
    ax.axis('tight')
    ax.axis('off')
    #table = ax.table(cellText=table_data.values, colLabels=table_data.columns, loc='center', cellLoc='center')
    table = ax.table(cellText=table_data.values, colLabels=table_data.columns, loc='center', cellLoc='center',
    fontsize=12,
                     colWidths=[0.1, 0.1, 0.1])
    return table

def make_table_altair(table_data):
    # not using this right now
    # Convert the DataFrame to a long format
    table_data_long = table_data.melt(var_name='Column', value_name='Text')

    # Creating a "table" using text marks
    text_chart = alt.Chart(table_data_long).mark_text(align='left', baseline='middle', dx=5).encode(
        y=alt.Y('Column:N', axis=alt.Axis(title='')),
        x=alt.X('row_number:O', axis=alt.Axis(title=''), sort=None),
        text='Text:N',
        detail='Column:N'
    ).transform_window(
        row_number='row_number()'
    ).transform_filter(
        alt.datum.Column != 'Column'  # Optional filter if you have a 'Column' column
    )
    return text_chart


import click

@click.command()
@click.option('--file', help='path to the target csv file that contains the survey scan number, fragment scan number, and spectrum file')
def main(file):
    # example
    file = "./Sites_AssigmentVerification.csv"
    df = pd.read_csv(file)
    check_cols(df)


    spec_file_basenames = df.SpectrumFile.unique()
    files = get_files(spec_file_basenames)

    # move out later
    list_of_dicts = df.to_dict(orient='records')


    scans = defaultdict(dict)
    # mzml_info = parse_mzml_files(files['53640_1_EXP_MCF7_EGFRa_LF_phos']['mzml_files'])
    # pepxml_info = parse_mzml_files(files['53640_1_EXP_MCF7_EGFRa_LF_phos']['pepxml_files'])

    for file in files.keys():
        scaninfo = df[ df.SpectrumFile == file ]
        mzml_files = files[file]['mzml_files']
        pepxml_files = files[file]['pepxml_files']
        # TODO fix for if there is more than one file # or maybe not necessary if basefilename is unique by fraction (believe it should be)

        targetscans = scaninfo.FragScanNumber.tolist()
        mzml_file = mzml_files[0]
        pepxml_file = pepxml_files[0]

        logger.info(f"processing {mzml_file} and {pepxml_file}")

        with mzml.MzML(mzml_file.__str__()) as reader:
            for spectrum in tqdm(reader):
                scan_number = spectrum['index'] + 1
                if scan_number in targetscans:
                    print(f"got scan {scan_number}")
                    scans[scan_number]['mzml'] = spectrum

        with pepxml.PepXML(pepxml_file.__str__()) as reader:
            for spectrum in tqdm(reader):
                scan_number = spectrum['start_scan']
                if scan_number in targetscans:
                    print(f"got scan {scan_number}")
                    scans[scan_number]['pepxml'] = spectrum

    for ix, scan in scans.items():
        #break
        handle_scan(scan)







if __name__ == '__main__':
    main()