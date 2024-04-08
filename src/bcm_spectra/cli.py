# main.py

import os
from collections import defaultdict
import argparse
import logging

# import glob
from pathlib import Path
import pandas as pd
from pyteomics import mzml, pepxml
from pyteomics import mass
import spectrum_utils
import click

import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup
import spectrum_utils.iplot as supi
from tqdm import tqdm
import altair as alt


from .utils import get_files, get_filescans
from . import db
from . import utils
from . import parser
from . import plot


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def check_cols(df):
    assert "SurveyScanNumber" in df.columns
    assert "FragScanNumber" in df.columns
    assert "SpectrumFile" in df.columns
    return df


annotation_settings = {
    "fragment_tol_mass": 0.05,
    "fragment_tol_mode": "Da",
    "ion_types": "by",
    "max_ion_charge": 2,
    "neutral_losses": {"NH3": -17.026549, "H2O": -18.010565, "H3PO4": -97.976896},
}


# def get_scan_info(df, scan_numbers)


def handle_scan(scan: dict):

    mzmlobj = scan.get("mzml")
    pepxmlobj = scan.get("pepxml")

    if not mzmlobj["ms level"] == 2:
        print(f"scan['id'] not msms")
        return

    mzarray = mzmlobj["m/z array"]
    intensityarray = mzmlobj["intensity array"]
    _name = mzmlobj["id"]  # concat other info to make more comprehensive
    scanno = mzmlobj["index"] + 1  # concat other info to make more comprehensive
    start_scan = pepxmlobj["start_scan"]

    assert scanno == start_scan
    assert mzmlobj["precursorList"]["count"] == 1

    precursor = mzmlobj["precursorList"]["precursor"][0]
    precursor_id = precursor["spectrumRef"]
    selected_ion = precursor["selectedIonList"]["selectedIon"][0]
    precursor_mz = selected_ion["selected ion m/z"]
    precursor_charge = selected_ion["charge state"]

    def get_massdiff():
        out = dict()
        for ix in modifications_aa:
            aa = modifications_aa[ix]
            if aa not in mass.std_aa_mass:
                ValueError(f"aa {aa} not in mass.std_aa_mass")
            massdiff = modifications[ix] - mass.std_aa_mass.get(aa, 0)
            out[ix] = massdiff
        return out

    search_hits = pepxmlobj["search_hit"]

    for search_hit in search_hits:
        # this has many many little ways to get tricked up. so many things to map

        # dict_keys(['peptide', 'num_missed_cleavages', 'num_tot_proteins', 'tot_num_ions', 'hit_rank', 'num_matched_ions', 'search_score', 'modified_peptide', 'massdiff', 'calc_neutral_pep_mass', 'is_rejected', 'proteins', 'modifications'])
        hit_rank = search_hit["hit_rank"]
        hit_massdiff = search_hit["massdiff"]

        peptide = search_hit["peptide"]
        modifications = {
            mod["position"]: mod["mass"] for mod in search_hit["modifications"]
        }
        modifications_aa = {
            k: peptide[k - 1] for k in modifications.keys()
        }  # 1 based index
        massdiff_dict = get_massdiff()

        proforma_sequence = ""
        for i, aa in enumerate(peptide, start=1):
            proforma_sequence += aa
            if i in massdiff_dict:
                # Add the modification in square brackets
                mass_diff = massdiff_dict[i]
                proforma_sequence += f"[+{mass_diff:.4f}]"

        msms = sus.MsmsSpectrum(
            identifier=_name,
            precursor_mz=float(precursor_mz),
            precursor_charge=precursor_charge,
            mz=mzarray,
            intensity=intensityarray,
            retention_time=1.0,
            # peptide = peptide,
            # modifications = modifications
        )

        msms.annotate_proforma(proforma_sequence, **annotation_settings)

        # fig, ax = plt.subplots()
        # sup.spectrum(msms, ax=ax, grid=False)
        msms_annot_table = extract_ion_annotation(msms)

        outpath = Path("./spectra")
        outpath.mkdir(exist_ok=True)

        # fig, axs = plt.subplots(1, 1, figsize=(10.5, 7))

        fig = plt.figure(figsize=(14, 6))
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[3, 1])
        ax1 = fig.add_subplot(gs[0, 0])  # First row, first column
        ax2 = fig.add_subplot(gs[1, 0])  # Second row, first column
        ax3 = fig.add_subplot(gs[:, 1])  # Both rows, second column
        ax3.axis("off")
        sup.spectrum(msms, ax=ax1, grid=False)
        sup.mass_errors(msms, ax=ax2, plot_unknown=False)
        table = pd.plotting.table(ax3, msms_annot_table, loc="center")
        table.auto_set_font_size(True)
        search_score = str(search_hit["search_score"])
        title = (
            f"scan {scanno} {proforma_sequence}\n"
            f"hit rank: {hit_rank} massdiff: {hit_massdiff:4f}\n"
            f"search score: {search_score}"
        )
        fig.suptitle(proforma_sequence, fontsize=16, weight="bold", color="black")
        fig.axes[0].title.set_text(title)
        table.set_fontsize(12)

        # fig = sup.facet(
        #    spec_top=msms,

        #    spec_mass_errors=msms,
        #    mass_errors_kws={"plot_unknown": False},
        #    height=7,
        #    width=10.5,
        # )

        plt.tight_layout()
        outname = f"scan_{scanno}_{peptide}_hit_{hit_rank}"
        fulloutname = outpath.joinpath(outname + ".png")
        if not os.path.exists(fulloutname):
            plt.savefig(fulloutname, dpi=300, bbox_inches="tight")
        plt.close(fig)

        # ============ altair chart

        achart = supi.spectrum(msms)
        # achart.properties(width=640, height=400)
        achart.width = 1200
        achart.height = 400
        achart.title = proforma_sequence
        # Define the footnote text
        footnote_text = f"Neutral losses:" + str(annotation_settings["neutral_losses"])

        # Create a text chart for the footnote
        afootnote = (
            alt.Chart({"values": [{"text": footnote_text}]})
            .mark_text(align="left", baseline="top", fontSize=10)
            .encode(text="text:N")
        )

        final_chart = alt.vconcat(achart, afootnote, spacing=5)
        # altair_table = make_table_altair(msms_annot_table)

        outname = outname
        fulloutname = outpath.joinpath(outname + ".html")
        if not os.path.exists(fulloutname):
            final_chart.interactive().save(fulloutname)




@click.command()
@click.option(
    "--data-dir", help="data directory where the mzML and pepXML files are located"
)
@click.option(
    "--file",
    help="path to the target csv file that contains the survey scan number, fragment scan number, and spectrum file",
)
@click.option("--sqlite-db", help="path to the sqlite database", default="scans.db")
def main(
    data_dir,
    file,
    sqlite_db,
):
    """
    :param file: path to the target csv file that contains the survey scan number, fragment scan number, and spectrum file
    percisely the columns are 'SurveyScanNumber', 'FragScanNumber', 'SpectrumFile'
    """
    # example
    # file = "./Sites_AssigmentVerification.csv"
    df = pd.read_csv(file)
    df = df.head(1)
    check_cols(df)

    # move out later
    list_of_dicts = df.to_dict(orient="records")

    spec_file_basenames = df.SpectrumFile.unique()

    if sqlite_db is None:
        raise NotImplementedError("sqlite db required")
    if sqlite_db:
        engine_url = f"sqlite:///{sqlite_db}"
        session = db.get_global_session(engine_url)

    files = utils.get_files(spec_file_basenames, basepath=data_dir)
    # scans = defaultdict(dict)
    filescans = utils.get_all_filescans(files, df, session=session)



    # this can be moved to test
    assert all( not isinstance(q, dict) for v in filescans.values() for q in v.values() )


    for ix, scan_mapping in filescans.items():
        for scanno, scan in scan_mapping.items():
            1+1
            # do something
            # handle_scan(scan)


            fragments = scan.fragments
            if len(fragments) == 0:
                logger.warning(f"no fragments found for scan {scanno}")
                continue
            fragment = fragments[0]
            mz_array = np.frombuffer(fragment.scan.mz_array, dtype=np.float64)
            intensity_array = np.frombuffer(fragment.scan.mz_array, dtype=np.float64)

            assert len(mz_array) == len(intensity_array)

            _name = f"{fragment.run.filename}_{fragment.scan.scan_number}_{fragment.precursor_charge}"

            rt = scan.rt_seconds / 60

            msms = sus.MsmsSpectrum( identifier=_name, precursor_mz=fragment.precursor_mz,
                                     precursor_charge=fragment.precursor_charge, mz=mz_array, intensity=intensity_array,
                 retention_time=rt,
            )

                # peptide = peptide,
                # modifications = modifications
            #)
            import ipdb; ipdb.set_trace()

            for search_result in fragment.search_results:
                proforma_sequence = search_result.proforma_sequence
                hit_rank = search_result.rank
                scoreinfo = search_result.score # dict
                # the "name" of the score may not be constant - not sure about others
                search_score = scoreinfo.get("hyperscore", )
                nextscore = scoreinfo.get("nextscore", )
                mass_error = search_result.mass_error

                info = (
                    f"scan {scanno} {proforma_sequence}\n"
                    f"hit rank: {hit_rank} massdiff: {mass_error:4f}\n"
                    f"search score: {search_score}"
                )
                msms.annotate_proforma(proforma_sequence, **annotation_settings)
                fig = plot.make_spectrum_plot(msms, proforma_sequence=proforma_sequence,
                                         additional_info=info)

    # scanobj = filescans[ list(filescans.keys())[0] ].get(43816)
    # f1 = scanobj.fragments[0]

    # all the scans are in this filescans dictionary object
    # the structure is
    # filescans = {
    #     '53640_1_EXP_MCF7_EGFRa_LF_phos':
    #         {'mzml': mzml_obj, 'pepxml': pepxml_obj},
    #         ...
    #     },
    #     '53640_2_EXP_MCF7_EGFRa_LF_phos':
    #         {'mzml': mzml_obj, 'pepxml': pepxml_obj},
    #         ...
    #     },
    #     ...
    # } # double check if this is right. need to make tests
    # for each scan in each file, we can now handle the scan
    # for example

    # for ix, scan in filescans.items():
    #     #break
    #     handle_scan(scan)


if __name__ == "__main__":
    main()
