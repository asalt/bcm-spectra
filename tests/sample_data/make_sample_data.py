from pyteomics import mzml
import os

# List of scan numbers you are interested in
interested_scans = ["scan=1", "scan=2"]  # Adjust this list to your scan numbers

# first scans
scans = [
    43816,
    19233,
    19233,
    108782,
    108375,
    107908,
    107933,
    115674,
    116353,
    116792,
    57186,
    58595,
    58221,
    57863,
    71341,
    71692,
    58339,
    57979,
    58717,
    57601,
    72067,
    71975,
    56882,
    67181,
    59752,
    59345,
    58623,
    58271,
    57895,
    60952,
    61508,
    58985,
    63738,
]

# Input and output file paths
input_mzml = "./53640_1_EXP_MCF7_EGFRa_LF_phos.mzML"
output_mzml = "short.mzML"


# rr = mzml.read(input_mzml)

# Reading from the original mzML and writing selected scans to a new file
with mzml.read(input_mzml) as reader, mzml.write(output_mzml, use_index=True) as writer:
    for spectrum in reader:
        if spectrum.get("index", 0) + 1 in interested_scans:
            writer.write(spectrum)

print(f"Abbreviated mzML file written to {os.path.abspath(output_mzml)}")
