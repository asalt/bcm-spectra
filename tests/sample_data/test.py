from pyteomics import mzml
import os

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

interested_scans = [f'controllerType=0 controllerNumber=1 scan={scan}' for scan in scans]

# Input and output file paths
input_mzml = "./sample_data/53640_1_EXP_MCF7_EGFRa_LF_phos.mzML"
output_mzml = 'short.mzML'

from psims.transform.mzml import MzMLTransformer, cvstr

def myfilter(spectrum):
    if spectrum.get('id') in interested_scans:
        return spectrum
    return None

with open(input_mzml, 'rb') as in_stream, open(output_mzml, 'wb') as out_stream:
    transformer = MzMLTransformer(in_stream, out_stream, myfilter).write()

# this isn't worty doing right now
# import psims
# writer = psims.mzml.MzMLWriter(output_mzml)
# writer.begin()
# writer.file_description([
#     # "MS1 spectrum",
#     "MS2 spectrum",
#     "centroid spectrum",
# ])
# writer.controlled_vocabularies()
# writer.run(id=1)




# os.setwd(os.path.dirname(__file__))

reader = mzml.read(input_mzml)
for spectrum in reader:
    if spectrum.get('id') in interested_scans:
        mz_array = spectrum['m/z array']
        intensity_array = spectrum['intensity array']
n = next(reader)


# Reading from the original mzML and writing selected scans to a new file
with mzml.read(input_mzml) as reader, mzml.write(output_mzml, use_index=True) as writer:

print(f"Abbreviated mzML file written to {os.path.abspath(output_mzml)}")
