# db.py

#ys(['index', 'id', 'defaultArrayLength', 'scanList', 'precursorList', 'MSn spectrum', 'ms level', 'positive scan', 'centroid spectrum', 'base peak m/z', 'base peak intensity', 'total ion current', 'lowest observed m/z', 'highest observed m/z', 'count', 'm/z array', 'intensity array'])

#   mzarray = mzmlobj['m/z array']
#  intensityarray = mzmlobj['intensity array']
#  _name = mzmlobj['id'] # concat other info to make more comprehensive
#  scanno = mzmlobj['index'] + 1 # concat other info to make more comprehensive
#  start_scan = pepxmlobj['start_scan']
#
#  assert scanno == start_scan
#  assert mzmlobj['precursorList']['count'] == 1
#
#  precursor = mzmlobj['precursorList']['precursor'][0]
#  precursor_id = precursor['spectrumRef']
#  selected_ion = precursor['selectedIonList']['selectedIon'][0]
#  precursor_mz = selected_ion['selected ion m/z']
#  precursor_charge = selected_ion['charge state']


# >>> spectrum['defaultArrayLength']
# 120
# >>> spectrum['scanList']
# {'count': 1, 'scan': [{'scanWindowList': {'count': 1, 'scanWindow': [{'scan window lower limit': 120.0, 'scan window upper limit': 1415.516357421875}]}, 'scan start time': 41.083011782583, 'mass resolving power': 15000.0, 'filter string': 'FTMS + p NSI d Full ms2 683.2767@hcd32.00 [120.0000-1415.5164]', 'preset scan configuration': 2.0, 'ion injection time': 29.999999329448, '[Thermo Trailer Extra]Monoisotopic M/Z:': 683.2766581365408}], 'no combination': ''}
# i need  ms level
# two tables. one for ms level 1 and one for ms level 2.
# ms level 2 will link forgein key to ms level 1
# both have mz array and intensity array
# these are binary blobs we can decode with numpy.frombuffer?? or something
# it is done within pyteomics library code somewhere perhaps

# then we have pepxml table. which is a 1 to many relationship from mzml level 2 and each pepxml record.
