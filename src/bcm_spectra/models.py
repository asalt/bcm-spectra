# ys(['index', 'id', 'defaultArrayLength', 'scanList', 'precursorList', 'MSn spectrum', 'ms level', 'positive scan', 'centroid spectrum', 'base peak m/z', 'base peak intensity', 'total ion current', 'lowest observed m/z', 'highest observed m/z', 'count', 'm/z array', 'intensity array'])

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


from sqlalchemy import (
    create_engine,
    Column,
    Integer,
    Float,
    String,
    ForeignKey,
    JSON,
    LargeBinary,
)
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import relationship

Base = declarative_base()


class Run(Base):
    __tablename__ = "runs"
    id = Column(Integer, primary_key=True)
    filename = Column(String, nullable=False)
    scans = relationship("Scan", back_populates="run")

    search_results = relationship("SearchResult", back_populates="run")
    precursors = relationship("Precursor", back_populates="run")
    fragments = relationship("Fragment", back_populates="run")


class Scan(Base):
    __tablename__ = "scans"
    id = Column(Integer, primary_key=True)
    scan_number = Column(Integer, nullable=False)

    rt_seconds = Column(Float)

    ms_level = Column(Integer)
    mz_array = Column(LargeBinary)  # this actually a blob we need to define encoding
    intensity_array = Column(
        LargeBinary
    )  # this actually a blob we need to define encoding

    other_info = Column(JSON)

    run_id = Column(Integer, ForeignKey("runs.id"))
    run = relationship("Run", back_populates="scans")

    precursors = relationship("Precursor", back_populates="scan")
    fragments = relationship("Fragment", back_populates="scan")

    # search_results = relationship("SearchResult", back_populates="scan")



class Precursor(Base):
    """
    TODO: this needs to be renamed to SurveyScan
    we are not making an entire table for Precursor selected ion info -
    # that information is stored in the fragments table
    # this is for survey scans
    """
    __tablename__ = "precursors"
    id = Column(Integer, primary_key=True)
    mz = Column(Float)
    intensity = Column(Float)
    charge = Column(Integer)

    scan_id = Column(Integer, ForeignKey("scans.id"))

    scan = relationship("Scan", back_populates="precursors")

    fragments = relationship("Fragment", back_populates="precursors")

    run_id = Column(Integer, ForeignKey("runs.id"))
    run = relationship("Run", back_populates="precursors")


class Fragment(Base):
    __tablename__ = "fragments"
    id = Column(Integer, primary_key=True)
    # mz = Column(Float)
    # intensity = Column(Float)
    precursor_mz = Column(Float)
    precursor_intensity = Column(Float)
    precursor_charge = Column(Integer)

    run_id = Column(Integer, ForeignKey("runs.id"))
    run = relationship("Run", back_populates="fragments")

    scan_id = Column(Integer, ForeignKey("scans.id"))
    scan = relationship("Scan", back_populates="fragments")

    precursor_id = Column(Integer, ForeignKey("precursors.id"))
    precursors = relationship("Precursor", back_populates="fragments")

    search_results = relationship("SearchResult", back_populates="fragment")


class SearchResult(Base):
    __tablename__ = "search_results"
    id = Column(Integer, primary_key=True)
    peptide_sequence = Column(String)
    proforma_sequence = Column(String)
    rank = Column(Integer)

    score = Column(JSON)
    mass_error = Column(Float)
    mass_shifts = Column(JSON)
    mass_diffs = Column(JSON)

    # fragment_id = relationship("Fragment", backref="search_results")

    run_id = Column(Integer, ForeignKey("runs.id"))
    run = relationship("Run", back_populates="search_results")

    # scan_id = Column(Integer, ForeignKey("scans.id"))
    # scan = relationship("Scan", back_populates="search_results")

    fragment_id = Column(Integer, ForeignKey("fragments.id"))
    fragment = relationship("Fragment", back_populates="search_results")
