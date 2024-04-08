import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import numpy as np

from bcm_spectra.models import Base, Run, Scan, Precursor, SearchResult


@pytest.fixture
def session():
    # engine = create_engine('sqlite:///:memory:')  # Use an in-memory SQLite database for tests
    engine = create_engine(
        "sqlite:///test_db.db"
    )  # Use an in-memory SQLite database for tests
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)()
    yield Session
    Session.close()


@pytest.fixture()
def run(session):
    run = Run(filename="test_run")
    session.add(run)
    session.commit()
    return run


def test_scan_insertion(session):
    run = Run(filename="test_run")
    scan = Scan(scan_number=123, ms_level=2, run=run)
    session.add(scan)
    session.commit()
    assert session.query(Run).filter_by(filename="test_run").first() is not None


def test_binary_store(session, run):

    # filename = "test_binary_store"
    # run = Run(filename=filename)
    filename = run.filename

    assert run is not None
    # we already have session and run
    # now we add scan to the run

    mz_array = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    intensity_array = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    scan = Scan(
        scan_number=123,
        ms_level=2,
        run=run,
        mz_array=mz_array.tobytes(),
        intensity_array=intensity_array.tobytes(),
    )
    session.add(scan)
    session.commit()

    # run = session.query(Run).filter_by(filename=filename).first()

    scans = run.scans

    assert len(scans) == 1

    scan = scans[0]

    mz_array_decoded = np.frombuffer(scan.mz_array, dtype=np.float32)

    with pytest.raises(ValueError):
        mz_array_decoded_wrong = np.frombuffer(scan.mz_array, dtype=np.float64)

    intensity_array_decoded = np.frombuffer(scan.intensity_array, dtype=np.float64)

    assert np.allclose(mz_array, mz_array_decoded)
