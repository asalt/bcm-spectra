
import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from bcm_spectra.models import Base, Run, Scan, Precursor, SearchResult

@pytest.fixture
def session():
    engine = create_engine('sqlite:///:memory:')  # Use an in-memory SQLite database for tests
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    return Session()

def test_scan_insertion(session):
    run = Run(filename="test_run")
    scan = Scan(scan_number=123, ms_level=2, run=run)
    session.add(scan)
    session.commit()
    assert session.query(Run).filter_by(filename="test_run").first() is not None
