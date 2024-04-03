import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from bcm_spectra.db import init_db, get_session

def test_init_db():
    engine = init_db('sqlite:///:memory:')
    assert engine is not None

def test_get_session():
    engine = create_engine('sqlite:///:memory:')
    session = get_session(engine)
    assert session is not None
    assert session.bind == engine
    session.close()


def test_get_session_with_engine():
    engine = create_engine('sqlite:///:memory:')
    session = get_session(engine)
    assert session is not None
    assert session.bind == engine
    session.close()
