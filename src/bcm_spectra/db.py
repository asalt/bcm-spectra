# db.py

from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session
from sqlalchemy.orm import sessionmaker

from .models import Base

# engine = create_engine('sqlite:///ms_data.db')
# Base.metadata.create_all(engine)

def init_db(db_url):
    engine = create_engine(db_url)
    Base.metadata.create_all(engine)
    return engine

# Function to get a new session
def get_session(engine):
    # Using scoped_session to ensure thread safety in a web app context or similar
    session_factory = sessionmaker(bind=engine)
    Session = scoped_session(session_factory)
    return Session

# Function to provide a global point of access to the session
def get_global_session(db_url):
    engine = init_db(db_url)
    return get_session(engine)