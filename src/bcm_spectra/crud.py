from sqlalchemy.orm import Session

from .db import init_db, get_session
from .models import Base, Run, Scan, Precursor, SearchResult


def add_and_commit(db: Session, obj):
    db.add(obj)
    db.commit()
    db.refresh(obj)
    return obj


def get_run_by_filename(session: Session, filename):
    run = session.query(Run).filter_by(filename=filename).first()
    session.close()
    return run


def get_or_create_run(session: Session, filename):
    run = get_run_by_filename(session, filename)
    if run is None:
        run = Run(filename=filename)
        run = add_and_commit(session, run)
    return run


def get_scan_from_run_by_scan_number(session: Session, run: Run, scan_number):
    scan = session.query(Scan).filter_by(run=run, scan_number=scan_number).first()
    return scan


def add_scan_to_run(
    session: Session,
    run: Run,
    scan_number,
    ms_level,
    mz_array,
    intensity_array,
    other_info,
):
    scan = Scan(
        scan_number=scan_number,
        ms_level=ms_level,
        mz_array=mz_array,
        intensity_array=intensity_array,
        other_info=other_info,
        run=run,
    )
    scan = add_and_commit(session, scan)
    return scan
