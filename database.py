#!/usr/bin/env python3
# coding: utf-8

'''
Functions for initializing and managing the cache databases.
'''

from __future__ import annotations

import time
import sqlite3

from pathlib import Path

from contextlib import contextmanager

from typing import Generator

import sqlite3
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def init_smiles_cache(db_path: Path) -> Generator[sqlite3.Connection]:
    '''
    Create or open a SQLite CID to SMILES cache, enforce the expected schema, and yield the connection.

    Parameters
    ----------
    db_path: Path
        Path to the sqlite database file.

    Returns
    -------
    conn: sqlite3.Connection
        Open SQLite connection to the cache database.
    '''
    db_path = Path(db_path)

    # Make the parent directories if they do not exist
    db_path.parent.mkdir(parents=True, exist_ok=True)

    # Open/create database
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row

    try:
        # Set pragmatic defaults for notebook usage
        conn.execute('PRAGMA foreign_keys = ON;')
        conn.execute('PRAGMA journal_mode = WAL;')
        conn.execute('PRAGMA synchronous = NORMAL;')

        # Create schema if missing
        conn.execute(
            '''
            CREATE TABLE IF NOT EXISTS cid_smiles (
                cid INTEGER PRIMARY KEY,
                smiles TEXT NOT NULL,
                retrieved_at TEXT NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ', 'now'))
            );
            '''
        )

        # Validate/enforce schema for existing DBs
        info = conn.execute("PRAGMA table_info('cid_smiles');").fetchall()
        if not info:
            raise RuntimeError("Schema error: table 'cid_smiles' was not created and does not exist.")

        col = {row['name']: row for row in info}
        expected_cols = {'cid', 'smiles', 'retrieved_at'}
        if set(col.keys()) != expected_cols:
            raise RuntimeError(
                f"Schema mismatch for 'cid_smiles'. Expected columns {sorted(expected_cols)}, "
                f"found {sorted(col.keys())}. Delete the DB or run a migration."
            )

        # cid: INTEGER PRIMARY KEY
        if col['cid']['type'].upper() != 'INTEGER' or col['cid']['pk'] != 1:
            raise RuntimeError("Schema mismatch: 'cid' must be INTEGER PRIMARY KEY.")

        # smiles: TEXT NOT NULL
        if col['smiles']['type'].upper() != 'TEXT' or col['smiles']['notnull'] != 1:
            raise RuntimeError("Schema mismatch: 'smiles' must be TEXT NOT NULL.")

        # retrieved_at: TEXT NOT NULL with default timestamp expression
        retrieved_default = col['retrieved_at']['dflt_value']
        if col['retrieved_at']['type'].upper() != 'TEXT' or col['retrieved_at']['notnull'] != 1:
            raise RuntimeError("Schema mismatch: 'retrieved_at' must be TEXT NOT NULL.")

        if not retrieved_default or 'strftime' not in str(retrieved_default).lower():
            raise RuntimeError(
                "Schema mismatch: 'retrieved_at' must have a DEFAULT strftime(...) timestamp."
            )

        conn.commit()
        yield conn
    finally:
        conn.close()


def get_cached_smiles(conn: sqlite3.Connection,
                      cid: int) -> tuple:
    '''
    Retrieve a cached SMILES entry.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open sqlite connection.

    cid: int
        PubChem CID.

    Returns
    -------
    result: tuple
        (smiles, status) if present, else (None, None).
    '''
    cur = conn.execute('SELECT smiles, status FROM smiles_cache WHERE cid=?', (cid,))
    row = cur.fetchone()
    if row is None:
        return None, None
    return row[0], row[1]


def upsert_cached_smiles(conn: sqlite3.Connection, cid: int, smiles: str, status: str) -> None:
    '''
    Insert or update a cached SMILES entry.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open sqlite connection.

    cid: int
        PubChem CID.

    smiles: str
        SMILES string (can be empty or None-like string, depending on policy).

    status: str
        One of: 'ok', 'missing', 'http_error', 'parse_error'.

    Returns
    -------
    None: None
        Writes to the cache table.
    '''
    conn.execute(
        '''
        INSERT INTO smiles_cache (cid, smiles, status, updated_at)
        VALUES (?, ?, ?, ?)
        ON CONFLICT(cid) DO UPDATE SET
            smiles=excluded.smiles,
            status=excluded.status,
            updated_at=excluded.updated_at
        ''',
        (cid, smiles, status, time.time())
    )