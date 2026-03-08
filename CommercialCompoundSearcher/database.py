#!/usr/bin/env python3
# coding: utf-8

'''
For initializing and managing a PubChem cache database
'''
import json
import sqlite3
import logging

from pathlib import Path

from .utils import PubchemVendor

VALID_IDENTIFIER_TYPES = {'smiles', 'cas', 'inchikey', 'name'}
VALID_STATUSES = {'ok', 'not_found', 'ambiguous', 'http_error', 'parse_error'}
VALID_ENDPOINTS = {'vendors'}


def connect_pubchem_cache(db_path: Path) -> sqlite3.Connection:
    '''
    Create or open the a PubChem cache database and ensure the expected schema exists.

    Parameters
    ----------
    db_path: Path
        Path to the sqlite database file.

    Returns
    -------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.
    '''
    db_path = Path(db_path)

    # Make sure the output directory exists before opening the database
    db_path.parent.mkdir(parents=True, exist_ok=True)

    # Open the database and configure rows to behave like dictionaries
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row

    # Apply pragmatic SQLite settings for lightweight local caching
    conn.execute('PRAGMA foreign_keys = ON;')
    conn.execute('PRAGMA journal_mode = WAL;')
    conn.execute('PRAGMA synchronous = NORMAL;')

    # Create the minimal unified schema.
    _create_schema(conn)

    # Make sure the expected tables were actually created.
    _validate_schema(conn)

    return conn


def _create_schema(conn: sqlite3.Connection) -> None:
    '''
    Create the unified database schema if it does not already exist.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    Returns
    -------
    None: None
        Creates the required tables if they are missing.
    '''
    # Store compound data by CID
    conn.execute(
        '''
        CREATE TABLE IF NOT EXISTS compounds (
            cid INTEGER PRIMARY KEY,
            canonical_smiles TEXT,
            isomeric_smiles TEXT,
            inchi TEXT,
            inchikey TEXT UNIQUE,
            molecular_formula TEXT,
            iupac_name TEXT,
            last_updated TEXT NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ', 'now'))
        );
        '''
    )

    # Store all non-CID identifiers in one lookup table.
    conn.execute(
        '''
        CREATE TABLE IF NOT EXISTS identifier_lookup (
            identifier_type TEXT NOT NULL,
            identifier_value TEXT NOT NULL,
            cid INTEGER,
            status TEXT NOT NULL CHECK (status IN ('ok', 'not_found', 'ambiguous', 'http_error', 'parse_error')),
            last_checked TEXT NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ', 'now')),
            PRIMARY KEY (identifier_type, identifier_value),
            FOREIGN KEY (cid) REFERENCES compounds(cid)
        );
        '''
    )

    # Store one-to-many vendor data separately from the core compound record.
    conn.execute(
        '''
        CREATE TABLE IF NOT EXISTS vendors (
            cid INTEGER NOT NULL,
            registry_id TEXT,
            sid INTEGER,
            source_detail TEXT,
            vendor_name TEXT NOT NULL,
            source_record_url TEXT,
            source_url TEXT,
            vendor_info_json TEXT,
            last_updated TEXT NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ', 'now')),
            PRIMARY KEY (cid, vendor_name),
            FOREIGN KEY (cid) REFERENCES compounds(cid) ON DELETE CASCADE
        );
        '''
    )

    # Track whether endpoint-specific data has already been fetched.
    conn.execute(
        '''
        CREATE TABLE IF NOT EXISTS endpoint_status (
            cid INTEGER NOT NULL,
            endpoint TEXT NOT NULL,
            status TEXT NOT NULL CHECK (status IN ('ok', 'not_found', 'ambiguous', 'http_error', 'parse_error')),
            last_checked TEXT NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ', 'now')),
            PRIMARY KEY (cid, endpoint),
            FOREIGN KEY (cid) REFERENCES compounds(cid) ON DELETE CASCADE
        );
        '''
    )

    conn.commit()


def _validate_schema(conn: sqlite3.Connection) -> None:
    '''
    Confirm that the expected tables exist in the unified cache database.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    Returns
    -------
    None: None
        Raises an error if a required table is missing.
    '''
    # Query SQLite for all existing table names.
    rows = conn.execute(
        '''
        SELECT name
        FROM sqlite_master
        WHERE type = 'table'
        '''
    ).fetchall()

    table_names = {row['name'] for row in rows}
    expected_table_names = {'compounds', 'identifier_lookup', 'vendors', 'endpoint_status'}

    missing_table_names = expected_table_names - table_names
    if missing_table_names:
        raise RuntimeError(f'Missing expected tables: {sorted(missing_table_names)}')


def _validate_identifier_type(identifier_type: str) -> None:
    '''
    Validate that the identifier type is supported by the lookup table.

    Parameters
    ----------
    identifier_type: str
        Identifier category such as smiles, cas, inchikey, or name.

    Returns
    -------
    None: None
        Raises an error if the identifier type is not supported.
    '''
    if identifier_type not in VALID_IDENTIFIER_TYPES:
        raise ValueError(f'Invalid identifier_type: {identifier_type!r}. Must be one of {sorted(VALID_IDENTIFIER_TYPES)}.')


def _validate_status(status: str) -> None:
    '''
    Validate that the cache status string is supported.

    Parameters
    ----------
    status: str
        Cache status string for an identifier lookup or endpoint fetch.

    Returns
    -------
    None: None
        Raises an error if the status is not supported.
    '''
    if status not in VALID_STATUSES:
        raise ValueError(
            f'Invalid status: {status!r}. Must be one of {sorted(VALID_STATUSES)}.'
        )


def _validate_endpoint(endpoint: str) -> None:
    '''
    Validate that the endpoint label is supported.

    Parameters
    ----------
    endpoint: str
        Endpoint label such as vendors.

    Returns
    -------
    None: None
        Raises an error if the endpoint label is not supported.
    '''
    if endpoint not in VALID_ENDPOINTS:
        raise ValueError(
            f'Invalid endpoint: {endpoint!r}. Must be one of {sorted(VALID_ENDPOINTS)}.'
        )


def upsert_compound(conn: sqlite3.Connection,
                    cid: int,
                    canonical_smiles: str | None = None,
                    isomeric_smiles: str | None = None,
                    inchi: str | None = None,
                    inchikey: str | None = None,
                    molecular_formula: str | None = None,
                    iupac_name: str | None = None) -> None:
    '''
    Insert or update the canonical PubChem compound record for a CID.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    cid: int
        PubChem CID.

    canonical_smiles: str | None
        Canonical SMILES string.

    isomeric_smiles: str | None
        Isomeric SMILES string.

    inchi: str | None
        InChI string.

    inchikey: str | None
        InChIKey string.

    molecular_formula: str | None
        Molecular formula string.

    iupac_name: str | None
        IUPAC name string.

    Returns
    -------
    None: None
        Writes the compound record to the database.
    '''
    # Insert a new compound or update only the fields that were provided.
    conn.execute(
        '''
        INSERT INTO compounds (
            cid,
            canonical_smiles,
            isomeric_smiles,
            inchi,
            inchikey,
            molecular_formula,
            iupac_name,
            last_updated
        )
        VALUES (
            ?, ?, ?, ?, ?, ?, ?, strftime('%Y-%m-%dT%H:%M:%fZ', 'now')
        )
        ON CONFLICT(cid) DO UPDATE SET
            canonical_smiles = COALESCE(excluded.canonical_smiles, compounds.canonical_smiles),
            isomeric_smiles = COALESCE(excluded.isomeric_smiles, compounds.isomeric_smiles),
            inchi = COALESCE(excluded.inchi, compounds.inchi),
            inchikey = COALESCE(excluded.inchikey, compounds.inchikey),
            molecular_formula = COALESCE(excluded.molecular_formula, compounds.molecular_formula),
            iupac_name = COALESCE(excluded.iupac_name, compounds.iupac_name),
            last_updated = strftime('%Y-%m-%dT%H:%M:%fZ', 'now')
        ''',
        (
            cid,
            canonical_smiles,
            isomeric_smiles,
            inchi,
            inchikey,
            molecular_formula,
            iupac_name
        )
    )

    conn.commit()


def get_compound(conn: sqlite3.Connection,
                 cid: int) -> dict | None:
    '''
    Retrieve the canonical compound record for a CID.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    cid: int
        PubChem CID.

    Returns
    -------
    compound: dict | None
        Compound record as a dictionary, or None if no record exists.
    '''
    # Read the compound row and convert it to a plain dictionary.
    row = conn.execute(
        'SELECT * FROM compounds WHERE cid = ?',
        (cid,)
    ).fetchone()

    if row is None:
        return None

    return dict(row)


def upsert_identifier_lookup(conn: sqlite3.Connection,
                             identifier_type: str,
                             identifier_value: str,
                             cid: int | None,
                             status: str) -> None:
    '''
    Insert or update a cached identifier to CID mapping.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    identifier_type: str
        Identifier category such as smiles, cas, inchikey, or name.

    identifier_value: str
        Input identifier value.

    cid: int | None
        PubChem CID if resolved, otherwise None.

    status: str
        Cache status string describing the lookup result.

    Returns
    -------
    None: None
        Writes the identifier lookup record to the database.
    '''
    _validate_identifier_type(identifier_type)
    _validate_status(status)

    # Store the mapping from an arbitrary identifier to a CID.
    conn.execute(
        '''
        INSERT INTO identifier_lookup (
            identifier_type,
            identifier_value,
            cid,
            status,
            last_checked
        )
        VALUES (
            ?, ?, ?, ?, strftime('%Y-%m-%dT%H:%M:%fZ', 'now')
        )
        ON CONFLICT(identifier_type, identifier_value) DO UPDATE SET
            cid = excluded.cid,
            status = excluded.status,
            last_checked = excluded.last_checked
        ''',
        (
            identifier_type,
            identifier_value,
            cid,
            status
        )
    )

    conn.commit()


def get_CID_from_identifier(conn: sqlite3.Connection,
                            identifier_type: str,
                            identifier_value: str) -> tuple[int | None, str | None]:
    '''
    Retrieve a chached CID from an identifier.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    identifier_type: str
        Identifier category such as smiles, cas, inchikey, or name.

    identifier_value: str
        Input identifier value.

    Returns
    -------
    result: tuple[int | None, str | None]
        Cached CID and status, or (None, None) if no record exists.
    '''
    _validate_identifier_type(identifier_type)

    # Read the mapping for the provided identifier.
    row = conn.execute(
        '''
        SELECT cid, status
        FROM identifier_lookup
        WHERE identifier_type = ?
        AND identifier_value = ?
        ''',
        (
            identifier_type,
            identifier_value
        )
    ).fetchone()

    if row is None:
        return None, None

    return row['cid'], row['status']


def get_identifiers_from_CID(conn: sqlite3.Connection,
                             identifier_type: str,
                             cid: int) -> tuple[list[str] | None, str | None]:
    '''
    Retrieve cached identifiers from a CID.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    identifier_type: str
        Identifier category such as smiles, cas, inchikey, or name.

    cid: int
        PubChem CID.

    Returns
    -------
    result: tuple[list[str] | None, str | None]
        Cached identifier values and status, or (None, None) if no record exists.
    '''
    _validate_identifier_type(identifier_type)

    # Read all mappings for the provided CID and identifier type.
    rows = conn.execute(
        '''
        SELECT identifier_value, status
        FROM identifier_lookup
        WHERE identifier_type = ?
        AND cid = ?
        ''',
        (
            identifier_type,
            cid
        )
    ).fetchall()

    if len(rows) == 0:
        return None, None

    identifier_values = [row['identifier_value'] for row in rows]

    # If multiple statuses somehow exist, use the first one.
    status = rows[0]['status']

    return identifier_values, status


def get_cid_from_smiles(conn: sqlite3.Connection,
                        smiles: str) -> tuple[int | None, str | None]:
    '''
    Retrieve a cached CID for a SMILES string.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    smiles: str
        Input SMILES string.

    Returns
    -------
    result: tuple[int | None, str | None]
        Cached CID and status, or (None, None) if no record exists.
    '''
    return get_identifier_lookup(conn=conn,
                                 identifier_type='smiles',
                                 identifier_value=smiles)


def cache_smiles_to_cid(conn: sqlite3.Connection,
                        smiles: str,
                        cid: int | None,
                        status: str) -> None:
    '''
    Cache a SMILES to CID mapping.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    smiles: str
        Input SMILES string.

    cid: int | None
        PubChem CID if resolved, otherwise None.

    status: str
        Cache status string describing the lookup result.

    Returns
    -------
    None: None
        Writes the SMILES lookup record to the database.
    '''
    upsert_identifier_lookup(conn=conn,
                             identifier_type='smiles',
                             identifier_value=smiles,
                             cid=cid,
                             status=status)


def get_cid_from_cas(conn: sqlite3.Connection,
                     cas: str) -> tuple[int | None, str | None]:
    '''
    Retrieve a cached CID for a CAS number.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    cas: str
        Input CAS number.

    Returns
    -------
    result: tuple[int | None, str | None]
        Cached CID and status, or (None, None) if no record exists.
    '''
    return get_identifier_lookup(conn=conn,
                                 identifier_type='cas',
                                 identifier_value=cas)


def cache_cas_to_cid(conn: sqlite3.Connection,
                     cas: str,
                     cid: int | None,
                     status: str) -> None:
    '''
    Cache a CAS to CID mapping.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    cas: str
        Input CAS number.

    cid: int | None
        PubChem CID if resolved, otherwise None.

    status: str
        Cache status string describing the lookup result.

    Returns
    -------
    None: None
        Writes the CAS lookup record to the database.
    '''
    upsert_identifier_lookup(conn=conn,
                             identifier_type='cas',
                             identifier_value=cas,
                             cid=cid,
                             status=status)


def get_smiles_from_cid(conn: sqlite3.Connection,
                        cid: int) -> str | None:
    '''
    Retrieve the cached canonical SMILES string for a CID.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    cid: int
        PubChem CID.

    Returns
    -------
    smiles: str | None
        Cached canonical SMILES string, or None if it is unavailable.
    '''
    # Read the canonical SMILES directly from the compound table.
    row = conn.execute(
        '''
        SELECT canonical_smiles
        FROM compounds
        WHERE cid = ?
        ''',
        (cid,)
    ).fetchone()

    if row is None:
        return None

    return row['canonical_smiles']


def get_smiles_from_cas(conn: sqlite3.Connection,
                        cas: str) -> tuple[str | None, str | None]:
    '''
    Retrieve a cached canonical SMILES string for a CAS number.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    cas: str
        Input CAS number.

    Returns
    -------
    result: tuple[str | None, str | None]
        Cached canonical SMILES and lookup status, or (None, None) if no record exists.
    '''
    # First resolve the CAS number to a CID using the identifier table.
    cid, status = get_cid_from_cas(conn=conn, cas=cas)

    if cid is None:
        return None, status

    # Then read the canonical SMILES from the compound table.
    smiles = get_smiles_from_cid(conn=conn, cid=cid)
    return smiles, status


def replace_vendors_for_cid(conn: sqlite3.Connection,
                            cid: int,
                            vendors: list[PubchemVendor],
                            status: str = 'ok') -> None:
    '''
    Replace the cached vendor list for a CID using PubchemVendor objects.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    cid: int
        PubChem CID.

    vendors: list[PubchemVendor]
        List of PubchemVendor objects returned for the CID.

    status: str
        Cache status string describing the vendor fetch result.

    Returns
    -------
    None: None
        Replaces the cached vendor list and updates endpoint fetch status.
    '''
    _validate_status(status)

    # Make sure the parent compound row exists before storing child rows.
    conn.execute(
        '''
        INSERT INTO compounds (cid, last_updated)
        VALUES (?, strftime('%Y-%m-%dT%H:%M:%fZ', 'now'))
        ON CONFLICT(cid) DO NOTHING
        ''',
        (cid,)
    )

    # Remove the old vendor list so the cache always reflects the latest fetch.
    conn.execute(
        '''
        DELETE FROM vendors
        WHERE cid = ?
        ''',
        (cid,)
    )

    # Insert the full vendor payload for each cached vendor.
    for vendor in vendors:
        conn.execute(
            '''
            INSERT INTO vendors (
                cid,
                vendor_name,
                registry_id,
                sid,
                source_detail,
                source_record_url,
                source_url,
                vendor_info_json,
                last_updated
            )
            VALUES (
                ?, ?, ?, ?, ?, ?, ?, ?, strftime('%Y-%m-%dT%H:%M:%fZ', 'now')
            )
            ''',
            (
                cid,
                vendor.SourceName,
                vendor.RegistryID,
                vendor.SID,
                vendor.SourceDetail,
                vendor.SourceRecordURL,
                vendor.SourceURL,
                json.dumps(vendor.vendor_info)
            )
        )

    # Record that the vendors endpoint has been checked for this CID.
    conn.execute(
        '''
        INSERT INTO endpoint_status (
            cid,
            endpoint,
            status,
            last_checked
        )
        VALUES (
            ?, 'vendors', ?, strftime('%Y-%m-%dT%H:%M:%fZ', 'now')
        )
        ON CONFLICT(cid, endpoint) DO UPDATE SET
            status = excluded.status,
            last_checked = excluded.last_checked
        ''',
        (
            cid,
            status
        )
    )

    conn.commit()


def get_vendors_from_cid(conn: sqlite3.Connection,
                         cid: int) -> tuple[list[PubchemVendor] | None, str | None]:
    '''
    Retrieve the cached vendor list for a CID as PubchemVendor objects.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    cid: int
        PubChem CID.

    Returns
    -------
    result: tuple[list | None, str | None]
        Cached list of PubchemVendor objects and endpoint status, or
        (None, None) if vendors were never fetched.
    '''
    _validate_endpoint('vendors')

    # Read the endpoint status first so empty vendor lists can still be cached.
    status_row = conn.execute(
        '''
        SELECT status
        FROM endpoint_status
        WHERE cid = ?
        AND endpoint = 'vendors'
        ''',
        (cid,)
    ).fetchone()

    if status_row is None:
        return None, None

    # Read all cached vendor fields for the CID.
    rows = conn.execute(
        '''
        SELECT
            vendor_name,
            registry_id,
            sid,
            source_detail,
            source_record_url,
            source_url,
            vendor_info_json
        FROM vendors
        WHERE cid = ?
        ORDER BY vendor_name
        ''',
        (cid,)
    ).fetchall()

    # Reconstruct PubchemVendor objects from the cached rows.
    vendor_objects = []
    for row in rows:
        if row['vendor_info_json'] is not None:
            vendor_info = json.loads(row['vendor_info_json'])
        else:
            vendor_info = {
                'RegistryID': row['registry_id'],
                'SID': row['sid'],
                'SourceDetail': row['source_detail'],
                'SourceName': row['vendor_name'],
                'SourceRecordURL': row['source_record_url'],
                'SourceURL': row['source_url']
            }

        vendor_objects.append(PubchemVendor(cid, vendor_info))

    return vendor_objects, status_row['status']


def close_pubchem_cache(conn: sqlite3.Connection) -> None:
    '''
    Commit any pending work and close the unified cache database.

    Parameters
    ----------
    conn: sqlite3.Connection
        Open SQLite connection to the unified cache database.

    Returns
    -------
    None: None
        Closes the SQLite connection.
    '''
    # Make sure pending writes are flushed before closing.
    conn.commit()
    conn.close()
