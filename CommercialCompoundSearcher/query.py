#!/usr/bin/env python3
# coding: utf-8

'''
For querying PubChem
'''

from __future__ import annotations

import json
import urllib
import logging
import urllib.request

from .utils import PubchemVendor

logger = logging.getLogger(__name__)


def query_cid_from_inchi_key(inchi_key: str):
    '''
    Gets the CID from Pubchem using the PUG REST API.

    Parameters
    ----------
    inchi_key: str
        WRITE DOC STRING

    Returns
    ----------
    str
    '''
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchi_key}/cids/JSON/'
    with urllib.request.urlopen(url) as u:
        data = u.read()
    data = data.decode('utf-8')
    j = json.loads(data)

    return j['IdentifierList']['CID'][0]


def query_CID_from_CAS(cas: str):
    '''
    Gets the CID from Pubchem using the PUG REST API.

    Parameters
    ----------
    cas: str
        WRITE DOC STRING

    Returns
    ----------
    str
    '''
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas}/cids/JSON/'
    with urllib.request.urlopen(url) as u:
        data = u.read()
    data = data.decode('utf-8')
    j = json.loads(data)

    return j['IdentifierList']['CID'][0]


def query_vendors_from_cid(cid: int) -> list[PubchemVendor]:
    '''
    Gets a list of vendors from Pubchem based on
    the CID using the PUG REST API

    Returns empty list if no vendors are found
    '''
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/{cid}/JSON/?response_type=view'
    with urllib.request.urlopen(url) as u:
        data = u.read()
    data = data.decode('utf-8')
    j = json.loads(data)

    chemical_vendors = None
    for category in j['SourceCategories']['Categories']:
        if category['Category'].casefold() == 'Chemical Vendors'.casefold():
            chemical_vendors = category
            break

    if chemical_vendors is None:
        print(f'CID {cid} does not have vendors')
        return []

    vendors = [PubchemVendor(cid, x) for x in chemical_vendors['Sources']]

    return vendors


def get_vendor_json(cid: int) -> list:
    '''
    Gets a list of vendors from Pubchem based on
    the CID using the PUG REST API

    Returns empty list if no vendors are found
    '''
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/{cid}/JSON/?response_type=view'
    with urllib.request.urlopen(url) as u:
        data = u.read()
    data = data.decode('utf-8')
    j = json.loads(data)

    chemical_vendors = None
    for category in j['SourceCategories']['Categories']:
        if category['Category'].casefold() == 'Chemical Vendors'.casefold():
            chemical_vendors = category
            break

    if chemical_vendors is None:
        print(f'CID {cid} does not have vendors')
        return []

    results = []
    for d in chemical_vendors['Sources']:
        results.append(PubchemVendor(d))
    return results


def query_CAS_from_cid(cid: int) -> str | None:
    '''
    Gets a CAS number from Pubchem based on
    the CID using the PUG REST API
    '''

    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON/'

    with urllib.request.urlopen(url) as u:
        data = u.read()
    data = data.decode('utf-8')
    j = json.loads(data)
    CAS = [x for x in j['InformationList']['Information'][0]['Synonym'] if len(''.join([z for z in x if z == '-'])) == 2]
    CAS = [x for x in CAS if all([z.isdigit() or z == '-' for z in x])]

    if len(CAS) == 0:
        return None

    if len(CAS) != 1:
        logger.warning('Found %d CAS number for CID %d. %s', len(CAS), cid, str(CAS))

    return str(CAS[0])


def query_SMILES_from_CID(cid: int) -> str | None:
    '''
    Gets the SMILES string from Pubchem based on the
    CID using the PUG REST API

    Returns None if no smiles was found
    '''

    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON/'
    with urllib.request.urlopen(url) as u:
        data = u.read()
    data = data.decode('utf-8')
    j = json.loads(data)

    try:
        SMILES = j['PropertyTable']['Properties'][0]['ConnectivitySMILES']
    except KeyError as e:
        print(f'Warning: ConnectivitySMILES keyerror. Using CanonicalSMILES')
        try:
            SMILES = j['PropertyTable']['Properties'][0]['ConnectivitySMILES']
        except KeyError as e:
            print(f'Could not get SMILES for {cid}')
            return None

    if len(SMILES) == 0:
        return None

    return str(SMILES)


def get_SMILES_from_CAS(cas: int) -> str | None:
    '''
    Gets a canonical SMILES string from Pubchem based
    on the CAS number using the PUG REST API

    Returns empty string if no SMILES is found
    '''
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas}/property/MolecularFormula,CanonicalSMILES/JSON'
    #print(url)
    with urllib.request.urlopen(url) as u:
        data = u.read()
    data = data.decode('utf-8')
    j = json.loads(data)

    try:
        SMILES = j['PropertyTable']['Properties'][0]['ConnectivitySMILES']
    except KeyError as e:
        print(f'Warning: ConnectivitySMILES keyerror. Using CanonicalSMILES')
        try:
            SMILES = j['PropertyTable']['Properties'][0]['ConnectivitySMILES']
        except KeyError as e:
            print(f'Could not get SMILES for {cas}')
            return None

    if len(SMILES) == 0:
        return None

    return str(SMILES)
