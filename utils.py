#!/usr/bin/env python3
# coding: utf-8

'''
Functions for processing
commercially available compounds
'''

from __future__ import annotations

import re
import ast
import json
import urllib
import urllib.request

from collections import Counter

import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from pathlib import Path

EXAMPLE_VENDORS_TO_REMOVE = ['Chemieliva Pharmaceutical Co., Ltd',
 'Elsa Biotechnology',
 'AmicBase - Antimicrobial Activities',
 'BOC Sciences',
 'Parchem',
 'Chemenu Inc.',
 'Cyclic PharmaTech',
 'Syntech Labs',
 'BLD Pharm',
 'Twinhang OU',
 'ChemShuttle',
 'Anward',
 'Shanghai Send Pharmaceutical Technology Co., Ltd',
 'Tenova Pharmaceuticals Inc',
 'Sun-shine Chemical',
 'MIC Scientific',
 '3A SpeedChemical Inc',
 'Syntechem',
 'B&C Chemical',
 'AbaChemScene',
 'Calbiochem',
 'Tractus',
 'CAPOT',
 'Vitas-M Laboratory',
 'Chemhere',
 'MolPort',
 'Excenen Pharmatech',
 'Paragos',
 'ApexBio Technology',
 'Hangzhou APIChem Technology',
 'Tetrahedron Scientific Inc',
 'Moldb',
 'AstaTech, Inc.',
 'T&J Chemicals (Singapore)',
 'Allbio Pharm Co., Ltd',
 'ChemExper Chemical Directory',
 'CEGChem',
 'Finetech Industry Limited',
 'Chembase.cn',
 'ChemMol',
 'Axon Medchem',
 'Aurora Fine Chemicals LLC',
 '3WAY PHARM INC',
 'DAOGE BIOPHARMA',
 'Hairui Chemical',
 'Alsachim',
 'eInhibitors',
 'SynHet - Synthetic Heterocycles',
 '001Chemical',
 'Watec Laboratories',
 'Specs',
 'Chemchart',
 'Win-Win Chemical',
 'Biocore',
 'Strem Chemicals, Inc.',
 'eNovation Chemicals',
 'Acemol',
 'Aurum Pharmatech LLC',
 'SPECIFIC POLYMERS (SP)',
 'ET Co.,Ltd.',
 'Alomone Labs',
 'Chiralblock Biosciences',
 'ISpharm',
 'Amatye',
 'Acadechem',
 'Mcule',
 'Sinfoo Biotech',
 'CymitQuimica',
 'Bic Biotech',
 'FondChemical Co., Ltd',
 'CoreSyn',
 'Nextpeptide',
 'SharyPharm',
 'Lan Pharmatech',
 'VladaChem',
 'Apexmol',
 'Suntto Chemical',
 'AvaChem Scientific',
 'Tocris Bioscience',
 'Carbott PharmTech Inc.',
 'MolCore BioPharmatech',
 'Amadis Chemical',
 'Cooke Chemical Co., Ltd',
 'Glentham Life Sciences Ltd.',
 'Assembly Blocks Pvt. Ltd.',
 'Founder Pharma',
 'Kingston Chemistry',
 'Luminescence Technology Corp. (Lumtec)',
 'Synthonix, Inc.',
 'Wolves R&D chemical',
 'DC Chemicals',
 'Parkway Scientific',
 'AHH Chemical co.,ltd',
 'AbovChem LLC',
 'Biopharma PEG Scientific Inc',
 'DSL Chemicals',
 'Activate Scientific',
 'CD Formulation',
 'Jamson Pharmachem Technology',
 'Ennopharm',
 'NovoSeek',
 'SMID',
 'Eximed Laboratory',
 'LGC Standards',
 'Ark Pharm, Inc.',
 'Cangzhou Enke Pharma Tech Co.,Ltd.',
 'Smolecule',
 'LabNetwork, a WuXi AppTec Company',
 'Suzhou GeAo New Materials Co.,Ltd',
 'Hello Bio',
 'ChemTik',
 'MicroCombiChem GmbH',
 'AEchem Scientific Corp., USA',
 'Frinton Laboratories',
 'Fluoropharma Co.,Ltd',
 'AlchemyPharm',
 'InFarmatik',
 'MuseChem',
 'Shandong Youbang Biochemical',
 'AN PharmaTech',
 'Molepedia',
 'LEAPCHEM',
 'Watanabe Chemical Ind.',
 'Alfa Chemistry',
 'ChemBlock',
 'BroadPharm',
 'Bangyong Technology  Co., Ltd.',
 'BePharm Ltd.',
 'Syntree',
 'Hunan Huateng Pharmaceutical Co., Ltd.',
 'Irvine Chemistry Lab',
 'Porsechemical',
 'Fragmenta',
 'EDASA Scientific',
 'Up-Fluorochem',
 'ChemDiv',
 'CHESS fine organics',
 'labseeker',
 'OtavaChemicals',
 'Matrix Scientific',
 'Aaron Chemicals LLC',
 'AbMole Bioscience',
 'Rosewachem',
 'A2B Chem',
 'iChemical Technology USA Inc',
 'King Scientific',
 'SpiroChem',
 'Exclusive Chemistry Ltd',
 'zealing chemical',
 'ACT Chemical',
 'Fluoropharm Co.,Ltd',
 'BioCrick',
 'ABBLIS Chemicals',
 'Aribo Reagent',
 'Enamine',
 'Creasyn Finechem',
 'PepTech Corp.',
 'Milestone Pharmtech USA Inc.',
 'MolMall',
 'AZEPINE',
 'Wubei-Biochem',
 'AOBChem USA',
 'Achemtek',
 'Yuhao Chemical',
 'Achemica',
 'Pi Chemicals',
 'Annker Organics',
 'Shenzhen Nexconn Pharmatechs. Ltd',
 'RSChem, LLC',
 'OChem',
 'Biopurify Phytochemicals',
 'KCS Online',
 'Adooq BioScience',
 'J&W PharmLab',
 'TripleBond',
 'True PharmaChem',
 'EMD Millipore',
 'AAA Chemistry',
 'Hangzhou Trylead Chemical Technology',
 'Axispharm',
 'Toref Standards',
 'Boerchem',
 'Selleck Chemicals',
 'BerrChem',
 'Wuhan Atomole Chemicals Co., Ltd.',
 'BenchChem',
 'Angene Chemical',
 'Biosynth',
 'Achemo Scientific Limited',
 'Career Henan Chemical Co',
 'Inhibitor 2',
 'Beijing Advanced Technology Co, Ltd',
 'Aromsyn catalogue',
 'Georganics',
 'MP Biomedicals',
 'RR Scientific',
 'AKos Consulting & Solutions',
 'Boronpharm',
 'Vesino Industrial Co., Ltd',
 'BioChemPartner',
 'ChangChem',
 'Rieke Metals, LLC',
 'CSNpharm',
 'Ark Pharma Scientific Limited',
 'Oakwood Products',
 'Aromalake Chemical',
 'Infinium PharmaChem Pvt Ltd',
 'A&J Pharmtech CO., LTD.',
 'Shangyu Catsyn Co., Ltd.',
 'R&D Chemicals',
 'Debye Scientific Co., Ltd',
 'Key Organics/BIONET',
 'ABI Chem',
 'Total TOSLab Building-Blocks',
 'Aronis',
 'China MainChem Co., Ltd',
 'Innovapharm',
 'Alichem',
 'Wilshire Technologies',
 'W&J PharmaChem',
 'OXCHEM CORPORATION',
 'Changzhou Highassay Chemical Co., Ltd',
 'Life Chemicals',
 'Changzhou Naide Chemical',
 'Iodochem',
 'Hoffman Fine Chemicals',
 'Zjartschem',
 'Santa Cruz Biotechnology, Inc.',
 'ChemFish Tokyo Co., Ltd.',
 'Phion Ltd',
 'Chem-Space.com Database',
 'SYNCHEM OHG',
 'abcr GmbH',
 'MedChemexpress MCE',
 'Ampyridine Co.,Ltd',
 'Chemodex Ltd.',
 '1st Scientific',
 'Angel Pharmatech Ltd.',
 'Yick-Vic Chemicals & Pharmaceuticals (HK) Ltd.',
 'Synblock Inc',
 'ChemFaces',
 'ZINC',
 'Macsen Labs',
 'THE BioTek',
 'Xlinebio',
 'TimTec',
 'ForeChem',
 'Hefei Hirisun Pharmatech Co., Ltd',
 'Fisher Chemical',
 'Vichem Chemie Ltd.',
 'HDH Pharma',
 '4C Pharma Scientific Inc',
 'Day Biochem',
 'A1 BioChem Labs',
 'J&H Chemical Co.,ltd',
 'Van Aroma',
 'Wutech',
 'Hunan Chemfish Pharmaceutical Co., Ltd.',
 'AK Scientific, Inc. (AKSCI)',
 'Aceschem Inc',
 'Shanghai Sinofluoro Chemicals Co., Ltd',
 'SelectLab Chemicals GmbH',
 '3B Scientific (Wuhan) Corp',
 'Thoreauchem',
 'Acorn PharmaTech Product List',
 'Shanghai Mathcon Pharmaceutical Co.,LTD.',
 'Active Biopharma',
 'Sarchem Laboratories, Inc.',
 'Halochem',
 'BioAustralis Fine Chemicals',
 'Broadpharm',
 'AA BLOCKS',
 'Alinda Chemical',
 'Advanced Technology & Industrial Co., Ltd.',
 'A2Z Chemical',
 'AOB Chem USA',
 'ChemBridge',
 'Clearsynth',
 'Tyger Scientific',
 'Shanghai FWD Chemicals Limited',
 'Heterocyclics Research Chemicals & Building blocks',
 'Bestdo Inc',
 'Zhuorui Chemical Technology Co.,Ltd',
 'Boroncore',
 'Accela ChemBio Inc.',
 'Ambinter',
 'TargetMol',
 'Clinivex',
 'Race Chemical',
 'Nantong Baihua Bio-Pharmaceutical Co., Ltd',
 'Starshine Chemical',
 'ASINEX',
 'Heteroz LLC',
 'SQUARIX GmbH',
 'Qingdao Truelight Functional Materials Technology Co., Ltd.',
 'Qingdao Wanyuan Mountain Biotech Co.,Ltd',
 'LKT Labs',
 'Herbest Bio-Tech',
 'Nature Science Technologies Ltd',
 'Uchem Meditech',
 ]

EXAMPLE_VENDORS_TO_KEEP = ['TCI (Tokyo Chemical Industry)',
 'NIH Clinical Collection',
 'Ambeed',
 'Combi-Blocks',
 'Shanghai Tauto Biotech Co., Ltd',
 'Thermo Fisher Scientific',
 'ChemProbes',
 'Sigma-Aldrich',
 'VWR, Part of Avantor']

def canonicalize_smiles(smiles: str):
    '''

    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToSmiles(mol, canonical=True)
    else:
        print(f'SMILES {smiles} made None mol. Returning SMILES.')
        return smiles

def smiles_to_inchi(smiles: str) -> str:
    '''

    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToInchi(mol)
    else:
        print(f'SMILES {smiles} made None mol.')

def smiles_to_inchi_key(smiles: str) -> str:
    '''

    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToInchiKey(mol)
    else:
        print(f'SMILES {smiles} made None mol.')

def remove_duplicate_inchi_keys(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    '''
    Removes indices of df that posess a duplicate
    InChI Key in the INCHI_KEY column.

    Must have INCHI_KEY column

    Returns two dataframes: one that is the filtered and one
    that is the set that was removed
    '''

    # Get the duplicate INCHI_KEYS
    c = Counter(df['INCHI_KEY'])
    dups = []
    for k, v in c.items():
        if v > 1:
            print(f'Found duplicate of {k}')
            dups.append(k)

    # Get the entries in the df that have a duplicated inchi key
    filtered_out = df[df['INCHI_KEY'].isin(dups)].copy(deep=True)

    # Get a new dataframe that only has the unique InChI keys
    filtered = df[~df['INCHI_KEY'].isin(dups)].copy(deep=True)

    # Drop the duplicates from the filtered out
    #TODO check if this is only exact duplicates
    filtered_out = filtered_out.drop_duplicates(subset='INCHI_KEY')

    # Add the filtered out that has the dropped
    # duplicates back into the filtered dataframe
    filtered = pd.concat([filtered, filtered_out], axis=0)

    return filtered, filtered_out

def get_cid_from_inchi_key(inchi_key: str):
    '''
    Gets the CID from Pubchem using the PUG REST API.
    '''
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchi_key}/cids/JSON/'
    with urllib.request.urlopen(url) as u:
        data = u.read()
    data = data.decode('utf-8')
    j = json.loads(data)

    return j['IdentifierList']['CID'][0]

def get_vendor_list_from_cid(cid: int) -> list:
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

    chemical_vendors = chemical_vendors['Sources']

    return [vendor["SourceName"] for vendor in chemical_vendors]

def convert_str_list(s: str) -> list:
    if isinstance(s, list):
        return s
    return list(ast.literal_eval(s))

def _remove_from_str_list(list: list[str], to_remove_list: list[str]) -> str:
    '''
    Helper function for df.apply to remove vendors from
    the list_of_vendors if theyare in vendors_to_remove
    '''
    return [x for x in list if x not in to_remove_list]

def remove_specific_vendors_from_dataframe(dataframe: pd.DataFrame, vendors: list):
    dataframe['VENDORS'] = dataframe['VENDORS'].apply(_remove_from_str_list, to_remove_list=vendors)
    return dataframe

def draw_molecules_to_grid_image(smiles: list[str],
                                 mols_per_row: tuple = 5,
                                 maxMols: int = 50,
                                 img_resolution: int = 500,
                                 save_path: Path = None,
                                 ):
    '''
    Draws the molecules in dataframe['SMILES'] to a file
    in filtered_molecules
    '''

    if isinstance(save_path, str):
        save_path = Path(save_path)

    # Specify grid shape
    grid_width = mols_per_row

    # Specify resolution of the molecule images
    image_size = (img_resolution, img_resolution)

    mols = [Chem.MolFromSmiles(s) for s in smiles]

    png = Draw.MolsToGridImage(mols,
                            molsPerRow=grid_width,
                            subImgSize=image_size,
                            highlightAtomLists=None,
                            highlightBondLists=None,
                            maxMols=maxMols,
                            returnPNG=False,
                            legends=[f"SMILES: {s}" for s in smiles])

    if save_path is not None:
        png.save(save_path)

    return png