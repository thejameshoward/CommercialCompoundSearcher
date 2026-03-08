#!/usr/bin/env python3
# coding: utf-8

'''
For reading and writing files.
'''

from __future__ import annotations

import pandas as pd
from pathlib import Path

def read_table_file(file: Path,
                    sheet_name: str | None) -> pd.DataFrame:
    '''
    Read a delimited table file into a pandas DataFrame. Supports .txt, .csv,
    and .xlsx files. The DataFrame must contain the column name "SMILES". The
    sheet_name must be specified if using a .xlsx file.

    Parameters
    ----------
    file: Path
        Path to the input file. Supported extensions are `.txt` and `.csv`.

    sheet_name: str | None
        The name of the sheet in a .xlsx file.

    Returns
    -------
    df: pd.DataFrame
        DataFrame containing the parsed table data.

    Raises
    ------
    ValueError
        If `file` does not have a supported extension.
    '''

    # Read in the file
    if file.suffix == '.txt':
        df = pd.read_table(file, header=0)
    elif file.suffix == '.csv':
        df = pd.read_csv(file, header=0)
    elif file.suffix == '.xlsx':
        if sheet_name is None:
            raise ValueError(f'Must specify sheet_name to read in .xlsx files.')
        df = pd.read_excel(file, sheet_name=sheet_name, engine='openpyxl')
    else:
        raise ValueError(f'{file.name} does not have a supported extension.')

    # Check that the file is formatted correctly
    if 'SMILES' not in df.columns:
        raise KeyError('The column "SMILES" must be in your provided spreadsheet or table.')

    return df
