# CommercialCompoundSearcher
A set of functions and complementary Jupyter Notebook for assessing the commercial availability of molecules. This notebook will "clean" a list of SMILES strings with RDKit and query Pubchem for [Pubchem CIDs](https://pubchem.ncbi.nlm.nih.gov/docs/compounds) and a list of vendors. The remaining cells of the notebook are for filtering vendors.

## Installation
A conda environment is provided in the `env.yml` file. To create the environment, run

```
conda env create --name commercialSearch --file=env.yml
```

Then activate the environment with

```
conda activate commercialSearch
```

## Requirements
If you'd like to use your own environment, here is a list of known dependencies.

1.  [Python 3.11.6](https://www.python.org/downloads/release/python-3116/)
2.  [RDKit 2022.9.5](https://pypi.org/project/rdkit/)
2.  [Pandas 2.1.2](https://pandas.pydata.org/docs/getting_started/install.html)

## Contributors

-   James Howard, PhD
