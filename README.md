# CommercialCompoundSearcher

This repository provides a set of functions and complementary Jupyter Notebook for assessing the commercial
availability of molecules. This notebook will "clean" a list of SMILES strings with RDKit and query PubChem for
[PubChem CIDs](https://pubchem.ncbi.nlm.nih.gov/docs/compounds) and chemical vendors. The remaining cells of the
notebook are for filtering vendors. Unfortunately, price queries are not supported at this time.

## Installation

A conda environment is provided in the `env.yml` file. To create the environment, run

```bash
conda env create --name commercialSearch --file=env.yml
```

Then activate the environment with

```bash
conda activate commercialSearch
```

## Additional Features

## Contributors

- James Howard, PhD
