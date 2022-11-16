# AUA4GRO

Automated routine for the creation of GROMACS input files for MD simulation of organic molecules.

## Instalation

### Pre-requisite

This instruction of installation is based on a computer with Ubuntu OS. It is also requeried haviong Anaconda installed.

### Create anaconda enviroment

Open the terminal on the folder project with the `conda init` activated and enter the following command:

```
conda env create -f aua4gro.yml
```

Activate the new enviroment:

```
conda activate aua4gro
```

## Usage

Run the Python script `AppAUA4GRO.py`:

```
python AppAUA4GRO.py
```

Introduce the SMILES codes and press enter. The script will generate coordinate file ('.gro') and topology file ('.top').

## Citation

Please cite this work as: Castro-Anaya, L., Gomez, S.Y., Orozco, G.. Comprehensive automated routine implementation, validation, and benchmark of the anisotropic force field (AUA4) using Python and GROMACS. XXXX. 

