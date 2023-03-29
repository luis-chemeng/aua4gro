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

## Description of files

- aua4gro.yml : Yaml file for the creation of conda enviroment.
- AppAUA4GRO.py : Automated routine (main python file).
- AUA4_parameters_2.py : Python file containing all AUA4 force field parameters.
- utilities.py : Python module containing classes and auxiliary functions used by the main code.
- opt/ : Folder containing files optimized for alkanolamines.

## Citation

Please cite this work as: Castro-Anaya, L., Gomez, S.Y., Orozco, G.. Comprehensive automated routine implementation, validation, and benchmark of the anisotropic force field (AUA4) using Python and GROMACS. J. Phys. Chem. A 2023, 127, 1555−1563.

Refer to the manuscript for more information: https://doi.org/10.1021/acs.jpca.2c08335

## License 

This project is licensed under the terms of the GNU General Public License v3.0 (GPU GPLv3).
