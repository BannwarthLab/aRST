# Installation

This page describes how to install *aRST* 2.3.0 and the external programs it calls.

**Contents**

[TOC]

---

## Requirements

### Python
*aRST* 2.3.0 supports Python 3.9-3.12

!!! warning "Python 3.13 is not supported yet"
    *aRST* needs [molbar](https://gitlab.git.nrw/rwth-bannwarthlab/molbar) 1.1.x, that has no build for Python 3.13, yet version 1.2 changed its API.

### External programs

| Program | Need for |
|---|---|
| [xtb](https://github.com/grimme-lab/xtb) >=6.7.1 | GFN*n*-xTB optimisations and scans |
| [g-xtb](https://github.com/grimme-lab/g-xtb) >=2.0.0, or [ORCA](https://www.faccts.de/orca/) >=6.0.0 | single-point calculations |
| [Modified sTDA](https://gitlab.git.nrw/rwth-bannwarthlab/arst) | AFO features from molden files |
| [FFiTS](https://gitlab.git.nrw/rwth-bannwarthlab/ffits) *(optional)* | force constants *k* for the mechanochemical analysis (`-m`) and force-field generation of TS guess |
| [Multiwfn](http://sobereva.com/multiwfn/download.html) 3.8.noGUI *(optional)* | AFO features with the ORCA backend |
| TeraChem *(optional)* | AFO features with the TeraChem backend |
| [CREST](https://github.com/crest-lab/crest) *(optional)* | protonation / deprotonation tasks |

---

## Step 1: Pull the source code and install

Clone the repository from Git:
```bash
git clone https://gitlab.git.nrw/rwth-bannwarthlab/arst.git
python3.12 -m venv arstenv          # or python3.11
source arstenv/bin/activate
pip install --upgrade pip
cd aRST
pip install -e .                    # installs aRST and the `aRST` command
```
---

## Step 2: Build modified sTDA
*aRST* requires *sTDA* for reading Cartesian `MOLDEN` type file.
A precompiled version of modified *sTDA* is provided in the repository root `stda-v1.0-intel-linux-x86_64.tar.gz`.
Unpack the archive in your working directory:
```
tar -xzf stda-v1.0-intel-linux-x86_64.tar.gz
```
This will extract the modified *sTDA* executable, which is required for AFO feature generation.

## Step 3: Set the binary paths

The template `set_environment_variables.sh` is in the top-level `aRST/` folder.
Replace every `path/to/...` placeholder in it with the real path on your machine:

```bash
# required
export XTB_BIN="/path/to/xtb"              
export GXTB_BIN="/path/to/gxtb"
export STDA_BIN="/path/to/stda" 

# If you need run DFT-level calculations using ORCA:
export ORCA_BIN="/path/to/orca/orca"      
export MWFN_BIN="/path/to/Multiwfn"      # reactivity.backend = "orca" 
#
# If you need run DFT-level calculations using TeraChem:
export TC_BIN="/path/to/terachem"        # reactivity.backend = "tc" 

# If you need to generate TS guess using FFiTS
export FFITS_VENV="/path/to/ffits/ffitsenv"   # the separate FFiTS venv
export FFITS_BIN="/path/to/ffits"        # alternative if you do not use a venv
# FFiTS also calls xtb: keep xtb on PATH as well
export PATH="$(dirname "$XTB_BIN"):$PATH"
```

---

### Configure Multiwfn settings
!!! Important "This setting only requied when you are generating AFO feature at DFT level."
Since *sTDA* requires Cartesian basis functions but *ORCA* supports only spherical basis functions as output, the use of *Multiwfn* is necessary.
To enable the conversion fuction in *Multiwfn*, you need to set `iloadasCart= 1` in its `settings.ini` setting file.


### FFiTS in its own virtual environment
!!! Important "This setting only requied when you are generating TS guess using Force-field method or when you are running `-m` command in aRST"
!!! danger "Do not install FFiTS into `arstenv`"
    FFiTS pins its own versions of numpy and other Python packages, and these conflict with the
    dependencies of *aRST* (`molbar` pins `networkx==3.1`, for example).
    *aRST* never imports FFiTS: it only runs the `ffits` command as a subprocess,
    so FFiTS can live in a **separate virtual environment**.
```bash
deactivate 2>/dev/null              # leave the aRST venv first
git clone git@gitlab.git.nrw:rwth-bannwarthlab/ffits.git
cd ffits
python3.12 -m venv ffitsenv
source ffitsenv/bin/activate
make install-dev
deactivate
```


## Step 4: Check the installation

```bash
source set_environment_variables.sh
source arstenv/bin/activate
aRST -v                             # prints the banner and the configured binary paths
```

!!! tip "Running on a cluster"
    `source` the environment script and activate `arstenv` inside the job script.
    *aRST* takes its core count from `#SBATCH --ntasks` × `#SBATCH --cpus-per-task`
    (a missing one counts as 1), or from `resources.cpu_budget` in the input file.

**Enjoy working with aRST!**
