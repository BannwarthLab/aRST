# Installation

This page contains detailed information about the installation of aRST.

## Prepare
- [Python](https://www.python.org/) 3.8
- [xtb](https://github.com/grimme-lab/xtb) 6.7.0 or later
- [ORCA](https://www.faccts.de/orca/) 6.0.0 or later or [g-xtb](https://github.com/grimme-lab/g-xtb)
- [Multiwfn](http://sobereva.com/multiwfn/download.html) 3.8.noGUI
- Modified sTDA (see below) 

## Step 1: Pull the source code
Clone the repository from Git:
```angular2html
git clone https://git.rwth-aachen.de/bannwarthlab/aRST.git
cd aRST
```

## Step 2 (Recommended): Create a virtual environment
```angular2html
python -m venv venv
source venv/bin/activate
```

## Step 3: Install the package
```angular2html
pip install -e .
```
You can check that the package is installed: `pip list | grep aRST`.

## Step 4: Build modified sTDA
A precompiled version of modified *sTDA* is provided in the repository root `stda-v1.0-intel-linux-x86_64.tar.gz`.
Unpack the archive in your working directory:
```
tar -xzf stda-v1.0-intel-linux-x86_64.tar.gz
```
This will extract the modified *sTDA* executable, which is required for AFO feature generation.

## Step 5: Set up environment variables
In the repository root, locate the file: `set_environment_variables.sh`.
Edit this file and set the paths to the required external programs:
```angular2html
export XTB_BIN="/path/to/xtb"              # Semiempirical QM methods (GFNn-xTB)
export STDA_BIN="/path/to/modified/stda"   # Modified sTDA executable

# Single-point calculations
export ORCA_BIN="/path/to/orca"            # ORCA
export GXTB_BIN="/path/to/gxtb"            # g-xTB
```
For parallel usage set the threads for OMP and the MKL linear algebra backend by
```angular2html
export OMP_NUM_THREADS=<ncores> MKL_NUM_THREADS=<ncores>
```
For larger systems please adjust the stack size accordingly, otherwise stack overflows will occur. 
Use something along the lines of this:
```angular2html
ulimit -s unlimited
export OMP_STACKSIZE=4G
```
After editing the file, activate the environment variables by
```angular2html
source set_environment_variables.sh
```

## Step 6: Activate the aRST virtual environment
After setting the environment variables, activate the aRST Python environment:
```angular2html
source venv/bin/activate
```

** Enjoy working with aRST! **