#!/bin/bash
# In this file, you can set the environment variables used in aRST
# next, in your current shell session, you need to source the script: source set_environment_variables.sh

################# necessary setting #################
export XTB_BIN="path/to/xtb" # for performing semiempirical quantum mechanical methods GFNn-xTB in xtb
export STDA_BIN="path/to/modified/stda" # to generate AFO features from molden files

# running single-point calculation in orca:
export ORCA_BIN="path/to/ORCA"
# running single-point calculation in gxtb:
export GXTB_BIN="path/to/gxtb"

################# optional setting #################
# generate molden files for AFO features in ORCA
# export ORCA_BIN="path/to/ORCA"
export MWFN_BIN="path/to/Multiwfn" # for converting spherical basis functions to Cartesian ones. Note: change settings.ini/iloadasCart= 1
# generate molden files for AFO features in TeraChem
export TC_BIN="path/to/TeraChem"
# perform potonation/depotonation tasks
export CREST_BIN="path/to/CREST"
