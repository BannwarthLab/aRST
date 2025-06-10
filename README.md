# aRST: an Automated Reaction Search Tool

This package provides an automated tool for modeling chemical reaction paths. 
It searchs thermodynamic-allowed products based on the given electronic structures. 

Other useful built-in functions include evaluating possible reactive site of given structure;
and aligning two structures at specified atomic sites to minimize the spatial resistance in between.

## Getting Started

Prepare a venv envirment:
```
python3 -m venv path/to/venv
source path/to/venv/bin/activate
```

Download source codes to target folder and unzip it.

Use pip install all required packages:
```
pip install -r /path/to/requirements.txt
```

Use pip install aRST package:
```
pip install -e /path/to/aRST
```

Before running a calculation, make sure the following quantum chemical programs have been loaded into your system. 
The software required specifically for calculations based on DFT-based reactivity features is labeled as "DFT-features required" in parentheses.
```
ORCA (version 6.0.0 for DFT-features required)
xtb
Multiwfn (DFT-features required)
sTDA (DFT-features required)
```

Now you are able to used aRST as commandline tool :)


## User Manual
To start a calculation, the essential input consists of the ```.xyz``` files of reactants, 
a ```in.toafo``` file with detailed settings.

### essential input for all type calculation in in.toafo

```
[geom.reactant1]
name = "coord1.xyz"                 # str, file name of coord1
charge = 0                          # int, chagrge of coord1
multiplicity = 1                    # int, multiplicity of coord1

[geom.reactant2]
name =
charge =
multiplicity =

...
```

For inputs with a stoichiometric number more than 1, you can add under each entry with

```
[geom.reactantX]
number =                            # int, stoichiometric number of coordX
```

### reactivity assessment: HOAO and LUAO
Print out reactive information of given structures.

With essential input, calling command:

```
aRST in.toafo --afo_feature > record.log 2>&1 &
```


### alignment

Aligning two structures with the given atom index.

Besides essential input, adding target atom index for each coord with

```
[geom.reactantX]
at =                                # int, specified atom index                      
```

Calling commandline:

```
aRST in.toafo --align > record.log 2>&1 &
```

### reaction pathes exploration

Explore reaction pathes based on given reactants.

Besides essential input, (optional) adding detailed searching settings with

```
[read.general]
dftopt = False                      # bool, doing dftopt when first reading input coord
gfn2opt = False                     # bool, doing gfn2opt when first reading input coord
jumprefine = False                  # bool, if jump all sp/dft calculations


[searching_setting.general]
cycle = 1                           # int, number of searching cycle
elim =  300                         # float, energy limitation to kill this reaction path, unit kcal/mol
hydrogen = True                     # bool, searching range includes/excludes H atoms
halogen = True                      # bool, searching range includes/excludes halogen atoms
spr_mode = 0                        # 0(default): globally compare both same spin and oppo spin HL combinations
                                    # 1: only compare same spin (HL_aa and HL_bb)
                                    # 2: only compare opposite spin (HL_ab and HL_ba)
scan_mode = 0                       # 0(default): use xtb-scan simulating reactions
                                    # 1: use molbar FF simulating reactions
sep_mode = 0                        # 0(default): directly separate fragments from scan-complex  
                                    # 1: do unconstrain opt to scan-complex and seperate fragments                   
                                                       
mb_level =                          # None or "loose" to only screen with topology spectrum
md = True                           # bool, if doing short md to test the stability of product structure

[searching_setting.md]
### if md=True, you will need: ###
time=10                             # int, total run time of MD simulation
step=1                              # int, time step for propagation
temp=300                            # int, electronic temp, by defalut using scan temp
shake=1                             # 0/1/2, 0: not using shake algorithm; 1: constrain H bonds; 2: constrain all bonds
dump=50                             # int, interval for trajectory printout

[searching_setting.reactivity]
afo_method ="gfn2"/"stda"            # str, method for afo feature calculations
##############################################
MO_method = "orca"/"terachem"       # str, program to get molden file
### if MO_method=terachem, you will need: ###
afo_tc_functional="BHandHLYP"        # str, if use stda, method for structure sp calculations in TC
afo_tc_basis="vdzp"                  # str, if use stda, basis for structure sp calculations in TC
afo_tc_epsilon=                      # float, if use stda, pcm solvent mode
##############################################
### if MO_method=orca, you will need: ###
afo_orca_functional="BHandHLYP"        # str, if use stda, method for structure sp calculations in ORCA
afo_orca_basis="vdzp"                  # str, if use stda, basis for structure sp calculations in ORCA
afo_orca_cpcm=                         # str, if use stda, cpcm solvent mode in ORCA
##############################################

[searching_setting.inter]
number = 1                          # int, number of intermolecular reaction paths to be search in one cycle
search_mode = 1                     # 0: only search association possibility between same structures (A+A -> and B+B ->);
                                    # 1(default): only search association possibility between different structures (A+B ->);
                                    # 2: search association possibility in all structures (A+A ->, A+B-> and B+B ->).

at1=                                # int (index starts with 0) manually define reactive atom
at2=

loose_bond =                        # bool(default as False), loose a bit the weakest bond in the structure
maxovlap =                          # bool(default as False), different alignment direction
aligning_distance = 5               # int, scale of aligning distance
concerted_pair= False               # bool, create concerted pair in simulations (good for DA reactions)

[searching_setting.inter]
number = 1
at1=                                # int (index starts with 0) manually define reactive atom
at2=
concerted_pair=False

[searching_setting.disso]
number = 1
at1=                                # int (index starts with 0) manually define reactive atom
at2=
mode = 0                            # 0: use mbo as bo
                                      1: use embo as bo
                                      
concerted_pair=False

pulling_distance = 2               # int, scale of pulling distance
depotonation = False               # boll, if detecting depotonated form of reactant structures


[xtb_setting]
strucopt=False                     # bool, do xtb preopt for structure
method = "gfn2"                    # str, method for xtb calculations
optlevel = ""                      # str, xtb acceptable optimization levels
alpb = "h2o"                       # str, alpb solvent mode
struc_etemp = 1000                 # int, electronic temp for structure opt
scan_etemp = 1000                  # int, electronic temp for scan
k = 1                              # geom constrain strength (used in scan)
step = 20                          # pushing/pulling steps
xtb_wd = "path/to/modified/xtb/version" 

[orca_setting]
functional = "r2scan-3c"            # str, method for dft calculations
basis = None
dftopt = False                      # bool, do dft opt for products
solmode = "water"                   # str, CPCM solvent mode
ncpu = 1                            # int, parallel calculations for orca

[stda_setting]
stda_wd = "path/to/modified/stda/version"

[mwfn_setting]
mwfn_wd = "path/to/mwfn/"    # change settings.ini/iloadasCart= 1
```

If products are known, you can add following structure infomation, aRST will compare each simulated reaction path results for products.

```
[geom.product1]
name =
charge =
multiplicity =

[geom.product2]
name =
charge =
multiplicity =
```

Calling commandline:

```
aRST in.toml --explore > record.log 2>&1 &
```
