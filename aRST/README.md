# aRST: an Automated Reaction Search Tool

This package provides as an automated tool for modeling chemical reaction paths, searching thermodynamic-allowed
products based on the given electronic structures. Other useful built-in functions include evaluating atomic reactive
sites and generating the reactive complex.

## Getting Started

for now (test version):
source my version of QC progromm and venv environment

please talk to me before you use this package

```
source "/mnt/group-bin/venv-python/ychen/aRST/bin/activate"
module load openmpi/4.1.1
module load orca
module load xtb
module load libraries/libgfortran/5
module load terachem/dev-3c

cmd for calling aRST
```

now you are able to used aRST as commandline tool or call its functions :)
To start the calculation, the essential input consists of the .xyz file of the structure, the .toml file with detailed
settings.

## User Manual for setting input file

### essential input

```
[geom.reactant1]
name = "coord1.xyz"                 # str, file name of coord1
charge = 0                          # int,chagrge of coord1
multiplicity = 1                    # int,multiplicity of coord1

[geom.reactant2]
name =
charge =
multiplicity =

...
```

For inputs with a stoichiometric number more than 1, you can add under each entry with

```
[geom.reactantX]
number =                            # stoichiometric number of coord
```

### HOAO and LUAO Infomation

Print out reactive information of given structures.

with essential input, calling command with

```
aRST in.toml --ml_feature > record.log 2>&1 &
```

### Alignment

Aligning two strucs with given atom index.

besides essential input, adding target atomic loation for each coord with

```
[geom.reactantX]
at =                                # reactive atom number                      
```

calling command with

```
aRST in.toml --align > record.log 2>&1 &
```

### Reaction Pathes Exploration

Explore reaction pathes based on given reactants.

besides essential input, (optional) adding detailed searching settings with

```
[read.general]
dftopt = False                      # bool, doing dftopt when first reading input coord
gfn2opt = False                     # bool, doing gfn2opt when first reading input coord

[searching_setting.general]
cycle = 1                           # int, number of searching cycle
elim =  300                         # float, energy limitation to kill this reaction path, unit kcal/mol
hydrogen = True                     # bool, searching range includs/excludes H atoms
halogen = True                      # bool, searching range includs/excludes halogen atoms
scan_mode = 0                       # 0(default): use xtb-scan simulating reactions
                                    # 1: use molbar FF simulating reactions
sep_mode = 0                        # 0(default): directly seperate fragments from scan-komplex  
                                    # 1: do unconstrain opt to scan-komplex and seperate fragments                                      
mb_level =                          # None or "loose" to only screen with topology spectrum


[searching_setting.reactivity]
ml_method ="gfn2"/"stda"            # str, method for ml feature calculations
##############################################
MO_method = "orca"/"terachem"       # str, program to get molden file
### if MO_method=terachem, you will need: ###
ml_tc_functional="BHandHLYP"            # str, if use stda, method for struc sp calculations in TC
ml_tc_basis="vdzp"                  # str, if use stda, basis for struc sp calculations in TC
ml_tc_epsilon=                      # float, if use stda, pcm solvent mode
##############################################
### if MO_method=orca, you will need: ###
ml_orca_functional="BHandHLYP"            # str, if use stda, method for struc sp calculations in ORCA
ml_orca_basis="vdzp"                  # str, if use stda, basis for struc sp calculations in ORCA
ml_orca_cpcm=                         # str, if use stda, cpcm solvent mode in ORCA
##############################################

[searching_setting.inter]
number = 1                          # int, number of inermolecular reaction paths to be search in one cycle
search_mode = 1                     # 0: only search association possibility between same structures (A+A -> and B+B ->);
                                    # 1(default): only search association possibility between different structures (A+B ->);
                                    # 2: search association possibility in all structures (A+A ->, A+B-> and B+B ->).


loose_bond =                        # bool (defalut as False), loose a bit the weakest bond in the struc
maxovlap =                          # bool (defalut as False), different alignment direction
aligning_distance = 5               # int, scale of aligning distance

[searching_setting.disso]
number = 1                          # int, number of inermolecular reaction paths to be search in one cycle
mode="bo/wboe"                      # str, searching mode              
pulling_distance =                  # int, scale of pulling distance

[searching_setting.intra]
number =                            # int, number of inermolecular reaction paths to be search in one cycle
loose_bond =                        # bool (defalut as False), loose a bit the weakest bond in the struc


[xtb_setting]
strucopt=False                     # bool, do xtb preopt for struc
method = "gfn2"                    # str, method for xtb calculations
optlevel = ""                      # str, xtb acceptable optimization levels
alpb = "h2o"                       # str, alpb solvent mode
struc_etemp = 1000                 # int, electronic temp for struc opt
scan_etemp = 1000                  # int, electronic temp for scan
k = 1                              # geom constrain strength (used in scan)
step = 20                          # pushing/pulling steps
xtb_wd = "/home/guests/ychen/public/xtbml/_build/xtb" 

[orca_setting]
functional = "r2scan-3c"                # str, method for dft calculations
basis = None
dftopt = False                      # bool, do dft opt
solmode = "water"                   # str, CPCM solvent mode
ncpu = 1                            # int, parallel calculations for orca

[stda_setting]
stda_wd = "/home/guests/ychen/public/stdaml/stda"
```

if products are known, you can add following struc infomation, aRST will compare each simulated reaction path results
for products

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

commandline

```
aRST in.toml --explore > record.log 2>&1 &
```
