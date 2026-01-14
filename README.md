# Welcome to the documentation site of aRST

## Introduction 
aRST is a Python package that automates chemical reaction exploration from a quantum chemical perspective.
By treating elementary reactions as bond formation and bond breaking events, aRST is able to identifies potential reactive sites, carries out bond changing simulations, and generates chemically reasonable product structures.

Chemical bonding usually occurs between nucleophilic and electrophilic sites.
To capture these bonding tendencies, aRST assigns each atomic site an average orbital energy value for both occupied (HOAO) and virtual (LUAO) orbitals, which indicate nucleophilicity and electrophilicity, respectively.
We set a weighting factor for each molecular orbital and atom by considering the contribution of that MO to each atom according to its Mulliken population and a denominator that amplifies the importance of orbitals near the Fermi level: 

$$w_{A,i} = \frac{p_{A,i}}{\left ( \varepsilon_{i} - \mu \right )^2} \ \ \& \ \ w_{A,a} = \frac{v_{A,a}}{\left ( \mu - \varepsilon_{a} \right )^2}$$
Here, $p_{A,q} = n_i\sum_{\kappa \in A} \sum_\lambda C_{\kappa i}C_{\lambda i}S_{ \kappa \lambda }$ is the Mulliken electron and  \linebreak[4]
$v_{A,a} = (1-n_a)\sum_{\kappa \in A} \sum_\lambda C_{\kappa a}C_{\lambda a}S_{ \kappa \lambda }$ is the Mulliken vacancy population.
$n_q$ is the orbital occupation number, 
$\varepsilon_q$ is orbital energy of orbital $q$, $\mu$ is the Fermi-level which we approximate as the arithmetic mean of the HOMO and LUMO.
It is obtained as the weighted value of occupied orbitals contributing to atom $A$:
$$\text{HOAO}\left(A\right) = \mu - \frac{\sum_i w_{A,i} \left(\mu - \varepsilon_i \right)}{Q_A\sum_j w_{A,j}}$$
$$\text{LUAO}\left(A\right) = \mu + \frac{\sum_a w_{A,a} \left( \varepsilon_a - \mu \right)}{V_A\sum_b w_{A,b}}$$
where $V_A=\sum_a p_{A,a}$ becomes the net Mulliken vacancy population.
Therefore, reactive sites can be identified by a numerical comparison of the HOAO-LUAO energy difference (AFO gaps).

In the current version of aRST, the bond breaking event locatiosn are predicted using Mayer Bond Order. 

After deteriming the atomic reacvtive sites, aRST surpports the use of either the Newton trajectory scan at the GFN2-xTB semi-empirical method or harmonic force field to force the breaking or formation of target bonds.
For bimolecular reactions, the two reactant structures are pre-aligned to form a reactive complex prior to the bond-changing process. 
During this bond-changing process, only the identified atom pair(s) are strictly constrained during these bond-changing events, while the remaining degrees of freedom in the system are allowed to relax. This approach enables other bonds in the system to form or break naturally.

To process the product structure, aRST first separates possible fragments from the primary product. 
Charges are then assigned to each fragment based on the atomic charges in the original structure. 
Finally, the geometry of each fragment is optimized at the GFN2-xTB level to obtain the final product structures.

- ChemRxiv Paper: [An automated intermolecular reaction discovery approach relying on heuristic atom-partitioned frontier orbital features](https://chemrxiv.org/engage/chemrxiv/article-details/6872432de1957b8c618377a6)
- Source code: [aRST_gitlab](https://git.rwth-aachen.de/bannwarthlab/aRST)
- Bug reports: [aRST_issuse](https://git.rwth-aachen.de/bannwarthlab/aRST/-/issues)
- Email contact：[ying.chen@pc.rwth-aachen](mailto:ying.chen@pc.rwth-aachen)

**Contents**

- [Installation](#installation)
- [Preparing Input Files](#preparing-input-files)
- [Running a Task](#running-a-task)
  - [Reaction Exploration](#reaction-exploration)
  - [Generating AFO Features](#generating-afo-features)
  - [Aligning reactive structures](#aligning-reactive-structures)


## Installation

1. Prepare a venv envirment:
```
python3 -m venv path/to/venv
source path/to/venv/bin/activate
```
2. Download source codes
   + a) from repository:
      ```
      git clone git@git.rwth-aachen.de:bannwarthlab/aRST.git
      ```
   + b) from pip:
      ```
      pip install arst-chem
      ```
3. Set up the environment variables necessary for executing QC programs: 

change env paths in `./set_environment_variables.sh` 
and then `source set_environment_variables.sh`.

## Preparing input files
1. XYZ file of reactants

Note: Only one reactant structure can exist in each xyz file. If you are exploring bimolecular reactions, you should prepare two xyz files.
<details>
  <summary>coord1.xyz</summary>
  <pre><code>
4

C            0.00000038603655       -0.00000010798198        0.00044709980789
H            1.03077378701329       -0.27758440553714       -0.00613760517443
H           -0.75578436932805       -0.75390194272838        0.00335566910181
H           -0.27498980372179        1.03148645624750        0.00233483626474
  </code></pre>
</details>
<details>
  <summary>coord2.xyz</summary>
  <pre><code>
1

Cl           0.00000000000000        0.00000000000000        0.00000000000000
  </code></pre>
</details>

2. `.toml` setting file
The syntax of TOML mainly consists of `key = value` pairs, allowing you to customize the details of the reaction exploration settings. 
The following section introduces the exploration parameters supported by aRST.

### Available parameter in aRST
+ **Processing of reactant structures**

```
[geom.reactant1]
name="coord1.xyz"                      # str, file name of coord1
charge=0                               # int, chagrge of coord1
multiplicity=2                         # int, multiplicity of coord1

[geom.reactant2]
name="coord2.xyz"                      # str, file name of coord2
charge=0                               # int, chagrge of coord2
multiplicity=2                         # int, multiplicity of coord2

[read.general] 
dftopt=False                           # bool, doing dftopt when first reading input coord
gfn2opt=False                          # bool, doing gfn2opt when first reading input coord
jumprefine=False                       # bool, jump all sp/dft calculations and read structure information from folder "struc"
```

+ **Generating GFN2-xTB-based AFO features in xtb type molden file**
```
[searching_setting.reactivity]
afo_method="gfn2"                      # str, define generating GFN2-xTB-based AFO features
```
+ **Generating DFT-based AFO features from ORCA type molden file**
```
[searching_setting.reactivity]
afo_method="stda"                      # str, define generating DFT or HF-based AFO features
MO_method="orca"                       # str, running sp calculation in ORCA to generate molden file
afo_orca_functional="BHandHLYP"        # str, functional method for sp calculations
afo_orca_basis="vdzp"                  # str, basis method for sp calculations
afo_orca_cpcm="ch2cl2"                 # str, keyword with the name of the solvent used in cpcm solvent model
```
+ **Generating DFT-based AFO features from TeraChem type molden file**
```
[searching_setting.reactivity]
afo_method="stda"                      # str, define generating DFT or HF-based AFO features
MO_method="terachem"                   # str, running sp calculation in TC to generate molden file
afo_tc_functional="BHandHLYP"          # str, if use stda, method for structure sp calculations in TC
afo_tc_basis="vdzp"                    # str, if use stda, basis for structure sp calculations in TC
afo_tc_epsilon=9.08                    # float, the dielectric constant of solvents used in pcm solvent model
```

+ **General settings for exploring reactions**
```
[searching_setting.general]
cycle=1                                # int, number of searching cycle
elim=300                               # float, energy threshold for $E_{\mathrm{el, products} - E_{\mathrm{el, reactants}}$ to keep this trial reaction path, unit in kcal/mol
hydrogen=True                          # bool, if include the H atoms when generating AFO reactive pair preditions
halogen=True                           # bool, if include the halogen atoms when generating AFO reactive pair preditions

spr_mode=0                             # int(0/1/2), there are up to eight different AFO feature combinations for each pair of open-shell reactants (different HOAO or LUAO providers and among different spin channels).
                                       # 0(default): consider only the lowest AFO gap from the eight possible combinations as a repre-sentative measure of this atom-atom pair
                                       # 1: consider only the lowest AFO gap from the same spin channels (HL_aa and HL_bb)
                                       # 2: consider only the lowest AFO gap from the opposite spin channels (HL_ab and HL_ba)
scan_mode=0                            # int(0/1), method for simulating bond-change events
                                       # 0: use scanning algorithm to perform a relaxed PES scan in xtb
                                       # 1: use harmonic force field (mb FF) optimization to define expected product topology using MolBar package, tuning parameter for this FF see below.
sep_mode=0                             # int(0/1), treatment style of primary product structure after processing bond-change events
                                       # 0(default): separate the passible fragments directly from the superstructure (For example, in SN2 reactions, multiple product structures may already exist within the superstructure. Separating the geometry before optimization helps push the reaction to the product site.)
                                       # 1:  do fully relaxed optimization on the superstructure and then separate the passible fragments
                                                       
mb_level=""                            # None(default) or "loose": Molbar identifier level for identifying newly generately product structures
                                       # None: use full Molbar identifier to distinguish newly generated structures
                                       # "loose": only use topology spectrum part of Molbar identifier to distinguish newly generated structures
```
+ **Search settings for intermolecular association reactions based on top reactive AFO gap values**


```
[searching_setting.inter]
number=1                           # int, number of top reactive intermolecular atom pairs to be searched in one cycle
aligning_distance=5                # int, separating distance between reactant structures based on the sum of covalent radii of reactive atom pair
maxovlap=False                     # bool, alignment style for reactant structures
                                   # False(default): Minimal geometric overlap between the two reactants while maintaining exposed reaction sites
                                   # True: Maxmial geometric overlap between the two reactants while maintaining exposed reaction sites
concerted_pair=False               # bool, if search for the second equivalent atom pair with degenerated AFO gap (useful for Diels–Alder reactions)
loose_bond=False                   # bool, if loose a bit the weakest bond in the preliminary product structure by pulling this bond to 1.2 times of its current length 
```
 
+ **Search settings for intramolecular association reactions based on top reactive AFO gap values**

```
[searching_setting.intra]
number=1                           # int, number of top reactive intramolecular atom pairs to be searched in one cycle
concerted_pair=False               # bool, if search for the second equivalent atom pair with degenerated AFO gap (useful for Diels–Alder reactions)
```

+ **Search settings for dissociation reactions based on the weakest mayer bond order**

```
[searching_setting.disso]
number=1                           # int, number of the weakest bonds to be searched in one cycle
concerted_pair=False               # bool, if search for the second equivalent atom pair with degenerated AFO gap (useful for Diels–Alder reactions)
pulling_distance=2                 # int, scale of pulling distance
depotonation=False                 # bool, if search for potential stable depotonated form of reactant structures in CREST
potonation=False                   # bool, if search for potential stable potonated form of reactant structures in CREST
```

+ **Manually define reaction sites**
```
[searching_setting.inter/intra/disso]
at1=0                              # int, the index of reactive atom in coord1.xyz
at2=0                              # int, the index of reactive atom in coord2.xyz

# for manually define a second bond formation:
at3=1                              # int, the index of reactive atom in coord1.xyz
at4=1                              # int, the index of reactive atom in coord2.xyz
```

+ **Perform molecular dynamics (MD) simulations on product structures**
```
[searching_setting.general]
md=True                             # bool, if performing short MD simulations to test the stability of the product structure

[searching_setting.md]
time=10                             # int, total run time of MD simulation, unit in ps
step=1                              # int, time step for propagation, unit in fs
temp=300                            # int, electronic temp, unit in K
shake=1                             # 0/1/2, 0: not using shake algorithm; 1: constrain H bonds; 2: constrain all bonds
dump=50                             # int, interval for trajectory printout, unit in fs
```

**Setting for computational methods**
```
[searching_setting.singlepoint]
method="gxtb"/"dft"                 # str, running either gxtb or DFT methods in ORCA for computing $\E_{\mathrm{el}}$

[xtb_setting]
method="gfn2"                       # str, method for GFN calculations
optlevel=""                         # str, xtb acceptable optimization levels
alpb="h2o"                          # str, alpb solvent mode
struc_etemp=300                     # int, electronic temp used in geometry optimization, unit in K
scan_etemp=1000                     # int, electronic temp used in bond-formation process, unit in K
k = 1                               # geom constrain strength used in bond-formation process, unit in Eh/Bohr2
step = 20                           # pushing/pulling steps

[orca_setting]
functional="r2scan-3c"              # str, functional method used in geometry optimization and sp calculations
basis=None                          # str, basis used in geometry optimization and sp calculations
dftopt=False                        # bool, if doing optimization for product geometries
solmode="ch2cl2"                    # str, keyword with the name of the solvent used in cpcm solvent model
ncpu = 1                            # int, running a calculation in parallel

[mbFF_setting]                      
repulsion_charge=100                # float/int, charge used in the Coulomb term in a harmonic potential E_coul = sum Q/r_ij  
```

## Running a Task

### Reaction Exploration
To generate product structures based on the given reactant structures,
execute the following command in the folder where the input files are located:
```
aRST in.toml --explore > record.log 2>&1 &
```

After the calculation completes, the following files and folders will be generated:

`record.log`: a log file recording the execution process of the task

`struc.json`: a JSON file containing structural details, including ID, XYZ geometry, charge, spin multiplicity (m), energy level, and MolBar code

`reaction.json`: a JSON file summarizing the reactions identified during the search

`scan/`: a folder containing detailed records of all trial reaction paths

`struc/`: a folder containing computational details for each structure used in the search

`react_l`: a JSON file that ranks the reactivity of all atomic pair combinations based on their AFO gap values


### Generating AFO Features
To obtain a preview of the predicted reactions based on the given reactant structures,
execute the following command in the folder where the input files are located:
```
aRST in.toml --afo_feature > record.log 2>&1 &
```
After the calculation completes, in addition to the `record.log`, and `struc/` files/folders, the following files and folders will also be generated:

`react_l_interA`: a JSON file that ranks the bimolecular association reactivity of all atomic pair combinations based on their AFO gap values

`react_l_intraA`: a JSON file that ranks the unimolecular association reactivity of all atomic pair combinations based on their AFO gap values

`react_l_disso`: a JSON file that ranks the dissociation reactivity of all atomic pair combinations based on their MBO values


### Aligning reactive structures
To obtain a pre-aligned reactive supramolecular structure based on the two given reactant structures,
execute the following command in the folder where the input files are located:
```
aRST in.toml --align > record.log 2>&1 &
```
After the calculation completes, in addition to the `record.log`, and `struc/` files/folders, the following files and folders will also be generated:

`alignment`: a folder containing supramolecular structure XYZ files aligned by both minimal and maximal overlap styles


