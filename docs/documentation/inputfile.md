# Input File Documentation

This page contains a guide to *aRST version 2.2.0* input files <INPUTFILE> <br>
Some keys may not be supported in older versions of *aRST*. For configuration details of earlier releases, please refer to the archived `README.md`.

**Contents**

[TOC]

---

A general `input.toml` file is structured as follows:

```toml

[geom.reactant1]
name="coord1.xyz"
charge=0                 
multiplicity=2

[read.general] 
gfn2opt=true

```

In TOML format:

- Square brackets `section.subsection` define a configuration block. The dot (`.`) denots that `subsection` belongs
  to `section`.
- Each entry inside a block is written as a key-value pair in the form `key=value`.

!!! note
    The input file should be placed in the **same folder** as the XYZ input geometry file(s).


---

## Input geometry settings `[geom.reactantX]`

`X` represents the reactant index (e.g., `reactant1`, `reactant2`).

| Key            | Default value, *Type*    | Description                                  |
|----------------|--------------------------|----------------------------------------------|
| `name`         | `"basename1.xyz"`, *str* | Specify the XYZ file of input geometry X     |
| `charge`       | `-`, *int*               | Specify the charge of input geometry X       |
| `multiplicity` | `-`, *int*               | Specify the multiplicity of input geometry X |

---

## Pre-processing of input geometry `[read.general]`

| Key          | Default value, Type | Description                                                                                                                                                                                                                                                               |
|--------------|---------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `dftopt`     | `false`, *bool*     | Specify whether to perform geometry optimization at the DFT level in *ORCA*. Detailed *ORCA* settings are defined in the `[orca_setting]` block                                                                                                                           |
| `gfn2opt`    | `false`, *bool*     | Specify whether to perform geometry optimization at the GFN2-xTB level. Detailed *xtb* settings are defined in the `[xtb_setting]` block                                                                                                                                  |
| `jumprefine` | `false`, *bool*     | Specify whether to skip all single-point/optimization pre-processing of inputfiles and read information from folder `struc`.<br/> If set to *true*, the reactant index must match the corresponding index in `struc` (e.g., read `reactant1` information from `struc/0`). |

!!! note
    *aRST*, by default, will not reoptimize the input geometries. 

---

## Reactivity settings `[searching_setting.reactivity]`

`afo_method` defines how AFO features (and MBO features) are generated.<br>
Allowed values: `"gfn2"`,`"stda"`,`"gxtb"` (`afo_method=gxtb` can currently only be used internally within the Bannwarth group)

### Option 1: GFN2-xTB AFO `afo_method="gfn2"`

Use this option when you want to generate GFN2-xTB-based AFO and MBO features.<br>
This enables these reactivity features to be generated from a Molden file produced via single-point calculations in *xtb*. Detailed *xtb* settings are specified in the `[xtb_setting]` block.

### Option 2: DFT/HF AFO `afo_method="stda"`

Use this option when you want to generate DFT/HF-based AFO and MBO features.<br>
This enables these reactivity features to be generated from a Molden file produced via *ORCA* or *TeraChem* single-point calculations.
Specifically, *aRST* invokes *ORCA* or *TeraChem*, as specified by the user to generate a Molden file. It then calls *Multiwfn* to convert spherical basis functions into Cartesian form. The converted Molden file is then used by the
modified *sTDA* procedure to generate AFO features.

**Example**

```toml
[searching_setting.reactivity]
afo_method="stda"
MO_method="orca" # or "terachem"
```

[//]: # (todo)
!!! important
    modified *sTDA* can be found in xx.<br>
    To enable the *Multiwfn*'s conversion fuction, you must set `iloadasCart= 1` in the `settings.ini` setting file of *Multiwfn*.

#### Option 2a: Generate molden file via ORCA `MO_method="orca"`

The level of theory for AFO features can be selected using the following keys:

| Key                   | Default value, Type  | Description                                                                                                                               |
|-----------------------|----------------------|-------------------------------------------------------------------------------------------------------------------------------------------|
| `afo_orca_functional` | `"BHandHLYP"`, *str* | Specify the DFT method for AFO feature generation. The *str* value will be passed to the *ORCA* input file.                               |
| `afo_orca_basis`      | `"vdzp"`, *str*      | Specify the basis set for AFO feature generation. The *str* value will be passed to the *ORCA* input file.                                |
| `afo_orca_cpcm`       | `"ch2cl2"`, *str*    | Specify the solvent used in the CPCM solvation model for AFO feature generation. The *str* value will be passed to the *ORCA* input file. |

#### Option 2b: Generate molden file via TeraChem `MO_method="terachem"`

The level of theory for AFO features can be selected using the following keys:

| Key                 | Default value, Type  | Description                                                                                                                                                    |
|---------------------|----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `afo_tc_functional` | `"BHandHLYP"`, *str* | Specify the DFT method for AFO feature generation. The *str* value will be passed to the *TeraChem* input file.                                                |
| `afo_tc_basis`      | `"vdzp"`, *str*      | Specify the basis set for AFO feature generation. The *str* value will be passed to the *TeraChem* input file.                                                 |
| `afo_tc_epsilon`    | `9.08`, *float*      | Specify the dielectric constant of solvents used in the PCM solvation model for AFO feature generation. The value will be passed to the *TeraChem* input file. |

---

## Reaction exploration

### General settings `[searching_setting.general]`

This block contains information on *how* elementary reaction-step searching is defined in *aRST*, including the number of reaction steps, number, the specification of the reactive atoms, the workflow for generating product guess structures, thermodynamic energy thresholds, and whether isomeric product structures are considered when accepting a new reaction.

| Key                     | Default value, Type                                 | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
|-------------------------|-----------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `wd`                    | `"."`, *str*                                        | Specify the working directory where the input files are located.                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `cycle`                 | `1`, *int*                                          | Specify the number of elementary steps to be automatically searched via *aRST*.                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `sp_method`             | `dft`, *str*                                        | Specify the theoretical method used for single-point energy calculations.<br/>**"dft"**=use DFT method in *ORCA* (detailed *ORCA* settings are defined in the `[orca_setting]` block);<br/>**"gxtb"**=use the g-xTB method.                                                                                                                                                                                                                                                                                                        |
| `elim`                  | `300`, *int*                                        | Specify the allowed electronic energy difference threshold for accepting a new reaction step, unit in *kcal/mol*.                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `exclude_pairs`         | `[[Z1,Z2],[Z3,Z4],...]`, <em>[[int, int], ...]</em> | Specify the element pair that should be excluded from the reaction search, e.g., `[[1,8]]` to exclude reactions between *H* and *O* atoms.                                                                                                                                                                                                                                                                                                                                                                                         |
| `exclude_elements`      | `[Z1,Z2,...]`, *[int,int,...]*                      | Specify the elements that should be excluded from the reaction search, e.g.,`["O"]` to exclude all reactions involes *O* atom.                                                                                                                                                                                                                                                                                                                                                                                                     |
| `equivalent_afogap_tol` | `0.01`, *float*                                     | Specify the tolerance for identifying the denegerated AFO gaps, unit in *eV*.                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `spr_mode`              | `0`, *int*                                          | When evaluating the reactivity of a given atom pair (e.g., A–B) based on AFO gaps, multiple HOAO–LUAO combinations may exist due to $\alpha$ and $\beta$ electrons This option specifies which spin combinations are considered:<br>**0**=all combinations ($\alpha$–$\alpha$, $\beta$–$\beta$, $\alpha$–$\beta$, $\beta$–$\alpha$);<br/>**1**=same-spin only ($\alpha$–$\alpha$, $\beta$–$\beta$);<br/>**2**=opposite-spin only ($\alpha$–$\beta$, $\beta$–$\alpha$).<br/>For closed-shell reactants, this setting has no effect. |
| `scan_mode`             | `0`, *int*                                          | Specify which workflow is used to generate products:<br>**0**=Scanning algorithm: performs an xTB PES scan along the interatomic distance of the reactive atom pair (detailed *xtb* settings are defined in the `[xtb_setting]` block);<br/>**1**=Harmonic force field algorithm: generates products via harmonic FF optimization using the *MolBar* package                                                                                                                                                                       |
| `sep_mode`              | `0`, *int*                                          | After the scan, multiple product fragments may exist within the precursor structure. This option specifies the order of separation and optimization of the supramolecule:<br>**0**=First separate the possible fragments, then fully relax each fragment geometry;<br/>**1**=First fully relax the supramolecule, then separate the possible fragments.                                                                                                                                                                            |
| `md`                    | `false`, *bool*                                     | Specify wether to perform a short xTB-MD simulation on supramolecule (after scan). If `sep_mode=1` and `md=true`, the supramolecule is first fully relaxed, then subjected to MD, and finally separated into possible fragments. Detailed MD settings are defined in the `[searching_setting.md]` block.                                                                                                                                                                                                                           |
| `mb_level`              | `0`, *int*                                          | *aRST* use MolBar identifier to distinguish structural uniqueness. This option specifies how a new structure is identified:<br/>**0**=full MolBar identifier; <br/>  **1**=topology spectrum only (encoding atomic connectivity).                                                                                                                                                                                                                                                                                                  |

### Molecular Dynamics Simulations `[searching_setting.md]`

This block contains the settings for performing MD simulations with xTB.

| Key     | Default value, Type | Description                                                     |
|---------|---------------------|-----------------------------------------------------------------|
| `time`  | `10`, *int*         | Specify the total run time of MD simulation, unit in *ps*       |
| `step`  | `1`, *int*          | Specify time step for propagation, unit in *fs*.                |
| `temp`  | `300`, *int*        | Specify electronic temp, unit in *K*.                           |
| `shake` | `1`, *int*          | Specify shake mode: **0**=off; **1**=X-H only, **2**=all bonds. |
| `dump`  | `50`, *int*         | Specify the interval for trajectory printout, unit in *fs*.     |

### Search settings for intermolecular reactions `[searching_setting.inter]`

This block specifies the selection of reactive atom pairs for intermolecular reactions and workflow used to aligning reactants.

| Key                 | Default value, Type                       | Description                                                                                                                                                                                                                                                                                                                                                                                               |
|---------------------|-------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `number`            | `-`, *int*                                | Specify how many atom pairs, ranked by the lowest AFO gaps, should be simulated for intermolecular association events in one cycle.<br/>If a number bigger than the available number of atom pairs is specified, all available atom pairs in the system will be simulated.                                                                                                                                |
| `concerted_pair`    | `false`, *bool*                           | Specify whether the searched reactions are concerted. When set to `true`, only concerted intermolecular association reactions (i.e., reactive pairs include two equivalent atom pairs) will be simulated in one cycle.                                                                                                                                                                                    |
| `atom_pair`         | `[[AT1,AT2]]`, <em>[[int, int], ...]</em> | Specify the atomic indices (starting from 0) of each reactant structure to be simulated for intermolecular association events. For example,`[[0,0]]` means that atom 0 of `coord1.xyz` and atom 0 of `coord2.xyz` will form a bond. <br/>When this option is used, `number` is overwritten by the number of provided atom pairs. If more than one pair is given, `concerted_pair` is set to `true`.       |
| `aligning_distance` | `5`, *int*                                | Specify the alignment distance between the two reactant structures. The given value places the reactants along the z-axis, with their reactive centers positioned at (0,0,0) and (0,0,z). The separtaed distance is set as the sum of the covalent radii of the reactive atom pair.<br/>For concerted atom pairs, the reactive center is set as the midpoint of the two reactive atoms in each structure. |
| `maxovlap`          | `false`, *bool*                           | Specify the alignment style of given reactant structures. By default, *aRST* aligns the reactants so that the principal axes of the two structures face opposite directions, minimizing steric hindrance at the reaction site. For reactions such as Diels–Alder, setting this option to `true` may help generate different stereoisomeric products.                                                      |

### Search settings for intramolecular reactions `[searching_setting.intra]`

This block specifies the selection of reactive atom pairs for intramolecular association.

| Key              | Default value, Type                        | Description                                                                                                                                                                                                                                                                                                                                                               |
|------------------|--------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `number`         | `-`, *int*                                 | Specify how many atom pairs, ranked by the lowest AFO gaps, should be simulated for intramolecular association events in one cycle.<br/>If a number bigger than the available number of atom pairs is specified, all available atom pairs in the system will be simulated.                                                                                                |
| `concerted_pair` | `false`, *bool*                            | Specify the atomic indices (starting from 0) of reactant structure to be simulated for intramolecular association events. For example,`[[0,1]]` means that atom 0 and 1 of `coord1.xyz` will form a bond. <br/>When this option is used, `number` is overwritten by the number of provided atom pairs. If more than one pair is given, `concerted_pair` is set to `true`. |
| `atom_pair`      | `[[AT1, AT2]]`, <em>[[int, int], ...]</em> | Specify the atomic indices (starting from 0) of the atom pair to be simulated for intramolecular association events. When this option is used,  `number` is overwritten by the number of provided atom pairs. If more than one pair is given, `concerted_pair` is set to `true`.                                                                                          |

### Search settings for dissociations `[searching_setting.disso]`

This block specifies the selection of reactive atom pairs for dissociations.

| Key                                                             | Default value, Type                       | Description                                                                                                                                                                                                                                                                                                                                                                |
|-----------------------------------------------------------------|-------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `number`                                                        | `-`, *int* or *str*                       | Specify how many bonded atom pairs, ranked by the lowest Mayer bond order (MBO > 0.3), should be simulated for dissociations in one cycle.<br/>If the string `"all"` is given, all available atom pairs in the system will be simulated.                                                                                                                                   |
| `concerted_pair`                                                | `false`, *bool*                           | Specify whether the searched reactions are concerted. When set to `true`, only concerted dissociations (i.e., reactive pairs include two equivalent atom pairs) will be simulated in one cycle.                                                                                                                                                                            |
| `atom_pair`                                                     | `[[AT1,AT2]]`, <em>[[int, int], ...]</em> | Specify the atomic indices (starting from 0) of the reactant structure to be simulated for dissociation events. For example,`[[0,1]]` means that the bond between atom 0 and 1 of `coord1.xyz` will break. <br/>When this option is used, `number` is overwritten by the number of provided atom pairs. If more than one pair is given, `concerted_pair` is set to `true`. |
| `pulling_distance`                                              | `5`, *int* or *float*,or `"auto"`         | Specify the scaling factor for the pulling distance. The separated distance is set as a multiple of the sum of the covalent radii of the reactive atom pair.                                                                                                                                                                                                               |
| `pulling_distance_add`                                          | `-`, *int* or *float*                     | Specify the separation distance as the current interatomic distance of the reactive atom pair plus the given value, unit in *Å*.                                                                                                                                                                                                                                           |
| `depotonation`                                                  | `false`, *bool*                           | Specify whether to perform deprotonation using CREST before breaking the selected weakest bond. Note that deprotonation by CREST may reorder the atom indices.                                                                                                                                                                                                             |
| `potonation`                                                    | `false`, *bool*                           | Specify whether to perform potonation using CREST before breaking the selected weakest bond. Note that potonation by CREST may reorder the atom indices.                                                                                                                                                                                                                   |
| **The following keywords are supported only in `--fmax` mode:** |
| `terminal_pair`                                                 | `[[AT1,AT2]]`, <em>[[int, int], ...]</em> | Specify the atomic indices (starting from 0) of the force application points in the reactant structure. <br/>If not specified, aRST will automatically predict potential force-transmission pathways.                                                                                                                                                                      |
| `stretch_scale`                                                 | `1.5`, *int* or *float*                   | Specify the stretching factor for the bond to be ruptured. <br/>The factor should be chosen such that only the target bond is broken, while the remaining bonding topology is preserved.                                                                                                                                                                                   |
| `use_d0`                                                        | `false`, *bool*                           | Specify whether consider ZPE contribuition in the dissociation energy                                                                                                                                                                                                                                                                                                      |

---

## Computational methodology settings

### *xtb* settings `[xtb_setting]`

This block provides detailed settings for electronic structure calculations at the xTB level using *xtb*.

| Key           | Default value, *Type* | Description                                                                                                                                                                                                 |
|---------------|-----------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `method`      | `"gfn2"`, *str*       | Specify the xTB-based method. All valid string keywords accepted by *xtb* are supported.                                                                                                                    |
| `optlevel`    | `"normal"`, *str*     | Specify optimizing level. All valid string keywords accepted by *xtb* are supported.                                                                                                                        |
| `alpb`        | `"ch2cl2"`, *str*     | Specify the solvent used in the ALPB solvation model. All valid string keywords accepted by *xtb* are supported.                                                                                            |
| `struc_etemp` | `300`, *int*          | Specify the electronic temperature used in equilibrium electronic structure calculations. Unit in *K*.                                                                                                      |
| `scan_etemp`  | `1000`, *int*         | Specify the electronic temperature used in constrained geometry optimizations, such as when generating preliminary product structures using either scanning or harmonic force field algorithm. Unit in *K*. |
| `k`           | `1`, *int*            | Specify force constant used in the scanning algorithm. Unit in *Eh/Bohr<sup>2</sup>*.                                                                                                                       |
| `step`        | `20`, *int*           | Specify step number used in the scanning algorithm.                                                                                                                                                         |

### *ORCA* settings `[orca_setting]`

This block provides detailed settings for electronic structure calculations at the DFT or HF level using *ORCA*.

| Key          | Default value, *Type* | Description                                                                                                                                                                                                          |
|--------------|-----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `functional` | `"r2scan-3c"`, *str*  | Specify the functional method. All valid string keywords accepted by *ORCA* are supported.                                                                                                                           |
| `basis`      | `""`, *str*           | Specify the basis set. If a composite functional is used, this option can be left empty; otherwise, all valid string keywords accepted by *ORCA* are supported.                                                      |
| `dftopt`     | `"false"`, *bool*     | Specify to perform geometry optimization for the newly generated product structure using the specified DFT method.                                                                                                   |
| `cpcm`       | `"ch2cl2"`, *str*     | !!!CHECK CODE!!!Specify the solvent used in the CPCM solvation model. All valid string keywords accepted by *ORCA* are supported.                                                                                    |
| `ncpu`       | `1`, *int*            | Specify the number of processors used to run a calculation in *ORCA*. For further setup instructions, refer to the [ORCA parallel tutorial](https://www.faccts.de/docs/orca/6.0/tutorials/first_steps/parallel.html) |

---