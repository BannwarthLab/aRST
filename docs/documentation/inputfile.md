# Input File Documentation

This page describes the `input.toml` file of *aRST version 2.3.0*.
!!! danger "Incompatible input format"
    The 2.3.0 input format is **not compatible** with 2.2.0: most sections were renamed. For the old format, see the [2.2.0 input file documentation](../archive-2.2.0/documentation/inputfile.md).

**Contents**

[TOC]

---

## General structure

A minimal `input.toml` only needs the reactant(s):

```toml
[reactants.reactant1]
name = "coord1.xyz"
charge = 0
multiplicity = 1
```

Every other section is optional and falls back to its defaults. In TOML format:

- Square brackets `[section.subsection]` start a configuration block. The dot (`.`) means that `subsection` belongs to `section`.
- Each entry inside a block is a key-value pair, `key = value`.

The input file has the following blocks:

| Block | Content |
|---|---|
| [`[reactants.reactantX]`](#input-structures-reactantsreactantx) | input structures |
| [`[read]`](#pre-processing-read) | pre-processing of the input structures |
| [`[working_dir]`](#working-directories-working_dir) | run folder and output folders |
| [`[reactivity]`](#reactivity-assessment-reactivity) | how reactive atom pairs are ranked |
| [`[search]`](#reaction-search-search) | which reactions are explored |
| [`[path_modeling]`](#path-modeling-path_modeling) | how product guesses are generated |
| [`[refinement]`](#product-refinement-refinement) | optimisation of new structures |
| [`[energy]`](#energy-evaluation-energy) | energies, thresholds and transition states |
| [`[structures]`](#structure-comparison-structures) | how identical structures are recognised |
| [`[calculators]`](#calculator-settings-calculators) | shared settings for *xtb*, *ORCA* and *TeraChem* |
| [`[resources]`](#computational-resources-resources) | cores and memory |


!!! note "Atom indices"
    All atom indices in the input file are **0-based**: atom `0` is the first atom in the XYZ file.

---

## Input structures `[reactants.reactantX]`

`X` is the reactant number (`reactant1`, `reactant2`, …, up to `reactant99`). At least one reactant is required.

| Key | Default value, *Type* | Description |
|---|---|---|
| `name` | –, *str* | XYZ file of the reactant, read from `working_dir.work_dir`. **Required.** |
| `charge` | –, *int* | Charge. **Required.** |
| `multiplicity` | –, *int* | Spin multiplicity (> 0). **Required.** |
| `number` | `1`, *int* | How many copies of this reactant take part in the search, e.g. `number = 2` for a dimerisation of one molecule. |

---

## Pre-processing `[read]`

| Key | Default value, *Type* | Description |
|---|---|---|
| `opt` | `false`, *bool* | Optimise the input structures before the search, at the level set in [`[refinement]`](#product-refinement-refinement). By default the input geometries are used as they are. |

---

## Working directories `[working_dir]`

| Key | Default value, *Type* | Description |
|---|---|---|
| `work_dir` | current directory, *str* | Run folder. The XYZ files are read from here and all output is written here. |

---

## Reactivity assessment `[reactivity]`

Reactive atom pairs for association reactions are ranked by AFO gaps (and bonds by MBO), computed at the level selected here.

| Key | Default value, *Type* | Description |
|---|---|---|
| `backend` | `"xtb"`, *str* | Program that computes the AFO/MBO features:<br>**`"xtb"`** = xTB molden file, method from `[calculators.xtb]`;<br>**`"orca"`** = DFT/HF molden file from *ORCA*, converted with *Multiwfn* and read by the modified *sTDA*. Defaults to BHandHLYP/vDZP;<br>**`"tc"`** = DFT/HF molden file from *TeraChem*, read by the modified *sTDA*. Defaults to BHandHLYP/vDZP. |
| `overrides` | `{}`, *table* | Calculator settings that apply to this step only, see [Overrides](#overrides). |
| `equivalent_afogap_tol` | `0.01`, *float* | Tolerance for treating two AFO gaps as equal (degenerate pairs), in *eV*. |
| `spr_mode` | `0`, *int* | Which spin combinations of HOAO–LUAO are considered for an atom pair:<br>**0** = all ($\alpha$–$\alpha$, $\beta$–$\beta$, $\alpha$–$\beta$, $\beta$–$\alpha$);<br>**1** = same spin only ($\alpha$–$\alpha$, $\beta$–$\beta$);<br>**2** = opposite spin only ($\alpha$–$\beta$, $\beta$–$\alpha$).<br>No effect for closed-shell reactants. |
| `exclude_elements` | `[]`, *[int, …]* | Atomic numbers of elements that never react, e.g. `[8]` excludes all reactions involving Oxygen element. |
| `exclude_pairs` | `[]`, *[[int, int], …]* | Element pairs (atomic numbers) that never react with each other, e.g. `[[1, 8]]` excludes H–O pairs. |

!!! important
    For `backend = "orca"`, set `iloadasCart= 1` in the `settings.ini` file of *Multiwfn* (see [Installation](../installation.md)).

**Example**

```toml
[reactivity]
backend = "orca"

[reactivity.overrides]
functional = "PBE0"
basis = "def2-SVP"
```

---

## Reaction search `[search]`

### General settings `[search]`

| Key | Default value, *Type* | Description |
|---|---|---|
| `cycle` | `1`, *int* | Number of elementary reaction steps searched in a row. In every cycle after the first, the products of the previous cycle are used as new reactants. |
| `keep_every_cycle` | `3`, *int* | From cycle 2 on: how many reactions of the previous cycle are explored further, the ones with the lowest reaction energy first. |

A reaction type is only searched when its `count` is not `0` or `atom_pairs` is given. All three types are off by default.

### Intermolecular association `[search.inter]`

Needs at least two reactants.

| Key | Default value, *Type* | Description |
|---|---|---|
| `count` | `0`, *int* or `"all"` | Number of atom pairs to simulate, ranked by the lowest AFO gap. `"all"` simulates every pair. |
| `atom_pairs` | `[]`, *[[int, int], …]* | Atom pairs to simulate instead of the ranking. `[[0, 3]]` means atom 0 of `reactant1` and atom 3 of `reactant2` form a bond. |
| `concerted_pair` | `false`, *bool* | Simulate concerted reactions: two bonds form at the same time. With `atom_pairs`, all given pairs are simulated together in one path. |
| `aligning_distance` | `5.0`, *float* | Distance at which the two reactants are placed before the scan, as a multiple of the sum of the covalent radii of the reactive pair. |

### Intramolecular association `[search.intra]`

| Key | Default value, *Type* | Description |
|---|---|---|
| `count` | `0`, *int* or `"all"` | Number of atom pairs to simulate, ranked by the lowest AFO gap. `"all"` simulates every pair. |
| `atom_pairs` | `[]`, *[[int, int], …]* | Atom pairs of `reactant1` that should form a bond, e.g. `[[0, 5]]`. |

### Dissociation `[search.dissociation]`

| Key | Default value, *Type* | Description |
|---|---|---|
| `count` | `0`, *int* or `"all"` | Number of bonds to break, ranked by the lowest MBO. `"all"` breaks every bond. |
| `atom_pairs` | `[]`, *[[int, int], …]* | Bonds to break instead of the ranking, e.g. `[[8, 9]]`. In mechanochemical mode (`aRST -m`) these bonds are pulled **in addition** to the `count` top-ranked bonds. |

#### Pulling distance `[search.dissociation.pulling]`

Sets how far a bond is stretched. Choose one `mode`:

| `mode` | Further keys | Target distance |
|---|---|---|
| `"auto"` *(default)* | – | Open-chain bond: 5 × the sum of the covalent radii. Ring bond: the length of the alternative (the shortest) molecular graph path around the ring + 1 Å. |
| `"scale"` | `distance`, *float* (required) | `distance` × the sum of the covalent radii. |
| `"increment"` | `distance`, *float* (required) | Current bond length + `distance`, in *Å*. |

```toml
[search.dissociation.pulling]
mode = "scale"
distance = 1.5
```

#### Mechanochemical analysis `[search.dissociation.mechanical_settings]`

Only used by `aRST -m`. See [Mechanochemical Analysis](../examples/mechano.md) for the workflow.

| Key | Default value, *Type* | Description |
|---|---|---|
| `pull_position` | `"weakest_bond"`, *str* | What is pulled:<br>**`"weakest_bond"`** = the top-ranked bonds of the mechanophore B-set;<br>**`"terminals"`** = the terminal atom pairs (force application points). |
| `rank_by` | `"mbo"`, *str* | How B-set bonds are ranked, smallest first: **`"mbo"`** = Mayer bond order, **`"k"`** = FFiTS force constant, **`"dist"`** = bond length. |
| `use_d0` | `true`, *bool* | Use D0 (including the ZPE change, needs `energy.mode = "g"`) instead of De for the rupture force. |
| `terminal_pairs` | detected automatically, *[[int, int], …]* | Your own force application points. If not given, *aRST* finds the force-transmission paths itself. |

---

## Path modeling `[path_modeling]`

Settings for generating the product guess of every reaction path.

| Key | Default value, *Type* | Description |
|---|---|---|
| `scan_mode` | `"xtb"`, *str* | **`"xtb"`** = constrained xTB scan along the reactive atom pair;<br>**`"mbff"`** = *MolBar* harmonic force field, followed by a constrained xTB optimisation. |
| `backend` | `"xtb"`, *str* | Only `"xtb"` is available. |
| `overrides` | `{}`, *table* | Calculator settings for the scan, see [Overrides](#overrides). Scan-specific keys: `etemp` (default **4000** K), `force_constant` (default `1.0` Eh/bohr²), `steps` (default `20`), `opt_lvl` (default `"normal"`). |
| `enable_md` | `false`, *bool* | Run a short xTB MD simulation on the end structure of the scan before it is relaxed. The lowest-energy frame is used. |

!!! note
    The scan always uses its own electronic temperature (4000 K by default). An `etemp` in `[calculators.xtb]` does not change it; set it with `[path_modeling.overrides]` instead.

```toml
[path_modeling]
scan_mode = "xtb"

[path_modeling.overrides]
etemp = 4000
steps = 30
```

### Molecular dynamics `[path_modeling.md]`

Used when `enable_md = true`. The values are passed to the `$md` block of *xtb*.

| Key | Default value, *Type* | Description |
|---|---|---|
| `time` | `2.0`, *float* | Total simulation time, in *ps*. |
| `step` | `0.5`, *float* | Time step, in *fs*. |
| `temp` | `300.0`, *float* | Thermostat temperature, in *K*. |
| `shake` | `0`, *int* | SHAKE constraints: **0** = off, **1** = X–H bonds only, **2** = all bonds. |
| `dump` | `10`, *int* | Interval for writing trajectory frames, in *fs*. |
| `overrides` | `{}`, *table* | Calculator settings for the MD, see [Overrides](#overrides). |

---

## Product refinement `[refinement]`

New structures (and the input structures with `read.opt = true`) are optimised at this level.

| Key | Default value, *Type* | Description |
|---|---|---|
| `backend` | `"xtb"`, *str* | **`"xtb"`** = method from `[calculators.xtb]` (GFN2-xTB by default). If an optimisation fails, it is retried at 4000 K;<br>**`"orca"`** = *ORCA* optimisation, r2SCAN-3c by default. |
| `overrides` | `{}`, *table* | Calculator settings for this step, see [Overrides](#overrides). |

---

## Energy evaluation `[energy]`

| Key | Default value, *Type* | Description |
|---|---|---|
| `sp_method` | `"gxtb"`, *str* | Single-point method for reaction energies: **`"gxtb"`** = g-xTB, **`"dft"`** = DFT. |
| `sp_backend` | `"xtb"`, *str* | Program for the single points. Use `"xtb"` with `sp_method = "gxtb"` and `"orca"` with `sp_method = "dft"` (r2SCAN-3c by default). |
| `overrides` | `{}`, *table* | Calculator settings for the single points, see [Overrides](#overrides). |
| `mode` | `"e"`, *str* | Energy used to judge reactions:<br>**`"e"`** = electronic energy;<br>**`"g"`** = free energy: single-point energy + thermal correction (xTB frequencies through *ORCA*) + solvation correction (xTB, only if a solvent is set). |
| `elim` | `300.0`, *float* | Highest accepted reaction energy, in *kcal/mol*. Reactions above it are recorded but not explored further. |
| `tseval` | `false`, *bool* | Search a transition state for every accepted reaction. |
| `ts_method` | `"nebts"`, *str* | Transition state method: **`"nebts"`** = NEB-TS, **`"ffits"`** = FFiTS. |
| `tslim` | `50.0`, *float* | Highest accepted barrier, in *kcal/mol*. Only used with `tseval = true`. |

```toml
[energy]
sp_method = "dft"
sp_backend = "orca"
mode = "g"

[energy.overrides]
functional = "wB97X-3c"
```

---

## Structure comparison `[structures]`

*aRST* uses the *MolBar* identifier to decide whether two structures are the same.

| Key | Default value, *Type* | Description |
|---|---|---|
| `mbcut_length` | `7`, *int* | Number of leading *MolBar* fields that are compared. Smaller values compare less of the identifier, so more structures count as identical. |

---

## Calculator settings `[calculators]`

These blocks hold the settings that are **shared by every step** that uses the program. Only keys you write here are passed on; everything else keeps the default of the step (see the table in [Overrides](#overrides)).

### *xtb* `[calculators.xtb]`

| Key | Default value, *Type* | Description |
|---|---|---|
| `method` | `"gfn2"`, *str* | xTB method: `"gfn0"`, `"gfn1"`, `"gfn2"`, `"gfnff"` (run with `XTB_BIN`) or `"gxtb"` (run with `GXTB_BIN`). |
| `solvent` | `"ch2cl2"`, *str* | Implicit solvent. Any solvent name accepted by *xtb*. |
| `solvent_model` | `"alpb"`, *str* | Solvation model. g-xTB always uses `"cosmo"`. |
| `solvent_reference` | `"bar1M"`, *str* | Reference state for the solvation free energy (only used with `energy.mode = "g"`). |
| `etemp` | `300`, *int* | Electronic temperature, in *K*. |

The optimisation level (`opt_lvl`, default `"normal"`) is not a shared setting. Set it per step in `[refinement.overrides]` or `[path_modeling.overrides]`.

### *ORCA* `[calculators.orca]`

| Key | Default value, *Type* | Description |
|---|---|---|
| `functional` | `"BHandHLYP"` (reactivity),<br>`"r2scan-3c"` (refinement, energy), *str* | Functional, written to the *ORCA* input as given. |
| `basis` | `"vdzp"` (reactivity),<br>`""` (refinement, energy), *str* | Basis set. Leave empty for composite methods such as r2SCAN-3c. |
| `solvent` | `"ch2cl2"`, *str* | Implicit solvent. |
| `solvent_model` | `"cpcm"`, *str* | Solvation model. |

### *TeraChem* `[calculators.tc]`

| Key | Default value, *Type* | Description |
|---|---|---|
| `functional` | `"BHandHLYP"`, *str* | Functional. |
| `basis` | `"vdzp"`, *str* | Basis set. |
| `epsilon` | `9.08`, *float* | Dielectric constant for the PCM solvation model. |

### Overrides

Each step uses one calculator. The settings of a step are built in this order, later ones win:

1. the defaults of the step (table below),
2. the keys you set in `[calculators.<program>]`,
3. the keys in the step's `overrides` table.

| Step | Block with `backend` / `overrides` | Programs | Step defaults |
|---|---|---|---|
| reactivity | `[reactivity]` | xtb, orca, tc | ORCA/TeraChem: BHandHLYP/vDZP |
| scan | `[path_modeling]` | xtb | `etemp = 4000`, `force_constant = 1.0`, `steps = 20` |
| MD | `[path_modeling.md]` | xtb | – |
| refinement | `[refinement]` | xtb, orca | ORCA: r2SCAN-3c |
| energy | `[energy]` (`sp_backend`) | xtb, orca | xtb: g-xTB with COSMO; ORCA: r2SCAN-3c |

The `backend` cannot be changed inside `overrides`.

**Example:** GFN2-xTB in water everywhere, but GFN1-xTB for the reactivity:

```toml
[calculators.xtb]
method = "gfn2"
solvent = "h2o"

[reactivity.overrides]
method = "gfn1"
```

---

## Computational resources `[resources]`

| Key | Default value, *Type* | Description |
|---|---|---|
| `cpu_budget` | SLURM or all cores, *int* | Total number of cores. Under SLURM: `--ntasks` × `--cpus-per-task` (a missing one counts as 1). Otherwise all cores available to the process. |
| `xtb_threads` | `cpu_budget`, *int* | Threads per *xtb* run. Cannot exceed `cpu_budget`. |
| `orca_nprocs` | `cpu_budget`, *int* | Processes per *ORCA* run. Cannot exceed `cpu_budget`. |
| `orca_maxcore_mb` | `1000`, *int* | Memory per *ORCA* process, in *MB*. |

---

## Full example

```toml
[reactants.reactant1]
name = "coord1.xyz"
charge = 0
multiplicity = 1

[reactants.reactant2]
name = "coord2.xyz"
charge = 0
multiplicity = 1

[read]
opt = true

[search]
cycle = 2
keep_every_cycle = 3

[search.inter]
count = 5

[search.dissociation]
count = 3

[path_modeling]
scan_mode = "xtb"

[path_modeling.overrides]
etemp = 4000
steps = 30

[energy]
sp_method = "gxtb"
sp_backend = "xtb"
mode = "e"
elim = 50.0

[calculators.xtb]
method = "gfn2"
solvent = "h2o"

[resources]
cpu_budget = 8
```