# Mechanochemical Analysis
Search for the optimal force application point, rupture force and the correspongding rupture product using the given reactant structure with *aRST version dev*.

**Contents**

[TOC]

---
## Standard analyse
In the following, standard runs are demonstrated using the Cyclobutane (No.**2** of [Klein dataset](https://pubs.acs.org/jacsat/article/142/38/16364/618132/Validation-of-the-CoGEF-Method-as-a-Predictive)) as an example.

=== "command"

    ```
    aRST input.toml --fmax ff1.csv
    ```
    This command must be executed from the command line. `ff1.csv` is the force field fitted force constant file coming from FFiTS
    You can save the terminal output by appending `> record.log 2>&1` to the end of the command.
    
    The resulting log file will look similar to the example shown in the `output` tab above. 
    In addition, the file `struc.json` records detailed electronic information for all structures encountered during the search, while file `reaction.json` records detailed reaction information.
=== "input.toml"

    ```toml
    [geom.reactant1]
    name="coord.xyz"
    charge = 0
    multiplicity = 1
    
    [read.general]
    gfn2opt = true
    jumprefine=true
    [searching_setting.general]
    sp_method = "gxtb"
    scan_mode = 1
    exclude_elements=[1]
    mb_level=1
    [searching_setting.reactivity]
    afo_method = "gfn2"
    
    [searching_setting.disso]
    number = 1
    stretch_scale=1.5
    pulling_distance='auto'
    
    [xtb_setting]
    method = "gfn2"
    alpb = "ch2cl2"
    scan_etemp = 1000
    k = 1
    ```
=== "coord.xyz"

    ```
    42
    
    C          0.75489       -2.56549        0.36897
    C          0.45793       -3.99841       -0.15675
    C         -0.95848       -3.52446       -0.56080
    C         -0.44605       -2.05344       -0.49823
    C         -1.35572       -1.04760        0.16123
    C          2.09959       -1.97479        0.02641
    O          2.62356       -2.04947       -1.06561
    O          2.65159       -1.33248        1.08597
    C          3.93279       -0.70920        0.86519
    C          3.73261        0.68968        0.29855
    O          5.04521        1.27347        0.21306
    C          5.09433        2.53358       -0.28959
    O          4.11208        3.14927       -0.63567
    C          6.51887        3.03155       -0.34736
    O         -1.78737       -1.15191        1.29048
    O         -1.64348       -0.01277       -0.66562
    C         -2.50839        1.01335       -0.13724
    C         -3.96480        0.61657       -0.33542
    O         -4.74347        1.73414        0.12953
    C         -6.08994        1.59075        0.03317
    O         -6.62149        0.60607       -0.42682
    C         -6.79845        2.80709        0.58017
    H          0.53145       -2.45330        1.43117
    H          1.08095       -4.22098       -1.02741
    H          0.53687       -4.81398        0.56724
    H         -1.36123       -3.86171       -1.51968
    H         -1.69329       -3.70799        0.22823
    H         -0.10416       -1.69923       -1.47254
    H          4.40915       -0.66622        1.84668
    H          4.53129       -1.31992        0.18596
    H          3.27557        0.65135       -0.69360
    H          3.09662        1.29450        0.95243
    H          7.12041        2.37718       -0.98690
    H          6.53151        4.04816       -0.74147
    H          6.96647        3.01361        0.65180
    H         -2.29183        1.17256        0.92125
    H         -2.26976        1.91453       -0.70556
    H         -4.18320        0.42088       -1.38958
    H         -4.21309       -0.27640        0.24402
    H         -6.42813        3.71555        0.09496
    H         -7.87172        2.70512        0.41686
    H         -6.59667        2.90339        1.65254
    ```
=== "output"

    ```

     ==========================================================================================
     |                                                                                        |
     |                                                                                        |
     |                      automated Reaction Search Tool (aRST) V2.2.0                      |
     |                                                                                        |
     |                                                                                        |
     ==========================================================================================
     # Author:             Y. Chen, C. Bannwarth
     # Institute:          Institute of Physical Chemistry
     #                     RWTH Aachen University
     # Version date:       03.2026
     # Job running at:     10-08-2026 11:03:06
     ==========================================================================================

    Checking input TOML...done.
    
    Checking envirment:
    XTB_BIN checked, PATH=/bins/group-bin/xtb/6.7.1/install/bin/xtb
    stda_bin checked, PATH=/mnt/group-bin/aRST/dependencies/stda
    GXTB_BIN checked, PATH=/bins/group-bin/g-xtb/1.1.0/g-xtb-1.1.0/bin/gxtb
    done.
    
                      ┌──────────────────────────────────────────────────────────────┐
                      │                     Theory Level Setting                     │
                      ├──────────────────────────────────────────────────────────────┤
                      ├──────────────────────────────────────────────────────────────┤
                      │                    Geometry Optimization                     │
                      ├──────────────────────────────────────────────────────────────┤
                      │ Method            : gfn2                                     │
                      │ Opt Level         : normal                                   │
                      │ Solvation         : ALPB(ch2cl2)                             │
                      │ Etemp             : 300 K                                    │
                      │ Re-opt input      : False                                    │
                      ├──────────────────────────────────────────────────────────────┤
                      │                     Reactivity Features                      │
                      ├──────────────────────────────────────────────────────────────┤
                      │ Method            : gfn2                                     │
                      │ Solvation         : ALPB(ch2cl2)                             │
                      │ AFO Gap tolerance : 0.01 eV                                  │
                      │ Exclude elements  : [1]                                      │
                      ├──────────────────────────────────────────────────────────────┤
                      │                     Reaction Simulation                      │
                      ├──────────────────────────────────────────────────────────────┤
                      │ Method            : harmonic FF                              │
                      │ Etemp             : 1000 K                                   │
                      ├──────────────────────────────────────────────────────────────┤
                      │                     Single-Point Energy                      │
                      ├──────────────────────────────────────────────────────────────┤
                      │ Method            : gxtb                                     │
                      └──────────────────────────────────────────────────────────────┘
    
    Start doing mechanical analysis with given structues:  0
    
    Following force tranmission paths are considered:
    [21, 19, 18, 17, 16, 15, 4, 3, 0, 5, 7, 8, 9, 10, 11, 13]
    
              --------------------------------------------------------------------------------
                                             Mechano Assessment
              --------------------------------------------------------------------------------
               #      atl       mbo    uhf    deltaE[eh]    Fmax[nN]     k       a[1/bohr]
              --------------------------------------------------------------------------------
               0     [3, 0]    0.9384   0       0.0819       1.7167    0.0424      0.5088
              --------------------------------------------------------------------------------
    
                                     *********************************
                                     * Exploration Summary of Node 0 *
                                     *********************************
    
      ΔE = Esp_gas(P)-Esp_gas(R), unit given in kcal/mol.
      ------------------------------------------------------------------------------------------------
       #   RXNID   Scan_wd       R        Atom Pair         P          ΔE‡         ΔE             Note
      ------------------------------------------------------------------------------------------------
       0     0       0_0        [0]        [[3, 0]]       [1, 1]       inf        24.2    new rxn
       1     0       0_1        [0]        [[3, 0]]       [1, 1]       inf        24.2    rxn exists
      ------------------------------------------------------------------------------------------
    
    ***** Job finished in 0 hour 0 mins 19.07 sec *****
                  

    ```
    After running, each indexed reaction will be written to the `reaction.json` file, and each indexed structure will be written to the `struc.json` file.

