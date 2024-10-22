# aRST: an Automated Reaction Search Tool

This package provides as an automated tool for modeling chemical reaction paths, searching thermodynamic-allowed products based on the given electronic structures. Other useful built-in functions include evaluating atomic reactive sites and generating the reactive complex.

## Getting Started
To clone and run this package you will need:

Clone
```bash
# Clone this repository
git clone https://github.com/xx
```

Create a virtual environment for python
```bash
python3 -m venv your_path
```

Activate this environment
```bash
source your_path/bin/activate
```

Install aRST
```bash
run setup.py
```

now you are able to used aRST as commandline tool or call its functions :)

## Commandlines for calling usages
Example of different usages are given below as well as can be found in test folder.

To start the calculation, the essential input consists of the .xyz file of the structure, the .toml file with detailed settings.

for my version of xtb and stda i can provide a complied binary address with
xtb
```
/home/guests/ychen/public/xtbml/_build/xtb
```
stda
```
/home/guests/ychen/public/stdaml/stda
```

### Reaction Pathes Exploration
Explore reaction pathes based on given reactants.
```
aRST in.toml --explore > record.log 2>&1 &
```

### HOAO and LUAO Infomation
Print out reactive information of given structures.
```
aRST in.toml -e -t > record.log 2>&1 &
```

### Alignment
Aligning two strucs with given atom index.
```
aRST in.toml --align > record.log 2>&1 &
```

