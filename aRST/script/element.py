element = {
    'H': 1,
    'He': 2,
    'Li': 3,
    'Be': 4,
    'B': 5,
    'C': 6,
    'N': 7,
    'O': 8,
    'F': 9,
    'Ne': 10,
    'Na': 11,
    'Mg': 12,
    'Al': 13,
    'Si': 14,
    'P': 15,
    'S': 16,
    'Cl': 17,
    'CL': 17,
    'Ar': 18,
    'K': 19,
    'Ca': 20,
    'Sc': 21,
    'Ti': 22,
    'V': 23,
    'Cr': 24,
    'Mn': 25,
    'Fe': 26,
    'Co': 27,
    'Ni': 28,
    'Cu': 29,
    'Zn': 30,
    'Ga': 31,
    'Ge': 32,
    'As': 33,
    'Se': 34,
    'Br': 35,
    'Kr': 36,
    'Rb': 37,
    'Sr': 38,
    'Y': 39,
    'Zr': 40,
    'Nb': 41,
    'Mo': 42,
    'Tc': 43,
    'Ru': 44,
    'Rh': 45,
    'Pd': 46,
    'Ag': 47,
    'Cd': 48,
    'In': 49,
    'Sn': 50,
    'Sb': 51,
    'Te': 52,
    'I': 53,
    'Xe': 54,
    'Cs': 55,
    'Ba': 56,
    'La': 57,
    'Pb': 82,
    'Bi': 83
}

element2 = {v: k for k, v in element.items()}

r0 = {'H': 0.32, 'D': 0.32, 'He': 0.46, 'Li': 1.20, 'Be': 0.94, 'B': 0.77, 'C': 0.75, 'N': 0.71, 'O': 0.63,
      'F': 0.64, 'Ne': 0.67, 'Na': 1.40, 'Mg': 1.25, 'Al': 1.13, 'Si': 1.04, 'P': 1.10, 'S': 1.02,
      'CL': 0.99, 'Ar': 0.96, 'K': 1.76, 'Ca': 1.54, 'Sc': 1.48, 'Ti': 1.36, 'V': 1.34, 'Cr': 1.22,
      'Mn': 1.19, 'Fe': 1.16, 'Co': 1.11, 'Ni': 1.10, 'Cu': 1.12, 'Zn': 1.18, 'Ga': 1.24, 'Ge': 1.21,
      'As': 1.21, 'Se': 1.16, 'Br': 1.14, 'Kr': 1.17, 'Rb': 2.10, 'Sr': 1.85, 'Y': 1.63, 'Zr': 1.54,
      'Nb': 1.47, 'Mo': 1.38, 'Tc': 1.28, 'Ru': 1.25, 'Rh': 1.25, 'Pd': 1.20, 'Ag': 1.28, 'Cd': 1.36,
      'In': 1.42, 'Sn': 1.40, 'Sb': 1.40, 'Te': 1.36, 'I': 1.33, 'Xe': 1.31, 'Cs': 2.32, 'Ba': 1.96,
      'La': 1.80, 'Ce': 1.63, 'Pr': 1.76, 'Nd': 1.74, 'Pm': 1.73, 'Sm': 1.72, 'Eu': 1.68, 'Gd': 1.69,
      'Tb': 1.68, 'Dy': 1.67, 'Ho': 1.66, 'Er': 1.65, 'Tm': 1.64, 'Yb': 1.70, 'Lu': 1.62, 'Hf': 1.52,
      'Ta': 1.46, 'W': 0.95 * 1.37, 'Re': 1.31, 'Os': 1.29, 'Ir': 1.22, 'Pt': 1.23, 'Au': 1.24, 'Hg': 1.33,
      'Tl': 1.44, 'Pb': 1.44, 'Bi': 1.51, 'Po': 1.45, 'At': 1.47, 'Rn': 1.42, 'Fr': 2.23, 'Ra': 2.01,
      'Ac': 1.86, 'Th': 1.75, 'Pa': 1.69, 'U': 1.70, 'Np': 1.71, 'Pu': 1.72, 'Am': 1.66, 'Cm': 1.66,
      'Bk': 1.68, 'Cf': 1.68, 'Es': 1.65, 'Fm': 1.67, 'Md': 1.73, 'No': 1.76, 'Lr': 1.61, 'Rf': 1.57,
      'Db': 1.49, 'Sg': 1.43, 'Bh': 1.41, 'Hs': 1.34, 'Mt': 1.29, 'Ds': 1.21, 'Rg': 1.21, 'Cn': 1.22,
      'Nh': 1.36, 'Fl': 1.43, 'Mc': 1.62, 'Lv': 1.75, 'Ts': 1.65, 'Og': 1.57, 'Cl': 0.99}
scaled_r0 = {k: v*4/3 for k, v in element.items()}

valence_electrons = {
    # Group 1 - Alkali Metals
    "H": 1, "Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1, "Fr": 1,

    # Group 2 - Alkaline Earth Metals
    "Be": 2, "Mg": 2, "Ca": 2, "Sr": 2, "Ba": 2, "Ra": 2,

    # Group 3 (commonly 3 valence electrons, though transition metals are more complex)
    "Sc": 3, "Y": 3, "La": 3, "Ac": 3,

    # Group 4
    "Ti": 4, "Zr": 4, "Hf": 4, "Rf": 4,

    # Group 5
    "V": 5, "Nb": 5, "Ta": 5, "Db": 5,

    # Group 6
    "Cr": 6, "Mo": 6, "W": 6, "Sg": 6,

    # Group 7
    "Mn": 7, "Tc": 7, "Re": 7, "Bh": 7,

    # Group 8
    "Fe": 8, "Ru": 8, "Os": 8, "Hs": 8,

    # Group 9
    "Co": 9, "Rh": 9, "Ir": 9, "Mt": 9,

    # Group 10
    "Ni": 10, "Pd": 10, "Pt": 10, "Ds": 10,

    # Group 11
    "Cu": 1, "Ag": 1, "Au": 1, "Rg": 1,

    # Group 12
    "Zn": 2, "Cd": 2, "Hg": 2, "Cn": 2,

    # Group 13 - Boron Group
    "B": 3, "Al": 3, "Ga": 3, "In": 3, "Tl": 3, "Nh": 3,

    # Group 14 - Carbon Group
    "C": 4, "Si": 4, "Ge": 4, "Sn": 4, "Pb": 4, "Fl": 4,

    # Group 15 - Nitrogen Group
    "N": 5, "P": 5, "As": 5, "Sb": 5, "Bi": 5, "Mc": 5,

    # Group 16 - Oxygen Group (Chalcogens)
    "O": 6, "S": 6, "Se": 6, "Te": 6, "Po": 6, "Lv": 6,

    # Group 17 - Halogens
    "F": 7, "Cl": 7, "Br": 7, "I": 7, "At": 7, "Ts": 7,

    # Group 18 - Noble Gases
    "He": 2, "Ne": 8, "Ar": 8, "Kr": 8, "Xe": 8, "Rn": 8, "Og": 8,

    # Lanthanides (Variable valence electrons, but typically 3)
    "Ce": 3, "Pr": 3, "Nd": 3, "Pm": 3, "Sm": 3, "Eu": 3, "Gd": 3,
    "Tb": 3, "Dy": 3, "Ho": 3, "Er": 3, "Tm": 3, "Yb": 3, "Lu": 3,

    # Actinides (Variable, but typically 3)
    "Th": 3, "Pa": 3, "U": 3, "Np": 3, "Pu": 3, "Am": 3, "Cm": 3,
    "Bk": 3, "Cf": 3, "Es": 3, "Fm": 3, "Md": 3, "No": 3, "Lr": 3
}
