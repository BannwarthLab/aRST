import sys

from rdkit import Chem
from rdkit.Chem import Draw
from os.path import join, exists
from rdkit.Chem import rdDetermineBonds
import json
from rdkit.Chem import AllChem

wd = "/home/ychen/aRST/test/tmp"
reactants = [0, 1, 2, 3]
list_1 = [8,9]
list_2 = [4,5,6,7,10,11]
list_3 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 23, 26, 44, 72]
# list_3 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 25, 32, 33, 34, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 71, 72, 73, 77, 80, 81, 82, 83, 84, 85, 87, 88, 89, 90, 92, 93, 94, 95, 97, 101, 104, 105, 106, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 121, 122, 123, 124, 125, 126, 127, 128, 130, 131, 132, 133, 134, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 151, 152, 154, 155, 156, 157, 158, 159, 160, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174]


def get_Erxn(reaction_data, target_sid):
    matches = []
    for item in reaction_data:
        if "info" in item:
            rxninfo = item["info"]
            rxnid = rxninfo["rxnid"]
            products = rxninfo.get("products", [])
            Erxn = rxninfo.get("E_rxn", None)

            if target_sid in products:
                matches.append({
                    "rxnid": rxnid,
                    "E_rxn": Erxn,
                })

    if not matches:
        return [], None

    lowest = min(matches, key=lambda x: x["E_rxn"])

    return matches, lowest

struc_list0 = reactants
struc_list = list_2
# struc_list = [_ for _ in list_1 if _ not in reactants] # and _ not in list_2


print(struc_list)
with open(join(wd, "struc.json"), "r") as f:
    struc_data = json.load(f)

struc_map = {item["id"]: item for item in struc_data}

with open(join(wd, "reaction.json"), "r") as f:
    reaction_data = json.load(f)



mols = []
legends = []

for struc in struc_list:
    xyzwd = join(wd, "struc", str(struc), "gxtb", "coord.xyz")

    if not exists(xyzwd):
        print(f"Missing xyz: {struc}")
        continue

    general = struc_map[struc]["struc_info"]["general"]
    charge = int(general["charge"])
    m = int(general["multiplicity"])

    mol = Chem.MolFromXYZFile(xyzwd)

    if mol is None:
        print(f"Failed to read xyz: {struc}")
        continue

    # print(struc, charge, m)
    match, lowest_entry = get_Erxn(reaction_data,struc)

    if m == 1:
        rdDetermineBonds.DetermineBonds(mol, charge=charge)
        Chem.SanitizeMol(mol)
        status = "Lewis"
    else:
        rdDetermineBonds.DetermineConnectivity(mol)
        status = "radical"

    AllChem.Compute2DCoords(mol)

    mols.append(mol)

    if struc in reactants:
        legends.append(f"ID{struc} ({charge},{m})")
        continue
    if lowest_entry:
        legends.append(f"ID{struc}, ({charge},{m}), Erxn={lowest_entry['E_rxn']:.1f}")
    else:
        legends.append(f"ID{struc} ({charge},{m})")

img = Draw.MolsToGridImage(
    mols,
    molsPerRow=10,
    subImgSize=(300, 220),
    legends=legends,
    useSVG=False
)

img.save(join(wd, "layer3.png"))