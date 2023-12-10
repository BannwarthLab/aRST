import os
import shutil
import sys
import time
from collections import Counter

import pandas as pd
from molbar.barcode import get_molbar_from_file
from toolbox import create_folder, write_xyz, seperate, creat_randomfolder


class StrucOperation:
    def __init__(self, coord, m, charge):
        self.coord = coord
        self.m = m
        self.charge = charge

    def check_if_fall_back(self, reactants_mb):
        pwd = os.getcwd()
        random_wd = creat_randomfolder(pwd)
        write_xyz(self.coord, 'coord')
        seperate()
        fragment_mb = []
        for k in os.listdir(f'{random_wd}'):
            if os.path.isdir(f'{random_wd}/{k}') and str.isdigit(k):
                self.coord = pd.read_table(f'{random_wd}/{k}/coord.xyz', header=None, skiprows=2, sep='\s+')
                fragment_mb.append(self.get_molbarcode())

        if Counter(reactants_mb) == Counter(fragment_mb):
            return True
        else:
            os.chdir(pwd)
            shutil.rmtree(f'{random_wd}', ignore_errors=True)
            return False

    def get_nfragement(self):
        pwd = os.getcwd()
        wd = create_folder(os.path.join(pwd,"sep"))
        write_xyz(self.coord, 'coord')
        seperate()
        l = [k for k in os.listdir('.') if os.path.isdir(f'.') and str.isdigit(k)]
        dir_l = [os.path.join(wd,k) for k in l]
        if len(l) == 1:
            os.chdir(pwd)
            shutil.rmtree(f'{wd}', ignore_errors=True)
            return dir_l # if nfragment =1, dir_l no meaning, only len counts
        else:
            return dir_l


    def compare_structs(self, Struc2):
        # compare charge,m,mb
        if self.charge != Struc2.charge:
            return False
        elif self.m != Struc2.m:
            return False
        elif self.get_molbarcode() != Struc2.get_molbarcode():
            return False
        else:
            return True

    def get_molbarcode(self):
        pwd = os.getcwd()
        random_wd = creat_randomfolder(pwd)
        write_xyz(self.coord, 'coord')
        mb = get_molbar_from_file("coord.xyz")
        os.chdir(pwd)
        shutil.rmtree(f'{random_wd}', ignore_errors=True)
        return mb