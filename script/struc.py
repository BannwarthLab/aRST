import sys

from callsys import call_orca_dft, call_xtb, call_xtb_ml
from operation import creat_folder, write_xyz



class Struc:
    def __init__(self, coord, m, charge):
        self.coord = coord
        # self.nat = nat
        # self.nel = nel
        self.m = m
        self.charge = charge

    def runxtb(self, wd, **xtb_setting):
        creat_folder(f'{wd}')
        write_xyz(self.coord, 'coord')
        nat = self.coord.shape[0]
        if nat > 1:
            rc = call_xtb(f'{wd}', self.charge, caltyp='--opt', uhf=self.m - 1, **xtb_setting)
        else:
            rc = call_xtb(f'{wd}', self.charge, uhf=self.m - 1, **xtb_setting)
        return rc


    def runorcasp(self, wd, caltyp, **orca_setting):
        creat_folder(f'{wd}')
        write_xyz(self.coord, 'coord')
        rc = call_orca_dft(f'{wd}', self.charge, self.m, caltyp=caltyp, **orca_setting)
        return rc

    def runxtbml(self, wd, **xtb_setting):
        creat_folder(f'{wd}')
        write_xyz(self.coord, 'coord')
        call_xtb_ml(f'{wd}', self.charge, self.m - 1, **xtb_setting)

