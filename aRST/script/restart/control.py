import pickle
import os
import aRST.script.initialization as storage
from os.path import join, exists

class ProgramRunner:
    # a contronal for saving the breaking point of each complete reaction searching simulations
    def __init__(self):
        self.rxnnetwork = None
        self.allstruc = None
        self.results = {}



    def checkpoint(self,fn="restart.pkl"):
        path = join(storage.setting.searching_setting.general.wd,fn)
        tmp = path + ".tmp"
        with open(tmp, "wb") as f:
            pickle.dump(self, f)
        os.replace(tmp, path)

    @classmethod
    def load(cls, fn="restart.pkl"):
        path = join(storage.setting.searching_setting.general.wd, fn)
        if os.path.exists(path):
            with open(path, "rb") as f:
                return pickle.load(f)
        return cls()

    def run(self,tasks):
        for i, task in enumerate(tasks):
            if i < self.step:
                continue

            result = self.run_one(task)

            self.results[i] = result
            self.step = i + 1

            self.checkpoint()
