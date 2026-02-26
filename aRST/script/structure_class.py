from typing import Optional
from pydantic import BaseModel,ConfigDict
import pandas as pd

class GeneralStrucInfo(BaseModel):
    coordinate: Optional[list] = None
    ele_l: Optional[list] = None
    charge: Optional[float] = None
    multiplicity: Optional[int] = None
    number:Optional[int] = None
    nel:Optional[int] = None
    atcharge:Optional[list] = None
    mb:Optional[str] = None
    atomic_connection:Optional[list] = None

class ReactivityStrucInfo(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    afo: Optional[pd.DataFrame] = None
    bo:Optional[pd.DataFrame] = None

class EnergyStrucInfo(BaseModel):
    sp: Optional[float] = None


class StrucInfo(BaseModel):
    general: GeneralStrucInfo = GeneralStrucInfo()
    energy: EnergyStrucInfo = EnergyStrucInfo()
    reactivity:ReactivityStrucInfo=ReactivityStrucInfo()


class Structure(BaseModel):
    id: int
    struc_info: StrucInfo = StrucInfo()