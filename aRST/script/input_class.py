import re
import os
from typing import Dict, Optional, Literal, Union, List, Tuple
from pydantic import BaseModel, ConfigDict, field_validator, model_validator, Field
from enum import Enum


class MyBaseModel(BaseModel):
    model_config = ConfigDict(extra="forbid")

    @field_validator("*", mode="before")
    @classmethod
    def normalize_strings(cls, v):
        if isinstance(v, str):
            return v.strip().lower()
        return v


class Reactant(MyBaseModel):
    model_config = ConfigDict(extra="forbid")
    name: str
    charge: int
    multiplicity: int
    number: int = 1


class GeneralRead(MyBaseModel):
    model_config = ConfigDict(extra="forbid")
    dftopt: bool = False
    gfn2opt: bool = False
    jumprefine: bool = False


class Read(MyBaseModel):
    model_config = ConfigDict(extra="forbid")
    general: GeneralRead = GeneralRead()


class AFOMethod(str, Enum):
    gfn2 = "gfn2"
    stda = "stda"
    gxtb = "gxtb"


class MOMethod(str, Enum):
    orca = "orca"
    terachem = "terachem"


class ReactivitySearchingSetting(MyBaseModel):
    model_config = ConfigDict(extra="forbid")
    afo_method: AFOMethod = AFOMethod.gfn2
    MO_method: Optional[MOMethod] = None

    # only meaningful when MO_method == "orca"
    afo_orca_functional: str = "BHandHLYP"
    afo_orca_basis: str = "vdzp"
    afo_orca_cpcm: str = "ch2cl2"

    # only meaningful when MO_method == "terachem"
    afo_tc_functional: str = "BHandHLYP"
    afo_tc_basis: str = "vdzp"
    afo_tc_epsilon: float = 9.08

    @model_validator(mode="after")
    def validate_conditional_settings(self):
        # Rule: if STDA, MO_method must be provided
        if self.afo_method == AFOMethod.stda and self.MO_method is None:
            raise ValueError("MO_method must be provided when afo_method is 'stda'")

        orca_fields = {"afo_orca_functional", "afo_orca_basis", "afo_orca_cpcm"}
        tc_fields = {"afo_tc_functional", "afo_tc_basis", "afo_tc_epsilon"}

        # user explicitly set
        set_fields = set(self.model_fields_set)

        if self.MO_method == MOMethod.orca:
            # Forbid user setting any TeraChem-only fields
            bad = sorted(tc_fields & set_fields)
            if bad:
                raise ValueError(f"{bad} can only be set when MO_method='terachem'")

        elif self.MO_method == MOMethod.terachem:
            # Forbid user setting any ORCA-only fields
            bad = sorted(orca_fields & set_fields)
            if bad:
                raise ValueError(f"{bad} can only be set when MO_method='orca'")

        else:
            # MO_method is None: forbid setting any MO-specific fields
            bad = sorted((orca_fields | tc_fields) & set_fields)
            if bad:
                raise ValueError(f"{bad} requires MO_method to be set (orca/terachem)")

        return self


class MDSetting(MyBaseModel):
    model_config = ConfigDict(extra="forbid")
    time: int = 10
    step: int = 1
    temp: int = 300
    shake: int = 1
    dump: int = 50


class GeneralSearchingSetting(MyBaseModel):
    model_config = ConfigDict(extra="forbid")
    wd: str = Field(default_factory=os.getcwd)
    scanwd: Optional[str] = None
    strucwd: Optional[str] = None
    bufferwd: Optional[str] = None

    cycle: int = 1
    sp_method: Literal["dft", "gxtb"] = "dft"
    elim: int = 300
    tslim: int=100
    tseval:bool = False
    equivalent_afogap_tol: float = 0.01
    spr_mode: Literal[0, 1, 2] = 0
    scan_mode: Literal[0, 1] = 0
    sep_mode: Literal[0, 1] = 0
    md: bool = False
    md_setting: Optional[MDSetting] = None
    mb_level: Literal[0, 1] = 0
    mbcut_length:Optional[int] = None

    # optional inputs
    exclude_elements: List[int] = Field(default_factory=list)
    exclude_pairs: List[Tuple[int, int]] = Field(default_factory=list)

    @model_validator(mode="after")
    def handle_md_logic(self):
        if self.md:
            # if md is True but no settings provided → use defaults
            if self.md_setting is None:
                self.md_setting = MDSetting()
        else:
            # if md is False but settings provided → forbid
            if self.md_setting is not None:
                raise ValueError(
                    "md_setting provided but md is False"
                )
        return self

    @model_validator(mode="after")
    def set_default_paths(self):
        if self.scanwd is None:
            self.scanwd = os.path.join(self.wd, "scan")
        if self.strucwd is None:
            self.strucwd = os.path.join(self.wd, "struc")
        if self.bufferwd is None:
            self.bufferwd = os.path.join(self.wd, "buffer")
        return self

    @model_validator(mode="after")
    def set_mbcut_length(self):
        # for compare new struc
        if self.mbcut_length is None:
            self.mbcut_length = 5 if self.mb_level == "loose" else 7
        return self

class InterSearchingSetting(MyBaseModel):
    model_config = ConfigDict(extra="forbid")
    number: Union[int, Literal["all"]] = 0
    concerted_pair: bool = False
    atom_pair: List[Tuple[int, int]] = Field(default_factory=list)
    aligning_distance: int = 5
    maxovlap: bool = False
    loose_bond: bool = False


class IntraSearchingSetting(MyBaseModel):
    model_config = ConfigDict(extra="forbid")
    number: Union[int, Literal["all"]] = 0
    concerted_pair: bool = False
    atom_pair: List[Tuple[int, int]] = Field(default_factory=list)
    loose_bond: bool = False

    depotonation: bool = False
    potonation: bool = False


class DissoSearchingSetting(MyBaseModel):
    model_config = ConfigDict(extra="forbid")
    number: Union[int, Literal["all"]] = 0
    concerted_pair: bool = False
    atom_pair: List[Tuple[int, int]] = Field(default_factory=list)

    pulling_distance: Union[int, Literal["auto"]] = "auto"
    pulling_distance_add: Optional[float] = None

    depotonation: bool = False
    potonation: bool = False

    @model_validator(mode="after")
    def check_exclusive(self):
        set_fields = self.model_fields_set

        # If both explicitly set → error
        if "pulling_distance" in set_fields and "pulling_distance_add" in set_fields:
            raise ValueError(
                "Only one of 'pulling_distance' or 'pulling_distance_add' can be provided."
            )

        return self


class SearchingSetting(MyBaseModel):
    model_config = ConfigDict(extra="forbid")
    reactivity: ReactivitySearchingSetting = ReactivitySearchingSetting()
    general: GeneralSearchingSetting = GeneralSearchingSetting()
    inter: InterSearchingSetting = InterSearchingSetting()
    intra: IntraSearchingSetting = IntraSearchingSetting()
    disso: DissoSearchingSetting = DissoSearchingSetting()


class XTBSetting(MyBaseModel):
    model_config = ConfigDict(extra="forbid")
    method: str = "gfn2"
    optlevel: str = ""
    alpb: str = "ch2cl2"
    struc_etemp: int = 300
    scan_etemp: int = 1000
    k: int = 1
    step: int = 20


class ORCASetting(MyBaseModel):
    model_config = ConfigDict(extra="forbid")
    functional: str = "r2scan-3c"
    basis: str = ""
    dftopt: bool = False
    cpcm: str = "ch2cl2"
    ncpu: int = 1
    maxcore: int=1000


class InputModel(MyBaseModel):
    model_config = ConfigDict(extra="forbid")

    # geom section
    geom: Dict[str, Reactant]

    # read section
    read: Read = Read()

    # searching_setting section
    searching_setting: SearchingSetting = SearchingSetting()

    # xtb setting
    xtb_setting: XTBSetting = XTBSetting()

    # orca setting
    orca_setting: ORCASetting = ORCASetting()

    # check if geom name is coorect
    @field_validator("geom")
    @classmethod
    def only_reactant_keys(cls, geom: Dict[str, Reactant]):
        bad = [k for k in geom.keys() if not re.fullmatch(r"reactant\d+", k)]
        if bad:
            raise ValueError(f"Invalid geom keys (expected reactantN): {bad}")
        if not geom:
            raise ValueError("At least one reactant (reactant1, reactant2, ...) is required.")
        return geom
