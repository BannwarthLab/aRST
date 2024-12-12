""""
"struc":{
    "struc_ID":{
        "struc_info":{
            "general":{
                "formula":str
                "mb":str
                "charge":int
                "multiplicity":int
                "nel":int
                "coord":dataframe
            }
            "energy":{
                "sp":float
                "sp_gfn":float
                "trv":float
                "solv":float
                "G:float
            }
            "reactivity":{
                "ml":dataframe
                "bo":dataframe
        }
        "reaction_path":{
            "as_reactant":list([reactionid])
            "as_product":list
        }
    }
}
"""
import pandas as pd

"""
"reaction":{
    "reaction_ID" :{
        "formula" :str
        "reaction_path" :list
        "struc_info" :{
             "reactants" :{
                list([reactionid])
            }
             "products":{
                list([reactionid])
            }
        }
        "energy_singlecycle":{ 
            sp_gfn2 :float
            sp_dft :float
            G :float
        }
        "energy_total":{ 
            sp_gfn2 :float
            sp_dft :float
            G :float
        }
    }
}
"""



class DictionaryInspector:
    def __init__(self, data):
        self.data = data

    def get_value(self, target_key, default_value=None):
        keys = target_key.split('.')
        current_data = self.data

        for key in keys:
            if isinstance(current_data, dict) and key in current_data:
                current_data = current_data[key]
            elif isinstance(current_data, list):
                try:
                    index = int(key)
                    current_data = current_data[index]
                except (ValueError, IndexError):
                    return current_data
            else:
                return default_value

        return current_data

    # def find_strucid_from_mb(self,mb):
    #     if mb is not None:
    #         out = []
    #         for key in self.data.keys():
    #             tmp_mb = self.get_value(f"{key}.struc_info.general.mb")
    #             if tmp_mb ==mb:
    #                 out.append(key)
    #         return out
    #     else:
    #         return None
    def set_value(self, target_key, new_value):
        keys = target_key.split('.')
        current_data = self.data

        for key in keys[:-1]:
            if isinstance(current_data, dict):
                if key not in current_data:
                    current_data[key] = {}
                current_data = current_data[key]
            else:
                raise KeyError(f"Invalid key '{key}' in the dictionary.")

        last_key = keys[-1]
        current_data[last_key] = new_value


    # def give_id(self):
    #     current_datanr = len(self.data.keys())
    #     return str(current_datanr + 1).zfill(4)

    # def update_strucrec_values(self, dict2,cal_G=False):
    #     struc_idl = dict2.keys()
    #     dict2_keys = get_deepest_keys(dict2)
    #     di2 = DictionaryInspector(dict2)
    #
    #
    #     def should_update(struc_id):
    #         # if lower energy, update all values
    #         egfn_old = self.get_value(f"{struc_id}.struc_info.energy.sp_gfn2")
    #         edft_old = self.get_value(f"{struc_id}.struc_info.energy.sp_dft")
    #
    #         egfn_new = di2.get_value(f"{struc_id}.struc_info.energy.sp_gfn2")
    #         edft_new = di2.get_value(f"{struc_id}.struc_info.energy.sp_dft")
    #
    #         if egfn_old is None and edft_old is None:
    #             return True
    #         elif edft_old is not None and edft_new is not None and edft_old > edft_new:
    #             return True
    #         elif egfn_old is not None and egfn_new is not None and egfn_old > egfn_new:
    #             return True
    #         else:
    #             return False
    #
    #     tmp_d = {i:should_update(i) for i in struc_idl}
    #     for i in dict2_keys:
    #         value = di2.get_value(i)
    #         strucdid = i.split(".")[0]
    #         if value is not None and tmp_d[strucdid]:
    #             oldvalue = self.get_value(i)
    #             if isinstance(value,list):
    #                 if oldvalue is not None:
    #                     value += oldvalue
    #                     value = list(set(value))
    #             self.set_value(i, value)
    #             if cal_G and self.get_value(f"{strucdid}.struc_info.energy.G") is None:
    #                 if self.get_value(f"{strucdid}.struc_info.energy.sp_dft") is None:
    #                     eg = self.get_value(f"{strucdid}.struc_info.energy.sp_dft")+\
    #                          self.get_value(f"{strucdid}.struc_info.energy.trv")+\
    #                          self.get_value(f"{strucdid}.struc_info.energy.solv")
    #                 else:
    #                     eg = self.get_value(f"{strucdid}.struc_info.energy.sp_gfn2") + \
    #                          self.get_value(f"{strucdid}.struc_info.energy.trv") + \
    #                          self.get_value(f"{strucdid}.struc_info.energy.solv")
    #                 self.set_value(f"{strucdid}.struc_info.energy.G", eg)
    #
    # def update_reactionrec_values(self, dict2):
    #
    #     dict2_keys = get_deepest_keys(dict2)
    #     di2 = DictionaryInspector(dict2)
    #     def should_update(key):
    #         if "energy" in key.split("."):
    #             if erec_gfn is None or erec_gfn > di2.get_value(f"{reaction_id}.energy.sp_gfn2"):
    #                 return True
    #             elif erec_dft is None or erec_dft > di2.get_value(f"{reaction_id}.energy.sp_dft"):
    #                 return True
    #             elif erec_g is None or erec_g > di2.get_value(f"{reaction_id}.energy.G"):
    #                 return True
    #         else:
    #             return True
    #
    #     for i in dict2_keys:
    #         reaction_id = i.split('.')[0]
    #         erec_gfn = self.get_value(f"{reaction_id}.energy.sp_gfn2")
    #         erec_dft = self.get_value(f"{reaction_id}.energy.sp_dft")
    #         erec_g = self.get_value(f"{reaction_id}.energy.G")
    #         if should_update(i):
    #             value = di2.get_value(i)
    #             if isinstance(value,list):
    #                 oldvalue = self.get_value(i)
    #                 oldvalue = [] if oldvalue is None else oldvalue
    #                 if "reaction_path" in i.split("."):
    #                     self.set_value(i, list(set(value+oldvalue)))
    #                 else:
    #                     self.set_value(i, value + oldvalue)
    #             else:
    #                 self.set_value(i, value)
    #
    #
    #
    # def set_reaction_values(self, STRUC_rec, reactants, products,unreactants,
    #                         reactionid, reactionpath,at_l,gap,onlyenergy=False):
    #     if not onlyenergy:
    #         name_reactants = [try_get_formular(STRUC_rec,i) for i in reactants if try_get_formular(STRUC_rec,i) is not None]
    #         name_products = [try_get_formular(STRUC_rec,i) for i in products if try_get_formular(STRUC_rec,i) is not None]
    #
    #         if all(name_reactants) is not None and all(name_products) is not None:
    #             str1 = ", ".join(name_reactants)
    #             str2 = ", ".join(name_products)
    #             self.set_value(f"{reactionid}.formula", f"{str1} to {str2}")
    #         tmp = self.get_value(f"{reactionid}.reaction_path")
    #         if tmp is None:
    #             if isinstance(reactionpath,list):
    #                 self.set_value(f"{reactionid}.reaction_path", list(set(reactionpath)))
    #             else:
    #                 self.set_value(f"{reactionid}.reaction_path", [reactionpath])
    #         else:
    #             if isinstance(reactionpath,list):
    #                 self.set_value(f"{reactionid}.reaction_path", tmp+reactionpath)
    #             else:
    #                 tmp.append(reactionpath)
    #                 tmp = list(set(tmp))
    #                 self.set_value(f"{reactionid}.reaction_path", tmp)
    #         self.set_value(f"{reactionid}.struc_info.reactants", reactants)
    #         self.set_value(f"{reactionid}.struc_info.products", products)
    #         if not unreactants:
    #             self.set_value(f"{reactionid}.struc_info.unreactants", None)
    #         else:
    #             self.set_value(f"{reactionid}.struc_info.unreactants", unreactants)
    #
    #         if at_l is not None:
    #             self.set_value(f"{reactionid}.atl", at_l)
    #         if at_l is not None:
    #             self.set_value(f"{reactionid}.gap", gap)
    #     for item in ['sp_gfn2', 'sp_dft', 'G']:
    #         item_values = [STRUC_rec.get_value(f"{i}.struc_info.energy.{item}") for i in reactants + products]
    #
    #         # Check for valid values
    #         if all(value is not None and value != float('inf') for value in item_values):
    #             results = sum([STRUC_rec.get_value(f"{i}.struc_info.energy.{item}") for i in products]) - \
    #                       sum([STRUC_rec.get_value(f"{i}.struc_info.energy.{item}") for i in reactants])
    #
    #             self.set_value(f"{reactionid}.energy.{item}", results*627.205)
    #         else:
    #             self.set_value(f"{reactionid}.energy.{item}", None)
    #
    #
    #
    # def find_sturcid_with_reactionpath(self, rpath, stype="as_reactant"):
    #     allkeys = get_deepest_keys(self.data)
    #     targetkeys = [i for i in allkeys if stype in i]
    #     struc_idl = [i.split(".")[0] for i in targetkeys if rpath in self.get_value(i)]
    #     return struc_idl
    #
    #
    # def update_struc_d_after_mp(self,dict2):
    #     di2 = DictionaryInspector(dict2)
    #     out_d = {}
    #     for strucid in dict2.keys():
    #         if strucid not in self.data.keys(): # new struc
    #             mb = dict2[strucid]["struc_info"]["general"]["mb"]
    #             if mb is None:
    #                 from aRST.script.struc import StrucOperation
    #                 coord = dict2[strucid]["struc_info"]["general"]["coord"]
    #                 mb = StrucOperation(coord, 1, 0).get_molbarcode(opt=False)
    #
    #             exist_id = self.find_strucid_from_mb(mb)
    #             exist_id_l = [k for k in exist_id if not "ref" in k.split(".")[0]]
    #
    #             if not exist_id_l:
    #                 exist_id = None
    #             else:
    #                 exist_id = exist_id_l[0].split(".")[0]
    #
    #
    #             if exist_id and exist_id is not None: # repeated struc
    #                 as_rl = di2.get_value(f"{strucid}.reaction_path.as_reactant")
    #                 as_pl = di2.get_value(f"{strucid}.reaction_path.as_product")
    #                 oldvalue1 = self.get_value(f"{exist_id}.reaction_path.as_reactant")
    #                 oldvalue2 = self.get_value(f"{exist_id}.reaction_path.as_product")
    #
    #                 if as_rl is not None:
    #                     if oldvalue1 is not None:
    #                         oldvalue1.extend(as_rl)
    #                     else:
    #                         oldvalue1 = as_rl
    #                     self.set_value(f"{exist_id}.reaction_path.as_reactant",oldvalue1)
    #
    #                 if as_pl is not None:
    #                     if oldvalue2 is not None:
    #                         oldvalue2.extend(as_pl)
    #                     else:
    #                         oldvalue2 = as_pl
    #                     self.set_value(f"{exist_id}.reaction_path.as_product",oldvalue2)
    #                 out_d[strucid] = exist_id
    #             else:
    #                 tmp_d = {strucid: dict2[strucid]}
    #                 self.update_strucrec_values(tmp_d)
    #     return out_d
    #
    #
    # def update_reaction_d_after_mp(self,dict2,out_d=None):
    #     di2 = DictionaryInspector(dict2)
    #     for reactionid in dict2.keys():
    #         current_reactants_idl = di2.get_value(f"{reactionid}.struc_info.reactants") or []
    #         current_unreactants_idl = di2.get_value(f"{reactionid}.struc_info.unreactants") or []
    #         current_products_idl = di2.get_value(f"{reactionid}.struc_info.products") or []
    #
    #
    #         tmp_l = [strucid if strucid not in out_d.keys() else out_d[strucid] for strucid in current_products_idl]
    #         di2.set_value(f"{reactionid}.struc_info.products", tmp_l)
    #
    #
    #         existid = self.test_if_reaction_exist(current_reactants_idl, current_unreactants_idl, tmp_l)
    #         if existid is not None:
    #             oldreactionpath = self.get_value(f"{existid}.reaction_path") or []
    #
    #             oldreactionpath.extend(di2.get_value(f"{reactionid}.reaction_path"))
    #             self.set_value(f"{existid}.reaction_path", oldreactionpath)
    #
    #
    #         else:
    #             self.update_reactionrec_values({reactionid: dict2[reactionid]})
    #
    #
    #
    #     return
    #
    #
    #
    # def test_if_reaction_loop(self,reactionid):
    #     # reactants for this cycle can't be identical to existing reactants
    #     products_idl = self.get_value(f"{reactionid}.struc_info.products") or []
    #     unreactants_idl = self.get_value(f"{reactionid}.struc_info.unreactants") or []
    #
    #
    #     reactants_idl = products_idl+unreactants_idl
    #     for current_reactionid, items in self.data.items():
    #         if items.get("struc_info", None) is not None:
    #             current_reactants = items.get("struc_info").get("reactants") or []
    #             # current_products = items.get("struc_info").get("products")
    #             current_unreactants = items.get("struc_info").get("unreactants") or []
    #             current_fullreactants = current_reactants+current_unreactants
    #
    #             if current_fullreactants and reactants_idl:
    #                 if Counter(reactants_idl) == Counter(current_fullreactants):
    #                     return True
    #     else:
    #         return False
    #
    #
    #
    #
    # def test_if_reaction_exist(self,reactants_idl, unreactants_idl, products_idl):
    #     outid = None
    #     for reactionid,items in self.data.items():
    #         current_reactants = self.get_value(f"{reactionid}.struc_info.reactants") or []
    #         current_unreactants_idl = self.get_value(f"{reactionid}.struc_info.unreactants") or []
    #         current_products_idl = self.get_value(f"{reactionid}.struc_info.products") or []
    #
    #         full_reactants = current_reactants+current_unreactants_idl
    #         target_reactants = reactants_idl+unreactants_idl if unreactants_idl is not None else reactants_idl
    #         if all([full_reactants,current_products_idl,target_reactants,products_idl]):
    #
    #             if Counter(full_reactants)==Counter(target_reactants) and \
    #                     Counter(current_products_idl) == Counter(products_idl):
    #                 outid = reactionid
    #                 continue
    #     return outid





def convert_df2json(d):
    for key, value in d.items():
        if isinstance(value, dict):
            convert_df2json(value)
        elif isinstance(value, pd.DataFrame):
            d[key] = value.to_json()
        else:
            pass
    return d




import json
class DataFrameEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, pd.DataFrame):
            return obj.to_json(orient='split', index=False)
        return super().default(obj)




def get_deepest_keys(dictionary, parent_key="", deepest_keys=None):
    if deepest_keys is None:
        deepest_keys = []

    for key, value in dictionary.items():
        current_key = f"{parent_key}.{key}" if parent_key else key

        if isinstance(value, dict):
            get_deepest_keys(value, current_key, deepest_keys)
        else:
            deepest_keys.append(current_key)

    return deepest_keys


def convert_json2dict(jsondata):
    d = {}
    for struc in jsondata.keys():
        d[struc] = {}
        for item in list(jsondata[struc].keys()):

            try:
                df = pd.read_json(jsondata[struc][item])
                d[struc][item] = df
            except:
                d[struc][item] = jsondata[struc][item]
    return d


class JSONInspector:
    def __init__(self, json_data):
        self.json_data = json_data

    def get_value(self, target_key, default_value=None):
        keys = target_key.split('.')
        current_data = self.json_data

        for key in keys:
            if isinstance(current_data, dict) and key in current_data:
                current_data = current_data[key]
            elif isinstance(current_data, list):
                try:
                    index = int(key)
                    current_data = current_data[index]
                except (ValueError, IndexError):
                    return default_value
            else:
                return default_value
        return current_data


