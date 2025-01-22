import pandas as pd


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