import json
import os


def identity(*args):
    if len(args) == 1:
        return args[0]
    return args


class config:

    def __init__(self, var, type=identity, ignore_case=False, default=None):
        self.var = var
        self.type = type
        self.ignore_case = ignore_case
        self.default = default


class json_config:

    def __init__(self):
        self.configs = {}
        self.config_set = set()
        self.defaults = {}

    def add_config(self, var, type=identity, ignore_case=True, synonyms=[], default=None):
        temp = config(var, type, ignore_case, default)
        self.config_set.add(temp)
        self.defaults[var] = default
        var = var.lower() if(ignore_case) else var
        self.configs[var] = temp
        for syn in synonyms:
            syn = syn.lower if(ignore_case) else syn
            self.configs[syn] = temp

    def load_config(self, json_file, no_fail=True):
        if(not os.path.isfile(json_file)):
            return self.defaults
        f = open(json_file)
        data = json.load(f)
        f.close()
        assert isinstance(data, dict)
        ret = {}
        for key in data:
            lower_key = key.lower()
            if(key in self.configs):
                ret[self.configs[key].var] = self.configs[key].type(data[key])
            elif(lower_key in self.configs and self.configs[lower_key].ignore_case):
                ret[self.configs[lower_key].var] = self.configs[lower_key].type(data[key])
            else:
                raise Exception('Unrecognized config option "{0!s}" found within {1!s}'.format(key, json_file))
        temp_defaults = dict(self.defaults)
        temp_defaults.update(ret)
        return temp_defaults


if(__name__ == '__main__'):

    def test_json_config():
        test_file_name = 'test_config.json'
        test_dict = {'Opt1': 1}
        f = open(test_file_name, 'w')
        json.dump(test_dict, f)
        f.close()
        conf = json_config()
        conf.add_config('opt1')
        conf.add_config('opt2')
        data = conf.load_config(test_file_name)
        print(data)
        os.remove(test_file_name)
