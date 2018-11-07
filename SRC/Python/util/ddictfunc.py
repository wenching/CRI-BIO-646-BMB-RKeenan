from collections import defaultdict
import json
import pprint
import sys

'''
[One-line Tree in Python](https://gist.github.com/hrldcpr/2012250)
'''
def ddicts(): return defaultdict(ddicts)


def subset(d, key_list, invert = False, join=['inner', 'left'][0]):
    if invert:
        return {k: d[k] for k in d if k not in key_list}
    else:
        if(join == "inner"): return {k: d[k] for k in key_list if k in d}
        if(join == "left"): return {k: d.get(k, None) for k in key_list}

'''
[How to convert defaultdict of defaultdicts to dict of dicts?](https://stackoverflow.com/questions/26496831/how-to-convert-defaultdict-of-defaultdicts-of-defaultdicts-to-dict-of-dicts-o)
'''
def ddicts_2_dict(d): 
    if isinstance(d, defaultdict):
        d = {k: ddicts_2_dict(v) for k, v in d.items()}
    return d


def print_ddicts(ddicts, test_exit = True):
    print(json.dumps(ddicts))

    if(test_exit): sys.exit()


def pprint_ddicts(ddicts, key_list = None, test_exit = True):
    if key_list is None:
        pprint.pprint(ddicts_2_dict(ddicts))
    else:
        pprint.pprint(subset(ddicts_2_dict(ddicts), key_list))

    if(test_exit): sys.exit()
