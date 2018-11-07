# different implemention of nested dict in Python
from collections import defaultdict


# from: http://stackoverflow.com/questions/635483/what-is-the-best-way-to-implement-nested-dictionaries-in-python
class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


# from: http://stackoverflow.com/questions/635483/what-is-the-best-way-to-implement-nested-dictionaries-in-python
class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value


# from: http://stackoverflow.com/questions/651794/whats-the-best-way-to-initialize-a-dict-of-dicts-in-python
def ddict():
    return defaultdict(ddict)


# from @lpp1985 https://github.com/lpp1985
class Ddict(defaultdict, dict):
    def __init__(self):
        defaultdict.__init__(self, Ddict)

    def __repr__(self):
        return dict.__repr__(self)

# test
'''
for c in [AutoVivification, Vividict, ddict, Ddict]:
    print("\n%s\n%s" % ('=' * 78, c))

    d = c()
    d[1][2][3] = 4
    d['a']['b']['c'] = 'd'

    print(d)
'''
