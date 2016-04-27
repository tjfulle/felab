from numpy import *
from pandas import MultiIndex, DataFrame

SCALAR = 'Scalar'
VECTOR = 'Vector'
SYMTENSOR = 'Symmetric Tensor'
TENSOR = 'Tensor'

class DataWharehouse(object):

    def __init__(self, dimensions, node_labels, node_vars, elem_blocks):

        # NODE DATA
        index = MultiIndex.from_product([['one','two'], node_labels])
        self.node_df = DataFrame(index=index)
        self.node_key_map = {}
        for node_var in node_vars:
            name, type = node_var[:2]
            idata = 0. if len(node_var) <= 2 else node_var[2]
            if type == SCALAR:
                self.node_df[name] = idata
            elif type == VECTOR:
                for c in 'XYZ'[:dimensions]:
                    key = '{0}_{1}'.format(name, c)
                    self.node_df[key] = idata
                    self.node_key_map.setdefault(name, []).append(key)
            else:
                raise TypeError('UNKNOWN NODE DATA TYPE')

        # ELEMENT DATA

    def _get_key_labels(self, arg):
        assert len(arg) == 2
        key, labels = arg[0], arg[1]
        if isinstance(key, (str,)):
            key = self._expand_key(key)
        return key, labels

    def _expand_key(self, arg):
        return self.node_key_map[arg]

    def __getitem__(self, arg):
        if arg in self.node_df.loc['one'].columns:
            return self.node_df.loc['one', arg]
        if isinstance(arg, (str,)):
            key = self._expand_key(arg)
            return self.node_df.loc['one', key]
        if isinstance(arg, (list, tuple, ndarray)):
            key, labels = self._get_key_labels(arg)
            if all(in1d(labels, self.node_df.loc['one'].index)):
                return self.node_df.loc['one'].ix[labels, key]
        raise IndexError(arg)

    def __setitem__(self, arg, value):
        if arg in self.node_df.loc['two'].columns:
            self.node_df.loc['two',arg] = value
            return
        if isinstance(arg, (str,)):
            for (i, key) in enumerate(self._expand_key(arg)):
                self.node_df.loc['two',key] = value[:,i]
            return
        raise IndexError(arg)

    def advance(self):
        for key in self.node_df.loc['one'].index:
            self.node_df.ix['one',key] = self.node_df.loc['two',key]

def test():
    dimensions = 3
    node_labels = [1, 2, 5]
    node_vars = (('U', VECTOR), ('T', SCALAR, 195))
    elem_blocks = None
    dw = DataWharehouse(dimensions, node_labels, node_vars, elem_blocks)

    print('='*80)
    print(dw.node_df)
    print('='*80)
    print(dw['U_X'])
    print('='*80)
    u = dw['U']
    print(u)
    print('='*80)
    print(dw['U', (2,5)])
    dw['U'] = arange(u.size).reshape(u.shape)
    dw['T'] = 500
    print(dw.node_df)
    dw.advance()
    print(dw.node_df)



if __name__ == '__main__':
    test()
