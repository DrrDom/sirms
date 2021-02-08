#-------------------------------------------------------------------------------
# Name:        files
# Purpose:     working with general purpose files
#
# Author:      Pavel Polishchuk
#
# Created:     06.06.2013
# Copyright:   (c) Pavel Polishchuk 2013-2015
# Licence:     BSD 3-clause
#-------------------------------------------------------------------------------

import os
from string import ascii_uppercase
from .ppgfunctions import GetFileNameNoExt, GetWorkDir
from collections import OrderedDict


class SvmSaver:

    def __init__(self, file_name):
        self.__fname = file_name
        self.__varnames_fname = os.path.splitext(file_name)[0] + '.colnames'
        self.__molnames_fname = os.path.splitext(file_name)[0] + '.rownames'
        self.__varnames = []
        # self.__molnames = []
        if os.path.isfile(self.__fname):
            os.remove(self.__fname)
        if os.path.isfile(self.__molnames_fname):
            os.remove(self.__molnames_fname)
        if os.path.isfile(self.__varnames_fname):
            os.remove(self.__varnames_fname)

    def __convert_value(self, value):
        if type(value) is int or value.is_integer():
            res = str(int(value))
        elif isinstance(value, float):
            res = "{:.5f}".format(value)
        else:
            res = str(value)
        return res

    def save_mol_descriptors(self, mol_name, mol_descr_dict):

        values = []

        for i, varname in enumerate(self.__varnames):
            if varname in mol_descr_dict:
                value = self.__convert_value(mol_descr_dict[varname])
                values.append((i, value))

        new_varnames = list(set(mol_descr_dict) - set(self.__varnames))
        for i, varname in enumerate(new_varnames):
            value = self.__convert_value(mol_descr_dict[varname])
            values.append((len(self.__varnames) + i, value))

        if values:  # values can be empty if all descriptors are zero

            self.__varnames.extend(new_varnames)

            with open(self.__molnames_fname, 'at') as f:
                f.write(mol_name + '\n')

            if new_varnames:
                with open(self.__varnames_fname, 'at') as f:
                    f.write('\n'.join(new_varnames) + '\n')

            with open(self.__fname, 'at') as f:
                values = sorted(values)
                values = ('%i:%s' % (i, v) for i, v in values)
                f.write(' '.join(values) + '\n')


def GetAtomPropertyFromSetup(setup_fname):
    with open(setup_fname) as f:
        output = [line.strip().split('=')[0] for line in f]
    return output


# def NotExistedPropertyFiles(opt_diff, input_fname):
#     """
#     Checks property files existence and returns names of not existed property files,
#     these properties will be loaded from the initial sdf file
#     """
#     output = []
#     for s_diff in opt_diff:
#         if s_diff not in ["elm", "none", "uffd", "uffe"]:   # built-in types are ignored
#             if not os.path.isfile(os.path.join(GetWorkDir(input_fname),
#                                                GetFileNameNoExt(input_fname) + '.' + s_diff)):
#                 output.append(s_diff)
#     return (output)


def ReadPropertyRange(file_setup_name, prop_name):
    res = None
    try:
        f = open(file_setup_name)
        for line in f:
            if prop_name == line.split('=')[0]:
                res = list(map(float, line.strip().split('=')[1].split('<')))
                break
    except FileNotFoundError:
        print('File setup.txt cannot be found in the folder containing input sdf file.')
        exit()
    else:
        f.close()
    return res


def RangedLetter(value, prop_range, none_label="NA"):
    if value is None:
        return (none_label)
    if value <= prop_range[0]:
        return ('A')
    for i in range(1, len(prop_range)):
        if prop_range[i - 1] < value <= prop_range[i]:
            return (ascii_uppercase[i])
        if i == len(prop_range) - 1 and value > prop_range[-1]:
            return (ascii_uppercase[i + 1])
    print('Property label was not assigned.')
    return (None)


def LoadRangedProperty(mols, setup_dir, prop_fname):

    def ReadPropertyFile(fname):
        f = open(fname)
        sep = '---'
        d = {}
        f.readline()  # header
        # read each mol as list of lists where inner list is ['atom_id', 'value']
        while True:
            line = f.readline()
            if not line: break
            if line.strip() == sep:
                title = f.readline().strip()
                d[title] = list()
                continue
            d[title].append(line.strip().split())
        f.close()
        # sort values for each mol by atom_id and keep values only
        for k in d.keys():
            d[k].sort(key=lambda x: int(x[0]))
            # if value is not numeric remain it as is (as string)
            try:
                d[k] = [float(i[1]) for i in d[k]]
            except ValueError:
                d[k] = [i[1] for i in d[k]]
        return (d)

    prop_name = os.path.splitext(os.path.basename(prop_fname))[1][1:]
    prop_range = ReadPropertyRange(os.path.join(setup_dir, 'setup.txt'), prop_name)
    # read property file
    d = ReadPropertyFile(os.path.join(setup_dir, prop_fname))
    # assign property values to atoms
    for m in mols.values():
        values = d[m.title]
        for i, a in enumerate(sorted(m.atoms.keys())):
            # if label is string use it as label, otherwise use ranges
            label = values[i] if type(values[i]) is str else RangedLetter(values[i], prop_range)
            m.atoms[a]['property'][prop_name] = {'value': values[i], 'label': label}


def LoadMixturesTxt(fname, mix_type):
    """
    Load mixture composition
    INPUT: tab-delimited text file
    mol_name_1  mol_name_2  first_component_ratio  second_component_ratio
    mol_name_1  mol_name_3  first_component_ratio  second_component_ratio
    ...
    NOTE: number of components can be any
    OUTPUT:
    { "mol_name_1_(ratio_1)+mol_name_2_(ratio_2)":
        {'names': ["mol_name_1", "mol_name_2"], 'ratios': [ratio_1, ratio_2]},
      "mol_name_1_(ratio_1)+mol_name_3_(ratio_2)":
        {'names': ["mol_name_1", "mol_name_3"], 'ratios': [ratio_1, ratio_2]},
      ...
    }
    OUTPUT EXAMPLE:
    { 'P-Ser-Ser_(1.0)+E-37_(1.0)': {'names': ['P-Ser-Ser', 'E-37'], 'ratios': [1.0, 1.0]},
      'K-Arg-Ser_(1.0)+MS-7_(1.0)': {'names': ['K-Arg-Ser', 'MS-7'], 'ratios': [1.0, 1.0]},
      'K-Ala-Glu_(1.0)+E-21_(1.0)': {'names': ['K-Ala-Glu', 'E-21'], 'ratios': [1.0, 1.0]} }
    """
    d = OrderedDict()

    if mix_type == 'abs':
        abs_ratio = True
    else:
        abs_ratio = False

    f = open(fname, 'rt')

    while True:
        line = f.readline()
        if not line: break
        tmp = line.strip().split('\t')
        mix_num = len(tmp) // 2
        mol_names = tmp[:mix_num]
        mol_ratios = list(map(float, tmp[(len(tmp) // 2):]))
        if not abs_ratio:
            mol_ratios = [v / sum(mol_ratios) for v in mol_ratios]
        mol_ratios = [round(v, 4) for v in mol_ratios]
        mix_name = '+'.join([mol_names[i] + '_(' + str(mol_ratios[i]) + ')' for i in range(mix_num)])
        d[mix_name] = {'names': mol_names, 'ratios': mol_ratios}
    return d


def LoadFragments(fname):
    """
    INPUT
    Input file format - tab delimited text file, each line contains
    mol_name  fragment_name  atom_index_1  atom_index_2 ... atom_index_n
    where atom ids of fragment atoms in the molecule are listed
    OUTPUT
    Output dict of molecules which are dict of fragments with corresponding atom ids
    { "mol_name_1":
        {
          "frag_name_1#0": [atom_id1, atom_id2, ...],
          "frag_name_1#1": [atom_id1, atom_id2, ...],
          "frag_name_2#0": [atom_id1, atom_id2, ...]
        }
      ...
    }
    """
    if fname is None:
        return (None)
    d = OrderedDict()
    f = open(fname, 'rt')
    for line in f:
        tmp = line.strip().split("\t")
        if tmp[0] not in d.keys():
            d[tmp[0]] = OrderedDict()
        d[tmp[0]][tmp[1] + "#" + str(len(d[tmp[0]]))] = list(map(int, tmp[2:]))
    f.close()
    return (d)