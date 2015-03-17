#-------------------------------------------------------------------------------
# Name:        files
# Purpose:     working with general purpose files
#
# Author:      Pavel Polishchuk
#
# Created:     06.06.2013
# Copyright:   (c) Pavel Polishchuk 2013-2015
# Licence:     GPLv3
#-------------------------------------------------------------------------------

import os
from string import ascii_uppercase
from ppgfunctions import GetFileNameNoExt, GetWorkDir


def NotExistedPropertyFiles(opt_diff, input_fname):
    """
    Checks property files existence and returns names of not existed property files,
    these properties will be loaded from the initial sdf file
    """
    output = []
    for s_diff in opt_diff:
        if s_diff != "elm":
            if not os.path.isfile(os.path.join(GetWorkDir(input_fname),
                                               GetFileNameNoExt(input_fname) + '.' + s_diff)):
                output.append(s_diff)
    return (output)


def ReadPropertyRange(file_setup_name, prop_name):
    try:
        f = open(file_setup_name)
    except:
        print('File setup.txt cannot be found in the folder containing input sdf file.')
        exit()
    res = None
    for line in f:
        if prop_name == line.split('=')[0]:
            res = list(map(float, line.strip().split('=')[1].split('<')))
            break
    f.close()
    return (res)


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
    if prop_range is None:
        print('Cannot find "%s" property in setup.txt' % prop_name)
        exit()
    # read property file
    d = ReadPropertyFile(os.path.join(setup_dir, prop_fname))
    # assign property values to atoms
    for m in mols.values():
        values = d[m.title]
        for i, a in enumerate(sorted(m.atoms.keys())):
            # if label is string use it as label, otherwise use ranges
            label = values[i] if type(values[i]) is str else RangedLetter(values[i], prop_range)
            m.atoms[a]['property'][prop_name] = {'value': values[i], 'label': label}


def LoadMixturesTxt(fname):
    """
    Load mixture composition
    INPUT: tab-delimited textfile with first line
    !absolute ratio or !relative ratio (which desognates how to process ratio of components)
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
    d = {}
    f = open(fname, 'rt')
    header = f.readline()
    if header.strip() == '!absolute ratio':
        abs_ratio = True
    elif header.strip() == '!relative ratio':
        abs_ratio = False
    else:
        raise Exception(
            'The first line of mixtures.txt file should contain the header "!absolute ratio" or "!relative ratio"')
        sys.exit()
    while True:
        line = f.readline()
        if not line: break
        tmp = line.strip().split('\t')
        mix_num = len(tmp) // 2
        mol_names = tmp[:mix_num]
        mol_ratios = list(map(float, tmp[(len(tmp) // 2):]))
        if not abs_ratio:
            mol_ratios = [v / sum(mol_ratios) for v in mol_ratios]
        mix_name = '+'.join([mol_names[i] + '_(' + str(mol_ratios[i]) + ')' for i in range(mix_num)])
        d[mix_name] = {'names': mol_names, 'ratios': mol_ratios}
    return (d)


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
    d = {}
    f = open(fname, 'rt')
    for line in f:
        tmp = line.strip().split("\t")
        if tmp[0] not in d.keys():
            d[tmp[0]] = {}
        d[tmp[0]][tmp[1] + "#" + str(len(d[tmp[0]]))] = list(map(int, tmp[2:]))
    f.close()
    return (d)