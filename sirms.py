#-------------------------------------------------------------------------------
# Name:        sirms.py --in mols.sdf --out descriptors.txt
#              sirms.py -h for help
#
#
# Purpose:     calculate simplex descriptors for molecules in sdf file
#
# Author:      Pavel Polishchuk
#
# Created:     01.01.2013
# Copyright:   (c) Pavel Polishchuk 2013
# Licence:     GPLv3
# Python Version: 3.2
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import sys, os, argparse
from itertools import combinations, chain, product
from sdf import ReadSDF
from files import LoadRangedProperty, LoadMixturesTxt, LoadFragments
from ppgfunctions import *
from canon import LoadSirmsDict, GetCanonNameByDict, GenCanonName, GetSirmsType2, GetSirmsType
from mols import MergeMols3
from multiprocessing import Pool, cpu_count
from datetime import datetime
import time

#===============================================================================
# Save simplexes

def SaveSimplexes(fname, sirms):
    """
    Save calculated decriptors in the tab-delimited text file format
    Descriptors   compound_1  compound_2  ...  compound_n
    descriptor_1           5           2                0
    descriptor_2           1           1                4
    ...
    """
    f_out = open(fname, 'wt')
    # get sorted unique simplex names
    sirms_names = sorted(list(set(list(chain.from_iterable([list(s.keys()) for s in sirms.values()])))))
    compound_names = sorted([k for k in sirms.keys()])
    f_out.write("Descriptors\t" + "\t".join(compound_names) + "\n")
    for n in sirms_names:
        line = []
        for m in compound_names:
            line.append(sirms[m].get(n, 0))
        f_out.write(n + "\t" + "\t".join(map(str, line)) + "\n")
    f_out.close()

#===============================================================================
# SIRMS

def SetLabels(mols, opt_diff, input_fname):
    if 'type' in opt_diff:
        for m in mols.values():
            m.SetSybylTypes()
    if 'chg' in opt_diff: LoadRangedProperty(mols, os.path.join(os.path.abspath(os.path.dirname(input_fname))), GetFileNameNoExt(input_fname) + '.chg')
    if 'lip' in opt_diff: LoadRangedProperty(mols, os.path.join(os.path.abspath(os.path.dirname(input_fname))), GetFileNameNoExt(input_fname) + '.lip')
    return(None)

def CalcSingleCompSirms(mol, sirms_dict, diff_list, sirms_types, noH, verbose, frags=None):
    """
    INPUT:
        mol: Mol3 object;
        sirms_dict: loaded precomputed dictionary;
        diff_list: list of labeling to compute;
        sirms_types: list of topological simplex types;
    OUTPUT:
        returns dict of dicts
        if frags are specified then it returns descriptors for fragment-depleted molecule
        { "mol_name":               { "descriptors_1": 0, "descriptor_2": 5, ... },
          "mol_name#frag_name_1":   { "descriptors_1": 0, "descriptor_2": 3, ... },
          ...
        }

    """

    if verbose:
        cur_time = time.time()
    nodict = len(sirms_dict) == 0

    d = {mol.title: {}}
    if frags is not None:
        local_frags = frags.get(mol.title, None)
        if local_frags is not None:
            for frag in local_frags.keys():
                d[frag] = {}
    else:
        local_frags = None

    if noH:
        atoms = [a for a in mol.atoms.keys() if mol.atoms[a]["label"] != 'H']
    else:
        atoms = mol.atoms.keys()

    for a in combinations(atoms, 4):
        # if simplex is of non-allowed type it will not be calculated
        if GetSirmsType2(mol, a) not in sirms_types: continue
        # start simplex calculation
        bonds = [mol.GetBondType(b[0], b[1]) for b in combinations(a, 2)]
        for s_diff in diff_list:
            if s_diff == 'elm':
                labels = [mol.atoms[a_id]["label"] for a_id in a]
            else:
                labels = [mol.atoms[a_id]['property'][s_diff]['label'] for a_id in a]
            if nodict:
                canon_name = GenCanonName(labels, bonds, a)
            else:
                canon_name = GetCanonNameByDict(labels, bonds, sirms_dict)
            sirms_name = 'S|A|' + s_diff + '|' + canon_name
            d[mol.title][sirms_name] = d[mol.title].get(sirms_name, 0) + 1
            if local_frags is not None:
                for frag in local_frags.keys():
                    # if there is no common atoms in simplex and fragment
                    if set(a).isdisjoint(local_frags[frag]):
                        d[frag][sirms_name] = d[frag].get(sirms_name, 0) + 1

    # rename frags in output dict
    output = {}
    for n in d.keys():
        if n == mol.title:
            output[mol.title] = d[n]
        else:
            output[mol.title+"#"+n] = d[n]
    del(d)

    if verbose:
        print(mol.title, (str(round(time.time()-cur_time, 1)) + ' s').rjust(78-len(mol.title)))

    return(output)

def CalcBinMixSirms(mol_list, id_list, sirms_dict, diff_list, sirms_types, noH, mix_ordered, verbose):
    """
    Calculate simplex descriptors which belong only to both components of a binary mixture
    INPUT:
        mol_list: list of two MOL objects
        id_list: list of index number of a compunent in ordered mixtures,
            for unordered mixtures this list can be of any value, since it
            doesn't use further
    OUTPUT: dict of descriptors
    {"descriptor_1": 2,
     "descriptor_2": 6,
     "descriptor_3": 8,
     ...
    }
    """

    def GetMixBonds2(a1, a2):
        b = [0,0,0,0,0,0]
        if len(a2) == 3:
            b[3] = mol_list[1].GetBondType(a2[0], a2[1])
            b[4] = mol_list[1].GetBondType(a2[0], a2[2])
            b[5] = mol_list[1].GetBondType(a2[1], a2[2])
        elif len(a1) == 3:
            b[0] = mol_list[0].GetBondType(a1[0], a1[1])
            b[1] = mol_list[0].GetBondType(a1[0], a1[2])
            b[3] = mol_list[0].GetBondType(a1[1], a1[2])
        elif len(a1) == 2:
            b[0] = mol_list[0].GetBondType(a1[0], a1[1])
            b[5] = mol_list[1].GetBondType(a2[0], a2[1])
        return(b)

    if verbose:
        cur_time = time.time()
    d = {}
    nodict = len(sirms_dict) == 0

    if noH:
        atoms1 = [a for a in mol_list[0].atoms.keys() if mol_list[0].atoms[a]["label"] != 'H']
        atoms2 = [a for a in mol_list[1].atoms.keys() if mol_list[1].atoms[a]["label"] != 'H']
    else:
        atoms1 = mol_list[0].atoms.keys()
        atoms2 = mol_list[1].atoms.keys()

    for j in range(1,4):
        for a1, a2 in product(combinations(atoms1, j), combinations(atoms2, 4-j)):
            bonds = GetMixBonds2(a1, a2)
            # if simplex is of non-allowed type it will not be calculated
            if GetSirmsType(bonds) not in sirms_types: continue
            for s_diff in diff_list:
                # get labels
                if s_diff == 'elm':
                    if mix_ordered:
                        labels = [id_list[0] + mol_list[0].atoms[a11]["label"] for a11 in a1] + [id_list[1] + mol_list[1].atoms[a22]["label"] for a22 in a2]
                    else:
                        labels = [mol_list[0].atoms[a11]["label"] for a11 in a1] + [mol_list[1].atoms[a22]["label"] for a22 in a2]
                else:
                    if mix_ordered:
                        labels = [id_list[0] + mol_list[0].atoms[a11]['property'][s_diff]['label'] for a11 in a1] + [id_list[1] + mol_list[1].atoms[a22]['property'][s_diff]['label'] for a22 in a2]
                    else:
                        labels = [mol_list[0].atoms[a11]['property'][s_diff]['label'] for a11 in a1] + [mol_list[1].atoms[a22]['property'][s_diff]['label'] for a22 in a2]
                canon_name = GenCanonName(labels, bonds, [1,2,3,4]) if nodict else GetCanonNameByDict(labels, bonds, sirms_dict)
                sirms_name = 'M|A|' + s_diff + '|' + canon_name
                d[sirms_name] = d.get(sirms_name, 0) + 1
    if verbose:
        print(mol_list[0].title, mol_list[1].title, (str(round(time.time()-cur_time, 1)) + ' s').rjust(77-len(mol_list[0].title)-len(mol_list[1].title)))
    return(d)

def CalcMixSirms(mol_list, ratio_list, single_sirms, base_bin_mix_sirms, mix_ordered):
    """
    Calculate simplexes for mixture of compounds based on pre-calculated
    simplexes for binary compositions and individual components
    INPUT:
        mol_list - list of molecules which are components of a mixture in order
                   given in a dataset text file
        ratio_list - list of ratio of each component given in the same order
    OUTPUT:
        dict of weighted descriptors for a mixture accordin to given ratio of components
    """
    d = {}
    # descriptors from separate components
    if mix_ordered:
        # rename descriptors (add component ids to labels) and calc weighted values
        for i, mol in enumerate(mol_list):
            for s_name, value in single_sirms[mol.title].items():
                tmp = s_name.split('|')
                tmp[3] = ','.join([str(i+1) + a for a in tmp[3].split(',')])
                new_name = '|'.join(tmp)
                d[new_name] = value * ratio_list[i]
    else:
        # combine descriptors of each pair of components and calc weighted values
        m_titles = [m.title for m in mol_list]
        for i1, i2 in combinations(range(len(mol_list)), 2):
            for s_name in set(list(single_sirms[m_titles[i1]].keys()) + list(single_sirms[m_titles[i2]].keys())):
                d[s_name] = single_sirms[m_titles[i1]].get(s_name, 0) * ratio_list[i1] + single_sirms[m_titles[i2]].get(s_name, 0) * ratio_list[i2]
    # descriptors from mixtures
    m_titles = [str(i+1) + '#' + m.title for i, m in enumerate(mol_list)] if mix_ordered else [m.title for m in mol_list]
    for i1, i2 in combinations(range(len(mol_list)), 2):
        min_ratio = min(ratio_list[i1], ratio_list[i2])
        mix_name = (m_titles[i1], m_titles[i2]) if m_titles[i1] < m_titles[i2] else (m_titles[i2], m_titles[i1])
        for s_name in base_bin_mix_sirms[mix_name].keys():
            d[s_name] = d.get(s_name, 0) + min_ratio * base_bin_mix_sirms[mix_name][s_name]
    return(d)

def GetBaseBinMixSirmsMP(uniq_bin_mix, mols, sirms_dict, sirms_diff, sirms_types, noH, mix_ordered, verbose):
    """
    return dictionary of descriptors of all binary mixtures in a dataset,
    it doesn't take into account the raio of components
    OUTPUT:
    {
      ("1#mol_name_1", "2#mol_name_2"): {"descriptor_1": 3, "descriptor_2": 6, ... },
      ...
    }

    """
    if verbose:
        print(' Basic binary mixtures calculation '.center(79, '-'))
    mix_set = sorted(uniq_bin_mix)
    p = Pool(processes = cpu_count())
    if not mix_ordered:
        res = [p.apply_async(CalcBinMixSirms, [[mols[n1], mols[n2]], [0, 1], sirms_dict, sirms_diff, sirms_types, noH, mix_ordered, verbose]) for n1, n2 in mix_set]
    else:
        res = []
        for n1, n2 in mix_set:
            id1, name1 = n1.split('#', 1)
            id2, name2 = n2.split('#', 1)
            res.append(p.apply_async(CalcBinMixSirms, [[mols[name1], mols[name2]], [id1, id2], sirms_dict, sirms_diff, sirms_types, noH, mix_ordered, verbose]))
    d = {mix_set[i]: r.get() for i, r in enumerate(res)}
    p.close()
    return(d)

def CalcSingleCompSirmsMP(mol_list, sirms_dict, opt_diff, opt_types, opt_noH, ncores, opt_verbose, frags=None):
    """
    Multi-processed calculation of descriptors of single compounds depending on ncores option
    INPUT:
        list of molecules and additional parameters for calculation of simplex descriptors
    OUTPUT:
        dict of compounds with calculated descriptors as a dict
        { "compound_1":
            { "descriptor_1: 2,
              "descriptor_2: 5,
              ...
            }
          ...
        }
    NOTE: To speed up calculation for small compounds or with enabled option noH
    it is better to choose single-core calculation (ncores = 1), since there is
    an overhead in process creation and management.

    """
    if opt_verbose:
        print(' Single compounds calculation '.center(79, '-'))
    d = {}
    if ncores > 1:
        p = Pool(processes = min(len(mol_list), ncores))
        res = [p.apply_async(CalcSingleCompSirms, [mol, sirms_dict, opt_diff, opt_types, opt_noH, opt_verbose, frags]) for mol in mol_list]
        for r in res:
            d.update(r.get())
        p.close()
    else:
        for mol in mol_list:
            d.update(CalcSingleCompSirms(mol, sirms_dict, opt_diff, opt_types, opt_noH, opt_verbose, frags))
    return(d)

#===============================================================================

def GetUniqBinMixNames(mix, ordered):
    """
    INPUT:
        mix - dict loaded from the text file mixtures.txt
    OUTPUT:
        set of tuples with components names of all binary mixtures
        for unordered mixtures:
            (("mol_name_1", "mol_name_2"), ("mol_name_1", "mol_name_3"), ...)
        for ordered mixtures:
            (("1#mol_name_1", "2#mol_name_2"), ("1#mol_name_1", "2#mol_name_3"), ...)
    """
    d = set()
    for m in mix.values():
        names = m['names']
        if ordered:
            names = [str(i+1) + '#' + n for i, n in enumerate(names)]
        for c1, c2 in combinations(names, 2):
            if c1 > c2: c1, c2 = c2, c1
            d.add((c1, c2))
    return(d)

#===============================================================================
# Main cycle

def main():
    parser = argparse.ArgumentParser(description='Calculate simplex descriptors for molecules in sdf-file.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
           help='input sdf file with standartized structures, molecules should have titles')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
           help='output txt file with calculated descriptors')
    parser.add_argument('-n', '--nodict', action='store_true', default=False,
           help='if set this flag the simplexes will be generated slower but this procedure can handle any bond types, while the other approach (which uses dictionary) can handle structures containing only 0-4 bond types')
    parser.add_argument('-d', '--diff', metavar='', default=['elm'], nargs='*',
           choices=['elm', 'type', 'lip', 'chg'],
           help='list of atom labeling schemes separated by space, allowed values: type, elm, chg, lip. Charges and lipophilicity should be calculated in HiTQSAR and saved with the same names as sdf-file and extensions .chg and .lip correspondingly along with sdf-file. In the same folder the "setup.txt" file should be stored. Default value = elm')
    parser.add_argument('-t', '--types', metavar='', default='all',
           choices=['all', 'bounded', 'extended'],
           help='list of simplex types which should be calculated. There three possible values: all, bounded=5,6,8-11, extended=3-11. Default value = all')
    parser.add_argument('-m', '--mixtures', metavar='mixtures.txt',
           help='text file containing list of mixtures of components and their ratios. Names of components should be the same as in input.sdf file. The header should contain the string "!absolute ratio" or "!relative ratio".')
    parser.add_argument('-r', '--mix_ordered', action='store_true', default=False,
           help='if set this flag the mixtures will be considered ordered, otherwise as unordered')
    parser.add_argument('-c', '--ncores', metavar='[all, 1, 2, ...]', default='1',
           help='number of cpu cores which will be used for calculation of sirms for single compounds. All cores are always used for calculation of mixtures. Default = all. Hint: for small single compounds the best choice is 1 for highest calculation speed')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
           help='if set this flag progress will be printed out (may cause decrease in speed).')
    parser.add_argument('-x', '--noH', action='store_true', default=False,
           help='if set this flag hydrogen atoms will be excluded from the simplexes calculation')
    parser.add_argument('-f', '--fragments', metavar='fragments.txt',
           help='text file containing list of names of single compounds, fragment names and atom indexes of fragment in the structure of corresponding compound')

    args = vars(parser.parse_args())
    opt_mix_ordered = None
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "nodict":
            if v: sirms_dict = {}
            else: sirms_dict = LoadSirmsDict()
        if o == "diff": opt_diff = v
        if o == "types":
            if v == "all": opt_types = list(range(1,12))
            if v == "bounded": opt_types = [5,6,8,9,10,11]
            if v == "extended": opt_types = list(range(3,12))
        if o == "mixtures": mix_fname = v
        if o == "mix_ordered": opt_mix_ordered = v
        if o == "ncores": opt_ncores = cpu_count() if v == "all" else int(v)
        if o == "verbose": opt_verbose = v
        if o == "noH": opt_noH = v
        if o == "fragments": frag_fname = v

    mols = ReadSDF(in_fname)
    SetLabels(mols, opt_diff, in_fname)
    if mix_fname is None:
        frags = LoadFragments(frag_fname)
        # calc simplex descriptors
        sirms = CalcSingleCompSirmsMP([m for m in mols.values()], sirms_dict, opt_diff, opt_types, opt_noH, opt_ncores, opt_verbose, frags)
        SaveSimplexes(out_fname, sirms)
    else:
        mix = LoadMixturesTxt(mix_fname)
        mols_used = set(chain.from_iterable([m['names'] for m in mix.values()]))
        base_single_sirms = CalcSingleCompSirmsMP([m for m in mols.values() if m.title in mols_used], sirms_dict, opt_diff, opt_types, opt_noH, opt_ncores, opt_verbose)
        # calc basic sirms for all binary mixtures (no weights, all mixtures are 1:1)
        uniq_bin_mix = GetUniqBinMixNames(mix, opt_mix_ordered)
        base_bin_mix_sirms = GetBaseBinMixSirmsMP(uniq_bin_mix, mols, sirms_dict, opt_diff, opt_types, opt_noH, opt_mix_ordered, opt_verbose)
        mix_sirms = {}
        for m_name, m in mix.items():
            mix_sirms[m_name] = CalcMixSirms([mols[mol_name] for mol_name in m['names']], m['ratios'], base_single_sirms, base_bin_mix_sirms, opt_mix_ordered)
        SaveSimplexes(out_fname, mix_sirms)

if __name__ == '__main__':
    main()
