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
# Copyright:   (c) Pavel Polishchuk 2013-2015
# Licence:     GPLv3
# Python Version: 3.2
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import argparse
import time
import files
from itertools import combinations, chain, product
from collections import OrderedDict

from sdf import ReadSDF, ReadRDF, ReadRXN
from ppgfunctions import *
from canon import LoadSirmsDict, GetCanonNameByDict, GenCanonName, GetSirmsType2, GetSirmsType
from multiprocessing import Pool, cpu_count


#===============================================================================
# Save simplexes

def SaveSimplexes(fname, sirms, ndigits=5):
    """
    Save calculated decriptors in the tab-delimited text file format
    Descriptors   descriptor_1  descriptor_2  ...  descriptor_n
    compound_1               5             2                  0
    compound_2               1             1                  4
    ...
    """
    f_out = open(fname, 'wt')
    # get sorted unique simplex names
    sirms_names = sorted(list(set(list(chain.from_iterable([list(s.keys()) for s in sirms.values()])))))
    # compound_names = sorted([k for k in sirms.keys()])
    f_out.write("Compounds\t" + "\t".join(sirms_names) + "\n")
    s = "{:." + str(ndigits) + "f}"
    for n in sirms.keys():
        line = []
        for m in sirms_names:
            value = sirms[n].get(m, 0)
            # round floating point values, all others will be remained unchanged
            if type(value) is int or value.is_integer():
                line.append(str(int(value)))
            elif isinstance(value, float):
                line.append(s.format(value))
            else:
                line.append(str(value))
        f_out.write(n + "\t" + "\t".join(map(str, line)) + "\n")
    f_out.close()


#===============================================================================
# create reaction sirms
def concat_reaction_sirms(sirms):

    output = OrderedDict()

    while sirms:

        # take first element
        k, v = sirms.popitem(last=False)
        rx_id, role = k.rsplit('_', 1)
        # first letter from role is used as additional description
        new_sirms_names = [role[0] + '|' + s for s in sorted(v.keys())]
        s1 = {new_sirms_names[i]: v[sname] for i, sname in enumerate(sorted(v.keys()))}

        # take complimentary element (reactants or products of the same reaction)
        if role == 'reactants':
            k = rx_id + '_products'
            new_sirms_names = ['p|' + s for s in sorted(sirms[k].keys())]
        elif role == 'products':
            k = rx_id + '_reactants'
            new_sirms_names = ['r|' + s for s in sorted(sirms[k].keys())]
        else:
            print('Impossible error')
        v = sirms[k]
        s2 = {new_sirms_names[i]: v[sname] for i, sname in enumerate(sorted(sirms[k].keys()))}
        del(sirms[k])

        # combine them
        s1.update(s2)
        output[rx_id] = s1

    return output


#===============================================================================
# SIRMS

def SetLabels(mols, opt_diff, input_fname):
    for s_diff in opt_diff:
        if s_diff not in ['elm', 'none']:
            files.LoadRangedProperty(mols, GetWorkDir(input_fname), GetFileNameNoExt(input_fname) + '.' + s_diff)
    return (None)


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
        if GetSirmsType2(mol, a) not in sirms_types:
            continue
        # start simplex calculation
        bonds = [mol.GetBondType(b[0], b[1]) for b in combinations(a, 2)]
        for s_diff in diff_list:
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
            output[mol.title + "#" + n] = d[n]
    del (d)

    if verbose:
        print(mol.title, (str(round(time.time() - cur_time, 1)) + ' s').rjust(78 - len(mol.title)))

    return (output)


def CalcBinMixSirms(mol_list, id_list, sirms_dict, diff_list, sirms_types, noH, mix_ordered, verbose):
    """
    Calculate simplex descriptors which belong only to both components of a binary mixture
    INPUT:
        mol_list: list of two MOL objects
        id_list: list of index number of a component in ordered mixtures,
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
        b = [0, 0, 0, 0, 0, 0]
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
        return (b)

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

    for j in range(1, 4):
        for a1, a2 in product(combinations(atoms1, j), combinations(atoms2, 4 - j)):
            bonds = GetMixBonds2(a1, a2)
            # if simplex is of non-allowed type it will not be calculated
            if GetSirmsType(bonds) not in sirms_types: continue
            for s_diff in diff_list:
                # get labels
                if mix_ordered:
                    labels = [id_list[0] + mol_list[0].atoms[a11]['property'][s_diff]['label'] for a11 in a1] + [
                        id_list[1] + mol_list[1].atoms[a22]['property'][s_diff]['label'] for a22 in a2]
                else:
                    labels = [mol_list[0].atoms[a11]['property'][s_diff]['label'] for a11 in a1] + [
                        mol_list[1].atoms[a22]['property'][s_diff]['label'] for a22 in a2]

                canon_name = GenCanonName(labels, bonds, [1, 2, 3, 4]) if nodict else GetCanonNameByDict(labels, bonds,
                                                                                                         sirms_dict)
                sirms_name = 'M|A|' + s_diff + '|' + canon_name
                d[sirms_name] = d.get(sirms_name, 0) + 1
    if verbose:
        print(mol_list[0].title, mol_list[1].title, (str(round(time.time() - cur_time, 1)) + ' s').rjust(
            77 - len(mol_list[0].title) - len(mol_list[1].title)))
    return (d)


def CalcMixSirms(mol_list, ratio_list, single_sirms, base_bin_mix_sirms, mix_ordered):
    """
    Calculate simplexes for mixture of compounds based on pre-calculated
    simplexes for binary compositions and individual components
    INPUT:
        mol_list - list of molecules which are components of a mixture in order
                   given in a dataset text file
        ratio_list - list of ratio of each component given in the same order
    OUTPUT:
        dict of weighted descriptors for a mixture according to given ratio of components
    """
    d = {}

    # descriptors from separate components
    if mix_ordered:
        # rename descriptors (add component ids to labels) and calc weighted values
        for i, mol in enumerate(mol_list):
            for s_name, value in single_sirms[mol.title].items():
                tmp = s_name.split('|')
                tmp[3] = ','.join([str(i + 1) + a for a in tmp[3].split(',')])
                new_name = '|'.join(tmp)
                d[new_name] = value * ratio_list[i]
    else:
        # combine descriptors of each pair of components and calc weighted values
        m_titles = [m.title for m in mol_list]
        if len(mol_list) == 1:
            for s_name, s_value in single_sirms[m_titles[0]].items():
                d[s_name] = s_value * ratio_list[0]
        else:
            for i1, i2 in combinations(range(len(mol_list)), 2):
                for s_name in set(list(single_sirms[m_titles[i1]].keys()) + list(single_sirms[m_titles[i2]].keys())):
                    d[s_name] = single_sirms[m_titles[i1]].get(s_name, 0) * ratio_list[i1] + single_sirms[m_titles[i2]].get(
                        s_name, 0) * ratio_list[i2]

    # descriptors from mixtures
    m_titles = [str(i + 1) + '#' + m.title for i, m in enumerate(mol_list)] if mix_ordered else [m.title for m in
                                                                                                 mol_list]
    for i1, i2 in combinations(range(len(mol_list)), 2):
        min_ratio = min(ratio_list[i1], ratio_list[i2])
        mix_name = (m_titles[i1], m_titles[i2]) if m_titles[i1] < m_titles[i2] else (m_titles[i2], m_titles[i1])
        for s_name in base_bin_mix_sirms[mix_name].keys():
            d[s_name] = d.get(s_name, 0) + min_ratio * base_bin_mix_sirms[mix_name][s_name]
    return (d)


def GetBaseBinMixSirmsMP(uniq_bin_mix, mols, sirms_dict, sirms_diff, sirms_types, noH, mix_ordered, ncores, verbose):
    """
    return dictionary of descriptors of all binary mixtures in a dataset,
    it doesn't take into account the ratio of components
    OUTPUT:
    {
      ("1#mol_name_1", "2#mol_name_2"): {"descriptor_1": 3, "descriptor_2": 6, ... },
      ...
    }

    """
    if verbose:
        print(' Basic binary mixtures calculation '.center(79, '-'))
    mix_set = sorted(uniq_bin_mix)
    res = []
    if ncores == -1:
        if not mix_ordered:
            res = [CalcBinMixSirms([mols[n1], mols[n2]], [0, 1], sirms_dict, sirms_diff, sirms_types, noH, mix_ordered,
                                   verbose) for n1, n2 in mix_set]
        else:
            for n1, n2 in mix_set:
                id1, name1 = n1.split('#', 1)
                id2, name2 = n2.split('#', 1)
                res.append(
                    CalcBinMixSirms([mols[name1], mols[name2]], [id1, id2], sirms_dict, sirms_diff, sirms_types, noH,
                                    mix_ordered, verbose))
        d = {mix_set[i]: r for i, r in enumerate(res)}
    else:
        p = Pool(processes=cpu_count()) if ncores > 0 else Pool(processes=abs(ncores))
        if not mix_ordered:
            res = [p.apply_async(CalcBinMixSirms,
                                 [[mols[n1], mols[n2]], [0, 1], sirms_dict, sirms_diff, sirms_types, noH, mix_ordered,
                                  verbose]) for n1, n2 in mix_set]
        else:
            for n1, n2 in mix_set:
                id1, name1 = n1.split('#', 1)
                id2, name2 = n2.split('#', 1)
                res.append(p.apply_async(CalcBinMixSirms,
                                         [[mols[name1], mols[name2]], [id1, id2], sirms_dict, sirms_diff, sirms_types,
                                          noH, mix_ordered, verbose]))
        d = {mix_set[i]: r.get() for i, r in enumerate(res)}
        p.close()
    return (d)


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
    d = OrderedDict()
    if abs(ncores) > 1:
        p = Pool(processes=min(len(mol_list), abs(ncores)))
        res = [p.apply_async(CalcSingleCompSirms, [mol, sirms_dict, opt_diff, opt_types, opt_noH, opt_verbose, frags])
               for mol in mol_list]
        for r in res:
            d.update(r.get())
        p.close()
    else:
        for mol in mol_list:
            d.update(CalcSingleCompSirms(mol, sirms_dict, opt_diff, opt_types, opt_noH, opt_verbose, frags))
    return (d)


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
            names = [str(i + 1) + '#' + n for i, n in enumerate(names)]
        for c1, c2 in combinations(names, 2):
            if c1 > c2: c1, c2 = c2, c1
            d.add((c1, c2))
    return (d)


#===============================================================================

def GenQuasiMix(mol_names):
    return {n: {'names': [n, n], 'ratios': [0.5, 0.5]} for n in mol_names}


#===============================================================================
# Main cycle

def main_params(in_fname, out_fname, opt_no_dict, opt_diff, opt_types, mix_fname, opt_mix_ordered, opt_ncores,
                opt_verbose, opt_noH, frag_fname, parse_stereo, quasimix, id_field_name):

    if opt_no_dict:
        sirms_dict = {}
    else:
        sirms_dict = LoadSirmsDict()

    # define which property will be loaded from external file or from sdf-file
    opt_diff_sdf = files.NotExistedPropertyFiles(opt_diff, in_fname)
    opt_diff_ext = [el for el in opt_diff if el not in opt_diff_sdf]

    # load sdf, rdf or rxn file depending on its extension
    input_file_extension = in_fname.strip().split(".")[-1].lower()
    setup_path = os.path.join(GetWorkDir(in_fname), "setup.txt")
    if input_file_extension == 'sdf':
        mols = ReadSDF(in_fname, id_field_name, opt_diff_sdf, setup_path, parse_stereo)
    elif input_file_extension == 'rdf':
        mols, mix = ReadRDF(in_fname, id_field_name)
    elif input_file_extension == 'rxn':
        mols, mix = ReadRXN(in_fname, id_field_name)
    else:
        print("Input file extension should be SDF, RDF or RXN. Current file has %s. Please check it." % input_file_extension.upper())
        return None

    # set property labels on atoms from external data files
    SetLabels(mols, opt_diff_ext, in_fname)

    if mix_fname is None and not quasimix and input_file_extension == 'sdf':

        frags = files.LoadFragments(frag_fname)
        # calc simplex descriptors
        sirms = CalcSingleCompSirmsMP([m for m in mols.values()], sirms_dict, opt_diff, opt_types, opt_noH, opt_ncores,
                                      opt_verbose, frags)
        SaveSimplexes(out_fname, sirms)

    else:

        # create mix data for sdf (for rdf/rxn mix is created during file loading)
        if input_file_extension == 'sdf':
            if quasimix:
                mix = GenQuasiMix(list(mols.keys()))
            else:
                mix = files.LoadMixturesTxt(mix_fname)
        mols_used = set(chain.from_iterable([m['names'] for m in mix.values()]))
        base_single_sirms = CalcSingleCompSirmsMP([m for m in mols.values() if m.title in mols_used], sirms_dict,
                                                  opt_diff, opt_types, opt_noH, opt_ncores, opt_verbose)
        # calc basic sirms for all binary mixtures (no weights, all mixtures are 1:1)
        uniq_bin_mix = GetUniqBinMixNames(mix, opt_mix_ordered)
        base_bin_mix_sirms = GetBaseBinMixSirmsMP(uniq_bin_mix, mols, sirms_dict, opt_diff, opt_types, opt_noH,
                                                  opt_mix_ordered, opt_ncores, opt_verbose)
        mix_sirms = OrderedDict()
        for m_name, m in mix.items():
            mix_sirms[m_name] = CalcMixSirms([mols[mol_name] for mol_name in m['names']], m['ratios'],
                                             base_single_sirms, base_bin_mix_sirms, opt_mix_ordered)

        if input_file_extension in ['rdf', 'rxn']:
            mix_sirms = concat_reaction_sirms(mix_sirms)

        SaveSimplexes(out_fname, mix_sirms)


def main():

    parser = argparse.ArgumentParser(description='Calculate simplex descriptors for single molecules, quasi-mixtures, '
                                                 'mixtures and reactions.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
                        help='input file (allowed formats: sdf, rdf, rxn) with standardized structures, '
                             'molecules or reactions should have titles')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='output txt file with calculated descriptors')
    parser.add_argument('-n', '--nodict', action='store_true', default=False,
                        help='if set this flag the simplexes will be generated slower but this procedure can handle '
                             'any bond types, while the other approach (which uses dictionary) can handle structures '
                             'containing only 0-4 bond types')
    parser.add_argument('-d', '--diff', metavar='', default=['elm'], nargs='*',
                        help='list of atom labeling schemes separated by space. Built-in scheme is element (elm) and '
                             'topology (none). '
                             'To include other schemes user should specify the name of the corresponding property '
                             'value identical to the name of SDF field, which contains calculated atomic properties. '
                             'Fields names are case-insensitive and converted to lowercase. For RDF/RXN input files '
                             'only built-in types can be used. '
                             'Default value = elm')
    parser.add_argument('-t', '--types', metavar='', default='extended',
                        choices=['all', 'bounded', 'extended'],
                        help='list of simplex types which should be calculated. There three possible values: all, '
                             'bounded=5,6,8-11, extended=3-11. Default value = extended')
    parser.add_argument('-s', '--stereo', action='store_true', default=False,
                        help='parse stereo information from the <stereoanalysis> field generated by Chemaxon. '
                             'Works only for double bonds in sdf.')
    parser.add_argument('-m', '--mixtures', metavar='mixtures.txt', default=None,
                        help='text file containing list of mixtures of components and their ratios. Names of components'
                             ' should be the same as in input.sdf file. The header should contain the string '
                             '"!absolute ratio" or "!relative ratio". Works only with sdf files.')
    parser.add_argument('-r', '--mix_ordered', action='store_true', default=False,
                        help='if set this flag the mixtures will be considered ordered, otherwise as unordered. '
                             'Used only in combination with -m key.')
    parser.add_argument('-q', '--quasi_mix', action='store_true', default=False,
                        help='calculate quasi-mixture descriptors for pure compounds. Works only with sdf files.')
    parser.add_argument('-c', '--ncores', metavar='[all, 1, 2, ..., -1, -2, ...]', default='1',
                        help='negative number defines number of cores which will be available for calculation of '
                             'single compounds and mixtures. Positive number defines number of cores which will be '
                             'available for calculation of single compounds only; all cores will be used for '
                             'calculation of mixtures. Default = 1. Hint: for small single compounds and mixtures '
                             'the best choice is -1 for highest calculation speed')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='if set this flag progress will be printed out (may cause decrease in speed).')
    parser.add_argument('-x', '--noH', action='store_true', default=False,
                        help='if set this flag hydrogen atoms will be excluded from the simplexes calculation')
    parser.add_argument('-f', '--fragments', metavar='fragments.txt', default=None,
                        help='text file containing list of names of single compounds, fragment names and atom '
                             'indexes of fragment in the structure of corresponding compound')
    parser.add_argument('-w', '--id_field_name', metavar='field_name', default=None,
                        help='field name of unique ID for compounds (sdf) or reactions (rdf/rxn). '
                             'If omitted for sdf molecule titles will be used or auto-generated names; '
                             'for rdf $RIREG/$REREG field if it is not empty or auto-generated names')

    args = vars(parser.parse_args())
    opt_mix_ordered = None
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "nodict": opt_no_dict = v
        if o == "diff": opt_diff = v
        if o == "types":
            if v == "all": opt_types = list(range(1, 12))
            if v == "bounded": opt_types = [5, 6, 8, 9, 10, 11]
            if v == "extended": opt_types = list(range(3, 12))
        if o == "mixtures": mix_fname = v
        if o == "mix_ordered": opt_mix_ordered = v
        if o == "ncores": opt_ncores = cpu_count() if v == "all" else int(v)
        if o == "verbose": opt_verbose = v
        if o == "noH": opt_noH = v
        if o == "fragments": frag_fname = v
        if o == "stereo": parse_stereo = v
        if o == "quasi_mix": quasimix = v
        if o == "id_field_name": id_field_name = v
    if quasimix:
        opt_mix_ordered = False
        mix_fname = None
    if in_fname.split('.')[-1] in ['rdf', 'rxn']:
        opt_mix_ordered = False
        mix_fname = None
    opt_diff = [s.lower() for s in opt_diff]

    main_params(in_fname, out_fname, opt_no_dict, opt_diff, opt_types, mix_fname, opt_mix_ordered, opt_ncores,
                opt_verbose, opt_noH, frag_fname, parse_stereo, quasimix, id_field_name)


if __name__ == '__main__':
    main()
