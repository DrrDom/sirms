#!/usr/bin/env python3
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

import argparse
import time
import files
import copy
from itertools import combinations, chain, product, combinations_with_replacement
from collections import OrderedDict

from sdf import ReadSDF, ReadRDF, ReadRXN
from labels import SetLabelsInternal, builtin_types
from ppgfunctions import *
import sirms_name

mol_frag_sep = "###"


#===============================================================================
# Save simplexes

def sort_lists_by(*lists, key_list=0, desc=False):
    return zip(*sorted(zip(*lists), reverse=desc, key=lambda x: x[key_list]))


def SaveSimplexes(fname, sirms, output_format, ndigits=5):
    """
    Save calculated decriptors in the tab-delimited text file format
    Descriptors   descriptor_1  descriptor_2  ...  descriptor_n
    compound_1               5             2                  0
    compound_2               1             1                  4
    ...

    or in svm sparse format with two additional files with extensions colnames and rownames
    zero based indices
    0:12 1:23 4:1 ...
    """

    if output_format == 'txt':
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
    elif output_format == 'svm':
        sirms_names = sorted(list(set(list(chain.from_iterable([list(s.keys()) for s in sirms.values()])))))
        colnames = sirms_names
        rownames = sirms.keys()
        sirms_names = dict(zip(sirms_names, range(len(sirms_names))))
        s = "{:." + str(ndigits) + "f}"
        missing_rows = []  # contains mols with all zero descriptors to omit in output
        with open(fname, 'wt') as f:
            for i, name in enumerate(sirms.keys()):
                if len(sirms[name]) == 0:
                    missing_rows.append(i)
                    continue
                # get new integer names of simplexes
                int_names = [sirms_names[v] for v in sirms[name].keys()]
                values = []
                for value in sirms[name].values():
                    # round floating point values, all others will be remained unchanged
                    if type(value) is int or value.is_integer():
                        values.append(str(int(value)))
                    elif isinstance(value, float):
                        values.append(s.format(value))
                    else:
                        values.append(str(value))
                int_names, values = sort_lists_by(int_names, values)
                f.write(' '.join([str(n) + ':' + v for n, v in zip(int_names, values)]) + '\n')
        open(os.path.splitext(fname)[0] + '.colnames', 'wt').write('\n'.join(colnames))
        open(os.path.splitext(fname)[0] + '.rownames', 'wt').write(
            '\n'.join([row for i, row in enumerate(rownames) if i not in missing_rows]))


#===============================================================================
# create reaction sirms
def concat_reaction_sirms(sirms):

    output = OrderedDict()

    while sirms:

        # take first element
        k, v = sirms.popitem(last=False)
        rx_id, role = k.rsplit('_', 1)
        # first letter from role is used as additional description
        new_sirms_names = [sirms_name.insert_reaction_info(s, role[0]) for s in sorted(v.keys())]
        s1 = {new_sirms_names[i]: v[sname] for i, sname in enumerate(sorted(v.keys()))}

        # take complimentary element (reactants or products of the same reaction)
        if role == 'reactants':
            k = rx_id + '_products'
            new_sirms_names = [sirms_name.insert_reaction_info(s, 'p') for s in sorted(sirms[k].keys())]
        elif role == 'products':
            k = rx_id + '_reactants'
            new_sirms_names = [sirms_name.insert_reaction_info(s, 'r') for s in sorted(sirms[k].keys())]
        else:
            print('Impossible error')
        v = sirms[k]
        s2 = {new_sirms_names[i]: v[sname] for i, sname in enumerate(sorted(sirms[k].keys()))}
        del (sirms[k])

        # combine them
        s1.update(s2)
        output[rx_id] = s1

    return output


#===============================================================================
# SIRMS

def SetLabelsExternal(mols, opt_diff, input_fname):
    for s_diff in opt_diff:
        files.LoadRangedProperty(mols, GetWorkDir(input_fname), GetFileNameNoExt(input_fname) + '.' + s_diff)


def CalcMolSingleSirms(mol, diff_list, min_num_atoms, max_num_atoms, min_num_components,
                       max_num_components, noH, verbose, frags=None):
    """
    INPUT:
        mol: Mol3 object;
        sirms_dict: loaded precomputed dictionary;
        diff_list: list of labeling to compute;
    OUTPUT:
        returns dict of dicts
        if frags are specified then it returns descriptors for fragment-depleted molecule
        { "mol_name":               { "descriptors_1": 0, "descriptor_2": 5, ... },
          "mol_name#frag_name_1":   { "descriptors_1": 0, "descriptor_2": 3, ... },
          ...
        }

        in case of for_mix=True
        returns additional dict with sirms names and atom counts
    """

    if verbose:
        cur_time = time.time()

    d = {mol.title: {}}
    if frags is not None:
        local_frags = frags.get(mol.title, None)
        if local_frags is not None:
            for frag in local_frags.keys():
                d[frag] = {}
    else:
        local_frags = None

    for a in mol.GetAtomsCombinations(min_num_components=min_num_components, max_num_components=max_num_components,
                                      min_num_atoms=min_num_atoms, max_num_atoms=max_num_atoms, noH=noH):

        for s_diff in diff_list:
            labels = [mol.atoms[a_id]['property'][s_diff]['label'] for a_id in a]
            for labels_set in product(*labels):
                s_name = sirms_name.create_full_name(prod_react='', single_mix='S', num_prob='n',
                                                     atom_count=len(a), atom_labeling=s_diff,
                                                     smiles=mol.get_name(a, labels_set))
                # sirms_name = 'S|' + s_diff + '|' + mol.get_name(a, labels_set)
                d[mol.title][s_name] = d[mol.title].get(s_name, 0) + 1
                if local_frags is not None:
                    # if there is no common atoms in simplex and fragment
                    for frag in local_frags.keys():
                        if set(a).isdisjoint(local_frags[frag]):
                            d[frag][s_name] = d[frag].get(s_name, 0) + 1

    # rename frags in output dict
    output = {}
    for n in d.keys():
        if n == mol.title:
            output[mol.title] = d[n]
        else:
            output[mol.title + mol_frag_sep + n] = d[n]
    del (d)

    if verbose:
        print(mol.title, (str(round(time.time() - cur_time, 1)) + ' s').rjust(78 - len(mol.title)))

    return output


#===============================================================================

def CalcSingleSirms(mol_list, opt_diff, min_num_atoms, max_num_atoms, min_num_components,
                    max_num_components, opt_noH, opt_verbose, frags):

    if opt_verbose:
        print(' Single compounds calculation '.center(79, '='))

    sirms = OrderedDict()
    for mol in mol_list:
        res = CalcMolSingleSirms(mol, opt_diff, min_num_atoms, max_num_atoms, min_num_components,
                                 max_num_components, opt_noH, opt_verbose, frags)
        sirms.update(res)

    return sirms


#===============================================================================

def CalcMixSirms(single_sirms, mix, atom_labeling, min_num_atoms=4, max_num_atoms=4,
                 min_num_mix_components=2, max_num_mix_components=2, verbose=False, ordered=False,
                 self_assembly_mix=False):

    def mult(values):
        p = 1
        for item in values:
            p *= item
        return p

    def gen_mix_sirms_name(sirms_names, ordered, ids=None):
        if not ordered:
            tmp = [sirms_name.get_smile(item) for item in sirms_names]
            tmp = [item.split('.') for item in tmp]
            # flatten list and sort to obtain canonical order
            tmp = sorted([item for sublist in tmp for item in sublist])
        else:
            # add component id before smile
            tmp = sorted([str(id + 1) + '_' + sirms_name.get_smile(item) for id, item in zip(ids, sirms_names)])
        # get atom count
        atom_count = sum(sirms_name.get_atomcount(item) for item in sirms_names)
        # get atomic property of descriptors (they all must have the same)
        return sirms_name.create_full_name(prod_react='', single_mix='M', num_prob='n', atom_count=atom_count,
                                           atom_labeling=sirms_name.get_atom_labeling(sirms_names[0]),
                                           smiles='.'.join(tmp))

    def select_sirms_by_labeling_one_component(sirms_names, labeling_name):
        substr = '|' + labeling_name + '|'
        return [item for item in sirms_names if item.find(substr) > -1 and item.split('|')[-1].find('.') == -1]

    if verbose:
        print(' Mixtures calculation '.center(79, '='))

    d = OrderedDict()

    for mix_name, mix_data in mix.items():

        if verbose:
            cur_time = time.time()

        if self_assembly_mix:
            iterator = combinations_with_replacement
        else:
            iterator = combinations

        d_mix = dict()

        # multiply descriptors on compounds ratios
        d_tmp = dict()
        for i, mol_name in enumerate(mix_data['names']):
            d_tmp[mol_name] = copy.deepcopy(single_sirms[mol_name])
            for name, value in d_tmp[mol_name].items():
                d_tmp[mol_name][name] = value * mix_data['ratios'][i]

        for n in range(min_num_mix_components, max_num_mix_components + 1):
            for ids in iterator(range(len(mix_data['names'])), n):   # ids - mol id in mixture
                comb = [mix_data['names'][i] for i in ids]
                for labeling in atom_labeling:
                    for p in product(*[select_sirms_by_labeling_one_component(d_tmp[mol_name].keys(), labeling) for mol_name in comb]):  # p - combination of sirms names from molecules
                        s = sum(sirms_name.get_atomcount(item) for item in p)
                        if min_num_atoms <= s <= max_num_atoms:
                            mix_sirs_name = gen_mix_sirms_name(p, ordered, ids)
                            d_mix[mix_sirs_name] = d_mix.get(mix_sirs_name, 0) + mult([d_tmp[mol_name].get(p[i], 0) for i, mol_name in enumerate(comb)])

        # add single sirms and filter them with given number of atoms
        sirms_names = set(list(chain.from_iterable(list(s.keys()) for s in single_sirms.values())))
        for name in sirms_names:
            if min_num_atoms <= sirms_name.get_atomcount(name) <= max_num_atoms:
                d_mix[name] = sum(single_sirms[mol_name].get(name, 0) for mol_name in mix_data['names'])

        d[mix_name] = d_mix

        if verbose:
            print(mix_name, (str(round(time.time() - cur_time, 1)) + ' s').rjust(78 - len(mix_name)))

    return d


#===============================================================================

def GenQuasiMix(mol_names):
    d = OrderedDict()
    for n in mol_names:
        d[n] = {'names': [n, n], 'ratios': [0.5, 0.5]}
    return d


#===============================================================================

def CalcReactionDiffSirms(sirms):

    def calc_diff(d):
        output = dict()
        checked_names = set()
        for name, value in d.items():
            r, common_part = sirms_name.split_by_reaction_part(name)
            if common_part not in checked_names:
                checked_names.add(common_part)
                if r == 'p':
                    res = d[name] - d.get(sirms_name.join_reaction_part('r', common_part), 0)
                elif r == 'r':
                    res = d.get(sirms_name.join_reaction_part('p', common_part), 0) - d[name]
                else:
                    print("Incorrect reaction_part label %s, can be only r, p or pr." % r)
                    continue
                if res != 0:
                    output[sirms_name.join_reaction_part('pr', common_part)] = res
        return output


    for mol_name, data in sirms.items():
        sirms[mol_name] = calc_diff(sirms[mol_name])
    return sirms


#===============================================================================

def CalcProbSirms(sirms, type):
    if type in ['prob', 'both']:
        for mol_name in sirms.keys():
            d = dict()
            s_mix = sum(v for k, v in sirms[mol_name].items() if sirms_name.get_mix_single(k) == 'M')
            s_single = sum(v for k, v in sirms[mol_name].items() if sirms_name.get_mix_single(k) == 'S')
            for k, v in sirms[mol_name].items():
                if sirms_name.get_mix_single(k) == 'M':
                    d[sirms_name.invert_num_prob_type(k)] = v/s_mix if s_mix != 0 else 0
                else:
                    d[sirms_name.invert_num_prob_type(k)] = v/s_single if s_single != 0 else 0
            if type == 'both':
                sirms[mol_name].update(d)
            elif type == 'prob':
                sirms[mol_name] = d
    return sirms


#===============================================================================
# Main cycle

def main_params(in_fname, out_fname, opt_diff, min_num_atoms, max_num_atoms, min_num_components,
                max_num_components, min_num_mix_components, max_num_mix_components, mix_fname,
                descriptors_transformation, mix_type, opt_mix_ordered, opt_verbose, opt_noH,
                frag_fname, self_association_mix, reaction_diff, quasimix, id_field_name, output_format):

    # define which property will be loaded from external file or from sdf-file
    opt_diff_builtin = [v for v in opt_diff if v in builtin_types]
    opt_diff_sdf = files.NotExistedPropertyFiles([v for v in opt_diff if v not in opt_diff_builtin], in_fname)
    opt_diff_ext = [el for el in opt_diff if el not in opt_diff_sdf and el not in opt_diff_builtin]

    # load sdf, rdf or rxn file depending on its extension
    input_file_extension = in_fname.strip().split(".")[-1].lower()
    setup_path = os.path.join(GetWorkDir(in_fname), "setup.txt")
    if input_file_extension == 'sdf':
        mols = ReadSDF(in_fname, id_field_name, opt_diff_sdf, setup_path)
    elif input_file_extension == 'rdf':
        mols, mix = ReadRDF(in_fname, id_field_name)
    elif input_file_extension == 'rxn':
        mols, mix = ReadRXN(in_fname, id_field_name, opt_diff_sdf, setup_path)
    else:
        print("Input file extension should be SDF, RDF or RXN. Current file has %s. Please check it." %
              input_file_extension.upper())
        return None

    # set property labels on atoms from external data files
    SetLabelsExternal(mols, opt_diff_ext, in_fname)
    # set labels of built-in types
    SetLabelsInternal(mols, opt_diff_builtin, setup_path)

    if mix_fname is None and not quasimix and input_file_extension == 'sdf':

        frags = files.LoadFragments(frag_fname)
        # calc simplex descriptors
        sirms = CalcSingleSirms(mols.values(), opt_diff, min_num_atoms, max_num_atoms, min_num_components,
                                max_num_components, opt_noH, opt_verbose, frags)

    else:

        # create mix data for sdf (for rdf/rxn mix is created during file loading)
        if input_file_extension == 'sdf':
            if quasimix:
                mix = GenQuasiMix(list(mols.keys()))
            elif mix_fname is not None:
                mix = files.LoadMixturesTxt(mix_fname, mix_type)
            else:
                print("Strange error occurred during mix preparation")
                exit()

        mols_used = set(chain.from_iterable([m['names'] for m in mix.values()]))

        # min_num_atoms set to 1 to be able to generate mixtures
        sirms = CalcSingleSirms([mols[mol_name] for mol_name in mols_used], opt_diff, 1,
                                max_num_atoms, min_num_components, max_num_components, opt_noH,
                                opt_verbose, None)

        sirms = CalcMixSirms(single_sirms=sirms,
                             mix=mix,
                             atom_labeling=opt_diff,
                             min_num_atoms=min_num_atoms,
                             max_num_atoms=max_num_atoms,
                             min_num_mix_components=min_num_mix_components,
                             max_num_mix_components=max_num_mix_components,
                             verbose=opt_verbose,
                             ordered=opt_mix_ordered,
                             self_assembly_mix=self_association_mix)

        # filter single sirms with number of atoms lower than min_num_atoms

        if input_file_extension in ['rdf', 'rxn']:
            sirms = concat_reaction_sirms(sirms)

    if descriptors_transformation in ['prob', 'both']:
        sirms = CalcProbSirms(sirms, descriptors_transformation)

    if reaction_diff:
        sirms = CalcReactionDiffSirms(sirms)

    SaveSimplexes(out_fname, sirms, output_format)


def main():
    parser = argparse.ArgumentParser(description='Calculate simplex descriptors for single molecules, quasi-mixtures, '
                                                 'mixtures and reactions.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
                        help='input file (allowed formats: sdf, rdf, rxn) with standardized structures, '
                             'molecules or reactions should have titles.')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='output file with calculated descriptors. Can be in text or sparse svm format.')
    parser.add_argument('-b', '--output_format', metavar='format_name', default='txt',
                        help='format of output file with calculated descriptors (txt|svm). '
                             'Txt is ordinary tab-separated text file. Svm is sparse format, two additional file will '
                             'be saved with extensions .colnames amd .rownames. Default: txt.')
    parser.add_argument('-a', '--atoms_labeling', metavar='', default=['elm'], nargs='*',
                        help='list of atom labeling schemes separated by space. Built-in scheme is element (elm) and '
                             'topology (none). '
                             'To include other schemes user should specify the name of the corresponding property '
                             'value identical to the name of SDF field, which contains calculated atomic properties. '
                             'Fields names are case-sensitive. For RDF/RXN input files only built-in types can be '
                             'used. Default: elm.')
    parser.add_argument('--min_atoms', metavar='', default=4,
                        help='The minimal number of atoms in the fragment. Default: 4')
    parser.add_argument('--max_atoms', metavar='', default=4,
                        help='The maximal number of atoms in the fragment. Default: 4')
    parser.add_argument('--min_components', metavar='', default=1,
                        help='The minimal number of disconnected groups of atoms in the fragment. '
                             'Default: 1 (mean fully connected fragments).')
    parser.add_argument('--max_components', metavar='', default=2,
                        help='The maximal number of disconnected groups of atoms in the fragment. '
                             'Default: 2.')
    parser.add_argument('-q', '--quasi_mix', action='store_true', default=False,
                        help='calculate quasi-mixture descriptors for pure compounds.')
    parser.add_argument('-m', '--mixtures', metavar='mixtures.txt', default=None,
                        help='text file containing list of mixtures of components and their ratios. Names of '
                             'components should be the same as in input.sdf file.')
    parser.add_argument('--min_mix_components', metavar='', default=2,
                        help='The minimal number of molecules which contribute to mixture fragments. '
                             'Default: 2 (and cannot be less).')
    parser.add_argument('--max_mix_components', metavar='', default=2,
                        help='The maximal number of molecules which contribute to mixture fragments. '
                             'Default: 2 (take into account only binary interactions in a mixture).')
    parser.add_argument('--descriptors_transformation', metavar='num|prob|both', default='num',
                        help='num: numbers of fragments (for single compounds) or number of fragments combinations '
                             'weighted by their molar ratios (for mixtures). '
                             'prob: final value of descriptors are divided on sum of all descriptors to '
                             'describe probability of each descriptor. Descriptors for single compounds and mixtures '
                             'are weighted separately. both: will generate both types of descriptors. Default: num.')
    parser.add_argument('--mix_type', metavar='abs|rel', default='abs',
                        help='abs: means that amount of components given in a mixture file will be considered as is. '
                             'rel: means that amount of components given in a mixture file will be taken as relative '
                             'amount and will be scaled to the sum of 1. Default: abs.')
    parser.add_argument('-r', '--mix_ordered', action='store_true', default=False,
                        help='if set this flag the mixtures will be considered ordered, otherwise as unordered. '
                             'In ordered mixtures role of each component is known and position of a component in a '
                             'mixture description file will be taken into account. In unordered mixtures all '
                             'components are equitable and their roles don''t depend on their positions in mixture '
                             'description. Used only in combination with -m key.')
    parser.add_argument('--mix_self_association', action='store_true', default=False,
                        help='calculates mixture descriptors between components with themselves in order to take into '
                             'account self-interaction of components. Default: false.')
    parser.add_argument('--reaction_diff', action='store_true', default=False,
                        help='if set this flag difference between product and reactant descriptors will be calculated. '
                             'By default these two feature vectors are concatenated.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='if set this flag progress will be printed out (may cause decrease in speed).')
    parser.add_argument('-x', '--noH', action='store_true', default=False,
                        help='if set this flag hydrogen atoms will be excluded from the simplexes calculation.')
    parser.add_argument('-f', '--fragments', metavar='fragments.txt', default=None,
                        help='text file containing list of names of single compounds, fragment names and atom '
                             'indexes of fragment to remove (all values are tab-separated).')
    parser.add_argument('-w', '--id_field_name', metavar='field_name', default=None,
                        help='field name of unique ID for compounds (sdf) or reactions (rdf/rxn). '
                             'If omitted for sdf molecule titles will be used or auto-generated names; '
                             'for rdf $RIREG/$REREG/$MIREG/$MEREG field if it is not empty or auto-generated names')

    args = vars(parser.parse_args())
    opt_mix_ordered = None
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "atoms_labeling": opt_diff = v
        if o == "mixtures": mix_fname = v
        if o == "mix_ordered": opt_mix_ordered = v
        if o == "verbose": opt_verbose = v
        if o == "noH": opt_noH = v
        if o == "fragments": frag_fname = v
        if o == "stereo": parse_stereo = v
        if o == "quasi_mix": quasimix = v
        if o == "id_field_name": id_field_name = v
        if o == "output_format": output_format = v
        if o == "min_atoms": min_num_atoms = int(v)
        if o == "max_atoms": max_num_atoms = int(v)
        if o == "min_components": min_num_components = int(v)
        if o == "max_components": max_num_components = int(v)
        if o == "min_mix_components": min_num_mix_components = int(v)
        if o == "max_mix_components": max_num_mix_components = int(v)
        if o == "mix_self_association": self_association_mix = v
        if o == "descriptors_transformation": descriptors_transformation = v
        if o == "mix_type": mix_type = v
        if o == "reaction_diff": reaction_diff = v
    if quasimix:
        opt_mix_ordered = False
        mix_fname = None
    if in_fname.split('.')[-1] in ['rdf', 'rxn']:
        opt_mix_ordered = False
        mix_fname = None
    if output_format not in ['svm', 'txt']:
        print("INPUT ERROR: wrong output format specified - %s. On;y txt or svm are allowed." % output_format)
        exit()
    if min_num_atoms > max_num_atoms:
        print("INPUT ERROR: min_num_atoms should not be greater than max_num_atoms.")
        exit()
    if min_num_components > max_num_components:
        print("INPUT ERROR: min_num_components should not be greater than max_num_components.")
        exit()
    if min_num_mix_components > max_num_mix_components or min_num_mix_components < 2:
        print("INPUT ERROR: min_num_mix_components should not be greater than max_num_mix_components and "
              "minimal value of min_num_mix_components should be 2 or greater.")
        exit()
    if min_num_mix_components > max_num_components:
        print("INPUT ERROR: minimal number of mixture components (min_num_mix_components) cannot be greater than "
              "maximal number of components in a fragment descriptor (min_num_components).")
        exit()
    if mix_type not in ['abs', 'rel']:
        print("INPUT ERROR: mixture type (mix_type) can be only abs or rel.")
        exit()
    if descriptors_transformation not in ['num', 'prob', 'both']:
        print("INPUT ERROR: type of mixture descriptors (types_mix_descriptors) can be only num, prob or both.")
        exit()
    if in_fname.split('.')[-1] not in ['rdf', 'rxn'] and reaction_diff:
        print("The option reaction diff can be used only with RDF or RXN input files.")
        exit()

    main_params(in_fname=in_fname,
                out_fname=out_fname,
                opt_diff=opt_diff,
                min_num_atoms=min_num_atoms,
                max_num_atoms=max_num_atoms,
                min_num_components=min_num_components,
                max_num_components=max_num_components,
                min_num_mix_components=min_num_mix_components,
                max_num_mix_components=max_num_mix_components,
                mix_fname=mix_fname,
                descriptors_transformation=descriptors_transformation,
                mix_type=mix_type,
                opt_mix_ordered=opt_mix_ordered,
                opt_verbose=opt_verbose,
                opt_noH=opt_noH,
                frag_fname=frag_fname,
                self_association_mix=self_association_mix,
                reaction_diff=reaction_diff,
                quasimix=quasimix,
                id_field_name=id_field_name,
                output_format=output_format)


if __name__ == '__main__':
    main()
