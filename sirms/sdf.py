#-------------------------------------------------------------------------------
# Name:        sdf
# Purpose:     operations with sdf files (reading)
#
# Author:      Pavel Polishchuk
#
# Created:     11.01.2013
# Copyright:   (c) Pavel Polishchuk 2013-2015
# Licence:     BSD 3-clause
#-------------------------------------------------------------------------------


import re
from collections import OrderedDict

from .mols import SmilesMol3 as Mol
from .files import ReadPropertyRange, RangedLetter


formal_charges_table = {'0': 0,
                        '1': 3,
                        '2': 2,
                        '3': 1,
                        '4': 0,
                        '5': -1,
                        '6': -2,
                        '7': -3}


def molstr_to_Mol(molstr):

    mol = Mol()

    mol.title = molstr[0]

    natoms = int(molstr[3][0:3])
    nbonds = int(molstr[3][3:6])

    # read atoms
    id = 0
    for line in molstr[4:4 + natoms]:
        x = float(line[0:10])
        y = float(line[10:20])
        z = float(line[20:30])
        label = line[30:33].strip()
        formal_charge = line[36:39].strip()
        id += 1
        mol.AddAtom(id, label, x, y, z, formal_charges_table.get(formal_charge, 0))

    # read bonds
    for line in molstr[4 + natoms:4 + natoms + nbonds]:
        id1 = int(line[0:3])
        id2 = int(line[3:6])
        bond_type = int(line[6:9])
        mol.AddBond(id1, id2, bond_type)

    return mol


def add_property_to_atoms(mol, data_dict, fsetup):
    for prop_name in data_dict.keys():
        prop_range = ReadPropertyRange(fsetup, prop_name)
        for i, a in enumerate(sorted(mol.atoms.keys())):
            # if value is numeric use ranges, if value is string use it as label itself
            label = data_dict[prop_name][i] if type(data_dict[prop_name][i]) is str else RangedLetter(data_dict[prop_name][i], prop_range)
            # if label is represented as A|D three labels should be created A, D and AD
            label_full = label.split('|')
            if len(label_full) > 1:
                label_full.insert(0, label.replace('|', ''))
            mol.atoms[a]['property'][prop_name] = {'value': data_dict[prop_name][i],
                                                   'label': label_full}


def get_sdf_field(sdf_lines, field_name):

    # case-sensitive field name search

    i = 0

    patt = "^> +<%s>" % field_name

    # find field
    for i, line in enumerate(sdf_lines):
        if re.match(patt, line):
            break

    if i == len(sdf_lines) - 1:
        print('Field %s was not found in the input file.' % field_name)
        return None

    # retrieve data if field was found
    i += 1
    data = []
    while i < len(sdf_lines):
        if sdf_lines[i].find('>  <') == 0 or sdf_lines[i].find('> <') == 0 or sdf_lines[i].rstrip() == '$$$$':
            break
        if sdf_lines[i]:
            data.append(sdf_lines[i])
        i += 1

    return data


def ReadSDF(fname, id_field_name, opt_diff, fsetup, parse_stereo=False):
    """
    INPUT: sdf-filename
    OUTPUT: dict of molecules, where key is the title of the moleculae taken from the first line of mol-record
    """

    def _MolstrToMol(molstr, id_field_name, opt_diff, parse_stereo, num_mol):

        mol = molstr_to_Mol(molstr)

        # get mol name from sdf field (multi-line data fields are concatenated via dot)
        if id_field_name is not None:
            mol.title = '.'.join(get_sdf_field(molstr, id_field_name))
        if not mol.title:
            mol.title = 'auto_generated_id_' + str(num_mol)
        mol.stereo = parse_stereo

        # read properties from sdf fields
        data_dict = dict()
        for diff in opt_diff:
            data = get_sdf_field(molstr, diff)
            if data is not None and len(data) == 1:
                s = data[0].split(';')
                try:
                    data_dict[diff] = [float(el.replace(",", ".")) if el != "" else None for el in s]
                except ValueError:
                    data_dict[diff] = [el if el != "" else None for el in s]

        # add labels
        add_property_to_atoms(mol, data_dict, fsetup)

        # parse stereo
        if parse_stereo:

            data = get_sdf_field(molstr, 'stereoanalysis'.upper())
            if data is not None:
                for line in data:
                    tmp = line.split(' ')
                    if tmp[0] == 'CISTRANS':
                        # atoms enumeration in Chemaxon output starts from 0
                        id1 = int(tmp[1][1:-1]) + 1
                        id2 = int(tmp[2][:-1]) + 1
                        bond_stereo = tmp[-1].lower()
                        if bond_stereo == 'wiggly':
                            bond_stereo = 0
                        elif bond_stereo == 'e':
                            bond_stereo = 2
                        elif bond_stereo == 'z':
                            bond_stereo = 1
                        mol.SetDoubleBondConfig(id1, id2, bond_stereo)

            mol.SetCyclicDoubleBondsCis()

        return mol

    molstr = []
    i = 0
    with open(fname) as f:
        for line in f:
            if line.rstrip() != "$$$$":
                molstr.append(line.rstrip())
            else:
                i += 1
                yield _MolstrToMol(molstr, id_field_name, opt_diff, parse_stereo, i)
                molstr = []


def get_rdf_field(rxn, field_name):
    i = 0
    while rxn[i].find("$DTYPE") != 0:
        i += 1
    while rxn[i].split()[0] != "$DTYPE" and rxn[i].split()[1].lower() != field_name:
        i += 1
    i += 1
    field_data = ''
    # data can be on multiple lines
    while i < len(rxn) and rxn[i].find("$DTYPE") == -1:
        field_data += rxn[i].strip()
        i += 1
    return field_data.split(' ', 1)[1] if field_data else None


def parse_rxn_mols(rxn_lines):

    n_reactants = int(rxn_lines[4][0:3])
    n_products = int(rxn_lines[4][3:6])

    start_pos = 6
    pos = 6

    # read mols
    mols = []
    while n_reactants + n_products > len(mols):
        if rxn_lines[pos] != 'M  END':
            pos += 1
        else:
            mols.append(molstr_to_Mol(rxn_lines[start_pos:pos + 1]))
            start_pos = pos + 2
            pos = start_pos

    return mols, n_reactants, n_products


def set_mol_titles(mols, rx_id, n_reactants, n_products):
    for i in range(n_reactants):
        mols[i].title = rx_id + '_reactant_' + str(i + 1)
    for i in range(n_products):
        mols[n_reactants + i].title = rx_id + '_product_' + str(i + 1)


def create_rx_mix(mols, rx_id, n_reactants, n_products):
    mix = OrderedDict()
    mix[rx_id + '_reactants'] = {'names': [mols[i].title for i in range(n_reactants)],
                                 'ratios': [1] * n_reactants}
    mix[rx_id + '_products'] = {'names': [mols[n_reactants + i].title for i in range(n_products)],
                                'ratios': [1] * n_products}
    return mix


def process_rdf_reaction(rxn_lines, id_field_name, auto_name):

    mols, n_reactants, n_products = parse_rxn_mols(rxn_lines[1:])

    # get reaction id
    rx_id = None
    if id_field_name is not None:
        rx_id = get_rdf_field(rxn_lines, id_field_name)
    if rx_id is None:
        # from $RFMT $RIREG NUMBER
        tmp = rxn_lines[0].split()
        if len(tmp) == 3 and tmp[1] in ['$RIREG', '$REREG', '$MIREG', '$MEREG']:
                rx_id = tmp[2]
        else:
            rx_id = auto_name

    set_mol_titles(mols, rx_id, n_reactants, n_products)

    mix = create_rx_mix(mols, rx_id, n_reactants, n_products)

    # for now atomic properties cannot be considered for weighting of atoms in reactions
    # chemaxon cannot calculate atomic properties for reactions

    # read atom properties
    # data_dict = dict()
    # for diff in opt_diff:
    #     data = get_rdf_field(rxn_lines, diff)
    #     if data is not None:
    #         try:
    #             data_dict[diff] = [float(el.replace(",", ".")) if el != "" else None for el in data]
    #         except ValueError:
    #             data_dict[diff] = [el if el != "" else None for el in data]

    # add atom labels
    # if data_dict:
    #     n_atoms = 0
    #     for mol in mols:
    #         mol_data_dict = {k: v[n_atoms:len(mol.atoms)] for k, v in data_dict.items()}
    #         add_property_to_atoms(mol, mol_data_dict, fsetup)
    #         n_atoms += len(mol.atoms)

    return mols, mix


def process_rxn_reaction(rxn_lines, id_field_name, auto_name, opt_diff, fsetup):

    mols, n_reactants, n_products = parse_rxn_mols(rxn_lines)

    rx_id = auto_name
    if id_field_name is not None:
        data = get_sdf_field(rxn_lines, id_field_name)
        if data is not None and len(data) == 1:
            rx_id = data[0]

    # read properties from sdf fields
    data_dict = dict()
    for diff in opt_diff:
        data = get_sdf_field(rxn_lines, diff)
        if data is not None and len(data) == 1:
            s = data[0].split(';')
            try:
                data_dict[diff] = [float(el.replace(",", ".")) if el != "" else None for el in s]
            except ValueError:
                data_dict[diff] = [el if el != "" else None for el in s]

    # add labels
    start = 0
    for mol in mols:
        end = start + len(mol.atoms)
        data = {k: v[start:end] for k, v in data_dict.items()}
        add_property_to_atoms(mol, data, fsetup)
        start = end

    set_mol_titles(mols, rx_id, n_reactants, n_products)
    mix = create_rx_mix(mols, rx_id, n_reactants, n_products)

    return mols, mix


def ReadRDF(fname, id_field_name):
    # reaction ID will be searched in
    # 1) named field if specified (not None)
    # 2) $RFMT $RIREG id / $RFMT $REREG id if present
    # 3) generated automatically

    mols = OrderedDict()
    mix = OrderedDict()

    with open(fname) as rdf:

        # read file header
        rdf.readline()
        rdf.readline()

        rxn = []
        i = 1

        for line in rdf:
            if line.find("$RFMT") == 0:
                if rxn:
                    rx_mols, rx_mix = process_rdf_reaction(rxn, id_field_name, 'rx_autogen_id_' + str(i))
                    for mol in rx_mols:
                        mols[mol.title] = mol
                    mix.update(rx_mix)
                    i += 1
                rxn = [line.rstrip()]
            else:
                rxn.append(line.rstrip())

        rx_mols, rx_mix = process_rdf_reaction(rxn, id_field_name, 'rx_autogen_id_' + str(i))
        for mol in rx_mols:
            mols[mol.title] = mol
        mix.update(rx_mix)

    return mols, mix


def ReadRXN(fname, id_field_name, opt_diff, fsetup):

    mols = OrderedDict()
    mix = OrderedDict()

    with(open(fname)) as f:

        rxn = []
        i = 1

        for line in f:
            if line.find("$RXN") == 0:
                if rxn:
                    rx_mols, rx_mix = process_rxn_reaction(rxn, id_field_name, 'rx_autogen_id_' + str(i), opt_diff, fsetup)
                    for mol in rx_mols:
                        mols[mol.title] = mol
                    mix.update(rx_mix)
                    i += 1
                rxn = [line.rstrip()]
            else:
                rxn.append(line.rstrip())

        rx_mols, rx_mix = process_rxn_reaction(rxn, id_field_name, 'rx_autogen_id_' + str(i), opt_diff, fsetup)
        for mol in rx_mols:
            mols[mol.title] = mol
        mix.update(rx_mix)

    return mols, mix