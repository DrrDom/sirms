__author__ = 'pavel'
# Basic functions to work with descriptors names


def get_smile(sirms_string):
    return sirms_string.rsplit('|', 1)[1]


def get_atomcount(sirms_string):
    return int(sirms_string.split('|', 6)[5])


def get_atom_labeling(sirms_string):
    return sirms_string.rsplit('|', 2)[1]


def get_mix_single(sirms_string):
    return sirms_string.split('|', 2)[1]


def invert_num_prob_type(sirms_string):
    tmp = sirms_string.split('|')
    if tmp[2] == 'n':
        tmp[2] = 'p'
    else:
        tmp[2] = 'n'
    return '|'.join(tmp)


def insert_reaction_info(sirms_string, prod_react):
    """
    prod_react: p or r or pr
    """
    return prod_react + '|' + sirms_string.split('|', 1)[1]


def split_by_reaction_part(sirms_name):
    return sirms_name.split('|', 1)


def join_reaction_part(react_part, common_part):
    return react_part + '|' + common_part


def create_full_name(prod_react, single_mix, num_prob, atom_count, atom_labeling, smiles):
    return prod_react + '|' + single_mix + '|' + num_prob + '|||' + str(atom_count) + '|||' + \
           atom_labeling + '|' + smiles


