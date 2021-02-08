#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 20-01-2016
# version         : 0.1
# python_version  : 3
# copyright       : Pavel Polishchuk 2016
# license         : BSD 3-clause
#==============================================================================

from .files import ReadPropertyRange, RangedLetter

elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
            'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
            'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',
            'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
            'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
            'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt',
            'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
            'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf',
            'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Uub', 'Uut', 'Uuq', 'Uup', 'Uuh',
            'Uus', 'Uuo']

uff_dist = {'H': 2.886, 'He': 2.362, 'Li': 2.451, 'Be': 2.745, 'B': 4.083, 'C': 3.851, 'N': 3.66, 'O': 3.5, 'F': 3.364,
       'Ne': 3.243, 'Na': 2.983, 'Mg': 3.021, 'Al': 4.499, 'Si': 4.295, 'P': 4.147, 'S': 4.035, 'Cl': 3.947,
       'Ar': 3.868, 'K': 3.812, 'Ca': 3.399, 'Sc': 3.295, 'Ti': 3.175, 'V': 3.144, 'Cr': 3.023, 'Mn': 2.961,
       'Fe': 2.912, 'Co': 2.872, 'Ni': 2.834, 'Cu': 3.495, 'Zn': 2.763, 'Ga': 4.383, 'Ge': 4.28, 'As': 4.23,
       'Se': 4.205, 'Br': 4.189, 'Kr': 4.141, 'Rb': 4.114, 'Sr': 3.641, 'Y': 3.345, 'Zr': 3.124, 'Nb': 3.165,
       'Mo': 3.052, 'Tc': 2.998, 'Ru': 2.963, 'Rh': 2.929, 'Pd': 2.899, 'Ag': 3.148, 'Cd': 2.848, 'In': 4.463,
       'Sn': 4.392, 'Sb': 4.42, 'Te': 4.47, 'I': 4.5, 'Xe': 4.404, 'Cs': 4.517, 'Ba': 3.703, 'La': 3.522, 'Ce': 3.556,
       'Pr': 3.606, 'Nd': 3.575, 'Pm': 3.547, 'Sm': 3.52, 'Eu': 3.493, 'Gd': 3.368, 'Tb': 3.451, 'Dy': 3.428,
       'Ho': 3.409, 'Er': 3.391, 'Tm': 3.374, 'Yb': 3.355, 'Lu': 3.64, 'Hf': 3.141, 'Ta': 3.17, 'W': 3.069,
       'Re': 2.954, 'Os': 3.12, 'Ir': 2.84, 'Pt': 2.754, 'Au': 3.293, 'Hg': 2.705, 'Tl': 4.347, 'Pb': 4.297,
       'Bi': 4.37, 'Po': 4.709, 'At': 4.75, 'Rn': 4.765, 'Fr': 4.9, 'Ra': 3.677, 'Ac': 3.478, 'Th': 3.396, 'Pa': 3.424,
       'U': 3.395, 'Np': 3.424, 'Pu': 3.424, 'Am': 3.381, 'Cm': 3.326, 'Bk': 3.339, 'Cf': 3.313, 'Es': 3.299,
       'Fm': 3.286, 'Md': 3.274, 'No': 3.248, 'Lr': 3.236}

uff_energy = {'H': 0.044, 'He': 0.056, 'Li': 0.025, 'Be': 0.085, 'B': 0.18, 'C': 0.105, 'N': 0.069, 'O': 0.06,
              'F': 0.05, 'Ne': 0.042, 'Na': 0.03, 'Mg': 0.111, 'Al': 0.505, 'Si': 0.402, 'P': 0.305, 'S': 0.247,
              'Cl': 0.227, 'Ar': 0.185, 'K': 0.035, 'Ca': 0.238, 'Sc': 0.019, 'Ti': 0.017, 'V': 0.016, 'Cr': 0.015,
              'Mn': 0.013, 'Fe': 0.013, 'Co': 0.014, 'Ni': 0.015, 'Cu': 0.005, 'Zn': 0.124, 'Ga': 0.415, 'Ge': 0.379,
              'As': 0.309, 'Se': 0.291, 'Br': 0.251, 'Kr': 0.22, 'Rb': 0.04, 'Sr': 0.235, 'Y': 0.072, 'Zr': 0.069,
              'Nb': 0.059, 'Mo': 0.056, 'Tc': 0.048, 'Ru': 0.056, 'Rh': 0.053, 'Pd': 0.048, 'Ag': 0.036, 'Cd': 0.228,
              'In': 0.599, 'Sn': 0.567, 'Sb': 0.449, 'Te': 0.398, 'I': 0.339, 'Xe': 0.332, 'Cs': 0.045, 'Ba': 0.364,
              'La': 0.017, 'Ce': 0.013, 'Pr': 0.01, 'Nd': 0.01, 'Pm': 0.009, 'Sm': 0.008, 'Eu': 0.008, 'Gd': 0.009,
              'Tb': 0.007, 'Dy': 0.007, 'Ho': 0.007, 'Er': 0.007, 'Tm': 0.006, 'Yb': 0.228, 'Lu': 0.041, 'Hf': 0.072,
              'Ta': 0.081, 'W': 0.067, 'Re': 0.066, 'Os': 0.037, 'Ir': 0.073, 'Pt': 0.08, 'Au': 0.039, 'Hg': 0.385,
              'Tl': 0.68, 'Pb': 0.663, 'Bi': 0.518, 'Po': 0.325, 'At': 0.284, 'Rn': 0.248, 'Fr': 0.05, 'Ra': 0.404,
              'Ac': 0.033, 'Th': 0.026, 'Pa': 0.022, 'U': 0.022, 'Np': 0.019, 'Pu': 0.016, 'Am': 0.014, 'Cm': 0.013,
              'Bk': 0.013, 'Cf': 0.013, 'Es': 0.012, 'Fm': 0.012, 'Md': 0.011, 'No': 0.011, 'Lr': 0.011}

builtin_types = ('elm', 'none', 'uffd', 'uffe')

builtin_continuous_types = ('uffd', 'uffe')


def GetSetupRanges(opt_diff, setup_path):
    ranges = dict()
    for diff in opt_diff:
        if diff in builtin_continuous_types:
            range = ReadPropertyRange(setup_path, diff)
            if range is not None:
                ranges[diff] = range
    return ranges


def SetLabelsInternalToMol(mol, opt_diff, setup_ranges):
    for atom in mol.atoms.values():
        if 'elm' in opt_diff:
            atom['property']['elm'] = {'label': [atom['label']], 'value': atom['label']}
        if 'none' in opt_diff:
            atom['property']['none'] = {'label': ['A'], 'value': 'A'}
        if 'uffd' in opt_diff:
            value = uff_dist[atom['label']]
            atom['property']['uffd'] = {'label': [RangedLetter(value, setup_ranges['uffd'])], 'value': value}
        if 'uffe' in opt_diff:
            value = uff_energy[atom['label']]
            atom['property']['uffe'] = {'label': [RangedLetter(value, setup_ranges['uffe'])], 'value': value}


def SetLabelsInternal(mols, opt_diff, setup_path):
    ranges = GetSetupRanges(opt_diff, setup_path)
    for mol in mols.values():
        SetLabelsInternalToMol(mol, opt_diff, ranges)
