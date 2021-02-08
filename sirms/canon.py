#-------------------------------------------------------------------------------
# Name:        canon
# Purpose:     generate canonic sirms names and stereodescriptors
#
# Author:      Pavel Polishchuk
#
# Created:     11.01.2013
# Copyright:   (c) Pavel Polishchuk 2013
# Licence:     BSD 3-clause
#-------------------------------------------------------------------------------

from itertools import combinations, permutations, chain, product
from collections import defaultdict
from math import acos, sqrt, degrees
from .ppgfunctions import SortTwoLists
import os
import json

# Global constants

dumb_label_set = ['A', 'B', 'C', 'D']

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

#base simplex types represented as sets of vertex neighbors
s_types = {(0, 0, 0, 0): 1, (0, 0, 1, 1): 2, (1, 1, 1, 1): 3, (0, 1, 1, 2): 4, (1, 1, 1, 3): 5,
           (1, 1, 2, 2): 6, (0, 2, 2, 2): 7, (1, 2, 2, 3): 8, (2, 2, 2, 2): 9, (2, 2, 3, 3): 10,
           (3, 3, 3, 3): 11}

bin_s_types = {(1, 0, 0, 1, 1, 1): 8, (0, 1, 0, 1, 1, 0): 6, (0, 1, 1, 1, 0, 1): 8, (1, 1, 0, 0, 1, 0): 6,
               (0, 0, 1, 0, 1, 0): 4, (0, 1, 1, 0, 0, 0): 4, (1, 0, 0, 0, 0, 0): 2, (0, 1, 1, 1, 1, 1): 10,
               (1, 0, 0, 1, 0, 0): 4, (1, 1, 1, 1, 1, 0): 10, (0, 0, 1, 0, 0, 1): 4, (0, 1, 1, 0, 1, 1): 8,
               (0, 0, 1, 1, 0, 0): 3, (0, 1, 1, 1, 0, 0): 6, (1, 0, 1, 1, 0, 1): 9, (0, 0, 1, 0, 1, 1): 5,
               (1, 1, 1, 0, 1, 0): 8, (0, 1, 1, 0, 0, 1): 7, (1, 0, 1, 0, 0, 0): 4, (0, 0, 1, 1, 1, 1): 8,
               (1, 1, 1, 1, 0, 0): 8, (1, 0, 1, 1, 1, 1): 10, (1, 0, 0, 1, 1, 0): 5, (0, 0, 1, 0, 0, 0): 2,
               (1, 1, 1, 0, 0, 1): 8, (0, 1, 1, 0, 1, 0): 6, (1, 0, 1, 0, 1, 1): 8, (0, 0, 1, 1, 0, 1): 6,
               (0, 0, 0, 1, 1, 0): 4, (1, 0, 1, 1, 0, 0): 6, (1, 1, 1, 0, 1, 1): 10, (1, 0, 1, 0, 0, 1): 6,
               (1, 1, 1, 1, 0, 1): 10, (0, 0, 1, 1, 1, 0): 6, (0, 0, 0, 1, 0, 1): 4, (0, 0, 0, 0, 0, 0): 1,
               (1, 1, 1, 0, 0, 0): 5, (1, 0, 1, 0, 1, 0): 7, (0, 1, 0, 0, 1, 0): 3, (0, 0, 0, 1, 1, 1): 7,
               (1, 1, 0, 1, 1, 0): 8, (0, 0, 0, 0, 1, 1): 4, (1, 0, 0, 1, 0, 1): 6, (0, 1, 0, 0, 0, 1): 4,
               (0, 0, 0, 1, 0, 0): 2, (1, 1, 0, 1, 0, 1): 8, (0, 1, 0, 1, 0, 0): 4, (0, 0, 0, 0, 0, 1): 2,
               (1, 1, 0, 0, 0, 0): 4, (0, 1, 0, 0, 1, 1): 6, (1, 0, 0, 0, 1, 0): 4, (1, 1, 0, 1, 1, 1): 10,
               (0, 1, 0, 1, 1, 1): 8, (0, 0, 0, 0, 1, 0): 2, (1, 1, 0, 0, 1, 1): 9, (0, 1, 0, 0, 0, 0): 2,
               (1, 0, 0, 0, 0, 1): 3, (1, 1, 0, 1, 0, 0): 7, (0, 1, 0, 1, 0, 1): 5, (0, 1, 1, 1, 1, 0): 9,
               (1, 1, 0, 0, 0, 1): 6, (1, 0, 1, 1, 1, 0): 8, (1, 1, 1, 1, 1, 1): 11, (1, 0, 0, 0, 1, 1): 6}


def LoadSirmsDict():
    f = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'short_sirms_dict.json'))
    tmp = json.load(f)
    f.close()
    return (tmp)


##    f = open(os.path.join(os.path.dirname(sys.argv[0]), 'short_sirms_dict.txt'), 'tr')
##    sirms_dict = {}
##    for line in f:
##        tmp = line.strip().split('\t')
##        sirms_dict[tmp[0]] = tmp[1]
##    f.close()
##    return(sirms_dict)

def det(a):
    return (a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2])
            - a[1][0] * (a[0][1] * a[2][2] - a[2][1] * a[0][2])
            + a[2][0] * (a[0][1] * a[1][2] - a[1][1] * a[0][2]))


def GetStereoRL(atom_elements, coord):
    atom_numbers = [elements.index(el) for el in atom_elements]
    # sort by atomic number, highest atom is last, lowest atom is first
    atom_numbers, coord = SortTwoLists(atom_numbers, coord)
    # calc diffefence between coord of the lowest atom and all other atoms
    b = [[j - coord[0][i] for i, j in enumerate(elm)] for elm in coord[1:]]
    # calc the signed volume of parallelepiped
    d = det(b)
    if d > 0: return ('R')
    if d < 0: return ('L')
    return (None)


def GetNormalCoord(c):
    # calc normal vector to the plane represented by three points (list of lists)
    # http://wiki.mirgames.ru/%D0%BD%D0%BE%D1%80%D0%BC%D0%B0%D0%BB%D1%8C#naxozhdenie_vektora_normali_k_ploskosti_postroennoj_po_trem_tochkam
    A = c[0][1] * (c[1][2] - c[2][2]) + c[1][1] * (c[2][2] - c[0][2]) + c[2][1] * (c[0][2] - c[1][2])
    B = c[0][2] * (c[1][0] - c[2][0]) + c[1][2] * (c[2][0] - c[0][0]) + c[2][2] * (c[0][0] - c[1][0])
    C = c[0][0] * (c[1][1] - c[2][1]) + c[1][0] * (c[2][1] - c[0][1]) + c[2][0] * (c[0][1] - c[1][1])
    return (A, B, C)


def GetDihedralAngle(coord1, coord2):
    # return dihedral angle in dergees
    A1, B1, C1 = GetNormalCoord(coord1)
    A2, B2, C2 = GetNormalCoord(coord2)
    return (degrees(
        acos((A1 * A2 + B1 * B2 + C1 * C2) / sqrt(A1 ** 2 + B1 ** 2 + C1 ** 2) / sqrt(A2 ** 2 + B2 ** 2 + C2 ** 2))))


def GetStereoZE(coord):
    # coord - coordinates of 4 consequentially bonded atoms
    # http://www.mathsisfun.com/geometry/dihedral-angles.html
    dihedral_angle = GetDihedralAngle(coord[0:3], coord[1:4])
    if 175 <= dihedral_angle <= 180:
        return ("E")
    if 0 <= dihedral_angle <= 5:
        return ("Z")
    return (None)


def GetStereoP(coord):
    # checks whether simplex is planar
    # coord - sorted coordinates of vertexes in order [1,1,2,4] (vertex degree)
    dihedral_angle = GetDihedralAngle(coord[0:2] + [coord[3]], coord[1:4])
    if (0 <= dihedral_angle <= 5) or (175 <= dihedral_angle <= 180):
        return (True)
    else:
        return (False)


def GetCanonNameByDict(labels, bonds, sirms_dict):
    uniq_labels = sorted(set(labels))
    dumb_labels = [dumb_label_set[uniq_labels.index(a)] for a in labels]
    s_name = ','.join(dumb_labels) + '|' + ','.join(map(str, bonds))
    canon_bonds = sirms_dict[s_name]
    return (','.join(sorted(labels)) + '|' + canon_bonds)


def GetSirmsType(bonds):
    ##    return(bin_s_types[tuple(b > 0 for b in bonds)])
    s_vertexes = [0, 0, 0, 0]
    for i, b in enumerate(combinations([0, 1, 2, 3], 2)):
        if bonds[i] > 0:
            s_vertexes[b[0]] += 1
            s_vertexes[b[1]] += 1
    return (s_types[tuple(sorted(s_vertexes))])


def GetSirmsType2(mol, atoms):
    b = []
    # atoms = set(atoms)
    for a in atoms:
        # b.append(len(set(mol.bonds[a].keys()).intersection(atoms)))
        tmp = 0
        for a1 in mol.bonds[a].keys():
            if a1 in atoms:
                tmp += 1
        b.append(tmp)
    idx = tuple(sorted(b))
    return (s_types[idx])


def GenCanonName(labels, bonds, a):
    # form simplex as dict
    s = {a_id: {"label": label} for label, a_id in zip(labels, a)}
    for i, b in enumerate(combinations(a, 2)):
        s[b[0]][b[1]] = s[b[1]][b[0]] = bonds[i]
    # get simplex type
    s_type = GetSirmsType(bonds)
    # create dict with labels and corresponding ids
    d = defaultdict(list)
    for k, v in s.items():
        d[v["label"]].append(k)
    # generate all combinations of ids of each label and combine generated sequences in specific order
    max_bond_value = [-1, -1, -1, -1, -1, -1]  # list which will store the maximum found value of bonds
    p = {k: permutations(v) for k, v in d.items()}
    for k in product(*map(p.get, sorted(p))):
        x = [s[i[0]][i[1]] for i in combinations(list(chain.from_iterable(k)), 2)]
        if x > max_bond_value:
            max_bond_value = x
    # return canonic name
    return (','.join(sorted([v["label"] for v in s.values()])) + '|' + ','.join(map(str, max_bond_value))) + '|' + str(
        s_type)
