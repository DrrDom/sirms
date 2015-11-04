#-------------------------------------------------------------------------------
# Name:        mols
# Purpose:     Mol class to operate with molecules
#
# Author:      Pavel Polishchuk
#
# Created:     11.01.2013
# Copyright:   (c) Pavel Polishchuk 2013-2015
# Licence:     GPLv3
#-------------------------------------------------------------------------------

from itertools import combinations

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


class Mol3:
    # double bonds can be:
    # 2 - non-steric or undefined double bond (2, 0)
    # 21 - Z double bond (2, 1)
    # 22 - E double bond (2, 2)
    # 23 - cyclic double bond in rings of size up to 7 (such bonds are always cis) (2, 3)

    def __init__(self):
        self.atoms = {}
        self.bonds = dict()
        self.title = ""
        self.stereo = False

    def AddAtom(self, id, label, x, y, z, formal_charge):
        self.atoms[id] = {'label': label, 'x': x, 'y': y, 'z': z,
                          'property': {'elm': {'label': label, 'value': label, 'weight': elements.index(label)},
                                       'none': {'label': 'A', 'value': 'A', 'weight': 0}},
                          'formal_charge': formal_charge}
        self.bonds[id] = dict()

    def AddBond(self, id1, id2, bond_type):
        if id1 not in self.bonds.keys():
            self.bonds[id1] = dict()
        if id2 not in self.bonds.keys():
            self.bonds[id2] = dict()
        # bond value is a tuple: 1 - bond order, 2 - double bond stereo type (0 - unspecified, 1 - cis, 2 -trans, 3 - cyclic)
        self.bonds[id1][id2] = self.bonds[id2][id1] = (bond_type, 0)

    def GetBondOrder(self, id1, id2):
        return self.bonds.get(id1, dict()).get(id2, (0, 0))[0]

    def GetBondType(self, id1, id2):
        bond_order = self.GetBondOrder(id1, id2)
        if self.stereo and bond_order == 2 and self.bonds[id1][id2][1] != 0:
            return self.bonds[id1][id2][0] * 10 + self.bonds[id1][id2][1]
        else:
            return bond_order

    def SetDoubleBondConfig(self, id1, id2, bond_stereo):
        if bond_stereo not in [0, 1, 2, 3]:
            raise Exception('Wrong double bond stereo!')
        self.bonds[id1][id2] = self.bonds[id2][id1] = (2, bond_stereo)

    def _Path(self, start_atom, list_atom, cycles_local, visited, size_range):
        for a in self.bonds[start_atom].keys():
            if len(list_atom) <= max(size_range) and a not in list_atom and a not in visited:
                self._Path(a, list_atom + [a], cycles_local, visited, size_range)
            elif len(list_atom) in size_range and a == list_atom[0]:
                if tuple(set(sorted(list_atom))) not in cycles_local:
                    cycles_local.add(tuple(set(sorted(list_atom))))

    def GetCycles(self, min_size, max_size):
        cycles = set()
        visited = set()  # atoms which already has been tested as a cycle member can be excluded from further cycles checks
        for a in sorted(self.atoms.keys()):
            visited.add(a)
            if len(self.bonds[a].keys()) > 1:
                self._Path(a, [a], cycles, visited, range(min_size, max_size + 1))
        return cycles

    def SetCyclicDoubleBondsCis(self, min_size=3, max_size=7):
        cycles = self.GetCycles(min_size, max_size)
        for cycle in cycles:
            for a1, a2 in combinations(cycle, 2):
                if self.GetBondOrder(a1, a2) == 2:
                    self.bonds[a1][a2] = self.bonds[a2][a1] = (2, 3)
