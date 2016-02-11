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
import copy


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
                          'property': {},    # will contain dicts of type 'elm': {'label': ['C'], 'value': 'C'}
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


class Mol4(Mol3):

    def __GetNumberConnectedComponents(self, atoms):

        res = 0
        visited = [False] * len(atoms)

        def dfs(i):
            visited[i] = True
            for a in self.bonds[atoms[i]].keys():
                if a in atoms and not visited[atoms.index(a)]:
                    dfs(atoms.index(a))

        for i in range(len(atoms)):
            if not visited[i]:
                dfs(i)
                res += 1

        return res

    def __GenConnectedSubgraphs(self, not_visited, min_num_atoms=4, max_num_atoms=4, curr_subset=set(), neighbors=set(), res=[]):

        if min_num_atoms <= len(curr_subset) <= max_num_atoms:
            res.append(curr_subset)

        if not curr_subset:
            candidates = set(not_visited)
        else:
            candidates = not_visited.intersection(neighbors)

        if candidates and len(curr_subset) < max_num_atoms:
            for a in candidates:
                not_visited.remove(a)
                tmp1 = set(not_visited)
                tmp2 = set(curr_subset)
                tmp2.add(a)
                tmp3 = not_visited.intersection(self.bonds[a].keys())
                tmp3 = neighbors.union(tmp3)
                self.__GenConnectedSubgraphs(tmp1, min_num_atoms=min_num_atoms, max_num_atoms=max_num_atoms, curr_subset=tmp2, neighbors=tmp3, res=res)

    def __GetAllNeighbours(self, atoms):
        output = set(atoms)
        for a in atoms:
            output = output.union(self.bonds[a].keys())
        return output


    def GetAtomsCombinations(self, min_num_components=1, max_num_components=2, min_num_atoms=4, max_num_atoms=4, noH=False):

        def CheckIntersection(sets, neighbour_sets, ids):
            if set.intersection(*[sets[i] for i in ids]):
                return True
            for i in ids:
                for j in ids:
                    if i != j:
                        if sets[i].intersection(neighbour_sets[j]):
                            return True
            return False

        if noH:
            atoms = set(a for a in self.atoms.keys() if self.atoms[a]["label"] != 'H')
        else:
            atoms = set(self.atoms.keys())

        # storage of results
        res = []
        # if only single component fragments are looking for then there is no need to search for fragments smaller than nim_num_atoms
        if max_num_components == 1:
            self.__GenConnectedSubgraphs(atoms, min_num_atoms=min_num_atoms, max_num_atoms=max_num_atoms, res=res)
        else:
            self.__GenConnectedSubgraphs(atoms, min_num_atoms=1, max_num_atoms=max_num_atoms, res=res)

        # get neighbours
        nb = [self.__GetAllNeighbours(v) for v in res]

        for n in range(min_num_components, max_num_components + 1):
            for comb in combinations(range(len(res)), n):
                # if min_num_atoms <= sum(len(res[i]) for i in comb) <= max_num_atoms and (len(comb) == 1 or not set.intersection(*[nb[i] for i in comb])):
                if min_num_atoms <= sum(len(res[i]) for i in comb) <= max_num_atoms and (len(comb) == 1 or not CheckIntersection(res, nb, comb)):
                    yield tuple(set.union(*[res[i] for i in comb]))


class SmilesMol3(Mol4):

    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
              109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
              233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359,
              367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,
              499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641,
              643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787,
              797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
              947, 953, 967, 971, 977, 983, 991, 997, 1009]

    tosmileskeys = {0: '.', 1: '-', 2: '=', 3: '#', 4: ':', 8: '~'}

    def __tosmiles(self, bond):
        return self.tosmileskeys[bond]

    def __getSmiles(self, trace, inter, labels_dict):
        # trace: atom ids
        # inter: selected atom id
        strace = set(trace)
        iterlist = set(self.bonds[inter].keys()).intersection(self.sub).difference(trace[-2:])

        # solution for unbounded fragments
        # if not iterlist and strace != self.sub:
        #     remained = {k: self.levels[k] for k in self.sub.difference(strace)}
        #     iterlist = {min(remained, key=remained.get)}

        # get atom label
        smi = [labels_dict[inter]]

        concat = []
        stoplist = []
        iterlen = len(iterlist) - 1
        for b, i in enumerate(sorted(list(iterlist), key=self.levels.get)):
            if i in strace:
                if i not in stoplist:
                    # костыль для циклов. чтоб не было 2х проходов.
                    cyc = next(self.nextnumb)
                    concat += [(i, cyc, inter)]
                    smi[0] += '%s%d' % (self.__tosmiles(self.GetBondOrder(inter, i)), cyc)
                if b == iterlen and len(smi) > 3:
                    smi[-1] = smi[-3] = ''
                continue
            deep = self.__getSmiles(copy.copy(trace + [i]), i, labels_dict)
            strace.update(deep[0])

            for j in deep[2]:
                if j[0] == inter:
                    stoplist += [j[2]]
                    smi[0] += '%s%d' % (self.__tosmiles(self.GetBondOrder(inter, j[2])), j[1])
                else:
                    concat.append(j)

            smi += ['(' if iterlen - b else '', '%s' % self.__tosmiles(self.GetBondOrder(inter, i)) + deep[1],
                    ')' if iterlen - b else '']

        return strace, ''.join(smi), concat

    def __getMorgan(self, sub, labels):
        """
        modified morgan algorithm
        init with prime numbers as weights according to labels
        iteratively update weight by sum up with weight of neighbours multiplied by bond order
        return: dict with atoms hash
        """
        a = sorted(set(labels))
        newweights = {x: self.primes[a.index(labels[i])] for i, x in enumerate(sub)}
        self.__buffer = 0
        while self.__accept(newweights) and len(set(newweights.values())) < len(sub):
            weights = copy.copy(newweights)
            for i in sub:
                l = sum([weights[x] * self.GetBondOrder(i, x) / 100 for x in self.bonds[i].keys() if x in sub])
                newweights[i] += l
        return newweights

    def __accept(self, numb):
        """
        stop-function to interrupt Morgan algorithm
        """
        numl = len(set(numb.values()))

        if self.__buffer >= numl:
            return False
        self.__buffer = numl
        return True

    # def __getWeininger(self, sub, labels):
    #     num = 0
    #     a = sorted(set(labels))
    #     newweights = ((x, a.index(labels[i]) * 100 + len(self.bonds[x])) for i, x in enumerate(sub))
    #     while len(set(newweights)) < len(sub):
    #
    #         newweights = {x: self.primes[a.index(labels[i])] for i, x in enumerate(sub)}
    #
    #         num = set()
    #
    #     pass

    def get_name(self, sub, labels):

        # if {10, 11, 12}.issubset(sub):
        #     # if 'H' in labels:
        #     if set(sub).intersection([24, 20]):
        #         print(sub)

        self.nextnumb = self.__numb()
        self.sub = set(sub)
        self.levels = self.__getMorgan(sub, labels)
        inter = min(self.levels, key=self.levels.get)

        res = [self.__getSmiles([inter], inter, dict(zip(sub, labels)))]
        visited = res[0][0]
        while visited != self.sub:
            remained = {k: self.levels[k] for k in self.sub.difference(visited)}
            inter = min(remained, key=remained.get)
            res.append(self.__getSmiles([inter], inter, dict(zip(sub, labels))))
            visited = visited.union(res[-1][0])

        # if {10, 11, 12}.issubset(sub):
        #     # if 'H' in labels:
        #     if set(sub).intersection([24, 20]):
        #         print(sub, '.'.join([i[1] for i in res]))

        return '.'.join([i[1] for i in res])

    def __numb(self):
        i = 1
        while True:
            yield i
            i += 1

