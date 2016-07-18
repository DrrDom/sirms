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
from collections import Counter
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

    def __getWeininger(self, sub, labels):
        """
        modified morgan algorithm
        init with prime numbers as weights according to labels
        iteratively update weight by sum up with weight of neighbours multiplied by bond order
        return: dict with atoms hash
        """
        a = sorted(set(labels))
        init_weights = [a.index(lab) * 100000 for lab in labels]
        for i, atom in enumerate(sub):
            bond_orders = [self.GetBondOrder(atom, nei) for nei in self.bonds[atom].keys() if nei in sub]
            c = Counter(bond_orders)
            init_weights[i] = init_weights[i] + len(bond_orders) * 10000 + c[4] * 1000 + c[3] * 100 + c[2] * 10 + c[1]
        # init_weights - label_rank, number_neighbours, number_bonds_1, number_bonds_2, number_bonds_3, number_bond_4 (aromatic)
        a = sorted(set(init_weights))
        ranks = [a.index(w) for w in init_weights]

        # -1 to be always TRUE on the first iteration
        previous_ranks_len = len(set(ranks)) - 1

        while len(set(ranks)) < len(ranks) and previous_ranks_len < len(set(ranks)):

            previous_ranks_len = len(set(ranks))

            primes = [self.primes[r] for r in ranks]

            primes_product = []
            for atom in sub:
                p = 1
                for nei in self.bonds[atom].keys():
                    try:
                        if nei in sub:
                            p *= primes[sub.index(nei)]
                    except ValueError:
                        continue
                primes_product.append(p)

            # re-rank
            new_ranks = [None] * len(sub)
            # d = defaultdict(list)
            d = dict()
            for i, r in enumerate(ranks):
                if r not in d.keys():
                    d[r] = []
                d[r].append(i)
            curr_rank = 0
            for i in sorted(d.keys()):
                if len(d[i]) == 1:
                    new_ranks[d[i][0]] = curr_rank
                    curr_rank += 1
                else:
                    pp = [primes_product[j] for j in d[i]]
                    a = sorted(set(pp))
                    local_ranks = [a.index(p) for p in pp]
                    for loc_r, j in zip(local_ranks, d[i]):
                        new_ranks[j] = curr_rank + loc_r
                    curr_rank += len(a)

            ranks = list(new_ranks)

        return {atom: rank for atom, rank in zip(sub, ranks)}

    def get_name(self, sub, labels):

        # if {10, 11, 12}.issubset(sub):
        #     # if 'H' in labels:
        #     if set(sub).intersection([24, 20]):
        #         print(sub)

        self.nextnumb = self.__numb()
        self.sub = set(sub)
        self.levels = self.__getWeininger(sub, labels)
        inter = min(self.levels, key=self.levels.get)

        res = [self.__getSmiles([inter], inter, dict(zip(sub, labels)))]
        visited = res[0][0]
        while visited != self.sub:
            remained = {k: self.levels[k] for k in self.sub.difference(visited)}
            inter = min(remained, key=remained.get)
            res.append(self.__getSmiles([inter], inter, dict(zip(sub, labels))))
            visited = visited.union(res[-1][0])

        # to get canonical multi-component SMILES
        res = sorted([r[1] for r in res])

        return '.'.join(res)

    def __numb(self):
        i = 1
        while True:
            yield i
            i += 1

