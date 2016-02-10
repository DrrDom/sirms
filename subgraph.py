#!/usr/bin/env python
# author          : Pavel
# date            : 08.02.16
# version         : 0.1
# python_version  : 3
# copyright       : Pavel 2016
# license         : GPL3
#==============================================================================

import itertools


# def _GenConnectedSubgraphs(mol, not_visited, min_num_atoms=4, max_num_atoms=4, curr_subset=set(), neighbors=set(), res=[]):
#
#     if min_num_atoms <= len(curr_subset) <= max_num_atoms:
#         res.append(curr_subset)
#
#     if not curr_subset:
#         candidates = set(not_visited)
#     else:
#         candidates = not_visited.intersection(neighbors)
#
#     if candidates and len(curr_subset) < max_num_atoms:
#         for a in candidates:
#             not_visited.remove(a)
#             tmp1 = set(not_visited)
#             tmp2 = set(curr_subset)
#             tmp2.add(a)
#             tmp3 = not_visited.intersection(mol["bonds"][a])
#             tmp3 = neighbors.union(tmp3)
#             _GenConnectedSubgraphs(mol, tmp1, min_num_atoms=min_num_atoms, max_num_atoms=max_num_atoms, curr_subset=tmp2, neighbors=tmp3, res=res)


def GetNumberConnectedComponents(mol, atoms):

    res = 0
    visited = [False] * len(atoms)

    def dfs(i):
        visited[i] = True
        for a in mol.bonds[atoms[i]].keys():
            if a in atoms and not visited[atoms.index(a)]:
                dfs(atoms.index(a))

    for i in range(len(atoms)):
        if not visited[i]:
            dfs(i)
            res += 1

    return res


def GenConnectedSubgraphs(mol, not_visited, min_num_atoms=4, max_num_atoms=4, curr_subset=set(), neighbors=set(), res=[]):

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
            tmp3 = not_visited.intersection(mol.bonds[a].keys())
            tmp3 = neighbors.union(tmp3)
            GenConnectedSubgraphs(mol, tmp1, min_num_atoms=min_num_atoms, max_num_atoms=max_num_atoms, curr_subset=tmp2, neighbors=tmp3, res=res)


def GetAllNeighbours(mol, atoms):
    output = set(atoms)
    for a in atoms:
        output = output.union(mol.bonds[a].keys())
    return output


def GetAtomsCombinations(mol, min_num_components=1, max_num_components=2, min_num_atoms=4, max_num_atoms=4, noH=False):

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
        atoms = set(a for a in mol.atoms.keys() if mol.atoms[a]["label"] != 'H')
    else:
        atoms = set(mol.atoms.keys())

    # storage of results
    res = []
    # if only single component fragments are looking for then there is no need to search for fragments smaller than nim_num_atoms
    if max_num_components == 1:
        GenConnectedSubgraphs(mol, atoms, min_num_atoms=min_num_atoms, max_num_atoms=max_num_atoms, res=res)
    else:
        GenConnectedSubgraphs(mol, atoms, min_num_atoms=1, max_num_atoms=max_num_atoms, res=res)

    # get neighbours
    nb = [GetAllNeighbours(mol, v) for v in res]

    for n in range(min_num_components, max_num_components + 1):
        for comb in itertools.combinations(range(len(res)), n):
            # if min_num_atoms <= sum(len(res[i]) for i in comb) <= max_num_atoms and (len(comb) == 1 or not set.intersection(*[nb[i] for i in comb])):
            if min_num_atoms <= sum(len(res[i]) for i in comb) <= max_num_atoms and (len(comb) == 1 or not CheckIntersection(res, nb, comb)):
                yield tuple(set.union(*[res[i] for i in comb]))


