# import networkx as nx
from itertools import chain, product, accumulate
from collections import Counter
from numpy import array_split
import numpy as np
# from math import ceil
# import time
# import random
# from itertools import chain

# N = 2
# max_val = 9
# max_vals = [max_val] * N
max_vals = [3,3,3]
# max_vals = [4,3,2]


def mysplit(ary, max_chunk_size):

    Ntotal = len(ary)
    Nsections = -(-Ntotal // max_chunk_size) #ceil
    Neach_section, extras = divmod(Ntotal, Nsections)
    section_sizes = ([0] +
                     extras * [Neach_section + 1] +
                     (Nsections - extras) * [Neach_section])
    div_points = list(accumulate(section_sizes))
    sublists = []
    for i in range(Nsections):
        st = div_points[i]
        end = div_points[i + 1]
        sublists.append(ary[st:end])

    return sublists


def testsplit(ary, Nsections):

    Ntotal = len(ary)
    Neach_section, extras = divmod(Ntotal, Nsections)
    print(divmod(Ntotal, Nsections))
    section_sizes = ([0] +
                     extras * [Neach_section+1] +
                     (Nsections-extras) * [Neach_section])
    div_points = np.array(section_sizes).cumsum()

    sub_arys = []
    for i in range(Nsections):
        st = div_points[i]
        end = div_points[i + 1]
        sub_arys.append(ary[st:end])

    return sub_arys

test = testsplit(np.arange(22), 16)
print(test)
beensplit = [mysplit(list(range(x + 1)), 2) for x in max_vals]
print(beensplit)
ranges = [[(min(sublist), max(sublist)) for sublist in lists] for lists in beensplit]

def square_gen(boundaries):
    job_id = 0
    for square_boundaries in product(*boundaries):
        square_sides = [list(range(x[0], x[1] + 1)) for x in square_boundaries]
        vertex = tuple([x[0] for x in square_boundaries])
        successors = []
        for i in range(len(ranges)):
            asdf = ranges[i].index(square_boundaries[i])
            if asdf < len(ranges):
                successors.append(ranges[asdf + 1][0])
        yield [job_id, vertex, list(product(*square_sides))]
        job_id += 1

print(ranges)
print(list(square_gen(ranges)))

vertices = [[min(sublist) for sublist in lists] for lists in beensplit]
print(vertices)
breakpoints = [[min(sublist) for sublist in lists] for lists in beensplit]
print(list(product(*vertices)))
import networkx as nx

graf = nx.DiGraph()

graf.add_nodes_from(all_nodes)
start = time.time()
for tup in all_nodes:
    for i in range(len(tup)):
        if tup[i] <= (max_vals[i] - 1):
            child = list(tup)
            child[i] += 1
            graf.add_edges_from([(tup, tuple(child))])

print(graf.successors((0,0,0,0)))
# n_possible_vals = max(max_vals) + 1
# all_nodes = list(product(range(n_possible_vals), repeat = N))
# print(len(all_nodes))

# layers = [sum(x) for x in all_nodes]



def get_predecessors(tup):
    parents = []
    for i in range(len(tup)):
        if tup[i] >= 1:
            parent = list(tup)
            parent[i] -= 1
            parents.append(tuple(parent))

    return parents

def get_successors(tup):
    children = []
    for i in range(len(tup)):
        if tup[i] <= (max_vals[i] - 1):
            child = list(tup)
            child[i] += 1
            children.append(tuple(child))

    return children


# graf = nx.DiGraph()

# graf.add_nodes_from(all_nodes)
# start = time.time()
# for tup in all_nodes:
#     for i in range(len(tup)):
#         if tup[i] <= (max_vals[i] - 1):
#             child = list(tup)
#             child[i] += 1
#             graf.add_edges_from([(tup, tuple(child))])

# print('elapsed', time.time() - start)

# toget_firstlayer = set()
# for node in all_nodes:
#     if (sum(node) == (N * max_val) // 2 - 1):
#         toget_firstlayer.add(node)

# toget_secondlayer = set()
# for node in all_nodes:
#     if (sum(node) == (N * max_val) // 2):
#         toget_secondlayer.add(node)




def all_perms_summing_to_x(length, total_sum, max_vals):
    if length == 1:
        yield (total_sum,)
    else:
        range_min = max(0, total_sum - sum(max_vals[len(max_vals) - length + 1:]))
        range_max = min(total_sum + 1, max_vals[len(max_vals) - length] + 1)
        # print(length, total_sum, max_vals, len(max_vals) - length + 1)
        # print(sum(max_vals[length - 1:]))
        # print(f'{range_min = }')
        # print(f'{range_max = }')
        for value in range(range_min, range_max):
            for permutation in all_perms_summing_to_x(length - 1, total_sum - value, max_vals):
                yield (value,) + permutation

# L = list(all_perms_summing_to_x(N, 14, max_val))
# print('total permutations:',len(L))
# otherway = set()
# for node in all_nodes:
#     if sum(node) == 14:
#         otherway.add(node)
# print(len(otherway))

toget_firstlayer = set(all_perms_summing_to_x(N, sum(max_vals) // 2 - 1, max_vals))
toget_secondlayer = set(all_perms_summing_to_x(N, sum(max_vals) // 2, max_vals))

# print(f'{toget_firstlayer = }')
# print(f'{toget_secondlayer = }')

print(len(toget_firstlayer))
print(len(toget_secondlayer))



class Batch:
    """This is a batch of jobs for HTCondor"""
    def __init__(self, max_vals, gotten_in_other_batches = set()):
        # self.graph_parent = graph_parent
        self.N = len(max_vals)
        self.firstlayer = set()
        self.secondlayer = set()
        self.max_vals = max_vals
        self.gotten_in_other_batches = gotten_in_other_batches #should be pointer

    def add_to_first_layer(self, nodes_to_add):
        firstlayer_additions = set(nodes_to_add).difference(self.firstlayer)
        self.firstlayer.update(firstlayer_additions)
        for node in firstlayer_additions:
            self.secondlayer.update(self.is_gotten(get_successors(node)).difference(self.gotten_in_other_batches))

    def add_to_second_layer(self, nodes_to_add):
        secondlayer_additions = set(nodes_to_add).difference(self.secondlayer).difference(self.gotten_in_other_batches)
        self.secondlayer.update(secondlayer_additions)
        for node in secondlayer_additions:
            self.firstlayer.update(get_predecessors(node))

    def get_siblings(self, start_node):
        sibs = set()
        for nodepred in get_predecessors(start_node):
            sibs.update(get_successors(nodepred))
        return sibs

    def add_cousins(self, start_node):
        for sib in self.get_siblings(start_node):
            self.add_to_second_layer(self.get_siblings(sib))


    # @staticmethod
    # def is_gotten(nodes_to_check_if_gotten, graph_parent, nodes_in_previous_layer, N):
    #     pass

    def is_gotten(self, nodes_to_check_if_gotten):
        gotten = set()
        for node in nodes_to_check_if_gotten:
            if not set(get_predecessors(node)).difference(self.firstlayer):
                gotten.add(node)
        return gotten

    def widen(self, n_to_add_per):
        firstlayerChildren = []
        for firstlayer_element in iter(self.firstlayer):
            firstlayerChildren.extend(list(get_successors(firstlayer_element)))
        countDict = Counter(firstlayerChildren)
        for node in countDict:
            countDict[node] += node.count(0) #each 0 means 1 fewer parent
        to_add = set()
        for node, n in countDict.items():
            if (n >= (self.N - n_to_add_per)) & (node not in self.gotten_in_other_batches):
                to_add.add(node)
        self.add_to_second_layer(to_add)

    
firstlayerofbatches = []
secondlayerofbatches = []
typical_batch_size = 1000
minimum_batch_size = 200


# think more about the first and last layers where there are no parents/children
all_done = 0
while True:
    gotten = set(chain.from_iterable(secondlayerofbatches))
    if len(gotten) == len(toget_secondlayer):
        break
    notgotten = toget_secondlayer.difference(gotten)
    
    thisbatch = Batch(max_vals, gotten)
    thisbatch.add_cousins(next(iter(notgotten)))
    while len(thisbatch.firstlayer) < typical_batch_size:
        n_to_add_per = 0
        old_batch_len = len(thisbatch.secondlayer)
        while old_batch_len == len(thisbatch.secondlayer):
            thisbatch.widen(n_to_add_per) #, gotten)
            # if len(gotten) == len(toget): #thisbatch.secondlayer.len
            #     break
            n_to_add_per += 1
            if n_to_add_per > N:
                diffs = set.difference(notgotten, thisbatch.secondlayer)
                if len(diffs) == 0:
                    all_done = 1
                    break
                else:
                    thisbatch.add_to_second_layer([next(iter(diffs))])
                    break
        if all_done == 1:
            break
    else: #once the batch size reaches the typical size, check if there are only a few remaining
        diffs = set.difference(notgotten, thisbatch.secondlayer)
        if len(diffs) < minimum_batch_size:
            thisbatch.add_to_second_layer(diffs)

    firstlayerofbatches.append(thisbatch.firstlayer)
    secondlayerofbatches.append(thisbatch.secondlayer)

    if all_done == 1:
        break




print([len(x) for x in secondlayerofbatches])

layer1lens = [len(x) for x in firstlayerofbatches]
print(layer1lens)
print(sum(layer1lens))
print(len(toget_firstlayer))






# jkl = set([1,2,3])
# print(jkl)
# asdf = set()
# asdf.update(jkl)
# print(asdf)

# x = Batch(graf, N)
# # print(len(x.firstlayer))
# # print(x.firstlayer)
# asdf = next(iter(toget_secondlayer))
# asdf = (4, 1, 2, 1)
# x.add_cousins(asdf)
# print(len(x.firstlayer))
# print(len(x.secondlayer))
# # x.widen(2)
# # x.widen(0)

# while len(x.secondlayer) < 55:
#     n_to_add_per = 0
#     old_batch_len = len(x.secondlayer)
#     while old_batch_len == len(x.secondlayer):
#         x.widen(n_to_add_per)
#         n_to_add_per += 1
#         if n_to_add_per > N:
#             x.add_to_second_layer(next(iter(set.difference(toget, x.secondlayer))))
#         print(n_to_add_per, old_batch_len, len(x.secondlayer))

# x.widen(0)


# print('**********************')
# print((1,4,2,1) in x.secondlayer)
# print(len(x.firstlayer))
# print(len(x.secondlayer))
# print(len(x.is_gotten(list())))
# print(len(toget_firstlayer))
# print(len(toget_secondlayer))