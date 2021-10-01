import networkx as nx
import itertools
from collections import Counter
import time
import random

N = 4
max_val = 4
n_possible_vals = max_val + 1

all_nodes = list(itertools.product(range(n_possible_vals), repeat = N))
max_vals = [max_val] * N

layers = [sum(x) for x in all_nodes]



print(len(all_nodes))

graf = nx.DiGraph()

graf.add_nodes_from(all_nodes)
start = time.time()
for tup in all_nodes:
    for i in range(len(tup)):
        if tup[i] <= (max_vals[i] - 1):
            child = list(tup)
            child[i] += 1
            graf.add_edges_from([(tup, tuple(child))])

print('elapsed', time.time() - start)




for node in graf.nodes:
    if (sum(node) == (N * max_val) // 2) & (sum([x == 0 for x in node]) == 0):
        mytrystartnode = node
        break

print(mytrystartnode)
mytrystartnode = (4, 1, 2, 1)

midlayer2 = []
for node in graf.nodes:
    if sum(node) == ((N * max_val) // 2):
        midlayer2.append(node)

# randomSS = 1609
# random_improved = []
# going_to_be_children = random.sample(midlayer2, randomSS // N)
# for child in going_to_be_children:
#     random_improved.extend(list(graf.predecessors(child)))
# while len(random_improved) < randomSS:
#     going_to_be_a_child = random.sample(midlayer2, 1)
#     for single_child in going_to_be_a_child:
#         print(single_child)
#         random_improved.extend(list(graf.predecessors(single_child)))


# random_improvedChildren = []
# for trying in iter(random_improved):
#     random_improvedChildren.extend(list(graf.successors(trying)))
# countDict_random_improved = Counter(random_improvedChildren)
# for node in countDict_random_improved:
#     countDict_random_improved[node] += node.count(0)
# print('random_improved')
# print(Counter(countDict_random_improved.values()))
# print(len(random_improved))

siblings = set()
cousins = set()
mytry = set()
nodepreds = graf.predecessors(mytrystartnode)

for nodepred in nodepreds:
    sibs = graf.successors(nodepred)
    siblings.update(list(sibs))


for sib in iter(siblings):
    sibpreds = graf.predecessors(sib)
    for sibpred in sibpreds:
        cousins.update(list(graf.successors(sibpred)))


for cousin in iter(cousins):
    mytry.update(list(graf.predecessors(cousin)))


print(len(mytry))

todo = [2] * 1
for n_to_add_per in todo:
    mytryChildren = []
    for trying in iter(mytry):
        mytryChildren.extend(list(graf.successors(trying)))
    countDict = Counter(mytryChildren)
    for node in countDict:
        countDict[node] += node.count(0)
    print(Counter(countDict.values()))
    for node, n in countDict.items():
        print(node, n)
        if (n >= (N - n_to_add_per)):
            for x in graf.predecessors(node):
                mytry.add(x)
    print(len(mytry))



mytryChildren = []
for trying in iter(mytry):
    mytryChildren.extend(list(graf.successors(trying)))
countDict = Counter(mytryChildren)
for node in countDict:
    countDict[node] += node.count(0)
print(Counter(countDict.values()))
gotten = set()
for node, n in countDict.items():
    if (n == N):
        gotten.add(node)
print('try1', len(mytry), len(gotten))




# for node in graf.nodes:
#     if (sum(node) == (N * max_val) // 2) & (node not in gotten):
#         newstartnode = node
#         break

# print(f'{newstartnode = }')


# siblings = set()
# cousins = set()
# mytry = set()
# nodepreds = graf.predecessors(newstartnode)

# for nodepred in nodepreds:
#     sibs = graf.successors(nodepred)
#     siblings.update(list(sibs))


# for sib in iter(siblings):
#     sibpreds = graf.predecessors(sib)
#     for sibpred in sibpreds:
#         cousins.update(list(graf.successors(sibpred)))


# for cousin in iter(cousins):
#     if cousin not in gotten:
#         mytry.update(list(graf.predecessors(cousin)))

# print(len(mytry))

# todo = [2] * 1
# for n_to_add_per in todo:
#     mytryChildren = []
#     for trying in iter(mytry):
#         mytryChildren.extend(list(graf.successors(trying)))
#     countDict = Counter(mytryChildren)
#     for node in countDict:
#         countDict[node] += node.count(0)
#     print(Counter(countDict.values()))
#     for node, n in countDict.items():
#         if n == 4:
#             print(node not in gotten)
#         if (n >= (N - n_to_add_per)) & (n < N) & (node not in gotten):
#             for x in graf.predecessors(node):
#                 mytry.add(x)
#     print(len(mytry))



# mytryChildren = []
# for trying in iter(mytry):
#     mytryChildren.extend(list(graf.successors(trying)))
# countDict = Counter(mytryChildren)
# for node in countDict:
#     countDict[node] += node.count(0)
# print(Counter(countDict.values()))
# print('try2', len(mytry))

# gotten2 = set()
# for node, n in countDict.items():
#     if (n == N):
#         gotten2.add(node)

# print('try1', len(gotten))
# print('try2', len(mytry), len(gotten2))

# toget = set()
# for node in graf.nodes:
#     if (sum(node) == (N * max_val) // 2):
#         toget.add(node)

# print(len(toget))
# print(len(set.union(gotten, gotten2)))

# class batch:
#     """docstring for ClassName"""
#     def __init__(self, graph_parent):
#         self.graph_parent = graph_parent
#         self.firstlayer = set()
#         self.secondlayer = set()

#     def add_to_first_layer(self, nodes_to_add):
#         self.firstlayer.update(list(nodes_to_add))
#         self.secondlayer.update(list(is_gotten(children(nodes_to_add))))

#     def add_to_second_layer(self, nodes_to_add):
#         self.secondlayer.update(list(nodes_to_add))
#         self.firstlayer.update(list(self.graph_parent.predecessors(nodes_to_add)))


#     def add_cousins(self, start_node):
#         self.arg = arg
        

# firstlayerofbatches = []
# secondlayerofbatches = []
# min_batch_size = 50
# # think more about the first and last layers where there are no parents/children
# while len(gotten) < len(toget):
#     notgotten = set.difference(toget, set(secondlayerofbatches))
#     thisbatch = batch(graf)
#     thisbatch.add_cousins(next(iter(notgotten)))
#     while thisbatch.firstlayer.len < min_batch_size:
#         n_to_add_per = 1
#         old_batch_len = thisbatch.secondlayer.len
#         while old_batch_len == thisbatch.secondlayer.len:
#             thisbatch.widen(gotten, n_to_add_per)
#             if thisbatch.secondlayer.len == len(toget):
#                 break
#             n_to_add_per += 1
#             if n_to_add_per > N:
#                 thisbatch.add_layer2(next(iter(set.difference(toget, thisbatch.secondlayer))))
#     firstlayerofbatches.append(thisbatch.firstlayer)
#     secondlayerofbatches.append(thisbatch.secondlayer)



