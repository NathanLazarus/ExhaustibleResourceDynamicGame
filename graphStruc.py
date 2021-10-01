import networkx as nx
from matplotlib import pyplot as plt
import itertools
from collections import Counter
import time
import random
from networkx.utils import pairwise

N = 4
max_val = 4
n_possible_vals = max_val + 1

all_nodes = list(itertools.product(range(n_possible_vals), repeat = N))
max_vals = [max_val] * N
# all_nodes = list(itertools.product([0,1,2,3,4], [0,1,2,3], [0,1,2]))
# all_nodes = list(itertools.product([0,1], [0,1], [0,1]))
# max_vals = [1, 1, 1]

layers = [sum(x) for x in all_nodes]



# subset_sizes = [5, 5, 4, 3, 2, 4, 4, 3]
subset_color = [
    'gold',
    'violet',
    'limegreen',
    'darkorange',
    'blue'
]

# subset_color = subset_color * 5 #(N * max_val) // len(subset_color) + 2
subset_color = list(itertools.chain.from_iterable(itertools.repeat(x, -(-(N * max_val) // len(subset_color))) for x in subset_color))
subset_color.extend(['blue', 'blue', 'blue', 'blue', 'blue'])
# -(-(N * max_val) // len(subset_color)) = ceil(N * max_val) / len(subset_color)

print(len(all_nodes))

graf = nx.DiGraph()
# for (ind, val) in enumerate(all_nodes):
    # graf.add_nodes_from([val]) #, layer = layers[ind])

graf.add_nodes_from(all_nodes)
start = time.time()
# edges = []
for tup in all_nodes:
    for i in range(len(tup)):
        if tup[i] <= (max_vals[i] - 1):
            # zeros = [0, 0, 0]
            # zeros[i] = 1
            # as_list = list(tup)
            # print(tup, tuple([x + y for x, y in zip(zeros, list(tup))]))
            child = list(tup)
            child[i] += 1
            # edges.append((tup, tuple(child)))
            graf.add_edges_from([(tup, tuple(child))])
# edges = [x for tup in all_nodes for i in range(N) if tup[i <= (max_val - 1)] ]
# graf.add_edges_from(edges)
print('elapsed', time.time() - start)



# intlist = [0] * len(all_nodes)
# for index, node in enumerate(all_nodes):
#     for ind, val in enumerate(node):
#         intlist[index] += val * n_possible_vals ** ind


# ints_to_tuples = {}
# for ind, val in enumerate(intlist):
#     ints_to_tuples[val] = all_nodes[ind]
# graf.add_nodes_from(intlist)
# print(time.time())
# start = time.time()
# for node in intlist:
#     for i in range(N):
#         if node // n_possible_vals ** i % n_possible_vals + 1 < n_possible_vals:
#             asdf = node + n_possible_vals ** i
#             graf.add_edges_from([(node, node + n_possible_vals ** i)])
# print(time.time() - start)






# # this takes forever to run
# for node1, node1layer in graf.nodes(data = True):
#     for node2, node2layer in graf.nodes(data = True):
#         if (node1layer['layer'] == node2layer['layer']) & (node1 != node2):
#             node1pred = set(graf.predecessors(node1))
#             node2pred = set(graf.predecessors(node2))
#             asdf = set.union(node1pred, node2pred)
#             if len(node1pred) + len(node2pred) - len(asdf) > 0:
#                 1+1#print(node1, node2, node1pred, node2pred)

# print(nx.shortest_path(graf))
# print(nx.single_source_shortest_path(graf, (0, 0, 0,)))
# print(nx.single_source_shortest_path_length(graf, (0, 0, 0,)))
# print(nx.ancestors(graf, (1, 1, 0)))
# print(nx.ancestors(graf, (1, 1, 1)))

# solved = set()
# for i in range(30):
#     newlySolved = set()
#     for node in graf.nodes:
#         ancestors = nx.ancestors(graf, node)
#         # if node == (1,1,1,1):
#         #     print(set.difference(ancestors, solved))
#         #     print(len(set.difference(ancestors, solved)))
#         # if node == (1,1,1,0):
#         #     print(set.difference(ancestors, solved))
#         #     print(len(set.difference(ancestors, solved)))
#         if len(set.difference(ancestors, solved)) == 0:
#             # print(node, len(set.difference(ancestors, solved)))
#             # print(f'{solved = }')
#             # print(node)
#             newlySolved.add(node)

#     solved = set.union(solved, newlySolved)
#     # print(newlySolved)
#     # print(i, len(solved))

# undirected = graf.to_undirected()

# for i in range(2):
    # for thislayer in range(5):
# nodeset = set()
# nodeset.update(graf.nodes)
# bucketlist = []
# while len(nodeset) > 0:
#     node = next(iter(nodeset))
#     siblings = set()
#     cousins = set()
#     nodepreds = graf.predecessors(node)
#     for nodepred in nodepreds:
#         sibs = graf.successors(nodepred)
#         for sib in sibs:
#             siblings.add(sib)
#             sibpreds = graf.predecessors(sib)
#             for sibpred in sibpreds:
#                 cous = graf.successors(sibpred)
#                 cousins.update(list(cous))
#                 # totry = set()
#                 # for cousin in cous:
#                 #     totry.update(list(graf.predecessors(cousin)))
#                 #     for trying in iter(totry):


#     # print(node, siblings, cousins)
#     inbucket = set.union(siblings, cousins)
#     inbucket.add(node)
#     bucketlist.append(inbucket)
#     nodeset = set.difference(nodeset, inbucket)


# print(bucketlist)
# print(len(bucketlist))
# print([len(x) for x in bucketlist])
# print([x for x in bucketlist if len(x) == 1])

# print(graf.nodes())
# print(graf.edges())

# plt.tight_layout()
# pos = nx.multipartite_layout(graf, subset_key='layer')
# color = [subset_color[data['layer']] for v, data in graf.nodes(data=True)]
# nx.draw(graf, pos, with_labels = True, node_color = color)
# plt.axis('equal')
# plt.show()
# # plt.savefig('graf.png', format='PNG')
# plt.clf()


for node in graf.nodes:
    if (sum(node) == (N * max_val) // 2) & (sum([x == 0 for x in node]) == 0):
        mytrystartnode = node
        break

# midlayer = []
# for node in graf.nodes:
#     if sum(node) == ((N * max_val) // 2) - 1:
#         midlayer.append(node)

# # node = max(intlist) // 2
# import random
# random3019 = random.sample(midlayer, 3019)
# random3019Children = []
# for trying in iter(random3019):
#     random3019Children.extend(list(graf.successors(trying)))
# countDict_random = Counter(random3019Children)
# for node in countDict_random:
#     countDict_random[node] += node.count(0)
# print('random')
# print(Counter(countDict_random.values()))
# print(len(random3019))


midlayer2 = []
for node in graf.nodes:
    if sum(node) == ((N * max_val) // 2):
        midlayer2.append(node)

randomSS = 1
# node = max(intlist) // 2
random_improved = []
going_to_be_children = random.sample(midlayer2, randomSS // N)
for child in going_to_be_children:
    random_improved.extend(list(graf.predecessors(child)))
while len(random_improved) < randomSS:
    going_to_be_a_child = random.sample(midlayer2, 1)
    for single_child in going_to_be_a_child:
        print(single_child)
        random_improved.extend(list(graf.predecessors(single_child)))


random_improvedChildren = []
for trying in iter(random_improved):
    random_improvedChildren.extend(list(graf.successors(trying)))
countDict_random_improved = Counter(random_improvedChildren)
for node in countDict_random_improved:
    countDict_random_improved[node] += node.count(0)
print('random_improved')
print(Counter(countDict_random_improved.values()))
print(len(random_improved))

# node = next(iter(nodeset))
siblings = set()
cousins = set()
mytry = set()
nodepreds = graf.predecessors(mytrystartnode)

# print(time.time())
for nodepred in nodepreds:
    sibs = graf.successors(nodepred)
    siblings.update(list(sibs))

# print(time.time())

for sib in iter(siblings):
    sibpreds = graf.predecessors(sib)
    for sibpred in sibpreds:
        cousins.update(list(graf.successors(sibpred)))

# print(time.time())

for cousin in iter(cousins):
    mytry.update(list(graf.predecessors(cousin)))


# mytryChildren = []
# for trying in iter(mytry):
#     mytryChildren.extend(list(graf.successors(trying)))

# countDict = Counter(mytryChildren)
# for node in countDict:
#     countDict[node] += node.count(0)
# print(Counter(countDict.values()))
# print(len(mytry))

# for node, n in countDict.items():
#     if (n == (N - 2)) | (n == (N - 1)):
#         print(node, n)
#         for x in graf.predecessors(node):
#             mytry.add(x)
#         # mytry.update(list(graf.predecessors(node)))

print(len(mytry))

# def grab_lowhanging_fruit(tryset, n_to_add_per):
todo = [2] * 3
# todo.append(1)
for n_to_add_per in todo:
    mytryChildren = []
    for trying in iter(mytry):
        mytryChildren.extend(list(graf.successors(trying)))
    countDict = Counter(mytryChildren)
    for node in countDict:
        countDict[node] += node.count(0)
    print(Counter(countDict.values()))
    for node, n in countDict.items():
        if (n >= (N - n_to_add_per)) | (n < N):
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
print(len(mytry))
gotten = set()
for node, n in countDict.items():
    if (n == N):
        gotten.add(node)




for node in graf.nodes:
    if (sum(node) == (N * max_val) // 2) & (sum([x == 0 for x in node]) == 0) & (node not in gotten):
        newstartnode = node
        break




siblings = set()
cousins = set()
mytry = set()
nodepreds = graf.predecessors(newstartnode)

for nodepred in nodepreds:
    sibs = graf.successors(nodepred)
    siblings.update(list(sibs))


for sib in iter(siblings):
    sibpreds = graf.predecessors(sib)
    for sibpred in sibpreds:
        cousins.update(list(graf.successors(sibpred)))


for cousin in iter(cousins):
    if cousin not in gotten:
        mytry.update(list(graf.predecessors(cousin)))




todo = [2] * 5
for n_to_add_per in todo:
    mytryChildren = []
    for trying in iter(mytry):
        mytryChildren.extend(list(graf.successors(trying)))
    countDict = Counter(mytryChildren)
    for node in countDict:
        countDict[node] += node.count(0)
    print(Counter(countDict.values()))
    for node, n in countDict.items():
        if (n >= (N - n_to_add_per)) | (n < N) & (node not in gotten):
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
print(len(mytry))

print(len(set(mytryChildren)))
print(len(gotten))