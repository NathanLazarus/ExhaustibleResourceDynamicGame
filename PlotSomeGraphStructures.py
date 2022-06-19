import networkx as nx
from matplotlib import pyplot as plt
import itertools
from collections import Counter
import time
import random
from networkx.utils import pairwise
import numpy as np

N = 3
max_val = 3
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
    'lightblue'
]

subset_color = subset_color * 5 #(N * max_val) // len(subset_color) + 2
subset_color = list(itertools.chain.from_iterable(itertools.repeat(x, -(-(N * max_val) // len(subset_color))) for x in subset_color))
# subset_color.extend(['lightblue', 'lightblue', 'lightblue', 'lightblue', 'lightblue'])
# -(-(N * max_val) // len(subset_color)) = ceil(N * max_val) / len(subset_color)

print(len(all_nodes))

graf = nx.DiGraph()
for (ind, val) in enumerate(all_nodes):
    graf.add_nodes_from([val], layer = layers[ind])

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

plt.tight_layout()
pos = nx.multipartite_layout(graf, subset_key='layer')
color = [subset_color[data['layer']] for v, data in graf.nodes(data=True)]
# color = [i for i in ['darkorange', 'cadetblue', 'violet'] * 2 + ['gold', 'chartreuse', 'maroon'] * 2 + ['indigo', 'olive', 'bisque'] * 2 for _ in range(2)]
# color = ['cadetblue'] * len(all_nodes)
# nx.draw(graf, pos, with_labels = True, node_color = color)
nx.draw(graf, pos, with_labels = False, node_color = color)
plt.axis('equal')
plt.show()
# plt.savefig('small2p.png', format='PNG')
plt.clf()