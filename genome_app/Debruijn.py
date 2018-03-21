#!/usr/bin/env python
import DNA_operations as RevComp
import networkx as nx
import matplotlib.pyplot as plt


with open('rosalind_dbru.txt') as input_data:
    k_mers = [line.strip() for line in input_data.readlines()]

# Get the edge elements.
DBG_edge_elmts = set()
for kmer in k_mers:
    DBG_edge_elmts.add(kmer)
    DBG_edge_elmts.add(RevComp.ReverseComplementDNA(kmer))

    
print(k_mers)

# Create the edges.

k = len(k_mers[0])
edge = lambda elmt: '('+elmt[0:k-1]+', '+elmt[1:k]+')'
DBG_edges = [edge(elmt) for elmt in DBG_edge_elmts]

# Write and save the adjacency list.
l='\n'.join(DBG_edges)
print(l)
map = str.maketrans('', '', '(),')

DBG_edges=[ element.translate(map) for element in DBG_edges]
with open('DBRU.txt', 'w') as output_file:
    output_file.write('\n'.join(DBG_edges))

ins = open( "DBRU.txt", "r" )
data = [tuple(str(n) for n in line.split()) for line in ins]

G = nx.MultiDiGraph(data)

pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos,  node_color = 'r',node_size=1000)
nx.draw_networkx_edges(G, pos,  arrows=True)
nx.draw_networkx_labels(G,pos)
plt.show()
