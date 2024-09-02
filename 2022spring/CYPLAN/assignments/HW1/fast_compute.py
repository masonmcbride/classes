import numpy as np
import networkx as nx
import pandas as pd

fl = pd.read_csv("internet.edgelist.txt", delimiter='\t')
print("read file fl")

g = nx.Graph()
for index, row in fl.iterrows():
    g.add_edge(row['0'], row['1'])
largest_cc = max(nx.connected_components(g), key=len)
G = g.subgraph(largest_cc)
print("created graph G")

print("computing shortest paths on G")
out = nx.shortest_path(G)
print("computed shortest paths on G")
print(out)
