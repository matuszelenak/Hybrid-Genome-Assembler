"""
==================
Spectral Embedding
==================

The spectral layout positions the nodes of the graph based on the
eigenvectors of the graph Laplacian $L = D - A$, where $A$ is the
adjacency matrix and $D$ is the degree matrix of the graph.
By default, the spectral layout will embed the graph in two
dimensions (you can embed your graph in other dimensions using the
``dim`` argument to either :func:`~drawing.nx_pylab.draw_spectral` or
:func:`~drawing.layout.spectral_layout`).

When the edges of the graph represent similarity between the incident
nodes, the spectral embedding will place highly similar nodes closer
to one another than nodes which are less similar.

This is particularly striking when you spectrally embed a grid
graph.  In the full grid graph, the nodes in the center of the
graph are pulled apart more than nodes on the periphery.
As you remove internal nodes, this effect increases.
"""

import matplotlib.pyplot as plt
import networkx as nx


options = {
    'node_color': 'C0',
    'node_size': 250,
    'font_color': "white"
}
G = nx.Graph()
G.add_weighted_edges_from([
    (0, 1, 2),
    (0, 4, 3),
    (1, 2, 2),
    (2, 3, 1.5),
    (3, 4, 1),
    (0, 1, 2),
    (4, 5, 0.3),
    (5, 6, 2),
    (6, 7, 1.5),
    (7, 8, 2.5),
    (8, 5, 2),
    (6, 9, 0.2),
    (7, 9, 0.3),
    (9, 11, 2),
    (9, 10, 3),
    (10, 11, 2.6)
])

edges = G.edges()
#
plt.subplot(121)
layout = nx.spring_layout(G, scale=2)
nx.draw(G, layout, **options, edges=edges, width=[G[u][v]['weight'] for u, v in edges], with_labels=True)
nx.draw_networkx_edge_labels(G, layout, edge_labels={(u, v): int(G[u][v]['weight'] * 100) for u, v in edges})

plt.subplot(122)
layout = nx.spectral_layout(G, scale=2)
nx.draw_networkx_nodes(G, layout, **options)
nx.draw_networkx_labels(G, layout, **options)

plt.show()
