import json
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import networkx as nx
import aRST.script.initialization as storage
from os.path import join


def plot_rxnnetwork(file=None):
    custom_colors = {
        'r': '#CC071E',
        'g': '#57ab27',
        'b': '#00549f',
        'lv': '#9E5C9B',
        't': '#0098a1',
        'o': '#f6a800',
        "lo": "#FFD27F",
        'y': '#ffed00',
        "dy": "#D9B600",
        'dgrey': '#646567',
        'lgrey': '#cfd1d2',
        'black': '#000000',
        'white': '#ffffff',
    }
    if not file:
        from aRST.script.exploration.network import rxnnetwork
        graph = rxnnetwork.graph
    else:
        # Load the reaction data
        with open(file, 'r') as f:
            data = json.load(f)

        # Create a directed graph
        graph = nx.DiGraph()

        # Add nodes with their attributes
        for node_data in data['nodes']:
            node = node_data['node']
            info = node_data['info']
            graph.add_node(node, info=info)

        # Add edges
        graph.add_edges_from(data['edges'])

    fig, ax = plt.subplots(dpi=200)

    node_labels = {node: graph.nodes[node].get("info", None).get('scanwd', '0') for node in graph.nodes()}
    E_rxn = [graph.nodes[node].get("info", None)["E_rxn"] for node in graph.nodes]


    norm = mcolors.Normalize(vmin=-50, vmax=300)
    cmap = plt.cm.get_cmap("RdYlBu")
    node_colors = [cmap(norm(e)) for e in E_rxn]

    pos = nx.spring_layout(graph, k=1, iterations=50)
    nx.draw_networkx_nodes(
        graph,
        pos,
        node_color=node_colors,
        node_size=100,
        edgecolors='none',
        linewidths=1
    )
    nx.draw_networkx_edges(graph,
                           pos,
                           arrowstyle='->',
                           arrowsize=10,
                           width=1.5,
                           edge_color=custom_colors["dgrey"],
                           connectionstyle="arc3,rad=0.2")



    label_pos = {node: (x, y + 0.07) for node, (x, y) in pos.items()}
    nx.draw_networkx_labels(graph, label_pos, labels=node_labels, font_size=8, font_weight='bold')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, label="E_rxn (kcal/mol)")


    plt.tight_layout()
    plt.savefig(join(storage.setting.wd0,"rxnnetwork"))
