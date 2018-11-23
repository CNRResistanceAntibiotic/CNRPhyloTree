#!/usr/bin/python3
import os

import numpy as np
import pandas as pd
import networkx as nx
import argparse
import json
from bokeh.io import show, output_file, save
from bokeh.plotting import figure
from bokeh.palettes import RdYlGn as RYG
from bokeh.palettes import Greys
from bokeh.palettes import Viridis
from bokeh.models import Plot, ColumnDataSource, Circle, MultiLine, Label, LabelSet, Legend, LegendItem


def remove_duplicate(names, array):
    dl_dic = {}
    df = pd.DataFrame(data=array, index=names, columns=names)
    for n, name in enumerate(names):
        n = n + 1
        cols = []
        for item in names:
            if item != name:
                cols.append(item)
        sd = df.loc[cols, name]
        sd = list(sd[sd == 0].index)
        if sd != []:
            dl_dic[name] = [name] + sd
            dl_dic[name] = list(set(dl_dic[name]))
            dl_dic[name].sort()
    # print(dl_dic)
    indexes = []
    for items in dl_dic.values():
        for item in items[1:]:
            if item not in indexes:
                indexes.append(item)
                df.drop(item, axis=0, inplace=True)
                df.drop(item, axis=1, inplace=True)
                df.rename(index={items[0]: items[0] + '\n' + item}, columns={items[0]: items[0] + '\n' + item},
                          inplace=True)
    print("\nRemove duplicates:", indexes)
    print(dl_dic)
    print('')
    arr = df.as_matrix()
    names = list(df.index)
    return names, arr


def load_matrix(mtx_file):
    """
    load the matrix of the file in a dictionary
    which associate strain name to all other strain name with the distance between them
    :param mtx_file: The file path of the matrix
    :return:
    """
    # matrix as numpy array
    arr = np.genfromtxt(mtx_file, delimiter='\t', dtype=str)
    names = arr[0, 1:].tolist()
    data = np.asarray(arr[1:, 1:], dtype=np.int)
    print(names, "\n\n", data)
    names, data = remove_duplicate(names, data)
    print(names, "\n\n", data)
    #    data as scipy sparce matrix
    #    data = csr_matrix(data)
    #    store column names, row names and data in mtx_dict dictionnary
    mtx_dict = {'names': names, 'matrix': data}
    return mtx_dict


def compute_network(mtx_dic, graph_name, st_file):
    """
    Compute a network of strain with a distance matrix of SNP difference and
     plot him in interactive format with html extension
    :param mtx_dic: This is the distance matrix which contain the distance between each strain in number of different SNP
    :param graph_name: The name and place of the graph
    :param st_file: The path of a configuration file that associate strain with their ST number
    :return:
    """

    e = []
    #    Create an empty graph
    G = nx.Graph()
    G.add_nodes_from(mtx_dic['names'])

    matrix = mtx_dic['matrix'].tolist()

    if st_file:
        conf = open("../configSTnumber.json", "r")
        config_dico = json.load(conf)
        config_st = config_dico['dico']

        # attribute ST number to strain which have ST number precise in csv configuration file, attribute also name
        for node in G.nodes:
            try:
                print(node, ": ", config_st[node])
                G.nodes[node]['st'] = config_st[node]
            except KeyError:
                G.nodes[node]['st'] = ""
                continue
    else:
        for node in G.nodes:
            G.nodes[node]['st'] = ""

    # create edges with the weight in attribute
    n = 0
    for parent in mtx_dic['names']:
        m = n + 1
        while m < len(mtx_dic['names']):
            if matrix[n][m] < 10:
                e.append((parent, mtx_dic['names'][m],
                          {'weight': matrix[n][m], "color": 'red', 'width': 8, 'alpha': 0.8}))
            elif matrix[n][m] < 25:
                e.append((parent, mtx_dic['names'][m],
                          {'weight': matrix[n][m], "color": 'orange', 'width': 8, 'alpha': 0.8}))
            elif matrix[n][m] < 50:
                e.append((parent, mtx_dic['names'][m],
                          {'weight': matrix[n][m], "color": 'yellow', 'width': 6, 'alpha': 0.8}))
            elif matrix[n][m] < 100:
                e.append((parent, mtx_dic['names'][m],
                          {'weight': matrix[n][m], "color": 'grey', 'width': 4, 'alpha': 0.8}))
            else:
                e.append((parent, mtx_dic['names'][m],
                          {'weight': matrix[n][m], "color": 'grey', 'width': 1, 'alpha': 0.8}))

            m += 1
        n += 1
        G.add_edges_from(e)

    # with weightless edge of the graph G create a tree with all node
    g = nx.minimum_spanning_tree(G)

    # take the position give by the Kamada-Kawai algorithm and assign them to corresponding strain
    pos = nx.kamada_kawai_layout(g)
    nx.set_node_attributes(g, pos, 'pos')
    # print(pos)

    # with the coordinate xy of pos attribute the coordinate x and y to each strain
    for i in pos:
        # xode[i] = pos[i][0]
        # yode[i] = pos[i][1]
        g.nodes[i]['x'] = pos[i][0]
        g.nodes[i]['y'] = pos[i][1]

    for start, end in g.edges():
        # print(g.nodes[start]['x'])
        g.edges[start, end]['xector'] = np.array([g.nodes[start]['x'], g.nodes[end]['x']])
        g.edges[start, end]['yector'] = np.array([g.nodes[start]['y'], g.nodes[end]['y']])
    # print(g.nodes.data('st'))

    # attribute a color in function of SNP number with the strain which are the closest
    for node in g.nodes():
        # print(g[node])
        for key in g[node]:
            if node == "Reference":
                g.nodes[node]['color'] = 'green'

            else:
                if g[node][key]['color'] == 'red':
                    g.nodes[node]['color'] = 'red'

                elif g[node][key]['color'] == 'orange':
                    try:
                        if g.nodes[node]['color'] != 'red':
                            g.nodes[node]['color'] = 'orange'
                    except KeyError:
                        g.nodes[node]['color'] = 'orange'

                elif g[node][key]['color'] == 'yellow':
                    try:
                        if g.nodes[node]['color'] != 'red' and g.nodes[node]['color'] != 'orange':
                            # print(g[node][key])
                            g.nodes[node]['color'] = 'yellow'
                    except KeyError:
                        # print(g[node][key])
                        g.nodes[node]['color'] = 'yellow'

                else:
                    try:
                        if g.nodes[node]['color'] != 'red' and g.nodes[node]['color'] != 'orange' \
                                and g.nodes[node]['color'] != 'yellow':
                            g.nodes[node]['color'] = 'grey'
                    except KeyError:
                        g.nodes[node]['color'] = 'grey'

    # print(g.nodes[node], " nd \n")

    bokeh_node_color = color_converter(nx.get_node_attributes(g, 'color'))
    # print(bokeh_node_color)
    nx.set_node_attributes(g, bokeh_node_color, 'bokeh_color')
    bokeh_edge_color = color_converter(nx.get_edge_attributes(g, 'color'))
    # print(bokeh_edge_color)
    nx.set_edge_attributes(g, bokeh_edge_color, 'bokeh_color')

    # create a dictionary with the st number of each strain
    node_st = nx.get_node_attributes(g, "st")
    dic_node_st = {}
    for node, st in node_st.items():
        if st:
            dic_node_st[node] = node + "\n" + st
        else:
            dic_node_st[node] = node

    # -------
    bokeh_edge = []
    for a, z, e in g.edges(data=True):
        e['start'] = a
        e['end'] = z
        bokeh_edge.append(e)
    # print(bokeh_edge)
    # print(pd.DataFrame.from_dict({k: v for k, v in g.nodes(data=True)}, orient='index'))

    source_node = ColumnDataSource(pd.DataFrame.from_dict({k: v for k, v in g.nodes(data=True)}, orient='index'))
    source_edge = ColumnDataSource(pd.DataFrame.from_dict(bokeh_edge))

    xedge = []
    yedge = []
    snp = []
    for edge in g.edges():
        xedge.append((g.node[edge[0]]['pos'][0] + g.node[edge[1]]['pos'][0]) / 2)
        yedge.append((g.node[edge[0]]['pos'][1] + g.node[edge[1]]['pos'][1]) / 2)
        snp.append(g.edges[edge]['weight'])
    source = ColumnDataSource({'x': [i for i in xedge],
                               'y': [i for i in yedge],
                               'dist': [i for i in snp]})

    plot = figure(plot_width=900, plot_height=900, x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
    try:
        plot.title.text = "Strain network in function of SNP difference \n Schema MLST :" + dic_node_st['MLST_schema']
    except KeyError:
        plot.title.text = "Strain network in function of SNP difference "

    plot.multi_line(xs='xector', ys='yector', line_color='bokeh_color',
                    line_width='width', line_alpha='alpha', source=source_edge)

    # option select legend
    legend_node = []
    for name, data in g.nodes(data=True):
        cercle = plot.circle(x=data['x'], y=data['y'], size=20, fill_color=data['bokeh_color'],
                             muted_color=Viridis[10][3])
        legend_node.append((name, [cercle]))

    strain_label = plot.text(x='x', y='y', text_baseline='middle', text_align='center', text='index',
                             muted_alpha=0, source=source_node)

    st_label = plot.text(x='x', y='y', text_baseline='bottom', text_align='center',
                         y_offset=-10, text='st',
                         muted_alpha=0, source=source_node)

    snp_label = plot.text(x='x', y='y', text_baseline='bottom', text_align='center', y_offset=10,
                          text='dist', text_font='times', text_font_style='italic', text_font_size='10pt',
                          muted_alpha=0, source=source)

    legend = Legend(
        items=[
                  ('Strain Name', [strain_label]),
                  ('ST Number', [st_label]),
                  ('SNP difference', [snp_label])
              ] + legend_node,
        location=(10, 0)
    )

    legend.click_policy = 'mute'
    plot.add_layout(legend, 'right')
    output_file(graph_name)
    # -------
    save(plot)


def color_converter(colors):
    """
    It convert the name of color in a format readable by the module bokeh
    :param colors: it is a dictionary of color associate with indices of nodes or edges or ...
    :return: A dictionary of color in 'bokeh format' associate with indices of nodes, edges or ...
    """
    print(type(colors), "\n\n")
    for color in colors:
        if colors[color] == 'green':
            colors[color] = RYG[10][1]
        elif colors[color] == 'red':
            colors[color] = RYG[9][8]
        elif colors[color] == 'orange':
            colors[color] = RYG[9][7]
        elif colors[color] == 'yellow':
            colors[color] = Viridis[3][2]
        elif colors[color] == 'grey':
            colors[color] = Greys[6][2]
    return colors


def pre_main(args):
    mtx_file = os.path.abspath(mtx_file)
    graph_name = args.graphName
    st_file = os.path.abspath(args.st_file)
    main(mtx_file, graph_name, st_file)


def main(mtx_file, graph_name, st_file):
    mtx_dic = load_matrix(mtx_file)
    compute_network(mtx_dic, graph_name, st_file)


def run():
    """
    Take the argument of the command and give a variable for this
    :return: args
    """

    parser = argparse.ArgumentParser(description="""This program convert a distance matrix in graphic representation of SNP's difference between strain.
    if you have the reference genome in your matrix please make it name \"Reference\" (it will give a different colour for him)""")
    parser.add_argument('-i', '--inputMatrix', dest='mtxFile', default='',
                        help="Enter you distance matrix with his path")
    parser.add_argument('-o', '--outputGraph', dest='graphName', default='',
                        help="Enter the directory with the name of the graph file that you want")
    parser.add_argument('-s', '--ST', dest='st_file', default='', help="""Enter the path of the file which contain 
    the ST number define in with his strain """)
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
