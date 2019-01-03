#!/usr/bin/python3
import os
from csv import DictReader

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
from networkx.drawing.nx_agraph import graphviz_layout


def remove_duplicate(names, array):
    dl_dic = {}
    legends = {}
    df = pd.DataFrame(data=array, index=names, columns=names)
    for n, name in enumerate(names):
        n = n + 1
        cols = []
        for item in names:
            if item != name:
                cols.append(item)
        sd = df.loc[cols, name]
        sd = list(sd[sd == 0].index)
        if sd:
            dl_dic[name] = [name] + sd
            dl_dic[name] = list(set(dl_dic[name]))
            dl_dic[name].sort()
        else:
            legends[name] = name

    indexes = []
    for items in dl_dic.values():
        for item in items[1:]:
            if item not in indexes:
                indexes.append(item)
                df.drop(item, axis=0, inplace=True)
                df.drop(item, axis=1, inplace=True)
                items_names = ""
                if len(items) == 2:
                    items_names = items[0]+","+items[1]
                else:
                    items_names = items[0]+","+items[1]+"...."

                df.rename(index={items[0]: items_names}, columns={items[0]: items_names},
                          inplace=True)
                legends[items_names] = ",".join(items)

    print("\nRemove duplicates: {0}\n".format(indexes))
    # print(dl_dic)
    arr = df.as_matrix()
    names = list(df.index)
    return names, arr, legends


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
    # print(names, "\n\n", data)
    names, data, legends = remove_duplicate(names, data)
    # print("\n")
    # print(names, "\n\n", data)
    #    data as scipy sparce matrix
    #    data = csr_matrix(data)
    #    store column names, row names and data in mtx_dict dictionnary
    mtx_dict = {'names': names, 'matrix': data}

    return mtx_dict, legends


def add_ST_label(config_file, legends, G, second_mlst):

    if config_file:
        config_list = load_config(config_file)

        # print(G.nodes['output_snippy_CNR2218_ecoli_cgMLST_ref_2018_v1_ecoli'])
        for name_node, value in legends.items():
            for row_dict in config_list:

                if 'MLST-2' in row_dict:
                    second_mlst = True

                if row_dict['strains'] == name_node:

                    # attribute ST number to strain which have ST number precise in csv configuration file, attribute also name
                    if name_node in G.nodes:

                        try:
                            G.nodes[name_node]['st-1'] = row_dict['MLST-1']

                            if second_mlst:
                                G.nodes[name_node]['st-2'] = row_dict['MLST-2']
                        except KeyError:
                            G.nodes[name_node]['st-1'] = ""
                            if second_mlst:
                                G.nodes[name_node]['st-2'] = ""
                            continue
                    else:
                        continue
                    break

                elif row_dict['strains'] in value and "," in name_node:
                    if 'st-1' in G.nodes[name_node] and row_dict['MLST-1'] not in G.nodes[name_node]['st-1'].values():

                        G.nodes[name_node]['st-1'] = "*"
                        if second_mlst:
                            G.nodes[name_node]['st-2'] = "*"
                    else:
                        G.nodes[name_node]['st-1'] = row_dict['MLST-1']
                        if second_mlst:
                            G.nodes[name_node]['st-2'] = row_dict['MLST-2']
                        break

    else:
        for node in G.nodes:
            G.nodes[node]['st-1'] = ""
            if second_mlst:
                G.nodes[node]['st-2'] = ""

    return G, second_mlst


def compute_network(mtx_dic, graph_name, config_file, legends, count_snp_keep):
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

    second_mlst = False

    # add ST Label
    G, second_mlst = add_ST_label(config_file, legends, G, second_mlst)

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
    # pos = nx.kamada_kawai_layout(g, scale=50)
    # pos = nx.spring_layout(g, scale=50)
    pos = nx.nx_agraph.graphviz_layout(g)

    nx.set_node_attributes(g, pos, 'pos')
    print(pos)

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

    x_edge = []
    y_edge = []
    snp = []
    for edge in g.edges():
        x_edge.append((g.node[edge[0]]['pos'][0] + g.node[edge[1]]['pos'][0]) / 2)
        y_edge.append((g.node[edge[0]]['pos'][1] + g.node[edge[1]]['pos'][1]) / 2)
        snp.append(g.edges[edge]['weight'])
    source = ColumnDataSource({'x': [i for i in x_edge],
                               'y': [i for i in y_edge],
                               'dist': [i for i in snp]})

    plot = figure(plot_width=1700, plot_height=900, x_range=(-1000.1, 1000.1), y_range=(-1000.1, 1000.1))
    try:
        plot.title.text = "Strain network in function of SNP difference (Total SNP : {0}) " \
                          "\n Schema MLST :".format(count_snp_keep) + dic_node_st['MLST_schema']
    except KeyError:
        plot.title.text = "Strain network in function of SNP difference (Total SNP : {0}) ".format(count_snp_keep)

    plot.multi_line(xs='xector', ys='yector', line_color='bokeh_color',
                    line_width='width', line_alpha='alpha', source=source_edge)

    # option select legend
    legend_node = []
    for name, data in g.nodes(data=True):
        circle = plot.circle(x=data['x'], y=data['y'], size=20, fill_color=data['bokeh_color'],
                             muted_color=Viridis[10][3])
        for key, value in legends.items():
            if name == key:
                legend_node.append((value, [circle]))
                break
            else:
                pass

    strain_label = plot.text(x='x', y='y', text_baseline='middle', text_align='center', text='index',
                             muted_alpha=0, source=source_node)

    if second_mlst:
        st_label_1 = plot.text(x='x', y='y', text_baseline='bottom', text_align='center',
                               y_offset=-25, text='st-1',
                               muted_alpha=0, source=source_node, text_color='orangered')
        st_label_2 = plot.text(x='x', y='y', text_baseline='bottom', text_align='center',
                         y_offset=-10, text='st-2',
                         muted_alpha=0, source=source_node, text_color='deepskyblue')
    else:
        st_label_1 = plot.text(x='x', y='y', text_baseline='bottom', text_align='center',
                               y_offset=-10, text='st-1',
                               muted_alpha=0, source=source_node, text_color='orangered')

    snp_label = plot.text(x='x', y='y', text_baseline='bottom', text_align='center', y_offset=10,
                          text='dist', text_font='times', text_font_style='italic', text_font_size='10pt',
                          muted_alpha=0, source=source)

    if second_mlst:
        legend = Legend(
            items=[
                      ('Strain Name', [strain_label]),
                      ('ST Number 1', [st_label_1]),
                      ('ST Number 2', [st_label_2]),
                      ('SNP difference', [snp_label])
                  ] + legend_node,
            location=(10, 0)
        )
    else:
        legend = Legend(
            items=[
                      ('Strain Name', [strain_label]),
                      ('ST Number', [st_label_1]),
                      ('SNP difference', [snp_label])
                  ] + legend_node,
            location=(10, 0)
        )

    legend.click_policy = 'mute'
    plot.add_layout(legend, 'right')
    output_file(graph_name)
    # -------
    save(plot)


def load_config(config_file):

    with open(config_file, "r") as conf:
        reader = DictReader(conf, delimiter=",")
        config_list = list(reader)

    return config_list


def color_converter(colors):
    """
    It convert the name of color in a format readable by the module bokeh
    :param colors: it is a dictionary of color associate with indices of nodes or edges or ...
    :return: A dictionary of color in 'bokeh format' associate with indices of nodes, edges or ...
    """

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
    mtx_file = os.path.abspath(args.mtxFile)
    graph_name = args.graphName
    st_file = os.path.abspath(args.st_file)
    main(mtx_file, graph_name, st_file)


def main(mtx_file, graph_name, config_file, count_snp_keep):
    mtx_dic, legends = load_matrix(mtx_file)
    compute_network(mtx_dic, graph_name, config_file, legends, count_snp_keep)


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
