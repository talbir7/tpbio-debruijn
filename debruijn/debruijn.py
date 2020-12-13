#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
#import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

# Création du graphe de de Bruijn

# 1) Identification des k-mer uniques

# Fonction qui prend un seul argument correspondant au fichier fastq et retourne un générateur de séquences

def read_fastq(fastq_f):

    with open(fastq_f, "rt") as f:
        for line in f:
        # Pour chaque ligne du fichier fastq, sauter la ligne en cours (commençant par "@"), retourner la ligne
        # suivante (contenant les séquences), et enfin passer les deux lignes suivantes (commençant par "+" et "J" respectivement)
            yield next(f).replace("\n", "")
            next(f)
            next(f)

# Fonction qui prend une séquence, une taille de k-mer et retourne un générateur de k-mer

def cut_kmer(read, kmer_size):

    count = len(read) #nombre de caractères de la séquence
    i = 0 #compteur
    s = kmer_size #variable qui conditionne la continuité de la boucle
    while(s <= count):
        kmer = read[i:i+kmer_size]
        i += 1
        s = i + kmer_size
        yield(kmer)

# Fonction qui prend un fichier fastq, une taille k- mer et retourne un dictionnaire ayant pour clé le k-mer et pour valeur le nombre d’occurrence de ce k-mer

def build_kmer_dict(fastq_file, kmer_size):

    # récupérer tous les k-mers (avec redondance) et les stocker dans une liste
    list_kmers = []
    for line in read_fastq(fastq_file):
        for kmer in cut_kmer(line, kmer_size):
            list_kmers.append(kmer)

    # construire le dictionnaire à partir de la liste contenant tous les k-mers
    dict = {kmer_key:list_kmers.count(kmer_key) for kmer_key in list_kmers}

    return dict

# Fonction qui prend en entrée un dictionnaire de k-mer et crée l’arbre de k-mers préfixes et suffixes

def build_graph(kmer_dict):

    G = nx.DiGraph()
    for kmer_key in kmer_dict.keys():
        # enlever le dernier caractère du k-mer
        prefix = kmer_key[:-1]
        # enlever le première caractère du k-mer
        suffix = kmer_key[1:]
        # itérer sur les k-mers et construire les arcs du graphe
        G.add_edge(prefix, suffix, weight= kmer_dict[kmer_key])
    return G

# Simplification du graphe de de Bruijn

# Résolution des bulles

# Fonction qui prend un graphe et une liste de chemin et retourne un graphe nettoyé des chemins indésirables

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):

    for path in path_list:
        for node in path:
            if path.index(node) == len(path)-1:
                if delete_sink_node:
                    graph.remove_node(node)
                else:
                    continue
            elif path.index(node) == 0:
                if delete_entry_node:
                    graph.remove_node(node)
                else:
                    continue
            else:
                graph.remove_node(node)
    return graph

# Fonction qui qui prend une liste de valeur, qui retourne l’écart type

def std(data):

    ecart_type = statistics.stdev(data)
    return ecart_type

# Fonction qui retourne un graphe nettoyé des chemins indésirables

def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):

    weight_idx = [path for path, weight in enumerate(weight_avg_list) if weight == max(weight_avg_list)]

    final_paths = []

    if len(weight_avg_list) > 1:
        final_paths = [path_list[i] for i in weight_idx]

    if len(path_length) > 1:
        path_length = [path_length[i] for i in weight_idx]
        idx_length = [path for path, length in enumerate(path_length) if length == max(path_length)]
        final_paths = [final_paths[i] for i in idx_length]

    if len(final_paths) > 1:
        random_index = random.randint(0, len(final_paths))
        final_paths = final_paths[random_index]
    path_list.remove(final_paths[0])

    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)

    return graph

# Fonction qui prend un graphe et un chemin et qui retourne un poids moyen

def path_average_weight(graph, path):

    s = sum([graph[path[i]][path[i+1]]['weight'] for i in range(len(path)-1)])

    return s/(len(path)-1)

# Fonction qui prend un graphe, un nœud ancêtre, un nœud descendant et retourne un graph nettoyé de la bulle se trouvant entre ces deux nœuds

def solve_bubble(graph, ancestor_node, descendant_node):

    paths = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    lengths, avg_weights = ([], [])

    for path in paths:

        lengths.append(len(path))
        avg_weights.append(path_average_weight(graph, path))

    graph = select_best_path(graph, paths, lengths, avg_weights)

    return graph

# Fonction qui prend un graphe et retourne un graphe sans bulle

def simplify_bubbles(graph):

    ancestors, descendant = ([], [])

    for node in graph.nodes:

        descs = [desc for desc in graph.predecessors(node)]
        succs = [succ for succ in graph.successors(node)]

        if len(descs) >= 2:
            descendant.append(node)

        if len(succs) >= 2:
            ancestors.append(node)

    for a, d in zip(ancestors, descendant):
        graph = solve_bubble(graph, a, d)

    return graph

# Détection des pointes

# Fonction qui prend un graphe et une liste de noeuds d’entrée et retourne graphe sans chemin d’entrée indésirable

def solve_entry_tips(graph, starting_nodes):

    graph = simplify_bubbles(graph)
    path_list, path_length, path_avg_weight = ([], [], [])
    for node in starting_nodes:
         descs = list(nx.descendants(graph, node))
         for d in descs:
             preds = list(graph.predecessors(d))
             if len(preds) >= 2:
                 paths = list(nx.all_simple_paths(graph, node, d))
                 for p in paths:
                     path_list.append(p)
                     path_length.append(len(p))
                     path_avg_weight.append(path_average_weight(graph, p))
    graph = select_best_path(graph, path_list, path_length, path_avg_weight,
                            delete_entry_node=True, delete_sink_node=False)

    return graph

# Fonction qui prend un graphe et une liste de noeuds de sortie et retourne graphe sans chemin de sortie indésirable

def solve_out_tips(graph, ending_nodes):

    graph = simplify_bubbles(graph)
    path_list, path_avg_weight, path_length = ([], [], [])
    final_b = []
    for node in ending_nodes:
        all_nodes = list(graph.nodes)
        for i in range(len(all_nodes)):
            succs = list(graph.successors(all_nodes[i]))
            if len(succs) > 1:
                s = all_nodes[i]
                final_b.append([s, node])
    for b in final_b:
        for path in nx.all_simple_paths(graph, source=b[0], target=b[1]):
            path_list.append(path)
            path_avg_weight.append(path_average_weight(graph, path))
            path_length.append(len(path))
    graph = select_best_path(graph, path_list, path_length, path_avg_weight,
                            delete_entry_node=False, delete_sink_node=True)

    return graph

# 2) Parcours du graphe de de Bruijn

# Fonction qui prend en entrée un graphe et retourne une liste de noeuds d’entrée

def get_starting_nodes(graph):

    # Si le degré d'entrée est égal à 0 alors le noeud représente un noeud d'entrée
    list_starting_nodes = [n for n in graph.nodes() if graph.in_degree(n)==0]

    return list_starting_nodes

# Fonction qui prend en entrée un graphe et retourne une liste de noeuds de sortie

def get_sink_nodes(graph):

    # Si le degré de sortie est égal à 0 alors le noeud représente un noeud de sortie
    list_sink_nodes = [n for n in graph.nodes() if not graph.out_degree(n)]

    return list_sink_nodes

# Fonction qui prend un graphe, une liste de noeuds d’entrée et une liste de sortie et retourne une liste de tuple(contig, taille du contig)

def get_contigs(graph, starting_nodes, ending_nodes):

    contigs = []
    all_length = dict(nx.all_pairs_shortest_path_length(graph))
    for starting_node in starting_nodes:
        for ending_node in ending_nodes:
            path = ''
            for paths in nx.all_simple_paths(graph, source=starting_node, target=ending_node):
                contig  = [v for k, v in enumerate(paths,1) if k%2!=0]
                for p in contig:
                    path  += p
            contigs.append( (path, all_length[starting_node][ending_node]+2) )

    return contigs


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

# Fonction qui qui prend une liste de tuple (contig, taille du contig) et un nom de fichier de sortie et écrit un fichier de sortie contenant les contigs selon le format fasta

def save_contigs(contigs_list, output_file):

    with open(output_file, "w+") as f:
        count = 0
        for i in range(0, len(contigs_list)):
            fill_ = fill(contigs_list[i][0])
            count += 1
            f.write(f">contig_" + str(count) + "len=" + str(contigs_list[i][1]) + "\n" )
            f.write(fill_)
            f.write("\n")

        f.close()

def draw_graph(graph, graphimg_file):

    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)



#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Lecture du fichier et construction du graphe
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    # Trouver les noeuds d'entrée et les noeuds de sortie
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    # Résolution des bulles
    graph = simplify_bubbles(graph)
    # Résolution des pointes d’entrée et de sortie
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)
    # Ecriture du/des contigs
    list_contigs = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(list_contigs, args.output_file)
    # Dessiner le graphe
    draw_graph(graph, 'mygraph')

if __name__ == '__main__':
    main()
