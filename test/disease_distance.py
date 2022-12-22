import networkx as nx
import numpy as np
import optparse
import sys
import os

print(os.getcwd())
os.chdir('data/')
print(os.getcwd())


def print_usage(option, opt, value, parser):
    usage_message = """
# ----------------------------------------------------------------------
    """
    print usage_message
    sys.exit()
    

# =============================================================================
def read_network(network_file):
    G = nx.Graph()
    for line in open(network_file,'r'):
        if line[0]=='#':
            continue
        line_data   = line.strip().split('\t')
        node1 = line_data[0]
        node2 = line_data[1]
        G.add_edge(node1,node2)
    print "\n> done loading network:"
    print "> network contains %s nodes and %s links" %(G.number_of_nodes(),G.number_of_edges())
    return G


# =============================================================================
def read_gene_list(gene_file):
    genes_set = set()
    for line in open(gene_file,'r'):
        if line[0]=='#':
            continue
        line_data = line.strip().split('\t')
        gene      = line_data[0]
        genes_set.add(gene)
    print "\n> done reading genes:"
    print "> %s genes found in %s" %(len(genes_set),gene_file)
    return genes_set


# =============================================================================
def remove_self_links(G):
    sl = G.selfloop_edges()
    G.remove_edges_from(sl)


# =============================================================================
def get_pathlengths_for_single_set(G,given_gene_set):
    # remove all nodes that are not in the network
    all_genes_in_network = set(G.nodes())
    gene_set = given_gene_set & all_genes_in_network
    all_path_lenghts = {}
    # calculate the distance of all possible pairs
    for gene1 in gene_set:
        if not all_path_lenghts.has_key(gene1):
            all_path_lenghts[gene1] = {}
        for gene2 in gene_set:
            if gene1 < gene2:
                try:
                    l = nx.shortest_path_length(G, source=gene1, target=gene2)
                    all_path_lenghts[gene1][gene2] = l
                except:
                    continue
    return all_path_lenghts



# =============================================================================
def get_pathlengths_for_two_sets(G,given_gene_set1,given_gene_set2):
    all_genes_in_network = set(G.nodes())
    gene_set1 = given_gene_set1 & all_genes_in_network
    gene_set2 = given_gene_set2 & all_genes_in_network
    all_path_lenghts = {}
    # calculate the distance of all possible pairs
    for gene1 in gene_set1:
        if not all_path_lenghts.has_key(gene1):
            all_path_lenghts[gene1] = {}
        for gene2 in gene_set2:
            if gene1 != gene2:
                try:
                    l = nx.shortest_path_length(G, source=gene1, target=gene2)
                    if gene1 < gene2:
                        all_path_lenghts[gene1][gene2] = l
                    else:
                        if not all_path_lenghts.has_key(gene2):
                            all_path_lenghts[gene2] = {}
                        all_path_lenghts[gene2][gene1] = l
                except:
                    continue
    return all_path_lenghts


# =============================================================================
def calc_single_set_distance(G,given_gene_set):
    all_genes_in_network = set(G.nodes())
    gene_set = given_gene_set & all_genes_in_network
    all_path_lenghts = get_pathlengths_for_single_set(G,gene_set)
    all_distances = []
    for geneA in gene_set:
        all_distances_A = []
        for geneB in gene_set:
            if geneA < geneB:
                if all_path_lenghts[geneA].has_key(geneB):
                    all_distances_A.append(all_path_lenghts[geneA][geneB])
            else:
                if all_path_lenghts[geneB].has_key(geneA):
                    all_distances_A.append(all_path_lenghts[geneB][geneA])
        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)
    mean_shortest_distance = np.mean(all_distances)
    return mean_shortest_distance


# =============================================================================
def calc_set_pair_distances(G,given_gene_set1,given_gene_set2):
    all_genes_in_network = set(G.nodes())
    gene_set1 = given_gene_set1 & all_genes_in_network
    gene_set2 = given_gene_set2 & all_genes_in_network
    all_path_lenghts = get_pathlengths_for_two_sets(G,gene_set1,gene_set2)
    all_distances = []
    for geneA in gene_set1:
        all_distances_A = []
        for geneB in gene_set2:
            if geneA == geneB:
                all_distances_A.append(0)
            else:
                if geneA < geneB:
                    try:
                        all_distances_A.append(all_path_lenghts[geneA][geneB])
                    except:
                        pass
                else:
                    try:
                        all_distances_A.append(all_path_lenghts[geneB][geneA])
                    except:
                        pass
        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)
    for geneA in gene_set2:
        all_distances_A = []
        for geneB in gene_set1:
            if geneA == geneB:
                all_distances_A.append(0)
            else:
                if geneA < geneB:
                    try:
                        all_distances_A.append(all_path_lenghts[geneA][geneB])
                    except:
                        pass
                else:
                    try:
                        all_distances_A.append(all_path_lenghts[geneB][geneA])
                    except:
                        pass
        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)
    mean_shortest_distance = np.mean(all_distances)
    return mean_shortest_distance

path = "disease/"
files = os.listdir(path)

#Alert add 2022-05
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('--g1',
                              help    ='file containing gene set 1',
                              dest    ='gene_file_1',
                              default ='ID.txt',
                              type    = "string")
    
    (opts, args) = parser.parse_args()
    gene_file_1  = opts.gene_file_1

    filname=gene_file_1.split('.')[0]+'_'+'disease_results.txt'
    if os.path.isfile(filname):
        try:
            import pandas as pd
            input_data = pd.read_table(filname,header =None)
            aa=input_data.iloc[:,2]
            aa2=aa.values
            files=[x for x in files if x not in aa2]
        except Exception:
            print(filname+' is Empty!')

for i in range(0,len(files)):
    if __name__ == '__main__':
        parser = optparse.OptionParser()
        parser.add_option('-u', '--usage',
                          help    ='print more info on how to use this script',
                          action="callback", callback=print_usage)
        parser.add_option('-n',
                          help    ='file containing the network edgelist [interactome.tsv]',
                          dest    ='network_file',
                          default ='interactome.tsv',
                          type    = "string")
        parser.add_option('--g1',
                          help    ='file containing gene set 1',
                          dest    ='gene_file_1',
                          default ='ID.txt',
                          type    = "string")
        parser.add_option('--g2',
                          help    ='file containing gene set 2',
                          dest    ='gene_file_2',
                          default =files[i],
                          type    = "string")
        parser.add_option('-o',
                          help    ='file for results [disease_results.txt]',
                          dest    ='results_file',
                          default ='disease_results.txt',
                          type    = "string")
        (opts, args) = parser.parse_args()
        network_file = opts.network_file
        gene_file_1  = opts.gene_file_1
        gene_file_2  = path+opts.gene_file_2
        results_file = opts.results_file
        results_file=gene_file_1.split('.')[0]+'_'+results_file

        # checking for input:
        if gene_file_1 == 'none' or gene_file_2 == 'none':
            error_message = """
            ERROR: you must specify two files with gene sets, for example:
            python disease_distance.py --g1 MS.txt --g2 PD.txt

            For more information, type
            python disease_distance.py --usage
            
            """
            print error_message
            sys.exit(0)
        if network_file == 'interactome.tsv':
            print '> default network from "interactome.tsv" will be used'
        # read network
        G  = read_network(network_file)
        # get all genes ad remove self links
        all_genes_in_network = set(G.nodes())
        remove_self_links(G)
        # read gene set 1
        genes_A_full = read_gene_list(gene_file_1)
        # removing genes that are not in the network:
        genes_A = genes_A_full & all_genes_in_network
        if len(genes_A_full) != len(genes_A):
            print "> ignoring %s genes that are not in the network" %(
                len(genes_A_full - all_genes_in_network))
            print "> remaining number of genes: %s" %(len(genes_A))
        # read gene set 1
        genes_B_full = read_gene_list(gene_file_2)
        # removing genes that are not in the network:
        genes_B = genes_B_full & all_genes_in_network
        if len(genes_B_full) != len(genes_B):
            print "> ignoring %s genes that are not in the network" %(
                len(genes_B_full - all_genes_in_network))
            print "> remaining number of genes: %s" %(len(genes_B))
        # distances WITHIN the two gene sets:
        d_A = calc_single_set_distance(G,genes_A)
        d_B = calc_single_set_distance(G,genes_B)
        # distances BETWEEN the two gene sets:
        d_AB = calc_set_pair_distances(G,genes_A,genes_B)
        # calculate separation
        s_AB = d_AB - (d_A + d_B)/2.
        # print and save results:
        results_message = """
\"%s\"\t%s\t\"%s\"\t%s\t%s\t%s\t%s\t%s"""%(gene_file_1,len(genes_A),gene_file_2,len(genes_B),d_A,d_B,d_AB,s_AB)
        print results_message
        fp = open(results_file,'a')
        fp.write(results_message)
        fp.close()
        print "> results have been saved to %s" % (results_file)
