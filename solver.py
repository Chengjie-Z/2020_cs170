import networkx as nx
import numpy as np
from parse import read_input_file, write_output_file
from utils import is_valid_network, average_pairwise_distance, average_pairwise_distance_fast
import sys
from itertools import combinations
import random


#Algo 1. go through all shortest path trees, very close to the optimal but still cost time, to obtain the quality of outputs we choose this to be our main algorithm
def solve(G):
    """
    Args:
        G: networkx.Graph

    Returns:
        T: networkx.Graph
    """
    opt = 0
    for i in range(G.number_of_nodes()):
        t = nx.single_source_dijkstra_path(G, i)
        t_edge = set()
        for j in range(G.number_of_nodes()):
            if len(t[j]) > 1:
                n = 1
                while n < len(t[j]):
                    t_edge.add((t[j][n-1], t[j][n], G[t[j][n-1]][t[j][n]]["weight"]))
                    n += 1 
        T = nx.Graph()
        T.add_nodes_from(range(G.number_of_nodes()))
        edge = list(t_edge)
        T.add_weighted_edges_from(edge)
        T = inherit(T)
        if i == 0:
            opt = average_pairwise_distance_fast(T)
            opt_T = T
        else:
            ave = average_pairwise_distance_fast(T)
            if ave < opt:
                opt = ave
                opt_T = T
    return opt_T

# genetic algorithm (a monte carlo algorithm to decide whether cut leaves in the ST)
def get_child(DNA_p1, DNA_p2):
    DNA_children = []
    abr = 0.3
    for i in range(10):
        DNA_children.append([
            random.choice([DNA_p1[j], DNA_p2[j]])
            if random.random() < abr else random.randint(0, 1)
            for j in range(len(DNA_p1))
        ])
    # DNA_children.append(DNA_p1)
    # DNA_children.append(DNA_p2)
    return DNA_children

def cut_leaf_DNA(G, leaf_set, DNA):
    Gm = nx.Graph(G)
    for i in range(len(leaf_set)):
        if DNA[i]:
            Gm.remove_node(leaf_set[i])
    return Gm

def calculate_pairwise_distance_with_DNA(G, leaf_set, DNA):
    Gm = cut_leaf_DNA(G, leaf_set, DNA)
    D = average_pairwise_distance_fast(Gm)
    return Gm, D

def inherit(G):
    # leaf_set, Gnl = cut_leaf(G)
    eps = 0.001 # float to compare the error
    leaf_set = find_leaves(G)
    Gm = nx.Graph(G)
    leaf_size = len(leaf_set)
    DNA1 = np.zeros(leaf_size)
    DNA2 = np.ones(leaf_size)
   
    win_list = [(calculate_pairwise_distance_with_DNA(Gm, leaf_set, DNA1)[1], DNA1)
        , (calculate_pairwise_distance_with_DNA(Gm, leaf_set, DNA2)[1], DNA2)]
    if win_list[0][0] > win_list[1][0]: win_list[0], win_list[1] = win_list[1], win_list[0]
    while True:
        old_D0, old_D1 = win_list[0][0], win_list[1][0]
        DNA_child = get_child(win_list[0][1], win_list[1][1])
        for DNA in DNA_child:
            G_DNA, D = calculate_pairwise_distance_with_DNA(Gm, leaf_set, DNA)
            
            if D < win_list[0][0] or D < win_list[1][0]:
                win_list[1] = (D, DNA)
                if win_list[0][0] > win_list[1][0]: win_list[0], win_list[1] = win_list[1], win_list[0]
        # if parents are not subsituted, we halt the algorithm
        if abs(old_D0 - win_list[0][0]) <= eps and abs(old_D1 - win_list[1][0]) <= eps:
            break
    return cut_leaf_DNA(Gm, leaf_set, win_list[0][1])

# to find all vertices with degree 1 in a graph
def find_leaves(G):
    leaf_set = []
    for degree in G.degree:
        if degree[1] == 1:
            leaf_set.append(degree[0])
    return leaf_set

def cut_leaf(G):
    leaf_set = find_leaves(G)
    Gm = nx.Graph(G)
    Gm.remove_nodes_from(leaf_set)
    return Gm, leaf_set





#Algo 2. Here we provided a very fast but not that optimal approximation algorithm, which is to find a mid point then to find a shortest path tree 
def modifed_find_leaves(G):
    leaf_set = []
    value_set = sorted(d for n, d in G.degree())
    value = value_set[int(random.randint(0, 1) * 0.5 * G.number_of_nodes())] # randomly choose leaves with a low degree
    for degree in G.degree:
        if degree[1] <= value:
            leaf_set.append(degree[0])
    return leaf_set

def mid_dijkstra(G):
    def find_mid_point(G):
        leaf_set = modifed_find_leaves(G)
        if len(leaf_set) == 0:
            n_d = (-1, 0)
            for degree in G.degree:
                if degree[1] > n_d[1]: n_d = degree
            return n_d[0]
        Gm = nx.Graph(G)
        Gm.add_node(-1)
        for leaf in leaf_set:
            Gm.add_edge(-1, leaf)
        mid_points = find_last_bfs(Gm, -1)
        return mid_points

    def find_last_bfs(G, src):
        Gm = nx.Graph(G)
        pq = [src]
        traversal = []
        while len(traversal) + len(pq) < len(G):
            n = pq.pop(0)
            traversal.append(n)
            nei = Gm.neighbors(n)
            [pq.append(n) for n in nei if (n not in traversal) and (n not in pq)]
            # Gm.remove_edges_from(nei)
        return pq[-1]
    
    mid_point = find_mid_point(G)
    Dip = nx.single_source_dijkstra(G, mid_point)[1]
    Gm = nx.Graph(G)
    Gm.remove_edges_from(G.edges)
    for i in range(len(Dip)):
        if len(Dip[i]) < 2: continue
        nx.add_path(Gm, Dip[i])
    for edge in Gm.edges:
        u, v = edge
        Gm[u][v]["weight"] = G.get_edge_data(u, v, "weight")["weight"]
    return Gm

def starting_point_algo(G):
    T = mid_dijkstra(G)
    opt = inherit(T)
    for _ in range(5):
        t = inherit(T)
        if average_pairwise_distance_fast(t) < average_pairwise_distance_fast(opt):
            opt = t
    return opt






#Algo 3. an MST based algorithm, but the answer is usually not as optimal as Algo 2
def mst_solution(G):
    T = nx.minimum_spanning_tree(G)
    T = inherit(T)
    return T





#Algo 4. here is a brute-force search algorithm to go through all spanning tree, but it's very slow
def min_ave_brute_force_search(G):
    e = list(G.edges())
    all_possible_edges = combinations(e, G.number_of_nodes() - 1) 

    # calculate the number of spanning tree for the sake of convience
    mat = nx.to_numpy_matrix(G)
    mat = np.where(mat > 0, 1, mat)
    np.fill_diagonal(mat, [G.degree[i] for i in range(G.number_of_nodes())])
    np.delete(mat, 0, 0)
    np.delete(mat, 0, 1)
    n = np.linalg.det(mat)

    #start to select one tree
    num = 0
    for t in all_possible_edges:
        lis = list(t)
        T = nx.Graph()
        T.add_nodes_from(range(G.number_of_nodes()))
        T.add_edges_from(lis)
        if nx.is_connected(T):
            num += 1
            for s in lis:
                T[s[0]][s[1]]["weight"] = G[s[0]][s[1]]["weight"]
            lis = []
            for k in T.degree:
                if k[1] == 1:
                    lis.append(k[0])
            for z in lis:
                T.remove_node(z) 
            # T = inherit(T)
            if num == 1:
                opt = average_pairwise_distance_fast(T)
                opt_T = T
            else:
                ave = average_pairwise_distance_fast(T)
                if ave < opt:
                    opt = ave
                    opt_T = T
        if num == n:
            break
    return opt_T
# Here's an example of how to run your solver.
# Usage: python3 solver.py test.in

if __name__ == '__main__':
    #assert len(sys.argv) == 2
    #path = sys.argv[1]
    G = read_input_file("medium-235.in")
    T = solve(G)
    assert is_valid_network(G, T)
    print("Average pairwise distance: {}".format(average_pairwise_distance(T)))
    #write_output_file(T, 'out/test.out')
