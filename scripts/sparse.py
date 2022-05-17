import os
import itertools
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix

from subprocess import Popen, STDOUT, PIPE

from grid import SquareGrid, HexagonGrid
from grid import generate_pattern, create_graph

def get_erc_format(graph, simple = False):
    n = len(graph.nodes)
    if simple:
        m = len(graph.edges)
    else:
        m = 2*len(graph.edges)
    

    nindex = np.zeros(n+1, dtype=np.int)
    nlist = np.zeros(m, dtype=np.int)
    eweight = np.zeros(m, dtype=np.int)
    
    nindex[0] = 0
    new_index = {i:i for i in range(len(graph.nodes))}
    #new_index = {n:i for i,n in enumerate(graph.nodes)}
    if simple:
        edges = sorted([(new_index[u], new_index[v], float(graph[u][v]["weight"])) \
            for u,v in graph.edges])
    else:
        edges = sorted([(new_index[u], new_index[v], float(graph[u][v]["weight"])) \
            for u,v in graph.edges]+[(new_index[v], new_index[u], float(graph[u][v]["weight"])) \
            for u,v in graph.edges])
    
    for i,(src,dst,wei) in enumerate(edges):
        nindex[src+1] = i+1
        nlist[i] = dst
        eweight[i] = wei
        
    for i in range(1,n+1):
        nindex[i] = max(nindex[i-1],nindex[i])
    
    return n, m, nindex, nlist, eweight

def write_CSR(graph, filename, simple):
    with open(filename, "wb") as f:
        n, m, nindex, nlist, eweight = get_erc_format(graph, simple)

        n = np.uint32(n)
        m = np.uint32(m)

        #print(n,m)
        #print(nindex, nlist, eweight)
        # === WRITE SIZE OF THE GRAPH
        f.write(n)
        f.write(m)

        # === WRITE EDGES ===
        for i in nindex:
            f.write(np.uint32(i))
        for l in nlist:
            f.write(np.uint32(l))
        for e in eweight:
            f.write(np.uint32(e))

def get_independent(graph):
    filename = os.path.dirname(os.path.abspath(__file__))+'/../data/graph.csr'
    write_CSR(graph, filename = filename, simple=False)
    p = Popen([os.path.dirname(os.path.abspath(__file__))+'/../src/bin/main', filename, '--silent', '1'],
          stdout=PIPE, stdin=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True)

    out,err = p.communicate()
    for line in out.split('\n')[:-1]:
        print(line)

def get_independent_std(graph):
    n, m, nindex, nlist, eweight = get_erc_format(graph, simple=False)

    stream = "{} {}\n".format(n, m)
        
    stream+=" ".join([str(i) for i in nindex])+"\n"
    stream+=" ".join([str(n) for n in nlist])+"\n"
    stream+=" ".join([str(e) for e in eweight])+"\n"
    
    print("Graph collected")
    p = Popen([os.path.dirname(os.path.abspath(__file__))+'/../src/bin/main', "nofile", '--silent', '2'],
          stdout=PIPE, stdin=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True)

    out,err = p.communicate(stream)
    for line in out.split('\n')[:-1]:
        print(line)
        
    return stream

def get_independent_patt(cells, filename):
    N,M = np.shape(cells)
    print("Patt begin")
    patts = generate_pattern(cells, N//2, M//2, r)
    print(f"Patt end [{len(patts)}]")

    stream = "{} {} {} ".format(N, M, len(patts))
    stream+=" ".join([f"{p[0]} {p[1]}" for p in patts])+"\n"

    #print(stream)
    
    print("Stream collected")
    p = Popen([os.path.dirname(os.path.abspath(__file__))+'/../src/bin/main', filename, '--silent', '0'],
          stdout=PIPE, stdin=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True)

    out,err = p.communicate(stream)
    for line in out.split('\n')[:-1]:
        print(line)

def launch_poligon_graph(type):
    graph = nx.Graph()
    # Cherry
    if type == "cherry":
        graph.add_edges_from([[0,1], [1,2]], weight=1)
    # Traingle
    elif type == "triangle":
        graph.add_edges_from([[0,1], [1,2], [2,0]], weight=1)
    # Square
    elif type == "square":
        graph.add_edges_from([[0,1], [1,2], [2,3], [3,0]], weight=1)

    get_independent(graph)


if __name__ == "__main__":
    mode = 3
    if mode == 0:
        r = 0.02
        #grid = HexagonGrid([-4,4], r, 400, 400)
        grid = SquareGrid([-4,4], r, 200, 200)
        G = create_graph(grid.cells, r)
        print(len(G.nodes), len(G.edges))
        print("Graph created")
        get_independent(G)
    elif mode == 1:
        r = 0.1
        grid = SquareGrid([-4,4], r, 80, 80)
        G = create_graph(grid.cells, r)
        
        print(len(G.nodes), len(G.edges))
        print("Graph created")
        get_independent(G)
        #get_independent_std(G)
    elif mode == 2:
        #r = 0.005
        #grid = SquareGrid([-4,4], r, 1000, 1000)
        r = 0.01
        grid = SquareGrid([-4,4], r, 1600, 1600)
        get_independent_patt(grid.cells)
    elif mode == 3:
        size = [500, 500]
        for r in [0.01]:
            print(f"r = {r}")
            grid = SquareGrid([-4,4], r, size[0], size[1])
            get_independent_patt(grid.cells, f"../data/r:{r}_size:{size[0]}x{size[1]}")
    else:
        print("Mode not found")
