import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix

from subprocess import Popen, STDOUT, PIPE

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
    p = Popen([os.path.dirname(os.path.abspath(__file__))+'/../src/bin/main', filename, '--silent'],
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
    p = Popen([os.path.dirname(os.path.abspath(__file__))+'/../src/bin/main', "nofile", '--silent'],
          stdout=PIPE, stdin=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True)

    out,err = p.communicate(stream)
    for line in out.split('\n')[:-1]:
        print(line)
        
    return stream

class Coord:
    def __init__(self, x,y):
        self.x,self.y = x,y
    
    def __repr__(self):
        return f"({self.x:.3f}, {self.y:.3f})"

class Cell:
    def __init__(self, center, n, length, base_dir=0):
        self.n = n
        self.node_coords = []
        self.center = Coord(*center)
        for i in range(n):
            x = self.center.x+length*np.cos(base_dir*np.pi + 2*i*np.pi/n)
            y = self.center.y+length*np.sin(base_dir*np.pi + 2*i*np.pi/n)
            self.node_coords.append(Coord(x,y))
    def xshift(self):
        return 2*self.node_coords[0].x
    def shift(self, k = 0):
        v1 = self.node_coords[k%self.n]
        v2 = self.node_coords[(k+1)%self.n]
        return [(v1.x+v2.x)-self.center.x, (v1.y+v2.y)-self.center.y]

class HexagonGrid:
    def __init__(self, center, length, n, m):
        self.cells = []
        self.base = Cell(center, 6, length)
        for col in range(m):
            self.cells.append(self.base)
            for row in range(n-1):
                self.cells.append(Cell(self.cells[-1].shift(-2), 6, length))
            self.base = Cell(self.base.shift(0 if(col%2==0) else -1), 6, length)

EPS = 0.000000001
def one_dist(c1, c2, r):
    dist_center = np.sqrt((c1.center.x-c2.center.x)**2 + (c1.center.y-c2.center.y)**2)
    # If the polygons are too close, or too far away
    if (dist_center + 2*r < 1.0-EPS) or (dist_center - 2*r > 1.0+EPS):
        return False
    else:
        smaller = False
        larger = False
        for v1 in c1.node_coords:
            for v2 in c2.node_coords:
                l2 = np.sqrt((v1.x-v2.x)**2 + (v1.y-v2.y)**2)
                #print(l2)
                if(l2 < 1.0-EPS): smaller = True
                elif(l2 >= 1.0+EPS): larger = True
                else: return True # 1-dist found

                # Bolzano theorem ==> there exist 1-dist
                if(larger and smaller):
                    return True
        return False
                
def create_graph(cells, r):
    N = len(cells)
    graph = nx.Graph()
    graph.add_nodes_from(range(N))

    for i in range(N):
        for j in range(i+1, N):
            if(one_dist(cells[i],cells[j], r)):
                graph.add_edge(i,j, weight=1)
    return graph

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
    r = 0.02
    grid = HexagonGrid([-4,4], r, 100, 100)
    G = create_graph(grid.cells, r)
    print(len(G.nodes), len(G.edges))
    print("Graph created")

    get_independent(G)
    #get_independent_std(G)
