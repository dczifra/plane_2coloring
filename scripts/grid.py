import itertools

import numpy as np
import networkx as nx

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
            col_of_cells = [self.base]
            for row in range(n-1):
                col_of_cells.append(Cell(col_of_cells[-1].shift(-2), 6, length))
            self.base = Cell(self.base.shift(0 if(col%2==0) else -1), 6, length)
            self.cells.append(col_of_cells)
    def get_cells(self):
        ind = 0
        for row in self.cells:
            for c in row:
                yield ind,c
                ind+=1

class SquareGrid:
    def __init__(self, center, length, n, m):
        self.cells = []
        self.base = Cell(center, 4, length)
        for col in range(m):
            col_of_cells = [self.base]
            for row in range(n-1):
                col_of_cells.append(Cell(col_of_cells[-1].shift(0), 4, length))
            self.base = Cell(self.base.shift(-1), 4, length)
            self.cells.append(col_of_cells)
    
    def get_cells(self):
        ind = 0
        for row in self.cells:
            for c in row:
                yield ind,c
                ind+=1

EPS = 1e-12
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

# TODO:
#      * Tets generate pattern
#      * Add directions
def generate_pattern(cells, x0, y0, r):
    N,M = np.shape(cells)
    patt = []
    for x in range(N//2):
        for y in range(M//2):
            if(one_dist(cells[x0][y0], cells[x0+x][y0+y], r)):
                patt.append([x,y])
    return patt

def create_graph(cells, r):
    N,M = np.shape(cells)
    print("Patt begin")
    patts = generate_pattern(cells, N//2, M//2, r)
    print("Patt end")

    graph = nx.Graph()
    graph.add_nodes_from(range(N*M))

    for i,j in itertools.product(range(N), range(M)):
        coord = i*M+j
        # Add edges to 1-distance:
        for x,y in patts:
            # Add 4 directions
            for a,b in itertools.product([1,-1], [1,-1]):
                new_coord = ((i+x*a)%N)*M+((j+y*b)%M)
                graph.add_edge(coord,new_coord, weight=1)
    return graph                
