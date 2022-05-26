#pragma once

#include <map>
#include <vector>
#include <iostream>

#include "ECLgraph.h"

struct Triple;
std::ostream& operator<<(std::ostream& os, Triple& t);

struct Triple{
    int x,y,z;

    Triple(){}
    Triple(int _x, int _y, int _z): x(_x), y(_y), z(_z) {}

    void add(Triple& t){
        x+=t.x;
        y+=t.y;
        z+=t.z;
    }

    bool pos(){
        return x+y+z > 0;
    }

    /**
     * [Description]:
     *   Returns false, if the given coordinates are out of the bounds
     *   of the hexagon grid
     * [Args]:
     *   N: size of the hexagon grid
     * [Return]: bool
     * */
    inline bool valid(int N){
        return (std::abs(x)+std::abs(y)+std::abs(z) <= 2*(N-1)+1) && ((x+y+z == -1) || (x+y+z == 1));
    }

    /**
     * [Description]:
     *   This Function is looking for the neighbouring 6 hexagon, and finds
     *   the one, where the coordinates are within the borders.
     *   After that we switch to the new (relative) coordinates.
     *  [Args]: N: size of the hexagon grid
     *  [Return]: None
     * */
    void rescale(const int N){
        if(valid(N)) return;
        else{
            std::vector<Triple> centers = {
                Triple(2*N, -N,-N), Triple(-2*N, N, N), Triple(-N, 2*N, -N),
                Triple(N, -2*N, N), Triple(-N,-N,2*N), Triple(N,N,-2*N)};
            
            Triple t;
            for(Triple& c: centers){
                t.x=x, t.y=y,t.z=z;
                t.add(c);
                if(t.valid(N)){
                    *this=t;
                    return;
                }
            }
            std::cout<<"\rCenter not found    ";
        }
    }

    // TODO N should be ull ???
    ull get_index(const ull N){
        //return x*N*N+y*N+z;
        ull NN = 2*N+1;
        return (x+N)*NN*NN+(y+N)*NN+(z+N);
    }

    static Triple get_coords(ull ind, ull N){
        ull NN = 2*N+1;
        int x = ind/(NN*NN);
        ull temp = ind%(NN*NN);
        int y = temp/NN;
        int z = temp % NN;
        return Triple(x-N, y-N, z-N);
    }
};

std::ostream& operator<<(std::ostream& os, Triple& t){
    os<<t.x<<" "<<t.y<<" "<<t.z;
    return os;
}

struct TripleLess{
    bool operator()(const Triple& a, const Triple& b) const{
        if(a.x<b.x) return true;
        else if((a.x == b.x) && (a.y < b.y)) return true;
        else if((a.x == b.x) && (a.y == b.y) && (a.z < b.z)) return true;
        else return false;
    }
};

class TriangleGrid{
public:
    TriangleGrid(){}

    ECLgraph create_hexagon(){
        ECLgraph g;
        // === Read grid parameters ===
        std::cin>>N>>K;

        // === Read pattern ===
        std::vector<Triple> patt_pos(K);
        std::vector<Triple> patt_neg(K);
        for(int i=0;i<K;i++) std::cin>>patt_pos[i].x>>patt_pos[i].y>>patt_pos[i].z;
        for(int i=0;i<K;i++) std::cin>>patt_neg[i].x>>patt_neg[i].y>>patt_neg[i].z;

        // === Init ECL graph ===
        ull node_ind = 0;
        ull bigN = 2*N*((2*N+1)*(2*N+1)+(2*N+1)+2*N);
        for(ull coord=0;coord<bigN;coord++){
            Triple node=Triple::get_coords(coord, N);
            if(node.valid(N)){
                ind_to_triple[node_ind] = node;
                triple_to_ind[node] = node_ind;
                node_ind++;
            }
        }
        g.nodes = triple_to_ind.size();
        g.edges = g.nodes*K;
        g.nindex = new ull[g.nodes + 1];
        g.nlist = new int[g.edges];
        g.eweight = NULL; // We dont use it
        for(int i=0;i<g.nodes+1;i++) g.nindex[i]=0;
        for(ull i=0;i<g.edges;i++) g.nlist[i]=0;
        std::cout<<"Read all right "<<N<<" "<<K<<" "<<g.nodes<<std::endl;

        // === INIT EDGES ===
        ull ind = 0;
        for(auto p_coord: triple_to_ind){
            Triple node=p_coord.first;
            //std::cout<<"["<<p_coord.second<<"]"<<node;

            ull edge_ind = 0;
            std::vector<int> edges(patt_pos.size());

            for(Triple& p: (node.pos()?patt_pos:patt_neg)){
                Triple neigh(node);
                neigh.add(p);
                neigh.rescale(N);

                //std::cout<<"    "<<neigh;

                //ull new_coord = neigh.get_index(N);
                ull new_coord = triple_to_ind[neigh];
                edges[edge_ind++] = new_coord;
            }
        
            std::sort(edges.begin(), edges.end());
            // === Update CSR representation ===
            if(edges.size() == 0) continue;
            int last = -1;
            for(size_t i=0;i<edges.size();i++){
                if(edges[i]!=last){
                    g.nlist[ind++] = edges[i];
                    g.nindex[p_coord.second+1] = ind;
                    last = edges[i];
                }
            }
        }

        std::cout<<"Nodes: "<<g.nodes<<" Edges: "<<g.edges<<std::endl;
        return g;
    }

    void log_node(std::ofstream& log, ull ind){
        log<<ind_to_triple[ind]<<",";
    }

private:
    ull N,K;
    std::map<Triple, ull, TripleLess> triple_to_ind;
    std::map<ull, Triple> ind_to_triple;
};