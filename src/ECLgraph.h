/*
Copyright (c) 2016-2020, Texas State University. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
   * Neither the name of Texas State University nor the names of its
     contributors may be used to endorse or promote products derived from
     this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL TEXAS STATE UNIVERSITY BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Author: Martin Burtscher
*/


#ifndef ECL_GRAPH
#define ECL_GRAPH

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <algorithm>

#define ull unsigned long

struct ECLgraph {
  ull nodes;
  ull edges;
  ull* nindex;
  int* nlist;
  int* eweight;
};

ECLgraph createECLgraph(){
  ECLgraph g;
  ull N,M, K;
  std::cin>>N>>M>>K;
  std::vector<std::pair<int, int>> patt(K);
  for(int i=0;i<K;i++) std::cin>>patt[i].first>>patt[i].second;

  g.nodes = N*M;
  g.edges = g.nodes*K*4;
  ull estimated = g.edges;

  g.nindex = new ull[g.nodes + 1];
  g.nlist = new int[g.edges];
  g.eweight = NULL; // We dont use it
  for(int i=0;i<g.nodes+1;i++) g.nindex[i]=0;
  for(ull i=0;i<g.edges;i++) g.nlist[i]=0;

  std::vector<std::pair<int,int>> signs = {{-1,-1}, {-1,1}, {1,-1}, {1,1}};
  ull ind = 0;
  for(int coord=0;coord<N*M;coord++){
    int i = coord%N;
    int j = coord/N;
    std::vector<int> edges(patt.size()*4);
    ull edge_ind = 0;
    for(std::pair<int,int>& p: patt){
      int x = p.first;
      int y = p.second;
      for(std::pair<int,int>& sign: signs){
        int new_x = i+x*sign.first;
        int new_y = j+y*sign.second;
        int new_coord = ((new_x<0?new_x+N:new_x)%N)+((new_y<0?new_y+M:new_y)%M)*N;
        edges[edge_ind++] = new_coord;
      }
    }

    std::sort(edges.begin(), edges.end());
    // === Update CSR representation ===
    if(edges.size() == 0) continue;
    int last = -1;
    for(int i=0;i<edges.size();i++){
      if(edges[i]!=last){
        g.nlist[ind++] = edges[i];
        g.nindex[coord+1] = ind;
        last = edges[i];
      }
    }
  }
  g.edges = ind;
  std::cout<<"True: "<<g.edges<<" estimated: "<<estimated<<std::endl;
  return g;
}

ECLgraph readECLgraph_std(){
  ECLgraph g;
  std::cin>>g.nodes>>g.edges;

  g.nindex = (ull*)malloc((g.nodes + 1) * sizeof(g.nindex[0]));
  g.nlist = (int*)malloc(g.edges * sizeof(g.nlist[0]));
  g.eweight = (int*)malloc(g.edges * sizeof(g.eweight[0]));
  if ((g.nindex == NULL) || (g.nlist == NULL) || (g.eweight == NULL)){
      fprintf(stderr, "ERROR: memory allocation failed\n\n");  exit(-1);
  }
    
  for(int i=0;i<g.nodes+1;i++){
      std::cin>>(g.nindex[i]);
  }
  for(ull i=0;i<g.edges;i++){
      std::cin>>(g.nlist[i]);
  }
  for(ull i=0;i<g.edges;i++){
      std::cin>>(g.eweight[i]);
  }
    
  //std::cout<<"Nodes2 "<<g.nodes<<g.edges<<std::endl;
  return g;
}

ECLgraph readECLgraph(const char* const fname)
{
  ECLgraph g;
  int cnt;

  FILE* f = fopen(fname, "rb");  if (f == NULL) {fprintf(stderr, "ERROR: could not open file %s\n\n", fname);  exit(-1);}
  cnt = fread(&g.nodes, sizeof(g.nodes), 1, f);  if (cnt != 1) {fprintf(stderr, "ERROR: failed to read nodes\n\n");  exit(-1);}
  cnt = fread(&g.edges, sizeof(g.edges), 1, f);  if (cnt != 1) {fprintf(stderr, "ERROR: failed to read edges\n\n");  exit(-1);}
  if ((g.nodes < 1) || (g.edges <= 0)) {fprintf(stderr, "ERROR: node or edge count too low\n\n");  exit(-1);}

  g.nindex = (ull*)malloc((g.nodes + 1) * sizeof(g.nindex[0]));
  g.nlist = (int*)malloc(g.edges * sizeof(g.nlist[0]));
  g.eweight = (int*)malloc(g.edges * sizeof(g.eweight[0]));
  if ((g.nindex == NULL) || (g.nlist == NULL) || (g.eweight == NULL)) {fprintf(stderr, "ERROR: memory allocation failed\n\n");  exit(-1);}

  cnt = fread(g.nindex, sizeof(g.nindex[0]), g.nodes + 1, f);  if (cnt != g.nodes + 1) {fprintf(stderr, "ERROR: failed to read neighbor index list\n\n");  exit(-1);}
  cnt = fread(g.nlist, sizeof(g.nlist[0]), g.edges, f);  if (cnt != g.edges) {fprintf(stderr, "ERROR: failed to read neighbor list\n\n");  exit(-1);}
  cnt = fread(g.eweight, sizeof(g.eweight[0]), g.edges, f);

  /*
  for(int i=0;i<g.nodes+1;i++) std::cout<<g.nindex[i]<<" ";
  std::cout<<std::endl;
  for(int i=0;i<g.edges;i++) std::cout<<g.nlist[i]<<" ";
  std::cout<<std::endl;
  for(int i=0;i<g.edges;i++) std::cout<<g.eweight[i]<<" ";
  std::cout<<std::endl;
  */
  
  if (cnt == 0) {
    free(g.eweight);
    g.eweight = NULL;
  } else {
    if (cnt != g.edges) {fprintf(stderr, "ERROR: failed to read edge weights\n\n");  exit(-1);}
  }
  fclose(f);

  return g;
}

void writeECLgraph(const ECLgraph g, const char* const fname)
{
  if ((g.nodes < 1) || (g.edges <= 0)) {fprintf(stderr, "ERROR: node or edge count too low\n\n");  exit(-1);}
  int cnt;
  FILE* f = fopen(fname, "wb");  if (f == NULL) {fprintf(stderr, "ERROR: could not open file %s\n\n", fname);  exit(-1);}
  cnt = fwrite(&g.nodes, sizeof(g.nodes), 1, f);  if (cnt != 1) {fprintf(stderr, "ERROR: failed to write nodes\n\n");  exit(-1);}
  cnt = fwrite(&g.edges, sizeof(g.edges), 1, f);  if (cnt != 1) {fprintf(stderr, "ERROR: failed to write edges\n\n");  exit(-1);}

  cnt = fwrite(g.nindex, sizeof(g.nindex[0]), g.nodes + 1, f);  if (cnt != g.nodes + 1) {fprintf(stderr, "ERROR: failed to write neighbor index list\n\n");  exit(-1);}
  cnt = fwrite(g.nlist, sizeof(g.nlist[0]), g.edges, f);  if (cnt != g.edges) {fprintf(stderr, "ERROR: failed to write neighbor list\n\n");  exit(-1);}
  if (g.eweight != NULL) {
    cnt = fwrite(g.eweight, sizeof(g.eweight[0]), g.edges, f);  if (cnt != g.edges) {fprintf(stderr, "ERROR: failed to write edge weights\n\n");  exit(-1);}
  }
  fclose(f);
}

void freeECLgraph(ECLgraph &g)
{
  if (g.nindex != NULL) free(g.nindex);
  if (g.nlist != NULL) free(g.nlist);
  if (g.eweight != NULL) free(g.eweight);
  g.nindex = NULL;
  g.nlist = NULL;
  g.eweight = NULL;
}

#endif
