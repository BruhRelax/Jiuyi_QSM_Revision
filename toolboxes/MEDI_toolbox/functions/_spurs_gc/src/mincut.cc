#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "graph.h"
#include <string.h>
#include <vector>
#include "mincut.h"

using namespace std;

//int mincut(vector<float>& sourcesink, vector<float>& remain,
//    float& flow, vector<float>& cutside)
int mincut(std::vector<sourcesink_t>& sourcesink, 
           std::vector<remain_t>& remain,
           float& flow, 
           std::vector<cutside_t>& cutside)
{
  //Declarations
  //float *sourcesinkVal;
  //float *remainVal;

  //int sourcerowLen, sourcecolLen;
  //int remainrowLen, remaincolLen;

  //Get matrix sourcesink
  //sourcesinkVal = &sourcesink[0];
  //sourcerowLen = 3;
  //sourcecolLen = sourcesink.size()/sourcerowLen;

  //Get matrix remain
  //remainVal    = &remain[0];
  //remainrowLen = 4;
  //remaincolLen = remain.size()/remainrowLen;	

  int num_nodes = -1;
  for(unsigned int i=0; i<sourcesink.size(); i++) {
    if (sourcesink[i].node > num_nodes) num_nodes = sourcesink[i].node;
  }
  num_nodes += 1;

  vector<Graph::node_id> nodes(num_nodes);	// dynamic allocation

  //Graph *g = new Graph();
  Graph g;	

  for(unsigned int t=0; t<sourcesink.size(); t++)
    //nodes[t] = g.add_node();
    nodes[(int)sourcesink[t].node] = g.add_node();

  /* capacity of arcs between terminals and nodes*/
  for(unsigned int k=0; k<sourcesink.size(); k++) // sourcesink
    g.set_tweights(nodes[(int)sourcesink[k].node], 
        sourcesink[k].weight_source, 
        sourcesink[k].weight_sink); 

  for(unsigned int h=0; h<remain.size(); h++)	
    g.add_edge(
        nodes[(int)remain[h].node_start],
        nodes[(int)remain[h].node_end],
        remain[h].weight,
        remain[h].inverseweight);

  //cout << "flow = " << g.get_flow() ;
  //Graph::flowtype flow = g.maxflow();
  flow = g.maxflow();
  //cout << " ==> " << g.get_flow() << endl; 


  //Allocate memory and assign output pointer
  //plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); //mxReal is our data-type
  //plhs[0] = mxCreateNumericArray(1,dims, mxSINGLE_CLASS, mxREAL);
  //int dims2[] = {sourcecolLen,2};
  //plhs[1] = mxCreateNumericArray(2,dims2, mxSINGLE_CLASS, mxREAL);
  cutside.resize(sourcesink.size());

  for(unsigned int i=0;i<cutside.size();i++)
  {
    //assigns the number of the nodes
    cutside[i].node = sourcesink[i].node;

    //If o nó becomes assigned to source we assign the value 0
    //Otherwise we assign the value 1.
    if (g.what_segment(nodes[(int)sourcesink[i].node]) == Graph::SOURCE)
      cutside[i].segment = 0;
    else
      cutside[i].segment = 1;
  }

  return 0;
} 
