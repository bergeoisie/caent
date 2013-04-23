#ifndef _caentHelper_h
#define _caentHelper_h

#include <boost/config.hpp> // put this first to suppress some VC++ warnings

//#include "tnt.h"
//#include "jama_eig.h"

#include <iostream>
#include <iterator>
#include <algorithm>
#include <time.h>
#include <fstream>
#include <tr1/unordered_map>
#include <stdlib.h>
#include <stdio.h>
#include <queue>
#include <sstream>
#include <ctime>


#include <armadillo>


#include <boost/utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>

using namespace std;
//using namespace TNT;
//using namespace JAMA;
using namespace boost;
using namespace arma;

typedef property<edge_name_t,string> lbl;
typedef property<vertex_name_t,string> vlbl;

typedef adjacency_list<vecS,vecS,bidirectionalS, vlbl, lbl> Graph;

typedef graph_traits<Graph>::vertex_iterator VI;
typedef graph_traits<Graph>::edge_iterator EI;
typedef graph_traits<Graph>::out_edge_iterator OEI;
typedef vector<graph_traits<Graph>::vertex_descriptor> VVec;


typedef std::tr1::unordered_map<string,graph_traits<Graph>::vertex_descriptor> VertexMap;


Graph HigherNBlock(Graph G, int n);
void PrintFullGraphInfo(Graph G, ostream& os = cout);
Graph Rename(Graph G, vector<string> names);
vector<string> RandomStringGenerator(int ,double);
int flip(double);
Graph RightResolvingRepn(Graph G);
Graph inducedRp(Graph &);
VVec compatibleSet(Graph&,
 VVec,
 string);
void printVVec(VVec);
bool VVecSubset(VVec&,VVec&);
string ssVVec(VVec);
//double PFEigenvalue(Graph &);
double MyPFEigenvalue(Graph &);
int oneCounter(vector<string> s);
Graph Trim(Graph &);
double standardDeviation(vector<double>&,double,int);

#endif
