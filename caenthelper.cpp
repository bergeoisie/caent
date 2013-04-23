#include <boost/config.hpp> // put this first to suppress some VC++ warnings

//#include "tnt.h"
//#include "jama_eig.h"

#include "caenthelper.h"

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

boost::random::mt19937 gen(time(0));


// Checks to see if one VVec (vector of vertex descriptors) is a subset of another
bool VVecSubset(VVec& a, VVec& b)
{
    VVec::iterator ait,bit;
    for(ait = a.begin(); ait < a.end(); ait++)
    {
        bool found=false;
        for(bit = b.begin(); bit < b.end(); bit++)
        {
            if(*ait == *bit)
            {
                found = true;
            }
        }
        if(!found)
        {
            // printVVec(a);
            //    cout << " is NOT a subset of ";
            //    printVVec(b); cout << endl;
            return false;
        }
    }
    
    if(a.size() == b.size())
    {
        return false;
    }
    //  printVVec(a); cout << "is a subset of "; printVVec(b); cout << endl;
    return true;
}

void printVVec(VVec a)
{
    VVec::iterator it;
    cout << "{";
    
    for(it = a.begin(); it < a.end()-1; it++)
    {
        cout << *it << " ";
    }
    cout << *it << "}";
}

// Makes a VVec (a vector of vertex descriptors) into a stringstream
string ssVVec(VVec a)
{
    stringstream ss;
    VVec::iterator it;
    ss << "{";
    
    for(it = a.begin(); it < a.end()-1; it++)
    {
        ss << *it << " ";
    }
    ss << *it << "}";
    
    return ss.str();
}

Graph HigherNBlock(Graph G, int n)
{

    Graph B;
    
    int i=0;
    
    string currentpWord,currentqWord,currenteWord,first,last,newpWord,newqWord,neweWord;
    
    // This map stores n strings and their terminal vertices. Before inserting a new
    // string, we will check to see if the string is already in the map. If so, we check
    // their terminal vertices. If they match, we do not add the string, if not, return false
    map<string,graph_traits<Graph>::vertex_descriptor> GVertexMap;
    //  map<string,string>::iterator tmi;
    
    // The OEI pairs allow us to keep both the current location of the iterator, its end
    // and the string it is carrying
    
    typedef tuple<graph_traits<Graph>::out_edge_iterator,
    graph_traits<Graph>::out_edge_iterator,
    string> GHOEITuple;
    
    
    property_map<Graph,edge_name_t>::type
    g_name = get(edge_name,G);
    
    property_map<Graph,vertex_name_t>::type
    b_vname = get(vertex_name,B);
    
    property_map<Graph,edge_name_t>::type
    gename = get(edge_name,G);
    
    property_map<Graph,vertex_name_t>::type
    g_vname = get(vertex_name,G);
    
    graph_traits<Graph>::vertex_iterator gvi,gvi_end;
    
    graph_traits<Graph>::out_edge_iterator goi,goi_end,gcurrenti,gcurrenti_end;
    
    stack<GHOEITuple> GHOEIStack;
    
    int ChunkSize,BlockSize,NameChunkSize,NameBlockSize,HomChunkSize,HomBlockSize;
    
    bool done= false;
    if(n <= 1)
    {
        cout << "n must be greater than 1, it is " << n << endl;
        return G;
    }
    
    for(tie(gvi,gvi_end)=vertices(G); gvi!=gvi_end; gvi++)
    {
        i=1;
        tie(goi,goi_end)=out_edges(*gvi,G);
        
        ChunkSize = (gename(*goi)).size();
        BlockSize = ChunkSize*(n-1);
        
        //    cout << "We've started looking at vertex " << g_vname(*gvi) << endl;
        
        GHOEIStack.push(GHOEITuple(goi,goi_end,""));
        
        while(!GHOEIStack.empty())
        {
            GHOEITuple current = GHOEIStack.top();
            
            gcurrenti = get<0>(current);
            gcurrenti_end = get<1>(current);
            currenteWord = get<2>(current);
            
            boost::graph_traits<Graph>::edge_descriptor e=*gcurrenti;
            
            if(i < n) // We want to go deeper
            {
                graph_traits<Graph>::out_edge_iterator ni,ni_end;
                
                tie(ni,ni_end)=out_edges(target(e,G),G);
                GHOEIStack.push(GHOEITuple(ni,ni_end,currenteWord+gename(e)));
                i++;
                //	      cout << "i value is " << i << endl;
                
                //	      cout << "Our new edge word is " << currenteWord+gename(e) << endl;
            }
            else // We're the deepest we want to go. 
            {
                graph_traits<Graph>::vertex_descriptor s,t;
                map<string,graph_traits<Graph>::vertex_descriptor>::iterator si,ti;
                GHOEIStack.pop();
                neweWord = currenteWord+gename(e);
                //    cout << "Our final length " << n << " word is " << neweWord << endl;
                first = string(neweWord,0,BlockSize);
                last = string(neweWord,ChunkSize,BlockSize);
                si = GVertexMap.find(first);
                
                if(si==GVertexMap.end())
                {
                    s = add_vertex(B);
                    put(b_vname,s,first);
                    GVertexMap[first]=s;
                    //   cout << "We have created the vertex " << first << endl;
                }
                else {
                    s = (*si).second; 
                }
                ti = GVertexMap.find(last);
                if(ti==GVertexMap.end())
                {
                    t = add_vertex(B);
                    put(b_vname,t,last);
                    GVertexMap[last]=t;
                    //  cout << "We have created the vertex " << last << endl;
                }
                else {
                    t = (*ti).second;
                }	      
                
                add_edge(s,t,neweWord, B);
                gcurrenti++;
                if(gcurrenti!=gcurrenti_end) {
                    GHOEIStack.push(GHOEITuple(gcurrenti,gcurrenti_end,currenteWord));
                }
                else { // Our top iterator has finished. We need to go down a level (or two..)
                    done = false;
                    while(!done && !GHOEIStack.empty()) {
                        current = GHOEIStack.top();
                        GHOEIStack.pop();
                        gcurrenti = get<0>(current);
                        gcurrenti_end = get<1>(current);
                        currenteWord = get<2>(current);
                        gcurrenti++;
                        i--;
                        if(gcurrenti!=gcurrenti_end)
                        {
                            done=true;
                            GHOEIStack.push(GHOEITuple(gcurrenti,gcurrenti_end,currenteWord));
                        }
                        
                    }
                }
            }
            
        }
    }

    cout << "Done with creating Higher N block" << endl;
    
    return B;

}

void PrintFullGraphInfo(Graph G,ostream& os)
{
    bool found;
    
    graph_traits<Graph>::vertex_iterator gvi,gvi_end,gwi,gwi_end;
    graph_traits<Graph>::out_edge_iterator goei,goei_end;

    property_map<Graph,edge_name_t>::type
    gename = get(edge_name,G);
    
    property_map<Graph,vertex_name_t>::type
    gvname = get(vertex_name,G);
    
    os << "Printing G Graph which has " << num_vertices(G) << " vertices and " << num_edges(G) << " edges." << endl;
    
    
    for(tie(gvi,gvi_end)=vertices(G); gvi!=gvi_end; gvi++)
    {
        os << gvname(*gvi) << " ";
        map<graph_traits<Graph>::vertex_descriptor,string> strStor;
        for(tie(goei,goei_end)=out_edges(*gvi,G);goei!=goei_end;goei++)
        {
            graph_traits<Graph>::vertex_descriptor t=target(*goei,G);
            if(strStor.find(t)==strStor.end()) {
                strStor[t]=gename(*goei);
            }
            else {
                strStor[t]=strStor[t]+ string(" + ") + gename(*goei);
            }
            
        }
        
        for(tie(gwi,gwi_end)=vertices(G); gwi!=gwi_end; gwi++)
        {

            if(strStor.find(*gwi)!=strStor.end()) {
                os << strStor[*gwi] << " ";
            }
            else {
                os << "e ";
            }
        }
        os << endl;
        
    }
    
}


Graph Rename(Graph G, vector<string> names)
{
    int i;
    Graph H = G;

    graph_traits<Graph>::edge_iterator ei,ei_end;

    property_map<Graph,edge_name_t>::type
    ename = get(edge_name,H);

    for(tie(ei,ei_end)=edges(H), i=0;ei!=ei_end;ei++,i++) {
        put(ename,*ei,names[i]);
    }

    return H;

}

vector<string> RandomStringGenerator(int n,double p)
{
    vector<string> rstring(n);
    int i;
    ostringstream ss;


    for(i=0;i<n;i++)
    {
        ss.str(std::string());
        ss << flip(p);
       // cout << flip();
        rstring[i] = ss.str();
    }
    return rstring;
}

int flip(double p) {
    double probs[] = {p, 1.0-p};
    boost::random::discrete_distribution<> dist(probs);
    return dist(gen);
}

VVec compatibleSet(Graph & G, VVec U, string s)
{
    VVec C;
    
    OEI oei,oei_end;
    
    vector<graph_traits<Graph>::vertex_descriptor>::iterator uit;
    
    property_map<Graph,edge_name_t>::type
    ename = get(edge_name,G);

    typedef tuple<graph_traits<Graph>::vertex_descriptor,string> CSTuple;
    
    stack<CSTuple> CSStack;
    bool found;
    
    //cout << "We are computing S+(U," << s << ")" << endl;
    
    for(uit = U.begin(); uit < U.end(); uit++)
    {
        CSStack.push(CSTuple(*uit,s));
    }
    
    while(!CSStack.empty())
    {
        CSTuple current = CSStack.top();
        CSStack.pop();
        
        graph_traits<Graph>::vertex_descriptor v = get<0>(current);
        string currs = get<1>(current);
        
        
        OEI ni,ni_end;
        
        for(tie(ni,ni_end)=out_edges(v,G);ni!=ni_end;ni++)
        {
            int len = (ename(*ni)).length(),lens = currs.length();
            
            //    cout << "p(edge) = " << pt(*ni) << "and our s is " << currs << endl;
            
            
            if(ename(*ni) == string(currs,0,len) && ename(*ni) != currs)
            {
                //   cout << "We are pushing (" << target(*ni,T.first) <<  "," << string(currs,len,lens) << ")" << endl;
                CSStack.push(CSTuple(target(*ni,G),string(currs,len,lens)));
            }
            else if(ename(*ni)==currs)
            {
                found = false;
                for(uit = C.begin(); uit < C.end(); uit++)
                {
                    if(*uit == target(*ni,G))
                    {
                        found = true;
                        
                    }
                }
                if(!found)
                {
                    C.push_back(target(*ni,G));
                    
                }
                
            }
        }
    }
    
    sort(C.begin(),C.end());
    
    return C;
    
}

Graph inducedRp(Graph & G)
{
    VI vi,vi_end;
    OEI oei,oei_end;
    EI ei,ei_end;
    
    Graph Gamma;
    

    property_map<Graph,edge_name_t>::type
    gename = get(edge_name,G);
    
    property_map<Graph,edge_name_t>::type
    gamma_ename = get(edge_name,Gamma);
    
    property_map<Graph,vertex_name_t>::type
    gamma_vname = get(vertex_name,Gamma);
    
    property_map<Graph,vertex_name_t>::type
    tfvname = get(vertex_name,G);
    
    stack<VVec> toInv;
    
    set<VVec> seen;
    set<string> codex;
    
    for(tie(vi,vi_end)=vertices(G);vi!=vi_end;vi++)
    {
        VVec v(1,*vi);
        toInv.push(v);
        seen.insert(v);
    }

    for(tie(ei,ei_end)=edges(G);ei!=ei_end;ei++)
    {
        string s = gename(*ei);
        codex.insert(s);
    }

    while(!toInv.empty())
    {
        set<string>::iterator sit;
        VVec::iterator it;
        VVec v = toInv.top();
        toInv.pop();

    /*  cout << "Look at {";

        for(it = v.begin(); it < v.end(); it++)
        {
            cout << *it << " ";
        }

        cout << "}" << endl; */

        for(sit=codex.begin(); sit != codex.end(); sit++)
        {
            VVec S = compatibleSet(G,v,*sit);
            if(S.size() > 0)
            {
    /*          cout << "S is of size " << S.size() << endl << "It is {";
                for(it = S.begin(); it < S.end(); it++)
                {
                    cout << *it << " ";
                }
                cout << "}" << endl; */

                if(seen.find(S)==seen.end())
                {
                    seen.insert(S);
                    toInv.push(S);
                }
                else{
        //          cout << "We've already seen it, not inserting." << endl;
                }
            }
        }

    }

    
    
    // Before we create the actual graph, we need to trim our seen set. We need to check if we have non-maximal
    // VVecs in there and if so, get rid of them.
    set<VVec>::iterator si,ti;
    set<VVec> toDelete;
    set<VVec>::iterator tdi;
    bool done=false,fdone=false,found=false;
    VertexMap vmap;
    
    for(si=seen.begin(); si!=seen.end(); si++)
    {
        for(ti=seen.begin(); ti!=seen.end(); ti++)
        {
            VVec a=*si,b=*ti;
            if(VVecSubset(a,b) && toDelete.find(a) == toDelete.end())
            {
             /**   
                cout << "We are planning to erase";
                printVVec(a); cout << endl;
                toDelete.insert(a);
                */
            }
        }
    }

  /*  cout << "We have " << toDelete.size() << " vectors to delete" << endl;
    
    for(tdi=toDelete.begin(); tdi!=toDelete.end(); tdi++)
    {
         //   cout << "We are actually erasing "; printVVec(*tdi); cout << endl;
        seen.erase(*tdi);
    }
   */ 
    
    cout << "Our final Seen consists of:";
    
    for(si=seen.begin(); si != seen.end(); si++)
    {

        graph_traits<Graph>::vertex_descriptor v = add_vertex(Gamma);
        VVec::iterator it;
        string vname;
        VVec curr=*si;
        
        printVVec(curr);
        vname = ssVVec(curr);
        put(gamma_vname,v,vname);
        
        //We need to create a VMap so we can add edges later
        vmap.insert(VertexMap::value_type(vname,v));
    }
    cout << endl;
    
    for(si=seen.begin(); si != seen.end(); si++)
    {
        string s=ssVVec(*si),t;
        graph_traits<Graph>::vertex_descriptor v=vmap[s],w;
        set<string>::iterator cit;
        VVec S;
        
        for(cit = codex.begin(); cit != codex.end(); cit++)
        {
            S = compatibleSet(G,(*si),(*cit));
            if(S.size() > 0) {
                t = ssVVec(S);
                w = vmap[t];
                add_edge(v,w,*cit,Gamma);
            }  
        }
        
        
    }
    
    /*  while(!toInv.empty())
     {
     VVec::iterator it;
     VVec v = toInv.top();
     for(it = v.begin(); it < v.end(); it++)
     {
     cout << *it;
     }
     cout << endl;
     toInv.pop();
     }
     */

     return Gamma;
 }



// This function finds the Perron Frobenius Eigenvalue of the input matrix. We use JAMA for this as
// I would definitely screw it up.
 /**double PFEigenvalue(Graph & G)
 {
    // We use the graph G to create an adjacency matrix, then we make a JAMA Eigen class for it
    int N = num_vertices(G),i;
    Array2D<double> Adj(N,N,0.0);
    double max = 0;
    EI gei,gei_end;
    
    for(tie(gei,gei_end)=edges(G); gei!=gei_end; gei++)
    {
        Adj[source(*gei,G)][target(*gei,G)]+=1.0;
    }
    
    Eigenvalue<double> E(Adj);
    
    Array1D<double> Eig; 
    E.getRealEigenvalues(Eig);
    
    for(i=0; i<N; i++)
    {
      //  cout << Eig[i] << endl;
        if(Eig[i] > max)
        {
            max = Eig[i];
        }
    }
    
    return max;
    
}*/


// My attempt at the PF eigenvalue calculation. 
 double MyPFEigenvalue(Graph & G)
 {
    // We use the graph G to create an adjacency matrix, then we make a JAMA Eigen class for it
    int N = num_vertices(G),i,j=1;
    SpMat<double> Adj(N,N);
    EI gei,gei_end;
    double tol=0.00001,eval;
    
    for(tie(gei,gei_end)=edges(G); gei!=gei_end; gei++)
    {
        Adj(source(*gei,G),target(*gei,G))+=1.0;
    }

    vec Eig(N);
    vec OldEig(N);
    vec Diff(N);
    Eig.fill(1.0); 
    
    do
    {
        OldEig = Eig;
        Eig = Adj*Eig;

        Eig = Eig / norm(Eig,2);

        Diff = OldEig-Eig;


        for(i=0; i<N; i++)
        {
//        cout << Eig(i) << " ";

        }
 //   cout << endl;
//        j++;    
    } while(norm(Diff,2)>tol);
    eval = norm(Adj*Eig,2) / norm(Eig,2);

    return eval;
    
}

int oneCounter(vector<string> s)
{
    int i,numOnes = 0,n=s.size();
    for(i=0; i<n;i++)
    {
        if(s[i]==string("1"))
            numOnes++;
    }
    return numOnes;
}

Graph Trim(Graph &G)
{
    Graph Trimmed=G;
    bool isEndIT=true,isEndTI=true,isEndPQ=true,isEndQP=true;
    unsigned int i,j=0,delcount=0,N;
    //  vector<graph_traits<GammaGraph>::edge_descriptor> toDelete;
    vector<int> vertToDelete;
    bool done,stable;
    time_t start,end;
    double dif,average=0;

    graph_traits<Graph>::vertex_iterator vi,vi_end;

    bool needToCheckVhoms=false,found=false,topStart=true;;
    
   // stack<VD> toCheck;  

    cout << "Our starting graph has " << num_edges(Trimmed) << " edges and " << num_vertices(Trimmed) << " vertices." << endl; 

    time(&start);
    do{
        j++;
        if(j % 100 == 0)
        {
         time(&end);
         dif = difftime(end,start);
         time(&start);
         N = j/100;
         average = (((N-1)*average)+dif)/N;
         cout << "On round " << j << " deleted at least " << delcount << " edges so far. The last round took " 
                << dif << " seconds, and the average time is " << average << " seconds." <<  '\r';
         cout.flush();
     } 
     stable = true;
     done = false;

        // Check for sinks 
     for(tie(vi,vi_end)=vertices(Trimmed); vi!=vi_end && !done; vi++)
     {
            //  cout << "We are checking " << vname(*vi) << endl;
        if(out_degree(*vi,Trimmed)==0 || in_degree(*vi,Trimmed)==0)
        {
                //   cout << "The vertex " << vname(*vi) << "is a sink or a source" << endl;
            clear_vertex(*vi,Trimmed);
            remove_vertex(*vi,Trimmed);
            stable = false;
            done = true;
            break;
        }
            //  cout << "It's not." << endl;
    }
    } while(!stable); // do while statement (don't panic)
    
    
    cout << "Original T.first had " << num_edges(G) << " edges and " << num_vertices(G)
    << " vertices, now has " << num_edges(Trimmed) << " edges and " << num_vertices(Trimmed) << endl;

    return Trimmed;
}

double standardDeviation(vector<double> &v,double mean,int n)
{
    double stdev=0;
    int i;

    for(i=0;i<v.size();i++)
    {
        stdev += pow((v[i]-mean),2.0);
    }

    stdev = stdev/n;

    return sqrt(stdev);

}
