#include "caenthelper.h"

using namespace std;

int main(void)
{

	Graph G(1);
	
	property_map<Graph,vertex_name_t>::type
	gvname = get(vertex_name,G);
    vector<string> s;
    int i,j,k=2;
    ofstream outfile("CAENToutput5.txt");
    double ceval,deval,seval,total=0,runs=10,stdev,p;
    bool equalZerosAndOnes = 0;
    clock_t start,end;
    vector<double> evals(runs);

    outfile << setprecision(15);
    cout << setprecision(15);
/**
	add_edge(0,0,string("00"),G);
	add_edge(0,1,string("01"),G);
	add_edge(1,0,string("10"),G);
	add_edge(1,1,string("11"),G);

	put(gvname,0,string("0"));
	put(gvname,1,string("1"));
**/
    add_edge(0,0,string("0"),G);
    add_edge(0,0,string("1"),G);

//    PrintFullGraphInfo(G);

    k=6;
    p=0.5;
    while(k<7)
    {
        Graph H = HigherNBlock(G,k);

//    PrintFullGraphInfo(H);
        for(j=0; j<runs; j++)
        {
            if(equalZerosAndOnes) {
                do{
                    s = RandomStringGenerator(num_edges(H),p);
           //     cout << oneCounter(s) << endl;
                } while(oneCounter(s)!=num_edges(H)/2);
            }
            else{
                s = RandomStringGenerator(num_edges(H),p);
            }

            for(i=0; i<s.size();i++)
            {
                cout << s[i];
            //    outfile << s[i];
            }
            cout <<  "," << oneCounter(s) << endl;
            //outfile << "," << oneCounter(s) << ",";

            Graph J = Rename(H,s);

    //        PrintFullGraphInfo(J);

   // PrintFullGraphInfo(J);

            start = time(0);
            Graph K = inducedRp(J);
    //        PrintFullGraphInfo(K);

            end = time(0);
            cout << "Right resolving generation took " << (double) end - start << " seconds" << endl;

/**
        start = time(0);
        ceval = PFEigenvalue(K);
        end = time(0);
        cout << "Eigenvalue calculation took " << (double) end - start << " seconds " << endl;

        cout << ceval << endl;
        outfile << ceval << endl;
**/

        start = time(0);
        Graph L = Trim(K);
//        deval = PFEigenvalue(L);
        seval = MyPFEigenvalue(L);
        end = time(0);
        cout << "Trim and eval calculation took " << (double) end - start << " seconds " << endl;

        cout << seval << endl;

        evals[j] = seval;
        total += seval;

     //   cout << deval << endl;
    //    outfile << deval << endl;

    }
    total = total / runs;
    outfile << k << ", " << p << ", " << runs << ", " <<  total << ", " << standardDeviation(evals,total,runs) << endl;
 //   p += 0.01;
    k++;
    total = 0;

}
//    PrintFullGraphInfo(K);


/**
Graph A(2);
    
    property_map<Graph,vertex_name_t>::type
    avname = get(vertex_name,A);

    add_edge(0,0,string("0"),A);
    add_edge(0,1,string("0"),A);
    add_edge(1,0,string("1"),A);

    put(avname,0,string("a"));
    put(avname,1,string("b"));

    PrintFullGraphInfo(A);

    Graph B = inducedRp(A);

    PrintFullGraphInfo(B);
*/
    
    outfile.close();

    return 0;

}