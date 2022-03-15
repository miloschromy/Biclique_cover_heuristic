// Autor Miloš chromý. Provede testy na biklikové pokrytí bipartitního grafu. 


#include "functions.h"
#include <iostream>

using namespace std;
int main(int argc, char** argv){
	int iterations_count = 100;
	vector<int> vararr; //různé počty vrcholů.
	vararr.push_back(100);
	vararr.push_back(200);
	vararr.push_back(500);
	vararr.push_back(1000);
	for (int var : vararr){
		int variables = var;
		int clauses = variables;
		for (; clauses <= variables*1.5; clauses+=variables/100){
			for (int kknf = 3 ;kknf <= variables; kknf+=1){
				int biclique_2_min = 0, biclique_2_rand = 0, biclique_2_max = 0, biclique_unbound_min = 0, biclique_unbound_rand = 0, biclique_unbound_max = 0;
							
				for (int iteration = 0; iteration < iterations_count;++iteration){
					vector<vector<int> >adjacency_matrix;
				
					//vyggenerujeme incidenční graf 
					generate_random_k_sat(adjacency_matrix, variables, clauses, kknf);
					
					//tety pro biklkiky
					biclique_2_max += test_bipartite_h(adjacency_matrix, variables, clauses, 2,-1);
					biclique_2_rand += test_bipartite_h(adjacency_matrix, variables, clauses, 2,0);
					biclique_2_min += test_bipartite_h(adjacency_matrix, variables, clauses, 2,1);
					biclique_unbound_max += test_bipartite_h(adjacency_matrix, variables, clauses, variables,-1);
					biclique_unbound_rand += test_bipartite_h(adjacency_matrix, variables, clauses, variables,0);
					biclique_unbound_min += test_bipartite_h(adjacency_matrix, variables, clauses, variables,1);
					
				}
			
				//ladící výpis
				cerr << "Variables: " << variables << ". Clauses: " << clauses << ". k-KNF: " << kknf << 
					". K23_SAT(max,rand,min): " << biclique_2_max << ", " << biclique_2_rand << ", "  << biclique_2_min <<
					". Kmax_SAT(max,rand,min): " << biclique_unbound_max << ", " << biclique_unbound_rand << ", "  << biclique_unbound_min << endl; 
				//výpis pro analýzu
				cout << variables << "," << clauses << "," << kknf << "," << 
					 biclique_2_max << "," << biclique_2_rand << ","  << biclique_2_min << "," << 
					 biclique_unbound_max << "," << biclique_unbound_rand << ","  << biclique_unbound_min << endl; 
			}
		}
	}
	
	return 0;
}
