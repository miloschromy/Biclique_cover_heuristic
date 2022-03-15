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
		int clauses = variables*0.6;
		for (; clauses <= variables; clauses+=variables/100){
			for (int kknf = 3 ;kknf <= 10; kknf+=1){
				int failed_forced_SAT = 0, matched_rand_SAT = 0, matched_kcnf_SAT = 0;
							
				for (int iteration = 0; iteration < iterations_count;++iteration){
					vector<vector<int> >adjacency_matrix_ksat, adjacency_matrix_prob;
				
					//vyggenerujeme incidenční grafy
					generate_random_k_sat(adjacency_matrix_ksat, variables, clauses, kknf);
					generate_random_bipartite_graph(adjacency_matrix_prob, variables, clauses, (double)kknf / (double) variables, 0.0);
					
					if (forced_matching(adjacency_matrix_prob, variables, clauses)){
						++failed_forced_SAT;
						matched_rand_SAT += test_matching(adjacency_matrix_prob, variables, clauses);
					}
					matched_kcnf_SAT+=test_matching(adjacency_matrix_ksat, variables, clauses);
					
				}
			
				//ladící výpis
				cerr << "Variables: " << variables << ". Clauses: " << clauses << ". Degree: " << kknf << 
					". Matched k-CNF: " << matched_kcnf_SAT << ". Failed unitGProp: " << failed_forced_SAT << ". Matched k random KNF "  << matched_rand_SAT <<
					"." << endl; 
				//výpis pro analýzu
				cout << variables << "," << clauses << "," << kknf << "," << 
					 matched_kcnf_SAT << "," << failed_forced_SAT << ","  << matched_rand_SAT << endl; 
			}
				
		}
	}
	
	return 0;
}
