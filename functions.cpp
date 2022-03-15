#include "functions.h"

#include <cassert>
#include <unordered_set>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <random>
#include <queue>
using namespace std;

/**
 * Struktura reprezentující bikliku
 */
struct biclique_t{
	biclique_t() {id = -1;}
	vector<int> vertex_partition; //velikost první partity
	vector<int> clauses_partition; //velikost druhé parity
	int id;
	int vertex_count() const {return vertex_partition.size();}
	int clauses_count() const {return clauses_partition.size();}
};

/**
 * Funkce pro odstranění hodnty z matice sousednosti. Pokud je seznam setřízený naruší jeho strukturu.
 */
int remove_value_from_neighbours(vector<vector<int> >& adjacency_matrix, const int& vertex, const int& value){
	for (size_t i = 0; i < adjacency_matrix[vertex].size(); ++i){
		if (adjacency_matrix[vertex][i] == value || adjacency_matrix[vertex][i] == -value){ 
			adjacency_matrix[vertex][i] = adjacency_matrix[vertex][adjacency_matrix[vertex].size()-1];
			adjacency_matrix[vertex].resize(adjacency_matrix[vertex].size()-1);
			return adjacency_matrix[vertex].size();
		}
	}
	return adjacency_matrix[vertex].size();
}

/**
 *Naleze zlepšující cestu v grafu.
 *Matching reprezentuje párování.
 */
bool find_augumenting_path(vector<vector<int>>& adjacency_matrix, const int& variables, const int& clauses, vector<int>& matching, list<int>& unpaired){
	queue<int> BFS_queue;
	vector<int> BFS_ancestor(variables+1, 0);
	int c, var, c_swap;
	for (int v : unpaired){
		BFS_ancestor[v] = -1;
		BFS_queue.push(v);
	}
	while(!BFS_queue.empty()){
		var = BFS_queue.front();
		BFS_queue.pop();
		for (int c_signed : adjacency_matrix[var]){
			c = c_signed > 0 ? c_signed : -c_signed;
			if (matching[c] != 0){
				if (BFS_ancestor[matching[c]] == 0){
					BFS_ancestor[matching[c]] = var;
					BFS_queue.push(matching[c]);
				}
			} else {
				matching[c] = var;
				while (BFS_ancestor[var] != -1){
					c_swap = c;
					c = matching[var];
					matching[c] = BFS_ancestor[var];
					matching[var] = c_swap;
					var = BFS_ancestor[var];
				}
				matching[var] = c;
				unpaired.remove(var);
				return true;
			}
		}
	}
	
return false;

}

/** 
 * Párování vytiskne na výstup.
 */
void print_matching(vector<int>& matching, int& variables){
	for (int i = 1 ; i <= variables; ++i){
		cout << i << " - " << matching[i] << "; ";
	}
	cout << endl;
}

/**
 * Smaže množinu R_{G,M} z grafu. Redukce R_{G,M}
 */
bool remove_RgmSet(vector<vector<int>>& adjacency_matrix,const int& variables, const int& clauses, vector<int>& matching, vector<biclique_t>& removed_matching, const bool& store_removed){
	queue<int> BFS_queue;
	vector<bool> visited_verticies(adjacency_matrix.size(), false);
	unordered_set<int> clauses_SAT;
	int cleared_v = 0, cleared_c = 0, uncovered_v = 0, uncovered_c = 0;
	//nepokryté klauzule dáme do fronty
	for (int c = variables+1; c < adjacency_matrix.size(); ++c){
		if (adjacency_matrix[c].size() > 0 &&  matching[c] == 0){
			BFS_queue.push(c);
			visited_verticies[c]=true;
		}
	}
	// všechny vrcholy dosažitelné alternující cestou vložíme označíme
	while (!BFS_queue.empty()){
		int clause = BFS_queue.front();
		BFS_queue.pop();
		for (int v : adjacency_matrix[clause]){
			v = v > 0 ? v : -v;
			if (visited_verticies[v]) continue;
			visited_verticies[v] = true;
			assert(matching[v] != 0);
			visited_verticies[matching[v]] = true;
			BFS_queue.push(matching[v]);
		}
	}
	removed_matching.push_back(biclique_t());
	// smažeme všechny vrcholy nedosažitelné alternující cestou z vrcholů nepokrytých párováním
	for (int v = 1; v <=variables; ++v){
		if (!visited_verticies[v] && adjacency_matrix[v].size() > 0){
			if (matching[v] == 0){ 
				++uncovered_v;
			} else if (matching[v]  != 0 && store_removed){
				biclique_t new_matching_edge = biclique_t();
				new_matching_edge.vertex_partition.push_back(v);
				new_matching_edge.clauses_partition.push_back(matching[v]);
			}
			if (store_removed){
				removed_matching.back().vertex_partition.push_back(v);
			}
			adjacency_matrix[v].clear();
			matching[v] = 0;
			++cleared_v;
		}
	}
	for (int c = variables+1; c < adjacency_matrix.size(); ++c){
		if (!visited_verticies[c] && adjacency_matrix[c].size() > 0){
			if (matching[c] == 0) ++uncovered_c;
			for (int v : adjacency_matrix[c]){
				v = v > 0 ? v : -v;
				if (adjacency_matrix[v].size() > 0){
					remove_value_from_neighbours(adjacency_matrix,v, c);
				}
			}
			if (store_removed){
				removed_matching.back().clauses_partition.push_back(c);
			}
			adjacency_matrix[c].clear();
			matching[c] = 0;
			++cleared_c;
		}
	}
	if (removed_matching.back().clauses_count() == 0) {
		removed_matching.pop_back();
	}
	return true;
}

/**
 * Funkce, která najde maximální párování. V heuristice je to funkce testMatched(G).
 */
size_t find_max_matching(vector<vector<int>>& adjacency_matrix,const int& variables, const int& clauses, vector<biclique_t>& removed_matching,const bool& store_removed){
	vector<int> matching(adjacency_matrix.size(), 0);
	list<int> unpaired;
	int nonempty_v = count_if(adjacency_matrix.begin(), adjacency_matrix.begin()+variables+1, [](vector<int>& v){return v.size() > 0; });
	int	nonempty_c = count_if(adjacency_matrix.begin()+variables+1, adjacency_matrix.end(), [](vector<int>& v){return v.size() > 0; });
	if (nonempty_v<nonempty_c) {return 0;}
	for (int i = 1; i <= variables; ++i){
		if (!adjacency_matrix[i].empty())
			unpaired.push_back(i);
	}
	while (find_augumenting_path(adjacency_matrix, variables, clauses, matching, unpaired)){
	
	}
	
	int matching_size = count_if(matching.begin()+1, matching.begin()+1+variables, [](int& m){ return m != 0; });
	int clauses_size = count_if(adjacency_matrix.begin()+variables+1, adjacency_matrix.end(), [](vector<int>& neighbours){ return neighbours.size() > 0;});
	remove_RgmSet(adjacency_matrix,variables,clauses,matching, removed_matching, store_removed);

	int reduced_clauses_size = count_if(adjacency_matrix.begin()+variables+1, adjacency_matrix.end(), [](vector<int>& neighbours){ return neighbours.size() > 0;});
	if (matching_size == clauses_size || reduced_clauses_size == 0){
		assert(matching_size == clauses_size && reduced_clauses_size == 0);
	}
	return matching_size;
}
size_t find_max_matching(vector<vector<int>>& adjacency_matrix,const int& variables, const int& clauses){
	vector<biclique_t> tmp;
	return find_max_matching(adjacency_matrix, variables, clauses, tmp, false);
}

/**
 * Funkce pro načtení knf ze souboru, jež odpovídá specifikaci knf souboru v sat comptetiton (http://baldur.iti.kit.edu/sat2014/competitions.html). 
 */
void load_formula(vector<vector<int> >& adjacency_matrix, int& variables, int& clauses, istream& input){
	char c;
	string line = "";
	input >> c;
	while (c=='c'){
		getline(input, line);
		input >> c;
	}
	input >> line;
	int cur_var;
	input >> variables;
	input >> clauses;
	adjacency_matrix.resize(variables+clauses+1);
	for(int i = variables+1; i <= variables+clauses; ++i){
		input >> cur_var;
		while (cur_var != 0){
			if (cur_var > 0){
				adjacency_matrix[cur_var].push_back(i);
				adjacency_matrix[i].push_back(cur_var);
			} else {
				adjacency_matrix[-cur_var].push_back(-i);
				adjacency_matrix[i].push_back(cur_var);
			}
			input >> cur_var;
			
		}

	}
}

/**
 * Vytiskne formuli, jejíž incidenční graf je adjacency_matrix na výstup output.
 */
void print_formula(const vector<vector<int> >& adjacency_matrix,const int& variables, const int&  clauses, ostream& output){
	for (int i = 1; i <= clauses+variables; ++i){
		output << i << ": ";
		for (int j : adjacency_matrix[i]) output << j << ", ";
		output << endl;
	}
	for (int i = variables+1; i <= clauses+variables; ++i){
		if (adjacency_matrix[i].size() > 0) continue;
		output << "(";
		if (adjacency_matrix[i][0] == 0) output << "FALSE";
		else for (int j : adjacency_matrix[i]) output << j << " ";
		output << ")";
		if (i != clauses+variables) output << " AND ";
		else output << endl;
	}

}

/**
 * Zapíše formuli do souboru s názvem outfile_name
 */
void print_formula_file(vector<vector<int> >& adjacency_matrix,const int& variables, const int&  clauses, const char* outfile_name){
	ofstream outfile;
	outfile.open(outfile_name);
	print_formula(adjacency_matrix, variables, clauses, outfile);
	outfile.close();
}

/**
 * Jendotková propagace na klauzulích.
 */
bool unit_propagation(vector<vector<int> >& adjacency_matrix,const int& variables, const int&  clauses){
	list<int> unit_clauses;
	int propagation_clause, propagation_variable, return_size;
	int polarity;
	//get init queue of unit clauses
	for (int i = variables+1; i <= clauses+variables; ++i){
		if (adjacency_matrix[i].size() == 1){
				unit_clauses.push_back(i);
		}
	}
	//while there exists unit clause
	while (!unit_clauses.empty()){
		propagation_clause = unit_clauses.front(); //unit clause
		propagation_variable = adjacency_matrix[propagation_clause][0]; //only var in unit clause and its polarity
		polarity = adjacency_matrix[propagation_clause][0] > 0 ? 1 : -1;
		adjacency_matrix[propagation_clause].clear(); // clear me
		unit_clauses.pop_front();
		for (int c : adjacency_matrix[polarity*propagation_variable]){
			if (c*polarity > 0){
				for (int var : adjacency_matrix[c*polarity]){
					if (var != propagation_variable)
						remove_value_from_neighbours(adjacency_matrix, var > 0 ? var : -var,c);
				}
				adjacency_matrix[c*polarity].clear();
				//clause is true, we can delete it and remove it from vars
			} else {
				return_size=remove_value_from_neighbours(adjacency_matrix, -polarity*c, -1*propagation_variable);
				if (return_size == 0) return false;
				else if (return_size == 1) unit_clauses.push_back(-polarity*c);
				//we have to delete just one literal and check if we have to add to queue
			}
		}
		adjacency_matrix[polarity*propagation_variable].clear();
	}
	return true;
}

/**
 * Kontrola, zdali je literál check_variable čistý.
 */
bool is_pure(vector<vector<int> >& adjacency_matrix, const int& check_variable){
	if (adjacency_matrix[check_variable].empty()) return false; // empty literal is true, but its useless to mark it as one
	int polarity = adjacency_matrix[check_variable][0] > 0 ? 1 : -1; // polarity of first occurence of literal
	for (int c : adjacency_matrix[check_variable]){
		if (c*polarity < 0 ){
			return false;
		}
	}
	return true;
}

/**
 * Aplikace redukce čistých literálů.
 */
bool pure_literals(vector<vector<int> >& adjacency_matrix,const int& variables, const int&  clauses){
	list<int> pure_literals;
	int polarity;
	//nejdriv asi spis sezbiram pure literaly a pak je projedu poradne
	for (int i = 1; i <= variables; ++i){
		if (adjacency_matrix[i].size() >= 1){
			if (is_pure(adjacency_matrix, i)){
				pure_literals.push_back(i);
			}
		}
	}
	while (!pure_literals.empty()){
		int variable = pure_literals.front();
		pure_literals.pop_front();
		if (adjacency_matrix[variable].empty()) continue;
		polarity = adjacency_matrix[variable][0] > 0 ? 1 : -1;
		for (int c : adjacency_matrix[variable]){
			for (int v : adjacency_matrix[c*polarity]){
				v = v > 0 ? v : -v;
				if (v != variable){
					remove_value_from_neighbours(adjacency_matrix, v, c);
					if (is_pure(adjacency_matrix, v)){
						pure_literals.push_back(v);
					}
				}
			}
			adjacency_matrix[c*polarity].clear();
		}
		adjacency_matrix[variable].clear();
	}
	return true;
}


/**
 * Funkce z heuristiky rozšiřBikliku. Pokud by rozšíření o největší vrchol vytvořilo nesplnitelnou formuli vrací false a nerozšíří ji.
 */
bool enhance_biclique(const vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, biclique_t& biclique){
	vector<bool> subsequence_mark(variables+clauses+1, false);
	bool only_small = true;
	int reduced_size, max_size = -1, max_vertex = -1;
	for (int i : biclique.vertex_partition){
		subsequence_mark[i] = true;	
	}
	for (int i : biclique.clauses_partition){
		subsequence_mark[i] = true;	
	}
	for (int i = 1; i <= variables; ++i){ // pokusíme se přidat formuli, zastupuje funkci zvolRozšiřujícíVrchol
		if (subsequence_mark[i]) continue; // proměnná i je v biklice
		reduced_size = 0;
		for (int c : adjacency_matrix[i]){ //zkusíme vypočítat velikost bikliky
			c = c > 0 ? c : -c;
			reduced_size+=subsequence_mark[c];
			if (adjacency_matrix[c].size() < biclique.vertex_count()+1){ // test, jestli nevytvoříme prázdnou klauzuli
				only_small = true;
				for (int v_c : adjacency_matrix[c]){
					v_c = v_c > 0 ? v_c : -v_c;
					if (!subsequence_mark[v_c]){
						only_small = false;
						break;
					}
				}
				if (only_small){ //vytvořili jsme prázdnou klauzuli
					break;
				}
			} else {only_small = false;}
		}
		if (reduced_size > max_size && !only_small){ //zatím nelepší rozšíření
			max_size = reduced_size;
			max_vertex=i;
		}
	}
	if (max_size >= pow(2,biclique.vertex_count())-1){ //pokud jsme našli proměnnou, jež dostatečně rozšiřuje použijeme jí
		biclique.clauses_partition.clear(); 
		biclique.vertex_partition.push_back(max_vertex); // vyměníme klauzule v biklice
		
		for (int c : adjacency_matrix[max_vertex]){ 
			c = c > 0 ? c : -c;
			if (subsequence_mark[c]) {
				biclique.clauses_partition.push_back(c);
			}
		}
		return true;
	}
	return false;

}


/**
 * Funkce vygenerujZárodky. Pro každou dvojici vygeneruje zárodek, tedy bikliku K_2,?
 */
bool make_biclique_seed(const vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, vector<biclique_t>& bicliques_seed){
	vector<int> subsequence_mark(variables+clauses+1,0);
	int pair_number = 1;
	for (int i = 1; i <= variables; ++i){
		for (int c : adjacency_matrix[i]){
			subsequence_mark[c > 0 ? c : -c] = i;
		}
		for (int j = i+1; j <= variables; ++j){ //every tuple i,j appears exeactly once i<j
			bicliques_seed.push_back(biclique_t());
			
			bicliques_seed.back().vertex_partition.push_back(i);
			bicliques_seed.back().vertex_partition.push_back(j);
			
			for (int c : adjacency_matrix[j]){
				if (subsequence_mark[c > 0 ? c : -c] == i){
					bicliques_seed.back().clauses_partition.push_back(c > 0 ? c : -c);
				}
			}
			if (bicliques_seed.back().clauses_count() < 3) bicliques_seed.pop_back(); // biclique seed is worse or equal matching
		}
	}
	return bicliques_seed.size() > 0;
}

/**
 * Funkce vyčistiZárodky
 */
bool clear_biclique_seeds(const vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, vector<biclique_t>& bicliques_seed){
	vector<int> biclique_incidence(variables+clauses+1, -1);
	for (int i = 0 ; i < bicliques_seed.size(); ++i){
		bicliques_seed[i].clauses_partition.clear();
		for (int j : adjacency_matrix[bicliques_seed[i].vertex_partition[0]]){
			biclique_incidence[j > 0 ? j : -j] = i;
		}
		for (int j : adjacency_matrix[bicliques_seed[i].vertex_partition[1]]){
			j = j > 0 ? j : -j;
			if (biclique_incidence[j] == i){
				bicliques_seed[i].clauses_partition.push_back(j);
			}
		}
	}
	bicliques_seed.erase(
	remove_if(bicliques_seed.begin(), bicliques_seed.end(), [&](const biclique_t& biclique){ // všechny zárodky s malou druhou partitou
		return (biclique.clauses_count() < 3);
	}), bicliques_seed.end());
	return false;
}

/*
 * Funkce jednotkováGPropagace. Pročistí graf od vrcholů z druhé partity stupně 1 tak, že hranu jim incidentní přidá do pokrytí.
 */
bool forced_matching(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, vector<biclique_t>& used_matching, bool store_matching){
	int count = count_if(adjacency_matrix.begin(), adjacency_matrix.end(), [](vector<int>& vec){ return vec.size() != 0;});
	vector<bool> variables_used(variables+1,false);
	list<int> queue;
	int removed_verticies = 0;
	for (int i = variables+1; i <= clauses+variables; ++i){
		if (adjacency_matrix[i].size() == 1){
			int v = adjacency_matrix[i][0] > 0 ? adjacency_matrix[i][0] : -adjacency_matrix[i][0];
			if (variables_used[v]){ return false; }
			variables_used[v] = true;
			queue.push_back(v);
			remove_value_from_neighbours(adjacency_matrix,v, i); 
			adjacency_matrix[i].clear();
			if (store_matching){
				biclique_t new_matching_edge = biclique_t();
				new_matching_edge.vertex_partition.push_back(v);
				new_matching_edge.clauses_partition.push_back(i);
			}
		}
	}
	while (!queue.empty()){
		int front = queue.front();
		queue.pop_front();
		++removed_verticies;
		for (int c : adjacency_matrix[front]){
			c = c > 0 ? c : -c;
			remove_value_from_neighbours(adjacency_matrix, c, front);
			if (adjacency_matrix[c].size() == 1){ 
				int v = adjacency_matrix[c][0] > 0 ? adjacency_matrix[c][0] : -adjacency_matrix[c][0];
				if (variables_used[v]){ return false; }
				variables_used[v] = true;
				queue.push_back(v);
				remove_value_from_neighbours(adjacency_matrix, v, c);
				adjacency_matrix[c].clear();
				if (store_matching){
					biclique_t new_matching_edge = biclique_t();
					new_matching_edge.vertex_partition.push_back(v);
					new_matching_edge.clauses_partition.push_back(c);
				}
			}
			else if (adjacency_matrix[c].empty()) return false;
		}
		adjacency_matrix[front].clear();
		
	}
	return true;
}
/**
 * Funkce jednotkovaGPropagace, která si nepamatuje pokrytí, takže se jen ptáme na existenci.
 */
bool forced_matching(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses){
	vector<biclique_t> tmp;
	return forced_matching(adjacency_matrix, variables, clauses, tmp, false);
}
/**
 * Funkce odstraňBikliku. Odstraní z grafu bikliku a všechny její incidentní hrany.
 */
bool remove_biclique(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, biclique_t& biclique){
	vector<bool> biclique_incidence(1+variables+clauses, false);
	for (int v : biclique.vertex_partition){
		biclique_incidence[v] = true;
	}
	for (int c : biclique.clauses_partition){
		biclique_incidence[c] = true;
	}
	
	for (int v : biclique.vertex_partition){
		for (int c : adjacency_matrix[v]){
			c = c > 0 ? c : -c;
			remove_value_from_neighbours(adjacency_matrix, c, v);
			if (!biclique_incidence[c] && adjacency_matrix[c].size() == 0) return false; // we create empty clause
		}
		adjacency_matrix[v].clear();
	}
	for (int c : biclique.clauses_partition){
		biclique_incidence[c] = true;
		for (int v : adjacency_matrix[c]){
			v = v > 0 ? v : -v;
			remove_value_from_neighbours(adjacency_matrix, v, c);
		}
		adjacency_matrix[c].clear();
	}
	return true;
}
/**
 * Funkce zvolZárodek. Vybere zárodek na základě strategie s z {-1 max, 0 rand, 1 min}
 */
void select_biclique(vector<biclique_t>& bicliques_seeds, biclique_t& biclique, const int& strategy){
	std::random_device rd;
    std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, bicliques_seeds.size()-1);
	
	
	switch(strategy){
		case -1:
			biclique = *std::max_element(bicliques_seeds.begin(), bicliques_seeds.end(), [](const biclique_t& b1, const biclique_t& b2) {return b1.clauses_count() < b2.clauses_count();});
			break;
		case 0:
			biclique = bicliques_seeds[dis(gen)];
		case 1:
			biclique = *std::min_element(bicliques_seeds.begin(), bicliques_seeds.end(), [](const biclique_t& b1, const biclique_t& b2) {return b1.clauses_count() < b2.clauses_count();});
			break;
	}
}


/**
 * Funkce omezZarodek. Omezí zárodek, tak aby byl omezenou biklikou.
 */
void reduce_biclique(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, biclique_t& biclique){
	sort(biclique.clauses_partition.begin(), biclique.clauses_partition.end(),[&](const int& c1,const int&c2){
		return adjacency_matrix[c1].size() < adjacency_matrix[c2].size();
	});
	biclique.clauses_partition.resize(pow(2,biclique.vertex_count())-1);
}


/**
 * Zkusí nalézt autarky na základě R_{G,B} množiny, tedy redukce R_{G,B} na základě biklikového pokrytí.
 */
bool find_autark(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, vector<biclique_t>& bicliques_cover, vector<bool>& uncovered_verticies, vector<int>& uncovered_clauses){
	bool changed = true;
	vector<int> vertex_biclique_incidence(adjacency_matrix.size()); // array for incidence 
	for (int i = 0; i < bicliques_cover.size(); ++i){
		for (int v : bicliques_cover[i].vertex_partition){
			vertex_biclique_incidence[v] = i;
		}
		for (int c : bicliques_cover[i].clauses_partition){
			vertex_biclique_incidence[c] = i;
		}
	}
	while (uncovered_clauses.size() > 0){
		int process_clause = uncovered_clauses.back();
		uncovered_clauses.pop_back();
		for (int v : adjacency_matrix[process_clause]){
			v = v > 0 ? v : -v;
			if (!uncovered_verticies[v]){
				uncovered_verticies[v] = true;
				for (int c : bicliques_cover[vertex_biclique_incidence[v]].clauses_partition){
					if (!uncovered_verticies[c]){
						uncovered_verticies[c] = true;
						uncovered_clauses.push_back(c);
					}
				}
				bicliques_cover[vertex_biclique_incidence[v]].id = -2;
			}
		}
	}

	
	bicliques_cover.erase(remove_if(bicliques_cover.begin(), bicliques_cover.end(), [](biclique_t& b){
			return (b.id == -2);
	}),bicliques_cover.end());
	
	return true;
}

/**
 * Funkce, která se pokusí udělat R_{G,B} redukci.
 */
bool reduce_graph_max_biclique_cover(vector<vector<int> > adjacency_matrix, const int& variables, const int& clauses, const int& k, const int& strategy, vector<biclique_t>& biclique_autarky){
	vector<biclique_t> bicliques_seed;
	biclique_t biclique;
	vector<int> uncovered_clauses;
	vector<bool> uncovered_verticies(adjacency_matrix.size(), false);
	vector<vector<int> > adjacency_matrix_backup = adjacency_matrix;
	//cerr << "biclique init" << endl;
	if (!forced_matching(adjacency_matrix, variables, clauses, biclique_autarky, true)){
		biclique_autarky.clear();	
		return false;
	}
	make_biclique_seed(adjacency_matrix, variables, clauses,bicliques_seed);	
	while (bicliques_seed.size() > 0){
		select_biclique(bicliques_seed, biclique, strategy);
			
		while (biclique.vertex_count() < k && enhance_biclique(adjacency_matrix,variables,clauses,biclique));
		reduce_biclique(adjacency_matrix, variables, clauses, biclique);
		//cerr << "Capitains log. Size of biclique is " << biclique[0] << ":" << biclique.size()-2-biclique[0] << endl;
		if (!remove_biclique(adjacency_matrix,variables, clauses, biclique) || !forced_matching(adjacency_matrix, variables, clauses, biclique_autarky, true)){
			while(biclique_autarky.end()->clauses_count() == 1){
				biclique_autarky.pop_back();
			}
			if (biclique_autarky.size() > 0) {
				biclique_autarky.pop_back();
			}
			for ( int v = 1; v <= variables; ++v ){
				if (adjacency_matrix[v].size() > 0){
					uncovered_verticies[v] = true;
				}
			}
			for ( int c = variables+1; c < adjacency_matrix.size(); ++c ){
				if (adjacency_matrix[c].size() > 0){
					uncovered_clauses.push_back(c);
					uncovered_verticies[c] = true;
				}
			}
			find_autark(adjacency_matrix_backup, variables, clauses, biclique_autarky, uncovered_verticies, uncovered_clauses);
			return false;
		}
		biclique_autarky.push_back(biclique);
		clear_biclique_seeds(adjacency_matrix, variables, clauses, bicliques_seed);		
	}
	
	find_max_matching(adjacency_matrix,variables,clauses, biclique_autarky, true);
	if (count_if(adjacency_matrix.begin()+variables+1, adjacency_matrix.end(), [](const vector<int>& n){return n.size() > 0; }) == 0){
		return true;
	}
	for ( int v = 1; v <= variables; ++v ){
		if (adjacency_matrix[v].size() > 0){
			uncovered_verticies[v] = true;
		}
	}
	for ( int c = variables+1; c < adjacency_matrix.size(); ++c ){
		if (adjacency_matrix[c].size() > 0){
			uncovered_clauses.push_back(c);
			uncovered_verticies[c] = true;
		}
	}
	find_autark(adjacency_matrix_backup, variables, clauses, biclique_autarky, uncovered_verticies, uncovered_clauses);
	return true;
}

/** 
 * Heuristika na biklikové pokrytí bipartitního grafu se strategií strategy.
 */
bool cover_bipartite_graph(vector<vector<int> > adjacency_matrix, const int& variables, const int& clauses, const int& k, const int& strategy){
	vector<biclique_t> bicliques_seed;
	biclique_t biclique;
	int nonempty_v = 0, nonempty_c = 0;
	
	if (!forced_matching(adjacency_matrix, variables, clauses)) return false;
	make_biclique_seed(adjacency_matrix, variables, clauses,bicliques_seed);	
	while (bicliques_seed.size() > 0){
		select_biclique(bicliques_seed, biclique, strategy);
			
		while (biclique.vertex_count() < k && enhance_biclique(adjacency_matrix,variables,clauses,biclique));
		reduce_biclique(adjacency_matrix, variables, clauses, biclique);
		
		if (!remove_biclique(adjacency_matrix,variables, clauses, biclique) || !forced_matching(adjacency_matrix, variables, clauses)){
			return false;
		}
		clear_biclique_seeds(adjacency_matrix, variables, clauses, bicliques_seed);
		if (find_max_matching(adjacency_matrix, variables, clauses) != 0){
			if (count_if(adjacency_matrix.begin()+variables+1, adjacency_matrix.end(), [](const vector<int>& n){return n.size() > 0; }) == 0)
				return true;
		}
		
	}
	
	find_max_matching(adjacency_matrix,variables,clauses);
	return count_if(adjacency_matrix.begin()+variables+1, adjacency_matrix.end(), [](const vector<int>& n){return n.size() > 0; }) == 0;
}

/**
 * Funkce generujZárodky pro backtrack biklikového pokrytí s K_{2,3}
 */
void backtrack_make_biclique_seed(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, vector<biclique_t>& bicliques_seed){
	vector<biclique_t> bad_bicliques_seed;
	make_biclique_seed(adjacency_matrix,variables,clauses,bad_bicliques_seed);
	for(biclique_t biclique : bad_bicliques_seed){
		for (int i = 0 ; i < biclique.clauses_count(); ++i){
			for (int j = i+1; j < biclique.clauses_count(); ++j){
				for (int k = j+1; k < biclique.clauses_count(); ++k){
					bicliques_seed.push_back(biclique_t());
					bicliques_seed.back().vertex_partition = biclique.vertex_partition;
					bicliques_seed.back().clauses_partition.push_back(biclique.clauses_partition[i]);
					bicliques_seed.back().clauses_partition.push_back(biclique.clauses_partition[j]);
					bicliques_seed.back().clauses_partition.push_back(biclique.clauses_partition[k]);
					bicliques_seed.back().id= bicliques_seed.size();
				}
			}
		}
	}
}

/**
 * Funkce vyčistiZárodky pro backtrack.
 */
void backtrack_clear_biclique_seeds(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, vector<biclique_t>& bicliques_seed){
	bicliques_seed.erase(remove_if(bicliques_seed.begin(),bicliques_seed.end(), [&](biclique_t& biclique){
			for( int c : biclique.clauses_partition){
				int vertex_hit = count_if(adjacency_matrix[c].begin(), adjacency_matrix[c].end(), [&](int v){
					v = v > 0 ? v : -v;
					return find(biclique.vertex_partition.begin(), biclique.vertex_partition.end(), v) != biclique.vertex_partition.end();
				});
				assert(vertex_hit <= biclique.vertex_count());
				if (vertex_hit < biclique.vertex_count()) return true;
			}
			return false;
	}),bicliques_seed.end());
}

/**
 * Backtrack přístup k pokrytí bipartitního grafu.
 */
bool backtrack_cover_bipartite_graph(vector<vector<int> > adjacency_matrix, const int& variables, const int& clauses, vector<biclique_t> bicliques_seed, biclique_t biclique){
	if (biclique.vertex_count()>0){
		if (!remove_biclique(adjacency_matrix,variables, clauses, biclique)) return false;
		if (!forced_matching(adjacency_matrix, variables, clauses)) return false;
		backtrack_clear_biclique_seeds(adjacency_matrix, variables, clauses, bicliques_seed);
	}
	int clauses_size = count_if(adjacency_matrix.begin()+variables+1, adjacency_matrix.end(), [](vector<int>& neighbours){ return neighbours.size() > 0;});
	int verticies_size = count_if(adjacency_matrix.begin(), adjacency_matrix.begin()+variables+1, [](vector<int>& neighbours){ return neighbours.size() > 0;});
	if (verticies_size > clauses_size){
		find_max_matching(adjacency_matrix, variables, clauses);
	}
	if (count_if(adjacency_matrix.begin()+1+variables, adjacency_matrix.end(), [](vector<int>& neig){ return neig.size() != 0;}) == 0){// if count of nonempty clauses is zero then sat
	 return true;	
	}
	for (biclique_t biclique_it : bicliques_seed){
		if (biclique_it.id > biclique.id)
			if (backtrack_cover_bipartite_graph(adjacency_matrix, variables, clauses, bicliques_seed, biclique_it)) return true;
	}
	return false;
}

/**
 * Backtrack algoritmus pro biklikové pokrytí K_{2,3}
 */
int backtrack(vector<vector<int> > adjacency_matrix, const int& variables, const int& clauses){
	vector<biclique_t> bicliques_seed;
	backtrack_make_biclique_seed(adjacency_matrix,variables,clauses,bicliques_seed);
	return backtrack_cover_bipartite_graph(adjacency_matrix, variables, clauses, bicliques_seed, biclique_t());
	
}

/**
 * Pro formuli provede jednotkovou propagaci a redukci čistých literálů
 **/
bool preprocess_formula(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses){
	bool preproces_value =  unit_propagation(adjacency_matrix, variables, clauses) && pure_literals(adjacency_matrix, variables, clauses);
	
	return preproces_value;
}

/** 
 * Vygeneruje incidenční graf náhodné formule, s pravděpodobností hrany p_positive. Pokud chceme signed graf pak pravděpodobnost negativní hrany je p_negative.
 */
void generate_random_bipartite_graph(vector<vector<int>>& adjacency_matrix,const int& variables,const int& clauses, const double& p_positive_edge, const double& p_negative_edge){
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
	adjacency_matrix.clear();
	adjacency_matrix.resize(1+clauses+variables);
	double p;
	vector<vector<int>> incidence_matrix(variables,vector<int>(clauses,0));
	for (int i = 0; i < variables; ++i){
		for (int j = 0; j < clauses; ++j){
			p = dis(gen);
			if (p < p_positive_edge){
				incidence_matrix[i][j] = 1;
			} else if (p < p_positive_edge+p_negative_edge){
				incidence_matrix[i][j] = -1;
			}
		}
	}
	for (int i = 0; i < variables; ++i){
		for (int j = 0; j < clauses; ++j){
			if (incidence_matrix[i][j] != 0){
				adjacency_matrix[i+1].push_back((variables+j+1)*incidence_matrix[i][j]);
				adjacency_matrix[variables+j+1].push_back((i+1)*incidence_matrix[i][j]);
			}
		}
	}
}
/**
  * Vygeneruje náhodný k-KNF formuli.
  */
void generate_random_k_sat(vector<vector<int>>& adjacency_matrix,const int& variables,const int& clauses, const int& k){
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, variables);
	adjacency_matrix.clear();
	adjacency_matrix.resize(1+clauses+variables);
	vector<bool> incidence_matrix(variables+1,0);
	int variable = 0, variables_count;
	for (int i = 1; i <= clauses; ++i){
		for (int j = 1; j < incidence_matrix.size(); ++j){incidence_matrix[j] = false;}
		variables_count = 0;
		while (variables_count < k){
			variable = dis(gen);
			while (incidence_matrix[variable]){
				variable = dis(gen);
			}
			incidence_matrix[variable] = true;
			adjacency_matrix[i+variables].push_back(variable);
			adjacency_matrix[variable].push_back(i+variables);
			++variables_count;
		}
	}

}
/**
 *  Vyzkouší pokrýt incidenční graf adjacency_matrix s omezením na bikliky z biklikového pokrytí biclique_bouna a strategií strategy.
 */
bool test_bipartite_h(vector<vector<int>> adjacency_matrix, const int& variables, const int& clauses, const int& biclique_bound,const int& strategy){
	
	return cover_bipartite_graph(adjacency_matrix, variables, clauses, biclique_bound, strategy);
}


/**
 * Experiment, který zjistí jestli exituje V_C perfektní párování.
 */
bool test_matching(vector<vector<int>> adjacency_matrix, int& variables, int& clauses){
	int matching_size = find_max_matching(adjacency_matrix, variables, clauses);
	int clauses_size = count_if(adjacency_matrix.begin()+variables+1, adjacency_matrix.end(), [](vector<int>& neighbours){ return neighbours.size() > 0;});
	int verticies_size = count_if(adjacency_matrix.begin(), adjacency_matrix.begin()+variables+1, [](vector<int>& neighbours){ return neighbours.size() > 0;});
	return (clauses_size == 0);
}

/**
  * Experiment, který zkouší vystavět autarky z R_{G,B} množiny.
  */
bool test_autartky_h(vector<vector<int>> adjacency_matrix, const int& variables, const int& clauses, const int& biclique_bound,const int& strategy, int& clauses_covered, int& variables_used){
	vector<biclique_t> biclique_autarky;
	vector<bool>counted(adjacency_matrix.size(), false);
	clauses_covered = 0;
	variables_used = 0;
	bool return_val = reduce_graph_max_biclique_cover(adjacency_matrix, variables, clauses, biclique_bound, strategy, biclique_autarky);

		
	for_each(biclique_autarky.begin(), biclique_autarky.end(),[&](const biclique_t& b){
		for (int v : b.vertex_partition){
			if (!counted[v]){
				counted[v] = true;
				++variables_used;
			}
		}
		for (int c : b.clauses_partition){
			if (!counted[c] && c > 0){
				counted[c] = true;
				++clauses_covered;
			}
		}
	});
	return return_val;
}
