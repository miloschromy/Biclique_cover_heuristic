#include <iostream>
#include <vector>
#include <list>
using namespace std;
/**
 * Struktura reprezentující bikliku
 */
struct biclique_t;

/**
 * Funkce pro odstranění hodnty z matice sousednosti. Pokud je seznam setřízený naruší jeho strukturu.
 */
int remove_value_from_neighbours(vector<vector<int> >& adjacency_matrix, const int& vertex, const int& value);

/**
 *Naleze zlepšující cestu v grafu.
 *Matching reprezentuje párování.
 */
bool find_augumenting_path(vector<vector<int>>& adjacency_matrix, const int& variables, const int& clauses, vector<int>& matching, list<int>& unpaired);

/** 
 * Párování vytiskne na výstup.
 */
void print_matching(vector<int>& matching, int& variables);

/**
 * Smaže množinu R_{G,M} z grafu. Redukce R_{G,M}
 */
bool remove_RgmSet(vector<vector<int>>& adjacency_matrix,const int& variables, const int& clauses, vector<int>& matching, vector<biclique_t>& removed_matching, const bool& store_removed);

/**
 * Funkce, která najde maximální párování. V heuristice je to funkce testMatched(G). Maximální párování si pamatuje.
 */
size_t find_max_matching(vector<vector<int>>& adjacency_matrix,const int& variables, const int& clauses, vector<biclique_t>& removed_matching,const bool& store_removed);

/**
 * Funkce, která najde maximální párování. V heuristice je to funkce testMatched(G). Bez ukládání maximálního párování.
 */
size_t find_max_matching(vector<vector<int>>& adjacency_matrix,const int& variables, const int& clauses);

/**
 * Funkce pro načtení knf ze souboru, jež odpovídá specifikaci knf souboru v sat comptetiton (http://baldur.iti.kit.edu/sat2014/competitions.html). 
 */
void load_formula(vector<vector<int> >& adjacency_matrix, int& variables, int& clauses, istream& input);

/**
 * Vytiskne formuli, jejíž incidenční graf je adjacency_matrix na výstup output.
 */
void print_formula(const vector<vector<int> >& adjacency_matrix,const int& variables, const int&  clauses, ostream& output);

/**
 * Zapíše formuli do souboru s názvem outfile_name
 */
void print_formula_file(vector<vector<int> >& adjacency_matrix,const int& variables, const int&  clauses, const char* outfile_name);

/**
 * Jendotková propagace na klauzulích.
 */
bool unit_propagation(vector<vector<int> >& adjacency_matrix,const int& variables, const int&  clauses);

/**
 * Kontrola, zdali je literál check_variable čistý.
 */
bool is_pure(vector<vector<int> >& adjacency_matrix, const int& check_variable);

/**
 * Aplikace redukce čistých literálů.
 */
bool pure_literals(vector<vector<int> >& adjacency_matrix,const int& variables, const int&  clauses);


/**
 * Funkce z heuristiky rozšiřBikliku. Pokud by rozšíření o největší vrchol vytvořilo nesplnitelnou formuli vrací false a nerozšíří ji.
 */
bool enhance_biclique(const vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, biclique_t& biclique);


/**
 * Funkce vygenerujZárodky. Pro každou dvojici vygeneruje zárodek, tedy bikliku K_2,?
 */
bool make_biclique_seed(const vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, vector<biclique_t>& bicliques_seed);

/**
 * Funkce vyčistiZárodky
 */
bool clear_biclique_seeds(const vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, vector<biclique_t>& bicliques_seed);

/*
 * Funkce jednotkováGPropagace. Pročistí graf od vrcholů z druhé partity stupně 1 tak, že hranu jim incidentní přidá do pokrytí.
 */
bool forced_matching(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, vector<biclique_t>& used_matching, bool store_matching);

/**
 * Funkce jednotkovaGPropagace, která si nepamatuje pokrytí, takže se jen ptáme na existenci.
 */
bool forced_matching(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses);

/**
 * Funkce odstraňBikliku. Odstraní z grafu bikliku a všechny její incidentní hrany.
 */
bool remove_biclique(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, biclique_t& biclique);

/**
 * Funkce zvolZárodek. Vybere zárodek na základě strategie s z {-1 max, 0 rand, 1 min}
 */
void select_biclique(vector<biclique_t>& bicliques_seeds, biclique_t& biclique, const int& strategy);


/**
 * Funkce omezZarodek. Omezí zárodek, tak aby byl omezenou biklikou.
 */
void reduce_biclique(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, biclique_t& biclique);


/**
 * Zkusí nalézt autarky na základě R_{G,B} množiny, tedy redukce R_{G,B} na základě biklikového pokrytí.
 */
bool find_autark(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, vector<biclique_t>& bicliques_cover, vector<bool>& uncovered_verticies, vector<int>& uncovered_clauses);

/**
 * Funkce, která se pokusí udělat R_{G,B} redukci.
 */
bool reduce_graph_max_biclique_cover(vector<vector<int> > adjacency_matrix, const int& variables, const int& clauses, const int& k, const int& strategy, vector<biclique_t>& biclique_autarky);

/** 
 * Heuristika na biklikové pokrytí bipartitního grafu se strategií strategy.
 */
bool cover_bipartite_graph(vector<vector<int> > adjacency_matrix, const int& variables, const int& clauses, const int& k, const int& strategy);

/**
 * Funkce generujZárodky pro backtrack biklikového pokrytí s K_{2,3}
 */
void backtrack_make_biclique_seed(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, vector<biclique_t>& bicliques_seed);

/**
 * Funkce vyčistiZárodky pro backtrack.
 */
void backtrack_clear_biclique_seeds(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses, vector<biclique_t>& bicliques_seed);

/**
 * Backtrack přístup k pokrytí bipartitního grafu.
 */
bool backtrack_cover_bipartite_graph(vector<vector<int> > adjacency_matrix, const int& variables, const int& clauses, vector<biclique_t> bicliques_seed, biclique_t biclique);

/**
 * Backtrack algoritmus pro biklikové pokrytí K_{2,3}
 */
int backtrack(vector<vector<int> > adjacency_matrix, const int& variables, const int& clauses);

/**
 * Pro formuli provede jednotkovou propagaci a redukci čistých literálů
 **/
bool preprocess_formula(vector<vector<int> >& adjacency_matrix, const int& variables, const int& clauses);

/** 
 * Vygeneruje incidenční graf náhodné formule, s pravděpodobností hrany p_positive. Pokud chceme signed graf pak pravděpodobnost negativní hrany je p_negative.
 */
void generate_random_bipartite_graph(vector<vector<int>>& adjacency_matrix,const int& variables,const int& clauses, const double& p_positive_edge, const double& p_negative_edge);

/**
  * Vygeneruje náhodný k-KNF formuli.
  */
void generate_random_k_sat(vector<vector<int>>& adjacency_matrix,const int& variables,const int& clauses, const int& k);

/**
 *  Vyzkouší pokrýt incidenční graf adjacency_matrix s omezením na bikliky z biklikového pokrytí biclique_bouna a strategií strategy.
 */
bool test_bipartite_h(vector<vector<int>> adjacency_matrix, const int& variables, const int& clauses, const int& biclique_bound,const int& strategy);

/**
 * Experiment, který zjistí jestli exituje V_C perfektní párování.
 */
bool test_matching(vector<vector<int>> adjacency_matrix, int& variables, int& clauses);

/**
  * Experiment, který zkouší vystavět autarky z R_{G,B} množiny.
  */
bool test_autartky_h(vector<vector<int>> adjacency_matrix, const int& variables, const int& clauses, const int& biclique_bound,const int& strategy, int& clauses_covered, int& variables_used);
