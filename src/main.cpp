#include <iostream>
#include <fstream>
#include <cstring>

#include "nex_parser.h"

#include "freqdiff.h"
#include "local_consensus.h"

using namespace std;

int main(int argc, char** argv) {

	if (argc != 3) {
		cout << "Usage: exec [freq|minrs|minis] filename" << endl;
		cout << "More info to come" << endl;
		return 0;
	}

	string algo(argv[1]);

	ifstream fin(argv[2]);
	vector<Tree*> trees = parse_nex(fin);

	Tree* consensus;
	if (algo == "freq") {
		consensus = freqdiff(trees, Tree::get_taxas_num() > 1000);
	} else if (algo == "minrlc_exact") {
		consensus = minRLC_exact(trees);
	} else if (algo == "minilc_exact") {
		consensus = minILC_exact(trees);
	} else if (algo == "aho-build") {
		consensus = ahoBuild(trees);
	} else {
		std::cerr << "Algorithm not supported" << std::endl;
		return 1;
	}

	if (consensus == NULL) {
		std::cout << "No valid consensus found." << std::endl;
	} else {
		cout << consensus->to_newick() << endl;
	}

	for (Tree* tree : trees) {
		delete tree;
	}
	delete consensus;
}
