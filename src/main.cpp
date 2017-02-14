#include <iostream>
#include <fstream>
#include <cstring>

#include "nex_parser.h"

#include "freqdiff.h"
#include "minrs.h"

using namespace std;

int main(int argc, char** argv) {

	if (argc != 2) {
		cout << "Usage: exec filename" << endl;
		cout << "More info to come" << endl;
		return 0;
	}

	ifstream fin(argv[1]);
	vector<Tree*> trees = parse_nex(fin);

	//Tree* fdt = freqdiff(trees, Tree::get_taxas_num() > 1000);
	Tree* minrs = minRS(trees);
	//cout << fdt->to_newick() << endl;
	cout << minrs << endl;

	//delete fdt;
}
