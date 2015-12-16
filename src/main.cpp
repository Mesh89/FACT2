#include <iostream>
#include <fstream>
#include <cstring>

#include "nex_parser.h"

#include "freqdiff.h"

using namespace std;

int main(int argc, char** argv) {

	if (argc != 3) {
		cout << "Usage: exec filename [use cp]" << endl;
		cout << "More info to come" << endl;
		return 0;
	}

	ifstream fin(argv[1]);
	vector<Tree*> trees = parse_nex(fin);
	Tree* fdt = freqdiff(trees, strcmp(argv[2], "1") == 0);
	cout << fdt->to_newick() << endl;
}
