#include <iostream>
#include <fstream>

#include "nex_parser.h"

#include "freqdiff.h"

using namespace std;

int main(int argc, char** argv) {

	if (argc == 1 || argc > 3) {
		cout << "Usage: exec filename [algorithms]" << endl;
		cout << "More info to come" << endl;
	}

	ifstream fin(argv[1]);
	vector<Tree*> trees = parse_nex(fin);
	std::cout << "TREES PARSED" << std::endl;
	freqdiff(trees);
}
