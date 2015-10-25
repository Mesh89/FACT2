#ifndef NEX_PARSER_H_
#define NEX_PARSER_H_

#include <istream>
#include <iostream>
#include <cassert>

#include "Tree.h"

std::vector<Tree*> parse_nex(std::istream& in) {
	std::string s;
	while (getline(in, s)) {
		if (s.find("BEGIN TREES;") != std::string::npos)
			break;
	}
	getline(in, s); // "translate"

	// Parse taxas
	while (getline(in, s)) {
		if (s.find(";") != std::string::npos)
			break;
	}
	// Parse trees
	std::vector<Tree*> trees;
	while (getline(in, s)) {
		if (s.find("END;") != std::string::npos)
			break;
		trees.push_back(new Tree(s));
	}
	return trees;
}


#endif
