#ifndef NEX_PARSER_H_
#define NEX_PARSER_H_

#include <istream>
#include <iostream>
#include <cassert>

#include "Tree.h"

void parse_nex(std::istream& in) {
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
	while (getline(in, s)) {
		if (s.find("END;") != std::string::npos)
			break;
		Tree* t = new Tree(s);
		std::cout << t->to_string() << std::endl;
	}
}


#endif
