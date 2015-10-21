#include <vector>
#include <unordered_map>
#include <cassert>
#include <sstream>
#include <iostream>

#include "Tree.h"

std::unordered_map<std::string,int> Tree::taxas;

Tree::Tree(std::string& newick_str) {
	const char* str = newick_str.c_str();
	while (*str != '(') str++;
	this->build_tree(str);
}

Tree::~Tree() {}

std::string Tree::to_string() {
	std::stringstream ss;
	for (Node* node : nodes) {
		ss << node->to_string() << std::endl;
	}
	return ss.str();
}

Tree::Node* Tree::build_tree(const char*& str) {
	assert(*str == '(');
	str++;

	int vecpos = nodes.size();
	nodes.push_back(new Tree::Node(vecpos));

	assert(*str == '(' || isdigit(*str));
	while (true) {
		Tree::Node* subtree;
		if (*str == '(') {
			subtree = build_tree(str);
		} else {
			std::string taxa;
			while (*str != ':') {
				taxa += *str;
				str++;
			}
			subtree = new Node(nodes.size(), get_taxa_id(taxa));
			nodes.push_back(subtree);
		}
		nodes[vecpos]->add_child(subtree);
		while (*str != ',' && *str != ')') str++;
		if (*str == ',') {
			str++;
		} else if (*str == ')') {
			str++;
			break;
		}
	}

	return nodes[vecpos];
}

int Tree::get_taxa_id(std::string& taxa) {
	auto taxa_entry = Tree::taxas.find(taxa);
	if (taxa_entry == Tree::taxas.end()) {
		int taxa_id = taxas.size();
		taxas.insert({taxa,taxa_id});
		return taxa_id;
	} else {
		return taxa_entry->second;
	}
}

void Tree::Node::add_child(Node* child) {
	children.push_back(child);
}

std::string Tree::Node::to_string() {
	std::stringstream ss;
	ss << id << " ";
	if (taxa == NONE) {
		ss << "[";
		for (size_t i = 0; i < children.size(); i++) {
			ss << (i > 0 ? "," : "") << children[i]->id;
		}
		ss << "] ";
		ss << children.size();
	} else {
		ss << taxa;
	}
	return ss.str();
}
