#include <vector>
#include <unordered_map>
#include <cassert>
#include <sstream>
#include <iostream>

#include "Tree.h"

std::unordered_map<std::string,int> Tree::taxa_ids;
std::unordered_map<int,std::string> Tree::taxa_names;

Tree::Tree(std::string& newick_str) {
	const char* str = newick_str.c_str();
	while (*str != '(') str++;
	this->build_tree(str);

	taxa_to_leaf.resize(Tree::get_taxas_num());
	for (Node* node : nodes) {
		if (node->is_leaf())
			taxa_to_leaf[node->get_taxa()] = node;
	}
}

Tree::~Tree() {}

std::string Tree::to_string() {
	std::stringstream ss;
	for (Node* node : nodes) {
		ss << node->to_string() << std::endl;
	}
	return ss.str();
}

size_t Tree::get_taxas_num() {
	return taxa_ids.size();
}

size_t Tree::get_nodes_num() {
	return nodes.size();
}

Tree::Node* Tree::get_node(int i) {
	return nodes[i];
}

Tree::Node* Tree::get_root() {
	return get_node(0);
}

Tree::Node* Tree::get_leaf(int taxa) {
	return taxa_to_leaf[taxa];
}

Tree::Node* Tree::Node::get_parent() {
	return parent;
}

void Tree::delete_nodes(bool* to_delete) {
	for (size_t i = 0; i < get_nodes_num(); i++) {
		nodes[i]->remove_children(to_delete);
	}

	int curr_pos = 0;
	for (size_t i = 0; i < get_nodes_num(); i++) {
		if (to_delete[i]) {
			delete nodes[i];
		} else {
			nodes[curr_pos++] = nodes[i];
		}
	}
	nodes.resize(curr_pos);
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
	auto taxa_entry = Tree::taxa_ids.find(taxa);
	if (taxa_entry == Tree::taxa_ids.end()) {
		int taxa_id = taxa_ids.size();
		taxa_ids.insert({taxa,taxa_id});
		taxa_names.insert({taxa_id,taxa});
		return taxa_id;
	} else {
		return taxa_entry->second;
	}
}

void Tree::Node::add_child(Node* child) {
	child->parent = this;
	children.push_back(child);
}

bool Tree::Node::is_leaf() {
	return taxa != NONE;
}

int Tree::Node::get_taxa() {
	return taxa;
}

int Tree::Node::get_id() {
	return id;
}

Tree::Node* Tree::Node::get_child(int i) {
	return children[i];
}

std::vector<Tree::Node*> Tree::Node::get_children() {
	return children;
}

bool Tree::Node::is_root() {
	return parent == NULL;
}

void Tree::Node::remove_children(bool* to_delete) {
	int curr_pos = 0;
	for (size_t i = 0; i < children.size(); i++) {
		if (!to_delete[children[i]->get_id()]) {
			children[curr_pos++] = children[i];
		}
	}
	children.resize(curr_pos);
}

void Tree::Node::clear_children() {
	children.clear();
}

std::string Tree::Node::to_string() {
	std::stringstream ss;
	ss << id << " ";
	if (taxa == NONE) {
		ss << "[";
		for (size_t i = 0; i < children.size(); i++) {
			ss << (i > 0 ? "," : "") << children[i]->id;
		}
		ss << "]";
	} else {
		ss << taxa_names[taxa];
	}
	return ss.str();
}
