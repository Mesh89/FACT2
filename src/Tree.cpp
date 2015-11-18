#include <vector>
#include <unordered_map>
#include <cassert>
#include <sstream>
#include <iostream>
#include <queue>

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

Tree* Tree::copy() {
	std::string t_newick = to_newick();
	Tree* other = new Tree(t_newick);
	for (Node* node : nodes) {
		other->get_node(node->get_id())->set_weight(node->get_weight());
	}
	return other;
}

std::string Tree::to_string() {
	std::stringstream ss;
	for (Node* node : nodes) {
		ss << node->to_string() << std::endl;
	}
	return ss.str();
}

std::string Tree::to_newick(Node* node) {
	if (node == NULL) {
		node = get_root();
	}

	if (node->is_leaf()) {
		return taxa_names[node->get_taxa()];
	}
	std::stringstream ss;
	ss << "(";
	for (size_t i = 0; i < node->get_children_num(); i++) {
		if (i > 0) ss << ",";
		ss << to_newick(node->get_child(i));
	}
	ss << ")";
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

Tree::Node* Tree::add_node() {
	Tree::Node* newnode = new Tree::Node(get_nodes_num());
	nodes.push_back(newnode);
	return newnode;
}

void Tree::delete_nodes(bool* to_delete) {
	for (size_t i = 1; i < get_nodes_num(); i++) {
		if (to_delete[i]) {
			for (Tree::Node* child : nodes[i]->get_children()) {
				nodes[i]->get_parent()->add_child(child);
			}
			nodes[i]->get_parent()->null_child(nodes[i]->get_pos_in_parent());
			delete nodes[i];
			nodes[i] = NULL;
		}
	}
	fix_tree();
}

void Tree::fix_tree() {
	Tree::Node* root = get_root();
	nodes.clear();
	fix_tree_supp(root);
}

void Tree::fix_tree_supp(Tree::Node* curr) {
	curr->set_id(nodes.size());
	nodes.push_back(curr);
	curr->fix_children();
	for (Tree::Node* child : curr->get_children()) {
		fix_tree_supp(child);
	}
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
			while (*str != ':' && *str != ',' && *str != ')') {
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

void Tree::Node::fix_children() {
	size_t curr_pos = 0;
	for (size_t i = 0; i < children.size(); i++) {
		if (children[i] != NULL) {
			children[curr_pos] = children[i];
			children[curr_pos]->pos_in_parent = curr_pos;
			curr_pos++;
		}
	}
	children.resize(curr_pos);
}



Tree::Node* Tree::Node::get_parent() {
	return parent;
}

size_t Tree::Node::get_pos_in_parent() {
	return pos_in_parent;
}

void Tree::Node::add_child(Node* child) {
	child->parent = this;
	child->pos_in_parent = children.size();
	children.push_back(child);
}
void Tree::Node::set_child(Node* child, size_t pos) {
	child->parent = this;
	child->pos_in_parent = pos;
	children[pos] = child;
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

int Tree::Node::set_id(int id) {
	return this->id = id;
}

Tree::Node* Tree::Node::get_child(int i) {
	return children[i];
}

std::vector<Tree::Node*> Tree::Node::get_children() {
	return children;
}

size_t Tree::Node::get_children_num() {
	return children.size();
}

bool Tree::Node::is_root() {
	return parent == NULL;
}

void Tree::Node::null_child(size_t pos) {
	children[pos] = NULL;
}

void Tree::Node::clear_children() {
	children.clear();
}

int Tree::Node::get_weight() {
	return weight;
}
void Tree::Node::set_weight(int weight) {
	this->weight = weight;
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
	if (!is_root()) {
		ss << " (weight: " << get_weight() << ")";
	}
	return ss.str();
}

