#ifndef TREE_H_
#define TREE_H_

#include <vector>
#include <unordered_map>

class Tree {
public:
	static std::unordered_map<int,std::string> taxa_names;

	Tree(std::string& newick_str);
	virtual ~Tree();

	static size_t get_taxas_num();

	class Node;

	Node* get_node(int i);
	Node* get_root();
	size_t get_nodes_num();

	Node* get_leaf(int taxa);

	void delete_nodes(bool* to_delete);

	std::string to_string();

private:
	std::vector<Node*> nodes;
	std::vector<Node*> taxa_to_leaf;

	static std::unordered_map<std::string,int> taxa_ids;
	static int get_taxa_id(std::string& taxa);

	Tree::Node* build_tree(const char*&);
};

class Tree::Node {
public:
	static const int NONE = -1;

	Node(int id) : parent(NULL), id(id), taxa(NONE) {};
	Node(int id, int taxa) : parent(NULL), id(id), taxa(taxa) {};

	int get_taxa();
	int get_id();

	Node* get_parent();
	void add_child(Node* child);
	Node* get_child(int i);
	std::vector<Node*> get_children();

	void remove_children(bool* to_delete);
	void clear_children();

	bool is_leaf();
	bool is_root();

	std::string to_string();

private:
	Node* parent;
	int id, taxa;
	std::vector<Node*> children;
};

#endif
