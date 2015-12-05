#ifndef TREE_H_
#define TREE_H_

#include <vector>
#include <unordered_map>

class Tree {
public:
	static std::unordered_map<int,std::string> taxa_names;

	Tree();
	Tree(std::string& newick_str);
	Tree(Tree* other);
	virtual ~Tree();

	static size_t get_taxas_num();

	class Node;

	Node* get_node(int i);
	Node* get_root();
	size_t get_nodes_num();

	Node* get_leaf(int taxa);

	Tree::Node* add_node(int taxa = -1);
	void delete_nodes(bool* to_delete);

	// Applies lazy children deletions and resorts the nodes in topological order
	void fix_tree(Node* root = NULL);

	std::string to_string();
	std::string to_newick(Node* node = NULL);

	void reorder(); // puts for each node, puts heaviest subtree as first child

private:
	std::vector<Node*> nodes;
	std::vector<Node*> taxa_to_leaf; // TODO: how much the map reduces performance?
	std::unordered_map<int, Node*> taxa_to_leaf_map;

	static std::unordered_map<std::string,int> taxa_ids;
	static int get_taxa_id(std::string& taxa);

	Tree::Node* build_tree(const char*&);
	void fix_tree_supp(Tree::Node* root);
};

class Tree::Node {
public:
	static const int NONE = -1;

	std::vector<Node*> children;

	Node* parent;
	size_t pos_in_parent;
	int id, secondary_id;
	int taxa, weight;
	size_t size; // number of leaves in subtree

	Node(int id);// : parent(NULL), pos_in_parent(NONE), id(id), taxa(NONE), weight(0) {};
	Node(int id, int taxa);// : parent(NULL), pos_in_parent(NONE), id(id), taxa(taxa), weight(0) {};

	size_t get_children_num();

	void add_child(Node* child);
	void set_child(Node* child, size_t pos);
	void null_child(size_t pos);
	void clear_children();
	void fix_children(); // remove null children and recompacts others

	bool is_leaf();
	bool is_root();

	std::string to_string();
	std::string to_newick();
};

#endif
