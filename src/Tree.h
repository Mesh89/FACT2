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

	Tree::Node* add_node();
	void delete_nodes(bool* to_delete);

	// Applies lazy children deletions and resorts the nodes in topological order
	void fix_tree();

	std::string to_string();

private:
	std::vector<Node*> nodes;
	std::vector<Node*> taxa_to_leaf;

	static std::unordered_map<std::string,int> taxa_ids;
	static int get_taxa_id(std::string& taxa);

	Tree::Node* build_tree(const char*&);
	void fix_tree_supp(Tree::Node* root);
};

class Tree::Node {
public:
	static const int NONE = -1;

	Node(int id) : parent(NULL), pos_in_parent(NONE), id(id), taxa(NONE), weight(0) {};
	Node(int id, int taxa) : parent(NULL), pos_in_parent(NONE), id(id), taxa(taxa), weight(0) {};

	int get_taxa();
	int get_id();
	int set_id(int id);

	Node* get_parent();
	size_t get_pos_in_parent();

	Node* get_child(int i);
	std::vector<Node*> get_children();
	size_t get_children_num();

	void add_child(Node* child);
	void set_child(Node* child, size_t pos);
	void null_child(size_t pos);
	void remove_children(bool* to_delete);
	void clear_children();
	void fix_children(); // remove null children and recompacts others

	bool is_leaf();
	bool is_root();

	int get_weight();
	void set_weight(int weight);

	std::string to_string();


private:
	Node* parent;
	size_t pos_in_parent;
	int id, taxa, weight;
	std::vector<Node*> children;
};

#endif
