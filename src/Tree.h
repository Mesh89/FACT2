#ifndef TREE_H_
#define TREE_H_

#include <vector>
#include <unordered_map>

class Tree {
public:
	Tree(std::string& newick_str);
	virtual ~Tree();

	std::string to_string();

	class Node {
	public:
		static const int NONE = -1;

		Node(int id) : id(id), taxa(NONE) {};
		Node(int id, int taxa) : id(id), taxa(taxa) {};

		void add_child(Node* child);
		std::string to_string();

	private:
		int id, taxa;
		std::vector<Node*> children;
	};

private:
	std::vector<Node*> nodes;
	static std::unordered_map<std::string,int> taxas;
	static int get_taxa_id(std::string& taxa);

	Tree::Node* build_tree(const char*&);
};

#endif
