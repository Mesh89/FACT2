/*
 * leaves_ranges.h
 *
 *  Created on: 21 Oct 2015
 *      Author: Mesh
 */

#ifndef TAXAS_RANGES_H_
#define TAXAS_RANGES_H_

// taxas contains the leaves (taxas ids) in left-to-right order in the tree
// intervals[node id].start/end contain the starting and ending point of the interval of taxas descendants of the node
struct taxas_ranges_t {
	size_t taxas_num = 0;
	int* taxas;
	struct interval_t {
		int start, end;
	}* intervals;

	taxas_ranges_t() : taxas(NULL), intervals(NULL) {}
	taxas_ranges_t(size_t taxas_num, size_t nodes_num) : taxas(new int[taxas_num]), intervals(new interval_t[nodes_num]) {}
	~taxas_ranges_t() { delete[] taxas; delete[] intervals; }
};

void build_taxas_ranges_supp(Tree::Node* node, taxas_ranges_t* tr) {
	if (node->is_leaf()) {
		tr->taxas[tr->taxas_num] = node->taxa;
		tr->intervals[node->id].start = tr->intervals[node->id].end = tr->taxas_num;
		tr->taxas_num++;
		return;
	}
	for (Tree::Node* child : node->children) {
		build_taxas_ranges_supp(child, tr);
	}
	tr->intervals[node->id].start = tr->intervals[node->children[0]->id].start;
	tr->intervals[node->id].end = tr->intervals[node->children[node->get_children_num()-1]->id].end;
}

taxas_ranges_t* build_taxas_ranges(Tree* tree) {
	taxas_ranges_t* tr = new taxas_ranges_t(tree->get_leaves_num(), tree->get_nodes_num());
	build_taxas_ranges_supp(tree->get_root(), tr);
	return tr;
}

int* get_taxas_ranks(taxas_ranges_t* tr) {
	int* ranks = new int[Tree::get_taxas_num()];
	for (size_t i = 0; i < tr->taxas_num; i++) {
		ranks[tr->taxas[i]] = i;
	}
	return ranks;
}


#endif /* TAXAS_RANGES_H_ */
