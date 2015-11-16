/*
 * freqdiff.h
 *
 *  Created on: 21 Oct 2015
 *      Author: Mesh
 */

#ifndef FREQDIFF_H_
#define FREQDIFF_H_

#include <iostream>
#include <cassert>

#include "taxas_ranges.h"
#include "lca_preprocessing.h"
#include "utils.h"

struct node_bitvec_t {
	int tree_id, node_id;
	bool* bitvec;

	node_bitvec_t() : tree_id(-1), node_id(-1), bitvec(NULL) {}
	node_bitvec_t(int tree_id, int node_id, bool* bitvec) : tree_id(tree_id), node_id(node_id), bitvec(bitvec) {}
};

inline bool equal(bool* bitv1, bool* bitv2, size_t size) {
	for (size_t i = 0; i < size; i++) {
		if (bitv1[i] != bitv2[i]) return false;
	}
	return true;
}

// calculate cluster weights using kn^2 method
void calc_w_kn2(std::vector<Tree*>& trees) {
	// Generate bit vectors
	std::vector<node_bitvec_t> bitvectors;
	size_t n = Tree::get_taxas_num();
	for (size_t t = 0; t < trees.size(); t++) {
		Tree* tree = trees[t];
		taxas_ranges_t* tr = build_taxas_ranges(tree);

		size_t nodes_num = tree->get_nodes_num();
		for (size_t i = 1; i < nodes_num; i++) {
			Tree::Node* node = tree->get_node(i);
			if (!node->is_leaf()) {
				// fill bit vector for current node
				bool* bitvector = new bool[n];
				std::fill(bitvector, bitvector+n, 0);
				for (int j = tr->intervals[node->get_id()].start; j <= tr->intervals[node->get_id()].end; j++) {
					bitvector[tr->taxas[j]] = 1;
				}
				node_bitvec_t node_bitvec(t, node->get_id(), bitvector);
				bitvectors.push_back(node_bitvec);
			}
		}
	}

	// Sort bit vectors
	int bitvc = bitvectors.size();
	node_bitvec_t* bit0 = new node_bitvec_t[bitvc];
	std::copy(bitvectors.begin(), bitvectors.end(), bit0);
	node_bitvec_t* bit1 = new node_bitvec_t[bitvc];
	int bit0c = 0, bit1c = 0;

	for (int i = n-1; i >= 0; i--) {
		for (int j = 0; j < bitvc; j++) {
			if (bit0[j].bitvec[i] == 1) {
				bit1[bit1c++] =  bit0[j];
				// we signal that this position is free as the occupant had the curr bit set to 1
				bit0[j].bitvec = NULL;
			}
		}
		for (int j = 0; j < bitvc; j++) {
			if (bit0[j].bitvec != NULL) {
				bit0[bit0c++] = bit0[j];
			}
		}
		for (int j = 0; j < bit1c; j++) {
			if (bit1[j].bitvec != NULL) {
				bit0[bit0c++] = bit1[j];
			}
		}
		assert(bit0c == bitvc);
		bit0c = bit1c = 0;
	}
	// bit0 contains sorted bitvectors

	// count adjacent equal bitvectors and set weights
	for (int i = 0; i < bitvc; ) {
		int w = 1;
		while (i+w < bitvc && equal(bit0[i].bitvec,bit0[i+w].bitvec,n)) w++;
		for (int j = i; j < i+w; j++) {
			trees[bit0[j].tree_id]->get_node(bit0[j].node_id)->set_weight(w);
		}
		i += w;
	}
}


// for node_id in T1:
// start[node_id] = smallest rank (i.e. leftmost leaf) in T2 among \Lambda(T1[node_id])
// stop[node_id] = biggest rank (i.e. rightmost leaf) in T2 among \Lambda(T1[node_id])
void compute_start_stop(Tree* tree1, Tree* tree2, int* t2_leaves_ranks, int*& start, int*& stop) {
	start = new int[tree1->get_nodes_num()];
	std::fill(start, start+tree1->get_nodes_num(), INT32_MAX);
	stop = new int[tree1->get_nodes_num()];
	std::fill(stop, stop+tree1->get_nodes_num(), 0);

	// for each cluster in t1, find its leftmost and rightmost leaves in t2
	for (int i = tree1->get_nodes_num()-1; i >= 0; i--) {
		Tree::Node* node = tree1->get_node(i);
		if (node->is_leaf()) {
			start[i] = stop[i] = t2_leaves_ranks[node->get_taxa()];
		} else {
			for (Tree::Node* child : node->get_children()) {
				start[i] = std::min(start[i], start[child->get_id()]);
				stop[i] = std::max(stop[i], stop[child->get_id()]);
			}
		}
	}
}


bool* filter_clusters(Tree* tree1, Tree* tree2, taxas_ranges_t* t1_tr, taxas_ranges_t* t2_tr, lca_t* t2_lcas) {

	int* t2_leaves_ranks = new int[Tree::get_taxas_num()];
	for (size_t j = 0; j < Tree::get_taxas_num(); j++) {
		t2_leaves_ranks[t2_tr->taxas[j]] = j;
	}

	int* min_rank_in_t2,* max_rank_in_t2;
	compute_start_stop(tree1, tree2, t2_leaves_ranks, min_rank_in_t2, max_rank_in_t2);

	// mark clusters in t1 to be deleted if a heavier incompatible cluster is in t2
	bool* to_del = new bool[tree1->get_nodes_num()];
	bool* marked = new bool[tree2->get_nodes_num()];
	std::fill(to_del, to_del+tree1->get_nodes_num(), false);
	for (size_t i = 1; i < tree1->get_nodes_num(); i++) {
		Tree::Node* node = tree1->get_node(i);
		if (node->is_leaf())
			continue;

		std::fill(marked, marked+tree2->get_nodes_num(), false);
		Tree::Node* leftmost_leaf = tree2->get_leaf(t2_tr->taxas[min_rank_in_t2[i]]);
		Tree::Node* rightmost_leaf = tree2->get_leaf(t2_tr->taxas[max_rank_in_t2[i]]);
		int lcaX = lca(t2_lcas, leftmost_leaf->get_id(), rightmost_leaf->get_id());
		marked[lcaX] = true;
		for (int j = t1_tr->intervals[i].start; j <= t1_tr->intervals[i].end; j++) {
			marked[tree2->get_leaf(t1_tr->taxas[j])->get_id()] = true;
		}

		for (int j = t2_tr->intervals[lcaX].start; j <= t2_tr->intervals[lcaX].end; j++) {
			Tree::Node* leaf = tree2->get_leaf(t2_tr->taxas[j]);
			if (marked[leaf->get_id()]) { // leaf in X, unmark and move on
				marked[leaf->get_id()] = false;
			} else { // leaf not in X, mark all ancestors up to lcaX
				Tree::Node* curr = leaf->get_parent();
				while (!marked[curr->get_id()]) {
					marked[curr->get_id()] = true;
					if (node->get_weight() <= tree2->get_node(curr->get_id())->get_weight()) {
						// weight of node i in t1 is lower than an incompatible cluster in t2
						// mark for deletion
						to_del[i] = true;
						break;
					}
					curr = curr->get_parent();

				}
			}
			if (to_del[i])
				break;
		}
	}

	return to_del;
}

void compute_m(Tree::Node* node, int* e, int* m, std::vector<Tree::Node*>* rsort_lists) {
	if (node->is_leaf()) {
		m[node->get_id()] = e[node->get_taxa()];
	}
	for (Tree::Node* child : node->get_children()) {
		compute_m(child, e, m, rsort_lists);
		if (m[node->get_id()] > m[child->get_id()]) {
			m[node->get_id()] = m[child->get_id()];
		}
	}
	if (!node->is_root()) {
		rsort_lists[m[node->get_id()]].push_back(node);
	}
}


// See Section 2.4 of paper [XXX]
void merge_trees(Tree* tree1, Tree* tree2, taxas_ranges_t* t1_tr, taxas_ranges_t* t2_tr, lca_t* t2_lcas) {
	//compute e
	int* e = new int[Tree::get_taxas_num()];
	for (size_t j = 0; j < Tree::get_taxas_num(); j++) {
		e[t1_tr->taxas[j]] = j;
	}
	int* m = new int[tree2->get_nodes_num()];
	std::fill(m, m+tree2->get_nodes_num(), INT32_MAX);
	std::vector<Tree::Node*>* rsort_lists = new std::vector<Tree::Node*>[Tree::get_taxas_num()];
	compute_m(tree2->get_root(), e, m, rsort_lists);

	// sort tree2
	for (size_t i = 0; i < tree2->get_nodes_num(); i++) {
		tree2->get_node(i)->clear_children();
	}
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) {
		for (auto it = rsort_lists[i].begin(); it != rsort_lists[i].end(); it++) {
			(*it)->get_parent()->add_child(*it);
		}
	}
	t2_tr = build_taxas_ranges(tree2);

	int* start,* stop;
	int* t2_leaves_ranks = new int[Tree::get_taxas_num()];
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) {
		t2_leaves_ranks[t2_tr->taxas[i]] = i;
	}
	compute_start_stop(tree1, tree2, t2_leaves_ranks, start, stop);

	// calc depth
	int* depths = new int[tree2->get_nodes_num()];
	depths[0] = 0; // root depth
	for (size_t i = 1; i < tree2->get_nodes_num(); i++) {
		depths[i] = 1 + depths[tree2->get_node(i)->get_parent()->get_id()];
	}

	// calc x_left and x_right
	Tree::Node** left = new Tree::Node*[Tree::get_taxas_num()];
	Tree::Node** right = new Tree::Node*[Tree::get_taxas_num()];
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) {
		Tree::Node* curr = tree2->get_leaf(i);
		Tree::Node* parent = curr->get_parent();
		while (parent != NULL && *(parent->get_children().begin()) == curr) {
			curr = parent;
			parent = curr->get_parent();
		}
		left[i] = curr;

		curr = tree2->get_leaf(i);
		parent = curr->get_parent();
		while (parent != NULL && *(parent->get_children().rbegin()) == curr) {
			curr = parent;
			parent = curr->get_parent();
		}
		right[i] = curr;
	}

	size_t* orig_pos_in_parent = new size_t[tree2->get_nodes_num()];
	for (size_t i = 0; i < tree2->get_nodes_num(); i++) {
		orig_pos_in_parent[i] = tree2->get_node(i)->get_pos_in_parent();
	}

	for (int i = tree1->get_nodes_num()-1; i >= 1; i--) {
		Tree::Node* a = tree2->get_leaf(t2_tr->taxas[start[i]]);
		Tree::Node* b = tree2->get_leaf(t2_tr->taxas[stop[i]]);
		if (a == b) continue; // tree1->node i is a leaf

		Tree::Node* ru = tree2->get_node(lca(t2_lcas, a->get_id(), b->get_id()));

		Tree::Node* a_left = left[a->get_taxa()];
		Tree::Node* b_right = right[b->get_taxa()];

		size_t du_pos = (depths[a_left->get_id()] > depths[ru->get_id()]) ?
				orig_pos_in_parent[a_left->get_id()] : 0;
		size_t eu_pos = (depths[b_right->get_id()] > depths[ru->get_id()]) ?
				orig_pos_in_parent[b_right->get_id()] : ru->get_children().size()-1;
		if (du_pos == 0 && eu_pos == ru->get_children().size()-1) continue;

		Tree::Node* newnode = tree2->add_node();
		for (size_t j = du_pos; j <= eu_pos; j++) {
			if (ru->get_child(j) != NULL) {
				newnode->add_child(ru->get_child(j));
				// lazy child deletion
				ru->null_child(j);
			}
		}
		ru->set_child(newnode, du_pos);
	}
	tree2->fix_tree();
}

Tree* freqdiff(std::vector<Tree*>& trees) {
	// weights[i][node id] = weights of cluster of node (with node id) in tree i
	// initialize all to "number of trees" because trivial clusters will have that value
	for (size_t i = 0; i < trees.size(); i++) {
		for (size_t j = 0; j < trees[i]->get_nodes_num(); j++) {
			trees[i]->get_node(j)->set_weight(trees.size());
		}
	}
	calc_w_kn2(trees);

	Tree* T = trees[0];
	for (size_t i = 1; i < trees.size(); i++) {
		taxas_ranges_t* tr_T = build_taxas_ranges(T);
		taxas_ranges_t* tr_Ti = build_taxas_ranges(trees[i]);

		lca_t* lca_T = lca_preprocess(T);
		lca_t* lca_Ti = lca_preprocess(trees[i]);

		// filter clusters
		bool* to_del_ti = filter_clusters(trees[i], T, tr_Ti, tr_T, lca_T);
		bool* to_del_t = filter_clusters(T, trees[i], tr_T, tr_Ti, lca_Ti);
		trees[i]->delete_nodes(to_del_ti);
		T->delete_nodes(to_del_t);

		lca_T = lca_preprocess(T);
		tr_T = build_taxas_ranges(T);
		tr_Ti = build_taxas_ranges(trees[i]);

		merge_trees(trees[i], T, tr_Ti, tr_T, lca_T);
	}

	return T;
}


#endif /* FREQDIFF_H_ */
