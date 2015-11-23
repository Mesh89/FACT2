/*
 * freqdiff.h
 *
 *  Created on: 21 Oct 2015
 *      Author: Mesh
 */

#ifndef FREQDIFF_H_
#define FREQDIFF_H_

#include <iostream>
#include <queue>
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

int* start,* stop;
int* e,* m;
std::vector<Tree::Node*>* rsort_lists;

bool* marked;

int* depths;
Tree::Node** left,** right;
size_t* orig_pos_in_parent;

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
				for (int j = tr->intervals[node->id].start; j <= tr->intervals[node->id].end; j++) {
					bitvector[tr->taxas[j]] = 1;
				}
				node_bitvec_t node_bitvec(t, node->id, bitvector);
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
		bit0c = bit1c = 0;
	}
	// bit0 contains sorted bitvectors

	// count adjacent equal bitvectors and set weights
	for (int i = 0; i < bitvc; ) {
		int w = 1;
		while (i+w < bitvc && equal(bit0[i].bitvec,bit0[i+w].bitvec,n)) w++;
		for (int j = i; j < i+w; j++) {
			trees[bit0[j].tree_id]->get_node(bit0[j].node_id)->weight = w;
		}
		i += w;
	}

	delete[] bit0;
	delete[] bit1;
}


// for node_id in T1:
// start[node_id] = smallest rank (i.e. leftmost leaf) in T2 among \Lambda(T1[node_id])
// stop[node_id] = biggest rank (i.e. rightmost leaf) in T2 among \Lambda(T1[node_id])
void compute_start_stop(Tree* tree1, Tree* tree2, int* t2_leaves_ranks) {
	// for each cluster in t1, find its leftmost and rightmost leaves in t2
	for (int i = tree1->get_nodes_num()-1; i >= 0; i--) {
		Tree::Node* node = tree1->get_node(i);
		if (node->is_leaf()) {
			start[i] = stop[i] = t2_leaves_ranks[node->taxa];
		} else {
			start[i] = INT32_MAX;
			stop[i] = 0;
			for (Tree::Node* child : node->children) {
				if (start[i] > start[child->id])
					start[i] = start[child->id];
				if (stop[i] < stop[child->id])
					stop[i] = stop[child->id];
			}
		}
	}
}


void filter_clusters_n2(Tree* tree1, Tree* tree2, taxas_ranges_t* t1_tr, taxas_ranges_t* t2_tr, lca_t* t2_lcas,
		bool* to_del) {

	compute_start_stop(tree1, tree2, t2_tr->ranks);

	// mark clusters in t1 to be deleted if a heavier incompatible cluster is in t2
	std::fill(to_del, to_del+tree1->get_nodes_num(), false);
	for (size_t i = 1; i < tree1->get_nodes_num(); i++) {
		Tree::Node* node = tree1->get_node(i);
		if (node->is_leaf())
			continue;

		std::fill(marked, marked+tree2->get_nodes_num(), false);
		Tree::Node* leftmost_leaf = tree2->get_leaf(t2_tr->taxas[start[i]]);
		Tree::Node* rightmost_leaf = tree2->get_leaf(t2_tr->taxas[stop[i]]);
		int lcaX = lca(t2_lcas, leftmost_leaf->id, rightmost_leaf->id);
		marked[lcaX] = true;

		for (int j = t1_tr->intervals[i].start; j <= t1_tr->intervals[i].end; j++) {
			//mark all ancestors
			Tree::Node* curr = tree2->get_leaf(t1_tr->taxas[j]);
			while (!marked[curr->id]) {
				marked[curr->id] = true;
				curr = curr->parent;
			}
		}

		for (int j = t2_tr->intervals[lcaX].start; j <= t2_tr->intervals[lcaX].end; j++) {
			Tree::Node* leaf = tree2->get_leaf(t2_tr->taxas[j]);
			if (!marked[leaf->id]) { // leaf not in X, mark all ancestors up to lcaX
				Tree::Node* curr = leaf->parent;
				while (curr->id != lcaX) {
					if (marked[curr->id] && node->weight <= curr->weight) {
						// weight of node i in t1 is lower than an incompatible cluster in t2
						// mark for deletion
						to_del[i] = true;
						break;
					}
					curr = curr->parent;

				}
			}
			if (to_del[i])
				break;
		}
	}
}


void filter_clusters_nlog2n(Tree* tree1, Tree* tree2, taxas_ranges_t* t1_tr, taxas_ranges_t* t2_tr, lca_t* t2_lcas,
		bool* to_del) {

	size_t* counter = new size_t[Tree::get_taxas_num()*2];
	std::fill(counter, counter+Tree::get_taxas_num()*2, 0);
	bool* BT = new bool[Tree::get_taxas_num()*2];
	std::fill(BT, BT+Tree::get_taxas_num()*2, false);
	std::priority_queue<std::pair<int, int> > BTw;

	Tree::Node* curr = tree1->get_root();
	while (!curr->is_leaf()) {
		curr = curr->children[0];
	}
	Tree::Node* rim1 = tree2->get_leaf(curr->taxa);
	curr = curr->parent;
	// curr is now p_2 (see [XXX])

	compute_start_stop(tree1, tree2, t2_tr->ranks);
	while (curr != NULL) {
		Tree::Node* ri = tree2->get_node(
				lca(t2_lcas, tree2->get_leaf(t2_tr->taxas[start[curr->id]])->id, tree2->get_leaf(t2_tr->taxas[stop[curr->id]])->id)
				);

		while (rim1 != ri) {
			BT[rim1->id] = true;
			BTw.push(std::make_pair(rim1->weight, rim1->id));
			rim1 = rim1->parent;
		}

		int Dstart = start[curr->id] + curr->children[0]->size,
			Dend = stop[curr->id];
		for (int i = Dstart; i <= Dend; i++) {
			Tree::Node* x = tree2->get_leaf(t2_tr->taxas[i]);
			while (!BT[x->id] && x != ri) {
				BT[x->id] = true;
				BTw.push(std::make_pair(x->weight, x->id));
				x = x->parent;
			}
		}

		for (int i = Dstart; i <= Dend; i++) {
			Tree::Node* x = tree2->get_leaf(t2_tr->taxas[i]);
			counter[x->id]++;
			while (x != NULL && counter[x->id] == x->size) {
				counter[x->parent->id] += x->size;
				BT[x->id] = false;
				x = x->parent;
			}
		}

		std::pair<int, int> top;
		while (!BTw.empty()) {
			top = BTw.top();
			BTw.pop();
			if (BT[top.second]) break;
		}
		int M = BTw.empty() ? 0 : top.first;

		// TODO: beta_i stuff

		if (curr->weight > M) {
			// TODO: insert
		}

		curr = curr->parent;
	}
}


void compute_m(Tree::Node* node, int* e, int* m, std::vector<Tree::Node*>* rsort_lists) {
	if (node->is_leaf()) {
		m[node->id] = e[node->taxa];
	}
	for (Tree::Node* child : node->children) {
		compute_m(child, e, m, rsort_lists);
		if (m[node->id] > m[child->id]) {
			m[node->id] = m[child->id];
		}
	}
	if (!node->is_root()) {
		rsort_lists[m[node->id]].push_back(node);
	}
}


// See Section 2.4 of paper [XXX]
void merge_trees(Tree* tree1, Tree* tree2, taxas_ranges_t* t1_tr, lca_t* t2_lcas) {
	//compute e
	for (size_t j = 0; j < Tree::get_taxas_num(); j++) {
		e[t1_tr->taxas[j]] = j;
	}
	std::fill(m, m+tree2->get_nodes_num(), INT32_MAX);
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) rsort_lists[i].clear();

	compute_m(tree2->get_root(), e, m, rsort_lists);

	// sort tree2
	for (size_t i = 0; i < tree2->get_nodes_num(); i++) {
		tree2->get_node(i)->clear_children();
	}
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) {
		for (auto it = rsort_lists[i].begin(); it != rsort_lists[i].end(); it++) {
			(*it)->parent->add_child(*it);
		}
	}

	taxas_ranges_t* t2_tr = build_taxas_ranges(tree2);
	compute_start_stop(tree1, tree2, t2_tr->ranks);

	// calc depth
	depths[0] = 0; // root depth
	for (size_t i = 1; i < tree2->get_nodes_num(); i++) {
		depths[i] = 1 + depths[tree2->get_node(i)->parent->id];
	}

	// calc x_left and x_right
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) {
		Tree::Node* curr = tree2->get_leaf(i);
		Tree::Node* parent = curr->parent;
		while (parent != NULL && *(parent->children.begin()) == curr) {
			curr = parent;
			parent = curr->parent;
		}
		left[i] = curr;

		curr = tree2->get_leaf(i);
		parent = curr->parent;
		while (parent != NULL && *(parent->children.rbegin()) == curr) {
			curr = parent;
			parent = curr->parent;
		}
		right[i] = curr;
	}

	for (size_t i = 0; i < tree2->get_nodes_num(); i++) {
		orig_pos_in_parent[i] = tree2->get_node(i)->pos_in_parent;
	}

	for (int i = tree1->get_nodes_num()-1; i >= 1; i--) {
		Tree::Node* a = tree2->get_leaf(t2_tr->taxas[start[i]]);
		Tree::Node* b = tree2->get_leaf(t2_tr->taxas[stop[i]]);
		if (a == b) continue; // tree1->node i is a leaf

		Tree::Node* ru = tree2->get_node(lca(t2_lcas, a->id, b->id));

		Tree::Node* a_left = left[a->taxa];
		Tree::Node* b_right = right[b->taxa];

		size_t du_pos = (depths[a_left->id] > depths[ru->id]) ?
				orig_pos_in_parent[a_left->id] : 0;
		size_t eu_pos = (depths[b_right->id] > depths[ru->id]) ?
				orig_pos_in_parent[b_right->id] : ru->get_children_num()-1;
		if (du_pos == 0 && eu_pos == ru->get_children_num()-1) continue;

		Tree::Node* newnode = tree2->add_node();
		newnode->weight = tree1->get_node(i)->weight;
		newnode->size = tree1->get_node(i)->size;
		for (size_t j = du_pos; j <= eu_pos; j++) {
			if (ru->children[j] != NULL) {
				newnode->add_child(ru->children[j]);
				// lazy child deletion
				ru->null_child(j);
			}
		}
		ru->set_child(newnode, du_pos);
	}
	tree2->fix_tree();
}

Tree* freqdiff(std::vector<Tree*>& trees) {

	start = new int[Tree::get_taxas_num()*2];
	stop = new int[Tree::get_taxas_num()*2];
	e = new int[Tree::get_taxas_num()];
	m = new int[Tree::get_taxas_num()*2];
	rsort_lists = new std::vector<Tree::Node*>[Tree::get_taxas_num()];

	marked = new bool[Tree::get_taxas_num()*2];
	bool* to_del_t = new bool[Tree::get_taxas_num()*2];
	bool* to_del_ti = new bool[Tree::get_taxas_num()*2];

	depths = new int[Tree::get_taxas_num()*2];
	left = new Tree::Node*[Tree::get_taxas_num()];
	right = new Tree::Node*[Tree::get_taxas_num()];
	orig_pos_in_parent = new size_t[Tree::get_taxas_num()*2];

	// weights[i][node id] = weights of cluster of node (with node id) in tree i
	// initialize all to "number of trees" because trivial clusters will have that value
	for (size_t i = 0; i < trees.size(); i++) {
		for (size_t j = 0; j < trees[i]->get_nodes_num(); j++) {
			trees[i]->get_node(j)->weight = trees.size();
		}
	}
	calc_w_kn2(trees);

	lca_t** lca_preps = new lca_t*[trees.size()];
	for (size_t i = 0; i < trees.size(); i++) {
		lca_preps[i] = lca_preprocess(trees[i]);
	}

	Tree* T = new Tree(trees[0]); //trees[0]->copy();
	for (size_t i = 1; i < trees.size(); i++) {
		Tree* Ti = new Tree(trees[i]); //trees[i]->copy();
		taxas_ranges_t* tr_T = build_taxas_ranges(T);
		taxas_ranges_t* tr_Ti = build_taxas_ranges(Ti);

		lca_t* lca_T = lca_preprocess(T);

		// filter clusters
		filter_clusters_n2(Ti, T, tr_Ti, tr_T, lca_T, to_del_ti);
		filter_clusters_nlog2n(Ti, T, tr_Ti, tr_T, lca_T, to_del_ti);
		filter_clusters_n2(T, Ti, tr_T, tr_Ti, lca_preps[i], to_del_t);
		filter_clusters_nlog2n(T, Ti, tr_T, tr_Ti, lca_preps[i], to_del_t);
		Ti->delete_nodes(to_del_ti);
		T->delete_nodes(to_del_t);

		lca_T = lca_preprocess(T);
		tr_Ti = build_taxas_ranges(Ti);

		merge_trees(Ti, T, tr_Ti, lca_T);

		// delete Ti;
	}

	for (size_t i = 0; i < trees.size(); i++) {
		taxas_ranges_t* tr_T = build_taxas_ranges(T);
		taxas_ranges_t* tr_Ti = build_taxas_ranges(trees[i]);

		filter_clusters_n2(T, trees[i], tr_T, tr_Ti, lca_preps[i], to_del_t);
		T->delete_nodes(to_del_t);
	}

	delete[] start;
	delete[] stop;
	delete[] e;
	delete[] m;
	delete[] rsort_lists;

	delete[] marked;
	delete[] to_del_ti;
	delete[] to_del_t;

	delete[] orig_pos_in_parent;
	delete[] left;
	delete[] right;

	return T;
}


#endif /* FREQDIFF_H_ */
