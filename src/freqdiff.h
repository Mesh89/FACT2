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
#include <boost/dynamic_bitset.hpp>

#include "taxas_ranges.h"
#include "lca_preprocessing.h"
#include "utils.h"

#include "Tree.h"

struct node_bitvec_t {
	Tree::Node* node;
	int bitvec_id;

	node_bitvec_t() : node(NULL), bitvec_id(-1) {}
	node_bitvec_t(Tree::Node* node, int bitvec_id) : node(node), bitvec_id(bitvec_id) {}
};

struct subpath_query_info_t {
	Tree::Node** cp_roots;
	gen_rmq_t** cp_rmqs;

	subpath_query_info_t(Tree::Node** cp_roots, gen_rmq_t** cp_rmqs) : cp_roots(cp_roots), cp_rmqs(cp_rmqs) {}
};

int* start,* stop;
int* e,* m;
std::vector<Tree::Node*>* rsort_lists;

bool* marked;

int* depths;
Tree::Node** left,** right;
size_t* orig_pos_in_parent;

size_t* counter;
bool* BT;
int* leaf_p_index;


// calculate cluster weights using kn^2 method
void calc_w_kn2(std::vector<Tree*>& trees) {

	int tot_int_nodes = 0;
	for (Tree* tree : trees) {
		tot_int_nodes += tree->get_nodes_num() - Tree::get_taxas_num() - 1;
	}

	// Generate bit vectors
	node_bitvec_t* node_bitvecs = new node_bitvec_t[tot_int_nodes];
	boost::dynamic_bitset<>* bitvectors = new boost::dynamic_bitset<>[tot_int_nodes];
	int bitvc = 0;
	size_t n = Tree::get_taxas_num();
	for (size_t t = 0; t < trees.size(); t++) {
		Tree* tree = trees[t];
		taxas_ranges_t* tr = build_taxas_ranges(tree);

		size_t nodes_num = tree->get_nodes_num();
		for (size_t i = 1; i < nodes_num; i++) {
			Tree::Node* node = tree->get_node(i);
			if (!node->is_leaf()) {
				// fill bit vector for current node
				node_bitvec_t node_bitvec(node, bitvc);
				//bitvectors.emplace_back(Tree::get_taxas_num());
				bitvectors[bitvc] = boost::dynamic_bitset<>(Tree::get_taxas_num());
				for (int j = tr->intervals[node->id].start; j <= tr->intervals[node->id].end; j++) {
					bitvectors[bitvc].set(tr->taxas[j]);
				}
				node_bitvecs[bitvc] = node_bitvec;
				bitvc++;
			}
		}
	}

	// Sort bit vectors
	node_bitvec_t* bit0 = node_bitvecs;//new node_bitvec_t[bitvc];
	node_bitvec_t* bit1 = new node_bitvec_t[bitvc];
	int bit0c = 0, bit1c = 0;

	for (int i = n-1; i >= 0; i--) {
		for (int j = 0; j < bitvc; j++) {
			if (bitvectors[bit0[j].bitvec_id].test(i)) {
				bit1[bit1c++] =  bit0[j];
				// we signal that this position is free as the occupant had the curr bit set to 1
				bit0[j].node = NULL;
			}
		}
		for (int j = 0; j < bitvc; j++) {
			if (bit0[j].node != NULL) {
				bit0[bit0c++] = bit0[j];
			}
		}
		for (int j = 0; j < bit1c; j++) {
			if (bit1[j].node != NULL) {
				bit0[bit0c++] = bit1[j];
			}
		}
		bit0c = bit1c = 0;
	}
	// bit0 contains sorted bitvectors

	// count adjacent equal bitvectors and set weights
	for (int i = 0; i < bitvc; ) {
		int w = 1;
		while (i+w < bitvc && bitvectors[bit0[i].bitvec_id] == bitvectors[bit0[i+w].bitvec_id]) w++;
		for (int j = i; j < i+w; j++) {
			bit0[j].node->weight = w;
		}
		i += w;
	}

	delete[] bitvectors;
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

	int* t2_ranks = get_taxas_ranks(t2_tr);
	compute_start_stop(tree1, tree2, t2_ranks);

	// mark clusters in t1 to be deleted if a heavier incompatible cluster is in t2
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


subpath_query_info_t* preprocess_subpaths_queries(Tree* tree) {
	// pointer from each node to the root of its cp
	Tree::Node** cp_roots = new Tree::Node*[tree->get_nodes_num()];
	cp_roots[0] = tree->get_root();
	for (size_t i = 1; i < tree->get_nodes_num(); i++) {
		Tree::Node* node = tree->get_node(i);
		if (node->pos_in_parent == 0) {
			cp_roots[i] = cp_roots[node->parent->id];
		} else {
			cp_roots[i] = node;
		}
	}

	gen_rmq_t** cp_rmqs = new gen_rmq_t*[tree->get_nodes_num()];
	for (size_t i = 0; i < tree->get_nodes_num(); i++) {
		cp_rmqs[i] = new gen_rmq_t;
		cp_rmqs[cp_roots[i]->id]->v.push_back(-tree->get_node(i)->weight);
	}
	for (size_t i = 0; i < tree->get_nodes_num(); i++) {
		if (cp_rmqs[i]->v.size() > 1) {
			general_rmq_preprocess(cp_rmqs[i]); // TODO: the rmq API needs some revisiting
		}
	}

	//int* depths = new int[tree->get_nodes_num()]; //TODO: think about externalizing depths
	depths[0] = 0; // root depth
	for (size_t i = 1; i < tree->get_nodes_num(); i++) {
		depths[i] = 1 + depths[tree->get_node(i)->parent->id];
	}

	/*rmq_t** leaf_rmqs = new rmq_t*[Tree::get_taxas_num()]; TODO
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) {
		Tree::Node* curr = tree->get_leaf(i)->parent;
		leaf_rmqs[i] = new rmq_t;
		while (curr != NULL) {
			gen_rmq_t* currpath_rmq = cp_rmqs[cp_roots[curr->id]->id];
			int query_endp = depths[curr->id] - depths[cp_roots[curr->id]->id];
			leaf_rmqs[i]->v.push_back(currpath_rmq->v[general_rmq(currpath_rmq, 0, query_endp)]);
			curr = cp_roots[curr->id]->parent;
		}
		if (leaf_rmqs[i]->v.size() > 1) { // FIXME: if size(v) is 1, rmq_preprocess crashes.
			rmq_preprocess(leaf_rmqs[i], leaf_rmqs[i]->v);
		}
	}*/

	return new subpath_query_info_t(cp_roots, cp_rmqs);
}

int max_subpath_query(subpath_query_info_t* subpq_info, Tree::Node* ancestor, Tree::Node* descendant) { // TODO: we might need only ids
	Tree::Node* curr = descendant;
	int res = 0;
	while (subpq_info->cp_roots[curr->secondary_id] != subpq_info->cp_roots[ancestor->secondary_id]) {
		gen_rmq_t* currpath_rmq = subpq_info->cp_rmqs[subpq_info->cp_roots[curr->secondary_id]->id];
		int query_endp = depths[curr->secondary_id] - depths[subpq_info->cp_roots[curr->secondary_id]->id];
		if (query_endp == 0) {
			res = std::min(res, currpath_rmq->v[0]);
		} else {
			res = std::min(res, currpath_rmq->v[general_rmq(currpath_rmq, 0, query_endp)]);
		}
		curr = subpq_info->cp_roots[curr->secondary_id]->parent;
	}
	gen_rmq_t* currpath_rmq = subpq_info->cp_rmqs[subpq_info->cp_roots[curr->secondary_id]->id];
	int query_startp = depths[ancestor->secondary_id] - depths[subpq_info->cp_roots[curr->secondary_id]->id];
	int query_endp = depths[curr->secondary_id] - depths[subpq_info->cp_roots[curr->secondary_id]->id];
	if (query_startp < query_endp) {
		res = std::min(res, currpath_rmq->v[general_rmq(currpath_rmq, query_startp+1, query_endp)]);
	}
	return -res;
}


subpath_query_info_t* subpq_info; // TODO: temporary
lca_t* lca_prep;
// leaves in marked must be sorted in left-to-right order
Tree* contract_tree_fast(Tree* tree, std::vector<int>& marked) {	//TODO: think about externalizing depths and lca_t
	if (marked.empty()) return NULL;

	depths[0] = 0; // root depth
	for (size_t i = 1; i < tree->get_nodes_num(); i++) {
		depths[i] = 1 + depths[tree->get_node(i)->parent->id];
	}

	std::vector<int> levels, ids;
	levels.push_back(depths[tree->get_leaf(marked[0])->id]);
	ids.push_back(tree->get_leaf(marked[0])->secondary_id);
	for (size_t i = 1; i < marked.size(); i++) {
		int lca_id = lca(lca_prep, tree->get_leaf(marked[i-1])->id, tree->get_leaf(marked[i])->id);
		levels.push_back(depths[lca_id]);
		ids.push_back(tree->get_node(lca_id)->secondary_id);
		levels.push_back(depths[tree->get_leaf(marked[i])->id]);
		ids.push_back(tree->get_leaf(marked[i])->secondary_id);
	}

	std::vector<int> vleft(levels.size(), -1), pointer(levels.size(), -1);
	for (size_t i = 0; i < levels.size(); i++) pointer[i] = i;

	std::stack<int> Sl;
	for (size_t i = 0; i < levels.size(); i++) {
		while (!Sl.empty() && levels[i] <= levels[Sl.top()]) {
			Sl.pop();
		}

		if (!Sl.empty()) {
			vleft[i] = Sl.top();
		}
		Sl.push(i);
	}

	std::vector<int> vright(levels.size(), -1);
	std::stack<int> Sr;
	for (int i = levels.size()-1; i >= 0; i--) {
		while (!Sr.empty() && levels[i] <= levels[Sr.top()]) {
			Sr.pop();
		}

		if (!Sr.empty()) {
			vright[i] = Sr.top();
		}
		Sr.push(i);
	}

	int root_pos = -1;
	std::vector<int> parent(levels.size(), -1);
	std::vector<bool> exists(levels.size(), false);
	for (size_t i = 0; i < levels.size(); i++) {
		if (root_pos == -1 && vleft[i] == -1 && vright[i] == -1) { //root
			root_pos = i;
		} else if (vleft[i] == -1) {
			parent[i] = pointer[vright[i]];
		} else if (vright[i] == -1) {
			parent[i] = pointer[vleft[i]];
		} else {
			if (levels[vleft[i]] >= levels[vright[i]]) {
				parent[i] = pointer[vleft[i]];
				if (levels[vleft[i]] == levels[vright[i]]) { // deals with non-binary nodes
					pointer[vright[i]] = pointer[vleft[i]];
				}
			} else {
				parent[i] = pointer[vright[i]];
			}
		}

		if (parent[i] >= 0) {
			exists[parent[i]] = true; // non-binary nodes create trash nodes, which should not actually be instantiated
		}
	}

	assert(root_pos >= 0);

	Tree* new_tree = new Tree();
	std::vector<Tree::Node*> tree_nodes(levels.size(), NULL);
	for (size_t i = 0; i < levels.size(); i++) {
		// create node
		if (i%2 == 0 || exists[i]) {
			if (tree_nodes[i] == NULL) {
				if (i%2 == 0) { //a leaf
					tree_nodes[i] = new_tree->add_node(marked[i/2]);
				} else {
					tree_nodes[i] = new_tree->add_node();
				}
			}
			tree_nodes[i]->weight = tree->get_node(ids[i])->weight;
			tree_nodes[i]->secondary_id = ids[i];
		}

		// connect to parent
		if (tree_nodes[i] != NULL && parent[i] != -1) {
			if (tree_nodes[parent[i]] == NULL) {
				tree_nodes[parent[i]] = new_tree->add_node();
			}
			tree_nodes[parent[i]]->add_child(tree_nodes[i]);
		}
	}

	new_tree->fix_tree(tree_nodes[root_pos]); // TODO: could it be done only once at the end?
	assert(tree_nodes[root_pos]->parent == NULL);

	// insert special nodes
	size_t newtree_nodes = new_tree->get_nodes_num();
	for (size_t i = 0; i < newtree_nodes; i++) {
		Tree::Node* curr_node = new_tree->get_node(i);

		if (curr_node->is_root()) continue;

		Tree::Node* desc_par = tree->get_node(curr_node->secondary_id)->parent;
		int msq = max_subpath_query(subpq_info, curr_node->parent, desc_par);
		if (msq > 0) {
			Tree::Node* sp_node = new_tree->add_node();
			sp_node->weight = msq;
			sp_node->secondary_id = curr_node->parent->secondary_id;
			curr_node->parent->null_child(curr_node->pos_in_parent);
			curr_node->parent->add_child(sp_node);
			sp_node->add_child(curr_node);
		}
	}
	new_tree->fix_tree(tree_nodes[root_pos]);

	return new_tree;
}


Tree* orig_t2 = NULL; // TODO: temporary
bool is_special(Tree::Node* node) {
	return node->size < orig_t2->get_node(node->secondary_id)->size;
}
void filter_clusters_nlog2n(Tree::Node* t1_root, Tree* tree2, taxas_ranges_t* t1_tr, lca_t* t2_lcas,
		bool* to_del) {

	//std::fill(leaf_p_index, leaf_p_index+Tree::get_taxas_num(), -1);

	std::priority_queue<std::pair<int, int> > BTw;

	Tree::Node* curr = t1_root;
	while (!curr->is_leaf()) {
		curr = curr->children[0];
	}

	int beta = 0;

	Tree::Node* rim1 = tree2->get_leaf(curr->taxa);
	curr = curr->parent;
	Tree::Node* p_2 = curr;

	/* Handle side trees */
	int p_index = 0;

	// build sorted sets of leaves
	leaf_p_index[t1_tr->taxas[t1_tr->intervals[t1_root->id].start]] = -1; // p_1 must be ignored

	std::vector<Tree::Node*> st_roots;
	while (curr != t1_root->parent) {
		for (size_t i = 1; i < curr->get_children_num(); i++) {
			Tree::Node* child = curr->children[i];
			for (int j = t1_tr->intervals[child->id].start; j <= t1_tr->intervals[child->id].end; j++) {
				leaf_p_index[t1_tr->taxas[j]] = p_index;
			}
			st_roots.push_back(child);

			p_index++;
		}
		curr = curr->parent;
	}
	taxas_ranges_t* t2_tr = build_taxas_ranges(tree2);
	std::vector<int>* leaves_sets = new std::vector<int>[p_index];
	for (size_t i = 0; i < t2_tr->taxas_num; i++) {
		if (leaf_p_index[t2_tr->taxas[i]] >= 0) {
			leaves_sets[leaf_p_index[t2_tr->taxas[i]]].push_back(t2_tr->taxas[i]);
		} else {
			assert(t2_tr->taxas[i] == t1_tr->taxas[t1_tr->intervals[t1_root->id].start]);
		}
	}

	// build restrictions and recurse
	for (int i = 0; i < p_index; i++) {
		if (leaves_sets[i].size() > 1) {
			Tree* contracted = contract_tree_fast(orig_t2, leaves_sets[i]);
			lca_t* lca_ct2 = lca_preprocess(contracted);
			filter_clusters_nlog2n(st_roots[i], contracted, t1_tr, lca_ct2, to_del);
		}
	}

	std::fill(BT, BT+tree2->get_nodes_num(), false);
	std::fill(counter, counter+tree2->get_nodes_num(), 0);

	// T_2[p_1] should is "initialized" separately
	counter[rim1->parent->id] = 1;
	rim1 = rim1->parent;

	/* Handle centroid path */
	int curr_t1_start = t1_tr->intervals[t1_root->id].start + 1;
	curr = p_2;
	while (curr != t1_root->parent) {
		int curr_t1_end = t1_tr->intervals[curr->id].end;

		Tree::Node* ri = rim1;
		for (int i = curr_t1_start; i <= curr_t1_end; i++) {
			ri = tree2->get_node(lca(t2_lcas, ri->id, tree2->get_leaf(t1_tr->taxas[i])->id));
		}

		if (counter[rim1->id] == rim1->size && ri != rim1 && rim1->get_children_num() != 1) {
			// if children of rim1 are all "marked", should not go into BT
			counter[rim1->parent->id] += rim1->size;
			if (rim1->weight > beta && is_special(rim1)) beta = rim1->weight;
			rim1 = rim1->parent;
		}
		while (rim1 != ri) {
			BT[rim1->id] = true;
			BTw.push(std::make_pair(rim1->weight, rim1->id));
			if (rim1->weight > beta && is_special(rim1)) beta = rim1->weight;
			rim1 = rim1->parent;
		}

		for (int i = curr_t1_start; i <= curr_t1_end; i++) {
			Tree::Node* x = tree2->get_leaf(t1_tr->taxas[i]);
			while (!BT[x->id] && x != ri) {
				BT[x->id] = true;
				BTw.push(std::make_pair(x->weight, x->id));
				if (x->weight > beta && is_special(x)) beta = x->weight; // TODO: and is special node
				x = x->parent;
			}
		}
		for (int i = curr_t1_start; i <= curr_t1_end; i++) {
			Tree::Node* x = tree2->get_leaf(t1_tr->taxas[i]);
			counter[x->id]++;
			while (x != ri && counter[x->id] == x->size && !is_special(x)) { // TODO: check
				counter[x->parent->id] += x->size;
				BT[x->id] = false;
				x = x->parent;
			}
		}

		std::pair<int, int> top;
		while (!BTw.empty()) {

			top = BTw.top();
			if (BT[top.second]) break;
			BTw.pop();
		}
		int M = BTw.empty() ? 0 : top.first;
		if (curr->weight <= std::max(M, beta)) {
			to_del[curr->id] = true;
		}

		curr = curr->parent;
		curr_t1_start = curr_t1_end + 1;
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
	int* t2_ranks = get_taxas_ranks(t2_tr);
	compute_start_stop(tree1, tree2, t2_ranks);

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

Tree* freqdiff(std::vector<Tree*>& trees, bool centroid_paths) {

	//if (centroid_paths)
	//	std::cout << "Using centroid paths version" << std::endl;

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

	BT = new bool[Tree::get_taxas_num()*2];
	counter = new size_t[Tree::get_taxas_num()*2];
	leaf_p_index = new int[Tree::get_taxas_num()];

	// weights[i][node id] = weights of cluster of node (with node id) in tree i
	// initialize all to "number of trees" because trivial clusters will have that value
	for (size_t i = 0; i < trees.size(); i++) {
		for (size_t j = 0; j < trees[i]->get_nodes_num(); j++) {
			trees[i]->get_node(j)->weight = trees.size();
		}
		trees[i]->reorder();
	}
	calc_w_kn2(trees);

	lca_t** lca_preps = new lca_t*[trees.size()];
	for (size_t i = 0; i < trees.size(); i++) {
		lca_preps[i] = lca_preprocess(trees[i]);
	}

	Tree* T = new Tree(trees[0]);
	for (size_t i = 1; i < trees.size(); i++) {
		//std::cout << i << std::endl;
		Tree* Ti = new Tree(trees[i]);
		taxas_ranges_t* tr_Ti = build_taxas_ranges(Ti);

		T->reorder();
		taxas_ranges_t* tr_T = build_taxas_ranges(T);

		lca_t* lca_T = lca_preprocess(T);

		// filter clusters
		std::fill(to_del_ti, to_del_ti+Ti->get_nodes_num(), false);
		if (centroid_paths) {
			orig_t2 = T; // TODO: temporary
			subpq_info = preprocess_subpaths_queries(T);
			lca_prep = lca_preprocess(T);
			filter_clusters_nlog2n(Ti->get_root(), T, tr_Ti, lca_T, to_del_ti);
		} else {
			filter_clusters_n2(Ti, T, tr_Ti, tr_T, lca_T, to_del_ti);
		}

		std::fill(to_del_t, to_del_t+T->get_nodes_num(), false);
		if (centroid_paths) {
			orig_t2 = Ti; // TODO: temporary
			lca_prep = lca_preprocess(Ti);
			subpq_info = preprocess_subpaths_queries(Ti);
			filter_clusters_nlog2n(T->get_root(), Ti, tr_T, lca_preps[i], to_del_t);
		} else {
			filter_clusters_n2(T, Ti, tr_T, tr_Ti, lca_preps[i], to_del_t);
		}

		Ti->delete_nodes(to_del_ti);
		T->delete_nodes(to_del_t);

		lca_T = lca_preprocess(T);
		tr_Ti = build_taxas_ranges(Ti);

		merge_trees(Ti, T, tr_Ti, lca_T);

		// reset secondary ids
		for (size_t j = 0; j < T->get_nodes_num(); j++) {
			T->get_node(j)->secondary_id = T->get_node(j)->id;
		}

		//delete Ti;
	}

	for (size_t i = 0; i < trees.size(); i++) {
		taxas_ranges_t* tr_T = build_taxas_ranges(T);
		taxas_ranges_t* tr_Ti = build_taxas_ranges(trees[i]);

		std::fill(to_del_t, to_del_t+T->get_nodes_num(), false);
		if (centroid_paths) {
			orig_t2 = trees[i];
			lca_prep = lca_preprocess(trees[i]);
			subpq_info = preprocess_subpaths_queries(trees[i]); // TODO: should not go here
			filter_clusters_nlog2n(T->get_root(), trees[i], tr_T, lca_preps[i], to_del_t);
		} else {
			filter_clusters_n2(T, trees[i], tr_T, tr_Ti, lca_preps[i], to_del_t);
		}
		T->delete_nodes(to_del_t); // TODO: could be moved outside?
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
