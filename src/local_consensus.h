/*
 * minrs.h
 *
 *  Created on: Feb 14, 2017
 *      Author: ramesh
 */

#ifndef SRC_LOCAL_CONSENSUS_H_
#define SRC_LOCAL_CONSENSUS_H_

#include <vector>
#include <iostream>
#include <bitset>

#include "Tree.h"

typedef unsigned int uint;


uint dfs(std::vector<int>* adjl, int v, bool* visited, std::vector<int>& indices) {
	uint comp = (1 << indices[v]);
	visited[v] = true;
	for (int i = 0; i < (int) adjl[v].size(); i++) {
		if (!visited[adjl[v][i]]) {
			comp |= dfs(adjl, adjl[v][i], visited, indices);
		}
	}
	return comp;
}

std::vector<uint> build_aho(std::vector<int>& indices, int*** triplets, bool minrs) {
	int n = indices.size();
	std::vector<int>* adjl = new std::vector<int>[n];
	for (int i = 0; i < n; i++) {
		int ii = indices[i];
		for (int j = i+1; j < n; j++) {
			int ij = indices[j];
			for (int k = 0; k < n; k++) {
				int ik = indices[k];
				if (triplets[ii][ij][ik]) {
					adjl[i].push_back(j);
					adjl[j].push_back(i);
					break;
				}
			}
		}
	}

	std::vector<uint> components;
	bool* visited = new bool[n];
	std::fill(visited, visited+n, false);
	for (int i = 0; i < n; i++) {
		if (!visited[i]) {
			uint comp = dfs(adjl, i, visited, indices);
			if (std::bitset<32>(comp).count() > minrs) { // if minRS do not consider singletons
				components.push_back(comp);
			}
		}
	}
	delete[] visited;
	delete[] adjl;

	return components;
}


void print_bitmask_as_indices(uint bitmask, int len) {
	std::cout << "( ";
	for (int i = 0; i < len; i++) {
		if (bitmask & (1<<i)) {
			std::cout << Tree::taxa_names[i] << " ";
		}
	}
	std::cout << ")" << std::endl;
}


void print_tree(uint Lbitmask, std::vector<uint>* components, int** dp_backtrack, Tree* tree, Tree::Node* node, bool minrs) {

	int m = components[Lbitmask].size();

	// print singleton elems
	if (minrs) {
		uint singletons = 0;
		for (int i = 0; i < m; i++) {
			singletons |= components[Lbitmask][i];
		}
		singletons = Lbitmask & ~singletons;
		for (uint i = 0; i < Tree::get_taxas_num(); i++) {
			if (singletons & (1<<i)) {
				node->add_child(tree->add_node(i));
			}
		}
	}

	if (m == 0) {
		return;
	} else if (m == 1) {
		Tree::Node* new_node = tree->add_node();
		node->add_child(new_node);
		print_tree(components[Lbitmask][0], components, dp_backtrack, tree, new_node, minrs);
		return;
	}

	uint Dbitmask = (1 << m)-1;
	int Xbitmask = dp_backtrack[Lbitmask][(1 << m)-1];
	int origXbitmask;
	do {
		origXbitmask = Xbitmask;
		Xbitmask = abs(Xbitmask);
		uint DmXbitmask = Dbitmask & ~Xbitmask;
		uint lambdaX = 0;
		for (int i = 0; i < m; i++) {
			if (Xbitmask & (1<<i)) {
				lambdaX |= components[Lbitmask][i];
			}
		}

		if (!minrs && std::bitset<32>(lambdaX).count() == 1) {
			for (uint i = 0; i < Tree::get_taxas_num(); i++) {
				if (lambdaX & (1<<i)) {
					node->add_child(tree->add_node(i));
				}
			}
		} else {
			Tree::Node* new_node = tree->add_node();
			node->add_child(new_node);
			print_tree(lambdaX, components, dp_backtrack, tree, new_node, minrs);
		}

		Dbitmask = DmXbitmask;
		Xbitmask = dp_backtrack[Lbitmask][DmXbitmask];
	} while (origXbitmask >= 0);

	Xbitmask = abs(Xbitmask);
	uint lambdaX = 0;
	for (int i = 0; i < m; i++) {
		if (Dbitmask & (1<<i)) {
			lambdaX |= components[Lbitmask][i];
		}
	}

	if (!minrs && std::bitset<32>(lambdaX).count() == 1) {
		for (uint i = 0; i < Tree::get_taxas_num(); i++) {
			if (lambdaX & (1<<i)) {
				node->add_child(tree->add_node(i));
			}
		}
	} else {
		Tree::Node* new_node = tree->add_node();
		node->add_child(new_node);
		print_tree(lambdaX, components, dp_backtrack, tree, new_node, minrs);
	}
}


int*** get_common_triplets(std::vector<Tree*>& trees) {
	int n = trees[0]->get_taxas_num();
	int*** triplets = new int**[n];
	for (int i = 0; i < n; i++) {
		triplets[i] = new int*[n];
		for (int j = 0; j < n; j++) {
			triplets[i][j] = new int[n];
			std::fill(triplets[i][j], triplets[i][j]+n, 0);
		}
	}

	for (Tree* tree : trees) {
		lca_t* lca_prep = lca_preprocess(tree);
		for (int i = 0; i < n; i++) {
			Tree::Node* leafi = tree->get_leaf(i);
			for (int j = i+1; j < n; j++) {
				Tree::Node* leafj = tree->get_leaf(j);
				int dlcaij = tree->get_node(lca(lca_prep, leafi->id, leafj->id))->depth;
				for (int k = j+1; k < n; k++) {
					Tree::Node* leafk = tree->get_leaf(k);
					int dlcaik = tree->get_node(lca(lca_prep, leafi->id, leafk->id))->depth;
					int dlcajk = tree->get_node(lca(lca_prep, leafj->id, leafk->id))->depth;


					if (dlcaij == dlcaik && dlcaij == dlcajk) continue; // fan

					if (dlcaik == dlcajk && dlcaij > dlcaik) { //ij|k
						triplets[i][j][k]++;
					} else if (dlcaij == dlcajk && dlcaik > dlcaij) { //ik|j
						triplets[i][k][j]++;
					} else if (dlcaij == dlcaik && dlcajk > dlcaij) { //jk|i
						triplets[j][k][i]++;
					} else {
						// should never enter here
						std::cout << dlcaij << " " << dlcaik << " " << dlcajk << std::endl;
						assert(false);
					}
				}
			}
		}
		delete lca_prep;
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				if (triplets[i][j][k] < (int) trees.size()) {
					triplets[i][j][k] = 0;
				}
			}
		}
	}

	return triplets;
}


Tree* minRILC(std::vector<Tree*>& trees, bool minrs) {

	int n = trees[0]->get_taxas_num();
	int*** triplets = get_common_triplets(trees);

	uint states = (1 << n);
	int* opt = new int[states];
	for (int i = 0; i < n; i++) { // init opt for single taxa
		opt[1 << i] = 0;
	}

	int* dp = new int[states];
	dp[0] = 0;

	std::vector<uint>* components = new std::vector<uint>[states];

	// dp_backtrack[Lbitmask][Dbitmask] is the bitmask of X (line 8 in the MCFS paper)
	// if positive, it means DP[D\X] < opt[Merge(D\X)], else the number if negative
	int** dp_backtrack = new int*[states];
	std::fill(dp_backtrack, dp_backtrack+states, (int*)NULL);

	for (int i = 2; i <= n; i++) {
		uint Lbitmask = (1 << i)-1; // L' bitmask
		while ((Lbitmask & (1 << n)) == false) {
			// indices contains the subset of taxa we are currently dealing with
			std::vector<int> indices;
			for (int j = 0; j < n; j++) {
				if (Lbitmask & (1<<j)) {
					indices.push_back(j);
				}
			}

			// aho graph components
			components[Lbitmask] = build_aho(indices, triplets, minrs);
			int m = components[Lbitmask].size();
			if (m == 1 && std::bitset<32>(components[Lbitmask][0]).count() == std::bitset<32>(Lbitmask).count()) {
				return NULL;
			}

			dp_backtrack[Lbitmask] = new int[1 << m];

			// init dp for single components
			for (int j = 0; j < m; j++) { // init opt for single taxa
				dp[1 << j] = opt[components[Lbitmask][j]];
				if (!minrs) {
					std::size_t countComp = std::bitset<32>(components[Lbitmask][j]).count();
					dp[1 << j] += comb2(countComp) * (i-countComp);
				}
			}

			for (int j = 2; j <= m; j++) {
				uint Dbitmask = (1 << j)-1; // D bitmask
				while ((Dbitmask & (1 << m)) == false) {
					std::vector<int> indices2; // D
					for (int k = 0; k < m; k++) {
						if (Dbitmask & (1<<k)) {
							indices2.push_back(k);
						}
					}
					dp[Dbitmask] = INT_MAX;

					int q = indices2.size();
					uint upper_bm = (1 << q)-1;
					for (uint subX = 1; subX < upper_bm; subX++) { // X bitmask
						uint lambdaX = 0; // \Lambda(X)
						uint lambdaDmX = 0; // \Lambda(D\X)
						uint Xbitmask = 0, DmXbitmask = 0;
						for (int k = 0; k < q; k++) {
							if (subX & (1<<k)) {
								lambdaX |= components[Lbitmask][indices2[k]];
								Xbitmask |= (1 << indices2[k]);

							} else {
								lambdaDmX |= components[Lbitmask][indices2[k]];
								DmXbitmask |= (1 << indices2[k]);
							}
						}

						int optX = opt[lambdaX];
						int optDmX = opt[lambdaDmX];
						if (!minrs) {
							std::size_t countX = std::bitset<32>(lambdaX).count();
							std::size_t countDmX = std::bitset<32>(lambdaDmX).count();
							optX += comb2(countX) * (i-countX);
							optDmX += comb2(countDmX) * (i-countDmX);
						}

						int min2 = std::min(dp[DmXbitmask], optDmX);
						if (dp[Dbitmask] > optX + min2) {
							dp[Dbitmask] = optX + min2;
							if (min2 == optDmX) {
								dp_backtrack[Lbitmask][Dbitmask] = -(int) Xbitmask;
							} else {
								dp_backtrack[Lbitmask][Dbitmask] = (int) Xbitmask;
							}
							// it is important that whenever DP and opt are the same, we store DmXbitmask as negative
							// (see print_tree)
						}
					}

					uint t = (Dbitmask | (Dbitmask - 1)) + 1;
					Dbitmask = t | ((((t & -t) / (Dbitmask & -Dbitmask)) >> 1) - 1);
				}
			}

			opt[Lbitmask] = dp[(1 << m)-1] + minrs;

			// http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
			uint t = (Lbitmask | (Lbitmask - 1)) + 1;
			Lbitmask = t | ((((t & -t) / (Lbitmask & -Lbitmask)) >> 1) - 1);
		}
	}


	uint Lbitmask = states-1;
	Tree* tree = new Tree();
	print_tree(Lbitmask, components, dp_backtrack, tree, tree->add_node(), minrs);

	delete[] opt;
	delete[] dp;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			delete[] triplets[i][j];
		}
		delete[] triplets[i];
	}
	delete[] triplets;
	for (int i = 0; i < (int) states; i++) {
		delete[] dp_backtrack[i];
	}
	delete[] dp_backtrack;

	delete[] components;

	return tree;
}


Tree* minRLC_exact(std::vector<Tree*>& trees) {
	return minRILC(trees, true);
}
Tree* minILC_exact(std::vector<Tree*>& trees) {
	return minRILC(trees, false);
}

void ahoRec(std::vector<int> indices, int*** triplets, int n, Tree* tree, Tree::Node* node) {
	std::vector<uint> components = build_aho(indices, triplets, false);
	for (uint comp : components) {
		std::vector<int> indices2;
		for (int j = 0; j < n; j++) {
			if (comp & (1<<j)) {
				indices2.push_back(j);
			}
		}

		if (indices2.size() == 1) {
			node->add_child(tree->add_node(indices2[0]));
		} else {
			Tree::Node* new_node = tree->add_node();
			node->add_child(new_node);
			ahoRec(indices2, triplets, n, tree, new_node);
		}
	}
}

Tree* ahoBuild(std::vector<Tree*>& trees) {
	int n = trees[0]->get_taxas_num();
	int*** triplets = get_common_triplets(trees);

	std::vector<int> indices;
	for (int i = 0; i < n; i++) {
		indices.push_back(i);
	}

	Tree* tree = new Tree();
	ahoRec(indices, triplets, n, tree, tree->add_node());
	return tree;
}


#endif /* SRC_LOCAL_CONSENSUS_H_ */
