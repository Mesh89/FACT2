/*
 * minrs.h
 *
 *  Created on: Feb 14, 2017
 *      Author: ramesh
 */

#ifndef SRC_MINRS_H_
#define SRC_MINRS_H_

#include <vector>
#include <iostream>
#include <bitset>

#include "Tree.h"

typedef unsigned int uint;


uint dfs(std::vector<int>* adjl, int v, bool* visited) {
	uint comp = (1 << v);
	visited[v] = true;
	for (int i = 0; i < (int) adjl[v].size(); i++) {
		if (!visited[adjl[v][i]]) {
			comp |= dfs(adjl, adjl[v][i], visited);
		}
	}
	return comp;
}

std::vector<uint> build_aho(std::vector<int>& indices, bool*** triplets) {
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
				}
			}
		}
	}

	std::vector<uint> components;
	bool* visited = new bool[n];
	for (int i = 0; i < n; i++) {
		if (!visited[i]) {
			uint comp = dfs(adjl, i, visited);
			if (std::bitset<32>(comp).count() > 1) {
				components.push_back(comp);
			}
		}
	}

	return components;
}


Tree* minRS(std::vector<Tree*>& trees) {
	int n = trees[0]->get_taxas_num();
	bool*** triplets = new bool**[n];
	for (int i = 0; i < n; i++) {
		triplets[i] = new bool*[n];
		for (int j = 0; j < n; j++) {
			triplets[i][j] = new bool[n];
			std::fill(triplets[i][j], triplets[i][j]+n, false);
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
						triplets[i][j][k] = true;
					} else if (dlcaij == dlcajk && dlcaik > dlcaij) { //ik|j
						triplets[i][k][j] = true;
					} else if (dlcaij == dlcaik && dlcajk > dlcaij) { //jk|i
						triplets[j][k][i] = true;
					} else {
						// should never enter here
						std::cout << dlcaij << " " << dlcaik << " " << dlcajk << std::endl;
						assert(false);
					}
				}
			}
		}
	}

	uint states = (1 << n);
	int* opt = new int[states];
	for (int i = 0; i < n; i++) { // init opt for single taxa
		opt[1 << i] = 0;
	}

	int* dp = new int[states];
	dp[0] = 0;

	for (int i = 2; i <= n; i++) {
		uint Lbitmask = (1 << i)-1; // L' bitmask
		while ((Lbitmask & (1 << i)) == false) {
			// indices contains the subset of taxa we are currently dealing with
			std::vector<int> indices;
			for (int j = 0; j < n; j++) {
				if (Lbitmask & (1<<j)) {
					indices.push_back(j);
				}
			}

			// aho graph components
			std::vector<uint> components = build_aho(indices, triplets);
			int m = components.size();

			// init dp for single components
			for (int j = 0; j < m; j++) { // init opt for single taxa
				dp[1 << j] = opt[components[j]];
			}

			for (int j = 2; j <= m; j++) {
				uint Dbitmask = (1 << j)-1; // D bitmask
				while ((Dbitmask & (1 << j)) == false) {
					dp[Dbitmask] = INT_MAX;

					std::vector<int> indices2; // D
					for (int k = 0; k < m; k++) {
						if (Dbitmask & (1<<k)) {
							indices2.push_back(k);
						}
					}

					int q = indices2.size();
					uint upper_bm = (1 << q)-1;
					for (uint Xbitmask = 1; Xbitmask < upper_bm; Xbitmask++) { // X bitmask
						uint lambdaX = 0; // \Lambda(X)
						uint lambdaDmX = 0; // \Lambda(D\X)
						for (int k = 0; k < q; k++) {
							if (Xbitmask & (1<<k)) {
								lambdaX |= components[indices2[k]];
							} else {
								lambdaDmX |= components[indices2[k]];
							}
						}

						uint DmXbitmask = ~Xbitmask & upper_bm;
						int min2 = std::min(dp[DmXbitmask], opt[lambdaDmX]);
						dp[Dbitmask] = std::min(dp[Dbitmask], opt[lambdaX]+min2);
					}

					uint t = (Dbitmask | (Dbitmask - 1)) + 1;
					Dbitmask = t | ((((t & -t) / (Dbitmask & -Dbitmask)) >> 1) - 1);
				}
			}

			opt[Lbitmask] = dp[(1 << m)-1] + 1;

			// http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
			uint t = (Lbitmask | (Lbitmask - 1)) + 1;
			Lbitmask = t | ((((t & -t) / (Lbitmask & -Lbitmask)) >> 1) - 1);
		}
	}

	std::cout << opt[states-1] << std::endl;

	return NULL;
}


#endif /* SRC_MINRS_H_ */
