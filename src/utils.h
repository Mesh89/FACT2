#ifndef UTILS_H_
#define UTILS_H_

#include <sstream>
#include <stack>

inline int** alloc_int_matrix(int n, int v) {
	int** m = new int*[n];
    m[0] = new int[n*n];
    for (int i = 1; i < n; i++)
        m[i] = m[0] + i*n;
    std::fill(m[0], m[0]+n*n, v);
    return m;
}
inline int** alloc_int_matrix(int n) {
	return alloc_int_matrix(n, 0);
}

template<class T> std::string vector_to_string(std::vector<T> v) {
	std::stringstream ss;
	for (unsigned int i = 0; i < v.size(); i++) {
		ss << v[i] << " ";
	}
	return ss.str();
}

inline std::string intarray_to_string(int* a, int n) {
	return vector_to_string(std::vector<int>(a, a+n));
}
inline std::string boolarray_to_string(bool* a, int n) {
	return vector_to_string(std::vector<bool>(a, a+n));
}

inline int int_log2(int n) {
	int targetlevel = 0;
	while (n >>= 1) ++targetlevel;
	return targetlevel;
}

inline int comb2(int n) {
	return n*(n-1)/2;
}

#endif /* UTILS_H_ */
