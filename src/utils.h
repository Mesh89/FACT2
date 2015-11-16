#ifndef UTILS_H_
#define UTILS_H_

#include <sstream>

int** alloc_int_matrix(int n, int v) {
	int* a = new int[n*n];
	std::fill(a, a+n*n, v);

	int** m = new int*[n];
	for (int i = 0; i < n; i++) {
		m[i] = a+(i*n);
	}
	return m;
}
int** alloc_int_matrix(int n) {
	return alloc_int_matrix(n, 0);
}

template<class T> std::string vector_to_string(std::vector<T> v) {
	std::stringstream ss;
	for (unsigned int i = 0; i < v.size(); i++) {
		ss << v[i] << " ";
	}
	return ss.str();
}

std::string intarray_to_string(int* a, int n) {
	return vector_to_string(std::vector<int>(a, a+n));
}
std::string boolarray_to_string(bool* a, int n) {
	return vector_to_string(std::vector<bool>(a, a+n));
}

inline int int_log2(int n) {
	int targetlevel = 0;
	while (n >>= 1) ++targetlevel;
	return targetlevel;
}


#endif /* UTILS_H_ */