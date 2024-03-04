#include "KalmanFilter.h"

int KalmanABC::initializate(
		const double* x, const double* A,
		const double* P, const double* Q, int n,
		const double* B, const double* u, int l) {
	int n2 = n * n;
	size_t N  = sizeof(double) * n;
	size_t N2 = sizeof(double) * n2;
	try {
		_x = new double[n];  memcpy(_x, x, N);
		_A = new double[n2]; memcpy(_A, A, N2);
		_P = new double[n2]; memcpy(_P, P, N2);
		_Q = new double[n2]; memcpy(_Q, Q, N2);
	}
	catch (std::bad_alloc) { return 1; }
	_n = n;

	if (l == 0) return 0;
	if (B == NULL || u == NULL) return 1;

	int ln = l * n;
	size_t L  = sizeof(double) * l;
	size_t LN = sizeof(double) * ln;

	try {
		_u = new double[l];  memcpy(_u, u, L);
		_B = new double[ln]; memcpy(_B, B, LN);
	}
	catch (std::bad_alloc) { return 1; }
	_l = l;

	return 0;
}
