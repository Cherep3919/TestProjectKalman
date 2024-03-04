#ifndef __KALMAN_FILTER__
#define __KALMAN_FILTER__

#include <string.h>
#include <stdexcept>
#include "linalg.h"

enum kf_type {
	KF_none,
	KF_classical,
	KF_extended, // EKF
};


// Abstract Kalman filter
class KalmanABC {
	// n - count of state parameters
	// m - count of measurements

	// x - state vector
	// A - state transition matrix
	// P - covariance of state vector
	// Q - noise covariance matrix
protected:
	double* _x, * _A, * _P, * _Q;
	double* _B, * _u;

	int _n, _m, _l;

	virtual void predict() = 0;
	virtual int  correct(const double* z, const double* H, const double* R, int m) = 0;

	// A = A + E, E - ones matrix
	virtual void inc_matrix(double* A) {
		for (int i = 0; i < _n; ++i)
			for (int j = 0; j < _n; ++j)
				A[j + i*_n] += j == i ? 1 : 0;
	};

public:
	int calculate(double* x, double* P, const double* z, const double* H, const double* R, int m) {
		if (_n <= 0) return 1;
		predict();
		if (correct(z, H, R, m)) return 1;

		size_t N = sizeof(double) * _n;
		memcpy(x, _x, N);
		memcpy(P, _P, N * _n);
		return 0;
	}


	virtual int initializate(
		const double* x, const double* A,
		const double* P, const double* Q, int n,
		const double* B = NULL, const double* u = NULL, int l = 0);

	virtual ~KalmanABC() {
		if (_n > 0) {
			delete []_x;
			delete []_A;
			delete []_P;
			delete []_Q;
		}
		if (_l > 0) {
			delete []_u;
			delete []_B;
		}
	};
};

class classical_KF : public KalmanABC {
protected:
	double* _x_f;		// forecast state vector
	double* _P_f;		// forecast error matrix
	double* _K;         // size NxM

	// intermediate matrix data
	double* _AP_temp;   // size NxN 
	double* _HP_temp;   // size MxN
	double* _R_temp;    // size MxM
	double* _KH_temp;   // size NxN
	double* _z_temp;    // size Mx1 (vector)

	int realoc(int m) {
		if (m < _m) return 0;;
		size_t nm = m * _n;
		delete []_KH_temp;
		delete []_HP_temp;
		delete []_R_temp;
		delete []_z_temp;
		delete []_K;
		try {
			_HP_temp = new double[nm];
			_KH_temp = new double[nm];
			_R_temp  = new double[m*m];
			_z_temp  = new double[m];
			_K       = new double[nm];
		}
		catch (std::bad_alloc) { return 1; }
		_m = m;
		return 0;
	};

	void predict() {
		// x_f = A*x
		matmul("NN", _n, 1, _n, 1.0, _A, _x, 0.0, _x_f);
		if (_B != NULL)
			matmul("NN", _n, 1, _n, 1.0, _B, _u, 1.0, _x_f); // x_f += B*u

		// P_f = A*P*A.T + Q
		matmul("NN", _n, _n, _n, 1.0, _A, _P, 0.0, _AP_temp);
		matmul("NT", _n, _n, _n, 1.0, _AP_temp, _A, 1.0, _P_f);
	};

	int correct(const double* z, const double* H, const double* R, int m) {
		if (realoc(m)) return 1;

		{	// K = P*H.T*((H*P*H.T + R)^-1)

			matmul("NN", m, _n, _n, 1.0, H, _P_f, 0.0, _HP_temp);		// HP_temp = H*P_f 
			matmul("NT", m, _n, m, 1.0, _HP_temp, H, 1.0, _R_temp);	// R_temp  = H*P_f*H.T + R

			if (inverse(_R_temp, m)) return 1;

			matmul("NT", _n, _n, m, 1.0, _P_f, H, 0.0, _HP_temp);		// HP_temp = P_f*H.T
			matmul("NN", _n, _n, _n, 1.0, _HP_temp, _R_temp, 0.0, _K); // K = HP_temp*R_temp
		}

		{	// P = (I - K*H)*P
			matmul("NN", _n, _n, m, -1.0, _K, H, 0.0, _KH_temp);
			inc_matrix(_KH_temp);
			matmul("NN", _n, _n, _n, 1.0, _KH_temp, _P_f, 0.0, _P);
		}

		{	// x = x_f + K*(z - H*x_f)
			int i;
			matmul("NN", m, _n, 1, -1.0, H, _x_f, 0.0, _z_temp);
			for (i = 0; i < m; i++) _z_temp[i] += z[i];

			matmul("NN", _n, m, 1, 1.0, _K, _z_temp, 0.0, _x);
			for (i = 0; i < _n; i++) _x[i] += _x_f[i];
		}
	};

public:
	virtual int initializate(
		const double* x, const double* A,
		const double* P, const double* Q, int n,
		const double* B = NULL, const double* u = NULL, int l = 0) {
		if (KalmanABC::initializate(x, A, P, Q, n, B, u, l)) return 1;
		try {
			_x_f = new double[_n];
			_P_f = new double[_n*_n];
			_AP_temp = new double[_n * _n];
		}
		catch (std::bad_alloc) { return 1; }
		return 0;
	}

	~classical_KF() {
		if (_n) {
			delete []_x_f;
			delete []_P_f;
			delete[]_AP_temp;
		}
		if (_m > 0) {
			delete []_K;     
			delete []_HP_temp;
			delete []_R_temp;
			delete []_KH_temp;
			delete []_z_temp;
		}
	}
};


class KalmanFilter {
private:
	KalmanABC* k;
	kf_type type;

	KalmanABC* set_kalman_type() {
		switch (type) {
			case KF_none:      return NULL;
			case KF_classical: return new classical_KF;
			case KF_extended:  return NULL;
		}
	};

public:
	int init(const double* x, const double* A,
		const double* P, const double* Q, int n,
		const double* B = NULL, const double* u = NULL, int l = 0) {
		k = set_kalman_type();
		return k->initializate(x, A, P, Q, n, B, u, l);
	};

	int calculate(double* x, double* P, const double* z, const double* H, const double* R, int m) {
		return k->calculate(x, P, z, H, R, m);
	}
	
	~KalmanFilter() {
		delete k;
	}

};

#endif
