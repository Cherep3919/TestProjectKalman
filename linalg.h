
/* ---- matrix multiplication ----- 
* s="NN" --->  C = a*A   * B   + c*C
* s="NT" --->  C = a*A   * B.T + c*C
* s="TN" --->  C = a*A.T * B   + c*C
* s="TT" --->  C = a*A.T * B.T + c*C
* -----------------------------------*/
void matmul(const char* s, int n, int m, int k,
	double a, const double* A, const double* B,
	double c, double* C);

int inverse(double* A, int n);
