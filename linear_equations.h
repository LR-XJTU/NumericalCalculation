#ifndef LINEAR_EQUATIONS
#define LINEAR_EQUATIONS

#include <vector>

class Ax_b {
protected:
	int m;
	int n;
	double** A;
	double* x;
	double* b;
	bool init_flag = false;
public:
	void delete_Ax_b();
	void A_init(int mn);
	void A_init(int mm, int nn);
	void enter_Ab();
	void init_x();
	void init_x(double x0);
	void in_A(int i, int j, double data);
	void in_b(int i, double data);
	void copy_A(double** AA);
	void copy_b(double* bb);
	void copy_A(std::vector<std::vector<double>> AA);
	void copy_b(std::vector<double> bb);
	void get_x(double *xx);
	void exchange_row(int a, int b);
	bool check_symmetry();
	bool check_tridiagonal();
	void construct_tri(double* d, double* ud, double* ld);
	void construct_sym_tri();
};
#endif
