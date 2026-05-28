#ifndef LINEAR_EQUATIONS
#define LINEAR_EQUATIONS

#include <vector>

class Ax_b {
protected:
	int m;
	int n;
	std::vector<std::vector<double>> A;
	std::vector<double> x;
	std::vector<double> b;
public:
	virtual ~Ax_b() = default;
	void A_init(int mn);
	void A_init(int mm, int nn);
	void enter_Ab();
	void init_x();
	void init_x(double x0);
	void in_A(int i, int j, double data);
	void in_b(int i, double data);
	void copy_A(const std::vector<std::vector<double>>& AA);
	void copy_b(const std::vector<double>& bb);
	void get_x(std::vector<double>& xx);
	void exchange_row(int a, int b);
	bool check_symmetry();
	bool check_tridiagonal();
	void construct_tri(const std::vector<double>& d, const std::vector<double>& ud, const std::vector<double>& ld);
	void construct_sym_tri();
};
#endif
