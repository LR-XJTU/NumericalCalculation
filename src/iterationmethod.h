#ifndef ITERATIONMETHOD_H
#define ITERATIONMETHOD_H

#include "linear_equations.h"
#include "common_fd.h"
#include "filelog.h"

class Iteration_method :protected Ax_b {
protected:
	std::vector<double> itrerr;
	double eps;
	int itcounter;
	int maxcounter;
	filelog fl;
public:
	void init();
	void save_Ab_sym_tri();
	void out_Ab();
	void out_x();
	void set_eps();
	void set_max();
	void resize_itrerr();
	virtual void calc() = 0;
	virtual void out_result();
	virtual ~Iteration_method() = default;
	void generate_m();
	void exchange_diag_no_0(int);
};


class Jacobi :public Iteration_method {
public:
	Jacobi();
	void calc();
	void out_result();
};

class Gauss_Seidel :public Iteration_method {
public:
	Gauss_Seidel();
	void calc();
	void out_result();
};

class SOR :public Iteration_method {
private:
	double omega;
public:
	SOR();
	void in_omega();
	void calc();
	void out_result();
};

class steepest_descent :public Iteration_method {
public:
	steepest_descent();
	void calc();
	void out_result();
};

class conjugate_gradient :public Iteration_method {
private:
	bool print_state;
public:
	conjugate_gradient();
	conjugate_gradient(const std::vector<std::vector<double>>& AA, const std::vector<double>& bb);
	void calc();
	void out_result();
	void get_result(std::vector<double>& xx);
};

class GMRES :public Iteration_method {
private:
	int m;
public:
	GMRES();
	void in_m();
	void calc();
	void out_result();
};

#endif
