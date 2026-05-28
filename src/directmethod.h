#ifndef DIRECTMETHOD_H
#define DIRECTMETHOD_H

#include <vector>
#include "linear_equations.h"
#include "common_fd.h"
#include "filelog.h"

class Direct_method:public Ax_b {
protected:
	filelog fl;
public:
	void init();
	void init(int mm, int nn, const std::vector<std::vector<double>>& AA, const std::vector<double>& bb);
	void init(int mm, int nn, const std::vector<double>& dd, const std::vector<double>& ud, const std::vector<double>& ld, const std::vector<double>& bb);
	void out_Ab();
	void out_x();
	virtual void calc() = 0;
	virtual void out_result() = 0;
	virtual ~Direct_method() = default;
};

class Gauss :public Direct_method {
public:
	Gauss();
	Gauss(bool no_init);
	void calc();
	void out_result();
};

class CP_Gauss :public Direct_method {
public:
	CP_Gauss();
	CP_Gauss(bool no_init);
	void calc();
	void out_result();
};

class Doolittle :public Direct_method {
private:
	bool print_flag;
public:
	Doolittle();
	Doolittle(bool no_init);
	void decomposition(std::vector<std::vector<double>> AA, std::vector<std::vector<double>>& L, std::vector<std::vector<double>>& R);
	void L_solve(std::vector<std::vector<double>> LL, std::vector<double>& yy, std::vector<double> bb);
	void U_solve(std::vector<std::vector<double>> UU, std::vector<double>& xx, std::vector<double> yy);
	void out_LU();
	void calc();
	void out_result();
};

class Cholesky :public Direct_method {
public:
	Cholesky();
	Cholesky(bool no_init);
	void out_G();
	void calc();
	void out_result();
};

class Improved_sqrt :public Direct_method {
public:
	Improved_sqrt();
	Improved_sqrt(bool no_init);
	void out_LD();
	void calc();
	void out_result();
};

class Chasing :public Direct_method {
public:
	Chasing();
	Chasing(bool no_init);
	void save_Ab_sym_tri();
	void out_LU();
	void calc();
	void out_result();
};

class Givens :public Direct_method {
private:
	std::vector<std::vector<double>> Q;
public:
	Givens();
	Givens(bool no_init);
	void out_QR();
	void calc();
	void out_result();
};

class Householder :public Direct_method {
private:
	std::vector<std::vector<double>> Q;
	std::vector<double> d;
	std::vector<double> alpha;
	bool print_flag;
public:
	Householder();
	Householder(bool no_init);
	Householder(int mm, int nn, std::vector<std::vector<double>> AA, std::vector<double> bb);
	void out_QR();
	void calc();
	std::vector<double> get_x();
	void out_result();
};

#endif
