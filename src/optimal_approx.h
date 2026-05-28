#ifndef OPTIMAL_APPROX_H
#define OPTIMAL_APPROX_H

#include "formula.h"
#include "filelog.h"

class optimal_approx {
protected:
	//continuous function
	formula fx;
	//the interval
	double interval[2];
	//weight function
	bool w_flag;
	std::vector<double> wl;
	formula wc;
	//basis function
	int np1;
	std::vector<std::vector<double>> c;
	std::vector<formula> basis;
	//result
	std::vector<double> coef;
	std::string resultstr;
	double err;
	filelog fl;
public:
	optimal_approx();
	virtual ~optimal_approx() = default;
	double inner_product(formula f, formula g, formula w, double a, double b);
	double inner_product(formula f, formula g, const std::vector<double>& x, const std::vector<double>& w);
	double inner_product(formula f, const std::vector<double>& y, const std::vector<double>& x, const std::vector<double>& w);
	double inner_product(const std::vector<double>& y1, const std::vector<double>& y2, const std::vector<double>& w);
	void Cheby_poly(std::vector<std::vector<double>>&, int);
	std::string polytostr(const std::vector<double>& p);
	virtual void calc() = 0;
	virtual void out_result() = 0;
};

class sqr_approx :public optimal_approx {
protected:
	bool cont;
	//list function
	int m;
	std::vector<double> x;
	std::vector<double> y;
	int basis_flag;
public:
	sqr_approx();
	void sort_xy();
	double innpro(formula f, formula g);
	void orth_poly(std::vector<formula>&, int);
	void calc();
	void out_result();
	void generate_m();
};

class uni_approx :public optimal_approx {
protected:
	int ua_method;
	formula fx_xt_cos;
public:
	uni_approx();
	void calc();
	void out_result();
	void generate_m();
	void Cheby_interp();
	void Trunc_Cheby();
};

#endif
