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
	double* wl;
	formula wc;
	//basis function
	int np1;
	double** c;
	formula* basis;
	//result
	double* coef;
	std::string resultstr;
	double err;
	filelog fl;
public:
	optimal_approx();
	double inner_product(formula f, formula g, formula w, double a, double b);
	double inner_product(formula f, formula g, double* x, int m, double *w);
	double inner_product(formula f, double* y, double* x, int m, double* w);
	double inner_product(double* y1, double* y2, int m, double* w);
	void Cheby_poly(double**, int);
	std::string polytostr(double* p, int l);
	virtual void calc() = 0;
	virtual void out_result() = 0;
};

class sqr_approx :public optimal_approx {
protected:
	bool cont;
	//list function
	int m;
	double* x;
	double* y;
	int basis_flag;
public:
	sqr_approx();
	~sqr_approx();
	void sort_xy();
	double innpro(formula f, formula g);
	void orth_poly(formula*, int);
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
	~uni_approx();
	void calc();
	void out_result();
	void generate_m();
	void Cheby_interp();
	void Trunc_Cheby();
};

#endif