#ifndef INT_DIFF_H
#define INT_DIFF_H

#include <vector>
#include "formula.h"
#include "filelog.h"

class Int_Diff {
protected:
	formula fx;
	double result;
	double eps;
	std::vector<double> err;
	filelog fl;
	bool print_state;
public:
	Int_Diff();
	Int_Diff(formula);
	Int_Diff(bool);
	void set_eps();
	void generate_m();
	double get_result();
	virtual void calc() = 0;
	virtual void out_result() = 0;
};

class Romberg:public Int_Diff {
private:
	double interval[2];
public:
	Romberg();
	Romberg(formula, double, double);
	void init();
	void calc();
	void out_result();
	double T(int n);
	double T_Vari_step(int k, double T2k);
	double S_ext(double Tn[2]);
	double C_ext(double Sn[2]);
	double R_ext(double Cn[2]);
};

class DerivExtra:public Int_Diff {
private:
	double x;
	double h;
	int order;
	bool err_flag;
public:
	DerivExtra();
	DerivExtra(formula);
	DerivExtra(bool);
	void init(formula);
	double deriv(double);
	std::vector<double> dxp(std::vector<double>, std::vector<int>, int, double);
	std::vector<double> dxn(std::vector<double>, std::vector<int>, int, double);
	double deriv(std::vector<double>, std::vector<int>, int);
	double orderderiv(int, double);
	void calc();
	void out_result();
	double get_result();
};

#endif