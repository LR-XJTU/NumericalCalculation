#ifndef NONLINEAR_EQUATIONS
#define NONLINEAR_EQUATIONS

#include <vector>
#include "int_diff.h"
#include "formula.h"
#include "filelog.h"

class nl_eq {
protected:
	formula fx;
	std::vector<double> x;
	double interval[2];
	filelog fl;
	double eps1;
	double eps2;
	int maxcounter;
	bool count_flag;
	int counter;
public:
	nl_eq();
	void input_fx();
	void input_x();
	void input_max_eps();
	void isolation();
	double bisection(double a, double b);
	virtual void calc() = 0;
	virtual void out_result() = 0;
};

class eq_Simple :public nl_eq {
private:
	std::string str;
public:
	eq_Simple();
	void input_fx();
	void input_x();
	void calc();
	void out_result();
};

class eq_Newton :public nl_eq {
public:
	eq_Newton();
	void calc();
	void out_result();
};

class eq_Secant :public nl_eq {
private:
	std::vector<double> x0;
public:
	eq_Secant();
	void input_x();
	void calc();
	void out_result();
};

class eq_Relaxation :public eq_Simple {
public:
	void calc();
	void out_result();
};

class eq_Aitken :public eq_Simple {
public:
	void calc();
	void out_result();
};

class nl_eqs {
protected:
	int n;
	formulae fx;
	std::vector<double> x;
	filelog fl;
	double eps1;
	double eps2;
	int maxcounter;
	int counter;
public:
	nl_eqs();
	void input_fx();
	void input_x();
	void input_max_eps();
	virtual void calc() = 0;
	virtual void out_result() = 0;
};

class eqs_Simple :public nl_eqs {
public:
	eqs_Simple();
	void input_fx();
	void calc();
	void out_result();
};

class eqs_Newton :public nl_eqs {
public:
	eqs_Newton();
	void calc();
	void out_result();
};

class eqs_Secant :public nl_eqs {
private:
	double h;
public:
	eqs_Secant();
	void input_h();
	void calc();
	void out_result();
};

class eqs_Broyden :public nl_eqs {
public:
	eqs_Broyden();
	void calc();
	void out_result();
};


#endif