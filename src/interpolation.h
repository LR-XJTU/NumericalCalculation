#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "formula.h"
#include "filelog.h"

class Interpolation {
protected:
	int np1;
	bool input_fx;
	double* x;
	double* y;
	formula fx;
	filelog fl;
public:
	void init();
	virtual void input_data();
	virtual void cout_xy();
	virtual void save_xy();
	virtual void sort_x();
	void in_x(int i,double t);
	void in_y(int i,double t);
	void cout_poly(double* p, int l);
	void save_poly(double* p, int l);
	std::string polytostr(double* p, int l);
	virtual void cout_polynomial() = 0;
	virtual void save_polynomial() = 0;
	virtual void calc() = 0;
	virtual void out_result() = 0;
};

class Newton_Ip :public Interpolation {
protected:
	double** f;
	double* N;
public:
	Newton_Ip();
	Newton_Ip(bool no_init);
	Newton_Ip(double* xx, double* yy, double nnp1);
	~Newton_Ip();
	void cout_polynomial();
	void save_polynomial();
	void calc();
	void pure_calc();
	std::string polytostr_nodot(double* p, int l);
	std::string get_str();
	void out_result();
	void generate_m();
	void difference_quotient(int i, int j);
	void pi_expansion(int ord);
	void add_point(double xx, double yy);
};

class Hermite_Ip :public Newton_Ip {
protected:
	int* d;
public:
	Hermite_Ip();
	~Hermite_Ip();
	void input_data();
	void cout_xy();
	void save_xy();
	void sort_x();
	void in_d(int i, int t);
	void cout_polynomial();
	void save_polynomial();
	void calc();
	void out_result();
	void difference_quotient(int i, int j);
	void add_point(double xx, double yy, int dd);
};

class cube_spline :public Interpolation {
private:
	int bdcd_flag;
	double bd[2];
	double** S;
public:
	cube_spline();
	~cube_spline();
	void input_bd();
	void cout_polynomial();
	void save_polynomial();
	void calc();
	void out_result();
	void generate_m();
};

#endif