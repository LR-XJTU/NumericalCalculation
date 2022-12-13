#ifndef FORMULA_H
#define FORMULA_H

#include <fstream>
#include <queue>
#include <vector>
#include <stack>

enum op { FNUM, FX, FADD, FSUB, FMUL, FDIV, FPOW, FFAC, FSQRT, FEXP, FLN, FSIN, FCOS, FTAN, FASIN, FACOS, FATAN, FABS, FBRAL, FBRAR };

class formula {
private:
	std::string fstr;
	std::vector<op> frpn;
	std::queue<double> fnum;
	std::vector<double> dfx;
	std::vector<double> dfy;
public:
	void init();
	void init(std::string);
	void init(std::string fx, std::string xt);
	void showtips();
	void showtips_xnum();
	bool x_flag();
	void define_xy();
	bool find_xy(double);
	double list_xy(double);
	void set_fstr(std::string);
	std::string get_fstr();
	std::string replace_x(std::string fx, std::string xt);
	std::vector<op> get_frpn();
	std::queue<double> get_fnum();
	std::vector<double> get_dfx();
	std::vector<double> get_dfy();
	void rf_str(std::string& str);
	void trans_rpn(std::string str);
	void matlab_format();
	void check_brackets(std::string& str, size_t anchor, int length);
	double f(double x);
	double f_xnum(std::vector<double> x);
	void check_stack(std::stack<double> temp, int num);
	void skipbracket(std::string str, int &pos);
	formula operator=(formula);
	formula operator*(formula);
};

class formulae {
private:
	int n;
	std::vector<formula> fx;
	std::vector<std::vector<int>> xnum;
public:
	void init(int n);
	void init(int, const char*);
	formula getformula(int);
	std::vector<int> get_xnum(int);
	std::string get_fstr(int);
	void recog_xnum();
	std::vector<double> f(std::vector<double> x);
};

#endif