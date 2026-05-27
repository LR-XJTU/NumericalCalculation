#ifndef COMMON_FD_H
#define COMMON_FD_H

#define sqr(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

#define CHINESE_VERSION

#include <sstream>
#include <vector>

double sign(double x);
int strtoint(std::string s);
int in_int();
int calc_fac(int k);
double calc_fac(double k);
double vecnorm1(double* xx, int nn);
double vecnorm2(double* xx, int nn);
double vecnorm2(std::vector<double> xx);
double vecnorminf(double* xx, int nn);
double matnorm1(double** mm, int nn);
double matnorminf(double** mm, int nn);
std::vector<std::vector<double>> mat_inverse(std::vector<std::vector<double>>);
std::vector<std::vector<double>> mat_multi(std::vector<std::vector<double>>, std::vector<std::vector<double>>);
std::vector<double> mat_multi_vec(std::vector<std::vector<double>>, std::vector<double>);


template <class T>
T int_pow(T a, int x) {
	T r;
	r = (T)1;
	for (int i = 0; i < x; i++) r *= a;
	return r;
}

template <class Type>
Type StringToNum(const std::string& str)
{
	std::istringstream iss(str);
	Type num;
	iss >> num;
	return num;
}

//std::string doubletostr(double d);

/*
class mat {
private:
	char name[10];
	int m;
	int n;
	double **data;
public:
	mat();
	mat(const char *name);
	mat(const char *name, int m, int n, double** data);
	void input_mat();
	void set_m();
	void set_n();
	void set_data();
	int get_m();
	int get_n();
	double get_data(int i, int j);
	double matnorm1(double** mm, int nn);
	double matnorminf(double** mm, int nn);
};
*/

#endif
