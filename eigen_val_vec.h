#ifndef EIGEN_VAL_VEC
#define EIGEN_VAL_VEC

#include <vector>
#include "filelog.h"

class Eigen_val_vec {
protected:
	int n;
	std::vector<std::vector<double>> A;
	double m;
	std::vector<double> z;
	double eps1;
	double eps2;
	int counter;
	filelog fl;
public:
	void init();
	virtual void calc() = 0;
	virtual void out_result() = 0;
};

class Power_method :public Eigen_val_vec {
private:
	double p;
public:
	Power_method();
	void calc();
	void out_result();
};

class Inverse_power :public Eigen_val_vec {
private:
	double p;
	double lambda;
public:
	Inverse_power();
	void calc();
	void out_result();
};

#endif
