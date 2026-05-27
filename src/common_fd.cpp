#include <iostream>
#include <cmath>
#include <string>
#include "common_fd.h"

using namespace std;

//sign function
double sign(double x) {
	if(x>0.0) return 1.0;
	else if(x==0.0) return 0.0;
	else return -1.0;
}

//transform string into integer
int strtoint(string s) {
	int r=0;
	int temp=0;
	for(int i=(((int)s.length())-1);i>0;i--) {
		r += (s[i]-'0')*int_pow(10,temp);
		temp++;
	}
	if(s[0]=='-') r *= -1;
	else r += (s[0]-'0')* int_pow(10,temp);
	return r;
}

//input an integer
int in_int() {
	string s;
	cin >> s;
	for (int i = 0; i < s.length(); i++)
	{
		if (s[i] - '0' < 0 || s[i] - '9' > 0)
		{
#ifdef CHINESE_VERSION
			cout << "错误：输入的不是整数" << endl;
#else
			cout << "Error: The input is not a integer." << endl;
#endif
			throw 0;
		}
	}
	int intnum = strtoint(s);
	return intnum;
}

//calculate the factorial
int calc_fac(int k) {
	int r = 1;
	for (int i = k; i > 0; i--) r *= i;
	return r;
}

double calc_fac(double k) {
	int x = (int)k;
	if ((k - x) != 0.0) {
#ifdef CHINESE_VERSION
		cout << "错误：输入的不是整数，无法阶乘" << endl;
#else
		cout << "Error: The input is not a integer for factorial." << endl;
#endif
		throw 0;
	}
	double r = 1.0;
	for (int i = x; i > 0; i--) r *= i;
	return r;
}

//1-norm of vector
double vecnorm1(double* xx, int nn) {
	double norm = 0.0;
	for (int i = 0; i < nn; i++) norm += fabs(xx[i]);
	return norm;
}

//2-norm of vector
double vecnorm2(double* xx, int nn) {
	double norm = 0.0;
	for (int i = 0; i < nn; i++) norm += sqr(xx[i]);
	norm = sqrt(norm);
	return norm;
}

//2-norm of vector
double vecnorm2(vector<double> xx) {
	double norm = 0.0;
	for (int i = 0; i < xx.size(); i++) norm += sqr(xx[i]);
	norm = sqrt(norm);
	return norm;
}

//inf-norm of vector
double vecnorminf(double* xx, int nn) {
	double norm = fabs(xx[0]);
	for (int i = 1; i < nn; i++) if (fabs(xx[i]) > norm) norm = fabs(xx[i]);
	return norm;
}

//1-norm of matrix
double matnorm1(double** mm, int nn) {
	double norm = 0.0;
	for (int i = 0; i < nn; i++) norm += fabs(mm[i][0]);
	for (int j = 1; j < nn; j++) {
		double temp = 0.0;
		for (int i = 0; i < nn; i++) temp += fabs(mm[i][j]);
		if (temp > norm) norm = temp;
	}
	return norm;
}

//inf-norm of matrix
double matnorminf(double** mm, int nn) {
	double norm = 0.0;
	for (int i = 0; i < nn; i++) norm += fabs(mm[0][i]);
	for (int i = 1; i < nn; i++) {
		double temp = 0.0;
		for (int j = 0; j < nn; j++) temp += fabs(mm[i][j]);
		if (temp > norm) norm = temp;
	}
	return norm;
}

vector<vector<double>> mat_inverse(vector<vector<double>> A0) {
	vector<vector<double>> A;
	int n = A0.size();
	for (int i = 0; i < n; i++) {
		if (A0[i].size() != n) {
#ifdef CHINESE_VERSION
			cout << "错误：矩阵不是方阵，无法求逆" << endl;
#else
			cout << "Error: The matrix calculated for inverse matrix is not square." << endl;
#endif
			throw 0;
		}
	}
	A.resize(n);
	for (int i = 0; i < n; i++) {
		A[i].resize(n);
		A[i][i] = 1.0;
	}
	double t;
	for (int i = 0; i < n; i++) {
		for (int ii = i + 1; ii < n; ii++) {
			t = A0[ii][i] / A0[i][i];
			for (int j = i; j < n; j++) A0[ii][j] -= t * A0[i][j];
			for (int j = 0; j < n; j++) A[ii][j] -= t * A[i][j];
		}
		t = 1.0 / A0[i][i];
		for (int j = i; j < n; j++) A0[i][j] *= t;
		for (int j = 0; j < n; j++) A[i][j] *= t;
		for (int ii = i - 1; ii >= 0; ii--) {
			t = A0[ii][i] / A0[i][i];
			for (int j = i; j < n; j++) A0[ii][j] -= t * A0[i][j];
			for (int j = 0; j < n; j++) A[ii][j] -= t * A[i][j];
		}
	}
	return A;
}

vector<vector<double>> mat_multi(vector<vector<double>> A, vector<vector<double>> B) {
	int m1, m2, n1, n2;
	m1 = A.size();
	n1 = A[0].size();
	m2 = B.size();
	n2 = B[0].size();
	for (int i = 0; i < m1; i++) {
		if (A[i].size() != n1) {
#ifdef CHINESE_VERSION
			cout << "错误：矩阵阶数不一致" << endl;
#else
			cout << "Error: The order of matrix is inconsistent" << endl;
#endif
			throw 0;
		}
	}
	for (int i = 0; i < m2; i++) {
		if (B[i].size() != n2) {
#ifdef CHINESE_VERSION
			cout << "错误：矩阵阶数不一致" << endl;
#else
			cout << "Error: The order of matrix is inconsistent" << endl;
#endif
			throw 0;
		}
	}
	if (n1 != m2) {
#ifdef CHINESE_VERSION
		cout << "错误：矩阵阶数不对应，无法相乘" << endl;
#else
		cout << "Error: Wrong order of matrix multiplication." << endl;
#endif
		throw 0;
	}
	vector<vector<double>> C(m1);
	for (int i = 0; i < m1; i++) C[i].resize(n2);
	for (int i = 0; i < m1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n1; k++) C[i][j] += A[i][k] * B[k][j];
	return C;
}

vector<double> mat_multi_vec(vector<vector<double>> A, vector<double> a) {
	int m, n;
	m = A.size();
	n = A[0].size();
	for (int i = 0; i < m; i++) if (A[i].size() != n) {
#ifdef CHINESE_VERSION
		cout << "错误：矩阵阶数不一致" << endl;
#else
		cout << "Error: The order of matrix is inconsistent" << endl;
#endif
		throw 0;
	}
	if (a.size() != n) {
#ifdef CHINESE_VERSION
		cout << "错误：矩阵与向量阶数不对应，无法相乘" << endl;
#else
		cout << "Error: Incosistent order between the matrix and the vector." << endl;
#endif
		throw 0;
	}
	vector<double> b(m);
	for (int i = 0; i < m; i++) {
		b[i] = 0.0;
		for (int j = 0; j < n; j++) b[i] += A[i][j] * a[j];
	}
	return b;
}

/*
string doubletostr(double d) {
	string t;
	if (d == 0.0) {
		t += "0";
		return t;
	}
	if (d < 0.0) {
		t += "-";
		d = fabs(d);
	}
	int eflag = 0;
	if (d >= 999999.5) eflag = 1;
	else if (d < 0.00009999995) eflag = -1;

	int ord = 0;
	while (d >= 10.0) {
		d /= 10.0;
		ord++;
	}
	while (d < 1.0) {
		d *= 10.0;
		ord--;
	}
	int tempint = (int)(d * 1000000.0);
	if (tempint - ((int)(tempint / 10)) * 10 >= 5) tempint += 10;
	tempint /= 10;
	if (tempint == 1000000) {
		tempint /= 10;
		ord++;
	}
	int zeros = 0;
	while (tempint % 10 == 0) {
		tempint /= 10;
		zeros++;
	}
	char *tempchar = new char[6 - zeros];
	for (int i = 6 - zeros - 1; i >= 0; i--) {
		tempchar[i] = tempint - ((int)(tempint / 10)) * 10 + '0';
		tempint /= 10;
	}
	if (tempint != 0) {
		cout << "Fail to convert double to char*" << endl;
		t += "0";
		return t;
	}
	for (int i = 0; i < 6 - zeros; i++) {
		t += tempchar[i];
		if (i == ord) t += '.';
	}
	delete [] tempchar;
	if (eflag > 0) {
		t += "e+";
		t += doubletostr((double)ord);
	}
	else if (eflag < 0) {
		ord = abs(ord);
		t += "e-";
		t += doubletostr((double)ord);
	}
	else {
		for (int i = 0; i < zeros; i++) t += '0';
	}
	return t;
}
*/
