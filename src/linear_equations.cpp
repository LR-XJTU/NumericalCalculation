#include <iostream>
#include "linear_equations.h"
#include "common_fd.h"
using namespace std;

void Ax_b::delete_Ax_b() {
	if (init_flag) {
		for (int i = 0; i < m; i++) delete[] A[i];
		delete[] A;
		delete[] b;
		delete[] x;
	}
}

void Ax_b::enter_Ab() {
#ifdef CHINESE_VERSION
	cout << "\n请输入 A (" << m << "*" << n << ")：" << endl;
#else
	cout << "\nPlease input A (" << m << "*" << n << "):" << endl;
#endif
	double data;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cin >> data;
			in_A(i, j, data);
		}
	}
#ifdef CHINESE_VERSION
	cout << "\n请输入 b (" << m << "*" << 1 << "):" << endl;
#else
	cout << "\nPlease input b (" << m << "*" << 1 << "):" << endl;
#endif
	for (int i = 0; i < m; i++)
	{
		cin >> data;
		in_b(i, data);
	}
#ifdef CHINESE_VERSION
	cout << "\n输入完成" << endl;
#else
	cout << "\nFinish inputting." << endl;
#endif
}

void Ax_b::init_x() {
#ifdef CHINESE_VERSION
	cout << "\n初始化 x ：" << endl;
#else
	cout << "\nPlease initialize x:" << endl;
#endif
	for (int i = 0; i < n; i++) cin >> x[i];
}

void Ax_b::init_x(double x0) {
	for (int i = 0; i < n; i++) x[i] = x0;
}

void Ax_b::A_init(int mn) {
	m = mn;
	n = mn;
	A = new double* [n];
	for (int i = 0; i < n; i++) A[i] = new double[n];
	b = new double[n];
	x = new double[n];
	init_flag = true;
}

void Ax_b::A_init(int mm, int nn) {
	m = mm;
	n = nn;
	A = new double* [m];
	for (int i = 0; i < m; i++) A[i] = new double[n];
	b = new double[m];
	x = new double[n];
	init_flag = true;
}

void Ax_b::in_A(int i, int j, double data) {
	A[i][j] = data;
}

void Ax_b::in_b(int i, double data) {
	b[i] = data;
}

void Ax_b::copy_A(double** AA) {
	for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) A[i][j] = AA[i][j];
}

void Ax_b::copy_b(double* bb) {
	for (int i = 0; i < m; i++) b[i] = bb[i];
}

void Ax_b::copy_A(vector<vector<double>> AA) {
	for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) A[i][j] = AA[i][j];
}

void Ax_b::copy_b(vector<double> bb) {
	for (int i = 0; i < m; i++) b[i] = bb[i];
}

void Ax_b::get_x(double* r) {
	for (int i = 0; i < n; i++) r[i] = x[i];
}

void Ax_b::exchange_row(int i1, int i2)
{
	double temp;
	for (int j = 0; j < n; j++)
	{
		temp = A[i1][j];
		A[i1][j] = A[i2][j];
		A[i2][j] = temp;
	}
	temp = b[i1];
	b[i1] = b[i2];
	b[i2] = temp;
}

bool Ax_b::check_symmetry()
{
	if (m == n)
	{
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < i; j++)
			{
				if (A[i][j] != A[j][i]) return false;
			}
		}
		return true;
	}
	return false;
}

bool Ax_b::check_tridiagonal()
{
	if (m == n)
	{
		for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) if (abs(i - j) > 1 && A[i][j] != 0.0) return false;
		return true;
	}
	return false;
}

void Ax_b::construct_tri(double* dd, double* ud, double* ld) {
	for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) {
		if (i == j) A[i][j] = dd[i];
		else if (i - j == 1) A[i][j] = ld[j];
		else if (j - i == 1) A[i][j] = ud[i];
		else A[i][j] = 0.0;
	}
}

void Ax_b::construct_sym_tri() {
	double x1, x2;
#ifdef CHINESE_VERSION
	cout << "\n构造对称三角矩阵 A" << endl;
	cout << "例如：n = 4\nA =" << endl;
#else
	cout << "\nConstructing a symmetric tridiagonal matrix A." << endl;
	cout << "e.g. n = 4\nA =" << endl;
#endif
	cout << "\td\ts\t0\t0" << endl;
	cout << "\ts\td\ts\t0" << endl;
	cout << "\t0\ts\td\ts" << endl;
	cout << "\t0\t0\ts\td" << endl;
#ifdef CHINESE_VERSION
	cout << "\n请输入对角元素：\nd = ";
#else
	cout << "\nPlease enter diagonal element:\nd = ";
#endif
	cin >> x1;
#ifdef CHINESE_VERSION
	cout << "\n请输入次对角元素：\ns = ";
#else
	cout << "\nPlease enter subdiagonal element:\ns = ";
#endif
	cin >> x2;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) A[i][j] = x1;
			else if (i - j == 1 || j - i == 1) A[i][j] = x2;
			else A[i][j] = 0.0;
		}
	}
#ifdef CHINESE_VERSION
	cout << "\n构造 b" << endl;
	cout << "例如：n = 4\nb =" << endl;
#else
	cout << "\nConstructing b." << endl;
	cout << "e.g. n = 4\nb =" << endl;
#endif
	cout << "\ta\n\tc\n\tc\n\ta" << endl;
#ifdef CHINESE_VERSION
	cout << "\n请输入a：\na = ";
#else
	cout << "\nPlease enter a:\na = ";
#endif
	cin >> x1;
#ifdef CHINESE_VERSION
	cout << "\n请输入c：\nc = ";
#else
	cout << "\nPlease enter c:\nc = ";
#endif
	cin >> x2;
	b[0] = x1;
	b[n - 1] = x1;
	for (int i = 1; i < n - 1; i++) b[i] = x2;
}