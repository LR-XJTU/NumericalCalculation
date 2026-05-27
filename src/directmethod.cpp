#include <iostream>
#include <cmath>
#include "directmethod.h"

using namespace std;

//initialization without parameters
void Direct_method::init() {
	fl.init("Direct_method.txt");
	cout.setf(ios::left);
	int mm, nn;
#ifdef CHINESE_VERSION
	cout << "\n请输入矩阵A(m*n)的行数和列数：\nm = ";
#else
	cout << "\nPlease enter the numbers of rows and columns of coefficient matrix A(m*n):\nm = ";
#endif
	mm = in_int();
	if(mm<=0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：请输入正整数" << endl;
#else
		cout << "Error: Please enter a positive integer." << endl;
#endif
		throw 0;
	}
	cout << "n = ";
	nn = in_int();
	if(nn<=0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：请输入正整数" << endl;
#else
		cout << "Error: Please enter a positive integer." << endl;
#endif
		throw 0;
	}
	if(mm<nn)
	{
#ifdef CHINESE_VERSION
		cout << "错误：请输入正整数" << endl;
#else
		cout << "Error: Please enter a positive integer." << endl;
#endif
		throw 0;
	}
	A_init(mm,nn);
}

//initialization with 4 parameters
void Direct_method::init(int mm,int nn,double **AA,double *bb) {
	if(mm<=0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：矩阵A的行数不是正整数" << endl;
#else
		cout << "Error: The number of rows of A is not a positive integer." << endl;
#endif
		throw 1;
	}
	if(nn<=0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：矩阵A的列数不是正整数" << endl;
#else
		cout << "Error: The number of columns of A is not a positive integer." << endl;
#endif
		throw 1;
	}
	if(mm<nn)
	{
#ifdef CHINESE_VERSION
		cout << "错误：方程组有无穷解" << endl;
#else
		cout << "Error: Equations have infinite solutions." << endl;
#endif
		throw 1;
	}
	A_init(mm,nn);
	copy_A(AA);
	copy_b(bb);
}

//initialization with 6 parameters
void Direct_method::init(int mm,int nn,double *dd,double *ud,double *ld,double *bb) {
	if(mm<=0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：矩阵A的行数不是正整数" << endl;
#else
		cout << "Error: The number of rows of A is not a positive integer." << endl;
#endif
		throw 1;
	}
	if(nn<=0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：矩阵A的列数不是正整数" << endl;
#else
		cout << "Error: The number of columns of A is not a positive integer." << endl;
#endif
		throw 1;
	}
	if(mm<nn)
	{
#ifdef CHINESE_VERSION
		cout << "错误：方程组有无穷解" << endl;
#else
		cout << "Error: Equations have infinite solutions." << endl;
#endif
		throw 1;
	}
	A_init(mm,nn);
	construct_tri(dd,ud,ld);
	copy_b(bb);
}

//cout and save matrices A, b
void Direct_method::out_Ab() {
	cout << "\nA (" << m << "*" << n << ") =" << endl;
	for (int i = 0; i < m; i++)
	{
		cout << "\t";
		for (int j = 0; j < n; j++)
		{

			cout.width(15);
			cout << A[i][j];
		}
		cout << endl;
	}
	cout << "\nb (" << m << "*" << 1 << ") =" << endl;
	for (int i = 0; i < m; i++) cout << "\t" << b[i] << endl;
	fl << "\nA =\n\n";
	fl.set_mat_size(m, n) << A;
	fl << "\nb = \n\n";
	fl.set_array_len(m) << b;
}

//cout and save matrix x
void Direct_method::out_x()
{
	cout << "\nx (" << n << "*" << 1 << ") =" << endl;
	for (int i = 0; i < n; i++) cout << "\t" << x[i] << endl;
	fl << "\nx = \n\n";
	fl.set_array_len(n) << x;
}

//1.Gauss elimination

Gauss::Gauss() {
	init();
	enter_Ab();
	out_Ab();
}

Gauss::Gauss(bool no_init) {
	if (no_init) return;
	else Gauss();
}

void Gauss::calc()
{
#ifdef CHINESE_VERSION
	cout << "\n正在使用高斯消去法求解 Ax=b ..." << endl;
#else
	cout << "\nSolving Ax=b by Gauss elimination..." << endl;
#endif
	for (int i = 0; i < n - 1; i++)
	{
		if (A[i][i] == 0.0)
		{
			int temp = 0;
			while (A[i + temp][i] == 0.0)
			{
				temp++;
				if (i + temp == m) break;
			}
			if (i + temp < m) exchange_row(i, i + temp);
			else
			{
#ifdef CHINESE_VERSION
				cout << "错误：方程组有无穷解" << endl;
#else
				cout << "Error: Equations have infinite solutions." << endl;
#endif
				throw 0;
			}
		}
		for (int j = i + 1; j < m; j++)
		{
			A[j][i] /= A[i][i];
			for (int k = i + 1; k < n; k++)
			{
				A[j][k] -= A[j][i] * A[i][k];
			}
			b[j] -= A[j][i] * b[i];
		}
	}
	if (A[n - 1][n - 1] == 0.0)
	{
		if (b[n - 1] == 0.0)
		{
#ifdef CHINESE_VERSION
			cout << "错误：方程组有无穷解" << endl;
#else
			cout << "Error: Equations have infinite solutions." << endl;
#endif
			throw 0;
		}
		else
		{
#ifdef CHINESE_VERSION
			cout << "错误：方程组无解" << endl;
#else
			cout << "Error: Equations have no solution." << endl;
#endif
			throw 0;
		}
	}
	if (m > n)
	{
		double temp = b[n - 1] / A[n - 1][n - 1];
		for (int i = n; i < m; i++) if (b[i] / A[i][n - 1] != temp)
		{
#ifdef CHINESE_VERSION
			cout << "错误：方程组无解" << endl;
#else
			cout << "Error: Equations have no solution." << endl;
#endif
			throw 0;
		}
	}
	x[n - 1] = b[n - 1] / A[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = b[i];
		for (int j = i + 1; j < n; j++)
		{
			x[i] -= A[i][j] * x[j];
		}
		x[i] /= A[i][i];
	}
#ifdef CHINESE_VERSION
	cout << "\n求解完成" << endl;
#else
	cout << "\nFinish solving." << endl;
#endif
}

void Gauss::out_result()
{
#ifdef CHINESE_VERSION
	cout << "\n通过高斯消去法解得" << endl;
	fl << "\n通过高斯消去法解得\n";
#else
	cout << "\nSolved by Gauss elimination." << endl;
	fl << "\nSolved by Gauss elimination.\n";
#endif
	out_x();
}

Gauss::~Gauss() {
	delete_Ax_b();
}

//2.Column principle Gaussian elimination

CP_Gauss::CP_Gauss() {
	init();
	enter_Ab();
	out_Ab();
}

CP_Gauss::CP_Gauss(bool no_init) {
	if (no_init) return;
	else CP_Gauss();
}

void CP_Gauss::calc()
{
#ifdef CHINESE_VERSION
	cout << "\n正在使用列主元高斯消去法求解 Ax=b ..." << endl;
#else
	cout << "\nSolving Ax=b by column principle Gaussian elimination..." << endl;
#endif
	for (int i = 0; i < n - 1; i++)
	{
		double temp_max = fabs(A[i][i]);
		int temp = i;
		for (int j = i + 1; j < m; j++) if (fabs(A[j][i]) > temp_max)
		{
			temp_max = fabs(A[j][i]);
			temp = j;
		}
		if (temp_max == 0.0)
		{
#ifdef CHINESE_VERSION
			cout << "错误：方程组有无穷解" << endl;
#else
			cout << "Error: Equations have infinite solutions." << endl;
#endif
			throw 0;
		}
		else exchange_row(i, temp);
		for (int j = i + 1; j < m; j++)
		{
			A[j][i] /= A[i][i];
			for (int k = i + 1; k < n; k++)
			{
				A[j][k] -= A[j][i] * A[i][k];
			}
			b[j] -= A[j][i] * b[i];
		}
	}
	if (A[n - 1][n - 1] == 0.0)
	{
		if (b[n - 1] == 0.0)
		{
#ifdef CHINESE_VERSION
			cout << "错误：方程组有无穷解" << endl;
#else
			cout << "Error: Equations have infinite solutions." << endl;
#endif
			throw 0;
		}
		else
		{
#ifdef CHINESE_VERSION
			cout << "错误：方程组无解" << endl;
#else
			cout << "Error: Equations have no solution." << endl;
#endif
			throw 0;
		}
	}
	if (m > n)
	{
		double temp = b[n - 1] / A[n - 1][n - 1];
		for (int i = n; i < m; i++) if (b[i] / A[i][n - 1] != temp)
		{
#ifdef CHINESE_VERSION
			cout << "错误：方程组无解" << endl;
#else
			cout << "Error: Equations have no solution." << endl;
#endif
			throw 0;
		}
	}
	x[n - 1] = b[n - 1] / A[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = b[i];
		for (int j = i + 1; j < n; j++)
		{
			x[i] -= A[i][j] * x[j];
		}
		x[i] /= A[i][i];
	}
#ifdef CHINESE_VERSION
	cout << "\n求解完成" << endl;
#else
	cout << "\nFinish solving." << endl;
#endif
}

void CP_Gauss::out_result()
{
#ifdef CHINESE_VERSION
	cout << "\n通过列主元高斯消去法解得" << endl;
	fl << "\n通过列主元高斯消去法解得\n";
#else
	cout << "\nSolved by column principle Gaussian elimination." << endl;
	fl << "\nSolved by column principle Gaussian elimination.\n";
#endif
	out_x();
}

CP_Gauss::~CP_Gauss() {
	delete_Ax_b();
}

//3.Doolittle decomposition (LU)

Doolittle::Doolittle() {
	print_flag = true;
	init();
	enter_Ab();
	out_Ab();
}

Doolittle::Doolittle(bool no_init) {
	if (no_init) {
		print_flag = false;
		return;
	}
	else Doolittle();
}

void Doolittle::calc()
{
#ifdef CHINESE_VERSION
	if (print_flag) cout << "\n正在通过杜利特尔分解求解 Ax=b ..." << endl;
#else
	if (print_flag) cout << "\nSolving Ax=b by Doolittle decomposition..." << endl;
#endif
	if (m != n)
	{
#ifdef CHINESE_VERSION
		cout << "错误：A不是方阵" << endl;
#else
		cout << "Error: A is not a square matrix." << endl;
#endif
		throw 0;
	}
	if (A[0][0] == 0.0)
	{
		int temp = 1;
		while (A[temp][0] == 0.0)
		{
			temp++;
			if (temp == n) break;
		}
		if (temp < n) exchange_row(0, temp);
		else
		{
#ifdef CHINESE_VERSION
			cout << "错误：方程组有无穷解" << endl;
#else
			cout << "Error: Equations have infinite solutions." << endl;
#endif
			throw 0;
		}
	}
	for (int i = 1; i < n; i++) A[i][0] /= A[0][0];
	for (int i = 1; i < n - 1; i++)
	{
		double temp_data = A[i][i];
		int temp = 1;
		for (int k = 0; k < i; k++) A[i][i] -= A[i][k] * A[k][i];
		while (A[i][i] == 0.0)
		{
			A[i][i] = temp_data;
			exchange_row(i, i + temp);
			temp_data = A[i][i];
			for (int k = 0; k < i; k++) A[i][i] -= A[i][k] * A[k][i];
			temp++;
			if (i + temp == n) break;
		}
		if (A[i][i] == 0.0)
		{
#ifdef CHINESE_VERSION
			cout << "错误：方程组有无穷解" << endl;
#else
			cout << "Error: Equations have infinite solutions." << endl;
#endif
			throw 0;
		}
		for (int j = i + 1; j < n; j++)
		{
			for (int k = 0; k < i; k++)
			{
				A[i][j] -= A[i][k] * A[k][j];
				A[j][i] -= A[j][k] * A[k][i];
			}
			A[j][i] /= A[i][i];
		}
		for (int k = 0; k < i; k++) b[i] -= A[i][k] * b[k];
	}
	for (int k = 0; k < n - 1; k++) A[n - 1][n - 1] -= A[n - 1][k] * A[k][n - 1];
	if (A[n - 1][n - 1] == 0.0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：方程组有无穷解" << endl;
#else
		cout << "Error: Equations have infinite solutions." << endl;
#endif
		throw 0;
	}
	if (!print_flag) return;  //Only decomposition
	for (int k = 0; k < n - 1; k++) b[n - 1] -= A[n - 1][k] * b[k];
	x[n - 1] = b[n - 1] / A[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = b[i];
		for (int j = i + 1; j < n; j++) x[i] -= A[i][j] * x[j];
		x[i] /= A[i][i];
	}
#ifdef CHINESE_VERSION
	if (print_flag) cout << "\n求解完成" << endl;
#else
	if (print_flag) cout << "\nFinish solving." << endl;
#endif
}

void Doolittle::decomposition(vector<vector<double>> AA, vector<vector<double>>& L, vector<vector<double>>& R) {
	int mm = AA.size();
	int nn = AA[0].size();
	for (int i = 0; i < mm; i++) if (AA[i].size() != nn) {
#ifdef CHINESE_VERSION
		cout << "错误：矩阵阶数不一致" << endl;
#else
		cout << "Error: The order of matrix is inconsistent" << endl;
#endif
		throw 0;
	}
	A_init(mm, nn);
	copy_A(AA);
	vector<double> bb(m);
	copy_b(bb);
	calc();
	vector<vector<double>> temp(m);
	for (int i = 0; i < m; i++) temp[i].resize(n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) temp[i][j] = 1.0;
			else if (i < j) temp[i][j] = 0.0;
			else temp[i][j] = A[i][j];
			L[i][j] = temp[i][j];
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i > j) temp[i][j] = 0.0;
			else temp[i][j] = A[i][j];
			R[i][j] = temp[i][j];
		}
	}
}

void Doolittle::L_solve(vector<vector<double>> LL, vector<double>& yy, std::vector<double> bb) {
	int nn = LL.size();
	yy[0] = bb[0];
	for (int i = 1; i < nn; i++) {
		for (int j = 0; j < i; j++) bb[i] -= LL[i][j] * yy[j];
		yy[i] = bb[i];
	}
}

void Doolittle::U_solve(vector<vector<double>> UU, vector<double>& xx, std::vector<double> yy) {
	int nn = UU.size();
	xx[nn - 1] = yy[nn - 1] / UU[nn - 1][nn - 1];
	for (int i = nn - 2; i >= 0; i--)
	{
		xx[i] = yy[i];
		for (int j = i + 1; j < nn; j++) xx[i] -= UU[i][j] * xx[j];
		xx[i] /= UU[i][i];
	}
}

void Doolittle::out_LU()
{
	double** temp = new double* [n];
	for (int i = 0; i < n; i++) temp[i] = new double[n];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) temp[i][j] = 1.0;
			else if (i < j) temp[i][j] = 0.0;
			else temp[i][j] = A[i][j];
		}
	}
	cout << "\nL (" << m << "*" << n << ") = " << endl;
	for (int i = 0; i < n; i++)
	{
		cout << "\t";
		for (int j = 0; j < n; j++)
		{
			cout.width(15);
			cout << temp[i][j];
		}
		cout << endl;
	}
	fl << "\nL (" << m << "*" << n << ") =\n\n";
	fl.set_mat_size(n, n) << temp;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i > j) temp[i][j] = 0.0;
			else temp[i][j] = A[i][j];
		}
	}
	cout << "\nU (" << m << "*" << n << ") = " << endl;
	for (int i = 0; i < n; i++)
	{
		cout << "\t";
		for (int j = 0; j < n; j++)
		{
			cout.width(15);
			cout << temp[i][j];
		}
		cout << endl;
	}
	fl << "\nU (" << m << "*" << n << ") =\n\n";
	fl.set_mat_size(n, n) << temp;
	for (int i = 0; i < n; i++) delete[] temp[i];
	delete[] temp;
}

void Doolittle::out_result()
{
#ifdef CHINESE_VERSION
	cout << "\n通过杜利特尔分解法解得" << endl;
	fl << "\n通过杜利特尔分解法解得\n";
#else
	cout << "\nSolved by Doolittle decomposition." << endl;
	fl << "\nSolved by Doolittle decomposition.\n";
#endif
	out_LU();
	out_x();
}

Doolittle::~Doolittle() {
	delete_Ax_b();
}

//4.Cholesky decomposition

Cholesky::Cholesky() {
	init();
	enter_Ab();
	out_Ab();
}

Cholesky::Cholesky(bool no_init) {
	if (no_init) return;
	else Cholesky();
}

void Cholesky::calc()
{
#ifdef CHINESE_VERSION
	cout << "\n正在通过楚列斯基分解求解 Ax=b ..." << endl;
#else
	cout << "\nSolving Ax=b by Cholesky decomposition..." << endl;
#endif
	if (!check_symmetry())
	{
#ifdef CHINESE_VERSION
		cout << "错误：A不是对称矩阵" << endl;
#else
		cout << "Error: A is not symmetrical." << endl;
#endif
		throw 0;
	}
	if (A[0][0] <= 0.0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：A不是对称正定矩阵" << endl;
#else
		cout << "Error: A is not spd matrix." << endl;
#endif
		throw 0;
	}
	A[0][0] = sqrt(A[0][0]);
	b[0] /= A[0][0];
	for (int i = 1; i < n; i++) A[i][0] /= A[0][0];
	for (int j = 1; j < n - 1; j++)
	{
		for (int k = 0; k < j; k++) A[j][j] -= sqr(A[j][k]);
		if (A[j][j] <= 0.0)
		{
#ifdef CHINESE_VERSION
			cout << "错误：A不是对称正定矩阵" << endl;
#else
			cout << "Error: A is not spd matrix." << endl;
#endif
			throw 0;
		}
		A[j][j] = sqrt(A[j][j]);
		for (int i = j + 1; i < n; i++)
		{
			for (int k = 0; k < j; k++) A[i][j] -= A[i][k] * A[j][k];
			A[i][j] /= A[j][j];
		}
		for (int k = 0; k < j; k++) b[j] -= b[k] * A[j][k];
		b[j] /= A[j][j];
	}
	for (int k = 0; k < n - 1; k++) A[n - 1][n - 1] -= sqr(A[n - 1][k]);
	if (A[n - 1][n - 1] <= 0.0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：A不是对称正定矩阵" << endl;
#else
		cout << "Error: A is not spd matrix." << endl;
#endif
		throw 0;
	}
	A[n - 1][n - 1] = sqrt(A[n - 1][n - 1]);
	for (int k = 0; k < n - 1; k++) b[n - 1] -= b[k] * A[n - 1][k];
	b[n - 1] /= A[n - 1][n - 1];
	x[n - 1] = b[n - 1] / A[n - 1][n - 1];
	for (int j = n - 2; j >= 0; j--)
	{
		x[j] = b[j];
		for (int i = j + 1; i < n; i++) x[j] -= A[i][j] * x[i];
		x[j] /= A[j][j];
	}
#ifdef CHINESE_VERSION
	cout << "\n求解完成" << endl;
#else
	cout << "\nFinish solving." << endl;
#endif
}

void Cholesky::out_G() {
	double** temp = new double* [n];
	for (int i = 0; i < n; i++) temp[i] = new double[n];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i < j) temp[i][j] = 0.0;
			else temp[i][j] = A[i][j];
		}
	}
	cout << "\nG (" << m << "*" << n << ") = " << endl;
	for (int i = 0; i < n; i++)
	{
		cout << "\t";
		for (int j = 0; j < n; j++)
		{
			cout.width(15);
			cout << temp[i][j];
		}
		cout << endl;
	}
	fl << "\nG (" << m << "*" << n << ") =\n\n";
	fl.set_mat_size(n, n) << temp;
	for (int i = 0; i < n; i++) delete[] temp[i];
	delete[] temp;
}

void Cholesky::out_result()
{
#ifdef CHINESE_VERSION
	cout << "\n通过楚列斯基分解法解得" << endl;
	fl << "\n通过楚列斯基分解法解得\n";
#else
	cout << "\nSolved by Cholesky decomposition." << endl;
	fl << "\nSolved by Cholesky decomposition.\n";
#endif
	out_G();
	out_x();
}

Cholesky::~Cholesky() {
	delete_Ax_b();
}

//5.Improved Square Root

Improved_sqrt::Improved_sqrt() {
	init();
	enter_Ab();
	out_Ab();
}

Improved_sqrt::Improved_sqrt(bool no_init) {
	if (no_init) return;
	else Improved_sqrt();
}

void Improved_sqrt::calc()
{
#ifdef CHINESE_VERSION
	cout << "\n正在通过改进平方根法求解 Ax=b ..." << endl;
#else
	cout << "\nSolving Ax=b by improved square root..." << endl;
#endif
	if (!check_symmetry())
	{
#ifdef CHINESE_VERSION
		cout << "错误：A不是对称矩阵" << endl;
#else
		cout << "Error: A is not symmetrical." << endl;
#endif
		throw 0;
	}
	if (A[0][0] == 0.0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：A的顺序主子式为0" << endl;
#else
		cout << "Error: Sequential principal minor of A is 0." << endl;
#endif
		throw 0;
	}
	for (int j = 1; j < n; j++) A[j][0] = A[0][j] / A[0][0];
	for (int i = 1; i < n - 1; i++)
	{
		for (int k = 0; k < i; k++) A[i][i] -= A[i][k] * A[k][i];
		if (A[i][i] == 0.0)
		{
#ifdef CHINESE_VERSION
			cout << "错误：A的顺序主子式为0" << endl;
#else
			cout << "Error: Sequential principal minor of A is 0." << endl;
#endif
			throw 0;
		}
		for (int j = i + 1; j < n; j++)
		{
			for (int k = 0; k < i; k++) A[i][j] -= A[i][k] * A[k][j];
			A[j][i] = A[i][j] / A[i][i];
		}
	}
	for (int k = 0; k < n - 1; k++) A[n - 1][n - 1] -= A[n - 1][k] * A[k][n - 1];
	if (A[n - 1][n - 1] == 0.0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：A的顺序主子式为0" << endl;
#else
		cout << "Error: Sequential principal minor of A is 0." << endl;
#endif
		throw 0;
	}
	for (int i = 1; i < n; i++) for (int k = 0; k < i; k++) b[i] -= A[i][k] * b[k];
	x[n - 1] = b[n - 1] / A[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = b[i];
		for (int k = i + 1; k < n; k++) x[i] -= A[i][k] * x[k];
		x[i] /= A[i][i];
	}
#ifdef CHINESE_VERSION
	cout << "\n求解完成" << endl;
#else
	cout << "\nFinish solving." << endl;
#endif
}

void Improved_sqrt::out_LD() {
	double** temp = new double* [n];
	for (int i = 0; i < n; i++) temp[i] = new double[n];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) temp[i][j] = 1.0;
			else if (i < j) temp[i][j] = 0.0;
			else temp[i][j] = A[i][j];
		}
	}
	cout << "\nL (" << m << "*" << n << ") = " << endl;
	for (int i = 0; i < n; i++)
	{
		cout << "\t";
		for (int j = 0; j < n; j++)
		{
			cout.width(15);
			cout << temp[i][j];
		}
		cout << endl;
	}
	fl << "\nL (" << m << "*" << n << ") =\n\n";
	fl.set_mat_size(n, n) << temp;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) temp[i][j] = A[i][j];
			else temp[i][j] = 0.0;
		}
	}
	cout << "\nD (" << m << "*" << n << ") = " << endl;
	for (int i = 0; i < n; i++)
	{
		cout << "\t";
		for (int j = 0; j < n; j++)
		{
			cout.width(15);
			cout << temp[i][j];
		}
		cout << endl;
	}
	fl << "\nD (" << m << "*" << n << ") =\n\n";
	fl.set_mat_size(n, n) << temp;
	for (int i = 0; i < n; i++) delete[] temp[i];
	delete[] temp;
}

void Improved_sqrt::out_result()
{
#ifdef CHINESE_VERSION
	cout << "\n通过改进平方根法解得" << endl;
	fl << "\n通过改进平方根法解得\n";
#else
	cout << "\nSolved by improved square root." << endl;
	fl << "\nSolved by improved square root.\n";
#endif
	out_LD();
	out_x();
}

Improved_sqrt::~Improved_sqrt() {
	delete_Ax_b();
}

//6.Chasing

Chasing::Chasing() {
	init();
#ifdef CHINESE_VERSION
	cout << "\n是否要构造对称三角矩阵？（1 = 是 , 0 = 否）" << endl;
#else
	cout << "\nWant to construct a symmetric tridiagonal matrix? (1 = yes, 0 = no)" << endl;
#endif
	int temp_flag;
	temp_flag = in_int();
	if (temp_flag == 1) {
		construct_sym_tri();
		save_Ab_sym_tri();
		return;
	}
	enter_Ab();
	out_Ab();
}

Chasing::Chasing(bool no_init) {
	if (no_init) return;
	else Chasing();
}

void Chasing::calc()
{
#ifdef CHINESE_VERSION
	cout << "\n正在通过追赶法求解 Ax=b ..." << endl;
#else
	cout << "\nSolving Ax=b by chasing..." << endl;
#endif
	if (!check_tridiagonal())
	{
#ifdef CHINESE_VERSION
		cout << "错误：A不是三角矩阵" << endl;
#else
		cout << "Error: A is not a tridiagonal matrix." << endl;
#endif
		throw 0;
	}
	for (int i = 1; i < n; i++)
	{
		A[i][i - 1] /= A[i - 1][i - 1];
		A[i][i] -= A[i][i - 1] * A[i - 1][i];
		b[i] -= A[i][i - 1] * b[i - 1];
	}
	for (int i = 0; i < n; i++) if (A[i][i] == 0.0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：A不能通过追赶法求解" << endl;
#else
		cout << "Error: A can't be solved by chasing method." << endl;
#endif
		throw 0;
	}
	x[n - 1] = b[n - 1] / A[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--) x[i] = (b[i] - A[i][i + 1] * x[i + 1]) / A[i][i];
#ifdef CHINESE_VERSION
	cout << "\n求解完成" << endl;
#else
	cout << "\nFinish solving." << endl;
#endif
}

void Chasing::save_Ab_sym_tri() {
	double x1, x2;
	x1 = A[0][0];
	x2 = A[0][1];
	fl.set_left_right(true);
	fl << "\nA (" << n << "*" << n << ") = \n\n\t";
	fl.set_width(15);
	fl << x1;
	fl.set_width(15);
	fl << x2;
	fl << "\n\t";
	fl.set_width(15);
	fl << x2;
	fl.set_width(15);
	fl << x1;
	fl.set_width(15);
	fl << x2;
	fl << "\n\t";
	fl.set_width(15);
	fl << " ";
	fl.set_width(15);
	fl << "...";
	fl.set_width(15);
	fl << "...";
	fl.set_width(15);
	fl << "...";
	fl << "\n\t";
	fl.set_width(15);
	fl << " ";
	fl.set_width(15);
	fl << " ";
	fl.set_width(15);
	fl << x2;
	fl.set_width(15);
	fl << x1;
	fl.set_width(15);
	fl << x2;
	fl << "\n\t";
	fl.set_width(15);
	fl << " ";
	fl.set_width(15);
	fl << " ";
	fl.set_width(15);
	fl << " ";
	fl.set_width(15);
	fl << x2;
	fl.set_width(15);
	fl << x1;
	fl << "\n";
	x1 = b[0];
	x2 = b[1];
	fl << "\nb (" << n << "*" << 1 << ") = \n\n\t";
	fl << x1 << "\n\t" << x2 << "\n\t.\n\t.\n\t.\n\t" << x2 << "\n\t" << x1 << "\n";
}

void Chasing::out_LU()
{
	double** temp = new double* [n];
	for (int i = 0; i < n; i++) temp[i] = new double[n];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) temp[i][j] = 1.0;
			else if (i < j) temp[i][j] = 0.0;
			else temp[i][j] = A[i][j];
		}
	}
	cout << "\nL (" << m << "*" << n << ") = " << endl;
	for (int i = 0; i < n; i++)
	{
		cout << "\t";
		for (int j = 0; j < n; j++)
		{
			cout.width(15);
			cout << temp[i][j];
		}
		cout << endl;
	}
	fl << "\nL (" << m << "*" << n << ") =\n\n";
	fl.set_mat_size(n, n) << temp;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i > j) temp[i][j] = 0.0;
			else temp[i][j] = A[i][j];
		}
	}
	cout << "\nU (" << m << "*" << n << ") = " << endl;
	for (int i = 0; i < n; i++)
	{
		cout << "\t";
		for (int j = 0; j < n; j++)
		{
			cout.width(15);
			cout << temp[i][j];
		}
		cout << endl;
	}
	fl << "\nU (" << m << "*" << n << ") =\n\n";
	fl.set_mat_size(n, n) << temp;
	for (int i = 0; i < n; i++) delete[] temp[i];
	delete[] temp;
}

void Chasing::out_result()
{
#ifdef CHINESE_VERSION
	cout << "\n通过追赶法解得" << endl;
	fl << "\n通过追赶法解得\n";
#else
	cout << "\nSolved by chasing." << endl;
	fl << "\nSolved by chasing.\n";
#endif
	out_LU();
	out_x();
}

Chasing::~Chasing() {
	delete_Ax_b();
}

//7.Givens transformation (QR)

Givens::Givens() {
	init();
	enter_Ab();
	out_Ab();
	Q = new double* [m];
	for (int i = 0; i < m; i++) Q[i] = new double[m];
	for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) {
		Q[i][j] = 0.0;
		if (i == j) Q[i][j] = 1.0;
	}
}

Givens::Givens(bool no_init) {
	if (no_init) return;
	else Givens();
}

void Givens::calc()
{
#ifdef CHINESE_VERSION
	cout << "\n正在通过吉文斯变换求解 Ax=b ..." << endl;
#else
	cout << "\nSolving Ax=b by Givens transformation..." << endl;
#endif
	double c, s, temp;
	for (int k = 0; k < min(m - 1, n); k++)
	{
		for (int i = k + 1; i < m; i++)
		{
			if (A[i][k] == 0) continue;
			c = A[k][k] / sqrt(sqr(A[k][k]) + sqr(A[i][k]));
			s = A[i][k] / sqrt(sqr(A[k][k]) + sqr(A[i][k]));
			for (int j = k; j < n; j++)
			{
				temp = A[k][j];
				A[k][j] = c * A[k][j] + s * A[i][j];
				A[i][j] = c * A[i][j] - s * temp;
			}
			for (int j = 0; j < m; j++)
			{
				temp = Q[k][j];
				Q[k][j] = c * Q[k][j] + s * Q[i][j];
				Q[i][j] = c * Q[i][j] - s * temp;
			}
		}
		if (A[k][k] == 0.0)
		{
#ifdef CHINESE_VERSION
			cout << "错误：rank(A)<n" << endl;
#else
			cout << "Error: rank(A)<n" << endl;
#endif
			throw 0;
		}
	}
	if (m == n && A[n - 1][n - 1] == 0.0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：rank(A)<n" << endl;
#else
		cout << "Error: rank(A)<n" << endl;
#endif
		throw 0;
	}
	for (int i = 0; i < m; i++)
	{
		for (int j = i + 1; j < m; j++)
		{
			temp = Q[i][j];
			Q[i][j] = Q[j][i];
			Q[j][i] = temp;
		}
	}
	double* y = new double[m];
	for (int i = 0; i < m; i++)
	{
		y[i] = Q[0][i] * b[0];
		for (int j = 1; j < m; j++)
		{
			y[i] += Q[j][i] * b[j];
		}
	}
	if (m > n)
	{
		for (int i = n; i < m; i++)
		{
			if (y[i] != 0.0)
			{
#ifdef CHINESE_VERSION
				cout << "提示：方程组是超定的，此处为最小二乘解" << endl;
#else
				cout << "Tips: Equations are overdetermined. Here's the least-squares solution." << endl;
#endif
				break;
			}
		}
	}
	x[n - 1] = y[n - 1] / A[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = y[i];
		for (int j = i + 1; j < n; j++)
		{
			x[i] -= A[i][j] * x[j];
		}
		x[i] /= A[i][i];
	}
	delete[]y;
#ifdef CHINESE_VERSION
	cout << "\n求解完成" << endl;
#else
	cout << "\nFinish solving." << endl;
#endif
}

void Givens::out_QR()
{
	cout << "\nQ (" << m << "*" << m << ") =" << endl;
	fl << "\nQ (" << m << "*" << m << ") = \n\n";
	for (int i = 0; i < m; i++)
	{
		cout << "\t";
		fl << "\t";
		for (int j = 0; j < n; j++)
		{
			cout.width(15);
			cout << Q[i][j];
			fl.set_width(15);
			fl << Q[i][j];
		}
		if (m > n)
		{
			cout.width(5);
			cout << "|";
			fl.set_width(5);
			fl << "|";
			for (int j = n; j < m; j++)
			{
				cout.width(15);
				cout << Q[i][j];
				fl.set_width(15);
				fl << Q[i][j];
			}
		}
		cout << endl;
		fl << "\n";
	}
	cout << "\nR (" << m << "*" << n << ") =" << endl;
	fl << "\nR (" << m << "*" << n << ") = \n\n";
	for (int i = 0; i < n; i++)
	{
		cout << "\t";
		fl << "\t";
		for (int j = 0; j < n; j++)
		{
			cout.width(15);
			cout << A[i][j];
			fl.set_width(15);
			fl << A[i][j];
		}
		cout << endl;
		fl << "\n";
	}
	if (m > n)
	{
		cout << "\t";
		fl << "\t";
		for (int i = 0; i < n; i++)
		{
			cout.width(15);
			cout << "---------------";
			fl.set_width(15);
			fl << "---------------";
		}
		cout << endl;
		fl << "\n";
		for (int i = n; i < m; i++)
		{
			cout << "\t";
			fl << "\t";
			for (int j = 0; j < n; j++)
			{
				cout.width(15);
				cout << A[i][j];
				fl.set_width(15);
				fl << A[i][j];
			}
			cout << endl;
			fl << "\n";
		}
	}
}

void Givens::out_result()
{
#ifdef CHINESE_VERSION
	cout << "\n通过吉文斯变换解得" << endl;
	fl << "\n通过吉文斯变换解得\n";
#else
	cout << "\nSolved by Givens transformation." << endl;
	fl << "\nSolved by Givens transformation.\n";
#endif
	out_QR();
	out_x();
}

Givens::~Givens() {
	delete_Ax_b();
	for (int i = 0; i < m; i++) delete[] Q[i];
	delete[] Q;
}

//8.Householder transformation

Householder::Householder() {
	init();
	enter_Ab();
	out_Ab();
	Q = new double* [m];
	for (int i = 0; i < m; i++) Q[i] = new double[m];
	for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) {
		Q[i][j] = 0.0;
		if (i == j) Q[i][j] = 1.0;
	}
	d = new double[n];
	alpha = new double[n];
	print_flag = true;
}

Householder::Householder(bool no_init) {
	if (no_init) return;
	else Householder();
}

Householder::Householder(int mm, int nn, std::vector<std::vector<double>> AA, std::vector<double> bb) {
	if (mm <= 0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：矩阵A的行数不是正整数" << endl;
#else
		cout << "Error: The number of rows of A is not a positive integer." << endl;
#endif
		throw 1;
	}
	if (nn <= 0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：矩阵A的列数不是正整数" << endl;
#else
		cout << "Error: The number of columns of A is not a positive integer." << endl;
#endif
		throw 1;
	}
	if (mm < nn)
	{
#ifdef CHINESE_VERSION
		cout << "错误：方程组有无穷解" << endl;
#else
		cout << "Error: Equations have infinite solutions." << endl;
#endif
		throw 1;
	}
	A_init(mm, nn);
	copy_A(AA);
	copy_b(bb);
	Q = new double* [m];
	for (int i = 0; i < m; i++) Q[i] = new double[m];
	for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) {
		Q[i][j] = 0.0;
		if (i == j) Q[i][j] = 1.0;
	}
	d = new double[n];
	alpha = new double[n];
	print_flag = false;
}

void Householder::calc()
{
#ifdef CHINESE_VERSION
	if (print_flag) cout << "\n正在通过豪斯霍尔德变换求解 Ax=b ..." << endl;
#else
	if (print_flag) cout << "\nSolving Ax=b by Householder transformation..." << endl;
#endif
	double sigma, beta;
	for (int k = 0; k < n; k++)
	{
		if (k == n - 1 && m == n)
		{
			d[n - 1] = A[n - 1][n - 1];
			break;
		}
		sigma = 0.0;
		for (int i = k; i < m; i++) sigma += sqr(A[i][k]);
		sigma = -sign(A[k][k]) * sqrt(sigma);
		d[k] = sigma;
		if (d[k] == 0.0)
		{
#ifdef CHINESE_VERSION
			cout << "错误：rank(A)<n" << endl;
#else
			cout << "Error: rank(A)<n" << endl;
#endif
			throw 0;
		}
		alpha[k] = sigma * (sigma - A[k][k]);
		A[k][k] -= sigma;
		for (int j = k + 1; j < n; j++)
		{
			beta = A[k][k] * A[k][j];
			for (int i = k + 1; i < m; i++) beta += A[i][k] * A[i][j];
			beta /= alpha[k];
			A[k][j] -= beta * A[k][k];
			for (int i = k + 1; i < m; i++)
			{
				A[i][j] -= beta * A[i][k];
			}
		}
	}
	if (d[n - 1] == 0.0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：rank(A)<n" << endl;
#else
		cout << "Error: rank(A)<n" << endl;
#endif
		throw 0;
	}
	double* temp1 = new double[m];
	double* temp2 = new double[m];
	for (int k = 0; k < min(m - 1, n); k++)
	{
		for (int j = 0; j < m; j++)
		{
			for (int i = 0; i < m; i++)
			{
				temp1[i] = Q[i][j];
			}
			for (int i = k; i < m; i++)
			{
				for (int ii = 0; ii < m; ii++) temp2[ii] = 0.0;
				for (int ii = k; ii < m; ii++)
				{
					if (ii == i) temp2[ii] += 1.0;
					temp2[ii] -= (A[ii][k] * A[i][k]) / alpha[k];
				}
				Q[i][j] = 0.0;
				for (int ii = 0; ii < m; ii++)
				{
					Q[i][j] += temp1[ii] * temp2[ii];
				}
			}
		}
	}
	delete[]temp1;
	delete[]temp2;
	double temp;
	double* y = new double[m];
	for (int i = 0; i < m; i++)
	{
		for (int j = i + 1; j < m; j++)
		{
			temp = Q[i][j];
			Q[i][j] = Q[j][i];
			Q[j][i] = temp;
		}
	}
	for (int i = 0; i < m; i++)
	{
		y[i] = Q[0][i] * b[0];
		for (int j = 1; j < m; j++)

		{
			y[i] += Q[j][i] * b[j];
		}
	}
	if (m > n)
	{
		for (int i = n; i < m; i++)
		{
			if (y[i] != 0.0)
			{
#ifdef CHINESE_VERSION
				cout << "提示：方程组是超定的，此处为最小二乘解" << endl;
#else
				cout << "Tips: Equations are overdetermined. Here's the least-squares solution." << endl;
#endif
				break;
			}
		}
	}
	x[n - 1] = y[n - 1] / d[n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = y[i];
		for (int j = i + 1; j < n; j++)
		{
			x[i] -= A[i][j] * x[j];
		}
		x[i] /= d[i];
	}
	delete[]y;
#ifdef CHINESE_VERSION
	if (print_flag) cout << "\n求解完成" << endl;
#else
	if (print_flag) cout << "\nFinish solving." << endl;
#endif
}

void Householder::out_QR()
{
	cout << "\nQ (" << m << "*" << m << ") =" << endl;
	fl << "\nQ (" << m << "*" << m << ") = \n\n";
	for (int i = 0; i < m; i++)
	{
		cout << "\t";
		fl << "\t";
		for (int j = 0; j < n; j++)
		{
			cout.width(15);
			cout << Q[i][j];
			fl.set_width(15);
			fl << Q[i][j];
		}
		if (m > n)
		{
			cout.width(5);
			cout << "|";
			fl.set_width(5);
			fl << "|";
			for (int j = n; j < m; j++)
			{
				cout.width(15);
				cout << Q[i][j];
				fl.set_width(15);
				fl << Q[i][j];
			}
		}
		cout << endl;
		fl << "\n";
	}
	cout << "\nR (" << m << "*" << n << ") =" << endl;
	fl << "\nR (" << m << "*" << n << ") = \n\n";
	for (int i = 0; i < n; i++)
	{
		cout << "\t";
		fl << "\t";
		for (int j = 0; j < i; j++)
		{
			cout.width(15);
			cout << 0;
			fl.set_width(15);
			fl << 0.0;
		}
		cout.width(15);
		cout << d[i];
		fl.set_width(15);
		fl << d[i];
		for (int j = i + 1; j < n; j++)
		{
			cout.width(15);
			cout << A[i][j];
			fl.set_width(15);
			fl << A[i][j];
		}
		cout << endl;
		fl << "\n";
	}
	if (m > n)
	{
		cout << "\t";
		fl << "\t";
		for (int i = 0; i < n; i++)
		{
			cout.width(15);
			cout << "---------------";
			fl.set_width(15);
			fl << "---------------";
		}
		cout << endl;
		fl << "\n";
		for (int i = n; i < m; i++)
		{
			cout << "\t";
			fl << "\t";
			for (int j = 0; j < n; j++)
			{
				cout.width(15);
				cout << 0;
				fl.set_width(15);
				fl << 0.0;
			}
			cout << endl;
			fl << "\n";
		}
	}
}

vector<double> Householder::get_x() {
	vector<double> xx;
	xx.resize(n);
	for (int i = 0; i < n; i++) xx[i] = x[i];
	return xx;
}

void Householder::out_result()
{
#ifdef CHINESE_VERSION
	cout << "\n通过豪斯霍尔德变换解得" << endl;
	fl << "\n通过豪斯霍尔德变换解得\n";
#else
	cout << "\nSolved by Householder transformation." << endl;
	fl << "\nSolved by Householder transformation.\n";
#endif
	out_QR();
	out_x();
}

Householder::~Householder() {
	delete_Ax_b();
	for (int i = 0; i < m; i++) delete[] Q[i];
	delete[] Q;
	delete[] d;
	delete[] alpha;
}
