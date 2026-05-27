#include <iostream>
#include <cmath>
#include "iterationmethod.h"

using namespace std;

//construction function
void Iteration_method::init() {
	fl.init("Iteration_method.txt");
	cout.setf(ios::left);
	int nn;
#ifdef CHINESE_VERSION
	cout << "\n请输入系数矩阵A的阶数(n*n)：\nn = ";
#else
	cout<<"\nPlease input the order of coefficient matrix A(n*n):\nn = ";
#endif
	nn = in_int();
	if(nn<=0)
	{
#ifdef CHINESE_VERSION
		cout << "错误：请输入一个正整数。" << endl;
#else
		cout<<"Error: Please enter a positive integer."<<endl;
#endif
		throw 0;
	}
	A_init(nn);
#ifdef CHINESE_VERSION
	cout << "\n是否构造对称三对角矩阵？(1 = 是, 0 = 否)" << endl;
#else
	cout << "\nWant to construct a symmetric tridiagonal matrix? (1 = yes, 0 = no)" << endl;
#endif
	int temp;
	temp = in_int();
	if (temp == 1) {
		construct_sym_tri();
		init_x(0.0);
		save_Ab_sym_tri();
	}
	else {
		enter_Ab();
		init_x();
		out_Ab();
	}
	itcounter = 0;
	set_eps();
	fl.set_precision(eps);
	set_max();
	itrerr = new double[maxcounter];
}

//save matrices A, b after constructing them as symmetric triangular matrices
void Iteration_method::save_Ab_sym_tri() {
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

//cout and save matrices A, b
void Iteration_method::out_Ab() {
	cout << "\nA (" << n << "*" << n << ") =" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << "\t";
		for (int j = 0; j < n; j++)
		{

			cout.width(15);
			cout << A[i][j];
		}
		cout << endl;
	}
	cout << "\nb (" << n << "*" << 1 << ") =" << endl;
	for (int i = 0; i < n; i++) cout << "\t" << b[i] << endl;
	fl << "\nA =\n\n";
	fl.set_mat_size(n, n) << A;
	fl << "\nb = \n\n";
	fl.set_array_len(n) << b;
}

//cout and save matrix x
void Iteration_method::out_x()
{
	cout << "\nx (" << n << "*" << 1 << ") =" << endl;
	for (int i = 0; i < n; i++) cout << "\t" << x[i] << endl;
	fl << "\nx = \n\n";
	fl.set_array_len(n) << x;
}

//set the error limit eps
void Iteration_method::set_eps() {
#ifdef CHINESE_VERSION
	cout << "\n请设置误差限 (如 0.001 , 1e-6 ): eps = ";
#else
	cout << "\nPlease set the error limit ( e.g. 0.001 , 1e-6 ): eps = ";
#endif
	double temp;
	cin >> temp;
	eps = temp;
	if (eps < 1.0) {
		int pre;
		cout.setf(iostream::fixed);
		pre = -1 * floor(log10(eps));
		cout.precision(pre);
	}
}

//set the maximum iteration times
void Iteration_method::set_max() {
#ifdef CHINESE_VERSION
	cout << "\n请设置最大迭代次数：maxcounter = ";
#else
	cout << "\nPlease set the maximum permission iterative number: maxcounter = ";
#endif
	int t = in_int();
	maxcounter = t;
}

//if the actual iteration times is lower than the maximum, then resize the array saving the iteration error
void Iteration_method::resize_itrerr() {
	double* temp = new double[itcounter];
	for (int i = 0; i < itcounter; i++) temp[i] = itrerr[i];
	delete[] itrerr;
	itrerr = new double[itcounter];
	for (int i = 0; i < itcounter; i++) itrerr[i] = temp[i];
	delete[] temp;
}

//cout and save results
void Iteration_method::out_result() {
	out_x();
	if (itcounter < maxcounter) {
		resize_itrerr();
#ifdef CHINESE_VERSION
		cout << "\n迭代次数 = " << itcounter << endl;
		fl << "\n迭代次数 = " << itcounter << "\n";
#else
		cout << "\nIteration times = " << itcounter << endl;
		fl << "\nIteration times = " << itcounter << "\n";
#endif
	}
	else {
#ifdef CHINESE_VERSION
		cout << "\n迭代次数 = " << itcounter << " (已达到最大迭代次数)" <<endl;
		fl << "\n迭代次数 = " << itcounter << " (已达到最大迭代次数)\n";
#else
		cout << "\nIteration times = " << itcounter << " (reach the maximum iteration times)" <<endl;
		fl << "\nIteration times = " << itcounter << " (reach the maximum iteration times)\n";
#endif
	}
#ifdef CHINESE_VERSION
	cout << "\n误差限为 eps = " << eps << endl;
	fl << "\n误差限为 eps = " << eps << "\n";
	cout << "\n误差表如下：\n\n迭代次数 k\t\t误差" << endl;
	for (int i = 0; i < itcounter; i++) cout << "\t" << i + 1 << "\t\t\t" << itrerr[i] << endl;
	fl << "\n误差表如下：\n\n迭代次数 k\t\t误差\n";
#else
	cout << "\nThe error limit is eps = " << eps << endl;
	fl << "\nThe error limit is eps = " << eps << "\n";
	cout << "\nThe error table is as follow:\n\niteration times k\t\terror" << endl;
	for (int i = 0; i < itcounter; i++) cout << "\t" << i + 1 << "\t\t\t" << itrerr[i] << endl;
	fl << "\nThe error table is as follow:\n\niteration times k\t\terror\n";
#endif
	for (int i = 0; i < itcounter; i++) fl << "\t" << i + 1 << "\t\t\t" << itrerr[i] << "\n";
	generate_m();
}

//generate .m file to get the figure of iteration error in MATLAB
void Iteration_method::generate_m() {
#ifdef CHINESE_VERSION
	cout << "\n是否生成 .m 文件用于在MATLAB中绘制迭代误差曲线？(1 = 是, 0 = 否)" << endl;
#else
	cout << "\nWant to generate .m file to get the figure of iteration error in MATLAB? (1 = Yes , 0 = No)" << endl;
#endif
	int flag = in_int();
	if (flag == 1) {
		filelog gm;
		gm.init("Iteration_method.m");
		gm << "\nfigure\n";
		gm << "x = 1:" << itcounter << ";\n";
		gm << "data = [";
		for (int i = 0; i < itcounter; i++) gm << itrerr[i] << " ";
		gm << "];\n";
		gm << "y=log(data);\n";
		gm << "plot(y,'rx-')\n";
		gm << "xlabel('Iterative times (times)');\n";
		gm << "ylabel('Logarithm of error');\n";
		gm << "title('The curve of error varying with iteration times');\n";
	}
}

void Iteration_method::exchange_diag_no_0(int i) {
	int temp = 0;
	while (A[i + temp][i] == 0.0)
	{
		temp++;
		if (i + temp == m) break;
	}
	if (i + temp < m) exchange_row(i, i + temp);
	else
	{
		temp = 0;
		while (1) {
			while (A[i + temp][i] == 0.0) {
				temp--;
				if (i + temp < 0) break;
			}
			if (A[i][i + temp] != 0.0 || i + temp < 0) break;
		}
		if (i + temp < 0) {
#ifdef CHINESE_VERSION
			cout << "错误：由于对角线有零元素，无法使用基本迭代法求解该方程组。" << endl;
#else
			cout << "Error: Equations cannot be solved by basic iteration methods because of the zero element of the diagonal." << endl;
#endif
			throw 0;
		}
		else exchange_row(i, i + temp);
	}
}

//1.Jacobi iteration

Jacobi::Jacobi() {
	init();
}

void Jacobi::calc() {
#ifdef CHINESE_VERSION
	cout << "\n正在用Jacobi迭代法求解 Ax=b ..." << endl;
#else
	cout << "\nSolving Ax=b by Jacobi iterative method..." << endl;
#endif
	double* x0 = new double[n];
	double err;
	for (int i = 0; i < n; i++) {
		if (A[i][i] == 0.0) exchange_diag_no_0(i);
		b[i] /= A[i][i];
		for (int j = 0; j < n; j++) {
			if (i == j) continue;
			A[i][j] *= -1.0 / A[i][i];
		}
		A[i][i] = 0.0;
	}
	do {
		for (int i = 0; i < n; i++) x0[i] = x[i];
		for (int i = 0; i < n; i++) {
			x[i] = 0.0;
			for (int j = 0; j < n; j++) x[i] += A[i][j] * x0[j];
			x[i] += b[i];
		}
		for (int i = 0; i < n; i++) x0[i] = x[i] - x0[i];
		err = vecnorm2(x0, n);
		itrerr[itcounter] = err;
		itcounter++;
	} while (err > eps && itcounter < maxcounter);
	delete[]x0;
#ifdef CHINESE_VERSION
	cout << "\n求解完成。" << endl;
#else
	cout << "\nFinish solving." << endl;
#endif
}

void Jacobi::out_result() {
#ifdef CHINESE_VERSION
	cout << "\n通过Jacobi迭代法求解。" << endl;
	fl << "\n通过Jacobi迭代法求解。\n";
#else
	cout << "\nSolved by Jacobi iteration method." << endl;
	fl << "\nSolved by Jacobi iteration method.\n";
#endif
	Iteration_method::out_result();
}

Jacobi::~Jacobi() {
	delete_Ax_b();
	delete[]itrerr;
}

//2.Gauss-Seidel iteration

Gauss_Seidel::Gauss_Seidel() {
	init();
}

void Gauss_Seidel::calc() {
#ifdef CHINESE_VERSION
	cout << "\n正在用Gauss-Seidel迭代法求解 Ax=b ..." << endl;
#else
	cout << "\nSolving Ax=b by Gauss-Samuel iterative method..." << endl;
#endif
	double* x0 = new double[n];
	double err;
	for (int i = 0; i < n; i++) {
		if (A[i][i] == 0.0) exchange_diag_no_0(i);
		b[i] /= A[i][i];
		for (int j = 0; j < n; j++) {
			if (i == j) continue;
			A[i][j] *= -1.0 / A[i][i];
		}
		A[i][i] = 0.0;
	}
	do {
		for (int i = 0; i < n; i++) x0[i] = x[i];
		for (int i = 0; i < n; i++) {
			x[i] = 0.0;
			for (int j = 0; j < n; j++) x[i] += A[i][j] * x[j];
			x[i] += b[i];
		}
		for (int i = 0; i < n; i++) x0[i] = x[i] - x0[i];
		err = vecnorm2(x0, n);
		itrerr[itcounter] = err;
		itcounter++;
	} while (err > eps && itcounter < maxcounter);
	delete[]x0;
#ifdef CHINESE_VERSION
	cout << "\n求解完成。" << endl;
#else
	cout << "\nFinish solving." << endl;
#endif
}

void Gauss_Seidel::out_result() {
#ifdef CHINESE_VERSION
	cout << "\n通过Gauss-Seidel迭代法求解。" << endl;
	fl << "\n通过Gauss-Seidel迭代法求解。\n";
#else
	cout << "\nSolved by Gauss-Seidel iteration method." << endl;
	fl << "\nSolved by Gauss-Seidel iteration method.\n";
#endif
	Iteration_method::out_result();
}

Gauss_Seidel::~Gauss_Seidel() {
	delete_Ax_b();
	delete[]itrerr;
}

//3.SOR iteration

SOR::SOR() {
	init();
	in_omega();
}

void SOR::in_omega() {
#ifdef CHINESE_VERSION
	cout << "\n请输入松弛因子 omega：" << endl;
#else
	cout << "\nPlease input omega:" << endl;
#endif
	cin >> omega;
}

void SOR::calc() {
#ifdef CHINESE_VERSION
	cout << "\n正在用SOR迭代法求解 Ax=b ..." << endl;
#else
	cout << "\nSolving Ax=b by SOR iterative method..." << endl;
#endif
	double* x0 = new double[n];
	double err;
	for (int i = 0; i < n; i++) {
		if (A[i][i] == 0.0) exchange_diag_no_0(i);
		b[i] *= omega / A[i][i];
		for (int j = 0; j < n; j++) {
			if (i == j) continue;
			A[i][j] *= -1.0 * omega / A[i][i];
		}
		A[i][i] = 0.0;
	}
	do {
		for (int i = 0; i < n; i++) x0[i] = x[i];
		for (int i = 0; i < n; i++) {
			x[i] *= 1.0 - omega;
			for (int j = 0; j < n; j++) x[i] += A[i][j] * x[j];
			x[i] += b[i];
		}
		for (int i = 0; i < n; i++) x0[i] = x[i] - x0[i];
		err = vecnorm2(x0, n);
		itrerr[itcounter] = err;
		itcounter++;
	} while (err > eps && itcounter < maxcounter);
	delete[]x0;
#ifdef CHINESE_VERSION
	cout << "\n求解完成。" << endl;
#else
	cout << "\nFinish solving." << endl;
#endif
}

void SOR::out_result() {
#ifdef CHINESE_VERSION
	cout << "\n通过SOR迭代法求解。" << endl;
	fl << "\n通过SOR迭代法求解。\n";
#else
	cout << "\nSolved by SOR iteration method." << endl;
	fl << "\nSolved by SOR iteration method.\n";
#endif
	Iteration_method::out_result();
}

SOR::~SOR() {
	delete_Ax_b();
	delete[]itrerr;
}

//4.steepest descent method

steepest_descent::steepest_descent() {
	init();
}

void steepest_descent::calc() {
#ifdef CHINESE_VERSION
	cout << "\n正在用最速下降法求解 Ax=b ..." << endl;
#else
	cout << "\nSolving Ax=b by steepest descent method..." << endl;
#endif
	if (!check_symmetry())
	{
#ifdef CHINESE_VERSION
		cout << "错误：A不是对称矩阵。" << endl;
#else
		cout << "Error: A is not symmetrical." << endl;
#endif
		throw 0;
	}
	double* x0 = new double[n];
	double* p = new double[n];
	double alpha;
	double temp;
	double err;
	do {
		for (int i = 0; i < n; i++) x0[i] = x[i];
		for (int i = 0; i < n; i++) {
			p[i] = b[i];
			for (int j = 0; j < n; j++) p[i] -= A[i][j] * x[j];
		}
		alpha = 0.0;
		for (int i = 0; i < n; i++) alpha += sqr(p[i]);
		temp = 0.0;
		for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) temp += A[i][j] * p[i] * p[j];
		alpha /= temp;
		for (int i = 0; i < n; i++) x[i] += alpha * p[i];
		for (int i = 0; i < n; i++) x0[i] = x[i] - x0[i];
		err = vecnorm2(x0, n);
		itrerr[itcounter] = err;
		itcounter++;
	} while (err > eps && itcounter < maxcounter);
	delete[]x0;
	delete[]p;
#ifdef CHINESE_VERSION
	cout << "\n求解完成。" << endl;
#else
	cout << "\nFinish solving." << endl;
#endif
}

void steepest_descent::out_result() {
#ifdef CHINESE_VERSION
	cout << "\n通过最速下降法求解。" << endl;
	fl << "\n通过最速下降法求解。\n";
#else
	cout << "\nSolved by steepest descent method." << endl;
	fl << "\nSolved by steepest descent method.\n";
#endif
	Iteration_method::out_result();
}

steepest_descent::~steepest_descent() {
	delete_Ax_b();
	delete[]itrerr;
}

//5.conjugate gradient method

conjugate_gradient::conjugate_gradient() {
	init();
	print_state = true;
}

conjugate_gradient::conjugate_gradient(double** AA, double* bb, int nn) {
	A_init(nn);
	copy_A(AA);
	copy_b(bb);
	itcounter = 0;
	eps = 0.5e-8;
	maxcounter = 5 * nn;
	init_x(0.0);
	itrerr = new double[maxcounter];
	print_state = false;
}

void conjugate_gradient::calc() {
#ifdef CHINESE_VERSION
	if (print_state) cout << "\n正在用共轭梯度法求解 Ax=b ..." << endl;
#else
	if (print_state) cout << "\nSolving Ax=b by conjugate gradient method..." << endl;
#endif
	if (!check_symmetry())
	{
#ifdef CHINESE_VERSION
		cout << "错误：A不是对称矩阵。" << endl;
#else
		cout << "Error: A is not symmetrical." << endl;
#endif
		throw 0;
	}
	double* r0 = new double[n];
	double* r = new double[n];
	double* d = new double[n];
	double alpha;
	double beta;
	double temp;
	double err;
	for (int i = 0; i < n; i++) {
		r0[i] = b[i];
		for (int j = 0; j < n; j++) r0[i] -= A[i][j] * x[j];
	}
	for (int i = 0; i < n; i++) d[i] = r0[i];
	alpha = 0.0;
	for (int i = 0; i < n; i++) alpha += sqr(r0[i]);
	temp = 0.0;
	for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) temp += A[i][j] * d[i] * d[j];
	alpha /= temp;
	for (int i = 0; i < n; i++) x[i] += alpha * d[i];
	for (int i = 0; i < n; i++) {
		r[i] = b[i];
		for (int j = 0; j < n; j++) r[i] -= A[i][j] * x[j];
	}
	err = vecnorm2(r, n);
	itrerr[itcounter] = err;
	itcounter++;
	while (err > eps && itcounter < 5 * n) {
		beta = sqr(vecnorm2(r, n) / vecnorm2(r0, n));
		for (int i = 0; i < n; i++) {
			d[i] *= beta;
			d[i] += r[i];
		}
		for (int i = 0; i < n; i++) r0[i] = r[i];
		alpha = 0.0;
		for (int i = 0; i < n; i++) alpha += sqr(r[i]);
		temp = 0.0;
		for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) temp += A[i][j] * d[i] * d[j];
		alpha /= temp;
		for (int i = 0; i < n; i++) x[i] += alpha * d[i];
		for (int i = 0; i < n; i++) {
			r[i] = b[i];
			for (int j = 0; j < n; j++) r[i] -= A[i][j] * x[j];
		}
		err = vecnorm2(r, n);
		itrerr[itcounter] = err;
		itcounter++;
	}
	delete[]r0;
	delete[]r;
	delete[]d;
#ifdef CHINESE_VERSION
	if (print_state) cout << "\n求解完成。" << endl;
#else
	if (print_state) cout << "\nFinish solving." << endl;
#endif
}

void conjugate_gradient::out_result() {
#ifdef CHINESE_VERSION
	cout << "\n通过共轭梯度法求解。" << endl;
	fl << "\n通过共轭梯度法求解。\n";
#else
	cout << "\nSolved by conjugate gradient method." << endl;
	fl << "\nSolved by conjugate gradient method.\n";
#endif
	Iteration_method::out_result();
}

conjugate_gradient::~conjugate_gradient() {
	delete_Ax_b();
	delete[]itrerr;
}

void conjugate_gradient::get_result(double* xx) {
	get_x(xx);
}

//6.GMRES
GMRES::GMRES() {
	init();
	in_m();
}

void GMRES::in_m() {
#ifdef CHINESE_VERSION
	cout << "\n请输入低维子空间维度 m = ";
#else
	cout << "\nPlease input a low dimension m = ";
#endif
	cin >> m;
}

void GMRES::calc() {
#ifdef CHINESE_VERSION
	cout << "\n正在用GMRES求解 Ax=b ..." << endl;
#else
	cout << "\nSolving Ax=b by GMRES..." << endl;
#endif
	double* r0 = new double[n];
	double** Vm = new double*[n];
	for(int i=0;i<n;i++) Vm[i] = new double[m+1];
	double** H_trans = new double*[m]; //transpose of H
	for(int i=0;i<m;i++) H_trans[i] = new double[i+2];
	double* e1 = new double[m+1];
	double c,s,tmpi,tmpj;
	double* y = new double[m];
	double err;
	for(int i=0;i<n;i++) {
		r0[i] = b[i];
		for(int j=0;j<n;j++) r0[i] -= A[i][j]*x[j];
	}
	err = vecnorm2(r0,n);
	while (err > eps && itcounter < maxcounter) {
		//Arnoldi process
		for(int i=0;i<n;i++) Vm[i][0] = r0[i]/err;
		for(int k=0;k<m;k++) {
			for(int i=0;i<n;i++) Vm[i][m] = 0.0;
			for(int i=0;i<n;i++) for(int j=0;j<n;j++) Vm[i][m] += A[i][j]*Vm[j][k];
			for(int i=0;i<=k;i++) {
				H_trans[k][i] = 0.0;
				for(int j=0;j<n;j++) H_trans[k][i] += Vm[j][m]*Vm[j][i];
			}
			for(int i=0;i<n;i++) {
				Vm[i][k+1] = Vm[i][m];
				for(int j=0;j<=k;j++) Vm[i][k+1] -= H_trans[k][j]*Vm[i][j];
			}
			for(int i=0;i<n;i++) r0[i] = Vm[i][k+1];
			H_trans[k][k+1] = vecnorm2(r0,n);
			if(k==(m-1) && H_trans[m-1][m]==0.0) for(int i=0;i<n;i++) Vm[i][m] = 0.0;
			else for(int i=0;i<n;i++) Vm[i][k+1] /= H_trans[k][k+1];
		}
		//Least squares
		e1[0] = err;
		for(int i=1;i<=m;i++) e1[i] = 0.0;
		for(int i=0;i<m;i++) {
			c = H_trans[i][i]/sqrt(H_trans[i][i]*H_trans[i][i]+H_trans[i][i+1]*H_trans[i][i+1]);
			s = H_trans[i][i+1]/sqrt(H_trans[i][i]*H_trans[i][i]+H_trans[i][i+1]*H_trans[i][i+1]);
			for(int j=i;j<m;j++) {
				tmpi = H_trans[j][i];
				tmpj = H_trans[j][i+1];
				H_trans[j][i] = c*tmpi + s*tmpj;
				H_trans[j][i+1] = c*tmpj - s*tmpi;
			}
			tmpi = e1[i];
			tmpj = e1[i+1];
			e1[i] = c*tmpi + s*tmpj;
			e1[i+1] = c*tmpj - s*tmpi;
		}
		y[m-1] = e1[m-1]/H_trans[m-1][m-1];
		for(int i=m-2;i>=0;i--) {
			for(int j=m-1;j>i;j--) e1[i] -= H_trans[j][i]*y[j];
			y[i] = e1[i]/H_trans[i][i];
		}
		//Update x
		for(int i=0;i<n;i++) for(int j=0;j<m;j++) x[i] += Vm[i][j]*y[j];
		for(int i=0;i<n;i++) {
			r0[i] = b[i];
			for(int j=0;j<n;j++) r0[i] -= A[i][j]*x[j];
		}
		err = fabs(e1[m]);
		itrerr[itcounter] = err;
		itcounter++;
	}
	delete[] r0;
	for(int i=0;i<n;i++) delete[] Vm[i];
	delete[] Vm;
	for(int i=0;i<m;i++) delete[] H_trans[i];
	delete[] H_trans;
	delete[] e1;
	delete[] y;
#ifdef CHINESE_VERSION
	cout << "\n求解完成。" << endl;
#else
	cout << "\nFinish solving." << endl;
#endif
}

void GMRES::out_result() {
#ifdef CHINESE_VERSION
	cout << "\n通过GMRES求解。" << endl;
	fl << "\n通过GMRES求解。\n";
#else
	cout << "\nSolved by GMRES." << endl;
	fl << "\nSolved by GMRES.\n";
#endif
	Iteration_method::out_result();
}

GMRES::~GMRES() {
	delete_Ax_b();
	delete[]itrerr;
}
