#include <iostream>
#include <string>
#include "common_fd.h"
#include "formula.h"
#include "directmethod.h"
#include "iterationmethod.h"
#include "int_diff.h"
#include "interpolation.h"
#include "optimal_approx.h"
#include "nonlinear_equations.h"
#include "eigen_val_vec.h"

using namespace std;

void Numerical_calculating() {
	int chapter;
#ifdef CHINESE_VERSION
	cout << "\n--- 数值计算 ---" << endl;
	cout << "请选择功能：" << endl;
	cout << "1 直接法求解线性方程组" << endl;
	cout << "2 迭代法求解线性方程组" << endl;
	cout << "3 插值法" << endl;
	cout << "4 函数最优逼近" << endl;
	cout << "5 数值积分与数值微分" << endl;
	cout << "6 非线性方程（组）" << endl;
	cout << "7 矩阵特征值与特征向量" << endl;
	cout << "附加功能:" << endl;
	cout << "0 科学计算器" << endl;
#else
	cout << "\n--- Numerical calculating ---" << endl;
	cout << "Please select the chapter:" << endl;
	cout << "1.Direct method for solving linear equations" << endl;
	cout << "2.Iterative method for solving linear equations" << endl;
	cout << "3.Interpolation" << endl;
	cout << "4.Optimal approximation of function" << endl;
	cout << "5.Numerical integration and differentiation" << endl;
	cout << "6.Nonlinear equation(s)" << endl;
	cout << "7.Eigenvalue and eigenvector of matrix" << endl;
	cout << "Addtional function:" << endl;
	cout << "0.Scientific calculator" << endl;
#endif
	chapter = in_int();
	if (chapter < 0 || chapter >7)
	{
#ifdef CHINESE_VERSION
		cout << "错误：超出索引" << endl;
#else
		cout << "Error: Out of index range." << endl;
#endif
		return;
	}
	switch (chapter) {
		case 0: {
#ifdef CHINESE_VERSION
			cout << "\n--- 科学计算器 ---" << endl;
			cout << "* 需要退出时请输入000" << endl;
#else
			cout << "\n--- Scientific calculator ---" << endl;
			cout << "* Input 000 if you want to exit." << endl;
#endif
			formula fml;
			string str;
			double x;
			fml.showtips();
			cout << "\nf(x) = ";
			cin >> str;
			if (!str.compare("000")) break;
			fml.init(str);
			do {
				if (fml.x_flag()) {
					cout << "\nf(x) = " << fml.get_fstr();
#ifdef CHINESE_VERSION
					cout << "   (Matlab格式)" << endl;
#else
					cout << "   (Matlab format)" << endl;
#endif
					cout << "  x  = ";
					cin >> str;
					if (!str.compare("000")) break;
					x = StringToNum<double>(str);
					x = fml.f(x);
					cout << "f(" << str << ") = " << x << endl;
				}
				else {
					x = fml.f(0.0);
					cout << "     = " << x << endl;
					cout << "\nf(x) = ";
					cin >> str;
					if (!str.compare("000")) break;
					fml.init(str);
				}
			} while (1);
			break;
		}
		case 1: {
			Direct_method* dm = NULL;
			int method;
#ifdef CHINESE_VERSION
			cout << "\n--- 求解线性方程组 Ax=b ---" << endl;
			cout << "请选择求解方法：" << endl;
			cout << "1 高斯消去法" << endl;
			cout << "2 列主元高斯消去法" << endl;
			cout << "3 杜利特尔分解（LU分解）" << endl;
			cout << "4 楚列斯基分解" << endl;
			cout << "5 改进平方根法" << endl;
			cout << "6 追赶法" << endl;
			cout << "7 吉文斯变换（QR分解）" << endl;
			cout << "8 豪斯霍尔德变换" << endl;
#else
			cout << "\n--- To solve linear equations Ax=b ---" << endl;
			cout << "Please select the solution method:" << endl;
			cout << "1.Gauss elimination" << endl;
			cout << "2.Column principle Gaussian elimination" << endl;
			cout << "3.Doolittle decomposition (LU)" << endl;
			cout << "4.Cholesky decomposition" << endl;
			cout << "5.Improved square root" << endl;
			cout << "6.Chasing" << endl;
			cout << "7.Givens transformation (QR)" << endl;
			cout << "8.Householder transformation" << endl;
#endif
			
			method = in_int();
			if (method < 1 || method >8)
			{
#ifdef CHINESE_VERSION
				cout << "错误：超出索引" << endl;
#else
				cout << "Error: Out of index range." << endl;
#endif
				return;
			}
			switch (method) {
				case 1:
					dm = new Gauss;
					break;
				case 2:
					dm = new CP_Gauss;
					break;
				case 3:
					dm = new Doolittle;
					break;
				case 4:
					dm = new Cholesky;
					break;
				case 5:
					dm = new Improved_sqrt;
					break;
				case 6:
					dm = new Chasing;
					break;
				case 7:
					dm = new Givens;
					break;
				case 8:
					dm = new Householder;
					break;
				default:
					break;
			}
			dm->calc();
			dm->out_result();
			delete dm;
			break;
		}
		case 2: {
			Iteration_method* im = NULL;
			int method;
#ifdef CHINESE_VERSION
			cout << "\n--- 求解线性方程组 Ax=b ---" << endl;
			cout << "请选择求解方法：" << endl;
			cout << "1 雅克比迭代法" << endl;
			cout << "2 高斯-赛德尔迭代法" << endl;
			cout << "3 逐次超松弛迭代法" << endl;
			cout << "4 最速下降法" << endl;
			cout << "5 共轭梯度法" << endl;
			cout << "6 广义极小残余算法" << endl;
#else
			cout << "\n--- To solve linear equations Ax=b ---" << endl;
			cout << "Please select the solution method:" << endl;
			cout << "1.Jacobi" << endl;
			cout << "2.Gauss-Seidel" << endl;
			cout << "3.SOR" << endl;
			cout << "4.Steepest descent method" << endl;
			cout << "5.Conjugate gradient method" << endl;
			cout << "6.GMRES" << endl;
#endif
			method = in_int();
			if (method < 1 || method >6)
			{
#ifdef CHINESE_VERSION
				cout << "错误：超出索引" << endl;
#else
				cout << "Error: Out of index range." << endl;
#endif
				return;
			}
			switch (method)
			{
				case 1:
					im = new Jacobi;
					break;
				case 2:
					im = new Gauss_Seidel;
					break;
				case 3:
					im = new SOR;
					break;
				case 4:
					im = new steepest_descent;
					break;
				case 5:
					im = new conjugate_gradient;
					break;
				case 6:
					im = new GMRES;
				default:
					break;
			}
			im->calc();
			im->out_result();
			delete im;
			break;
		}
		case 3: {
			Interpolation* ipn = NULL;
			int method;
#ifdef CHINESE_VERSION
			cout << "\n--- 插值法 ---" << endl;
			cout << "请选择插值方法：" << endl;
			cout << "1 牛顿插值" << endl;
			cout << "2 埃尔米特插值（牛顿型）" << endl;
			cout << "3 三次样条插值" << endl;
#else
			cout << "\n--- Interpolation ---" << endl;
			cout << "Please select the interpolation method:" << endl;
			cout << "1.Newton interpolation" << endl;
			cout << "2.Hermite interpolation (Newton type)" << endl;
			cout << "3.Cubic spline interpolation" << endl;
#endif
			method = in_int();
			if (method < 1 || method >3)
			{
#ifdef CHINESE_VERSION
				cout << "错误：超出索引" << endl;
#else
				cout << "Error: Out of index range." << endl;
#endif
				return;
			}
			switch (method)
			{
				case 1:
					ipn = new Newton_Ip;
					break;
				case 2:
					ipn = new Hermite_Ip;
					break;
				case 3:
					ipn = new cube_spline;
					break;
				default:
					break;
			}
			ipn->calc();
			ipn->out_result();
			delete ipn;
			break;
		}
		case 4: {
			optimal_approx* oa = NULL;
			int method;
#ifdef CHINESE_VERSION
			cout << "\n--- 函数最优逼近 ---" << endl;
			cout << "请选择逼近方法：" << endl;
			cout << "1 最优平方逼近" << endl;
			cout << "2 近似最优一致逼近" << endl;
#else
			cout << "\n--- Optimal approximation of function ---" << endl;
			cout << "Please select the method:" << endl;
			cout << "1.Optimal square approximation" << endl;
			cout << "2.Approximate optimal uniform approximation" << endl;
#endif
			method = in_int();
			if (method < 1 || method>2) {
#ifdef CHINESE_VERSION
				cout << "错误：超出索引" << endl;
#else
				cout << "Error: Out of index range." << endl;
#endif
				return;
			}
			switch (method) {
			case 1:
				oa = new sqr_approx;
				break;
			case 2:
				oa = new uni_approx;
				break;
			default:
				break;
			}
			oa->calc();
			oa->out_result();
			delete oa;
			break;
		}
		case 5: {
			Int_Diff* ig = NULL;
			int method;
#ifdef CHINESE_VERSION
			cout << "\n--- 数值积分与数值微分 ---" << endl;
			cout << "请选择积分/微分方法：" << endl;
			cout << "1 龙贝格积分" << endl;
			cout << "2 外推法求导" << endl;
#else
			cout << "\n--- Numerical integration and differentiation ---" << endl;
			cout << "Please select the method:" << endl;
			cout << "1.Romberg integration" << endl;
			cout << "2.Derivation by extrapolation" << endl;
#endif
			method = in_int();
			if (method < 1 || method >2)
			{
#ifdef CHINESE_VERSION
				cout << "错误：超出索引" << endl;
#else
				cout << "Error: Out of index range." << endl;
#endif
				return;
			}
			switch (method)
			{
			case 1:
				ig = new Romberg;
				break;
			case 2:
				ig = new DerivExtra;
				break;
			default:
				break;
			}
			ig->calc();
			ig->out_result();
			delete ig;
			break;
		}
		case 6: {
			nl_eq* neq = NULL;
			nl_eqs* neqs = NULL;
			int method;
#ifdef CHINESE_VERSION
			cout << "\n--- 非线性方程（组） ---" << endl;
			cout << "请选择求解方法：" << endl;
			cout << "--- 非线性方程 ---" << endl;
			cout << "1 简单迭代法" << endl;
			cout << "2 牛顿法" << endl;
			cout << "3 弦割法" << endl;
			cout << "4 松弛加速法" << endl;
			cout << "5 艾特肯加速法（斯特芬森）" << endl;
			cout << "--- 非线性方程组 ---" << endl;
			cout << "6 简单迭代法" << endl;
			cout << "7 牛顿法" << endl;
			cout << "8 弦割法" << endl;
			cout << "9 布洛伊登法" << endl;
#else
			cout << "\n--- Nonlinear equation(s) ---" << endl;
			cout << "Please select the method:" << endl;
			cout << "--- Nonlinear equation ---" << endl;
			cout << "1.Simple iteration method" << endl;
			cout << "2.Newton method" << endl;
			cout << "3.Secant method" << endl;
			cout << "4.Relaxation acceleration method" << endl;
			cout << "5.Aitken acceleration method (Steffensen)" << endl;
			cout << "--- Nonlinear equations ---" << endl;
			cout << "6.Simple iteration method" << endl;
			cout << "7.Newton method" << endl;
			cout << "8.Secant method" << endl;
			cout << "9.Broyden method" << endl;
#endif
			method = in_int();
			if (method < 1 || method>9) {
#ifdef CHINESE_VERSION
				cout << "错误：超出索引" << endl;
#else
				cout << "Error: Out of index range." << endl;
#endif
				return;
			}
			switch (method) {
			case 1:
				neq = new eq_Simple;
				break;
			case 2:
				neq = new eq_Newton;
				break;
			case 3:
				neq = new eq_Secant;
				break;
			case 4:
				neq = new eq_Relaxation;
				break;
			case 5:
				neq = new eq_Aitken;
				break;
			case 6:
				neqs = new eqs_Simple;
				break;
			case 7:
				neqs = new eqs_Newton;
				break;
			case 8:
				neqs = new eqs_Secant;
				break;
			case 9:
				neqs = new eqs_Broyden;
				break;
			default:
				break;
			}
			if (method < 6) {
				neq->calc();
				neq->out_result();
				delete neq;
			}
			else {
				neqs->calc();
				neqs->out_result();
				delete neqs;
			}
			break;
		}
		case 7: {
			Eigen_val_vec* evv = NULL;
			int method;
#ifdef CHINESE_VERSION
			cout << "\n--- 矩阵特征值与特征向量 ---" << endl;
			cout << "1 乘幂法" << endl;
			cout << "2 反幂法" << endl;
#else
			cout << "\n--- Eigenvalue and eigenvector of matrix ---" << endl;
			cout << "1.Power method" << endl;
			cout << "2.Inverse power method" << endl;
#endif
			method = in_int();
			if (method < 1 || method>2) {
#ifdef CHINESE_VERSION
				cout << "错误：超出索引" << endl;
#else
				cout << "Error: Out of index range." << endl;
#endif
				return;
			}
			switch (method) {
			case 1:
				evv = new Power_method;
				break;
			case 2:
				evv = new Inverse_power;
				break;
			default:
				break;
			}
			evv->calc();
			evv->out_result();
			delete evv;
			break;
		}
		default:
			break;
		}
}

int main() {
	int c_flag = 1;
	do{
		cout.precision(15);
		try {
			Numerical_calculating();
		}
		catch (int err) {
#ifdef CHINESE_VERSION
			if (err == 0) cout << "\n程序运行中遇到了一处错误" << endl;
			else if (err == 1) cout << "\n输入数据出错" << endl;
			cout << "\n是否重新开始？（1 = 是 , 0 = 否）" << endl;
#else
			if (err == 0) cout << "\nAn error occurred in the program running." << endl;
			else if (err == 1) cout << "\nAn error occurred about the given data." << endl;
			cout << "\nWant to restart? (1 = yes , 0 = no)" << endl;
#endif
			c_flag = in_int();
			continue;
		}
#ifdef CHINESE_VERSION
		cout << "\n是否继续？（1 = 是 , 0 = 否）" << endl;
#else
		cout << "\nWant to continue? (1 = yes , 0 = no)" << endl;
#endif
		c_flag = in_int();
	}while(c_flag);
#ifdef CHINESE_VERSION
	cout << "\n感谢使用" << endl;
#else
	cout << "\nThanks for using." << endl;
#endif
	return 0;
}
