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
	cout << "\n--- ��ֵ���� ---" << endl;
	cout << "��ѡ���ܣ�" << endl;
	cout << "1 ֱ�ӷ�������Է�����" << endl;
	cout << "2 ������������Է�����" << endl;
	cout << "3 ��ֵ��" << endl;
	cout << "4 �������űƽ�" << endl;
	cout << "5 ��ֵ��������ֵ΢��" << endl;
	cout << "6 �����Է��̣��飩" << endl;
	cout << "7 ��������ֵ����������" << endl;
	cout << "���ӹ���:" << endl;
	cout << "0 ��ѧ������" << endl;
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
		cout << "���󣺳�������" << endl;
#else
		cout << "Error: Out of index range." << endl;
#endif
		return;
	}
	switch (chapter) {
		case 0: {
#ifdef CHINESE_VERSION
			cout << "\n--- ��ѧ������ ---" << endl;
			cout << "* ��Ҫ�˳�ʱ������000" << endl;
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
					cout << "   (Matlab��ʽ)" << endl;
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
			cout << "\n--- ������Է����� Ax=b ---" << endl;
			cout << "��ѡ����ⷽ����" << endl;
			cout << "1 ��˹��ȥ��" << endl;
			cout << "2 ����Ԫ��˹��ȥ��" << endl;
			cout << "3 �����ض��ֽ⣨LU�ֽ⣩" << endl;
			cout << "4 ����˹���ֽ�" << endl;
			cout << "5 �Ľ�ƽ������" << endl;
			cout << "6 ׷�Ϸ�" << endl;
			cout << "7 ����˹�任��QR�ֽ⣩" << endl;
			cout << "8 ��˹�����±任" << endl;
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
				cout << "���󣺳�������" << endl;
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
			cout << "\n--- ������Է����� Ax=b ---" << endl;
			cout << "��ѡ����ⷽ����" << endl;
			cout << "1 �ſ˱ȵ�����" << endl;
			cout << "2 ��˹-���¶�������" << endl;
			cout << "3 ��γ��ɳڵ�����" << endl;
			cout << "4 �����½���" << endl;
			cout << "5 �����ݶȷ�" << endl;
			cout << "6 ���弫С�����㷨" << endl;
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
				cout << "���󣺳�������" << endl;
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
			cout << "\n--- ��ֵ�� ---" << endl;
			cout << "��ѡ���ֵ������" << endl;
			cout << "1 ţ�ٲ�ֵ" << endl;
			cout << "2 �������ز�ֵ��ţ���ͣ�" << endl;
			cout << "3 ����������ֵ" << endl;
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
				cout << "���󣺳�������" << endl;
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
			cout << "\n--- �������űƽ� ---" << endl;
			cout << "��ѡ��ƽ�������" << endl;
			cout << "1 ����ƽ���ƽ�" << endl;
			cout << "2 ��������һ�±ƽ�" << endl;
#else
			cout << "\n--- Optimal approximation of function ---" << endl;
			cout << "Please select the method:" << endl;
			cout << "1.Optimal square approximation" << endl;
			cout << "2.Approximate optimal uniform approximation" << endl;
#endif
			method = in_int();
			if (method < 1 || method>2) {
#ifdef CHINESE_VERSION
				cout << "���󣺳�������" << endl;
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
			cout << "\n--- ��ֵ��������ֵ΢�� ---" << endl;
			cout << "��ѡ�����/΢�ַ�����" << endl;
			cout << "1 ���������" << endl;
			cout << "2 ���Ʒ���" << endl;
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
				cout << "���󣺳�������" << endl;
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
			cout << "\n--- �����Է��̣��飩 ---" << endl;
			cout << "��ѡ����ⷽ����" << endl;
			cout << "--- �����Է��� ---" << endl;
			cout << "1 �򵥵�����" << endl;
			cout << "2 ţ�ٷ�" << endl;
			cout << "3 �Ҹ" << endl;
			cout << "4 �ɳڼ��ٷ�" << endl;
			cout << "5 ���ؿϼ��ٷ���˹�ط�ɭ��" << endl;
			cout << "--- �����Է����� ---" << endl;
			cout << "6 �򵥵�����" << endl;
			cout << "7 ţ�ٷ�" << endl;
			cout << "8 �Ҹ" << endl;
			cout << "9 �������Ƿ�" << endl;
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
				cout << "���󣺳�������" << endl;
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
			cout << "\n--- ��������ֵ���������� ---" << endl;
			cout << "1 ���ݷ�" << endl;
			cout << "2 ���ݷ�" << endl;
#else
			cout << "\n--- Eigenvalue and eigenvector of matrix ---" << endl;
			cout << "1.Power method" << endl;
			cout << "2.Inverse power method" << endl;
#endif
			method = in_int();
			if (method < 1 || method>2) {
#ifdef CHINESE_VERSION
				cout << "���󣺳�������" << endl;
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
			if (err == 0) cout << "\n����������������һ������" << endl;
			else if (err == 1) cout << "\n�������ݳ���" << endl;
			cout << "\n�Ƿ����¿�ʼ����1 = �� , 0 = ��" << endl;
#else
			if (err == 0) cout << "\nAn error occurred in the program running." << endl;
			else if (err == 1) cout << "\nAn error occurred about the given data." << endl;
			cout << "\nWant to restart? (1 = yes , 0 = no)" << endl;
#endif
			c_flag = in_int();
			continue;
		}
#ifdef CHINESE_VERSION
		cout << "\n�Ƿ��������1 = �� , 0 = ��" << endl;
#else
		cout << "\nWant to continue? (1 = yes , 0 = no)" << endl;
#endif
		c_flag = in_int();
	}while(c_flag);
#ifdef CHINESE_VERSION
	cout << "\n��лʹ��" << endl;
#else
	cout << "\nThanks for using." << endl;
#endif
	return 0;
}
