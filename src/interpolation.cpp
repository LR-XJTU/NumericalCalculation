#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "interpolation.h"
#include "common_fd.h"
#include "directmethod.h"
#include "int_diff.h"

using namespace std;

//initialization
void Interpolation::init() {
	fl.init("Interpolation.txt");
	int cc;
	cout.setf(ios::left);
	cout << "\nWant to input a continuous function or a list funtion? (1 = continuous , 2 = list)" << endl;
	input_fx = false;
	cc = in_int();
	if (cc == 1) {
		fx.init();
		input_fx = true;
	}
}

//input the interpolation nodes
void Interpolation::input_data() {
	double temp;
	if (input_fx) {
		int tf;
		cout << "\nWant to input arbitrary x or get uniform sampling? (1 = arbitrary x , 2 = uniform sampling)" << endl;
		tf = in_int();
		if (tf == 1) {
			cout << "\nPlease input x of the points (x,f(x)):" << endl;
			for (int i = 0; i < np1; i++)
			{
				cout << "x = ";
				cin >> temp;
				in_x(i, temp);
				in_y(i, fx.f(temp));
			}
		}
		else if (tf == 2) {
			double a, b, h;
			cout << "\nPlease input the interval [a , b]:\na = ";
			cin >> a;
			cout << "b = ";
			cin >> b;
			h = (b - a) / (np1 - 1);
			for (int i = 0; i < np1; i++) {
				in_x(i, a + i * h);
				in_y(i, fx.f(a + i * h));
			}
		}
	}
	else {
		cout << "\nPlease input the points (x,y):" << endl;
		cout << "x\ty" << endl;
		for (int i = 0; i < np1 * 2; i++)
		{
			cin >> temp;
			if (i % 2 == 0) in_x(i / 2, temp);
			else in_y((i - 1) / 2, temp);
		}
	}
	cout << "\nFinish inputting." << endl;
	sort_x();
	cout_xy();
}

void Interpolation::cout_xy() {
	cout << "\nThe points are:" << endl;
	for (int i = 0; i < np1; i++) {
		cout << "(" << x[i] << " , " << y[i] << ")" << "\n";
	}
}

void Interpolation::save_xy() {
	fl << "\nThe points are:\n";
	for (int i = 0; i < np1; i++) {
		fl << "(" << x[i] << " , " << y[i] << ")\n";
	}
}

//sort the points by x
void Interpolation::sort_x() {
	double t;
	for (int i = 0; i < np1 - 1; i++) {
		for (int j = 0; j < np1 - 1 - i; j++) {
			if (x[j] > x[j + 1]) {
				t = x[j + 1];
				x[j + 1] = x[j];
				x[j] = t;
				t = y[j + 1];
				y[j + 1] = y[j];
				y[j] = t;
			}
		}
	}
}

void Interpolation::in_x(int i, double t) {
	x[i] = t;
}

void Interpolation::in_y(int i, double t) {
	y[i] = t;
}

//cout a polynomial
void Interpolation::cout_poly(double* p, int l) {
	bool zero_flag = 1;
	for (int i = 0; i < l; i++) {
		if (p[i] == 0.0 && (i < l - 1 || !zero_flag)) continue;
		if (p[i] > 0.0 && !zero_flag) cout << " + ";
		else if (p[i] < 0.0) cout << (zero_flag ? "- " : " - ");
		if (fabs(p[i]) != 1.0 || i == l - 1) cout << fabs(p[i]);
		if (i < l - 2) cout << "x^" << l - 1 - i;
		else if (i == l - 2) cout << "x";
		if (p[i] != 0.0) zero_flag = 0;
	}
}

//save a polynomial
void Interpolation::save_poly(double* p, int l) {
	bool zero_flag = 1;
	for (int i = 0; i < l; i++) {
		if (p[i] == 0.0 && (i < l - 1 || !zero_flag)) continue;
		if (p[i] > 0.0 && !zero_flag) fl << " + ";
		else if (p[i] < 0.0) fl << (zero_flag ? "- " : " - ");
		if (fabs(p[i]) != 1.0 || i == l - 1) fl << fabs(p[i]);
		if (fabs(p[i]) != 1.0 && i < l - 1) fl << "*";
		if (i < l - 2) fl << "x^" << l - 1 - i;
		else if (i == l - 2) fl << "x";
		if (p[i] != 0.0) zero_flag = 0;
	}
}

//transform a polynomial into string
std::string Interpolation::polytostr(double* p, int l) {
	std::string ts;
	bool zero_flag = 1;
	for (int i = 0; i < l; i++) {
		if (p[i] == 0.0 && (i < l - 1 || !zero_flag)) continue;
		if (p[i] > 0.0 && !zero_flag) ts += "+";
		else if (p[i] < 0.0) ts += "-";
		if (fabs(p[i]) != 1.0 || i == l - 1) {
			stringstream ss;
			ss << setprecision(15) << fabs(p[i]);
			ts += ss.str();
		}
		if (fabs(p[i]) != 1.0 && i < l - 1) ts += ".*";
		if (i < l - 2) {
			ts += "x.^";
			ts += to_string(l - 1 - i);
		}
		else if (i == l - 2) ts += "x";
		if (p[i] != 0.0) zero_flag = 0;
	}
	return ts;
}

//1.Newton interpolation

Newton_Ip::Newton_Ip() {
	init();
	cout << "\nPlease enter the number of points:" << endl;
	int nn = in_int();
	if (nn <= 0)
	{
		cout << "Error: Please enter a positive integer." << endl;
		throw 0;
	}
	np1 = nn;
	x = new double[np1];
	y = new double[np1];
	f = new double* [np1];
	for (int i = 0; i < np1; i++) f[i] = new double[np1];
	N = new double[np1];
	for (int i = 0; i < np1; i++) N[i] = 0.0;
	input_data();
}

Newton_Ip::Newton_Ip(bool no_init) {
	if (no_init) return;
	else Newton_Ip();
}

Newton_Ip::Newton_Ip(double* xx, double* yy, double nnp1) {
	np1 = nnp1;
	x = new double[np1];
	y = new double[np1];
	for (int i = 0; i < np1; i++) {
		x[i] = xx[i];
		y[i] = yy[i];
	}
	f = new double* [np1];
	for (int i = 0; i < np1; i++) f[i] = new double[np1];
	N = new double[np1];
	for (int i = 0; i < np1; i++) N[i] = 0.0;
}

void Newton_Ip::cout_polynomial() {
	cout << "\nN = ";
	cout_poly(N, np1);
	cout << endl;
}

void Newton_Ip::save_polynomial() {
	fl << "\nN = ";
	save_poly(N, np1);
	fl << "\n";
}

void Newton_Ip::calc() {
	cout << "\nCalculating the interpolation polynomial..." << endl;
	for (int j = 0; j < np1; j++) {
		for (int i = j; i < np1; i++) {
			difference_quotient(i, j);
		}
	}
	for (int i = 0; i < np1; i++) {
		pi_expansion(i);
	}
	cout_polynomial();
	int add_p = 1;
	while (1) {
		cout << "\nWant to add a point? (1 = yes , 0 = no)" << endl;
		cin >> add_p;
		if (!add_p) break;
		double xx, yy;
		cout << "\nPlease enter the new point:\nx = ";
		cin >> xx;
		if (input_fx) yy = fx.f(xx);
		else {
			cout << "y = ";
			cin >> yy;
		}
		add_point(xx, yy);
		for (int j = 0; j < np1; j++) difference_quotient(np1 - 1, j);
		pi_expansion(np1 - 1);
		cout_polynomial();
	}
	cout << "\nFinish interpolating." << endl;
}

void Newton_Ip::pure_calc() {
	for (int j = 0; j < np1; j++) {
		for (int i = j; i < np1; i++) {
			difference_quotient(i, j);
		}
	}
	for (int i = 0; i < np1; i++) {
		pi_expansion(i);
	}
}

//transform a polynomial into string
std::string Newton_Ip::polytostr_nodot(double* p, int l) {
	std::string ts;
	bool zero_flag = 1;
	for (int i = 0; i < l; i++) {
		if (p[i] == 0.0 && (i < l - 1 || !zero_flag)) continue;
		if (p[i] > 0.0 && !zero_flag) ts += "+";
		else if (p[i] < 0.0) ts += "-";
		if (fabs(p[i]) != 1.0 || i == l - 1) {
			stringstream ss;
			ss << setprecision(15) << fabs(p[i]);
			ts += ss.str();
		}
		if (fabs(p[i]) != 1.0 && i < l - 1) ts += "*";
		if (i < l - 2) {
			ts += "x^";
			ts += to_string(l - 1 - i);
		}
		else if (i == l - 2) ts += "x";
		if (p[i] != 0.0) zero_flag = 0;
	}
	return ts;
}

string Newton_Ip::get_str() {
	string s;
	s = polytostr_nodot(N, np1);
	return s;
}

void Newton_Ip::out_result() {
	sort_x();
	cout_xy();
	save_xy();
	cout << "\nThe Newton interpolation function is" << endl;
	fl << "\nThe Newton interpolation function is\n";
	cout_polynomial();
	save_polynomial();
	generate_m();
}

void Newton_Ip::generate_m() {
	cout << "\nWant to generate .m file to get figure of iteration error in MATLAB? (1 = Yes , 0 = No)" << endl;
	int flag = in_int();
	if (flag == 1) {
		filelog gm;
		gm.init("Interpolation.m");
		gm << "\nfigure\n";
		gm << "x = [";
		for (int i = 0; i < np1; i++) gm << x[i] << " ";
		gm << "];\n";
		gm << "y = [";
		for (int i = 0; i < np1; i++) gm << y[i] << " ";
		gm << "];\n";
		gm << "plot(x,y,'x')\n";
		gm << "legend('Sample points')\n";
		gm << "hold on\n";
		if (input_fx) {
			string s = fx.get_fstr();
			size_t an = 0;
			while (s.find("*", an) != s.npos) {
				an = s.find("*", an);
				s.replace(an, 1, ".*");
				an += 2;
				if (an >= s.size()) break;
			}
			an = 0;
			while (s.find("/", an) != s.npos) {
				an = s.find("/", an);
				s.replace(an, 1, "./");
				an += 2;
				if (an >= s.size()) break;
			}
			an = 0;
			while (s.find("^", an) != s.npos) {
				an = s.find("^", an);
				s.replace(an, 1, ".^");
				an += 2;
				if (an >= s.size()) break;
			}
			gm << "fplot(@(x) " << s << ",[" << x[0] << "," << x[np1 - 1] << "])\n" << "hold on\n";
		}
		gm << "fplot(@(x) " << polytostr(N, np1) << ",[" << x[0] << "," << x[np1 - 1] << "])\n";
		gm << "xlabel('x');\n";
		gm << "ylabel('y');\n";
		gm << "title('The graph of interpolation polynomial');\n";
	}
}

void Newton_Ip::difference_quotient(int i, int j) {
	if (j == 0) f[i][j] = y[i];
	else f[i][j] = (f[i][j - 1] - f[i - 1][j - 1]) / (x[i] - x[i - j]);
}

void Newton_Ip::pi_expansion(int ord) {
	int maxnum = int_pow(2, ord);
	int cnt, temp;
	double p;
	for (int i = 0; i < maxnum; i++) {
		p = 1.0;
		cnt = 0;
		for (int j = 0; j < ord; j++) {
			temp = int_pow(2, j);
			if (i & temp) {
				cnt++;
				p *= x[j];
			}
		}
		p *= int_pow(-1, cnt) * f[ord][ord];
		N[np1 - 1 - ord + cnt] += p;
	}
}

void Newton_Ip::add_point(double xx, double yy) {
	double* temp = new double[np1];
	for (int i = 0; i < np1; i++) temp[i] = x[i];
	delete[] x;
	np1++;
	x = new double[np1];
	for (int i = 0; i < np1 - 1; i++) x[i] = temp[i];
	x[np1 - 1] = xx;
	np1--;
	for (int i = 0; i < np1; i++) temp[i] = y[i];
	delete[] y;
	np1++;
	y = new double[np1];
	for (int i = 0; i < np1 - 1; i++) y[i] = temp[i];
	y[np1 - 1] = yy;
	np1--;
	for (int i = 0; i < np1; i++) temp[i] = N[i];
	delete[] N;
	np1++;
	N = new double[np1];
	N[0] = 0.0;
	for (int i = 1; i < np1; i++) N[i] = temp[i-1];
	delete[] temp;
	np1--;
	double** ff = new double* [np1];
	for (int i = 0; i < np1; i++) ff[i] = new double[np1];
	for (int i = 0; i < np1; i++) for (int j = 0; j <= i; j++) ff[i][j] = f[i][j];
	for (int i = 0; i < np1; i++) delete[] f[i];
	delete[] f;
	np1++;
	f = new double* [np1];
	for (int i = 0; i < np1; i++) f[i] = new double[np1];
	for (int i = 0; i < np1 - 1; i++) for (int j = 0; j <= i; j++) f[i][j] = ff[i][j];
	for (int i = 0; i < np1 - 1; i++) delete[] ff[i];
	delete[] ff;
}

Newton_Ip::~Newton_Ip() {
	delete[] x;
	delete[] y;
	for (int i = 0; i < np1; i++) delete[] f[i];
	delete[] f;
	delete[] N;
}

//2.Hermite interpolation

Hermite_Ip::Hermite_Ip() :Newton_Ip(true) {
	init();
	cout << "\nPlease enter the number of points: (including the multiple nodes)" << endl;
	int nn = in_int();
	if (nn <= 0)
	{
		cout << "Error: Please enter a positive integer." << endl;
		throw 0;
	}
	np1 = nn;
	x = new double[np1];
	y = new double[np1];
	f = new double* [np1];
	for (int i = 0; i < np1; i++) f[i] = new double[np1];
	N = new double[np1];
	for (int i = 0; i < np1; i++) N[i] = 0.0;
	d = new int[np1];
	for (int i = 0; i < np1; i++) d[i] = 0;
	input_data();
}

void Hermite_Ip::input_data() {
	double temp;
	int t;
	if (input_fx) {
		DerivExtra Deriv(fx);
		cout << "\nPlease input x and the orders of derivatives:" << endl;
		cout << "x\torder of derivative (0 represents the function value)" << endl;
		for (int i = 0; i < np1; i++)
		{
			cin >> temp;
			t = in_int();
			in_x(i, temp);
			in_d(i, t);
			if (t == 0) {
				in_y(i, fx.f(temp));
			}
			else {
				temp = Deriv.orderderiv(t, temp);
				in_y(i, temp);
			}
		}
	}
	else {
		cout << "\nPlease input the points and the orders of derivatives:" << endl;
		cout << "x\ty\torder of derivative (0 represents the function value)" << endl;
		for (int i = 0; i < np1 * 3; i++)
		{
			cin >> temp;
			if (i % 3 == 0) in_x(i / 3, temp);
			else if (i % 3 == 1) in_y(i / 3, temp);
			else in_d(i / 3, (int)temp);
		}
	}
	cout << "\nFinish inputting." << endl;
	sort_x();
	cout_xy();
}

void Hermite_Ip::cout_xy() {
	int t = 0;
	double tt = x[0];
	cout << "\nThe points are:" << endl;
	for (int i = 0; i < np1; i++) {
		if (i > 0 && tt == x[i]) t++;
		else tt = x[i];
		cout << "(x" << i - t << ",f";
		for (int j = 0; j < d[i]; j++) cout << "\47";
		cout << "(x" << i - t << ")) = ";
		cout << "(" << x[i] << " , " << y[i] << ")" << "\n";
	}
}

void Hermite_Ip::save_xy() {
	int t = 0;
	double tt = x[0];
	fl << "\nThe points are:\n";
	for (int i = 0; i < np1; i++) {
		if (i > 0 && tt == x[i]) t++;
		else tt = x[i];
		fl << "(x" << i - t << ",f";
		for (int j = 0; j < d[i]; j++) fl << "\47";
		fl << "(x" << i - t << ")) = ";
		fl << "(" << x[i] << " , " << y[i] << ")\n";
	}
}

void Hermite_Ip::sort_x() {
	double t;
	int tt;
	for (int i = 0; i < np1 - 1; i++) {
		for (int j = 0; j < np1 - 1 - i; j++) {
			if (x[j] > x[j + 1]) {
				t = x[j + 1];
				x[j + 1] = x[j];
				x[j] = t;
				t = y[j + 1];
				y[j + 1] = y[j];
				y[j] = t;
				tt = d[j + 1];
				d[j + 1] = d[j];
				d[j] = tt;
			}
		}
	}
	int cnt = 1;
	t = x[0];
	for (int i = 1; i < np1; i++) {
		if (x[i] == t) cnt++;
		else {
			if (cnt > 1) {
				for (int k = 0; k < cnt - 1; k++) {
					for (int j = i - cnt; j < i - k - 1; j++) {
						if (d[j] > d[j + 1]) {
							t = x[j + 1];
							x[j + 1] = x[j];
							x[j] = t;
							t = y[j + 1];
							y[j + 1] = y[j];
							y[j] = t;
							tt = d[j + 1];
							d[j + 1] = d[j];
							d[j] = tt;
						}
					}
				}
				cnt = 1;
			}
			t = x[i];
		}
	}
	if (cnt > 1) {
		for (int k = 0; k < cnt - 1; k++) {
			for (int j = np1 - cnt; j < np1 - k - 1; j++) {
				if (d[j] > d[j + 1]) {
					t = x[j + 1];
					x[j + 1] = x[j];
					x[j] = t;
					t = y[j + 1];
					y[j + 1] = y[j];
					y[j] = t;
					tt = d[j + 1];
					d[j + 1] = d[j];
					d[j] = tt;
				}
			}
		}
	}
}

void Hermite_Ip::in_d(int i, int t) {
	d[i] = t;
}

void Hermite_Ip::cout_polynomial() {
	cout << "\nH = ";
	cout_poly(N, np1);
	cout << endl;
}

void Hermite_Ip::save_polynomial() {
	fl << "\nH = ";
	save_poly(N, np1);
	fl << "\n";
}

void Hermite_Ip::calc() {
	cout << "\nCalculating the interpolation polynomial..." << endl;
	int add_p = 1;
	do {
		for (int j = 0; j < np1; j++) {
			for (int i = j; i < np1; i++) {
				difference_quotient(i, j);
			}
		}
		for (int i = 0; i < np1; i++) {
			pi_expansion(i);
		}
		cout_polynomial();
		cout << "\nWant to add a point? (1 = yes , 0 = no)" << endl;
		cin >> add_p;
		if (!add_p) break;
		double xx, yy;
		int dd;
		cout << "\nPlease enter the new point:\nx = ";
		cin >> xx;
		cout << "y = ";
		cin >> yy;
		cout << "order of derivative = ";
		cin >> dd;
		add_point(xx, yy, dd);
		sort_x();
	} while (add_p);
	cout << "\nFinish interpolating." << endl;
}

void Hermite_Ip::out_result() {
	cout_xy();
	save_xy();
	cout << "\nThe Hermite interpolation function is" << endl;
	fl << "\nThe Hermite interpolation function is\n";
	cout_polynomial();
	save_polynomial();
	generate_m();
}

void Hermite_Ip::difference_quotient(int i, int j) {
	if (j == 0) {
		if (d[i] > 0) {
			for (int k = 0; k < np1; k++) if (x[k] == x[i] && d[k] == 0) {
				f[i][j] = y[k];
				break;
			}
		}
		else f[i][j] = y[i];
	}
	else if (d[i] >= j) {
		for (int k = 0; k < np1; k++) if (x[k] == x[i] && d[k] == j) {
			f[i][j] = y[k] / calc_fac(j);
			break;
		}
	}
	else f[i][j] = (f[i][j - 1] - f[i - 1][j - 1]) / (x[i] - x[i - j]);
}

void Hermite_Ip::add_point(double xx, double yy, int dd) {
	int* temp = new int[np1];
	for (int i = 0; i < np1; i++) temp[i] = d[i];
	delete[] d;
	np1++;
	d = new int[np1];
	for (int i = 0; i < np1 - 1; i++) d[i] = temp[i];
	d[np1 - 1] = dd;
	delete[] temp;
	np1--;
	Newton_Ip::add_point(xx, yy);
	for (int i = 0; i < np1; i++) N[i] = 0.0;
}

Hermite_Ip::~Hermite_Ip() {
	delete[] x;
	delete[] y;
	for (int i = 0; i < np1; i++) delete[] f[i];
	delete[] f;
	delete[] N;
	delete[] d;
}

//3.cube spline interpolation

cube_spline::cube_spline() {
	init();
	cout << "\nPlease enter the number of points:" << endl;
	int nn = in_int();
	if (nn <= 0)
	{
		cout << "Error: Please enter a positive integer." << endl;
		throw 0;
	}
	np1 = nn;
	x = new double[np1];
	y = new double[np1];
	S = new double* [np1 - 1];
	for (int i = 0; i < np1 - 1; i++) S[i] = new double[4];
	input_data();
	input_bd();
}

void cube_spline::input_bd() {
	cout << "\nPlease select the type of boundary condition:" << endl;
	cout << "1 = Second derivatives of interval boundary are given." << endl;
	cout << "2 = First derivatives of interval boundary are given." << endl;
	cout << "3 = Periodic boundary condition" << endl;
	cin >> bdcd_flag;
	if (bdcd_flag < 3) {
		cout << "\nPlease input the boundary condition:" << endl;
	}
	if (bdcd_flag == 1) {
		cout << "The left end point: f\42(x_0) = ";
		cin >> bd[0];
		cout << "The right end point: f\42(x_n) = ";
		cin >> bd[1];
	}
	else if (bdcd_flag == 2) {
		cout << "The left end point: f\47(x_0) = ";
		cin >> bd[0];
		cout << "The right end point: f\47(x_n) = ";
		cin >> bd[1];
	}
}

void cube_spline::cout_polynomial() {
	cout << "\nS = ";
	for (int i = 0; i < np1 - 1; i++) {
		cout << "\t";
		cout_poly(S[i], 4);
		cout << " , " << x[i] << " <= x <= " << x[i + 1] << endl;
	}
}

void cube_spline::save_polynomial() {
	fl << "\nS = ";
	for (int i = 0; i < np1 - 1; i++) {
		fl << "\t";
		save_poly(S[i], 4);
		fl << " , " << x[i] << " <= x <= " << x[i + 1] << "\n";
	}
}

void cube_spline::calc() {
	cout << "\nCalculating the interpolation polynomial..." << endl;
	double* h = new double[np1 - 1];
	for (int i = 0; i < np1 - 1; i++) h[i] = x[i + 1] - x[i];
	double* M_solve, * M, * mu, * lambda, * d;
	double extra_lm[2];
	int siz;
	switch (bdcd_flag) {
	case 1:
		siz = np1 - 2;
		break;
	case 2:
		siz = np1;
		break;
	case 3:
		siz = np1 - 1;
		break;
	default:
		break;
	}
	M = new double[np1];
	M_solve = new double[siz];
	d = new double[siz];
	lambda = new double[siz - 1];
	mu = new double[siz - 1];
	switch (bdcd_flag) {
	case 1:
		for (int i = 0; i < siz - 1; i++) lambda[i] = h[i + 1] / (h[i] + h[i + 1]);
		extra_lm[0] = h[siz] / (h[siz - 1] + h[siz]);
		extra_lm[1] = h[0] / (h[0] + h[1]);
		for (int i = 0; i < siz - 1; i++) mu[i] = h[i + 1] / (h[i + 1] + h[i + 2]);
		for (int i = 0; i < siz; i++) d[i] = 6.0 / (h[i] + h[i + 1]) * ((y[i + 2] - y[i + 1]) / h[i + 1] - (y[i + 1] - y[i]) / h[i]);
		d[0] -= extra_lm[1] * bd[0];
		d[siz - 1] -= extra_lm[0] * bd[1];
		break;
	case 2:
		lambda[0] = 1.0;
		for (int i = 1; i < siz - 1; i++) lambda[i] = h[i] / (h[i - 1] + h[i]);
		for (int i = 0; i < siz - 2; i++) mu[i] = h[i] / (h[i] + h[i + 1]);
		mu[siz - 2] = 1.0;
		d[0] = 6.0 / h[0] * ((y[1] - y[0]) / h[0] - bd[0]);
		for (int i = 1; i < siz - 1; i++) d[i] = 6.0 / (h[i - 1] + h[i]) * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
		d[siz - 1] = 6.0 / h[np1 - 2] * (bd[1] - (y[np1 - 1] - y[np1 - 2]) / h[np1 - 2]);
		break;
	case 3:
		for (int i = 0; i < siz - 1; i++) lambda[i] = h[i + 1] / (h[i] + h[i + 1]);
		extra_lm[0] = h[0] / (h[np1 - 2] + h[0]);
		extra_lm[1] = h[0] / (h[0] + h[1]);
		for (int i = 0; i < siz - 2; i++) mu[i] = h[i + 1] / (h[i + 1] + h[i + 2]);
		mu[siz - 2] = h[np1 - 2] / (h[np1 - 2] + h[0]);
		for (int i = 0; i < siz - 1; i++) d[i] = 6.0 / (h[i] + h[i + 1]) * ((y[i + 2] - y[i + 1]) / h[i + 1] - (y[i + 1] - y[i]) / h[i]);
		d[siz - 1] = 6.0 / (h[np1 - 2] + h[0]) * ((y[1] - y[np1 - 1]) / h[0] - (y[np1 - 1] - y[np1 - 2]) / h[np1 - 2]);
		break;
	default:
		break;
	}
	double* a = new double[siz];
	for (int i = 0; i < siz; i++) a[i] = 2.0;
	Direct_method* mat;
	if (bdcd_flag == 3) mat = new Doolittle(true);
	else mat = new Chasing(true);
	mat->init(siz, siz, a, lambda, mu, d);
	if (bdcd_flag == 3) {
		mat->in_A(siz - 1, 0, extra_lm[0]);
		mat->in_A(0, siz - 1, extra_lm[1]);
	}
	mat->calc();
	mat->get_x(M_solve);
	delete mat;
	switch (bdcd_flag) {
	case 1:
		M[0] = bd[0];
		for (int i = 1; i < np1 - 1; i++) M[i] = M_solve[i - 1];
		M[np1 - 1] = bd[1];
		break;
	case 2:
		for (int i = 0; i < np1; i++) M[i] = M_solve[i - 1];
		break;
	case 3:
		M[0] = M_solve[np1 - 2];
		for (int i = 1; i < np1; i++) M[i] = M_solve[i - 1];
		break;
	default:
		break;
	}
	for (int i = 0; i < np1 - 1; i++) S[i][0] = (M[i + 1] - M[i]) / 6.0 / h[i];
	for (int i = 0; i < np1 - 1; i++) S[i][1] = (M[i] * x[i + 1] - M[i + 1] * x[i]) / 2.0 / h[i];
	for (int i = 0; i < np1 - 1; i++) S[i][2] = ((3.0 * sqr(x[i]) - sqr(h[i])) * M[i + 1] - (3.0 * sqr(x[i + 1]) - sqr(h[i])) * M[i] + 6.0 * (y[i + 1] - y[i])) / 6.0 / h[i];
	for (int i = 0; i < np1 - 1; i++) S[i][3] = ((sqr(h[i]) - sqr(x[i])) * x[i] * M[i + 1] - (sqr(h[i]) - sqr(x[i + 1])) * x[i + 1] * M[i] + 6.0 * (x[i + 1] * y[i] - x[i] * y[i + 1])) / 6.0 / h[i];
	delete[] h;
	delete[] M;
	delete[] M_solve;
	delete[] mu;
	delete[] lambda;
	delete[] d;
	cout << "\nFinish interpolating." << endl;
}

void cube_spline::out_result() {
	cout_xy();
	save_xy();
	cout << "\nThe cubic spline interpolation function is" << endl;
	fl << "\nThe cubic spline interpolation function is\n";
	cout_polynomial();
	save_polynomial();
	generate_m();
}

void cube_spline::generate_m() {
	cout << "\nWant to generate .m file to get the graph of interpolation in MATLAB? (1 = Yes , 0 = No)" << endl;
	int flag = in_int();
	if (flag == 1) {
		filelog gm;
		gm.init("Interpolation.m");
		gm << "\nfigure\n";
		gm << "x = [";
		for (int i = 0; i < np1; i++) gm << x[i] << " ";
		gm << "];\n";
		gm << "y = [";
		for (int i = 0; i < np1; i++) gm << y[i] << " ";
		gm << "];\n";
		gm << "plot(x,y,'x')\n";
		gm << "legend('Sample points')\n";
		gm << "hold on\n";
		if (input_fx) {
			string s = fx.get_fstr();
			size_t an = 0;
			while (s.find("*", an) != s.npos) {
				an = s.find("*", an);
				s.replace(an, 1, ".*");
				an += 2;
				if (an >= s.size()) break;
			}
			an = 0;
			while (s.find("/", an) != s.npos) {
				an = s.find("/", an);
				s.replace(an, 1, "./");
				an += 2;
				if (an >= s.size()) break;
			}
			an = 0;
			while (s.find("^", an) != s.npos) {
				an = s.find("^", an);
				s.replace(an, 1, ".^");
				an += 2;
				if (an >= s.size()) break;
			}
			gm << "fplot(@(x) " << s << ",[" << x[0] << "," << x[np1 - 1] << "])\n" << "hold on\n";
		}
		for (int i = 0; i < np1 - 1; i++) {
			gm << "fplot(@(x) " << polytostr(S[i], 4) << ",[" << x[i] << "," << x[i + 1] << "])\n";
			if (i < np1 - 2) gm << "hold on\n";
		}
		gm << "xlabel('x');\n";
		gm << "ylabel('y');\n";
		gm << "title('The graph of cubic spline interpolation function');\n";
	}
}

cube_spline::~cube_spline() {
	init();
	delete[] x;
	delete[] y;
	for (int i = 0; i < np1 - 1; i++) delete[] S[i];
	delete[] S;
}
