#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "optimal_approx.h"

#include "int_diff.h"
#include "iterationmethod.h"
#include "interpolation.h"

using namespace std;

optimal_approx::optimal_approx() {
	err = 0.0;
	fl.init("Optimal_approx.txt");
}

double optimal_approx::inner_product(formula f, formula g, formula w, double a, double b) {
	formula fg = w * f * g;
	Romberg integ(fg, a, b);
	integ.calc();
	return integ.get_result();
}

double optimal_approx::inner_product(formula f, formula g, double* x, int m, double* w) {
	double* fy = new double[m];
	double* gy = new double[m];
	for (int i = 0; i < m; i++) {
		fy[i] = f.f(x[i]);
		gy[i] = g.f(x[i]);
	}
	double r = 0.0;
	for (int i = 0; i < m; i++) r += w[i] * fy[i] * gy[i];
	delete [] fy, gy;
	return r;
}

double optimal_approx::inner_product(formula f, double* y, double* x, int m, double* w) {
	double* fy = new double[m];
	for (int i = 0; i < m; i++) fy[i] = f.f(x[i]);
	double r = 0.0;
	for (int i = 0; i < m; i++) r += w[i] * y[i] * fy[i];
	delete[] fy;
	return r;
}

double optimal_approx::inner_product(double* y1, double* y2, int m, double* w) {
	double r = 0.0;
	for (int i = 0; i < m; i++) r += w[i] * y1[i] * y2[i];
	return r;
}

void optimal_approx::Cheby_poly(double** c, int nnp1) {
	c[0][0] = 1.0;
	if (nnp1 == 1) return;
	c[1][0] = 1.0;
	c[1][1] = 0.0;
	for (int i = 2; i < nnp1; i++) {
		for (int j = 0; j < i + 1; j++) c[i][j] = 0.0;
		for (int j = 0; j < i; j++) c[i][j] += 2.0 * c[i - 1][j];
		for (int j = 2; j < i + 1; j++) c[i][j] -= c[i - 2][j - 2];
	}
}

sqr_approx::sqr_approx() {
	cout << "\nWant to approximate a continuous function or a list function? (1 = continuous, 2 = list)" << endl;
	int flag = in_int();
	if (flag == 1) cont = true;
	else if (flag == 2) cont = false;
	else {
		cout << "Error: Illegal input." << endl;
		throw 0;
	}
	if (cont) {
		fx.init();
		cout << "\nPlease input the interval [a,b]:" << endl;
		cout << "a = ";
		cin >> interval[0];
		cout << "b = ";
		cin >> interval[1];
	}
	else {
		cout << "\nPlease enter the number of points: m = ";
		m = in_int();
		x = new double[m];
		y = new double[m];
		wl = new double[m];
		cout << "\nPlease input (x,y):\nx\ty" << endl;
		for (int i = 0; i < m; i++) {
			cin >> x[i] >> y[i];
		}
		sort_xy();
		interval[0] = x[0];
		interval[1] = x[m - 1];
	}
	cout << "\nThe weight w uses the default value 1? (1 = Yes, 0 = No)" << endl;
	int wf = in_int();
	if (wf == 1) {
		w_flag = false;
		if (cont) {
			string s = "1";
			wc.init(s);
		}
		else {
			for (int i = 0; i < m; i++) {
				wl[i] = 1.0;
			}
		}
	}
	else if (wf == 0) {
		w_flag = true;
		if (cont) {
			cout << "\nPlease input weight function w:\nw = ";
			string s;
			cin >> s;
			wc.init(s);
		}
		else {
			cout << "\nPlease enter w_i: (i=0,1,...," << m << ")" << endl;
			for (int i = 0; i < m; i++) cin >> wl[i];
		}
	}
	else {
		cout << "Error: Illegal input." << endl;
		throw 0;
	}
	cout << "\nPlease select the basis functions:" << endl;
	cout << "1.Recursive orthogonal polynomials" << endl;
	cout << "2.self-defined basis functions" << endl;
	basis_flag = in_int();
	if (basis_flag == 1) {
		cout << "\nPlease enter the order of the polynomial: n = ";
		np1 = in_int() + 1;
		basis = new formula[np1];
		coef = new double[np1];
		for (int i = 0; i < np1; i++) coef[i] = 0.0;
		c = new double* [np1];
		orth_poly(basis, np1);
	}
	else if (basis_flag == 2) {
		cout << "\nPlease enter the number of basis functions: n = ";
		np1 = in_int();
		basis = new formula[np1];
		coef = new double[np1];
		for (int i = 0; i < np1; i++) coef[i] = 0.0;
		string s;
		for (int i = 0; i < np1; i++) {
			cout << "f" << i << "(x) = ";
			cin >> s;
			basis[i].init(s);
		}
		cout << "\nThe basis functions are:" << endl;
		for (int i = 0; i < np1; i++) cout << "g" << i << "(x) = " << basis[i].get_fstr() << endl;
		fl << "\nThe basis functions are:\n";
		for (int i = 0; i < np1; i++) fl << "\ng" << i << "(x) = " << basis[i].get_fstr();
		fl << "\n";
	}
	else {
		cout << "Error: Illegal input." << endl;
		throw 0;
	}
}

sqr_approx::~sqr_approx() {
	if (!cont) {
		delete[] x;
		delete[] y;
		delete[] wl;
	}
	delete[] basis;
	delete[] coef;
	for (int i = 0; i < np1; i++) delete[] c[i];
	delete[] c;
}

//sort the points by x
void sqr_approx::sort_xy() {
	double t;
	for (int i = 0; i < m - 1; i++) {
		for (int j = 0; j < m - 1 - i; j++) {
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

double sqr_approx::innpro(formula f, formula g) {
	double r;
	if (cont) r = inner_product(f, g, wc, interval[0], interval[1]);
	else r = inner_product(f, g, x, m, wl);
	return r;
}

//transform a polynomial into string
string optimal_approx::polytostr(double* p, int l) {
	string ts;
	bool zero_flag = true;
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
		if (p[i] != 0.0) zero_flag = false;
	}
	return ts;
}

void sqr_approx::orth_poly(formula* ff, int n) {
	string s;
	formula xx;
	s = "x";
	xx.init(s);
	c[0] = new double[1];
	c[0][0] = 1.0;
	s = polytostr(c[0], 1);
	ff[0].init(s);
	if (n == 1) return;
	double bk, ck, beta, gamma;
	beta = innpro(xx * ff[0], ff[0]);
	gamma = innpro(ff[0], ff[0]);
	bk = beta / gamma;
	c[1] = new double[2];
	c[1][0] = 1.0;
	c[1][1] = -1.0 * bk;
	s = polytostr(c[1], 2);
	ff[1].init(s);
	for (int i = 2; i < n; i++) {
		ck = 1.0 / gamma;
		beta = innpro(xx * ff[i-1], ff[i-1]);
		gamma = innpro(ff[i - 1], ff[i - 1]);
		bk = beta / gamma;
		ck *= gamma;
		c[i] = new double[i + 1];
		for (int j = 0; j < i + 1; j++) c[i][j] = 0.0;
		for (int j = 0; j < i; j++) c[i][j] += c[i - 1][j];
		for (int j = 1; j < i + 1; j++) c[i][j] += -1.0 * bk * c[i - 1][j - 1];
		for (int j = 2; j < i + 1; j++) c[i][j] += -1.0 * ck * c[i - 2][j - 2];
		s = polytostr(c[i], i + 1);
		ff[i].init(s);
	}
	cout << "\nThe basis functions are:" << endl;
	for (int i = 0; i < n; i++) cout << "g" << i << "(x) = " << ff[i].get_fstr() << endl;
}

void sqr_approx::calc() {
	cout << "\nCalculating optimal square approximation function..." << endl;
	double* cc = new double[np1];
	if (basis_flag == 1) {
		if (cont) {
			for (int i = 0; i < np1; i++) cc[i] = innpro(basis[i], fx) / innpro(basis[i], basis[i]);
		}
		else {
			for (int i = 0; i < np1; i++) cc[i] = inner_product(basis[i], y, x, m, wl) / innpro(basis[i], basis[i]);
		}
		for (int i = 0; i < np1; i++) {
			for (int j = np1 - i - 1; j < np1; j++) coef[j] += cc[i] * c[i][j - np1 + i + 1];
		}
		resultstr = polytostr(coef, np1);
	}
	else if (basis_flag == 2) {
		double** A = new double* [np1];
		for (int i = 0; i < np1; i++) A[i] = new double[np1];
		double* b = new double[np1];
		for (int i = 0; i < np1; i++) {
			for (int j = 0; j < np1; j++) {
				A[i][j] = innpro(basis[i], basis[j]);
			}
		}
		if (cont) {
			for (int i = 0; i < np1; i++) b[i] = innpro(basis[i], fx);
		}
		else {
			for (int i = 0; i < np1; i++) b[i] = inner_product(basis[i], y, x, m, wl);
		}
		conjugate_gradient conj(A, b, np1);
		conj.calc();
		conj.get_result(cc);
		for (int i = 0; i < np1; i++) delete[] A[i];
		delete[] A;
		delete[] b;
		for (int i = 0; i < np1; i++) coef[i] = cc[i];
		stringstream ss1;
		ss1 << setprecision(15) << coef[0];
		resultstr = ss1.str();
		resultstr += "*" + basis[0].get_fstr();
		for (int i = 1; i < np1; i++) {
			resultstr += (coef[i] >= 0) ? "+" : "-";
			stringstream ss2;
			ss2 << setprecision(15) << fabs(coef[i]);
			resultstr += ss2.str() + "*" + basis[i].get_fstr();
		}
	}
	if (cont) {
		for (int i = 0; i < np1; i++) {
			for (int j = 0; j < np1; j++) {
				err += cc[i] * cc[j] * innpro(basis[i], basis[j]);
			}
			err -= 2.0 * cc[i] * innpro(basis[i], fx);
		}
		err += innpro(fx, fx);
	}
	else {
		for (int i = 0; i < np1; i++) {
			for (int j = 0; j < np1; j++) {
				err += cc[i] * cc[j] * innpro(basis[i], basis[j]);
			}
			err -= 2.0 * cc[i] * inner_product(basis[i], y, x, m, wl);
		}
		err += inner_product(y, y, m, wl);
	}
	err = sqrt(err);
	delete [] cc;
	cout << "\nFinish calculating." << endl;
}

void sqr_approx::out_result() {
	if (cont) {
		cout << "\nf(x) = " << fx.get_fstr() << " , " << interval[0] << " <= x <= " << interval[1] << endl;
		fl << "\nf(x) = " << fx.get_fstr() << " , " << interval[0] << " <= x <= " << interval[1] << "\n";
		cout << "w(x) = " << wc.get_fstr() << endl;
		fl << "w(x) = " << wc.get_fstr() << "\n";
		fl << "\nThe basis functions are:\n";
		for (int i = 0; i < np1; i++) fl << "\ng" << i << "(x) = " << basis[i].get_fstr();
		fl << "\n";
		cout << "\nThe optimal square approximation of f(x) on [" << interval[0] << "," << interval[1] << "] is\n\np(x) = ";
		fl << "\nThe optimal square approximation of f(x) on [" << interval[0] << "," << interval[1] << "] is\n\np(x) = ";
	}
	else {
		cout << "\nThe sample points and the weight coefficients are:" << endl;
		fl << "\nThe sample points and the weight coefficients are:\n";
		for (int i = 0; i < m; i++) {
			cout << "(" << x[i] << " , " << y[i] << ") , w_" << i + 1 << " = " << wl[i] << endl;
			fl << "(" << x[i] << " , " << y[i] << ") , w_" << i + 1 << " = " << wl[i] << "\n";
		}
		fl << "\nThe basis functions are:\n";
		for (int i = 0; i < np1; i++) fl << "\ng" << i << "(x) = " << basis[i].get_fstr();
		fl << "\n";
		cout << "\nThe optimal square approximation of the sample points is\n\np(x) = ";
		fl << "\nThe optimal square approximation of the sample points is\n\np(x) = ";
	}
	cout << resultstr << endl;
	fl << resultstr << "\n";
	cout << "\nThe mean square error = " << err << endl;
	fl << "\nThe mean square error = " << err << "\n";
	generate_m();
}

void sqr_approx::generate_m() {
	cout << "\nWant to generate .m file to get figure of the optimal approximation function in MATLAB? (1 = Yes , 0 = No)" << endl;
	int flag = in_int();
	if (flag == 1) {
		string s;
		size_t an;
		filelog gm;
		gm.init("Optimal_approx.m");
		gm << "\nfigure\n";
		if (cont) {
			s = fx.get_fstr();
			an = 0;
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
			gm << "fplot(@(x) " << s << ",[" << interval[0] << "," << interval[1] << "])\n";
			gm << "legend; \n";
			gm << "hold on\n";
		}
		else {
			gm << "x = [";
			for (int i = 0; i < m; i++) gm << x[i] << " ";
			gm << "];\n";
			gm << "y = [";
			for (int i = 0; i < m; i++) gm << y[i] << " ";
			gm << "];\n";
			gm << "plot(x,y,'x')\n";
			gm << "legend('Sample points')\n";
			gm << "hold on\n";
		}

		s = resultstr;
		an = 0;
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
		gm << "fplot(@(x) " << s << ",[" << interval[0] << "," << interval[1] << "])\n";
		gm << "xlabel('x');\n";
		gm << "ylabel('y');\n";
		gm << "title('The graph of the optimal square approximation');\n";
	}
}

uni_approx::uni_approx() {
	cout << "\nPlease select the method of approximate optimal uniform approximation:" << endl;
	cout << "1.Chebyshev interpolation polynomials" << endl;
	cout << "2.Truncated Chebyshev series method" << endl;
	ua_method = in_int();
	if (ua_method < 1 || ua_method>2) {
		cout << "Error: Out of index range." << endl;
		throw 0;
	}
	string temp_fx;
	if (ua_method == 1) {
		fx.init();
	}
	else if (ua_method == 2) {
		fx.showtips();
		cout << "\nf(x) = ";
		cin >> temp_fx;
		fx.init(temp_fx);
	}
	cout << "\nPlease input the interval [a,b]:" << endl;
	cout << "a = ";
	cin >> interval[0];
	cout << "b = ";
	cin >> interval[1];
	cout << "\nPlease enter the order of the polynomial of approximate optimal uniform approximation: n = ";
	np1 = in_int() + 1;
	coef = new double[np1];
	for (int i = 0; i < np1; i++) coef[i] = 0.0;
	if (ua_method == 2) {
		double a1, a2;
		a1 = (interval[1] - interval[0]) / 2.0;
		a2 = (interval[0] + interval[1]) / 2.0;
		stringstream ss1, ss2;
		ss1 << setprecision(15) << a1;
		ss2 << setprecision(15) << fabs(a2);
		string s = "(" + ss1.str() + "x"; //x=f(t)
		s += (a2 >= 0) ? "+" : "-";
		s += ss2.str() + ")";
		string cosx = "cos ( x )"; //t=cos(theta)
		formula xt;
		xt.init(s, cosx);
		fx_xt_cos.init(temp_fx, xt.get_fstr());
	}
}

uni_approx::~uni_approx() {
	delete [] coef;
}

void uni_approx::calc() {
	if (ua_method == 1) Cheby_interp();
	else if (ua_method == 2) Trunc_Cheby();
}

void uni_approx::out_result() {
	cout << "\nf(x) = " << fx.get_fstr() << " , " << interval[0] << " <= x <= " << interval[1] << endl;
	fl << "\nf(x) = " << fx.get_fstr() << " , " << interval[0] << " <= x <= " << interval[1] << "\n";
	cout << "The optimal uniform approximation of f(x) on [" << interval[0] << "," << interval[1] << "] is\n\np(x) = ";
	fl << "The optimal uniform approximation of f(x) on [" << interval[0] << "," << interval[1] << "] is\n\np(x) = ";
	cout << resultstr << endl;
	fl << resultstr << "\n";
	cout << "\nSolved by ";
	fl << "\nSolved by ";
	if (ua_method == 1) {
		cout << "Chebyshev interpolation polynomials." << endl;
		fl << "Chebyshev interpolation polynomials.\n";
	}
	else if (ua_method == 2) {
		cout << "truncated Chebyshev series method." << endl;
		fl << "truncated Chebyshev series method.\n";
	}
	generate_m();
}

void uni_approx::generate_m() {
	cout << "\nWant to generate .m file to get figure of the optimal approximation function in MATLAB? (1 = Yes , 0 = No)" << endl;
	int flag = in_int();
	if (flag == 1) {
		string s;
		size_t an;
		filelog gm;
		gm.init("Optimal_approx.m");
		gm << "\nfigure\n";
		s = fx.get_fstr();
		an = 0;
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
		gm << "fplot(@(x) " << s << ",[" << interval[0] << "," << interval[1] << "])\n" << "hold on\n";

		s = resultstr;
		an = 0;
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
		gm << "fplot(@(x) " << s << ",[" << interval[0] << "," << interval[1] << "])\n";
		gm << "legend;\n";
		gm << "xlabel('x');\n";
		gm << "ylabel('y');\n";
		gm << "title('The graph of the optimal uniform approximation');\n";
	}
}

void uni_approx::Cheby_interp() {
	double* x = new double[np1];
	double* y = new double[np1];
	for (int i = 0; i < np1; i++) {
		x[i] = cos((2.0 * i + 1.0) / 2.0 / np1 * 3.14159265358979323846);
		x[i] = (interval[1] - interval[0]) / 2.0 * x[i] + (interval[1] + interval[0]) / 2.0;
		y[i] = fx.f(x[i]);
	}
	Newton_Ip newton(x, y, np1);
	newton.pure_calc();
	resultstr = newton.get_str();
	delete[]x;
	delete[]y;
}

void uni_approx::Trunc_Cheby() {
	c = new double* [np1]; //coefficients of Chebyshev polynomials
	for (int i = 0; i < np1; i++) c[i] = new double[i + 1];
	Cheby_poly(c, np1);

	double* cc = new double[np1]; //coefficients of Chebyshev series
	Romberg rg(fx_xt_cos, 0.0, 3.14159265358979323846);
	rg.calc();
	cc[0] = rg.get_result() / 3.14159265358979323846;
	for (int i = 1; i < np1; i++) {
		string s = "cos(" + to_string(i) + "*x)";
		formula coskx;
		coskx.init(s);
		Romberg rgt(fx_xt_cos * coskx, 0.0, 3.14159265358979323846);
		rgt.calc();
		cc[i] = rgt.get_result() * 2.0 / 3.14159265358979323846;
	}

	double ctx[2]; //t=f^(-1)(x)
	ctx[0] = 2.0 / (interval[1] - interval[0]);
	ctx[1] = -1.0 * (interval[1] + interval[0]) / (interval[1] - interval[0]);
	double** yanghui_tri = new double* [np1]; //Yanghui triangle
	yanghui_tri[0] = new double[1];
	yanghui_tri[0][0] = 1.0;
	if (np1 > 1) {
		yanghui_tri[1] = new double[2];
		yanghui_tri[1][0] = 1.0;
		yanghui_tri[1][1] = 1.0;
	}
	if (np1 > 2) {
		for (int i = 2; i < np1; i++) {
			yanghui_tri[i] = new double[i + 1];
			yanghui_tri[i][0] = 1.0;
			yanghui_tri[i][i] = 1.0;
			for (int j = 1; j < i; j++) yanghui_tri[i][j] = yanghui_tri[i - 1][j - 1] + yanghui_tri[i - 1][j];
		}
	}
	double** ccc = new double* [np1]; //coefficients of Chebyshev polynomials when the variable is x
	for (int i = 0; i < np1; i++) ccc[i] = new double[i + 1];
	for (int i = 0; i < np1; i++) for (int j = 0; j < i + 1; j++) ccc[i][j] = 0.0;
	for (int i = 0; i < np1; i++) {
		for (int j = 0; j < i + 1; j++) {
			for (int k = j; k < i + 1; k++) {
				ccc[i][k] += c[i][j] * yanghui_tri[i - j][k - j] * pow(ctx[0], i - k) * pow(ctx[1], k - j);
			}
		}
	}
	for (int i = 0; i < np1; i++) {
		for (int j = np1 - i - 1; j < np1; j++) coef[j] += cc[i] * ccc[i][j - np1 + i + 1];
	}

	resultstr = polytostr(coef, np1);
	for (int i = 0; i < np1; i++) delete [] c[i];
	delete [] c;
	delete [] cc;
	delete [] ccc;
	for (int i = 0; i < np1; i++) delete[] yanghui_tri[i];
	delete[] yanghui_tri;
}