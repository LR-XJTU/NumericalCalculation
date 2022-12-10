#include <iostream>
#include <string>
#include <cmath>
#include "common_fd.h"
#include "directmethod.h"
#include "nonlinear_equations.h"

using namespace std;

nl_eq::nl_eq() {
	fl.init("Nonlinear_equations.txt");
}

void nl_eq::input_fx() {
	cout << "\nNonlinear equation: f(x) = 0" << endl;
	fx.showtips();
	cout << "\nf(x) = ";
	string s;
	cin >> s;
	fx.init(s);
	cout << "\nThe equation is " << fx.get_fstr() << " = 0 (MATLAB format)" << endl;
	fl << "\nThe equation is " << fx.get_fstr() << " = 0 (MATLAB format)\n";
}

void nl_eq::input_x() {
	cout << "\nPlease enter the initial value x0 or an interval [a,b]: (1 = x0, 2 = [a,b]) ";
	int flag = in_int();
	if (flag == 1) {
		double x0;
		cout << "x0 = ";
		cin >> x0;
		x.push_back(x0);
		count_flag = true;
	}
	else if (flag == 2) {
		double a, b;
		cout << "a = ";
		cin >> a;
		cout << "b = ";
		cin >> b;
		interval[0] = a;
		interval[1] = b;
		isolation();
		count_flag = false;
	}
	else {
		cout << "Error: Illegal input." << endl;
		throw 0;
	}
}

void nl_eq::input_max_eps() {
	cout << "\nSet the errors and the maximum iteration times:" << endl;
	cout << "  |f(x_k+1)| < eps1 , or |x_k+1 - x_k| < eps2" << endl;
	cout << "  Default values: eps1 = 1.0e-6, eps2 = 1.0e-6, maximum iteration times = 10000" << endl;
	cout << "Use default values? (1 = Yes, 0 = No) ";
	int flag = in_int();
	if (flag == 1) {
		eps1 = 1.0e-6;
		eps2 = 1.0e-6;
		maxcounter = 10000;
	}
	else if (flag == 0) {
		cout << "eps1 = ";
		cin >> eps1;
		cout << "eps2 = ";
		cin >> eps2;
		cout << "maximum iteration times = ";
		maxcounter = in_int();
	}
	else {
		cout << "Error: Illegal input." << endl;
		throw 0;
	}
	if (eps1 < 1.0) {
		int pre;
		cout.setf(iostream::fixed);
		pre = -1 * floor(log10(eps1));
		cout.precision(pre);
	}
	fl.set_precision(eps1);
	fl << "\neps1 = " << eps1 << ", eps2 = " << eps2 << " (|f(x_k+1)| < eps1 , or |x_k+1 - x_k| < eps2)";
	fl << "\nmaximum iteration times = " << maxcounter << "\n";
}

void nl_eq::isolation() {
	double xk = interval[0], h = 0.01 * (interval[1] - interval[0]);
	cout << "\nNote: The distance of two adjacent roots should be above " << h << ", otherwise one of them may not be detected." << endl;
	double t1 = sign(fx.f(xk)), t2;
	if (t1 == 0.0) x.push_back(interval[0]);
	for (int i = 0; i < 100; i++) {
		xk += h;
		t2 = sign(fx.f(xk));
		if (t2 == 0.0) x.push_back(xk);
		else {
			if (t1 * t2 < 0.0) x.push_back(bisection(xk - h, xk));
		}
		t1 = t2;
	}
	cout << "Found " << x.size() << " roots on the interval." << endl;
	for (int i = 0; i < x.size(); i++) cout << "x = " << x[i] << endl;
}

double nl_eq::bisection(double a, double b) {
	double t, ft, fa, fb;
	fa = sign(fx.f(a));
	fb = sign(fx.f(b));
	if (fa == 0.0) return a;
	if (fb == 0.0) return b;
	for (int i = 0; i < 10; i++) {
		t = 0.5 * (a + b);
		ft = sign(fx.f(t));
		if (ft == 0.0) return t;
		if (ft * fa < 0) b = t;
		else if (ft * fb < 0) a = t;
		else cout << "\nWarning: An exception happened in bisection method. There may exist no root or more than one root in the interval." << endl;
	}
	t = 0.5 * (a + b);
	return t;
}

eq_Simple::eq_Simple() {
	input_fx();
	input_x();
	input_max_eps();
}

void eq_Simple::input_fx() {
	cout << "\nEquivalent transformation of nonlinear equation: x = g(x)" << endl;
	fx.showtips();
	cout << "\ng(x) = ";
	cin >> str;
	fx.init(str);
	cout << "\nThe equivalent transformation is x = " << fx.get_fstr() << " (MATLAB format)" << endl;
	fl << "\nThe equivalent transformation is x = " << fx.get_fstr() << " (MATLAB format)\n";
}

void eq_Simple::input_x() {
	cout << "\nPlease enter the initial value x0 or an interval [a,b]: (1 = x0 , 2 = [a,b]) ";
	int flag = in_int();
	if (flag == 1) {
		double x0;
		cout << "x0 = ";
		cin >> x0;
		x.push_back(x0);
		count_flag = true;
	}
	else if (flag == 2) {
		double a, b;
		cout << "a = ";
		cin >> a;
		cout << "b = ";
		cin >> b;
		interval[0] = a;
		interval[1] = b;
		string s = str + "-x";
		fx.init(s);
		isolation();
		fx.init(str);
		count_flag = false;
	}
	else {
		cout << "Error: Illegal input." << endl;
		throw 0;
	}
}

void eq_Simple::calc() {
	for (int i = 0; i < x.size(); i++) {
		counter = 0;
		double xk;
		while (counter < maxcounter) {
			xk = x[i];
			x[i] = fx.f(xk);
			if (fabs(fx.f(x[i]) - x[i]) < eps1) break;
			if (fabs(x[i] - xk) < eps2) break;
			counter++;
		}
		if (counter == maxcounter) cout << "\nWarning: Reach the maximum iteration number when calculate x = " << x[i] << endl;
	}
}

void eq_Simple::out_result() {
	cout << "\nCalculated by simple iteration method." << endl;
	fl << "\nCalculated by simple iteration method.";
	if (count_flag) {
		cout << "Iteration times = " << counter << endl;
		fl << "\nIteration times = " << counter << "\n";
	}
	for (int i = 0; i < x.size(); i++) {
		cout << "x = " << x[i] << endl;
		fl << "\nx = " << x[i];
	}
	fl << "\n";
}

eq_Newton::eq_Newton() {
	input_fx();
	input_x();
	input_max_eps();
}

void eq_Newton::calc() {
	DerivExtra df(fx);
	for (int i = 0; i < x.size(); i++) {
		counter = 0;
		double xk;
		while (counter < maxcounter) {
			xk = x[i];
			x[i] = xk - fx.f(xk) / df.orderderiv(1, xk);
			if (fabs(fx.f(x[i])) < eps1) break;
			if (fabs(x[i] - xk) < eps2) break;
			counter++;
		}
		if (counter == maxcounter) cout << "\nWarning: Reach the maximum iteration number when calculate x = " << x[i] << endl;
	}
}

void eq_Newton::out_result() {
	cout << "\nCalculated by Newton method." << endl;
	fl << "\nCalculated by Newton method.";
	if (count_flag) {
		cout << "Iteration times = " << counter << endl;
		fl << "\nIteration times = " << counter << "\n";
	}
	for (int i = 0; i < x.size(); i++) {
		cout << "x = " << x[i] << endl;
		fl << "\nx = " << x[i];
	}
	fl << "\n";
}

eq_Secant::eq_Secant() {
	input_fx();
	input_x();
	input_max_eps();
}

void eq_Secant::input_x() {
	cout << "\nPlease enter the initial value or an interval [a,b]: (1 = initial value , 2 = [a,b]) ";
	int flag = in_int();
	if (flag == 1) {
		double temp;
		cout << "x0 = ";
		cin >> temp;
		x0.push_back(temp);
		cout << "x1 = ";
		cin >> temp;
		x.push_back(temp);
		count_flag = true;
	}
	else if (flag == 2) {
		double a, b;
		cout << "a = ";
		cin >> a;
		cout << "b = ";
		cin >> b;
		interval[0] = a;
		interval[1] = b;
		isolation();
		double dt = (interval[1] - interval[0]) * 0.00001;
		for (int i = 0; i < x.size(); i++) x0.push_back(x[i] - dt);
		count_flag = false;
	}
	else {
		cout << "Error: Illegal input." << endl;
		throw 0;
	}
}

void eq_Secant::calc() {
	for (int i = 0; i < x.size(); i++) {
		counter = 0;
		double xk, xk1, t;
		xk1 = x0[i];
		xk = x[i];
		t = fx.f(xk) / fx.f(xk1);
		x[i] = xk - t / (t - 1.0) * (xk - xk1);
		if (fabs(fx.f(x[i])) < eps1) continue;
		if (fabs(x[i] - xk) < eps2) continue;
		counter++;
		while (counter < maxcounter) {
			xk1 = xk;
			xk = x[i];
			t = fx.f(xk) / fx.f(xk1);
			x[i] = xk - t / (t - 1.0) * (xk - xk1);
			if (fabs(fx.f(x[i])) < eps1) break;
			if (fabs(x[i] - xk) < eps2) break;
			counter++;
		}
		if (counter == maxcounter) cout << "\nWarning: Reach the maximum iteration number when calculate x = " << x[i] << endl;
	}
}

void eq_Secant::out_result() {
	cout << "\nCalculated by Secant method." << endl;
	fl << "\nCalculated by Secant method.";
	if (count_flag) {
		cout << "Iteration times = " << counter << endl;
		fl << "\nIteration times = " << counter << "\n";
	}
	for (int i = 0; i < x.size(); i++) {
		cout << "x = " << x[i] << endl;
		fl << "\nx = " << x[i];
	}
	fl << "\n";
}

void eq_Relaxation::calc() {
	DerivExtra dphi(fx);
	for (int i = 0; i < x.size(); i++) {
		double omega = dphi.orderderiv(1, x[i]);
		counter = 0;
		double xk;
		while (counter < maxcounter) {
			xk = x[i];
			x[i] = (fx.f(xk) - omega * xk) / (1 - omega);
			if (fabs(fx.f(x[i])) < eps1) break;
			if (fabs(x[i] - xk) < eps2) break;
			counter++;
		}
		if (counter == maxcounter) cout << "\nWarning: Reach the maximum iteration number when calculate x = " << x[i] << endl;
	}
}

void eq_Relaxation::out_result() {
	cout << "\nCalculated by relaxation acceleration method." << endl;
	fl << "\nCalculated by relaxation acceleration method.";
	if (count_flag) {
		cout << "Iteration times = " << counter << endl;
		fl << "\nIteration times = " << counter << "\n";
	}
	for (int i = 0; i < x.size(); i++) {
		cout << "x = " << x[i] << endl;
		fl << "\nx = " << x[i];
	}
	fl << "\n";
}

void eq_Aitken::calc() {
	for (int i = 0; i < x.size(); i++) {
		counter = 0;
		double x1, x2;
		while (counter < maxcounter) {
			x1 = fx.f(x[i]);
			x2 = fx.f(x1);
			x[i] = x2 - sqr(x2 - x1) / (x2 - 2 * x1 + x[i]);
			if (fabs(fx.f(x[i])) < eps1) break;
			if (fabs(x[i] - x2) < eps2) break;
			counter++;
		}
		if (counter == maxcounter) cout << "\nWarning: Reach the maximum iteration number when calculate x = " << x[i] << endl;
	}
}

void eq_Aitken::out_result() {
	cout << "\nCalculated by Aitken acceleration method." << endl;
	fl << "\nCalculated by Aitken acceleration method.";
	if (count_flag) {
		cout << "Iteration times = " << counter << endl;
		fl << "\nIteration times = " << counter << "\n";
	}
	for (int i = 0; i < x.size(); i++) {
		cout << "x = " << x[i] << endl;
		fl << "\nx = " << x[i];
	}
	fl << "\n";
}

nl_eqs::nl_eqs() {
	fl.init("Nonlinear_equations.txt");
}

void nl_eqs::input_fx() {
	cout << "\nNonlinear equations: f(x) = 0" << endl;
	cout << "\nPlease enter the number of functions: n = ";
	n = in_int();
	fx.init(n);
	cout << "\nThe equations are: (MATLAB format)" << endl;
	for (int i = 0; i < n; i++) {
		cout << fx.get_fstr(i) << " = 0" << endl;
	}
	fl << "\nThe equations are: (MATLAB format)\n\n";
	for (int i = 0; i < n; i++) {
		fl << fx.get_fstr(i) << " = 0\n";
	}
}

void nl_eqs::input_x() {
	cout << "\nPlease enter the initial vector x0:\nx0(" << n << "*" << 1 << ") = ";
	double temp;
	for (int i = 0; i < n; i++) {
		cin >> temp;
		cout << "\t  ";
		x.push_back(temp);
	}
}

void nl_eqs::input_max_eps() {
	cout << "\nSet the errors and the maximum iteration times:" << endl;
	cout << "  |f(x_k+1)| < eps1 , or |x_k+1 - x_k| < eps2" << endl;
	cout << "  Default values: eps1 = 1.0e-6, eps2 = 1.0e-6, maximum iteration times = 10000" << endl;
	cout << "Use default values? (1 = Yes, 0 = No) ";
	int flag = in_int();
	if (flag == 1) {
		eps1 = 1.0e-6;
		eps2 = 1.0e-6;
		maxcounter = 10000;
	}
	else if (flag == 0) {
		cout << "eps1 = ";
		cin >> eps1;
		cout << "eps2 = ";
		cin >> eps2;
		cout << "maximum iteration times = ";
		maxcounter = in_int();
	}
	else {
		cout << "Error: Illegal input." << endl;
		throw 0;
	}
	if (eps1 < 1.0) {
		int pre;
		cout.setf(iostream::fixed);
		pre = -1 * floor(log10(eps1));
		cout.precision(pre);
	}
	fl.set_precision(eps1);
	fl << "\neps1 = " << eps1 << ", eps2 = " << eps2 << " (|f(x_k+1)| < eps1 , or |x_k+1 - x_k| < eps2)";
	fl << "\nmaximum iteration times = " << maxcounter << "\n";
}

eqs_Simple::eqs_Simple() {
	input_fx();
	input_x();
	input_max_eps();
}

void eqs_Simple::input_fx() {
	cout << "\nEquivalent transformation of nonlinear equation: x = g(x)" << endl;
	cout << "\nPlease enter the number of functions: n = ";
	n = in_int();
	fx.init(n, "g");
	cout << "\nThe equivalent transformations are: (MATLAB format)" << endl;
	for (int i = 0; i < n; i++) {
		cout << "x" << i + 1 << " = " << fx.get_fstr(i) << endl;
	}
	fl << "\nThe equivalent transformations are: (MATLAB format)\n\n";
	for (int i = 0; i < n; i++) {
		fl << "x" << i + 1 << " = " << fx.get_fstr(i) << "\n";
	}
}

void eqs_Simple::calc() {
	counter = 0;
	vector<double> xk;
	while (counter < maxcounter) {
		xk = x;
		x = fx.f(xk);
		for (int i = 0; i < n; i++) xk[i] -= x[i];
		if (vecnorm2(xk) < eps2) break;
		xk = fx.f(x);
		for (int i = 0; i < n; i++) xk[i] -= x[i];
		if (vecnorm2(xk) < eps1) break;
		counter++;
	}
	if (counter == maxcounter) cout << "\nWarning: Reach the maximum iteration number" << endl;
}

void eqs_Simple::out_result() {
	cout << "\nCalculated by simple iteration method." << endl;
	fl << "\nCalculated by simple iteration method.";
	cout << "Iteration times = " << counter << endl;
	fl << "\nIteration times = " << counter << "\n";
	cout << "\nx (" << n << "*" << 1 << ") =" << endl;
	for (int i = 0; i < n; i++) cout << "\t" << x[i] << endl;
	fl << "\nx = \n\n";
	fl.set_array_len(n) << x;
}

eqs_Newton::eqs_Newton() {
	input_fx();
	input_x();
	input_max_eps();
}

void eqs_Newton::calc() {
	vector<DerivExtra> Deriv;
	vector<vector<double>> J;
	J.resize(n);
	for (int i = 0; i < n; i++) J[i].resize(n);
	DerivExtra de(true);
	for (int i = 0; i < n; i++) {
		de.init(fx.getformula(i));
		Deriv.push_back(de);
	}
	counter = 0;
	vector<double> xk;
	vector<double> fxk;
	vector<double> dxk;
	xk.resize(n);
	fxk.resize(n);
	dxk.resize(n);
	while (counter < maxcounter) {
		xk = x;
		for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) J[i][j] = Deriv[i].deriv(xk, fx.get_xnum(i), j);
		fxk = fx.f(xk);
		for (int i = 0; i < n; i++) fxk[i] *= -1.0;
		Householder equ(n, n, J, fxk);
		equ.calc();
		dxk = equ.get_x();
		for (int i = 0; i < n; i++) x[i] += dxk[i];
		if (vecnorm2(dxk) / vecnorm2(xk) < eps1) break;
		if (vecnorm2(fx.f(x)) < eps2) break;
		counter++;
	}
	if (counter == maxcounter) cout << "\nWarning: Reach the maximum iteration number" << endl;
}

void eqs_Newton::out_result() {
	cout << "\nCalculated by Newton method." << endl;
	fl << "\nCalculated by Newton method.";
	cout << "Iteration times = " << counter << endl;
	fl << "\nIteration times = " << counter << "\n";
	cout << "\nx (" << n << "*" << 1 << ") =" << endl;
	for (int i = 0; i < n; i++) cout << "\t" << x[i] << endl;
	fl << "\nx = \n\n";
	fl.set_array_len(n) << x;
}

eqs_Secant::eqs_Secant() {
	input_fx();
	input_x();
	input_max_eps();
	input_h();
}

void eqs_Secant::input_h() {
	cout << "\nUse the default value of step size h = 0.001 ? (1 = Yes, 0 = No) ";
	int flag = in_int();
	if (flag == 1) h = 0.001;
	else if (flag == 0) {
		cout << "h = ";
		cin >> h;
	}
	else {
		cout << "Error: Illegal input." << endl;
		throw 0;
	}
	fl << "step size h = " << h << "\n";
}

void eqs_Secant::calc() {
	double temp;
	vector<vector<double>> A;
	vector<double> b;
	A.resize(n);
	for (int i = 0; i < n; i++) A[i].resize(n);
	b.resize(n);
	counter = 0;
	while (counter < maxcounter) {
		for (int j = 0; j < n; j++) {
			b = x;
			b[j] += h;
			b = fx.f(b);
			for (int i = 0; i < n; i++) A[i][j] = b[i];
		}
		b = fx.f(x);
		Householder equ(n, n, A, b);
		equ.calc();
		b = equ.get_x();
		temp = 0.0;
		for (int i = 0; i < n; i++) temp += b[i];
		temp -= 1.0;
		temp = h / temp;
		for (int i = 0; i < n; i++) {
			A[0][i] = temp * b[i];
			x[i] += A[0][i];
		}
		if (vecnorm2(A[0]) < eps1) break;
		if (vecnorm2(fx.f(x)) < eps2) break;
		counter++;
	}
}

void eqs_Secant::out_result() {
	cout << "\nCalculated by Secant method." << endl;
	fl << "\nCalculated by Secant method.";
	cout << "Iteration times = " << counter << endl;
	fl << "\nIteration times = " << counter << "\n";
	cout << "\nx (" << n << "*" << 1 << ") =" << endl;
	for (int i = 0; i < n; i++) cout << "\t" << x[i] << endl;
	fl << "\nx = \n\n";
	fl.set_array_len(n) << x;
}

eqs_Broyden::eqs_Broyden() {
	input_fx();
	input_x();
	input_max_eps();
}

void eqs_Broyden::calc() {
	double temp;
	vector<vector<double>> A(n);
	vector<vector<double>> t(n);
	for (int i = 0; i < n; i++) A[i].resize(n);
	for (int i = 0; i < n; i++) t[i].resize(n);
	vector<DerivExtra> Deriv;
	DerivExtra de(true);
	for (int i = 0; i < n; i++) {
		de.init(fx.getformula(i));
		Deriv.push_back(de);
	}
	for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] = Deriv[i].deriv(x, fx.get_xnum(i), j);
	A = mat_inverse(A);
	vector<double> s(n), y(n), xx(n), fy(n), fyy(n);
	xx = x;           //x_0
	fyy = fx.f(xx);   //fx_0
	for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) x[i] -= A[i][j] * fyy[j];   //x_1
	fy = fx.f(x);     //fx_1
	counter = 0;
	while (counter < maxcounter) {
		for (int i = 0; i < n; i++) s[i] = x[i] - xx[i];    //s_k
		for (int i = 0; i < n; i++) y[i] = fy[i] - fyy[i];  //y_k
		xx = x;   //x_k
		fyy = fy; //fx_k
		//propagate A. Start
		for (int i = 0; i < n; i++) {
			fy[i] = 0.0;
			for (int j = 0; j < n; j++) fy[i] += s[j] * A[j][i];
		}
		temp = 0.0;
		for (int i = 0; i < n; i++) temp += fy[i] * y[i];
		fy = s;
		for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) fy[i] -= A[i][j] * y[j];
		for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) t[i][j] = fy[i] * s[j];
		t = mat_multi(t, A);
		for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] += t[i][j] / temp;
		//propagate A. Finish
		for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) x[i] -= A[i][j] * fyy[j];
		for (int i = 0; i < n; i++) fy[i] = x[i] - xx[i];
		if (vecnorm2(fy) < eps1) break;
		fy = fx.f(x);
		if (vecnorm2(fy) < eps2) break;
		counter++;
	}
}

void eqs_Broyden::out_result() {
	cout << "\nCalculated by Broyden method." << endl;
	fl << "\nCalculated by Broyden method.";
	cout << "Iteration times = " << counter << endl;
	fl << "\nIteration times = " << counter << "\n";
	cout << "\nx (" << n << "*" << 1 << ") =" << endl;
	for (int i = 0; i < n; i++) cout << "\t" << x[i] << endl;
	fl << "\nx = \n\n";
	fl.set_array_len(n) << x;
}