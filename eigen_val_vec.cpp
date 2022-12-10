#include <iostream>
#include <cmath>
#include "common_fd.h"
#include "directmethod.h"
#include "eigen_val_vec.h"

using namespace std;

void Eigen_val_vec::init() {
	fl.init("Eigen_val_vec.txt");
	m = 0.0;
	counter = 0;
	cout << "\nPlease enter the order of the matrix: n = ";
	n = in_int();
	A.resize(n);
	for (int i = 0; i < n; i++) A[i].resize(n);
	z.resize(n);
	for (int i = 0; i < n; i++) z[i] = 1.0;
	cout << "\nPlease input the matrix A(" << n << "*" << n << "):" << endl;
	for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) cin >> A[i][j];
	cout << "\nA (" << n << "*" << n << ") =" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) cout << A[i][j] << "\t";
		cout << endl;
	}
	fl << "\nA (" << n << "*" << n << ") =\n\n";
	for (int i = 0; i < n; i++) {
		fl << "  ";
		for (int j = 0; j < n; j++) fl << A[i][j] << "\t";
		fl << "\n";
	}
	cout << "\nSet the error limit:" << endl;
	cout << "  |m_k - m_k-1| <= eps1 , or ||z_k - z_k-1|| <= eps2" << endl;
	cout << "  Default values: eps1 = 1.0e-6, eps2 = 1.0e-6" << endl;
	cout << "Use default values? (1 = Yes, 0 = No) ";
	int flag = in_int();
	if (flag == 1) {
		eps1 = 1.0e-6;
		eps2 = 1.0e-6;
	}
	else if (flag == 0) {
		cout << "eps1 = ";
		cin >> eps1;
		cout << "eps2 = ";
		cin >> eps2;
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
	fl << "\neps1 = " << eps1 << ", eps2 = " << eps2 << " (|m_k - m_k-1| <= eps1 , or ||z_k - z_k-1|| <= eps2)\n";
}

Power_method::Power_method() {
	init();
	cout << "\nWant to move the origin? (1 = Yes, 0 = No) ";
	int flag = in_int();
	if (flag == 1) {
		cout << "Please input the moving distance: p = ";
		cin >> p;
	}
	else if (flag == 0) {
		p = 0.0;
	}
	else {
		cout << "Error: Illegal input." << endl;
		throw 0;
	}
}

void Power_method::calc() {
	for (int i = 0; i < n; i++) A[i][i] -= p;
	vector<double> y(n);
	double temp;
	while (counter < 100000) {
		y = mat_multi_vec(A, z);
		temp = 0.0;
		for (int i = 0; i < n; i++) if (fabs(y[i]) > fabs(temp)) temp = y[i];
		if (temp == 0.0) {
			cout << "Error: The vector is 0." << endl;
			throw 0;
		}
		for (int i = 0; i < n; i++) y[i] /= temp;
		if (fabs(temp - m) <= eps1) {
			m = temp;
			z = y;
			break;
		}
		for (int i = 0; i < n; i++) z[i] -= y[i];
		if (vecnorm2(z) <= eps2) {
			m = temp;
			z = y;
			break;
		}
		m = temp;
		z = y;
		counter++;
	}
	if (counter == 100000) {
		cout << "Warning: The iteration times reach the default maximum 100000." << endl;
		fl << "\nThe iteration times reach the default maximum 100000.\n";
	}
	m += p;
}

void Power_method::out_result() {
	cout << "\nCalculated by power method." << endl;
	fl << "\nCalculated by power method.";
	cout << "\nIteration times = " << counter << endl;
	fl << "\nIteration times = " << counter << "\n";
	cout << "\nThe eigenvalue with the largest absolute value is" << endl;
	cout << "lambda = " << m << endl;
	fl << "\nThe eigenvalue with the largest absolute value is";
	fl << "\nlambda = " << m << "\n";
	if (p != 0.0) {
		cout << "( p = " << p << " )" << endl;
		fl << "( p = " << p << " )\n";
	}
	cout << "\nThe corresponding eigenvector is" << endl;
	cout << "xi =" << endl;
	fl << "\nThe corresponding eigenvector is";
	fl << "\nxi =";
	for (int i = 0; i < z.size(); i++) {
		cout << "  " << z[i] << endl;
		fl << "\n  " << z[i];
	}
	fl << "\n";
}

Inverse_power::Inverse_power() {
	init();
	cout << "\nPlease input the approximate eigenvalue: lambda ~= ";
	cin >> p;
}

void Inverse_power::calc() {
	for (int i = 0; i < n; i++) A[i][i] -= p;
	vector<vector<double>> L(n), R(n);
	for (int i = 0; i < n; i++) {
		L[i].resize(n);
		R[i].resize(n);
	}
	Doolittle LR(true);
	LR.decomposition(A, L, R);
	vector<double> u(n);
	vector<double> y(n);
	double temp;
	while (counter < 100000) {
		if (counter == 0) {
			for (int i = 0; i < n; i++) u[i] = 1.0;
		}
		else {
			LR.L_solve(L, u, z);
		}
		LR.U_solve(R, y, u);
		temp = 0.0;
		for (int i = 0; i < n; i++) if (fabs(y[i]) > fabs(temp)) temp = y[i];
		if (temp == 0.0) {
			cout << "Error: The vector is 0." << endl;
			throw 0;
		}
		for (int i = 0; i < n; i++) y[i] /= temp;
		if (fabs(temp - m) <= eps1) {
			m = temp;
			z = y;
			lambda = p + 1.0 / m;
			break;
		}
		for (int i = 0; i < n; i++) z[i] -= y[i];
		if (vecnorm2(z) <= eps2) {
			m = temp;
			z = y;
			lambda = p + 1.0 / m;
			break;
		}
		m = temp;
		z = y;
		lambda = p + 1.0 / m;
		counter++;
	}
	if (counter == 100000) {
		cout << "Warning: The iteration times reach the default maximum 100000." << endl;
		fl << "\nThe iteration times reach the default maximum 100000.\n";
	}
}

void Inverse_power::out_result() {
	cout << "\nCalculated by inverse power method." << endl;
	fl << "\nCalculated by inverse power method.";
	cout << "\nIteration times = " << counter << endl;
	fl << "\nIteration times = " << counter << "\n";
	cout << "\nThe approximate eigenvalue is" << endl;
	cout << "lambda ~= " << p << endl;
	fl << "\nThe approximate eigenvalue is";
	fl << "\nlambda ~= " << p << "\n";
	cout << "\nThe more accurate eigenvalue is" << endl;
	cout << "lambda = " << lambda << endl;
	fl << "\nThe more accurate eigenvalue is";
	fl << "\nlambda = " << lambda << "\n";
	cout << "\nThe corresponding eigenvector is" << endl;
	cout << "xi =" << endl;
	fl << "\nThe corresponding eigenvector is";
	fl << "\nxi =";
	for (int i = 0; i < z.size(); i++) {
		cout << "  " << z[i] << endl;
		fl << "\n  " << z[i];
	}
	fl << "\n";
}