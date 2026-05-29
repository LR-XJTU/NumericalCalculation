#include <iostream>
#include <cmath>
#include "int_diff.h"
#include "common_fd.h"

using namespace std;

Int_Diff::Int_Diff() {
	fl.init("Int_Diff.txt");
	fx.init();
	print_state = true;
}

Int_Diff::Int_Diff(formula fx) {
	this->fx = fx;
	print_state = false;
}

Int_Diff::Int_Diff(bool) {
	print_state = false;
}

//set the error limit eps
void Int_Diff::set_eps() {
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

double Int_Diff::get_result() {
	return result;
}

//construction function
Romberg::Romberg() {
	double a, b;
#ifdef CHINESE_VERSION
	cout << "\n请设置积分区间 [a , b]:\na = ";
#else
	cout << "\nPlease set the integral interval [a , b]:\na = ";
#endif
	cin >> a;
	cout << "b = ";
	cin >> b;
	interval[0] = a;
	interval[1] = b;
#ifdef CHINESE_VERSION
	cout << "\n在 [" << a << " , " << b << "] 内有未定义的点吗？(例如 f(x)=1/x 在 x=0 处无定义。) (1 = 有, 0 = 无)" << endl;
#else
	cout << "\nIs there any undefined point in [" << a << " , " << b << "]? (e.g. f(x)=1/x has no definition at x=0.) (1 = Yes , 0 = No)" << endl;
#endif
	int flag = in_int();
	if (flag == 1) fx.define_xy();
	set_eps();
	fl.set_precision(eps);
}

//construction function
Romberg::Romberg(formula fx, double a, double b) :Int_Diff(fx) {
	interval[0] = a;
	interval[1] = b;
	this->eps = 0.5e-8;
}

void Romberg::init() {
	//empty
}

//Romberg integration
void Romberg::calc() {
#ifdef CHINESE_VERSION
	if (print_state) cout << "\n正在使用Romberg积分法进行计算..." << endl;
#else
	if (print_state) cout << "\nUsing Romberg to calculate the integral..." << endl;
#endif

	double T2k[2], S2k[2], C2k[2], R2k[2];
	double e;
	T2k[0] = T(1);

	T2k[1] = T_Vari_step(0, T2k[0]);
	S2k[0] = S_ext(T2k);

	T2k[0] = T2k[1];
	T2k[1] = T_Vari_step(1, T2k[0]);
	S2k[1] = S_ext(T2k);
	C2k[0] = C_ext(S2k);

	T2k[0] = T2k[1];
	T2k[1] = T_Vari_step(2, T2k[0]);
	S2k[0] = S2k[1];
	S2k[1] = S_ext(T2k);
	C2k[1] = C_ext(S2k);
	R2k[0] = R_ext(C2k);

	T2k[0] = T2k[1];
	T2k[1] = T_Vari_step(3, T2k[0]);
	S2k[0] = S2k[1];
	S2k[1] = S_ext(T2k);
	C2k[0] = C2k[1];
	C2k[1] = C_ext(S2k);
	R2k[1] = R_ext(C2k);

	e = fabs(R2k[1] - R2k[0]);
	err.push_back(e);
	if (e < eps) {
		result = R2k[1];
		return;
	}

	int k = 4;
	do {
		T2k[0] = T2k[1];
		T2k[1] = T_Vari_step(k, T2k[0]);
		S2k[0] = S2k[1];
		S2k[1] = S_ext(T2k);
		C2k[0] = C2k[1];
		C2k[1] = C_ext(S2k);
		R2k[0] = R2k[1];
		R2k[1] = R_ext(C2k);
		e = fabs(R2k[1] - R2k[0]);
		err.push_back(e);
		k++;
	} while (e >= eps);
	result = R2k[1];

#ifdef CHINESE_VERSION
	if (print_state) cout << "计算完成。" << endl;
#else
	if (print_state) cout << "Finish calculating." << endl;
#endif
}

//cout and save the results
void Romberg::out_result() {
#ifdef CHINESE_VERSION
	cout << "\nf(x) 在 [a,b] 上的积分为 \n\nR = " << result << endl;
	fl << "\nf(x) = " << fx.get_fstr() << "\n";
	fl << "\nf(x) 在 [" << interval[0] << "," << interval[1] << "] 上的积分等于 \n\nR = " << result << "\n";
#else
	cout << "\nThe integral of f(x) on [a,b] is \n\nR = " << result << endl;
	fl << "\nf(x) = " << fx.get_fstr() << "\n";
	fl << "\nThe integral of f(x) on [" << interval[0] << "," << interval[1] << "] equals to \n\nR = " << result << "\n";
#endif
#ifdef CHINESE_VERSION
	cout << "\n误差表如下：\n\n计算次数 k\t\t误差" << endl;
	for (int i = 0; i < (int)err.size(); i++) cout << "\t" << i + 1 << "\t\t\t" << err[i] << endl;
	fl << "\n误差限为 eps = " << eps << "\n";
	fl << "\n误差表如下：\n\n计算次数 k\t\t误差\n";
#else
	cout << "\nThe error table is as follow:\n\ncalculation number k\t\terror" << endl;
	for (int i = 0; i < (int)err.size(); i++) cout << "\t" << i + 1 << "\t\t\t" << err[i] << endl;
	fl << "\nThe error limit is eps = " << eps << "\n";
	fl << "\nThe error table is as follow:\n\ncalculation number k\t\terror\n";
#endif
	for (int i = 0; i < (int)err.size(); i++) fl << "\t" << i + 1 << "\t\t\t" << err[i] << "\n";
	generate_m();
}

//generate .m file to get the figure of calculation error in MATLAB
void Int_Diff::generate_m() {
#ifdef WEB_VISUALIZATION
	filelog gm;
	gm.init("Integration.html");
	gm << "<!DOCTYPE html>" << "\n";
	gm << "<html><head>" << "\n";
	gm << "<meta charset='UTF-8'>" << "\n";
	gm << "<title>Integration Error</title>" << "\n";
	gm << "<script src='https://cdn.plot.ly/plotly-3.0.0.min.js'></script>" << "\n";
	gm << "</head><body>" << "\n";
	gm << "<div id='plot' style='width:100%;height:100vh;'></div>" << "\n";
	gm << "<script>" << "\n";
		int n=(int)err.size();
		gm << "var x=[";
		for(int i=0;i<n;i++) gm << (i+1) << ((i<n-1)?",":"");
		gm << "];" << "\n";
		gm << "var y=[";
		for(int i=0;i<n;i++) gm << log2(err[i]) << ((i<n-1)?",":"");
		gm << "];" << "\n";
		gm << "Plotly.newPlot(" << "\n";
			gm << "'plot'," << "\n";
			gm << "[{x:x,y:y,mode:'lines+markers',name:'log2(error)'}]," << "\n";
			gm << "{title:'Integration Error Curve',xaxis:{title:'Calculating times'},yaxis:{title:'Log2 of error'}});" << "\n";
	gm << "</script>" << "\n";
	gm << "</body></html>" << "\n";
#else
#ifdef CHINESE_VERSION
	cout << "\n是否生成 .m 文件用于在MATLAB中绘制误差曲线？(1 = 是, 0 = 否)" << endl;
#else
	cout << "\nWant to generate .m file to get the figure of calcultion error in MATLAB? (1 = Yes , 0 = No)" << endl;
#endif
	int flag = in_int();
	if (flag == 1) {
		filelog gm;
		gm.init("Integration.m");
		gm << "\nfigure\n";
		gm << "x = 1:" << (int)err.size() << ";\n";
		gm << "data = [";
		for (int i = 0; i < (int)err.size(); i++) gm << err[i] << " ";
		gm << "];\n";
		gm << "y=log2(data);\n";
		gm << "plot(y,'rx-')\n";
		gm << "xlabel('Calculating times (times)');\n";
		gm << "ylabel('Logarithm of error');\n";
		gm << "title('The curve of error varying with calculating times');\n";
	}

#endif
}

//trapezoidal formula
double Romberg::T(int n) {
	double h = (interval[1] - interval[0]) / n;
	double xi;
	double Tn = fx.f(interval[0]) + fx.f(interval[1]);
	for (int i = 1; i <= n - 1; i++) {
		xi = interval[0] + i * h;
		Tn += 2.0 * fx.f(xi);
	}
	Tn *= h / 2.0;
	return Tn;
}

//trapezoidal formula with variable step size
double Romberg::T_Vari_step(int k, double T2k) {
	double T2k1 = 0.0;
	double xi;
	for (int i = 1; i <= int_pow(2, k); i++) {
		xi = interval[0] + (2.0 * i - 1.0) * (interval[1] - interval[0]) / int_pow(2, k + 1);
		T2k1 += fx.f(xi);
	}
	T2k1 *= (interval[1] - interval[0]) / int_pow(2, k + 1);
	T2k1 += T2k / 2.0;
	return T2k1;
}

//calculate S by extrapolation
double Romberg::S_ext(double Tn[2]) {
	double Sn = Tn[1] + (Tn[1] - Tn[0]) / 3.0;
	return Sn;
}

//calculate C by extrapolation
double Romberg::C_ext(double Sn[2]) {
	double Cn = Sn[1] + (Sn[1] - Sn[0]) / 15.0;
	return Cn;
}

//calculate R by extrapolation
double Romberg::R_ext(double Cn[2]) {
	double Rn = Cn[1] + (Cn[1] - Cn[0]) / 63.0;
	return Rn;
}

DerivExtra::DerivExtra() {
	int order;
	double x, h;
#ifdef CHINESE_VERSION
	cout << "\n请输入求导阶数：order = ";
#else
	cout << "\nPlease enter the order of the derivative: order = ";
#endif
	order = in_int();
	this->order = order;
#ifdef CHINESE_VERSION
	cout << "\n请输入 x：\nx = ";
#else
	cout << "\nPlease enter x:\nx = ";
#endif
	cin >> x;
	this->x = x;
#ifdef CHINESE_VERSION
	cout << "\n请设置初始步长 h：(区间 [x-h,x+h] 不应包含未定义的点)\nh = ";
#else
	cout << "\nPlease set the initial step h:(The interval [x-h,x+h] should not contain undefined points)\nh = ";
#endif
	cin >> h;
	this->h = h;
	set_eps();
	fl.set_precision(eps);
}

DerivExtra::DerivExtra(formula fx) :Int_Diff(fx) {
	h = 0.1;
	eps = 0.5e-8;
}

DerivExtra::DerivExtra(bool) : Int_Diff(true) {}

void DerivExtra::init(formula fx) {
	this->fx = fx;
	print_state = false;
	h = 0.1;
	eps = 0.5e-8;
}

double DerivExtra::deriv(double xx) {
	double hi = h;
	int n = 2;
	double e;
	vector<vector<double>> T;
	T.resize(n);
	for (int i = 0; i < n; i++) T[i].resize(n + 1);

	T[0][0] = (fx.f(xx + hi) - fx.f(xx - hi)) / 2.0 / hi;
	hi /= 2.0;
	T[0][1] = (fx.f(xx + hi) - fx.f(xx - hi)) / 2.0 / hi;
	hi /= 2.0;
	T[1][0] = T[0][1] + (T[0][1] - T[0][0]) / (int_pow(4, 1) - 1.0);
	T[0][2] = (fx.f(xx + hi) - fx.f(xx - hi)) / 2.0 / hi;
	hi /= 2.0;
	T[1][1] = T[0][2] + (T[0][2] - T[0][1]) / (int_pow(4, 1) - 1.0);
	e = fabs(T[1][1] - T[1][0]);
	if (err_flag) err.push_back(e);
	while (e >= eps) {
		n++;
		T.resize(n);
		for (int i = 0; i < n; i++) T[i].resize(n + 1);
		T[n - 1][0] = T[n - 2][1] + (T[n - 2][1] - T[n - 2][0]) / (int_pow(4, n - 1) - 1.0);
		T[0][n] = (fx.f(xx + hi) - fx.f(xx - hi)) / 2.0 / hi;
		hi /= 2.0;
		for (int i = 1; i < n; i++) T[i][n - i] = T[i - 1][n + 1 - i] + (T[i - 1][n + 1 - i] - T[i - 1][n - i]) / (int_pow(4, i) - 1.0);
		e = fabs(T[n - 1][1] - T[n - 1][0]);
		if (err_flag) err.push_back(e);
	}
	return T[n - 1][1];
}

vector<double> DerivExtra::dxp(vector<double> xx, vector<int> xnum, int j, double hi) {
	vector<double> xxx;
	xxx.resize(xnum.size());
	xx[j] += hi;
	for (size_t i = 0; i < xnum.size(); i++) xxx[i] = xx[xnum[i] - 1];
	return xxx;
}

vector<double> DerivExtra::dxn(vector<double> xx, vector<int> xnum, int j, double hi) {
	vector<double> xxx;
	xxx.resize(xnum.size());
	xx[j] -= hi;
	for (size_t i = 0; i < xnum.size(); i++) xxx[i] = xx[xnum[i] - 1];
	return xxx;
}

double DerivExtra::deriv(vector<double> xx, vector<int> xnum, int j) {
	double hi = h;
	int n = 2;
	double e;
	vector<vector<double>> T;
	T.resize(n);
	for (int i = 0; i < n; i++) T[i].resize(n + 1);
	vector<double> xp, xn;
	xp.resize(xnum.size());
	xn.resize(xnum.size());
	xp = dxp(xx, xnum, j, hi);
	xn = dxn(xx, xnum, j, hi);
	T[0][0] = (fx.f_xnum(xp) - fx.f_xnum(xn)) / 2.0 / hi;
	hi /= 2.0;

	xp = dxp(xx, xnum, j, hi);
	xn = dxn(xx, xnum, j, hi);
	T[0][1] = (fx.f_xnum(xp) - fx.f_xnum(xn)) / 2.0 / hi;
	hi /= 2.0;

	T[1][0] = T[0][1] + (T[0][1] - T[0][0]) / (int_pow(4, 1) - 1.0);
	xp = dxp(xx, xnum, j, hi);
	xn = dxn(xx, xnum, j, hi);
	T[0][2] = (fx.f_xnum(xp) - fx.f_xnum(xn)) / 2.0 / hi;
	hi /= 2.0;

	T[1][1] = T[0][2] + (T[0][2] - T[0][1]) / (int_pow(4, 1) - 1.0);
	e = fabs(T[1][1] - T[1][0]);
	if (err_flag) err.push_back(e);
	while (e >= eps) {
		n++;
		T.resize(n);
		for (int i = 0; i < n; i++) T[i].resize(n + 1);
		T[n - 1][0] = T[n - 2][1] + (T[n - 2][1] - T[n - 2][0]) / (int_pow(4, n - 1) - 1.0);
		xp = dxp(xx, xnum, j, hi);
		xn = dxn(xx, xnum, j, hi);
		T[0][n] = (fx.f_xnum(xp) - fx.f_xnum(xn)) / 2.0 / hi;
		hi /= 2.0;
		for (int i = 1; i < n; i++) T[i][n - i] = T[i - 1][n + 1 - i] + (T[i - 1][n + 1 - i] - T[i - 1][n - i]) / (int_pow(4, i) - 1.0);
		e = fabs(T[n - 1][1] - T[n - 1][0]);
		if (err_flag) err.push_back(e);
	}
	return T[n - 1][1];
}

double DerivExtra::orderderiv(int od, double xx) {
	if (od == 1) return deriv(xx);
	double hi = h;
	int n = 2;
	double e;
	vector<vector<double>> T;
	T.resize(n);
	for (int i = 0; i < n; i++) T[i].resize(n + 1);

	T[0][0] = (orderderiv(od - 1, xx + hi) - orderderiv(od - 1, xx - hi)) / 2.0 / hi;
	hi /= 2.0;
	T[0][1] = (orderderiv(od - 1, xx + hi) - orderderiv(od - 1, xx - hi)) / 2.0 / hi;
	hi /= 2.0;
	T[1][0] = T[0][1] + (T[0][1] - T[0][0]) / (int_pow(4, 1) - 1.0);
	T[0][2] = (orderderiv(od - 1, xx + hi) - orderderiv(od - 1, xx - hi)) / 2.0 / hi;
	hi /= 2.0;
	T[1][1] = T[0][2] + (T[0][2] - T[0][1]) / (int_pow(4, 1) - 1.0);
	e = fabs(T[1][1] - T[1][0]);
	if (od == order) err.push_back(e);
	while (e >= eps) {
		n++;
		T.resize(n);
		for (int i = 0; i < n; i++) T[i].resize(n + 1);
		T[n - 1][0] = T[n - 2][1] + (T[n - 2][1] - T[n - 2][0]) / (int_pow(4, n - 1) - 1.0);
		T[0][n] = (orderderiv(od - 1, xx + hi) - orderderiv(od - 1, xx - hi)) / 2.0 / hi;
		hi /= 2.0;
		for (int i = 1; i < n; i++) T[i][n - i] = T[i - 1][n + 1 - i] + (T[i - 1][n + 1 - i] - T[i - 1][n - i]) / (int_pow(4, i) - 1.0);
		e = fabs(T[n - 1][1] - T[n - 1][0]);
		if (od == order) err.push_back(e);
	}
	return T[n - 1][1];
}

void DerivExtra::calc() {
#ifdef CHINESE_VERSION
	if (print_state) cout << "\n正在使用外推法计算导数..." << endl;
#else
	if (print_state) cout << "\nUsing extrapolation to calculate the derivative..." << endl;
#endif
	if (order == 1) {
		err_flag = true;
		result = deriv(x);
	}
	else if (order >= 2) {
		err_flag = false;
		result = orderderiv(order,x);
	}
#ifdef CHINESE_VERSION
	if (print_state) cout << "计算完成。" << endl;
#else
	if (print_state) cout << "Finish calculating." << endl;
#endif
}

void DerivExtra::out_result() {
	if (order > 3) {
#ifdef CHINESE_VERSION
		cout << "\nf[" << order << "](" << x << ") = " << result << "   ([" << order << "] 表示导数阶数)" << endl;
		fl << "\nf(x) = " << fx.get_fstr() << "\n";
		fl << "\nf[" << order << "](" << x << ") = " << result << "   ([" << order << "] 是导数的阶数)\n";
#else
		cout << "\nf[" << order << "](" << x << ") = " << result << "   ([" << order << "] is the order of the derivative)" << endl;
		fl << "\nf(x) = " << fx.get_fstr() << "\n";
		fl << "\nf[" << order << "](" << x << ") = " << result << "   ([" << order << "] is the order of the derivative)\n";
#endif
	}
	else {
		cout << "\nf";
		for (int i = 1; i <= order; i++) cout << "'";
#ifdef CHINESE_VERSION
		cout << "(" << x << ") = " << result << endl;
		fl << "\nf(x) = " << fx.get_fstr() << "\n";
		fl << "\nf";
		for (int i = 1; i <= order; i++) fl << "'";
		fl << "(" << x << ") = " << result << "\n";
#else
		cout << "(" << x << ") = " << result << endl;
		fl << "\nf(x) = " << fx.get_fstr() << "\n";
		fl << "\nf";
		for (int i = 1; i <= order; i++) fl << "'";
		fl << "(" << x << ") = " << result << "\n";
#endif
	}
#ifdef CHINESE_VERSION
	cout << "\n误差表如下：\n\n计算次数 k\t\t误差" << endl;
	for (int i = 0; i < (int)err.size(); i++) cout << "\t" << i + 1 << "\t\t\t" << err[i] << endl;
	fl << "\n误差限为 eps = " << eps << "\n";
	fl << "\n误差表如下：\n\n计算次数 k\t\t误差\n";
#else
	cout << "\nThe error table is as follow:\n\ncalculation number k\t\terror" << endl;
	for (int i = 0; i < (int)err.size(); i++) cout << "\t" << i + 1 << "\t\t\t" << err[i] << endl;
	fl << "\nThe error limit is eps = " << eps << "\n";
	fl << "\nThe error table is as follow:\n\ncalculation number k\t\terror\n";
#endif
	for (int i = 0; i < (int)err.size(); i++) fl << "\t" << i + 1 << "\t\t\t" << err[i] << "\n";
	generate_m();
}

double DerivExtra::get_result() {
	return result;
}
