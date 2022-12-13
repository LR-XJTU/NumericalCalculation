#include <iostream>
#include <cmath>
#include <cstdlib>
#include "formula.h"
#include "common_fd.h"

using namespace std;

//Initialization of class formula
void formula::init() {
	string str;
	showtips();
	cout << "\nf(x) = ";
	cin >> str;
	rf_str(str);
	fstr = str;
	trans_rpn(str);
	matlab_format();
	cout << "\nThe formula is f(x) = " << fstr << "  (MATLAB format)" << endl;
}

void formula::init(string fxstr) {
	rf_str(fxstr);
	fstr = fxstr;
	trans_rpn(fxstr);
	matlab_format();
}

void formula::init(string fx, string xt) {
	rf_str(fx);
	rf_str(xt);
	fx = replace_x(fx, xt);
	fstr = fx;
	trans_rpn(fx);
	matlab_format();
}

void formula::showtips() {
	cout << "\nPlease input the function. The following elements are allowed:\n" << endl;
	cout << "   81 , -3.67 , 9.425E+3 ---- Rational number" << endl;
	cout << "                  e , pi ---- Constant: natural logarithm , circumference ratio" << endl;
	cout << "                       x ---- variable" << endl;
	cout << "           + , - , * , / ---- Four basic operations" << endl;
	cout << "                       ^ ---- Power" << endl;
	//cout << "                       ! ---- Factorial" << endl;
	cout << "                    sqrt ---- Square root" << endl;
	cout << "                     exp ---- Exponent function based on e" << endl;
	cout << "                      ln ---- Nature logarithmic function" << endl;
	cout << "                     log ---- Logarithmic function -> Usage examples: log2(x), log10(x), log3.14(x)" << endl;
	cout << "         sin , cos , tan ---- Trigonometric function (Radians)" << endl;
	cout << "arcsin , arccos , arctan ---- Inverse trigonometric function (Radians)" << endl;
	cout << "                       | ---- Absolute value" << endl;
	cout << "                     ( ) ---- Brackets" << endl;
	cout << "* Note: DO NOT use Chinese keyboard. *" << endl;
}

void formula::showtips_xnum() {
	cout << "\nPlease input the function. The following elements are allowed:\n" << endl;
	cout << "   81 , -3.67 , 9.425E+3 ---- Rational number" << endl;
	cout << "                  e , pi ---- Constant: natural logarithm , circumference ratio" << endl;
	cout << "           x1 , x2 , ... ---- variables" << endl;
	cout << "           + , - , * , / ---- Four basic operations" << endl;
	cout << "                       ^ ---- Power" << endl;
	//cout << "                       ! ---- Factorial" << endl;
	cout << "                    sqrt ---- Square root" << endl;
	cout << "                     exp ---- Exponent function based on e" << endl;
	cout << "                      ln ---- Nature logarithmic function" << endl;
	cout << "                     log ---- Logarithmic function -> Usage examples: log2(x), log10(x), log3.14(x)" << endl;
	cout << "         sin , cos , tan ---- Trigonometric function (Radians)" << endl;
	cout << "arcsin , arccos , arctan ---- Inverse trigonometric function (Radians)" << endl;
	cout << "                       | ---- Absolute value" << endl;
	cout << "                     ( ) ---- Brackets" << endl;
	cout << "* Note: DO NOT use Chinese keyboard. *" << endl;
}

bool formula::x_flag() {
	for (int i = 0; i < frpn.size(); i++) if (frpn[i] == FX) return true;
	return false;
}

//define points which have no definition
void formula::define_xy() {
	int flag;
	do {
		double x, y;
		cout << "\nDefining a single point for f(x):" << endl;
		cout << "x = ";
		cin >> x;
		cout << "y = ";
		cin >> y;
		dfx.push_back(x);
		dfy.push_back(y);
		cout << "Want to continue to define? (1 = Yes , 0 = No)" << endl;
		flag = in_int();
	} while (flag == 1);
}

//If the point has been defined, return true
bool formula::find_xy(double x) {
	for (int i = 0; i < dfx.size(); i++) {
		if (dfx[i] == x) return true;
	}
	return false;
}

//Return f(x) of the defined point
double formula::list_xy(double x) {
	for (int i = 0; i < dfx.size(); i++) {
		if (dfx[i] == x) return dfy[i];
	}
	cout << "\nError: Fail to find the definition of (x,f(x))." << endl;
	throw 0;
}

void formula::set_fstr(string str) {
	fstr = str;
}

//Get fstr of formula
std::string formula::get_fstr() {
	return fstr;
}

string formula::replace_x(string fx, string xt) {
	size_t an = 0;
	while (fx.find(" x", an) != fx.npos) {
		an = fx.find(" x", an);
		an++;
		fx.replace(an, 1, xt);
		an += xt.size();
		if (an >= fx.size()) break;
	}
	return fx;
}

vector<op> formula::get_frpn() {
	return frpn;
}

queue<double> formula::get_fnum() {
	return fnum;
}

vector<double> formula::get_dfx() {
	return dfx;
}

vector<double> formula::get_dfy() {
	return dfx;
}

//reform the string to facilitate the transformation into RPN
void formula::rf_str(string& str) {
	int t;
	string s;
	//judge |
	int c_al = 0, c_ar = 0;
	for (int i = 0; i < str.size(); i++) {
		if (str.at(i) == '|') {
			if (i == 0) {
				if (str.size() > 1) {
					str.insert(1, "(");
					c_al++;
				}
				else {
					cout << "Error: Nonstandard input of | | in f(x)." << endl;
					throw 0;
				}
			}
			else {
				if ((str.at(i - 1) - '0' >= 0 && str.at(i - 1) - '9' <= 0) || str.at(i - 1) == '.' || str.at(i - 1) == 'x' || str.at(i - 1) == 'e' || str.at(i - 1) == 'i' || str.at(i - 1) == ')') {
					str.replace(i, 1, ")");
					c_ar++;
				}
				else {
					if (i + 1 < str.size()) {
						str.insert(i + 1, "(");
						c_al++;
					}
					else {
						cout << "Error: Nonstandard input of | | in f(x)." << endl;
						throw 0;
					}
				}
			}
		}
	}
	if (c_al != c_ar) {
		cout << "Error: Nonstandard input of | | in f(x)." << endl;
		throw 0;
	}
	//erase ' ' and '\t', and supplement ')'
	int c_l = 0, c_r = 0;
	for (int i = 0; i < str.size(); i++) {
		if (str.at(i) == '(') c_l++;
		if (str.at(i) == ')') c_r++;
		if (str.at(i) == ' ' || str.at(i) == '\t') {
			str.erase(i, 1);
			i--;
		}
	}
	if (c_l < c_r) {
		cout << "Error: Nonstandard input of brackets in f(x). There're "<<c_r<<" ) but "<<c_l<<" (." << endl;
		throw 0;
	}
	for (int i = 0; i < c_l - c_r; i++) str += ")";
	//standardadize log
	while (str.find("log") != str.npos) {
		int tp = str.find("log");
		str.replace(tp, 3, "ln");
		tp += 2;
		int len = str.find("(", tp) - tp;
		string ts = "/ln(" + str.substr(tp, len) + ")";
		str.erase(tp, len);
		skipbracket(str, tp);
		tp++;
		if (tp < str.size()) str.insert(tp, ts);
		else str += ts;
	}
	//handle with negative number
	if (str.at(0) == '-') str.insert(0, "0");
	for (int i = 0; i < str.size(); i++) {
		if (str.at(i) == '(' && i + 1 < str.size()) {
			if (str.at(i + 1) == '-') str.insert(i + 1, "0");
		}
	}
	//segmentation
	for (int i = 0; i < str.size(); i++) {
		if ((str.at(i) - '0' >= 0 && str.at(i) - '9' <= 0) || str.at(i) == '.') {
			str.insert(i, " ");
			i++;
			while (i < str.size()) {
				if ((str.at(i) - '0' >= 0 && str.at(i) - '9' <= 0) || str.at(i) == '.') i++;
				else if (str.at(i) == 'E') {
					i++;
					if (str.at(i) == '+' || str.at(i) == '-') i++;
					else {
						cout << "Error: Nonstandard input of scientific notation." << endl;
						throw 0;
					}
				}
				else break;
			}
			i--;
			continue;
		}
		if (str.at(i) == 'x') {
			str.insert(i, " ");
			i++;
			continue;
		}
		if (str.at(i) == '+' || str.at(i) == '-' || str.at(i) == '*' || str.at(i) == '/' || str.at(i) == '^' || str.at(i)=='!' || str.at(i) == '|' || str.at(i) == '(' || str.at(i) == ')') {
			str.insert(i, " ");
			i++;
			continue;
		}
		if (str.at(i) == 'e') {
			str.insert(i, " ");
			i++;
			if (i + 2 < str.size()) {
				std::string tempstr = str.substr(i, 3);
				if (tempstr.find("exp") != tempstr.npos) i += 2;
			}
			continue;
		}
		s = str.substr(i, 2);
		if (s.find("ln") != s.npos || s.find("pi") != s.npos) {
			str.insert(i, " ");
			i += 2;
			continue;
		}
		s = str.substr(i, 3);
		if (s.find("sin") != s.npos || s.find("cos") != s.npos || s.find("tan") != s.npos) {
			str.insert(i, " ");
			i += 3;
			continue;
		}
		s = str.substr(i, 4);
		if (s.find("sqrt") != s.npos) {
			str.insert(i, " ");
			i += 4;
			continue;
		}
		s = str.substr(i, 6);
		if (s.find("arcsin") != s.npos || s.find("arccos") != s.npos || s.find("arctan") != s.npos) {
			str.insert(i, " ");
			i += 6;
			continue;
		}
	}
	//check brackets
	size_t anchor;
	anchor = str.find("sqrt");
	while (anchor != str.npos) {
		check_brackets(str, anchor, 4);
		anchor++;
		anchor = str.find("sqrt", anchor);
	}
	anchor = str.find("exp");
	while (anchor != str.npos) {
		check_brackets(str, anchor, 3);
		anchor++;
		anchor = str.find("exp", anchor);
	}
	anchor = str.find("ln");
	while (anchor != str.npos) {
		check_brackets(str, anchor, 2);
		anchor++;
		anchor = str.find("ln", anchor);
	}
	anchor = str.find("sin");
	while (anchor != str.npos) {
		check_brackets(str, anchor, 3);
		anchor++;
		anchor = str.find("sin", anchor);
	}
	anchor = str.find("cos");
	while (anchor != str.npos) {
		check_brackets(str, anchor, 3);
		anchor++;
		anchor = str.find("cos", anchor);
	}
	anchor = str.find("tan");
	while (anchor != str.npos) {
		check_brackets(str, anchor, 3);
		anchor++;
		anchor = str.find("tan", anchor);
	}
	anchor = str.find("arcsin");
	while (anchor != str.npos) {
		check_brackets(str, anchor, 6);
		anchor++;
		anchor = str.find("arcsin", anchor);
	}
	anchor = str.find("arccos");
	while (anchor != str.npos) {
		check_brackets(str, anchor, 6);
		anchor++;
		anchor = str.find("arccos", anchor);
	}
	anchor = str.find("arctan");
	while (anchor != str.npos) {
		check_brackets(str, anchor, 6);
		anchor++;
		anchor = str.find("arctan", anchor);
	}
	//standardize the product terms
	t = 0;
	while (t < str.size()) {
		while (str.at(t) == ' ') t++;
		if ((str.at(t) - '0' >= 0 && str.at(t) - '9' <= 0) || str.at(t) == '.' || str.at(t) == 'x') {
			if (t - 2 >= 0) {
				if (str.at(t - 2) != '+' && str.at(t - 2) != '-' && str.at(t - 2) != '*' && str.at(t - 2) != '/' && str.at(t - 2) != '^' && str.at(t - 2) != '(') {
					str.insert(t, "* ");
					t += 2;
				}
			}
			while (str.at(t) != ' ') {
				t++;
				if (t == str.size()) break;
			}
			if (t + 1 < str.size()) {
				if (str.at(t + 1) != '+' && str.at(t + 1) != '-' && str.at(t + 1) != '*' && str.at(t + 1) != '/' && str.at(t + 1) != '^' && str.at(t + 1) != ')') {
					str.insert(t + 1, "* ");
					t += 2;
				}
			}
		}
		else if (str.at(t) == '(') {
			if (t - 2 >= 0) {
				if (str.at(t - 2) == ')') {
					str.insert(t, "* ");
					t += 2;
				}
			}
			t++;
		}
		else {
			while (str.at(t) != ' ') {
				t++;
				if (t == str.size()) break;
			}
		}
	}
}

//Transform into Reverse Polish
void formula::trans_rpn(std::string str) {
	stack <op> sta;
	std::string sdata;
	int ord1, ord2;
	op temp;
	for (int i = 0; i < str.size(); i++) {
		while (str.at(i) != ' ') {
			i++;
			if (i >= str.size()) break;
		}
		if (i >= str.size()) break;
		while (str.at(i) == ' ') i++;
		switch (str.at(i)) {
		case '+':
			temp = FADD;
			break;
		case '-':
			temp = FSUB;
			break;
		case '*':
			temp = FMUL;
			break;
		case '/':
			temp = FDIV;
			break;
		case '^':
			temp = FPOW;
			break;
		case '!':
			temp = FFAC;
			break;
		case '|':
			temp = FABS;
			break;
		case '(':
			temp = FBRAL;
			break;
		case ')':
			temp = FBRAR;
			break;
		case 'p':
			temp = FNUM;
			break;
		case 'e':
			if (i + 2 < str.size()) {
				std::string tempstr = str.substr(i, 3);
				if (tempstr.find("exp") != tempstr.npos) temp = FEXP;
				else temp = FNUM;
			}
			else temp = FNUM;
			break;
		case 's':
			if (i + 1 < str.size()) {
				if (str.at(i + 1) == 'q') temp = FSQRT;
				else if (str.at(i + 1) == 'i') temp = FSIN;
			}
			break;
		case 'c':
			temp = FCOS;
			break;
		case 't':
			temp = FTAN;
			break;
		case 'l':
			temp = FLN;
			break;
		case 'a':
			if (i + 3 < str.size()) {
				if (str.at(i + 3) == 's') temp = FASIN;
				else if (str.at(i + 3) == 'c') temp = FACOS;
				else if (str.at(i + 3) == 't') temp = FATAN;
			}
			break;
		default:
			sdata = str.substr(i, str.size() - i);
			size_t tp = sdata.find(" ");
			if (tp != sdata.npos) sdata.erase(tp, sdata.size() - tp);
			if (sdata.at(0) == 'x') temp = FX;
			else temp = FNUM;
			break;
		}
		if (temp == FNUM) {
			frpn.push_back(temp);
			if (str.at(i) == 'p')  fnum.push(3.14159265358979323846);
			else if (str.at(i) == 'e') fnum.push(2.718281828459045);
			else fnum.push(atof(sdata.c_str()));
			continue;
		}
		else if (temp == FX) {
			frpn.push_back(temp);
			continue;
		}
		else if (temp == FBRAL) {
			sta.push(temp);
			continue;
		}
		else if (temp == FBRAR) {
			while (sta.top() != FBRAL) {
				frpn.push_back(sta.top());
				sta.pop();
			}
			sta.pop();
			continue;
		}
		else {
			switch (temp) {
			case FADD:
			case FSUB:
				ord1 = 1;
				break;
			case FMUL:
			case FDIV:
				ord1 = 2;
				break;
			default:
				ord1 = 3;
				break;
			}
			if (sta.empty()) {
				sta.push(temp);
				continue;
			}
			switch (sta.top()) {
			case FBRAL:
				ord2 = 0;
				break;
			case FADD:
			case FSUB:
				ord2 = 1;
				break;
			case FMUL:
			case FDIV:
				ord2 = 2;
				break;
			default:
				ord2 = 3;
				break;
			}
			if (ord1 > ord2) sta.push(temp);
			else {
				frpn.push_back(sta.top());
				sta.pop();
				while (!sta.empty()) {
					switch (sta.top()) {
					case FBRAL:
						ord2 = 0;
						break;
					case FADD:
					case FSUB:
						ord2 = 1;
						break;
					case FMUL:
					case FDIV:
						ord2 = 2;
						break;
					default:
						ord2 = 3;
						break;
					}
					if (ord1 > ord2) break;
					else {
						frpn.push_back(sta.top());
						sta.pop();
					}
				}
				sta.push(temp);
			}
		}
	}
	while (!sta.empty()) {
		frpn.push_back(sta.top());
		sta.pop();
	}
}

//Transform into Matlab format
void formula::matlab_format() {
	while (fstr.find("|") != fstr.npos) fstr.replace(fstr.find("|"), 1, "abs");
	while (fstr.find("arcsin") != fstr.npos) fstr.replace(fstr.find("arcsin"), 6, "asin");
	while (fstr.find("arccos") != fstr.npos) fstr.replace(fstr.find("arccos"), 6, "acos");
	while (fstr.find("arctan") != fstr.npos) fstr.replace(fstr.find("arctan"), 6, "atan");
	while(fstr.find("ln")!=fstr.npos) fstr.replace(fstr.find("ln"), 2, "log");
	for (int i = 0; i < fstr.size(); i++) {
		if (fstr.at(i) == ' ') {
			fstr.erase(i, 1);
			i--;
		}
	}
	if (fstr.size() > 1) if (fstr.at(0) == '0' && fstr.at(1) == '-') fstr.erase(0, 1);
	for (int i = 0; i < fstr.size(); i++) {
		if (fstr.at(i) == '(' && i + 2 < fstr.size()) {
			if (fstr.at(i + 1) == '0' && fstr.at(i + 2) == '-') fstr.erase(i + 1, 1);
		}
	}
}

void formula::check_brackets(std::string &str, size_t anchor, int length) {
	string temp_s = str.substr(anchor, length + 2);
	if (temp_s.find("(") == temp_s.npos) {
		size_t temp_p = anchor + length + 1;
		if (temp_p >= str.size()) {
			cout << "Error: The input of formula is not standard." << endl;
			throw 0;
		}
		str.insert(temp_p, "( ");
		temp_p += 2;
		while (temp_p < str.size()) {
			if (str.at(temp_p) != ' ') temp_p++;
			else break;
		}
		if (temp_p < str.size()) str.insert(temp_p, " )");
		else str += " )";
	}
}

double formula::f(double x) {
	double t1, t2, t3;
	stack<double> rs;
	queue<double> data = fnum;
	vector<op>::iterator it;
	for (it = frpn.begin(); it != frpn.end(); it++) {
		switch (*it) {
		case FNUM:
			rs.push(data.front());
			data.pop();
			break;
		case FX:
			rs.push(x);
			break;
		case FADD:
			check_stack(rs, 2);
			t2 = rs.top();
			rs.pop();
			t1 = rs.top();
			rs.pop();
			t3 = t1 + t2;
			rs.push(t3);
			break;
		case FSUB:
			check_stack(rs, 2);
			t2 = rs.top();
			rs.pop();
			t1 = rs.top();
			rs.pop();
			t3 = t1 - t2;
			rs.push(t3);
			break;
		case FMUL:
			check_stack(rs, 2);
			t2 = rs.top();
			rs.pop();
			t1 = rs.top();
			rs.pop();
			t3 = t1 * t2;
			rs.push(t3);
			break;
		case FDIV:
			check_stack(rs, 2);
			t2 = rs.top();
			rs.pop();
			if (t2 == 0.0) {
				if (find_xy(t2)) {
					rs.push(list_xy(t2));
					break;
				}
				else {
					cout << "Error: The calculation meets the operation of division by zero." << endl;
					throw 0;
				}
			}
			t1 = rs.top();
			rs.pop();
			t3 = t1 / t2;
			rs.push(t3);
			break;
		case FPOW:
			check_stack(rs, 2);
			t2 = rs.top();
			rs.pop();
			t1 = rs.top();
			rs.pop();
			t3 = pow(t1, t2);
			rs.push(t3);
			break;
		case FFAC:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			t3 = calc_fac(t1);
			rs.push(t3);
			break;
		case FSQRT:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			if (t1 < 0.0) {
				if (find_xy(t1)) {
					rs.push(list_xy(t1));
					break;
				}
				else {
					cout << "Error: The calculation meets the operation of square root of a negative number." << endl;
					throw 0;
				}
			}
			t3 = sqrt(t1);
			rs.push(t3);
			break;
		case FEXP:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			t3 = exp(t1);
			rs.push(t3);
			break;
		case FLN:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			if (t1 <= 0.0) {
				if (find_xy(t1)) {
					rs.push(list_xy(t1));
					break;
				}
				else {
					cout << "Error: The calculation meets the operation of logarithm of zero or a negative number." << endl;
					throw 0;
				}
			}
			t3 = log(t1);
			rs.push(t3);
			break;
		case FSIN:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			t3 = sin(t1);
			rs.push(t3);
			break;
		case FCOS:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			t3 = cos(t1);
			rs.push(t3);
			break;
		case FTAN:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			t3 = tan(t1);
			rs.push(t3);
			break;
		case FASIN:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			if (fabs(t1) > 1.0) {
				if (find_xy(t1)) {
					rs.push(list_xy(t1));
					break;
				}
				else {
					cout << "Error: The calculation meets the operation of arcsin of a number whose absolute value is greater than 1." << endl;
					throw 0;
				}
			}
			t3 = asin(t1);
			rs.push(t3);
			break;
		case FACOS:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			if (fabs(t1) > 1.0) {
				if (find_xy(t1)) {
					rs.push(list_xy(t1));
					break;
				}
				else {
					cout << "Error: The calculation meets the operation of arccos of a number whose absolute value is greater than 1." << endl;
					throw 0;
				}
			}
			t3 = acos(t1);
			rs.push(t3);
			break;
		case FATAN:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			t3 = atan(t1);
			rs.push(t3);
			break;
		case FABS:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			t3 = fabs(t1);
			rs.push(t3);
			break;
		default:
			break;
		}
	}
	return rs.top();
}

double formula::f_xnum(vector<double> x) {
	double t1, t2, t3;
	stack<double> rs;
	queue<double> data = fnum;
	vector<op>::iterator it;
	int xnum = 0;
	for (it = frpn.begin(); it != frpn.end(); it++) {
		switch (*it) {
		case FNUM:
			rs.push(data.front());
			data.pop();
			break;
		case FX:
			rs.push(x[xnum]);
			xnum++;
			break;
		case FADD:
			check_stack(rs, 2);
			t2 = rs.top();
			rs.pop();
			t1 = rs.top();
			rs.pop();
			t3 = t1 + t2;
			rs.push(t3);
			break;
		case FSUB:
			check_stack(rs, 2);
			t2 = rs.top();
			rs.pop();
			t1 = rs.top();
			rs.pop();
			t3 = t1 - t2;
			rs.push(t3);
			break;
		case FMUL:
			check_stack(rs, 2);
			t2 = rs.top();
			rs.pop();
			t1 = rs.top();
			rs.pop();
			t3 = t1 * t2;
			rs.push(t3);
			break;
		case FDIV:
			check_stack(rs, 2);
			t2 = rs.top();
			rs.pop();
			if (t2 == 0.0) {
				if (find_xy(t2)) {
					rs.push(list_xy(t2));
					break;
				}
				else {
					cout << "Error: The calculation meets the operation of division by zero." << endl;
					throw 0;
				}
			}
			t1 = rs.top();
			rs.pop();
			t3 = t1 / t2;
			rs.push(t3);
			break;
		case FPOW:
			check_stack(rs, 2);
			t2 = rs.top();
			rs.pop();
			t1 = rs.top();
			rs.pop();
			t3 = pow(t1, t2);
			rs.push(t3);
			break;
		case FFAC:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			t3 = calc_fac(t1);
			rs.push(t3);
			break;
		case FSQRT:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			if (t1 < 0.0) {
				if (find_xy(t1)) {
					rs.push(list_xy(t1));
					break;
				}
				else {
					cout << "Error: The calculation meets the operation of square root of a negative number." << endl;
					throw 0;
				}
			}
			t3 = sqrt(t1);
			rs.push(t3);
			break;
		case FEXP:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			t3 = exp(t1);
			rs.push(t3);
			break;
		case FLN:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			if (t1 <= 0.0) {
				if (find_xy(t1)) {
					rs.push(list_xy(t1));
					break;
				}
				else {
					cout << "Error: The calculation meets the operation of logarithm of zero or a negative number." << endl;
					throw 0;
				}
			}
			t3 = log(t1);
			rs.push(t3);
			break;
		case FSIN:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			t3 = sin(t1);
			rs.push(t3);
			break;
		case FCOS:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			t3 = cos(t1);
			rs.push(t3);
			break;
		case FTAN:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			t3 = tan(t1);
			rs.push(t3);
			break;
		case FASIN:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			if (fabs(t1) > 1.0) {
				if (find_xy(t1)) {
					rs.push(list_xy(t1));
					break;
				}
				else {
					cout << "Error: The calculation meets the operation of arcsin of a number whose absolute value is greater than 1." << endl;
					throw 0;
				}
			}
			t3 = asin(t1);
			rs.push(t3);
			break;
		case FACOS:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			if (fabs(t1) > 1.0) {
				if (find_xy(t1)) {
					rs.push(list_xy(t1));
					break;
				}
				else {
					cout << "Error: The calculation meets the operation of arccos of a number whose absolute value is greater than 1." << endl;
					throw 0;
				}
			}
			t3 = acos(t1);
			rs.push(t3);
			break;
		case FATAN:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			t3 = atan(t1);
			rs.push(t3);
			break;
		case FABS:
			check_stack(rs, 1);
			t1 = rs.top();
			rs.pop();
			t3 = fabs(t1);
			rs.push(t3);
			break;
		default:
			break;
		}
	}
	return rs.top();
}

void formula::check_stack(std::stack<double> temp, int num) {
	double* t = new double[num];
	for (int i = 0; i < num; i++) {
		if (temp.empty()) {
			cout << "Error: Nonstandard input of formula causes the stack error." << endl;
			throw 0;
		}
		t[i] = temp.top();
		temp.pop();
	}
	for (int i = num - 1; i >= 0; i--) temp.push(t[i]);
	delete [] t;
}

void formula::skipbracket(string str, int &pos) {
	int temp = 0;
	while (pos < str.size()) {
		if (str.at(pos) == '(') temp++;
		else if (str.at(pos) == ')') temp--;
		if (temp == 0) break;
		pos++;
	}
	if (pos == str.size()) {
		cout << "Error: Fail to skip the brackets." << endl;
		throw 0;
	}
}

formula formula::operator=(formula fx) {
	fstr = fx.get_fstr();
	frpn = fx.get_frpn();
	fnum = fx.get_fnum();
	dfx = fx.get_dfx();
	dfy = fx.get_dfy();
	return *this;
}

formula formula::operator*(formula fx) {
	formula temp = *this;
	temp.fstr.insert(0, "(");
	temp.fstr += ")*(" + fx.fstr + ")";
	temp.frpn.insert(temp.frpn.end(), fx.frpn.begin(), fx.frpn.end());
	temp.frpn.push_back(FMUL);
	while (!fx.fnum.empty()) {
		temp.fnum.push(fx.fnum.front());
		fx.fnum.pop();
	}
	if (!this->dfx.empty()) {
		double tx, ty;
		for (int i = 0; i < this->dfx.size(); i++) {
			tx = this->dfx[i];
			if (fx.find_xy(tx)) ty = this->dfy[i] * fx.list_xy(tx);
			else ty = this->dfy[i] * fx.f(tx);
			temp.dfx.push_back(tx);
			temp.dfy.push_back(ty);
		}
	}
	if (!fx.dfx.empty()) {
		double tx, ty;
		for (int i = 0; i < fx.dfx.size(); i++) {
			tx = fx.dfx[i];
			if (this->find_xy(tx)) ty = fx.dfy[i] * this->list_xy(tx);
			else ty = fx.dfy[i] * this->f(tx);
			temp.dfx.push_back(tx);
			temp.dfy.push_back(ty);
		}
	}
	return temp;
}

void formulae::init(int n) {
	this->n = n;
	fx.resize(n);
	xnum.resize(n);
	string str;
	fx[0].showtips_xnum();
	cout << endl;
	for (int i = 0; i < n; i++) {
		cout << "f" << i + 1 << "(x1,...,x" << n << ") = ";
		cin >> str;
		fx[i].set_fstr(str);
	}
	recog_xnum();
}

void formulae::init(int n, const char* ch) {
	this->n = n;
	fx.resize(n);
	xnum.resize(n);
	string str;
	fx[0].showtips_xnum();
	cout << endl;
	for (int i = 0; i < n; i++) {
		cout << ch << i + 1 << "(x1,...,x" << n << ") = ";
		cin >> str;
		fx[i].set_fstr(str);
	}
	recog_xnum();
}

formula formulae::getformula(int i) {
	return fx[i];
}

vector<int> formulae::get_xnum(int i) {
	return xnum[i];
}

string formulae::get_fstr(int i) {
	return fx[i].get_fstr();
}

void formulae::recog_xnum() {
	string str, num;
	size_t anchor;
	int k;
	for (int i = 0; i < n; i++) {
		str = fx[i].get_fstr();
		anchor = 0;
		while (str.find("x", anchor) != str.npos) {
			anchor = str.find("x", anchor);
			k = 1;
			num.erase();
			if (anchor + k >= str.size()) {
				cout << "Error: Illegal input of xi" << endl;
				throw 0;
			}
			if (str.at(anchor + k) == 'p') {
				anchor++;
				continue;
			}
			while (str.at(anchor + k) >= '0' && str.at(anchor + k) <= '9') {
				num += str.at(anchor + k);
				k++;
				if (anchor + k >= str.size()) break;
			}
			if (k == 1) {
				cout << "Error: Illegal input of xi" << endl;
				throw 0;
			}
			str.replace(anchor, k, "x");
			k = StringToNum<int>(num);
			xnum[i].push_back(k);
			anchor++;
			if (anchor >= str.size()) break;
		}
		fx[i].init(str);
		str = fx[i].get_fstr();
		anchor = 0;
		k = 0;
		while (str.find("x", anchor) != str.npos) {
			anchor = str.find("x", anchor);
			if (anchor + 1 >= str.size()) {
				str += to_string(xnum[i][k]);
				break;
			}
			if (str.at(anchor + 1) == 'p') {
				anchor++;
				continue;
			}
			str.insert(anchor + 1, to_string(xnum[i][k]));
			k++;
			anchor++;
			if (anchor >= str.size()) break;
		}
		fx[i].set_fstr(str);
	}
}

vector<double> formulae::f(vector<double> x) {
	vector<vector<double>> xx;
	vector<double> fy;
	xx.resize(n);
	for (int i = 0; i < n; i++) xx[i].resize(xnum[i].size());
	fy.resize(n);
	for (int i = 0; i < n; i++) for (int j = 0; j < xx[i].size(); j++) xx[i][j] = x[xnum[i][j] - 1];
	for (int i = 0; i < n; i++) fy[i] = fx[i].f_xnum(xx[i]);
	return fy;
}