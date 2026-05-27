#include <iostream>
#include <cmath>
#include "filelog.h"
#include "common_fd.h"

using namespace std;

filelog::filelog() {
	init_flag = false;
}

filelog::~filelog() {
	if (init_flag) {
		flog->close();
		delete flog;
	}
}

void filelog::init(const char* fn) {
	filename = fn;
	judge_type();
	ifstream test(filename);
	if (test) {
		test.close();
		cout << '\n' << filename;
#ifdef CHINESE_VERSION
		cout << " 已经存在，覆盖还是继续写入？ (1 = 覆盖 , 2 = 继续写入)" << endl;
#else
		cout << " already exists. Overwrite it or continue at the end? (1 = overwrite , 2 = continue)" << endl;
#endif
		int temp;
		temp = in_int();
		if (temp == 1) {
			flog = new fstream(filename, ios::out);
		}
		else if (temp == 2) {
			flog = new fstream(filename, ios::out | ios::app);
		}
	}
	else {
		flog = new fstream(filename, ios::out);
	}
	flog->precision(15);
	if (filetype == FILE_TXT) (*flog) << "\n===================================================================\n";
	init_flag = true;
}

void filelog::judge_type() {
	string str = filename;
	if (str.find(".txt") != str.npos) filetype = FILE_TXT;
	else if (str.find(".m") != str.npos) filetype = FILE_M;
}

void filelog::set_width(int wid) {
	flog->width(wid);
}

void filelog::set_left_right(bool left) {
	if (left) flog->setf(ios::left);
	else flog->setf(ios::right);
}

void filelog::set_precision(double eps) {
	if (eps < 1.0) {
		int pre;
		flog->setf(iostream::fixed);
		pre = -1 * floor(log10(eps));
		flog->precision(pre);
	}
}

filelog& filelog::set_array_len(int len) {
	m = len;
	return *this;
}

filelog& filelog::set_mat_size(int mm, int nn) {
	m = mm;
	n = nn;
	return *this;
}

filelog& filelog::operator<<(const std::string str) {
	(*flog) << str;
	(*flog) << flush;
	return *this;
}

filelog& filelog::operator<<(const char ch) {
	(*flog) << ch;
	(*flog) << flush;
	return *this;
}

filelog& filelog::operator<<(const char* chs) {
	(*flog) << chs;
	(*flog) << flush;
	return *this;
}

filelog& filelog::operator<<(const int x) {
	(*flog) << x;
	(*flog) << flush;
	return *this;
}

filelog& filelog::operator<<(const int* arr) {
	for (int i = 0; i < m; i++) (*flog) << "\t" << arr[i] << "\n";
	(*flog) << flush;
	return *this;
}

filelog& filelog::operator<<(const double x) {
	(*flog) << x;
	(*flog) << flush;
	return *this;
}

filelog& filelog::operator<<(const double* arr) {
	for (int i = 0; i < m; i++) (*flog) << "\t" << arr[i] << "\n";
	(*flog) << flush;
	return *this;
}

filelog& filelog::operator<<(const vector<double> arr) {
	for (int i = 0; i < m; i++) (*flog) << "\t" << arr[i] << "\n";
	(*flog) << flush;
	return *this;
}

filelog& filelog::operator<<(const double* const* mat) {
	set_left_right(true);
	for (int i = 0; i < m; i++) {
		(*flog) << "\t";
		for (int j = 0; j < n; j++) {
			set_width(15);
			(*flog) << mat[i][j];
		}
		(*flog) << "\n";
	}
	(*flog) << flush;
	return *this;
}