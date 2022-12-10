#ifndef FILELOG_H
#define FILELOG_H

#include <string>
#include <fstream>
#include <vector>

class filelog {
private:
	const char* filename;
	enum file_type { FILE_TXT, FILE_M };
	file_type filetype;
	std::fstream* flog;
	int m, n;
	bool init_flag;
public:
	filelog();
	~filelog();
	void init(const char*);
	void judge_type();
	void set_width(int);
	void set_left_right(bool);
	void set_precision(double);
	filelog& set_array_len(int);
	filelog& set_mat_size(int, int);
	filelog& operator<<(const std::string);
	filelog& operator<<(const char);
	filelog& operator<<(const char*);
	filelog& operator<<(const int);
	filelog& operator<<(const int*);
	filelog& operator<<(const double);
	filelog& operator<<(const double*);
	filelog& operator<<(const std::vector<double>);
	filelog& operator<<(const double* const*);
};

#endif