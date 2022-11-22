#pragma once
#include <vector>

using namespace std;

constexpr double PI = 3.14159265358979323846;

typedef double (*functionPointer)(const vector<double>& values);
typedef double (*fitness)(const double& value, unsigned dimensions);

class Function
{
public:
	double lowerBound, upperBound;
	functionPointer func;
	fitness fitFunc;
	int dimensions;
public:
	Function(const double lowerBound, const double upperBound, const functionPointer function, const fitness fitFunc, int dimensions) {
		this->lowerBound = lowerBound;
		this->upperBound = upperBound;
		this->func = function;
		this->fitFunc = fitFunc;
		this->dimensions = dimensions;
	}

};

double deJongsFunction(const vector<double>& values);
double deJongsFitness(const double& value, unsigned dimensions);
double schwefelsFunction(const vector<double>& values);
double schwefelsFitness(const double& value, unsigned dimensions);
double rastringsFunction(const vector<double>& values);
double rastringsFitness(const double& value, unsigned dimensions);
double michalewiczsFunction(const vector<double>& values);
double michalewiczsFitness(const double& value, unsigned dimensions);


double deJongsFunction(const vector<double>& values) {
	double sum = 0;

	for (const auto var : values)
		sum += var * var;

	return sum;
}
double deJongsFitness(const double& value, unsigned dimensions) {
	return 1/value;
}
double schwefelsFunction(const vector<double>& values) {
	double sum = 0;

	for (const auto var : values)
		sum += (-var) * (sin(sqrt(abs(var))));

	return sum;
}
double schwefelsFitness(const double& value, unsigned dimensions) {
	return	-(value - 420 * dimensions);
}

double rastringsFunction(const vector<double>& values) {
	double sum = 0;
	const double n = values.size();

	for (const auto var : values)
		sum += var * var - 10 * cos(2 * PI * var);
	sum += 10 * n;

	return sum;
}
double rastringsFitness(const double& value, unsigned dimensions) {
	return 1 / value;
}
double michalewiczsFunction(const vector<double>& values) {
	double sum = 0;
	unsigned i = 1;

	for (const auto var : values) {
		sum -= sin(var) * pow(sin((i * var * var) / PI), 20);
		i++;
	}

	return sum;
}
double michalewiczsFitness(const double& value, unsigned dimensions) {
	return -(value - dimensions);
}

