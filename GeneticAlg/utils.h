#pragma once
#include "functions.h"
#include <iostream>
#include <vector>


using namespace std;


double convertToDomain(const unsigned& value, const double& lowerBound, const double& a, const double& average);
double evaluate(const vector<bool>& candidate, const unsigned& nodeLength, const Function& function, const double& a, const double&
	average);
unsigned nodeToNumber(const vector<bool>& node, const unsigned& currentNode, const unsigned& nodeLength);
template<class T> void printVector(vector<T> array);


double convertToDomain(const unsigned& value, const double& lowerBound, const double& a, const double& average) {
	return (value / a) * average + lowerBound;
}

double evaluate(const vector<bool>& candidate, const unsigned& nodeLength, const Function& function, const double& a, const double&
	average) {

	const unsigned length = candidate.size() / nodeLength;
	vector<double> values;

	for (unsigned i = 0; i < length; i++) {
		values.push_back(convertToDomain(nodeToNumber(candidate, i, nodeLength), function.lowerBound, a, average));
	}

	return function.func(values);
}

unsigned nodeToNumber(const vector<bool>& node, const unsigned& currentNode, const unsigned& nodeLength) {
	unsigned power = 1, number = 0;
	const unsigned temp = currentNode * nodeLength;

	for (unsigned i = temp; i < temp + nodeLength; i++) {
		number += node[i] * power;
		power *= 2;
	}

	return number;
}

template <class T>
void printVector(vector<T> array) {
	for (const auto var : array)
		cout << var << ' ';
	cout << '\n';
}
