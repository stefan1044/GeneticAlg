#pragma once

#include "geneticAlgorithm.h"



enum paramNum { CM, CC, AC, ES };


void writeResults(ModParams parameters, double average, string nameOfFunc, unsigned dimensions) {
	string fileName = "Meta" + nameOfFunc + to_string(dimensions) + ".out";
	ofstream out(fileName);

	out << parameters.chanceToMutate << '\n';
	out << parameters.chanceToCrossOver << '\n';
	out << parameters.averageChromosomesToBeSelected << '\n';
	out << parameters.numberOfElites << '\n';
	out << average;
}

ModParams readResults(string fileName) {
	ifstream in(fileName);
	ModParams parameters;

	in >> parameters.chanceToMutate >> parameters.chanceToCrossOver >> parameters.averageChromosomesToBeSelected >> parameters.numberOfElites;

	return parameters;
}

void coutResults(ModParams parameters, double average, unsigned iteration) {
	cout << "Iteration: " << iteration << '\n';
	cout << "Chance to mutate: " << parameters.chanceToMutate << '\n';
	cout << "Chance to crossover: " << parameters.chanceToCrossOver << '\n';
	cout << "Average selected: " << parameters.averageChromosomesToBeSelected << '\n';
	cout << "Elites Selected: " << parameters.numberOfElites << '\n';
	cout << "With average result: " << average << "\n\n";
}

void metaAlg(Function func, unsigned dimension, string nameOfFunc) {

	const double CMRLower = 0, CMRUpper = 1;
	const double CCRLower = 0.1, CCRUpper = 1;
	const int ACRLower = 50, ACRUpper = 120;
	const int ESRLower = 1, ESRUpper = 30;

	const unsigned repetitions = 20;


	vector<double> stepsCM = { 0.01, 0.03, 0.05, 0.07 };
	vector<double> stepsCC = { 0.01, 0.03, 0.05, 0.07 };
	vector<int> stepsAC = { 1,3,5,8 };
	vector<int> stepsES = { 1,2,3,4 };



	auto parameters = readResults("Meta" + nameOfFunc + to_string(dimension) + ".out");

	unsigned iterationsRan = 0;
	paramNum paramNumber = CM;

	vector<double> results(8);

	double initAvg = 0;
	for (unsigned j = 0; j < repetitions; j++) {
		initAvg += runAlgorithm(2000, 200, dimension, func, nameOfFunc, j, parameters);
	}
	initAvg /= repetitions;
	double globalBest = initAvg;
	ModParams bestParams = parameters;


	cout << "Running " << nameOfFunc << " in " << dimension << " dimensions\n";
	cout << "Initial results: " << '\n';
	cout << "Chance to mutate: " << bestParams.chanceToMutate << '\n';
	cout << "Chance to crossover: " << bestParams.chanceToCrossOver << '\n';
	cout << "Average selected: " << bestParams.averageChromosomesToBeSelected << '\n';
	cout << "Elites Selected: " << bestParams.numberOfElites << '\n';
	cout << "With average result: " << globalBest << "\n\n";

	while (iterationsRan < 16) {
		double bestResult = globalBest;
		if (paramNumber == CM) {
			vector<double> values;
			double initParam = parameters.chanceToMutate;

			for (auto var : stepsCM) {
				double up = initParam + var;
				double down = initParam - var;
				if (CMRLower < up && up < CMRUpper) {
					values.push_back(up);
				}
				if (CMRLower < down && down < CMRUpper) {
					values.push_back(down);
				}
			}

			for (unsigned i = 0; i < values.size(); i++) {
				parameters.chanceToMutate = values[i];
				double med = 0;
				for (unsigned j = 0; j < repetitions; j++) {
					med += runAlgorithm(2000, 200, dimension, func, nameOfFunc, j, parameters);
				}
				med /= repetitions;
				results[i] = med;
			}

			bool changed = false;
			for (unsigned i = 0; i < values.size(); i++) {
				if (results[i] < bestResult) {
					bestResult = results[i];
					changed = true;
					parameters.chanceToMutate = values[i];
				}
			}
			if (changed == false)
				parameters.chanceToMutate = initParam;


			paramNumber = CC;
		}
		else if (paramNumber == CC) {
			vector<double> values;
			double range = 0.9;
			double initParam = parameters.chanceToCrossOver;

			for (auto var : stepsCC) {

				double up = initParam + var;
				double down = initParam - var;
				if (CCRLower < up && up <= CCRUpper) {
					values.push_back(up);
				}
				if (CCRLower < down && down <= CCRUpper) {
					values.push_back(down);
				}
			}

			for (unsigned i = 0; i < values.size(); i++) {
				parameters.chanceToCrossOver = values[i];
				double med = 0;
				for (unsigned j = 0; j < repetitions; j++) {
					med += runAlgorithm(2000, 200, dimension, func, nameOfFunc, j, parameters);
				}
				med /= repetitions;
				results[i] = med;
			}

			bool changed = false;
			for (unsigned i = 0; i < values.size(); i++) {
				if (results[i] < bestResult) {
					bestResult = results[i];
					changed = true;
					parameters.chanceToCrossOver = values[i];
				}
			}
			if (changed == false)
				parameters.chanceToCrossOver = initParam;

			paramNumber = AC;
		}
		else if (paramNumber == AC) {
			vector<int> values;
			int initParam = parameters.averageChromosomesToBeSelected;

			for (auto var : stepsAC) {
				int up = initParam + var;
				int down = initParam - var;
				if (ACRLower <= up && up <= ACRUpper) {
					values.push_back(up);
				}
				if (ACRLower <= down && down <= ACRUpper) {
					values.push_back(down);
				}
			}

			for (unsigned i = 0; i < values.size(); i++) {
				parameters.averageChromosomesToBeSelected = values[i];
				double med = 0;
				for (unsigned j = 0; j < repetitions; j++) {
					med += runAlgorithm(2000, 200, dimension, func, nameOfFunc, j, parameters);
				}
				med /= repetitions;
				results[i] = med;
			}

			bool changed = false;
			for (unsigned i = 0; i < values.size(); i++) {
				if (results[i] < bestResult) {
					bestResult = results[i];
					changed = true;
					parameters.averageChromosomesToBeSelected = values[i];
				}
			}
			if (changed == false)
				parameters.averageChromosomesToBeSelected = initParam;

			paramNumber = ES;
		}
		else {
			vector<int> values;
			int range = 80;
			int initParam = parameters.numberOfElites;

			for (auto var : stepsES) {

				int up = initParam + var;
				int down = initParam - var;
				if (ESRLower <= up && up <= ESRUpper) {
					values.push_back(up);
				}
				if (ESRLower <= down && down <= ESRUpper) {
					values.push_back(down);
				}
			}

			for (unsigned i = 0; i < values.size(); i++) {
				parameters.numberOfElites = values[i];
				double med = 0;
				for (unsigned j = 0; j < repetitions; j++) {
					med += runAlgorithm(2000, 200, dimension, func, nameOfFunc, j, parameters);
				}
				med /= repetitions;
				results[i] = med;
			}

			bool changed = false;
			for (unsigned i = 0; i < values.size(); i++) {
				if (results[i] < bestResult) {
					bestResult = results[i];
					changed = true;
					parameters.numberOfElites = values[i];
				}
			}
			if (changed == false)
				parameters.numberOfElites = initParam;

			paramNumber = CM;
		}


		coutResults(parameters, bestResult, iterationsRan);

		if (bestResult < globalBest) {
			globalBest = bestResult;
			bestParams = parameters;
		}
		iterationsRan++;
		writeResults(parameters, globalBest, nameOfFunc, dimension);
	}


	cout << "Final results: " << '\n';
	cout << "Chance to mutate: " << bestParams.chanceToMutate << '\n';
	cout << "Chance to crossover: " << bestParams.chanceToCrossOver << '\n';
	cout << "Average selected: " << bestParams.averageChromosomesToBeSelected << '\n';
	cout << "Elites Selected: " << bestParams.numberOfElites << '\n';
	cout << "With average result: " << globalBest << "\n\n";


}
