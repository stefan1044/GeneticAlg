#pragma once
#include <cfloat>
#include <iomanip>

#include "generators.h"
#include "utils.h"

#include <chrono>
#include <string>
#include <iostream>

double runAlgorithm(const unsigned maxGenerations, const unsigned populationSize, const unsigned dimensions,
	const Function& function, const string functionName, const int sample);

void crossOver(vector<vector<bool>>& startingPopulation, const vector<bool>& chromosome1, const vector<bool>& chromosome2,
	uniform_int_distribution<>& distribution, mt19937& gen, const unsigned& nodeLength, int dimensions);
void mutate(vector<bool>& chromosome, const double& chanceToMutate, uniform_real_distribution<double>& distribution, mt19937& gen);


double runAlgorithm(const unsigned maxGenerations, const unsigned populationSize, const unsigned dimensions,
	const Function& function, const string functionName, const int sample) {

	auto start = chrono::system_clock::now();

	constexpr unsigned precision = 5;
	const auto N = abs((function.lowerBound - function.upperBound) * pow(10, precision));
	const auto nodeLength = static_cast<unsigned>(ceil(log2(N)));
	const double a = pow(2, nodeLength) - 1;
	const double average = function.upperBound - function.lowerBound;
	int L = nodeLength * dimensions;

	uniform_real_distribution<double> random01(0, 1.0);
	uniform_int_distribution<> randomGenerator(0, dimensions);
	uniform_int_distribution<> randomBool(0, 1);

	vector<vector<bool>> startingPopulation(populationSize, vector<bool>(nodeLength * dimensions));


	unsigned seed = chrono::system_clock::now().time_since_epoch().count() * sample;
	mt19937 gen;
	gen.seed(seed);
	
	for (unsigned i = 0; i < populationSize; i++) {
		for (unsigned j = 0; j < L; j++) {
			startingPopulation[i][j] = randomBool(gen);
		}
		//cout << "Debug: "; printVector(startingPopulation[i]); cout << '\n';
	}

	vector<double> eval(populationSize), results(populationSize);

	//Count Variables
	unsigned t = 0;
	unsigned generationsSinceLastImprovement = 0;
	double averageSelectedChromosomes = 0;

	const double chanceToMutate = 0.3 / nodeLength;
	double bestResult = DBL_MAX;
	//int randomGenerator = rand() % dimensions;
	
	constexpr int averageChromosomesToBeSelected = 8;

	while (t < maxGenerations) {
		t++;
		//cout << "Started generation " << t << '\n';

		double totalFitness = 0;

		for (unsigned i = 0; i < populationSize; i++) {
			results[i] = evaluate(startingPopulation[i], nodeLength, function, a, average);
			if (results[i] < bestResult) {
				bestResult = results[i];
				generationsSinceLastImprovement = 0;
			}
			eval[i] = function.fitFunc(results[i]);
			totalFitness += eval[i];
			//cout <<"results["<<results[i]<<"] "<< "eval[" << i << "] " << eval[i] << ' ';
		}
		//cout << "TotalFitness " << totalFitness << '\n';
		generationsSinceLastImprovement++;
		if (generationsSinceLastImprovement > 250) {
			break;
		}

		vector<vector<bool>> selectedChromosomes;
		vector<double> p(populationSize), q(populationSize + 1, 0);

		vector<double> rVector(averageChromosomesToBeSelected);

		for (int i = 0; i < rVector.size(); i++) {
			rVector[i] = random01(gen);
		}
		//Individual sel.prob. + Accumulated sel.prob.
		for (unsigned i = 0; i < populationSize; i++) {
			p[i] = eval[i] / totalFitness;
			q[i + 1] = q[i] + p[i];
			//cout << "p[" << i << "]: " << p[i] << " q[" << i << "] " << q[i] << ' ';


			for (int j = 0; j < rVector.size(); j++) {
				if (q[i] < rVector[j] && rVector[j] <= q[i + 1]) {
					//cout << "Selected " << " for cross-over";
					selectedChromosomes.push_back(startingPopulation[i]);
					rVector.erase(rVector.begin() + j);
					averageSelectedChromosomes++;
					break;
				}
			}
			//cout << '\n';
		}

		//Generate new population
		unsigned chromosomesAdded = 0;

		if (selectedChromosomes.empty() == true || selectedChromosomes.size() == 1) {
			t--;
			continue;
		}

		startingPopulation.clear();
		const unsigned numberOfSelectedChromosomes = selectedChromosomes.size();

		for (unsigned i = 0; i < numberOfSelectedChromosomes - numberOfSelectedChromosomes % 2; i++) {
			startingPopulation.push_back(selectedChromosomes[i]);
			chromosomesAdded++;
		}
		while (chromosomesAdded < populationSize) {
			for (unsigned i = 0; i < numberOfSelectedChromosomes - 1; i++) {
				for (unsigned j = i + 1; j < numberOfSelectedChromosomes; j++) {

					crossOver(startingPopulation, selectedChromosomes[i],
						selectedChromosomes[j], randomGenerator, gen, nodeLength, dimensions);
					mutate(startingPopulation[chromosomesAdded], chanceToMutate, random01, gen);
					mutate(startingPopulation[chromosomesAdded + 1], chanceToMutate, random01, gen);
					chromosomesAdded += 2;

					if (chromosomesAdded == populationSize)
						break;
				}
				if (chromosomesAdded == populationSize)
					break;
			}
		}

		//cout << "\nBest Result: " <<fixed<<setprecision(5)<< bestResult<<'\n';

	}

	auto end = chrono::system_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();

	string path = functionName + "_" + to_string(dimensions) + '_' + to_string(sample);

	cout << path << '\n';

	cout << "Best Result: " << fixed << setprecision(5) << bestResult << " Average Selected " << averageSelectedChromosomes / t << '\n';
	cout << "Ended at generation: " << t << "\n";
	cout << "Duration: " << duration << '\n' << '\n';

	return bestResult;
}


/*double runAlgorithm2(const unsigned maxGenerations, const unsigned populationSize, const unsigned dimensions, const Function& function)
{
	int t = 0;

}*/

void crossOver(vector<vector<bool>>& startingPopulation, const vector<bool>& chromosome1, const vector<bool>& chromosome2,
	uniform_int_distribution<>& distribution, mt19937& gen, const unsigned& nodeLength, int dimensions) {

	const unsigned cutPoint = distribution(gen) * nodeLength;
	vector<bool> newChromosome1, newChromosome2;

	for (unsigned i = 0; i < cutPoint; i++) {
		newChromosome1.push_back(chromosome1[i]);
		newChromosome2.push_back(chromosome2[i]);
	}
	for (unsigned i = cutPoint; i < chromosome1.size(); i++) {
		newChromosome1.push_back(chromosome2[i]);
		newChromosome2.push_back(chromosome1[i]);
	}
	startingPopulation.push_back(newChromosome1);
	startingPopulation.push_back(newChromosome2);

}

void mutate(vector<bool>& chromosome, const double& chanceToMutate, uniform_real_distribution<double>& distribution, mt19937& gen){
	for (auto var : chromosome) {
		if (distribution(gen) < chanceToMutate)
			var = !var;
	}
}
