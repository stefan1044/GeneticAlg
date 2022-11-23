#pragma once
#include <algorithm>
#include <cfloat>
#include <chrono>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "generators.h"
#include "utils.h"



struct prob
{
	double probability;
	unsigned position;

};

bool probComparison(const prob& a, const prob& b) {
	if( a.probability < b.probability) {
		return true;
	}
	if(a.probability == b.probability) {
		return a.position < b.position;
	}
	return false;
}

double runAlgorithm(const unsigned maxGenerations, unsigned startingPopulationSize, const unsigned dimensions,
	const Function& function, const string& functionName, const int sample);

void crossOver(vector<vector<bool>>& startingPopulation, const vector<bool>& chromosome1, const vector<bool>& chromosome2,
	uniform_int_distribution<>& distribution, mt19937& gen, const unsigned& nodeLength);
void mutate(vector<bool>& chromosome, const double& chanceToMutate, uniform_real_distribution<double>& distribution, mt19937& gen);


double runAlgorithm(const unsigned maxGenerations, unsigned startingPopulationSize, const unsigned dimensions,
	const Function& function, const string& functionName, const int sample) {

	auto start = chrono::system_clock::now();

	constexpr unsigned precision = 5;
	const auto N = abs((function.lowerBound - function.upperBound) * pow(10, precision));
	const auto nodeLength = static_cast<unsigned>(ceil(log2(N)));
	const double a = pow(2, nodeLength) - 1;
	const double average = function.upperBound - function.lowerBound;
	unsigned L = nodeLength * dimensions;

	uniform_real_distribution<double> random01(0, 1.0);
	uniform_int_distribution<> randomGenerator(0, dimensions);
	uniform_int_distribution<> randomBool(0, 1);


	unsigned populationSize = startingPopulationSize;
	vector<vector<bool>> startingPopulation(populationSize, vector<bool>(nodeLength * dimensions));


	unsigned seed = chrono::system_clock::now().time_since_epoch().count() * sample;
	mt19937 gen(seed);

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
	double averagePopulationSize = 0;

	//Params
	double chanceToMutate = 1.5;
	constexpr double chanceToCrossOver = 0.8;
	//Very rough estimate
	constexpr int averageChromosomesToBeSelected = 120;
	

	chanceToMutate /= (nodeLength * dimensions);
	double bestResult = DBL_MAX;


	while (t < maxGenerations) {
		t++;
		//cout << "Started generation " << t << '\n';

		double totalFitness = 0;
		populationSize = startingPopulation.size();
		averagePopulationSize += populationSize;

		for (unsigned i = 0; i < populationSize; i++) {
			results[i] = evaluate(startingPopulation[i], nodeLength, function, a, average);
			if (results[i] < bestResult) {
				bestResult = results[i];
				generationsSinceLastImprovement = 0;
			}
			eval[i] = function.fitFunc(results[i], dimensions);
			totalFitness += eval[i];
			//cout <<"results["<<results[i]<<"] "<< "eval[" << i << "] " << eval[i] << ' ';
		}
		//cout << "TotalFitness " << totalFitness << '\n';
		generationsSinceLastImprovement++;
		if (generationsSinceLastImprovement > 250) {
			break;
		}

		vector<vector<bool>> selectedChromosomes;
		vector<prob> p(populationSize);
		vector<double> q(populationSize + 1, 0);

		deque<double> rVector;

		for (int i = 0; i < averageChromosomesToBeSelected; i++) {
			rVector.push_front(random01(gen));
		}

		//Individual sel.prob. + Accumulated sel.prob.
		for (unsigned i = 0; i < populationSize; i++) {
			p[i].probability = eval[i] / totalFitness;
			p[i].position = i;
		}

		sort(p.begin(), p.end(),probComparison );
		sort(rVector.begin(), rVector.end());
		bool finishedRandoms = false;
		for(unsigned i = 0;i<populationSize;i++) {
			if (finishedRandoms)
				break;

			q[i + 1] = q[i] + p[i].probability;
			if (q[i] < rVector[0] && rVector[0] <= q[i + 1]) {
				selectedChromosomes.push_back(startingPopulation[p[i].position]);
				rVector.pop_front();
				if (rVector.empty() == true)
					break;
				while (rVector[0] < q[i + 1]) {
					rVector.pop_front();
					if (rVector.empty() == true) {
						finishedRandoms = true;
						break;
					}
				}
				
				averageSelectedChromosomes++;
			}
		}

		//Generate new population
		unsigned chromosomesAdded = 0;

		//Kill-switch
		if (selectedChromosomes.empty() == true || selectedChromosomes.size() == 1) {
			t--;
			cout << "Selected none or 1";
			continue;
		}

		startingPopulation.clear();
		const unsigned numberOfSelectedChromosomes = selectedChromosomes.size();

		for (unsigned i = 0; i < numberOfSelectedChromosomes; i++) {
			startingPopulation.push_back(selectedChromosomes[i]);
			chromosomesAdded++;
		}

		for (unsigned i = 0; i < numberOfSelectedChromosomes - numberOfSelectedChromosomes % 2 - 1; i += 2) {
			if (random01(gen) < chanceToCrossOver) {
				crossOver(startingPopulation, selectedChromosomes[i],
					selectedChromosomes[i + 1], randomGenerator, gen, nodeLength);
				mutate(startingPopulation[chromosomesAdded], chanceToMutate, random01, gen);
				mutate(startingPopulation[chromosomesAdded + 1], chanceToMutate, random01, gen);
				chromosomesAdded += 2;
			}

		}


		//cout << "\nBest Result: " <<fixed<<setprecision(5)<< bestResult<<'\n';

	}

	auto end = chrono::system_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();

	string path = functionName + "_" + to_string(dimensions) + '_' + to_string(sample);
	ofstream fout(path);

	cout << path << '\n';

	fout << /*"Best Result: " << */fixed << setprecision(5) << bestResult << '\n';
	fout << /*"Average Selected " << */averageSelectedChromosomes / t << '\n';
	fout << averagePopulationSize / t << '\n';
	fout << /*"Ended at generation: " << */t << "\n";
	fout << /*"Duration: " << */duration << '\n';

	return bestResult;
}

void crossOver(vector<vector<bool>>& startingPopulation, const vector<bool>& chromosome1, const vector<bool>& chromosome2,
	uniform_int_distribution<>& distribution, mt19937& gen, const unsigned& nodeLength) {

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

void mutate(vector<bool>& chromosome, const double& chanceToMutate, uniform_real_distribution<double>& distribution, mt19937& gen) {
	for (auto var : chromosome) {
		if (distribution(gen) < chanceToMutate)
			var = !var;
	}
}
