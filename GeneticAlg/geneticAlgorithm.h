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



struct Prob
{
	double probability;
	unsigned position;

};

bool probComparison(const Prob& a, const Prob& b) {
	if (a.probability < b.probability) {
		return true;
	}
	if (a.probability == b.probability) {
		return a.position < b.position;
	}
	return false;
}

struct ModParams
{
	double chanceToMutate = 0.5; //		(0,1)		https://www.mathsisfun.com/numbers/sigma-calculator.html 1 to 10000 (chance)^n
	double chanceToCrossOver = 0.8;	//	(0.1,1]
	int averageChromosomesToBeSelected = 100;//		[5,100] :avgCTBS + avgCTBS*crossOverChance + Elites <=205
	unsigned numberOfElites = 5;	//  [0,100]
	double coolingConstant = 0.999; //  (0,1)
};

struct NonModParams
{
	const int sample;
	const unsigned maxGenerations = 2000;	// [10,2000]
	const unsigned startingPopulationSize = 200;	//	[50,200]
	const unsigned dimensions;	//{5,10,30}
	const string functionName;

};

//double runAlgorithm(const unsigned maxGenerations, unsigned startingPopulationSize, const unsigned dimensions,
//	const Function& function, const string& functionName, const int sample);

void crossOver(vector<vector<bool>>& startingPopulation, vector<bool>& chromosome1, vector<bool>& chromosome2,
               uniform_int_distribution<>& distribution, mt19937& gen, const unsigned& nodeLength);
void mutate(vector<bool>& chromosome, const double& chanceToMutate, uniform_real_distribution<double>& distribution, uniform_int_distribution
            <>& distribution2, mt19937& gen);

//unsigned averageMutated = 0;
//unsigned mutatedIndividuals = 0;


//double runAlgorithm(const unsigned maxGenerations, unsigned startingPopulationSize, const unsigned dimensions,
//	const Function& function, const string& functionName, const int sample) {
//
//
//	constexpr unsigned precision = 5;
//	const auto N = abs((function.lowerBound - function.upperBound) * pow(10, precision));
//	const auto nodeLength = static_cast<unsigned>(ceil(log2(N)));
//	const double a = pow(2, nodeLength) - 1;
//	const double average = function.upperBound - function.lowerBound;
//	unsigned L = nodeLength * dimensions;
//
//	uniform_real_distribution<double> random01(0, 1.0);
//	uniform_int_distribution<> randomGenerator(0, dimensions);
//	uniform_int_distribution<> randomGenerator2(0, nodeLength * dimensions - 1);
//	uniform_int_distribution<> randomBool(0, 1);
//
//
//	unsigned populationSize = startingPopulationSize;
//	vector<vector<bool>> startingPopulation(populationSize, vector<bool>(nodeLength * dimensions));
//
//
//	unsigned seed = chrono::system_clock::now().time_since_epoch().count() * sample;
//	mt19937 gen(seed);
//
//	for (unsigned i = 0; i < populationSize; i++) {
//		for (unsigned j = 0; j < L; j++) {
//			startingPopulation[i][j] = randomBool(gen);
//		}
//	}
//
//	vector<double> eval, results;
//
//	//Count Variables
//	unsigned t = 0;
//	unsigned generationsSinceLastImprovement = 0;
//	double averageSelectedChromosomes = 0;
//	double averagePopulationSize = 0;
//
//	//ModParams
//	double chanceToMutate = 0.99; // https://www.mathsisfun.com/numbers/sigma-calculator.html 1 to 10000 (1/chance)^n
//	double chanceToCrossOver = 0.5;
//	//Very rough estimate
//	int averageChromosomesToBeSelected = 100;
//	unsigned numberOfElites = 5;
//
//
//	//chanceToMutate /= (nodeLength * dimensions);
//	double bestResult = DBL_MAX;
//	double temperature = 1;
//	double coolingConstant = 0.9999;
//
//
//	auto start = chrono::system_clock::now();
//	while (t < maxGenerations) {
//		t++;
//
//		chanceToMutate *= temperature;
//		chanceToCrossOver *= (1.0 / temperature);
//
//		double totalFitness = 0;
//		populationSize = startingPopulation.size();
//		averagePopulationSize += populationSize;
//
//		//Evaluate
//		eval.resize(populationSize);
//		results.resize(populationSize);
//		for (unsigned i = 0; i < populationSize; i++) {
//			results[i] = evaluate(startingPopulation[i], nodeLength, function, a, average);
//			if (results[i] < bestResult) {
//				bestResult = results[i];
//				generationsSinceLastImprovement = 0;
//			}
//			eval[i] = function.fitFunc(results[i], dimensions);
//			totalFitness += eval[i];
//		}
//
//		//Shorter version
//		/*generationsSinceLastImprovement++;
//		if (generationsSinceLastImprovement > 250) {
//			break;
//		}*/
//
//		vector<vector<bool>> selectedChromosomes(averageChromosomesToBeSelected);
//		vector<Prob> p(populationSize);
//		vector<double> q(populationSize + 1, 0);
//
//		//Generate random numbers for wheel-of-fortune
//		deque<Prob> rVector;
//		for (int i = 0; i < averageChromosomesToBeSelected; i++) {
//			Prob temp;
//			temp.probability = random01(gen);
//			temp.position = i;
//			rVector.push_front(temp);
//		}
//
//		//Individual sel.prob.
//		for (unsigned i = 0; i < populationSize; i++) {
//			p[i].probability = eval[i] / totalFitness;
//			p[i].position = i;
//		}
//
//		sort(p.begin(), p.end(), probComparison);
//		sort(rVector.begin(), rVector.end(), probComparison);
//		bool finishedRandoms = false;
//		for (unsigned i = 0; i < populationSize - numberOfElites; i++) {
//			if (finishedRandoms)
//				break;
//
//			q[i + 1] = q[i] + p[i].probability;
//			if (q[i] < rVector[0].probability && rVector[0].probability <= q[i + 1]) {
//				selectedChromosomes[rVector[0].position] = startingPopulation[p[i].position];
//				rVector.pop_front();
//				if (rVector.empty() == true) {
//					finishedRandoms = true;
//					break;
//				}
//				while (rVector[0].probability <= q[i + 1]) {
//					selectedChromosomes[rVector[0].position] = startingPopulation[p[i].position];
//					rVector.pop_front();
//					if (rVector.empty() == true) {
//						finishedRandoms = true;
//						break;
//					}
//				}
//
//				//averageSelectedChromosomes++;
//			}
//		}
//
//		vector<vector<bool>> elites;
//		for (unsigned i = populationSize - numberOfElites; i < populationSize ; i++) {
//			if (finishedRandoms) {
//				elites.push_back(startingPopulation[p[i].position]);
//				continue;
//			}
//			q[i + 1] = q[i] + p[i].probability;
//			if (q[i] < rVector[0].probability && rVector[0].probability <= q[i + 1]) {
//				selectedChromosomes[rVector[0].position] = startingPopulation[p[i].position];
//				rVector.pop_front();
//				if (rVector.empty() == true) {
//					finishedRandoms = true;
//					continue;
//				}
//				while (rVector[0].probability <= q[i + 1]) {
//					selectedChromosomes[rVector[0].position] = startingPopulation[p[i].position];
//					rVector.pop_front();
//					if (rVector.empty() == true) {
//						finishedRandoms = true;
//						break;
//					}
//				}
//			}
//			else {
//				elites.push_back(startingPopulation[p[i].position]);
//			}
//		}
//
//		//Generate new population
//		unsigned chromosomesAdded = 0;
//
//		//Kill-switch
//		/*if (selectedChromosomes.empty() == true || selectedChromosomes.size() == 1) {
//			t--;
//			cout << "Selected none or 1";
//			continue;
//		}*/
//
//		startingPopulation.clear();
//		const unsigned numberOfSelectedChromosomes = selectedChromosomes.size();
//
//		/*for (unsigned i = 0; i < numberOfSelectedChromosomes; i++) {
//			startingPopulation.push_back(selectedChromosomes[i]);
//			chromosomesAdded++;
//		}*/
//		startingPopulation = selectedChromosomes;
//		for (unsigned i = 0; i < numberOfSelectedChromosomes - numberOfSelectedChromosomes % 2 - 1; i += 2) {
//			if (random01(gen) < chanceToCrossOver) {
//				crossOver(startingPopulation, selectedChromosomes[i],
//				          selectedChromosomes[i + 1], randomGenerator, gen, nodeLength);
//				mutate(startingPopulation[chromosomesAdded], chanceToMutate, random01,randomGenerator2, gen);
//				mutate(startingPopulation[chromosomesAdded + 1], chanceToMutate, random01,randomGenerator2, gen);
//				//mutatedIndividuals+=2;
//				chromosomesAdded += 2;
//			}
//		}
//
//		startingPopulation.insert(startingPopulation.begin(), elites.begin(), elites.end());
//		temperature *= coolingConstant;
//	}
//
//	auto end = chrono::system_clock::now();
//	auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
//
//	string path = functionName + "_" + to_string(dimensions) + '_' + to_string(sample);
//	ofstream fout(path);
//
//	cout << path << '\n';
//
//	fout << /*"Best Result: " << */fixed << setprecision(5) << bestResult << '\n';
//	fout << /*"Average Selected " << */averageSelectedChromosomes / t << '\n';
//	fout << averagePopulationSize / t << '\n';
//	fout << /*"Ended at generation: " << */t << "\n";
//	fout << /*"Duration: " << */duration << '\n';
//	/*cout << "Average mutations per individual: " << (double)averageMutated / (double)mutatedIndividuals << '\n';
//	averageMutated = 0;
//	mutatedIndividuals = 0;*/
//
//	return bestResult;
//}

void evolution(vector<vector<bool>>& startingPopulation, const unsigned dimensions, const Function& function,
	unsigned populationSize)
{
	constexpr unsigned precision = 5;
	const auto N = abs((function.lowerBound - function.upperBound) * pow(10, precision));
	const auto nodeLength = static_cast<unsigned>(ceil(log2(N)));
	const double a = pow(2, nodeLength) - 1;
	const double average = function.upperBound - function.lowerBound;
	unsigned L = nodeLength * dimensions;

	uniform_real_distribution<double> random01_real(0, 1.0);
	uniform_int_distribution<> randomGenerator(0, dimensions);
	uniform_int_distribution<> randomGenerator2(0, nodeLength * dimensions - 1);
	uniform_int_distribution<> randomBool(0, 1);

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 gen(seed);
	
	//for (unsigned i = 0; i < populationSize; i++) {
	//	for (unsigned j = 0; j < L; j++) {
	//		startingPopulation[i][j] = randomBool(gen);
	//	}
	//}

	vector<double> eval, results;
	double totalFitness = 0;
	

	//Wheel of Fortune
	eval.resize(populationSize);
	results.resize(populationSize);
	for (unsigned i = 0; i < populationSize; i++) {
		results[i] = evaluate(startingPopulation[i], nodeLength, function, a, average);
		eval[i] = function.fitFunc(results[i], dimensions);
		totalFitness += eval[i];
	}

	vector<double> prob;
	vector<double> accProb;

	prob.reserve(populationSize);
	accProb.reserve(populationSize + 1);
	accProb.push_back(0);

	for (int i = 0; i < populationSize; i++)
	{
		prob.push_back(eval[i] / totalFitness);
		accProb.push_back(prob[i]);
		accProb[i + 1] += accProb[i];
	}

	vector<vector<bool>> nextPopulation;

	for (int i = 0; i < populationSize; i++)
	{
		double r = random01_real(gen);
		for (int j = 0; j < populationSize; j++)
		{	
			if (accProb[j] < r && r <= accProb[j + 1])
			{
				nextPopulation.push_back(startingPopulation[j]);
				break;
			}
		}
	}

	//Mutation
	double changeToMutate = 0.015;
	for (int i = 0; i < populationSize; i++)
	{
		for (int j = 0; j < L; j++)
		{
			double random = random01_real(gen);
			if (random < changeToMutate)
				nextPopulation[i][j] = !nextPopulation[i][j];
		}
	}

	//CrossOver
	vector<pair<double, int>> chanceToCrossOver;
	chanceToCrossOver.reserve(populationSize);
	double crossOverProb = 0.25;
	for (int i = 0; i < populationSize; i++)
	{
		
		double random = random01_real(gen);
		chanceToCrossOver.push_back({ random, i });
	}
	sort(chanceToCrossOver.begin(), chanceToCrossOver.end());
	for (int i = 1; i < populationSize; i += 2)
	{
		if (chanceToCrossOver[i].first > crossOverProb)
			break;
		int cutPoint = randomGenerator2(gen);
		for (int j = 0; j <= cutPoint; j++)
		{
			swap(nextPopulation[chanceToCrossOver[i].second][j],
				nextPopulation[chanceToCrossOver[i - 1].second][j]);
		}
	}
	startingPopulation = nextPopulation;
}

double runAlgorithm(const unsigned dimensions, const Function& function, const string& functionName,
	const int sample)
{
	const int populationSize = 200;
	const int maxGenerations = 2000;

	constexpr unsigned precision = 5;
	const auto N = abs((function.lowerBound - function.upperBound) * pow(10, precision));
	const auto nodeLength = static_cast<unsigned>(ceil(log2(N)));
	unsigned L = nodeLength * dimensions;
	const double average = function.upperBound - function.lowerBound;
	const double a = pow(2, nodeLength) - 1;


	unsigned seed = chrono::system_clock::now().time_since_epoch().count() * sample;
	mt19937 gen(seed);
	uniform_int_distribution<> randomBool(0, 1);

	vector<vector<bool>> startingPopulation(populationSize, vector<bool>(nodeLength * dimensions));
	startingPopulation.reserve(populationSize);

	for (int i = 0; i < populationSize; i++)
	{
		startingPopulation[i].reserve(L);
		for (int j = 0; j < L; j++)
		{
			startingPopulation[i][j] = randomBool(gen);
		}
	}
	
	double bestResult = 0;
	vector<bool> bestPopulation;

	for (int g = 0; g < maxGenerations; g++)
	{
		for (int i = 0; i < populationSize; i++)
		{
			double current = function.fitFunc(evaluate(startingPopulation[i], nodeLength, function, a, average), dimensions);
			if (current > bestResult)
			{
				bestPopulation = startingPopulation[i];
				bestResult = current;
			}
		}
		//cout << bestResult << '\n';
		evolution(startingPopulation, dimensions, function, populationSize);
	}
	double result = evaluate(bestPopulation, nodeLength, function, a, average);
	cout << result << '\n';
	return result;
}

void crossOver(vector<vector<bool>>& startingPopulation, vector<bool>& chromosome1, vector<bool>& chromosome2,
               uniform_int_distribution<>& distribution, mt19937& gen, const unsigned& nodeLength) {

	const unsigned cutPoint = distribution(gen) * nodeLength;
	/*vector<bool> newChromosome1, newChromosome2;

	for (unsigned i = 0; i < cutPoint; i++) {
		newChromosome1.push_back(chromosome1[i]);
		newChromosome2.push_back(chromosome2[i]);
	}
	for (unsigned i = cutPoint; i < chromosome1.size(); i++) {
		newChromosome1.push_back(chromosome2[i]);
		newChromosome2.push_back(chromosome1[i]);
	}
	startingPopulation.push_back(newChromosome1);
	startingPopulation.push_back(newChromosome2);*/

	swap_ranges(chromosome1.begin(), chromosome1.begin() + cutPoint, chromosome2.begin());
	startingPopulation.push_back(chromosome1);
	startingPopulation.push_back(chromosome2);

}

void mutate(vector<bool>& chromosome, const double& chanceToMutate, uniform_real_distribution<double>& distribution, uniform_int_distribution
            <>& distribution2, mt19937& gen) {

	/*for (auto var : chromosome) {
		if (distribution(gen) < chanceToMutate) {
			var = !var;
			averageMutated++;
		}
	}*/

	while(distribution(gen) < chanceToMutate) {
		const int point = distribution2(gen);
		chromosome[point] = !chromosome[point];
		//averageMutated++;
	}
}