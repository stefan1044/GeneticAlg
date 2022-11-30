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
	double chanceToMutate = 0.75; //		(0,1)		https://www.mathsisfun.com/numbers/sigma-calculator.html 1 to 10000 (chance)^n
	double chanceToCrossOver = 0.9;	//	(0.1,1]
	int averageChromosomesToBeSelected = 100;//		[5,100] :avgCTBS + avgCTBS*crossOverChance + Elites <=205
	unsigned numberOfElites = 7;	//  [0,100]
};

struct NonModParams
{
	const int sample;
	const unsigned maxGenerations = 2000;	// [10,2000]
	const unsigned startingPopulationSize = 200;	//	[50,200]
	const unsigned dimensions;	//{5,10,30}
	const string functionName;

};

double runAlgorithm(const unsigned maxGenerations, unsigned startingPopulationSize, const unsigned dimensions,
                    const Function& function, const string& functionName, const int sample, const ModParams& parameters);

void crossOver(vector<vector<bool>>& startingPopulation, vector<bool>& chromosome1, vector<bool>& chromosome2,
	uniform_int_distribution<>& distribution, mt19937& gen, const unsigned& nodeLength);
void mutate(vector<bool>& chromosome, const double& chanceToMutate, uniform_real_distribution<double>& distribution, uniform_int_distribution
	<>& distribution2, mt19937& gen);

//unsigned averageMutated = 0;
//unsigned mutatedIndividuals = 0;
//unsigned genFoundBest = 0;


double runAlgorithm(const unsigned maxGenerations, unsigned startingPopulationSize, const unsigned dimensions,
                    const Function& function, const string& functionName, const int sample, const ModParams& parameters) {

	auto start = chrono::system_clock::now();

	constexpr unsigned precision = 5;
	const auto N = abs((function.lowerBound - function.upperBound) * pow(10, precision));
	const auto nodeLength = static_cast<unsigned>(ceil(log2(N)));
	const double a = pow(2, nodeLength) - 1;
	const double average = function.upperBound - function.lowerBound;
	unsigned L = nodeLength * dimensions;

	uniform_real_distribution<double> random01(0, 1.0);
	uniform_int_distribution<> randomGenerator(0, dimensions);
	uniform_int_distribution<> randomGenerator2(0, nodeLength * dimensions - 1);
	uniform_int_distribution<> randomBool(0, 1);


	unsigned populationSize = startingPopulationSize;
	vector<vector<bool>> startingPopulation(populationSize, vector<bool>(nodeLength * dimensions));


	unsigned seed = chrono::system_clock::now().time_since_epoch().count() * sample;
	mt19937 gen(seed);

	for (unsigned i = 0; i < populationSize; i++) {
		for (unsigned j = 0; j < L; j++) {
			startingPopulation[i][j] = randomBool(gen);
		}
	}

	vector<double> eval, results;

	//Count Variables
	unsigned t = 0;
	unsigned generationsSinceLastImprovement = 0;
	double averageSelectedChromosomes = 0;
	double averagePopulationSize = 0;

	//ModParams
	double chanceToMutate = parameters.chanceToMutate; // https://www.mathsisfun.com/numbers/sigma-calculator.html 1 to 10000 (1/chance)^n
	double chanceToCrossOver = parameters.chanceToCrossOver;
	int averageChromosomesToBeSelected = parameters.averageChromosomesToBeSelected;
	unsigned numberOfElites = parameters.numberOfElites;


	//chanceToMutate /= (nodeLength * dimensions);
	double bestResult = DBL_MAX;
	//double temperature = 1;
	//double coolingConstant = 0.9999;

	vector<bool> bestChromosome;

	while (t < maxGenerations) {
		t++;


		double totalFitness = 0;
		populationSize = startingPopulation.size();
		averagePopulationSize += populationSize;

		//Evaluate
		eval.resize(populationSize);
		results.resize(populationSize);
		for (unsigned i = 0; i < populationSize; i++) {
			results[i] = evaluate(startingPopulation[i], nodeLength, function, a, average);
			if (results[i] < bestResult) {
				bestResult = results[i];
				generationsSinceLastImprovement = 0;
				bestChromosome = startingPopulation[i];
			}
			eval[i] = function.fitFunc(results[i], dimensions);
			totalFitness += eval[i];
		}

		//Shorter version
		/*generationsSinceLastImprovement++;
		if (generationsSinceLastImprovement > 250) {
			break;
		}*/

		vector<vector<bool>> selectedChromosomes(averageChromosomesToBeSelected);
		vector<Prob> p(populationSize);
		vector<double> q(populationSize + 1, 0);

		//Generate random numbers for wheel-of-fortune
		deque<Prob> rVector;
		for (int i = 0; i < averageChromosomesToBeSelected; i++) {
			Prob temp;
			temp.probability = random01(gen);
			temp.position = i;
			rVector.push_front(temp);
		}

		//Individual sel.prob.
		for (unsigned i = 0; i < populationSize; i++) {
			p[i].probability = eval[i] / totalFitness;
			p[i].position = i;
		}

		sort(p.begin(), p.end(), probComparison);
		sort(rVector.begin(), rVector.end(), probComparison);
		bool finishedRandoms = false;
		for (unsigned i = 0; i < populationSize - numberOfElites; i++) {
			if (finishedRandoms)
				break;

			q[i + 1] = q[i] + p[i].probability;
			if (q[i] < rVector[0].probability && rVector[0].probability <= q[i + 1]) {
				selectedChromosomes[rVector[0].position] = startingPopulation[p[i].position];
				rVector.pop_front();
				if (rVector.empty() == true) {
					finishedRandoms = true;
					break;
				}
				while (rVector[0].probability <= q[i + 1]) {
					selectedChromosomes[rVector[0].position] = startingPopulation[p[i].position];
					rVector.pop_front();
					if (rVector.empty() == true) {
						finishedRandoms = true;
						break;
					}
				}

				//averageSelectedChromosomes++;
			}
		}

		vector<vector<bool>> elites;
		for (unsigned i = populationSize - numberOfElites; i < populationSize; i++) {
			if (finishedRandoms) {
				elites.push_back(startingPopulation[p[i].position]);
				continue;
			}
			q[i + 1] = q[i] + p[i].probability;
			if (q[i] < rVector[0].probability && rVector[0].probability <= q[i + 1]) {
				selectedChromosomes[rVector[0].position] = startingPopulation[p[i].position];
				rVector.pop_front();
				if (rVector.empty() == true) {
					finishedRandoms = true;
					continue;
				}
				while (rVector[0].probability <= q[i + 1]) {
					selectedChromosomes[rVector[0].position] = startingPopulation[p[i].position];
					rVector.pop_front();
					if (rVector.empty() == true) {
						finishedRandoms = true;
						break;
					}
				}
			}
			else {
				elites.push_back(startingPopulation[p[i].position]);
			}
		}

		//Generate new population
		unsigned chromosomesAdded = 0;

		//Kill-switch
		/*if (selectedChromosomes.empty() == true || selectedChromosomes.size() == 1) {
			t--;
			cout << "Selected none or 1";
			continue;
		}*/

		startingPopulation.clear();
		const unsigned numberOfSelectedChromosomes = selectedChromosomes.size();

		/*for (unsigned i = 0; i < numberOfSelectedChromosomes; i++) {
			startingPopulation.push_back(selectedChromosomes[i]);
			chromosomesAdded++;
		}*/
		startingPopulation = selectedChromosomes;
		for (unsigned i = 0; i < numberOfSelectedChromosomes - numberOfSelectedChromosomes % 2 - 1; i += 2) {
			if (random01(gen) < chanceToCrossOver) {
				crossOver(startingPopulation, selectedChromosomes[i],
					selectedChromosomes[i + 1], randomGenerator, gen, nodeLength);
				mutate(startingPopulation[chromosomesAdded], chanceToMutate, random01, randomGenerator2, gen);
				mutate(startingPopulation[chromosomesAdded + 1], chanceToMutate, random01, randomGenerator2, gen);
				//mutatedIndividuals+=2;
				chromosomesAdded += 2;
			}
		}

		startingPopulation.insert(startingPopulation.begin(), elites.begin(), elites.end());

		
		

	}

	//Hill-Climbing

	double current = bestResult;

	bool local = false;
	do
	{
		double minim = DBL_MAX;
		vector<bool> bitsminim;
		for (unsigned i = 0; i < L; i++)
		{
			bestChromosome[i] = !bestChromosome[i];
			double aux = evaluate(bestChromosome, nodeLength, function, a, average);
			if (aux < current && aux < minim)
			{
				minim = aux;
				bitsminim = bestChromosome;
			}
			bestChromosome[i] = !bestChromosome[i];
		}
		if (minim == DBL_MAX)
			local = true;
		else
		{
			bestChromosome = bitsminim;
			current = minim;
		}
	} while (!local);

	if (current < bestResult)
		bestResult = current;

	auto end = chrono::system_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();

	/*string path = functionName + "_" + to_string(dimensions) + '_' + to_string(sample);
	ofstream fout(path);*/

	//cout << path << '\n';

	//fout << /*"Best Result: " << */fixed << setprecision(5) << bestResult << '\n';
	//fout << /*"Average Selected " << */averageSelectedChromosomes / t << '\n';
	//fout << averagePopulationSize / t << '\n';
	//fout << /*"Ended at generation: " << */t << "\n";
	//fout << /*"Duration: " << */duration << '\n';


	/*cout << "Average mutations per individual: " << (double)averageMutated / (double)mutatedIndividuals << '\n';
	averageMutated = 0;
	mutatedIndividuals = 0;*/
	//cout << "Gen found best: " << genFoundBest << "\n";

	return bestResult;
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

	while (distribution(gen) < chanceToMutate) {
		const int point = distribution2(gen);
		chromosome[point] = !chromosome[point];
		//averageMutated++;
	}
}
