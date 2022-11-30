#include "geneticAlgorithm.h"
#include "metaAlg.h"
#include "prelucrare.h"
using namespace std;


int main() {

	/*prelucrare();
	exit(0);*/

	int dimensions[] = { 5, 10, 30 };
	string func[] = { "deJong", "schwefel", "rastrigin", "michalewicz" };
	unsigned startingPopulationSize = 200;
	unsigned maxGenerations = 2000;

	ModParams params;
	/*const Function rastrigin(-5.12, 5.12, rastringsFunction, rastringsFitness, 10);
	metaAlg(params, rastrigin, 10, "rastrigin");*/



	const Function deJong5(-5.12, 5.12, deJongsFunction, deJongsFitness, 5);
	const Function schwefel5(-500, 500, schwefelsFunction, schwefelsFitness, 5);
	const Function rastrigin5(-5.12, 5.12, rastringsFunction, rastringsFitness, 5);
	const Function michalewicz5(0, PI, michalewiczsFunction, michalewiczsFitness, 5);

	const Function deJong10(-5.12, 5.12, deJongsFunction, deJongsFitness, 10);
	const Function schwefel10(-500, 500, schwefelsFunction, schwefelsFitness, 10);
	const Function rastrigin10(-5.12, 5.12, rastringsFunction, rastringsFitness, 10);
	const Function michalewicz10(0, PI, michalewiczsFunction, michalewiczsFitness, 10);

	const Function rastrigin30(-5.12, 5.12, rastringsFunction, rastringsFitness, 30);
	const Function michalewicz30(0, PI, michalewiczsFunction, michalewiczsFitness, 30);
	const Function deJong30(-5.12, 5.12, deJongsFunction, deJongsFitness, 30);
	const Function schwefel30(-500, 500, schwefelsFunction, schwefelsFitness, 30);

	unsigned iterations = 1;

	while(true) {
		metaAlg(michalewicz5, 5, "Michalewiczs");
	}

	


	//for (int d : dimensions)
	//{
	//	const Function deJong(-5.12, 5.12, deJongsFunction, deJongsFitness, d);
	//	const Function schwefel(-500, 500, schwefelsFunction, schwefelsFitness, d);
	//	const Function rastrigin(-5.12, 5.12, rastringsFunction, rastringsFitness, d);
	//	const Function michalewicz(0, PI, michalewiczsFunction, michalewiczsFitness, d);



	//	metaAlg(deJong, d, "deJongs");
	//	metaAlg(schwefel, d, "Schwefels");
	//	metaAlg(rastrigin, d, "Rastrigins");
	//	metaAlg(michalewicz, d, "Michalewiczs");

	//	/*double averageTotal = 0;
	//	for (int i = 1; i <= 30; i++) {
	//		averageTotal += runAlgorithm(maxGenerations, startingPopulationSize, d, deJong, "deJong", i,params);
	//	}
	//	cout << "Average: " << averageTotal / 30 << '\n' << '\n';

	//	averageTotal = 0;
	//	for (int i = 1; i <= 30; i++) {
	//		averageTotal += runAlgorithm(maxGenerations, startingPopulationSize, d, schwefel, "Schwefel", i,params);
	//	}
	//	cout << "Average: " << averageTotal / 30 << '\n' << '\n';

	//	averageTotal = 0;
	//	for (int i = 1; i <= 30; i++) {
	//		averageTotal += runAlgorithm(maxGenerations, startingPopulationSize, d, rastrigin, "Rastrigin", i,params);
	//	}
	//	cout << "Average: " << averageTotal / 30 << '\n' << '\n';

	//	averageTotal = 0;
	//	for (int i = 1; i <= 30; i++) {
	//		averageTotal += runAlgorithm(maxGenerations, startingPopulationSize, d, michalewicz, "Michalewicz", i,params);
	//	}
	//	cout << "Average: " << averageTotal / 30 << '\n' << '\n';*/
	//}


	return 0;

}