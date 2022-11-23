#include "geneticAlgorithm.h"
#include "prelucrare.h"
using namespace std;


int main() {

	/*prelucrare();
	exit(0);*/

	int dimensions[] = {5, 10, 30};
	string func[] = { "deJong", "schwefel", "rastrigin", "michalewicz" };
	unsigned startingPopulationSize = 200;
	unsigned maxGenerations = 2000;
	//for (string s : func)
		for (int d : dimensions)
		{
			const Function deJong(-5.12, 5.12, deJongsFunction, deJongsFitness, d);
			const Function schwefel(-500, 500, schwefelsFunction, schwefelsFitness, d);
			const Function rastrigin(-5.12, 5.12, rastringsFunction, rastringsFitness, d);
			const Function michalewicz(0, PI, michalewiczsFunction, michalewiczsFitness, d);

			double averageTotal = 0;
			for (int i = 1; i <= 30; i++) {
				averageTotal += runAlgorithm(maxGenerations, startingPopulationSize, d, deJong, "deJong", i);
			}
			cout << "Average: " << averageTotal / 30 << '\n' << '\n';

			averageTotal = 0;
			for (int i = 1; i <= 30; i++) {
				averageTotal += runAlgorithm(maxGenerations, startingPopulationSize, d, schwefel, "Schwefel", i);
			}
			cout << "Average: " << averageTotal / 30 << '\n' << '\n';

			averageTotal = 0;
			for (int i = 1; i <= 30; i++) {
				averageTotal += runAlgorithm(maxGenerations, startingPopulationSize, d, rastrigin, "Rastrigin", i);
			}
			cout << "Average: " << averageTotal / 30 << '\n' << '\n';

			averageTotal = 0;
			for (int i = 1; i <= 30; i++) {
				averageTotal += runAlgorithm(maxGenerations, startingPopulationSize, d, michalewicz, "Michalewicz", i);
			}
			cout << "Average: " << averageTotal / 30 << '\n' << '\n';
		}

	
	return 0;

}
