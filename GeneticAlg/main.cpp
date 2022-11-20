#include "geneticAlgorithm.h"

using namespace std;


int main() {

	const Function deJong(-5.12, 5.12, deJongsFunction, deJongsFitness);
	const Function schwefel(-500, 500, schwefelsFunction, schwefelsFitness);
	const Function rastrigin(-5.12, 5.12, rastringsFunction, rastringsFitness);
	const Function michalewicz(0, PI, michalewiczsFunction, michalewiczsFitness);

	double averageTotal = 0;
	for (int i = 0; i < 30; i++) {
		averageTotal+=runAlgorithm(2000, 100, 30, rastrigin);
	}

	cout << averageTotal / 30;
	return 0;

}