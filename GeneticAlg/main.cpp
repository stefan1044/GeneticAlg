#include "geneticAlgorithm.h"

using namespace std;


int main() {

	int dimensions[] = {5, 10, 30};
	string func[] = { "deJong", "schwefel", "rastrigin", "michalewicz" };

	//for (string s : func)
		for (int d : dimensions)
		{
			const Function deJong(-5.12, 5.12, deJongsFunction, deJongsFitness, d);
			const Function schwefel(-500, 500, schwefelsFunction, schwefelsFitness, d);
			const Function rastrigin(-5.12, 5.12, rastringsFunction, rastringsFitness, d);
			const Function michalewicz(0, PI, michalewiczsFunction, michalewiczsFitness, d);

			double averageTotal = 0;
			for (int i = 1; i <= 30; i++) {
				averageTotal += runAlgorithm(2000, 50, d, deJong, "deJong", i);
			}
			cout << "Average: " << averageTotal / 30 << '\n' << '\n';

			averageTotal = 0;
			for (int i = 1; i <= 30; i++) {
				averageTotal += runAlgorithm(2000, 50, d, schwefel, "Schwefel", i);
			}
			cout << "Average: " << averageTotal / 30 << '\n' << '\n';

			averageTotal = 0;
			for (int i = 1; i <= 30; i++) {
				averageTotal += runAlgorithm(2000, 50, d, rastrigin, "Rastrigin", i);
			}
			cout << "Average: " << averageTotal / 30 << '\n' << '\n';

			averageTotal = 0;
			for (int i = 1; i <= 30; i++) {
				averageTotal += runAlgorithm(2000, 50, d, michalewicz, "Michalewicz", i);
			}
			cout << "Average: " << averageTotal / 30 << '\n' << '\n';
		}

	//const Function deJong(-5.12, 5.12, deJongsFunction, deJongsFitness, 5);
	////const Function schwefel(-500, 500, schwefelsFunction, schwefelsFitness);
	////const Function rastrigin(-5.12, 5.12, rastringsFunction, rastringsFitness);
	////const Function michalewicz(0, PI, michalewiczsFunction, michalewiczsFitness);

	//double averageTotal = 0;
	//for (int i = 1; i <= 30; i++) {
	//	averageTotal += runAlgorithm(2000, 100, 5, deJong, "deJong", i);
	//}
	//
	//cout << averageTotal / 30;
	//return 0;

	
	return 0;

}
