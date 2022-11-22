#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>

using namespace std;

int main()
{
	int dim[] = { 5, 10/*, 30*/ };
	string func[] = { "deJong", "schwefel", "michalewicz", "rastrigin" };

	for (string f : func)
		for (int d : dim)
		{
			string tip = f + "_" + to_string(d);
			cout << tip << '\n';

			double medie = 0;
			double stdev = 0;
			double valoare[350];
			double avgsel;
			double endgen;
			double minim = 1e9;
			double maxim = -1e9;
			double durata_total = 0;
			double durata_minima = 1e9, durata_maxima = 0;

			double avgselTotal = 0;
			double endgenTotal = 0;

			ofstream fout("aa_" + tip);

			for (int i = 1; i <= 30; i++)
			{
				string path = f + "_" + to_string(d) + "_" + to_string(i);
				ifstream fin(path);
				double durata = 0;
				fin >> valoare[i];
				fin >> avgsel;
				fin >> endgen;
				fin >> durata;

				medie += valoare[i];

				if (valoare[i] < minim)
					minim = valoare[i];
				if (valoare[i] > maxim)
					maxim = valoare[i];

				durata /= 1000;

				durata_total += durata;
				if (durata > durata_maxima)
					durata_maxima = durata;
				if (durata < durata_minima)
					durata_minima = durata;

				avgselTotal += avgsel;
				endgenTotal += endgen;

				//fout << durata << '\n';
			}

			double durata_medie = durata_total / 30;
			fout << "Durata medie: " << durata_medie << '\n';
			fout << "Durata minima:" << durata_minima << '\n';
			fout << "Durata maxima:" << durata_maxima << '\n';

			medie /= 30;

			fout << "Medie: " << setprecision(5) << medie << '\n';

			for (int i = 1; i <= 30; i++)
				stdev += (valoare[i] - medie) * (valoare[i] - medie);

			stdev = sqrt(stdev / 30);
			fout << "St Dev: " << setprecision(5) << stdev << '\n';

			fout << "Minim: " << minim << '\n';
			fout << "Maxim: " << maxim << '\n';

			double avgselMedie = avgselTotal / 30;
			double endgenMedie = endgenTotal / 30;

			fout << "Avg. Sel. average: " << avgselMedie << '\n';
			fout << "End Gen average: " << endgenMedie << '\n';
		}
	return 0;
}