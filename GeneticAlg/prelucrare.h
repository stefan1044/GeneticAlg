#pragma once
#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>


void prelucrare() {
	int dim[] = { 5, 10, 30 };
	std::string func[] = { "deJong", "schwefel", "michalewicz", "rastrigin" };

	for (std::string f : func)
		for (int d : dim)
		{
			std::string tip = f + "_" + std::to_string(d);
			std::cout << tip << '\n';

			double medie = 0;
			double stdev = 0;
			double valoare[350];
			double avgsel;
			double averagePopulationSize;
			double endgen;
			double minim = 1e9;
			double maxim = -1e9;
			double durata_total = 0;
			double durata_minima = 1e9, durata_maxima = 0;

			double avgselTotal = 0;
			double averagePopulationSizeTotal = 0;
			double endgenTotal = 0;

			std::ofstream fout("aa_" + tip);

			for (int i = 1; i <= 30; i++)
			{
				std::string path = f + "_" + std::to_string(d) + "_" + std::to_string(i);
				std::ifstream fin(path);
				double durata = 0;
				fin >> valoare[i];
				fin >> avgsel;
				fin >> averagePopulationSize;
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
				averagePopulationSizeTotal += averagePopulationSize;
				endgenTotal += endgen;

				//fout << durata << '\n';
			}

			double durata_medie = durata_total / 30;
			fout << "Durata medie: " << durata_medie << '\n';
			fout << "Durata minima:" << durata_minima << '\n';
			fout << "Durata maxima:" << durata_maxima << '\n';

			medie /= 30;

			fout << "Medie: " << std::setprecision(5) << medie << '\n';

			for (int i = 1; i <= 30; i++)
				stdev += (valoare[i] - medie) * (valoare[i] - medie);

			stdev = sqrt(stdev / 30);
			fout << "St Dev: " << std::setprecision(5) << stdev << '\n';

			fout << "Minim: " << minim << '\n';
			fout << "Maxim: " << maxim << '\n';

			double avgselMedie = avgselTotal / 30;
			double endgenMedie = endgenTotal / 30;
			double averagePopulationSizeMedie = averagePopulationSizeTotal / 30;

			fout << "Avg. Sel. average: " << avgselMedie << '\n';
			fout << "End Gen average: " << endgenMedie << '\n';
			fout << "Average Population Size: " << averagePopulationSizeMedie;
		}

}
