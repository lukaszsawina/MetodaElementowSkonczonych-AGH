#include <iostream>
#include <vector>

#include "calkowanie.h"
#include "MES.h"



//Zadanie domowe
double funkcja_testowa_1(double x) {
	// Definicja funkcji, kt�r� chcemy ca�kowa�
	return 5 * pow(x,2) + 3 * x + 6; // Przyk�adowa funkcja kwadratowa
}

double funkcja_testowa_2(double x, double y) {
	// Definicja funkcji, kt�r� chcemy ca�kowa�
	return 5 * pow(x,2) * pow(y,2) + 3 * x * y + 6; // Przyk�adowa funkcja kwadratowa
}


int main()
{

	//Mesh siatka;
	//siatka.readMeshFile("Test1_4_4.txt");
	//siatka.showGlobalData();
	//siatka.showNodes();
	//siatka.showElements();

	//std::cout << "przestrzen 1d 2 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss1d(funkcja_testowa_1, -1, 1, 2) << std::endl;

	//std::cout << "przestrzen 1d 3 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss1d(funkcja_testowa_1, -1, 1, 3) << std::endl;

	//std::cout << "przestrzen 2d 2 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss2d(funkcja_testowa_2, -1, 1, 2) << std::endl;

	//std::cout << "przestrzen 2d 3 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss2d(funkcja_testowa_2, -1, 1, 3) << std::endl;

	ElementUniwersalny elementUniwersalny;

	double x[4] = { 0, 0.025, 0.025, 0 };
	double y[4] = { 0, 0, 0.025, 0.025 };

	elementUniwersalny.init(2);

	elementUniwersalny.H(x, y);

	return 0;
}

