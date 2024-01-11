#include <iostream>
#include <vector>

#include "calkowanie.h"
#include "Mesh.h"
#include "ElementUniwersalny.h"
#include "RunType.h"


//Zadanie domowe
double funkcja_testowa_1(double x) {
	// Definicja funkcji, któr¹ chcemy ca³kowaæ
	return 5 * pow(x,2) + 3 * x + 6; // Przyk³adowa funkcja kwadratowa
}

double funkcja_testowa_2(double x, double y) {
	// Definicja funkcji, któr¹ chcemy ca³kowaæ
	return 5 * pow(x,2) * pow(y,2) + 3 * x * y + 6; // Przyk³adowa funkcja kwadratowa
}



int main()
{
#ifdef DEBUG
	std::cout << "Tryb Debug." << std::endl;
#endif

	Mesh siatka;
	siatka.readMeshFile("Test4_31_31_trapez.txt");

	ElementUniwersalny elementUniwersalny;
	elementUniwersalny.init(3);

	Mesh::calcTemperature(siatka, elementUniwersalny);

	//siatka.calcHForElements(elementUniwersalny);

	//siatka.calcHBCForElements(elementUniwersalny);

	//siatka.calcVectorPForElements(elementUniwersalny);

	//siatka.calcCForElements(elementUniwersalny);


	//std::cout << "przestrzen 1d 2 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss1d(funkcja_testowa_1, -1, 1, 2) << std::endl;

	//std::cout << "przestrzen 1d 3 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss1d(funkcja_testowa_1, -1, 1, 3) << std::endl;

	//std::cout << "przestrzen 2d 2 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss2d(funkcja_testowa_2, -1, 1, 2) << std::endl;

	//std::cout << "przestrzen 2d 3 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss2d(funkcja_testowa_2, -1, 1, 3) << std::endl;



	return 0;
}

