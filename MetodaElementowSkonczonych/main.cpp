#include <iostream>
#include <vector>

#include "calkowanie.h"
#include "MES.h"



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
	//Mesh siatka;
	//siatka.readMeshFile("Test1_4_4.txt");
	//siatka.showGlobalData();
	//siatka.showNodes();
	//siatka.showElements();

	ElementUniwersalny elementUniwersalny;
	elementUniwersalny.init(2);

	//siatka.calcHForNodes(elementUniwersalny);

	//siatka.calcHBCForElements(elementUniwersalny);


	double x[4] = { 0.1, 0.05469, 0.06239, 0.1 };
	double y[4] = { 0.005, 0.005, -0.03261, -0.04031 };

	int bc[4] = { 1,1,0,1 };

	calcVectorP(x, y, bc, elementUniwersalny);

	//calcHBC(x, y, bc, elementUniwersalny);

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

