#pragma once
#include <iostream>
#include <fstream>
#include <string>



struct GlobalData {
	int SimulationTime;
	int SimulationStepTime;
	int Conductivity;
	int Alfa;
	int Tot;
	int InitialTemp;
	int Density;
	int SpecificHeat;
	int NodesNumber;
	int ElementsNumber;

	void showGlobalData();
};

struct ElementUniwersalny {
	double** matdEta;
	double** matdKsi;
	int nPkt;
	void init(int n);
};

struct node {
	double x, y;
};

struct element {
	int ID_wezlow[4];
};

struct grid
{
	node nodes[2];
	element elements[2];
};

GlobalData* readMesh(std::string fileSrc);
