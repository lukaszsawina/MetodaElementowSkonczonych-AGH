#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>

#include "calkowanie.h"

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
};

struct ElementUniwersalny {
	double** matdEta;
	double** matdKsi;
	std::vector<double**> matNPktForEdges;
	
	int nPkt;
	void init(int n);
};

struct Node {
	int ID;
	double x, y;
	int BC = 0; //Czy jest warunkiem brzegowym i jakim
};

struct Element {
	int ID;
	int ID_wezlow[4];
	double** H;
	double** HBC;

	void calcH(double* x, double* y, ElementUniwersalny elUni);
	void calcHBC(double* x, double* y, int* BC, ElementUniwersalny elUni);
};

void calcVectorP(double* x, double* y, int* BC, ElementUniwersalny elUni);

struct Grid
{
	Node nodes[2];
	Element elements[2];
};

struct Mesh
{
	GlobalData* globalData;
	Node* nodes;
	Element* elements;

	void readMeshFile(std::string fileSrc);
	void showGlobalData();
	void showNodes();
	void showElements();

	void calcHForElements(ElementUniwersalny elUni);
	void calcHBCForElements(ElementUniwersalny elUni);

private:
	GlobalData* readMeshGlobalData(std::string fileSrc);
	Node* readMeshNodes(std::string fileSrc);
	Element* readMeshElements(std::string fileSrc);
};

