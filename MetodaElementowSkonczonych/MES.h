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

struct Node {
	int ID;
	double x, y;
};

struct Element {
	int ID;
	int ID_wezlow[4];
};

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

private:
	GlobalData* readMeshGlobalData(std::string fileSrc);
	Node* readMeshNodes(std::string fileSrc);
	Element* readMeshElements(std::string fileSrc);
};

struct ElementUniwersalny {
	double** matdEta;
	double** matdKsi;
	int nPkt;
	void init(int n);
	double** H(double* x, double* y);
};

