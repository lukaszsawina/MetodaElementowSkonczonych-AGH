#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>

#include "calkowanie.h"
#include "UkladyRownan.h"
#include "GlobalData.h"
#include "Element.h"
#include "ElementUniwersalny.h"
#include "Node.h"

struct Mesh
{
	GlobalData* globalData;
	Node* nodes;
	Element* elements;

	void readMeshFile(std::string fileSrc);
	void showGlobalData();
	void showNodes();
	void showElements();

	void calcHForElements(const ElementUniwersalny& elUni);
	void calcHBCForElements(const ElementUniwersalny& elUni);
	void calcVectorPForElements(const ElementUniwersalny& elUni);

	static double* calcTemperatureForElements(Mesh& mesh, const ElementUniwersalny& elUni);

	~Mesh();

private:
	GlobalData* readMeshGlobalData(std::string fileSrc);
	Node* readMeshNodes(std::string fileSrc);
	Element* readMeshElements(std::string fileSrc);
};