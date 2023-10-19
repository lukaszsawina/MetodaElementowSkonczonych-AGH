#pragma once
#include <iostream>
#include "Interpolacja.h"

double Calkowanie_metoda_prostokatow(double (*f)(double x), double a, double b, int n);
double Gauss1d(double (*f)(double x), double up, double down, int n);
double Gauss2d(double (*f)(double x, double y), double up, double down, int n);