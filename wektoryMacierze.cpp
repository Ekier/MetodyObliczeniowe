#include <iostream>
#include "wektoryMacierze.h"

double * utworzWektor(int n){
	double * w = new double[n];
	return w;
}

void usunWektor(double * w){
	delete [] w;
}

void wypiszWektor(double * w, int n){
	for(int i = 0; i < n; ++i)
		std::cout << *(w + i) << "\t";
	std::cout << std::endl;
}

double ** utworzMacierz(int n){
	double ** macierz = new double *[n];
	for (int i = 0; i < n; ++i)
		macierz[i] = new double[n];
	return macierz;
}

void usunMacierz(double ** A, int n){
	for (int i = 0; i < n; ++i)
		delete A[i];
	delete [] A;
}
