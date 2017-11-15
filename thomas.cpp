// Program realizuje algorytm Thomasa dla macierzy trójdiagonalnej NxN

#include <iostream>
#include "wektoryMacierze.h"
#include "thomas.h"

double * dekompozycjaLUmacierz(double * l, double * d, double * u, int n){
	for(int i = 1; i < n; ++i)
		*(d + i) -= (*(l + i) * *(u + i - 1))/ *(d + i - 1);
	return d;
}

double * dekompozycjaLUwektor(double * b, double * eta, double * l, int n){
	for(int i = 1; i < n; ++i)
		*(b + i) -= (*(l + i) * *(b + i - 1))/ *(eta + i - 1);
	return b;
}

double * rozwiazUklad(double * eta, double * r, double * u, int n){
	*(r + n - 1) /= *(eta + n - 1);
	for(int i = n - 2; i >= 0; --i)
		*(r + i) = (*(r + i) - *(u + i) * *(r + i + 1)) / *(eta + i);
	return r;
}
