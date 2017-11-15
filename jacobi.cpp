// Rozwiazywanie ukladu algebraicznych rownan liniowych postaci Ax = b metodą Jacobiego

#include <iostream>
#include <cmath>
#include "wektoryMacierze.h"
#include "jacobi.h"

double * jacobi(double ** A, double * b, int n, double * xn, double tolx, double tolf, double nmax){
	double * un = utworzWektor(n);
	double * estymator = utworzWektor(n);
	double * residuum = utworzWektor(n);
		
	double d;
	for(int i = 0; i < n; ++i){
		d = 1.0/A[i][i];
			b[i] *= d;
			for(int j = 0; j < n; ++j)
					A[i][j] *= d;
	}
	
	double s;
	for(int k = 1; k <= nmax; ++k){
		for(int i = 0; i < n; ++i){
			s = 0.0;
			for(int j = 0; j < n; ++j){
				if( j != i ){
						s += A[i][j]*xn[j];
				}
			}
			un[i] = b[i] - s;
		}
		for(int i = 0; i < n; ++i){
			estymator[i] = xn[i] - un[i];
			xn[i] = un[i];
			s = 0.0;
			for(int j = 0; j < n; ++j){
				s += A[i][j]*xn[j];
			}
			residuum[i] = b[i] - s;
		}
		/*std::cout << "*************(" << k << ")*************\n";
		std::cout << "x: " ; wypiszWektor(xn, n);
		std::cout << "e: " ; wypiszWektor(estymator, n);
		std::cout << "r: " ; wypiszWektor(residuum, n);*/
		
		bool flag = true;
		for(int i = 0; i < n; ++i){
			if(std::fabs(estymator[i]) >= tolx || std::fabs(residuum[i]) >= tolf)
				flag = false;
		}
		if(flag)
			break;
		
	}
	usunWektor(un);
	usunWektor(residuum);
	usunWektor(estymator);
	return xn;
}
