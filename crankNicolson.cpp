// Rozwiązanie zagadnienia opisującego transport ciepła z warunkiem początkowym i brzegowym

#include <iostream>
#include <fstream>
#include <cmath>
#include "wektoryMacierze.h"
#include "thomas.h"
#include "jacobi.h"
#include "calerf.h"

// Parametry
double t0 = 0.0, tmax = 1.0, b = 0.1, D = 1.0;
double a = 6.0*std::sqrt(D*tmax); // wartość do wyznaczenia przedziału skończonego [-a, a] do obliczeń numerycznych
double lambda = 1.0;
int ileH = 200; // liczba krokow przestrzennych
double h = 2.0*a/(ileH-1); // krok przestrzenny
double dt = lambda*h*h; // krok czasowy
int ileDt = static_cast<int>((tmax-t0)/dt + 1); // liczba poziomow czasowych

// zwraca wartosc analityczna rozwiazania zagadnienia w punktcie (x,t)
// korzysta z funkcji ERFCL z pakietu calerf
double rozwiazanieAnalityczne(double x, double t){
	return 0.5*std::exp(D*t/(b*b) - x/b) * ERFCL((2*D*t/b - x)/(2*std::sqrt(D*t)));
}

// Dyskretyzacja metodą Cranka-Nicholson z zastosowaniem algorytmu Thomasa do rozwiazania ukladu rownan
// n - liczba kroków przestrzennych
// m - liczba poziomow czasowych
// zwraca maksymalna wartosc bledu dla ostatniego poziomu czasowego
double crankNicolsonThomas(const int n, const int m){
	// utworzenie wektorów diagonalnych macierzy A w ukladzie rownan Ax=b
	double * l = utworzWektor(n);
	double * d = utworzWektor(n);
	double * u = utworzWektor(n-1);

	double ** XXX = new double *[m]; // macierz rozwiazan zagadnienia, rozwiazania na poszczegolnych poziomach czasowych
	for (int i = 0; i < m; i++)
		XXX[i] = new double[n];


	for(int i = 0; i < n; ++i){ // uzupelnienie wektorow diagonalnych macierzy A w ukladzie rownan Ax=b
		l[i]= lambda/2.0;
		d[i]= -(1.0+lambda);
		u[i]= lambda/2.0;
	}

	// warunki brzegowe
	l[n-1]= 0.0;
	d[0]= 1.0;
	d[n-1]= 1.0;
	u[0]= 0.0;

	dekompozycjaLUmacierz(l, d, u, n); // pierwsza część algorytmu Thomasa
	
	double * x = XXX[0]; // wektor x w ukladzie rownan Ax=b, rozwiazanie na poczatkowym poziomie czasowym
	double xx = -a; // poczatek zakresu dla zmiennej przestrzennej
	for(int i = 0; i < n; ++i){ // inicjacja wektora rozwiazan zgodnie z warunkiem poczatkowym
		if(xx < 0)
			x[i] = 0;
		else
			x[i] = std::exp(-xx/b); 
		xx += h;
	}

	double tt = t0; // aktualny czas
	double maxBlad;
	FILE * plik1;
	plik1 = fopen("crankNicolsonThomasCzasBlad.txt", "w");
	double * bb = utworzWektor(n); // wektor b w ukladzie rownan Ax=b
	for(int k = 1; k < m; ++k){ // kolejne poziomy czasowe
		for(int i = 1; i < n-1; ++i){ // aktualizacja wektora b w ukladzie rownan Ax=b
			bb[i] = -((lambda/2.0)*x[i-1] + (1.0-lambda)*x[i] + (lambda/2.0)*x[i+1]);
		}

		// warunki brzegowe
		bb[0] = 0.0;
		bb[n-1] = 0.0;
		
		FILE * plik11;
		int krok = 270;
		if(k == krok){
			plik11= fopen("crankNicolsonThomasRozw270.txt", "w");
		}
		double xx = -a + h;
		maxBlad = 0.0;
		for(int i = 1; i < n-1; ++i){ // wyznaczenie maksymalnego bezwzglednego bledu
			double aktualBlad = std::fabs(x[i]-rozwiazanieAnalityczne(xx, tt));
			if(k == krok){
				fprintf(plik11, "%15.8f %15.8f %15.8f %15.8f %15.8f\n", xx, tt, x[i], rozwiazanieAnalityczne(xx, tt), aktualBlad);
			}
			if(aktualBlad > maxBlad){
				maxBlad = aktualBlad;
			}
			xx += h;
		}

		if(k == krok){
			fclose(plik11);
		}
		
		fprintf(plik1, "%d %15.8f %15.8f\n", k, tt, maxBlad);

		// rozwiazanie ukladu rownan Ax=b na aktualnym poziomie czasowym, wynik zapisany w wektorze bb
		dekompozycjaLUwektor(bb, d, l, n);
		rozwiazUklad(d, bb, u ,n);

		x = XXX[k]; // wskaż na kolejny poziom czasowy w macierzy rozwiazan
		for(int j = 0; j < n; ++j){ // kopiowanie aktualnego rozwiazania do macierzy rozwiazan
			x[j] = bb[j];
		}
		tt += dt;
	}
	fclose(plik1);
	
	FILE *plik2;
	plik2 = fopen("crankNicolsonThomas.txt", "w");
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < m; ++j){
			fprintf(plik2, "%15.8f", XXX[j][i]);
		}
		fprintf(plik2, "\n");
	}
	fclose(plik2);

	usunWektor(l);
	usunWektor(d);
	usunWektor(u);
	usunWektor(bb);
	usunMacierz(XXX, m);
	return maxBlad;
}

// Dyskretyzacja metodą Cranka-Nicholson z zastosowaniem metody iteracyjnej Jacobiego do rozwiazania ukladu rownan
// n - liczba kroków przestrzennych
// m - liczba poziomow czasowych
// zwraca maksymalna wartosc bledu dla ostatniego poziomu czasowego
double crankNicolsonJacobi(const int n, const int m){
	double ** A = utworzMacierz(n); // macierz A w ukladzie rownan Ax=b
	double * bb = utworzWektor(n); // wektor b w ukladzie rownan Ax=b
	
	double ** XXX = new double *[m]; // macierz rozwiazan zagadnienia, rozwiazania na poszczegolnych poziomach czasowych
	for (int i = 0; i < m; i++)
		XXX[i] = new double[n];

	double * x = XXX[0]; // wektor x w ukladzie rownan Ax=b, rozwiazanie na poczatkowym poziomie czasowym
	double xx = -a; // poczatek zakresu dla zmiennej przestrzennej
	for(int i = 0; i < n; ++i){ // inicjacja wektora rozwiazan zgodnie z warunkiem poczatkowym
		if(xx < 0)
			x[i] = 0.0;
		else
			x[i] = std::exp(-xx/b); 
		xx += h;
	}

	double * x0 = utworzWektor(n); // wektor poczatkowego przyblizenia rozwiazania dla metody jacobiego
	for(int i = 0; i < n; ++i){
			x0[i] = 1.0;
	}
	double tt = t0; // aktualny czas
	double maxBlad;
	FILE * plik1;
	plik1 = fopen("crankNicolsonJacobiCzasBlad.txt", "w");
	for(int k = 1; k < m; ++k){ // kolejne poziomy czasowe
		// wypelnienie macierzy A
		for(int i = 0; i < n; ++i){
			for(int j = 0; j < n; ++j){
				A[i][j] = 0.0;
			}
		}
		for(int i = 1; i < n-1; ++i){
			A[i][i] = -(1.0 + lambda);
			A[i][i-1] = A[i][i+1] = lambda/2.0;
		}
		// warunki brzegowe
		A[0][0] = 1.0;
		A[n-1][n-1] = 1.0;

		for(int i = 1; i < n-1; ++i){ // aktualizacja wektora b w ukladzie rownan Ax=b
			bb[i] = -((lambda/2.0)*x[i-1] + (1.0-lambda)*x[i] + (lambda/2.0)*x[i+1]);
		}
		// warunki brzegowe
		bb[0] = 0.0;
		bb[n-1] = 0.0;

		FILE * plik11;
		int krok = 270;
		if(k == krok){
			plik11= fopen("crankNicolsonJacobiRozw270.txt", "w");
		}
		
		double xx = -a + h;
		maxBlad = 0.0;
		for(int i = 1; i < n-1; ++i){ // wyznaczenie maksymalnego bezwzglednego bledu
			double aktualBlad = std::fabs(x[i]-rozwiazanieAnalityczne(xx, tt));
			if(k == krok){
				fprintf(plik11, "%15.8f %15.8f %15.8f %15.8f %15.8f\n", xx, tt, x[i], rozwiazanieAnalityczne(xx, tt), aktualBlad);
			}
			if(aktualBlad > maxBlad){
				maxBlad = aktualBlad;
			}
			xx += h;
		}
		if(k == krok){
			fclose(plik11);
		}

		fprintf(plik1, "%d %15.8f %15.8f\n", k, tt, maxBlad);

		// rozwiazanie ukladu rownan Ax=b na aktualnym poziomie czasowym, wynik zapisany w wektorze x0
		x0 = jacobi(A, bb, n, x0, 1.0e-12, 1.0e-12, 69);

		x = XXX[k]; // wskaż na kolejny poziom czasowy w macierzy rozwiazan
		for(int j = 0; j < n; ++j){ // kopiowanie aktualnego rozwiazania do macierzy rozwiazan
			x[j] = x0[j];
		}
		tt += dt;
	}
	fclose(plik1);

	FILE * plik2;
	plik2 = fopen("crankNicolsonJacobi.txt", "w");
	for(int i = 0; i < m; ++i){
		fprintf(plik2, "%d: ", i);
		for(int j = 0; j < n; ++j){
			fprintf(plik2, "%15.8f", XXX[i][j]);
		}
		fprintf(plik2, "\n");
	}
	fclose(plik2);

	usunWektor(x0);
	usunWektor(bb);
	usunMacierz(A, n);
	usunMacierz(XXX, m);
	return maxBlad;
}

int main(){

	std::cout << "ileH: "<<ileH << ", h: " << h << std::endl;
	std::cout << "ileDt: " << ileDt << ", dt: " << dt << std::endl;

	FILE * maxBledy = fopen("maksymalneBledy.txt", "w");
	for(ileH = 20; ileH < 777; ileH += 20){
		h = 2.0*a/(ileH-1);
		lambda = 1.0;
 		dt = lambda*h*h; // krok czasowy
 		ileDt = static_cast<int>((tmax-t0)/dt + 1); // liczba poziomow czasowych
 	
		double maxBladT = crankNicolsonThomas(ileH, ileDt);
		double maxBladJ = crankNicolsonJacobi(ileH, ileDt);

		fprintf(maxBledy, "%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n", h, std::log10(h), maxBladT, std::log10(maxBladT), maxBladJ, std::log10(maxBladJ));
	}
	fclose(maxBledy);//*/

	//crankNicolsonThomas(ileH, ileDt);
	//crankNicolsonJacobi(ileH, ileDt);
}
