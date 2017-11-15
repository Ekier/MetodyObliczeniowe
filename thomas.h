#pragma once

// Metody do rozwiazywania ukladu rownan z uzyciem algorytmu Thomasa
// implementacja w pliku thomas.cpp


// dekompozycja macierzy, d - wektor elementow diagonalnych macierzy, u - wektor elementow nad diagonala, l - wektor elementow pod diagonala
// n - rozmiar macierzy
double * dekompozycjaLUmacierz(double * l, double * d, double * u, int n);

// dekompozycja wektora, korzysta z wektora eta zwroconego przez funkcje dekompozycjaLUmacierz
// n - dlugosc wektora
double * dekompozycjaLUwektor(double * b, double * eta, double * l, int n);

// rozwiazuje uklad rownan Ax = b, korzysta z wektorow eta i r zwroconych odpowiednio przez funkcje dekompozycjaLUmacierz i dekompozycjaLUwektor
// u - wektor elementow nad przekatna macierzy A
// n - rozmiar macierzy
double * rozwiazUklad(double * eta, double * r, double * u, int n);
