#pragma once

// Metoda iteracyjna Jacobiego
// implementacja w pliku jacobi.cpp

// METODA JACOBIEGO
// rozwiazuje uklad n rownan liniowych postaci Ax = b
// xn - poczatkowe przyblizenie wektora rozwiazan
// tolx - tolerancja błędu wartości pierwiastka
// tolf - tolerancja błędu wartości funkcji
// nmax - maksymalna liczba iteracji
double * jacobi(double ** A, double * b, int n, double * xn, double tolx, double tolf, double nmax);
