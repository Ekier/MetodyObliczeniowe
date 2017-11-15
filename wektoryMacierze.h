#pragma once

// Metody dotyczące wektorów i macierzy
// Implementacja w pliku wektoryMacierze.cpp

// alokuje pamiec na wektor o dlugosci n
double * utworzWektor(int n);

// delokuje pamiec z wektora
void usunWektor(double * w);

// wypisuje elementy wektora o dlugosci n poprzedzielone tabulacja
void wypiszWektor(double * w, int n);

// alokuje pamięć na macierz o wymiarach n x n
double ** utworzMacierz(int n);

// dealokuje pamiec z macierzy o liczbie wierszy n
void usunMacierz(double ** A, int n);
