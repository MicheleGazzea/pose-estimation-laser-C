#pragma once

#ifndef ACCESSORI_H
#define ACCESSORI_H

#include <armadillo>
#include <ctime> //clock()

double diffclock(clock_t clock1, clock_t clock2);


//Calcola il vettore normale del piano dai valori di elevazione (theta) e orientamento (phi)
arma::vec computeNormalVector(double theta, double phi);

//Prende un numero casuale a distribuzione uniforme centrata in x0 e intervallo (x0 - L; x0 + L)
double myUniformRand(double x0, double L);

//limita i valori di un vettore entro il range richiesto
arma::vec myBound(arma::vec X, arma::vec x0, arma::vec L);


//ROTAZIONI
arma::mat rotx(double angolo_deg);
arma::mat roty(double angolo_deg);
arma::mat rotz(double angolo_deg);


//Campiona una matrice con valore M
arma::mat mySample(arma::mat &A, int M);






#endif


