
#include "stdafx.h"

#include <math.h>
#include <armadillo>
#include <random>	//uniform random generator
#include <chrono>	//seed
#include <ctime>;

using namespace std;

const double PI = 3.14159265;


double diffclock(clock_t clock1, clock_t clock2)
{
	double diffticks = clock2 - clock1;
	return  (diffticks) / (CLOCKS_PER_SEC);
}


//************************************************************************************
//Calcola il vettore normale del piano dai valori di elevazione (theta) e orientamento (phi)
arma::vec computeNormalVector(double theta, double phi) {

	arma::vec n(3);

	n(0) = cos(theta * PI / 180) * cos(phi * PI / 180);
	n(1) = cos(theta * PI / 180) * sin(phi * PI / 180);
	n(2) = sin(theta * PI / 180);

	return n;
}

//Prende un numero casuale a distribuzione uniforme centrata in x0 e intervallo (x0 - L; x0 + L)
double myUniformRand(double x0, double L) {

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return x0 - L + 2 * L * number;
}

//limita i valori di un vettore entro il range richiesto
arma::vec myBound(arma::vec X, arma::vec x0, arma::vec L) {

	for (int i = 0; i < X.n_elem; ++i) {
		if (X(i) > x0(i) + L(i))
			X(i) = x0(i) + L(i);

		if (X(i) < x0(i) - L(i))
			X(i) = x0(i) - L(i);
	}
	return X;
}

//ROTAZIONI
arma::mat rotx(double angolo_deg) {
	double ang = angolo_deg * PI / 180;

	arma::mat Rx = { { 1,0,0 },
	{ 0, cos(ang), -sin(ang) },
	{ 0, sin(ang), cos(ang) } };
	return Rx;
}

arma::mat roty(double angolo_deg) {
	double ang = angolo_deg * PI / 180;

	arma::mat Rx = { { cos(ang), 0, sin(ang) },
					{ 0, 1, 0 },
					{ -sin(ang), 0, cos(ang) } };
	return Rx;
}

arma::mat rotz(double angolo_deg) {
	double ang = angolo_deg * PI / 180;

	arma::mat Rx = { {cos(ang), -sin(ang), 0 },
					{ sin(ang), cos(ang), 0},
					{ 0, 0, 1} };
	return Rx;
}




arma::mat mySample(arma::mat &A, int M) {

	//estremi
	int a = ceil((double)A.n_rows / (double)M);
	int b = ceil((double)A.n_cols / (double)M);

	arma::vec tmp = arma::zeros(a * b);
	arma::mat B = arma::zeros(a, b);

	int index = 0;
	for (int i = 0; i < A.n_rows; i = i + M) {
		for (int j = 0; j < A.n_cols; j = j + M) {
			tmp(index) = A(i, j);
			++index;
		}
	}

	index = 0;
	for (int i = 0; i < B.n_rows; ++i) {
		for (int j = 0; j < B.n_cols; ++j) {
			B(i, j) = tmp(index);
			++index;
		}
	}

	return B;
}