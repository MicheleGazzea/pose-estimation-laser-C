/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
/	Autore:				Michele Gazzea                              \
/	Data:				02/2018										\
/	Descrizione:		Stima la trasformazione nei 6 DoF			\
/						di una superficie cartesiana simulata		\
/vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/


#include "stdafx.h"

#include <iostream>
#include <armadillo>
#include <random>
#include <chrono>

#include "funcTools.h"
#include "accessori.h"


using namespace std;



int main() {
	

	//inizializzazione modello - posa originale
	spazioLavoro modello;
	modello.Tx = 0.1;	modello.Ty = 0.1;	modello.Dx = 5;	modello.Dy = 5;  	modello.P0 = { 15,0,0 };	 modello.theta = 0;	modello.psi = 0;
	myMeshgrid(modello, 0);

	//generazione funzione (in forma esplicita):
	modello.Zs = 2 * exp(  -0.1* (pow(modello.u, 2) + pow(modello.v, 2) )    ) +
		2 * exp(   -0.4* (pow(modello.u + 3, 2) + pow(modello.v, 2) )    ) +
		3 * exp(   -0.5* (pow(modello.u - 1, 2) + pow(modello.v + 1, 2) )    ) +
		1.5 * exp(   -0.1* (pow(modello.u + 4, 2) + pow(modello.v - 3, 2) )    );

	spazioLavoro scansione = modello; //copia adesso ma poi sciverò sempre su scansione

	// CREAZIONE TEMPLATE
	arma::vec T = { 10,0,90,15,-4,0 }; //trasformazione vera
	T.print("Trasformazione vera");
	myTrasform(T, modello, scansione);
	arma::mat templ = mySlice(scansione);
	myProfileProj(templ, scansione);


	//RICERCA PROFILO
	arma::vec T_stimata = {8,1,88,16,-1,3}; //guess iniziale
	double min = 1000;
	T_stimata.print("Guess iniziale:");		cout << endl;

	opzioniPSO opzioniSciame;
	opzioniSciame.S = 40;
	opzioniSciame.L = {3,3,3,3,3,3};
	opzioniSciame.max_iterazions = 200;		opzioniSciame.lowest_threshold = 0.02;		opzioniSciame.not_changed = 50;

	clock_t start = clock();
	ParticleSwarmOpt(T_stimata, min, templ, modello, opzioniSciame);
	T_stimata.print("Trasformazione stimata1");
	cout << endl << "errore: " << min;
	clock_t end = clock();


	
	cout << "\n  Tempo di esecuzione main(): " << diffclock(start, end) << "s" << endl;
} //fine main


