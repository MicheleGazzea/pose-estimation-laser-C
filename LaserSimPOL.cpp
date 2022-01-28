/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
/	Autore:				Michele Gazzea                              \
/	Data:				02/2018										\
/	Descrizione:		Stima la trasformazione nei 6 DoF			\
/						di un modello POL							\
/vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/

#include "stdafx.h"

#include <iostream>
#include <armadillo>
#include <random>
#include <chrono>
#include <fstream>

#include "funcTools.h" 
#include "accessori.h"

#include "funzioniPOL.h"


using namespace std;

int main() {

	spazioLavoro modello; //creo lo spazio per modello
	
	//inizializzo il modello
	FILE_POL* file = POL_fopen("Model2_originale.pol", "rb");
	modello.Tx = file->pol.stepscanx;		modello.Ty = file->pol.stepscany;
	modello.Ox = file->pol.orgx;			modello.Oy = file->pol.orgy;
	modello.righe = file->pol.nrigs;		modello.colonne = file->pol.ncols;
	modello.P0 = { -5,0,0 };

	myMeshgrid(modello, 1);

	float* tmp = POL_getsurface(file);
	int index = 0;
	modello.Zs.set_size(modello.colonne, modello.righe);
	for (int j = 0; j < modello.Zs.n_cols; ++j) {
		for (int i = 0; i < modello.Zs.n_rows; ++i) {
			modello.Zs(i, j) = *(tmp+index);
			++index;
		}
	}

	spazioLavoro scansione;
	//inizializzo la scansione
	file = POL_fopen("Model2_ruotato.pol", "rb");
	scansione.Tx = file->pol.stepscanx;		scansione.Ty = file->pol.stepscany;
	scansione.Ox = file->pol.orgx;			scansione.Oy = file->pol.orgy;
	scansione.righe = file->pol.nrigs;		scansione.colonne = file->pol.ncols;
	myMeshgrid(scansione, 1);
	tmp = POL_getsurface(file);
	index = 0;
	scansione.Zs.set_size(scansione.colonne, scansione.righe);
	for (int j = 0; j < scansione.Zs.n_cols; ++j) {
		for (int i = 0; i < scansione.Zs.n_rows; ++i) {
			scansione.Zs(i, j) = *(tmp + index);
			++index;
		}
	}


	//CAMPIONAMENTO
	int M = 6;		modello.Tx = 0.02*M;	scansione.Tx = 0.02*M;

	modello.u = mySample(modello.u, M);		modello.v = mySample(modello.v, M);
	modello.Zs = mySample(modello.Zs, M);

	scansione.u = mySample(scansione.u, M);		scansione.v = mySample(scansione.v, M);
	scansione.Zs = mySample(scansione.Zs, M);


	//GENERO UN TEMPLATE
	scansione.P0 = { -5,0,0 };
	arma::mat templ = mySlice(scansione);

	
	//RICERCA PROFILO
	arma::vec T = { 0,0,5,0,0,5 };
	T.print("Guess iniziale:");		cout << endl;
	double minimo = 10000;

	opzioniPSO opzioniSciame;
	opzioniSciame.S = 30;
	opzioniSciame.L = { 0.5,0.5,3,5,5,3 };
	opzioniSciame.max_iterazions = 200;		opzioniSciame.lowest_threshold = 0.1;		opzioniSciame.not_changed = 50;

	clock_t start = clock();
	ParticleSwarmOpt(T, minimo, templ, modello, opzioniSciame);
	clock_t end = clock();

	T.print("Trasformazione stimata:");
	cout << "\n" << "minimo: " << minimo;

	

	cout << "\n  Tempo di esecuzione main(): " << diffclock(start, end) << "s" << endl;


} //fine main


