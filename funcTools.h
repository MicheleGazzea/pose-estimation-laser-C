/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
/	Autore:				Michele Gazzea                              \
/	Data:				02/2018										\
/	Descrizione:		Funzioni per la stima						\
/						automatica della posa						\
/vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/

#pragma once

#ifndef FUNC_TOOLS_H
#define FUNC_TOOLS_H

#include <armadillo>

struct spazioLavoro {
	double Tx, Ty; //campionamento lungo i due assi (di solito uguale = 0.02)
	double Dx, Dy; //bordi della griglia
	double Ox = 0; //origine x
	double Oy = 0; //origine y
	double righe, colonne;
	arma::mat u;
	arma::mat v;
	arma::mat Zs;
	arma::vec P0; //punto passante per il piano intersecante
	double theta = 0; //elevazione piano
	double psi = 0; //orientazione piano
};

struct opzioniPSO {
	//numero di particelle
	int S = 0; 

	//parametri di comportamento
	double omega = 1.0;
	double PHIp = 0.4; 
	double PHIg = 0.2; 
	
	//spazio di ricerca
	arma::vec L = {2,2,2,2,2,2}; 

	//criteri di stop: ferma l'algoritmo se:
	double lowest_threshold; //raggiungi il minimo richiesto
	int max_iterazions;	//raggiungi le iterazioni massime ammesse
	int not_changed; //il minimo non cambia dopo # iterazioni
};

struct Particella {
	arma::mat stato = arma::zeros( 6, 3 ); //[position | velocity | best_position]
};


void myMeshgrid(spazioLavoro&, int mode);
//inizializza la griglia di lavoro u,v,Zs
//mode 0: simulazione con superifcie cartesiana
//mode 1: simulazione con POL



// ******************************************************************
// ParticleSwarmOpt: Usa il metodo PSO per cercare la trasformazione che minimizza l'errore
// tra template e profilo
void ParticleSwarmOpt(arma::vec &T, double &min, arma::mat templ, spazioLavoro  &modello, opzioniPSO &opzioniSciame);



// ******************************************************************
// Funzione obiettivo: da minimizzare
double errorFun(arma::vec T, arma::mat templ, spazioLavoro &modello, spazioLavoro &scansione);




// FUNZIONI INTERNE AD errorFun
//****************************

//Roto-traslo la superficie Zs di Trasf
void myTrasform(arma::vec T, spazioLavoro &modello, spazioLavoro &scansione);

//Taglio la superficie con il piano (fisso)
arma::mat mySlice(spazioLavoro &scansione);

//Proietta sul piano intersecante sezione e template
void myProfileProj(arma::mat &profilo, spazioLavoro &scansione);

//Confronto tra template e profilo estratto
double myProfileComparison(arma::mat templ, arma::mat &profilo);

#endif