/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
/	Autore:				Michele Gazzea                              \
/	Data:				02/2018										\
/	Descrizione:		Funzioni per la stima						\
/						automatica della posa						\
/vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/

#include "stdafx.h"
#include "funcTools.h"

#include <math.h>	//pow, sin cos
#include <armadillo>
#include <random>	//uniform random generator
#include <chrono>	//seed
#include <algorithm> //max,min

#include "accessori.h"

double INF = 1000000;
using namespace std;


void myMeshgrid(spazioLavoro &modello, int mode) {
	
	//con le superfici simulate
	if (mode == 0) {

		modello.u.set_size(2 * modello.Dy / modello.Ty + 1, 2 * modello.Dx / modello.Tx + 1);
		modello.v.set_size(2 * modello.Dy / modello.Ty + 1, 2 * modello.Dx / modello.Tx + 1);

		arma::rowvec Xaxis = arma::regspace(-modello.Dx, modello.Tx, modello.Dx).t();
		arma::vec Yaxis = arma::regspace(-modello.Dy, modello.Ty, modello.Dy);

		for (int i = 0; i <= (2 * modello.Dy) / modello.Ty; ++i) {
			modello.u.row(i) = Xaxis;
		}
		for (int i = 0; i <= (2 * modello.Dx) / modello.Tx; ++i) {
			modello.v.col(i) = Yaxis;
		}
	}

	//con i file .POL
	if (mode == 1) {
		modello.u.set_size(modello.colonne, modello.righe);
		modello.v.set_size(modello.colonne, modello.righe);

		arma::rowvec Xaxis = arma::regspace(modello.Ox, modello.Tx, modello.Ox + (modello.Tx * (modello.righe - 1))).t();
		arma::vec Yaxis = arma::regspace(modello.Oy, modello.Ty, modello.Oy + (modello.Ty * (modello.colonne - 1)));

		for (int i = 0; i < modello.colonne; ++i) 
			modello.u.row(i) = Xaxis;

		for (int i = 0; i < modello.righe; ++i)
			modello.v.col(i) = Yaxis;
	}
}




// Usa il metodo PSO per cercare la trasformazione che minimizza l'errore tra template e profilo

void ParticleSwarmOpt(arma::vec &T,double &minimo, arma::mat templ, spazioLavoro  &modello, opzioniPSO &opzioniSciame) {
	
	//creo lo spazio di lavoro per mettere i "tentativi" delle varie trasformazioni
	spazioLavoro scansione = modello;

	// Crea uno sciame (array di N particelle):
	int Ns = opzioniSciame.S;
	double omega = opzioniSciame.omega;
	double PHIp = opzioniSciame.PHIp;
	double PHIg = opzioniSciame.PHIg;
	arma::vec L = opzioniSciame.L;
	Particella* swarm = new Particella[Ns];
	arma::vec T0 = T; //guess iniziale

	//INIZIALIZZAZIONE
	for (int i = 0; i < Ns; ++i) {
		//posizioni:
		swarm[i].stato(0, 0) = myUniformRand(T[0], L[0]);
		swarm[i].stato(1, 0) = myUniformRand(T[1], L[1]);
		swarm[i].stato(2, 0) = myUniformRand(T[2], L[2]);
		swarm[i].stato(3, 0) = myUniformRand(T[3], L[3]);
		swarm[i].stato(4, 0) = myUniformRand(T[4], L[4]);
		swarm[i].stato(5, 0) = myUniformRand(T[5], L[5]);

		//velocità: initializzata di default a 0
		
		//best_position:
		swarm[i].stato(0, 2) = swarm[i].stato(0, 0);
		swarm[i].stato(1, 2) = swarm[i].stato(1, 0);
		swarm[i].stato(2, 2) = swarm[i].stato(2, 0);
		swarm[i].stato(3, 2) = swarm[i].stato(3, 0);
		swarm[i].stato(4, 2) = swarm[i].stato(4, 0);
		swarm[i].stato(5, 2) = swarm[i].stato(5, 0);
	}

	// INIZIO CICLO

	int iterazione = 1;
	bool flag = false;
	int counter = 0;
	arma::vec goal = T; //output

	// random settings
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	while (iterazione<opzioniSciame.max_iterazions && minimo>opzioniSciame.lowest_threshold && !flag) {

		for (int i = 0; i < Ns; ++i) { //per ogni particella
			for (int j = 0; j < 6; ++j) { //oer ogni dimensione
				
				double Rp = distribution(generator);
				double Rg = distribution(generator);

				//aggiorna velocità:
				swarm[i].stato(j,1) =  omega * swarm[i].stato(j, 1) + 
										PHIp * Rp * (swarm[i].stato(j, 2) - swarm[i].stato(j, 0)) + 
										PHIg * Rg *(goal(j) - swarm[i].stato(j, 0));
			}
			// aggiorna posizione: p(t+1)=p(t)+v(t) ed eventualmente limita entro lo spazio di ricerca
			swarm[i].stato.col(0) = myBound( swarm[i].stato.col(0) + swarm[i].stato.col(1), T0, L ) ;

			//aggiorna eventualmente la best_position
			if (errorFun(swarm[i].stato.col(0), templ, modello, scansione) < errorFun(swarm[i].stato.col(2), templ, modello, scansione)) {
				swarm[i].stato.col(2) = swarm[i].stato.col(0);
				if (errorFun(swarm[i].stato.col(2), templ, modello, scansione) < errorFun(goal, templ, modello, scansione)) {
					goal = swarm[i].stato.col(2);
					counter = 0;
				}
			}
		}

		//se dopo #iterazioni il minimo non è cambiato -> esci dal ciclo
		if (counter == opzioniSciame.not_changed) {
			flag = true;
		}

		++iterazione;
		++counter;
		minimo = errorFun(goal, templ, modello, scansione);
		T = goal;

		if (flag == true)
			cout << endl << "*** uscita per flag ***" << endl;
		else if (iterazione == opzioniSciame.max_iterazions)
			cout << endl << "*** uscita per numero massimo di iterazioni ***" << endl;
		else if (minimo < opzioniSciame.lowest_threshold)
			cout << endl << "*** uscita per minimo raggiunto ***"  << endl;

	}
	
	delete [] swarm; //When done, free memory pointed to by swarm.
	swarm = NULL; // Clear swarm to prevent using invalid memory reference.
	
}


// FUNZIONE PER CALCOLARE L'ERRORE

double errorFun(arma::vec T, arma::mat templ, spazioLavoro &modello, spazioLavoro &scansione) {

	//Roto-traslo la superficie Zs di Trasf, interseco col piano e confronto con il template
	myTrasform(T, modello, scansione);

	//Taglio la superficie con il piano (fisso)
	arma::mat profilo = mySlice(scansione);

	//Proietta sul piano intersecante sezione e template
	myProfileProj(profilo, scansione);

	//Confronto tra template e profilo estratto
	double errore = myProfileComparison(templ, profilo);


	return errore;
}



// FUNZIONI INTERNE AD errorFun


  //***********  <trasform>  ***********
void myTrasform(arma::vec T, spazioLavoro &modello, spazioLavoro &scansione) {
	//trasformo la superficie per estrarre il template

	//Prende la superficie di partenza e la roto-trasla del vettore T
	arma::mat R = rotz(T(2)) * roty(T(1)) * rotx(T(0));

	arma::vec punto = { 0,0,0 }; //punto generico temporaneo

	for (int i = 0; i < modello.u.n_rows; ++i) {
		for (int j = 0; j < modello.u.n_cols; ++j) {

			punto = { modello.u(i,j), modello.v(i,j), modello.Zs(i,j) };
			scansione.u(i, j) = as_scalar(R.row(0) * punto) + T(3);
			scansione.v(i, j) = as_scalar(R.row(1) * punto) + T(4);
			scansione.Zs(i, j) = as_scalar(R.row(2) * punto) + T(5);
		}
	}
}

//***********  <sliceSimple>  ***********
arma::mat mySlice(spazioLavoro &scansione) {

	//calcolo per ogni punto (=riga della sezione) la distanza dal piano intersecante: 
	//considero solo i punti la cui distanza è inferiore a Th (=threshold)

	double Th = 0.6 * scansione.Tx;

	arma::vec n = computeNormalVector(scansione.theta, scansione.psi);
	arma::mat distanza = abs(n(0)*(scansione.u - scansione.P0(0)) + n(1)*(scansione.v - scansione.P0(1)) + n(2)*(scansione.Zs - scansione.P0(2)));

	arma::mat profilo = arma::zeros(4 * max(scansione.u.n_rows, scansione.u.n_cols), 3);

	//scansiona la matrice "distanza" e prende solo i valori in [u,v,Zs] la cui distanza è minore di Th  
	int index = 0;
	for (int i = 0; i < distanza.n_rows; ++i) {
		for (int j = 0; j < distanza.n_cols; ++j) {

			if (distanza(i, j) < Th) {
				profilo(index, 0) = scansione.u(i, j);
				profilo(index, 1) = scansione.v(i, j);
				profilo(index, 2) = scansione.Zs(i, j);
				++index;
			}
		}
	}
	profilo.resize(index, 3);

	return profilo;
}


//***********  <profiloProjection>  ***********
void myProfileProj(arma::mat &profilo, spazioLavoro &scansione) {
	//proietto i punti appena scoperti sul piano : in questo modo sono tutti collineati in un certo s.d.r
	arma::vec n = computeNormalVector(scansione.theta, scansione.psi);

	for (int i = 0; i < profilo.n_rows; ++i) {
		arma::rowvec tmp = { scansione.P0(0), scansione.P0(1), scansione.P0(2) };
		profilo.row(i) = profilo.row(i) - tmp - (as_scalar((profilo.row(i) - tmp) * n) * n).t();
		
		//% ruoto i punti al contrario rispetto a theta e psi:
		// in questo modo sono tutti paralleli al piano XZ e posso e confontarli guardando solo le cordinate x e z
		arma::mat R = rotz(-scansione.psi)*rotx(-scansione.theta);
		profilo.row(i) = ( profilo.row(i) * R.t() ) + tmp;
	}
}


//***********  <profileComparison>  ***********
//Confronto tra template e profilo estratto
double myProfileComparison(arma::mat templ, arma::mat &profilo) {

	double errore;

	if (templ.n_rows < 20 || profilo.n_rows < 20) {
		errore = INF; //profilo acquisito troppo corto per una stima affidabile, evita tutto il calcolo di errorFun
	}
	else {

		//campionamento del template: 
		//---> (si può spostare sul main o su PSO per non dover campionare ad ogni chiamata di errorFun il template che è sempre uguale)
		int M = 5;
		int index = 0;
		for (int i = 0; i < templ.n_rows; i = i + M) {
			templ.row(index) = templ.row(i);
			++index;
		}
		templ.resize(index, 3);

		//campionamento del profilo
		index = 0;
		for (int i = 0; i < profilo.n_rows; i = i + M) {
			profilo.row(index) = profilo.row(i);
			++index;
		}
		profilo.resize(index, 3);

		double Sa = min(profilo.col(1));   double Sb = max(profilo.col(1));
		double Ta = min(templ.col(1));   double Tb = max(templ.col(1));

		if ((Sb - Sa) < (Tb - Ta) / 2)
			errore = INF;  //evita il confronto se il profilo è troppo corto rispetto al template 

		else {

			arma::vec Xaxis = arma::regspace(max(Sa, Ta), 0.2, min(Sb, Tb));

			int limite = Xaxis.n_rows;

			arma::vec templateI(Xaxis.n_rows); //non inizializzati
			arma::vec profiloI(Xaxis.n_rows);

			if (Xaxis.n_elem > 20) {
				//interpolo il template
				for (int i = 0; i < templ.n_rows - 1; ++i) {
					double P1x = templ.row(i)(1);			double P1y = templ.row(i)(2);
					double P2x = templ.row(i + 1)(1);		double P2y = templ.row(i + 1)(2);
					//cerco sul Xaxis le ascisse tra P1x e P2x
					index = 0;
					while (Xaxis(index) < P1x) {
						++index;
						if (index == Xaxis.n_elem - 1)
							break;
					}


					while (Xaxis(index) >= P1x && Xaxis(index) <= P2x) {	//range giusto [P1x,P2x]
						templateI(index) = P1y + ((P2y - P1y) / (P2x - P1x)) * (Xaxis(index) - P1x);
						++index;
						if (index == limite)
							break;
					}
				}

				//interpolo il profilo
				for (int i = 0; i < profilo.n_rows - 1; ++i) {
					double P1x = profilo.row(i)(1);			double P1y = profilo.row(i)(2);
					double P2x = profilo.row(i + 1)(1);		double P2y = profilo.row(i + 1)(2);
					//cerco sul Xaxis le ascisse tra P1x e P2x
					index = 0;
					while (Xaxis(index) < P1x) {
						++index;
						if (index == Xaxis.n_elem - 1)
							break;
					}

					while (Xaxis(index) >= P1x && Xaxis(index) <= P2x) {	//range giusto [P1x,P2x]
						profiloI(index) = P1y + ((P2y - P1y) / (P2x - P1x)) * (Xaxis(index) - P1x);
						++index;
						if (index == limite)
							break;
					}
				}

				//errore = arma::max(pow((templateI - profiloI), 2));
				errore = arma::accu(pow((templateI - profiloI), 2)) / templateI.n_elem; //media
			}
			else {
				errore = INF;
			}

		}//fine else

	}//fine else

	return errore;
}
