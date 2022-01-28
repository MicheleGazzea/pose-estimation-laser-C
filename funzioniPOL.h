#pragma warning(disable:4996)


//                Funzioni per l'utilizzo di superfici POL 

//#include "tchar.h" // for TCHAR
//#include "Winnt.h" //for LPCTSTR

enum costantiTipoFigura {
	SCAN_NULL, SCAN_GRECA, SCAN_SPIRALE, \
	SCAN_SPIRALE_INVERSA, SCAN_GRECA_INVERSA, \
	SCAN_RADIALE
};

enum costantiTipoModello { NULL_MODE, MASCHIO, FEMMINA };

typedef struct _DATI_POL
{
	short tipofigura;       // tipo di superfice
	short tipomodello;      // modello superfice
	long  nrigs;            // numero di righe
	long  ncols;             // numero di colonne
	long  semicerchi;       // ???
	float stepscanx;        // passo superfice in x
	float diametro;         // ???
	float alt;              // ???
	float larg_pezzo;       // larghezza superfice
	float lung_pezzo;       // lunghezza pezzo
	float quota_max;        // quota massima raggiunta
	float quota_min;        // quota minima raggiunta
	float initquotasuperf;  // quota di inizializzazione
	float stepscany;        // passo superfice in y
	float orgx;             // origine X  ----> 0
	float orgy;             // origine Y	 ----> 0

} tDATI_POL;

typedef struct _FILE_POL
{
	char *nomefile;
	int errore;
	int sizeHeader;
	tDATI_POL pol;

} FILE_POL;


// apertura e inizializzazione dei campi della struttura FILE_POL
FILE_POL* POL_fopen(const char *nomefile, const char *mode);

// libera la memoria utilizzata da *header
int POL_fclose(FILE_POL *f);

// restituisce le quote z del POL riferito a *header
float* POL_getsurface(FILE_POL *f);

// scrittura della superfice un su file tipo POL
void POL_putsurface(float *surface, FILE_POL *f);

// calcola la dimensione effettiva dell'header di un file POL
int POL_CalcolaDimensioneHeaderPOL(int ntipofigura, int nsemicerchi);