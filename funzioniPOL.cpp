//                Funzioni per l'utilizzo di superfici POL 

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <string.h>
#include <sys/stat.h>
#include <malloc.h>
#include <float.h>
#include <math.h>

#include "funzioniPOL.h"
#include "CLMATH.H"



//Funzione che calcola la dimensione effettiva dell'header del file POL
int POL_CalcolaDimensioneHeaderPOL(int ntipofigura, int nsemicerchi)
{
	int nSizeH = sizeof(tDATI_POL);
	if (ntipofigura == 0x101)
	{
		int nExtraSize = sizeof(short) + sizeof(double) * 3 + sizeof(double) + sizeof(C_ZONE3D);
		int nExtraSize2 = nExtraSize + sizeof(double) * 3 * 2;
		if (nsemicerchi>0)
			nSizeH += nExtraSize;

		if (nsemicerchi>nExtraSize)
			nSizeH += (nExtraSize2 - nExtraSize);

		if (nsemicerchi>nExtraSize2)
			nSizeH += sizeof(short) + sizeof(double) * 2;
	}
	return nSizeH;
}


// apertura e inizializzazione dei campi della struttura FILE_POL
FILE_POL*
POL_fopen(const char *nomefile, const char *mode)
{
	FILE_POL *ret;                        // uscita della funzione 

	ret = (FILE_POL*)malloc(sizeof(FILE_POL));
	if (ret)
	{
		size_t letti;
		FILE  *pfile;                      // puntatore al file POL 

		ret->nomefile = 0;

		// apertura del file POL  
		pfile = fopen(nomefile, mode);
		if (pfile == NULL)
		{
			ret->errore = 2;  /* errore nell'apertura del file */
			return ret;
		}

		ret->errore = 0;
		ret->sizeHeader = sizeof(tDATI_POL);

		ret->nomefile = (char*)malloc((strlen(nomefile) + 1) * sizeof(char));
		strcpy(ret->nomefile, nomefile);

		if (strcmp(mode, "rb") == 0)
		{
			letti = fread(&(ret->pol), 1, sizeof(tDATI_POL), pfile);
			if (letti != sizeof(tDATI_POL))
				ret->errore = 3; // errore nella lettura dell'header
			else
				ret->sizeHeader = POL_CalcolaDimensioneHeaderPOL(ret->pol.tipofigura, ret->pol.semicerchi);
		}

		fclose(pfile);
	}
	return ret;
}

// libera la memoria utilizzata da *header
int
POL_fclose(FILE_POL *f)
{
	if (f->nomefile)
		free(f->nomefile);
	free(f);
	return 0;
}

// scrittura della superfice un su file tipo POL
void
POL_putsurface(float *surface, FILE_POL *f)
{
	FILE *pfile;                      // puntatore al file POL 

	pfile = fopen(f->nomefile, "wb");
	if (pfile == NULL)
	{
		f->errore = 2;
		return;
	}

	fwrite(&(f->pol), sizeof(tDATI_POL), 1, pfile);
	fwrite(surface, sizeof(float), f->pol.nrigs*f->pol.ncols, pfile);
	fclose(pfile);
}

// restituisce le quote z del POL riferito a *header
float*
POL_getsurface(FILE_POL *f)
{
	size_t letti;
	float *ret;
	FILE  *pfile;                      // puntatore al file POL 

	if (f == NULL)
		return NULL;           // Errore nell' apertura del file sorgente

	pfile = fopen(f->nomefile, "rb");
	if (pfile == NULL)
	{
		f->errore = 2;
		return NULL;
	}

	fseek(pfile, sizeof(tDATI_POL), SEEK_SET);
	ret = (float*)malloc(f->pol.nrigs*f->pol.ncols * sizeof(float));
	if (ret == NULL)
	{
		f->errore = 1;             // Errore: memoria insufficente
		return NULL;
	}
	letti = fread(ret, sizeof(float), f->pol.nrigs*f->pol.ncols, pfile);
	if (letti != f->pol.nrigs*f->pol.ncols)
		f->errore = 4; // errore: header e superfice non combaciano

	fclose(pfile);
	return ret;
}
