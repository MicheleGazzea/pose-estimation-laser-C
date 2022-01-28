/****************************************************************************/
/*                                                                          */
/* Progetto:    Incicam CiElle.                                             */
/* Modulo:      Pih.h.                                                      */
/* Descrizione: Costanti e strutture dati retative al formato PIH.          */
/* Autore:      Michele Dekleva                                             */
/* Data:        8 Dicembre 1997                                             */
/* Versione:    3.0.                                                        */
/*                                                                          */
/* Revisione:   Daniele De Giorgio                                          */
/* Data:        24 Luglio 2000                                              */
/* Versione:    4.0.                                                        */
/* Commento:    Aggiunta gestione delle entità in formato double            */
/*                                                                          */
/****************************************************************************/
#ifndef PIH_H
#define PIH_H


/* allineamento al byte per efficienza nella memoriazzazione su disco */
#pragma pack(1)


/* versione corrente del file */
#define PIH_VERSION     4000    /* versione 4.0                             */


/* tipo di file */
#define OLD_TEC_FILE    10      /* file PIH 5 words chiavi                  */
#define FONT_FILE       11      /* file FON 25 words chiavi                 */
#define TECNICAL_FILE   12      /* file PIH 25 words chiavi                 */
#define FON_PREC        13      /* file FON 25 words chiavi, ARC_C+ARC_E    */
#define PIH_PREC        14      /* file PIH 25 words chiavi, ARC_C+ARC_E    */
#define PIH_PREC2       15      /* file PIH formato double,  ARC_C          */
#define PIH_PREC2_ELL   16      /* file PIH formato double,  ARC_C+ARC_E    */


/* numero di caratteri per FON (ASCII) */
#define FON_CHARS       256


/* unita` di misura utilizzabili */
#define UM_NDEF         0      /* unita` di misura non definita             */
#define UM_MM           1      /* millimetri                                */
#define UM_INCH         2      /* pollici                                   */

#define UM_DEG          0      /* gradi                                     */
#define UM_RAD          1      /* radianti                                  */
#define UM_DECDEG       2      /* gradi centesimali                         */


/* definizioni di uso interno */
typedef short              K_KEY;
typedef unsigned short     K_ETYPE;
typedef unsigned short     K_FTYPE;
typedef unsigned short     K_COUNT;
typedef unsigned short     K_U_MIS;
typedef unsigned short     K_UM_ANG;
typedef unsigned long      K_FLAGS;
typedef unsigned long      K_VERS;
typedef float              K_COORD;
typedef float              K_MARK;
typedef double             K_DIST;
typedef double             K2_COORD;


/* entita` gestite correntemente */
#define E_TMFILE              65535 /* di termine file (-1)              */
#define E_TMLIST              0     /* di termine lista                  */
#define E_START               1     /* di inizio disegno                 */
#define E_PENDOWN             2     /* comando PENDOWN                   */
#define E_PENUP               3     /* comando PENUP                     */
#define E_LINE                4     /* linea                             */
#define E_CURVE               5     /* curva (con mark.....)             */
#define E_POLYLINE            6     /* spezzata                          */
#define E_ARCH                7     /* arco di ellisse (obsoleto)        */
#define E_FULL                8     /* di regione da non svuotare        */
#define E_EMPTY               9     /* di regione da svuotare            */
#define E_ENDBORD             10    /* di fine regione                   */
#define E_FONTPAR             11    /* parametri generali di font        */
#define E_CHARPAR             12    /* parametri del singolo carattere   */
#define E_INDEX               13    /* offset inizio disegni nel file    */
#define E_ARCH_C              14    /* arco di cerchio                   */
#define E_ARCH_E              15    /* arco di ellisse                   */
#define E_ORIGINATOR          16    /* definizione del formato originale */
#define E_MACRO_ENTITY_START  17    /* inizio di una entita` multipla    */
#define E_MACRO_ENTITY_END    18    /* fine di una entita` multipla      */
#define E_CLOSED_PATH_START   19    /* inizio di un percorso chiuso      */
#define E_CLOSED_PATH_END     20    /* fine di un percorso chiuso        */
#define E_REGION_START        21    /* inizio di una zona                */
#define E_REGION_END          22    /* fine di una zona                  */
#define E_BEZIER              23    /* curva di bezier (n-punti)         */
#define E_LINE2               24    /* linea double                      */
#define E_POLYLINE2           25    /* polilinea double                  */
#define E_ARCH_C2             26    /* arco di cerchio double            */
#define E_ARCH_E2             27    /* arco di ellisse double            */


/* possibili formati originali di un file PIH */
#define ORIGIN_UNKNOW         0     /* derivazione sconosciuta          */
#define ORIGIN_CIELLE         1     /* derivazione da programmi CiElle  */
#define ORIGIN_DXF            2     /* derivazione da DXF (AutoCad)     */
#define ORIGIN_EPS            3     /* derivazione da EPS               */
#define ORIGIN_SCV            4     /* derivazione da SCV               */


/* definizione di un punto geometrico */
typedef struct _K_POINT  *PK_POINT;
typedef struct _K_POINT {

   K_COORD      x;
   K_COORD      y;

} K_POINT;


/* definizione di un punto geometrico 3D */
typedef struct _K_POINT3D  *PK_POINT3D;
typedef struct _K_POINT3D {

   K_COORD      x;
   K_COORD      y;
   K_COORD      z;

} K_POINT3D;


/* definizione rettangolo di contenimento */
typedef struct _K_ZONE     *PK_ZONE;
typedef struct _K_ZONE {

   K_COORD      minx;
   K_COORD      miny;
   K_COORD      maxx;
   K_COORD      maxy;

} K_ZONE;


/* struttura del record di testa files FON/PIH (5 words chiavi) */
typedef struct _OLDPIHHEAD    *POLDPIHHEAD;
typedef struct _OLDPIHHEAD {

   K_FTYPE  descriptor;       /* descrittore tipo file            */
   K_ZONE   Sheet;            /* dimensioni del foglio di disegno */
   K_ZONE   Draw;             /* dimensioni reali del disegno     */
   K_U_MIS  unit;             /* unita' di misura coordinate      */
   K_UM_ANG angle;            /* unita' di misura angoli          */
   K_COUNT  num_sheets;       /* numero di disegni contenuti      */
   K_FLAGS  flags;            /* flags vari associati al file     */
   K_VERS   version;          /* versione del file                */
   K_KEY    pal[5];           /* parametri per pal protezione     */

} OLDPIHHEAD;


/* struttura del record di testa files FON/PIH (25 words chiavi) */
typedef struct _PIHHEAD       *PPIHHEAD;
typedef struct _PIHHEAD {

   K_FTYPE  descriptor;       /* descrittore tipo file            */
   K_ZONE   Sheet;            /* dimensioni del foglio di disegno */
   K_ZONE   Draw;             /* dimensioni reali del disegno     */
   K_U_MIS  unit;             /* unita' di misura coordinate      */
   K_UM_ANG angle;            /* unita' di misura angoli          */
   K_COUNT  num_sheets;       /* numero di disegni contenuti      */
   K_FLAGS  flags;            /* flags vari associati al file     */
   K_VERS   version;          /* versione del file                */
   K_KEY    pal[25];          /* parametri per pal protezione     */

} PIHHEAD;


/* definizione record entita' files FON */
typedef struct _K_ENTITY      *PK_ENTITY;
typedef struct _K_ENTITY {

   K_ETYPE     type;           /* tipo entita' (seg, curva, etc.)  */
   K_COUNT     num_bytes;      /* bytes (oltre a questa testa)     */

} K_ENTITY;


/* definizione entita' di inizio carattere */
typedef struct _K_START       *PK_START;
typedef struct _K_START {

   K_ETYPE     type;             /* tipo entita' (seg, curva, etc.)  */
   K_COUNT     num_bytes;        /* bytes (oltre a questa testa)     */
   K_ZONE      Sheet;            /* dimensioni del foglio di disegno */
   K_ZONE      Draw;             /* dimensioni reali del disegno     */

} K_START;


/* struttura definizione entita` arco di ellisse (obsoleto) */
typedef struct _ARC_OBS       *PARC_OBS;
typedef struct _ARC_OBS {

   K_COORD      cx;            /* ascissa del centro                        */
   K_COORD      cy;            /* ordinata del centro                       */
   K_COORD      ra;            /* semiasse principale                       */
   K_COORD      rb;            /* semiasse secondario                       */
   K_COORD      inclinazione;  /* inclinazione semiasse A rispetto asse X   */
   K_COORD      start;         /* angolo iniziale dell'arco                 */
   K_COORD      stop;          /* angolo finale                             */

} ARC_OBS;


/* struttura definizione entita` arco di cerchio */
typedef struct _ARC_C       *PARC_C;
typedef struct _ARC_C {

   K_COORD     cx;            /* ascissa del centro                        */
   K_COORD     cy;            /* ordinata del centro                       */
   K_COORD     r;             /* raggio dell cerchio associato             */
   K_COORD     start;         /* angolo iniziale dell'arco                 */
   K_COORD     ampiezza;      /* ampiezza dell'arco (>0 => antiorario)     */

} ARC_C;


/* struttura definizione entita` arco di cerchio in formato double */
typedef struct _ARC_C2       *PARC_C2;
typedef struct _ARC_C2 {

   K2_COORD     cx;            /* ascissa del centro                        */
   K2_COORD     cy;            /* ordinata del centro                       */
   K2_COORD     r;             /* raggio dell cerchio associato             */
   K2_COORD     start;         /* angolo iniziale dell'arco                 */
   K2_COORD     ampiezza;      /* ampiezza dell'arco (>0 => antiorario)     */

} ARC_C2;


/* struttura definizione entita` arco di ellisse */
typedef struct _ARC_E       *PARC_E;
typedef struct _ARC_E {

   K_COORD      cx;            /* ascissa del centro                        */
   K_COORD      cy;            /* ordinata del centro                       */
   K_COORD      ra;            /* semiasse principale                       */
   K_COORD      rb;            /* semiasse secondario                       */
   K_COORD      inclinazione;  /* inclinazione semiasse A rispetto asse X   */
   K_COORD      start;         /* angolo iniziale dell'arco                 */
   K_COORD      ampiezza;      /* ampiezza dell'arco (>0 => antiorario)     */

} ARC_E;


/* struttura definizione entita` arco di ellisse */
typedef struct _ARC_E2       *PARC_E2;
typedef struct _ARC_E2 {

   K2_COORD      cx;            /* ascissa del centro                        */
   K2_COORD      cy;            /* ordinata del centro                       */
   K2_COORD      ra;            /* semiasse principale                       */
   K2_COORD      rb;            /* semiasse secondario                       */
   K2_COORD      inclinazione;  /* inclinazione semiasse A rispetto asse X   */
   K2_COORD      start;         /* angolo iniziale dell'arco                 */
   K2_COORD      ampiezza;      /* ampiezza dell'arco (>0 => antiorario)     */

} ARC_E2;


/* rilascio vincoli sull'allineamento delle strutture */
#pragma pack()
#endif
