#ifndef HWXCHROM_H_INCLUDED
#define HWXCHROM_H_INCLUDED

/* #define COUNTTYPE unsigned short */
typedef unsigned short COUNTTYPE;

void xChrom(int * rm,
	int * mf,
	int * rk,
	double * robservedVals, // observed stats: LLR, Prob, U, X2
	double * rPvals, // computed P values: LLR, Prob, U, X2
	int * rstatID, // which statistic to use for histogram (1-4)
	int * rhistobins, // number of bins for histogram. (no histogram if 0)
	double * rhistobounds, // Two values indicating the range for histogram
	double * rhistoData, // histogram data. length = histobounds.
	int * rsafeSecs, // abort calculation after this many seconds
	double * tables // the number of tables examined
);
static void twoAlleleSpecialCaseX();

static void heterozygoteX(int r, int c, double probl, double constMales, COUNTTYPE * R);

static void homozygoteX(int r, double probl, double constMales, COUNTTYPE * R);

static int sum();
static int enoughFemale();
static void recursiveEnumeration(int index);

static bool nearlyEqual(double a, double b, double epsilon = 1e-10);


#endif // HWXCHROM_H_INCLUDED
