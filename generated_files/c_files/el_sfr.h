#ifndef EL_SFR_H
#define EL_SFR_H

/* T [internal_units] * T_MYR = T [Myr] */
#define T_MYR (All.UnitTime_in_s / All.HubbleParam / SEC_PER_MEGAYEAR)
/* RHO [internal_units] * RHO_COSMO = RHO [cm^(-3)] */
#define RHO_COSMO (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv / PROTONMASS)
/* M [internal_units] * M_COSMO = M [Mₒ] */
#define M_COSMO (All.UnitMass_in_g / SOLAR_MASS)

/* Interpolation tables */
#define ETA_NROWS 107  // Number of rows in the η tables
#define ETA_NCOLS 7  // Number of columns in the η tables
#define R_ZSN_NROWS 7  // Number of rows in the R and Zsn table
#define R_ZSN_NCOLS 2  // Number of columns in the R and Zsn table

/* Paths */
static char *ETA_D_TABLE_PATH = "../code/src/el_sfr/tables/eta_d.txt";
static char *ETA_I_TABLE_PATH = "../code/src/el_sfr/tables/eta_i.txt";
static char *R_TABLE_PATH     = "../code/src/el_sfr/tables/R.txt";
static char *ZSN_TABLE_PATH   = "../code/src/el_sfr/tables/Zsn.txt";

/* ODE constants */

/* ϵff  = 1.0000 (stellar formation efficiency) */
/* Zsun = 0.0127 (solar metallicity) */
/* Cρ   = 100.0000 (clumping factor) */

#define N_EQU 6 /* Number of equations */
#define ODE_CS 5.1470081e+01  /* [Myr * cm^(-3/2)] */
#define ODE_CR 1.2187726e-01   /* [Myr * cm^(-3)] */
#define ODE_CC 5.7491245e-02  /* [Myr * cm^(-3)] */
#define ODE_CD 1.7736986e+03    /* [Myr * mp * cm^(-3)] */
#define TAU_DD 2.2953846e+03    /* [Myr] */
#define ZEFF 1.2700e-05      /* 1e-3 Zₒ */
#define CXD 2.8360e-01       /* [dimensionless] */

typedef struct DataTable
{
	double *data;  // Values of the table
	int n_rows;    // Number of rows in the table
	int n_cols;    // Number of columns in the table
} data_table;

#ifdef RHO_PDF

/*
 * Density PDF according to Burkhart (2018)
 * https://doi.org/10.3847/1538-4357/aad002
 *
 * We used the following parameters (all dimensionless):
 *
 * divisions = 20
 * range of ln(rho/rho_0) = (-6.0, 6.0)
 * α (power law slope) = 2.0
 * b (turbulent forcing parameter) = 0.5
 * Ms (mach number) = 10.0
 */

#define DIVISIONS 20

/* Integrated PDF of the interstellar gas density */
static const double PDF[] = {
	0.0107028600,
	0.0212984957,
	0.0379876559,
	0.0607272041,
	0.0870108777,
	0.1117416340,
	0.1286201183,
	0.1326954080,
	0.1227032705,
	0.1016972610,
	0.0755464976,
	0.0503002261,
	0.0300174715,
	0.0160555202,
	0.0076969348,
	0.0033071276,
	0.0012735615,
	0.0004395623,
	0.0001370353,
	0.0000412781,

};

/* Density factor: ρ = ρ₀ * F_RHO */
static const double F_RHO[] = {
	0.0033459655,
	0.0060967466,
	0.0111089965,
	0.0202419114,
	0.0368831674,
	0.0672055127,
	0.1224564283,
	0.2231301601,
	0.4065696597,
	0.7408182207,
	1.3498588076,
	2.4596031112,
	4.4816890703,
	8.1661699126,
	14.8797317249,
	27.1126389207,
	49.4024491055,
	90.0171313005,
	164.0219072999,
	298.8674009671,

};

#else /* #ifdef RHO_PDF */

#define DIVISIONS 1

static const double PDF[] = {
    1.0,
};

static const double F_RHO[] = {
    1.0,
};

#endif /* #ifdef RHO_PDF */

void *read_ftable(const char *file_path, const int n_rows, const int n_cols);
double rate_of_star_formation(const int index, double x);

#endif /* #ifdef EL_SFR_H */
	