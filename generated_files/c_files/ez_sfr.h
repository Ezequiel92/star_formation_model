#ifndef EZ_SFR_H
#define EZ_SFR_H

/* T [internal_units] * T_MYR = T [Myr] */
#define T_MYR (All.UnitTime_in_s / All.HubbleParam / SEC_PER_YEAR / 1000000)
/* RHO [internal_units] * RHO_COSMO = RHO [cm^(-3)] */
#define RHO_COSMO (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv / PROTONMASS)
/* M [internal_units] * M_COSMO = M [Mₒ] */
#define M_COSMO (All.UnitMass_in_g / SOLAR_MASS)

/* ODE error constants */
#define ABS_TOL 1.0e-8 /* Absolute tolerance */
#define REL_TOL 0.0    /* Relative tolerance */

/* η table constants */
#define ETA_NROWS 107  // Number of rows in the η table
#define ETA_NCOLS 7  // Number of columns in the η table
#define R_NROWS 7  // Number of rows in the R table
#define R_NCOLS 3  // Number of columns in the R table

/* ODE constants */
#define ODE_CS 2.5735041e+03  /* [Myr * cm^(-3/2)] */
#define ODE_CR 1.2234783e-01  /* [Myr * cm^(-3)] */
#define ODE_CC 3.5385031e+00  /* [Myr * cm^(-3)] */
#define ZEFF 1.3400e-05  /* 1e-3 Zₒ */
#define AW 0.00  /* Weight of the atomic fraction in the computation of the SFR */         
#define MW 1.00  /* Weight of the molecular fraction in the computation of the SFR */

/* Paths to the interpolation tables */
static char *ETA_D_TABLE = "../code/src/ez_sfr/tables/eta_d.txt";
static char *ETA_I_TABLE = "../code/src/ez_sfr/tables/eta_i.txt";
static char *R_TABLE = "../code/src/ez_sfr/tables/R_Zsn.txt";

#ifdef RHO_PDF

/*
 * Density PDF from Burkhart (2018) https://doi.org/10.3847/1538-4357/aad002
 *
 * We used the following parameters (all dimensionless):
 *
 * divisions = 20
 * range of log10(rho/rho_0) = (-6.0, 6.0)
 * α (power law slope) = 2.0
 * b (turbulent forcing parameter) = 0.5
 * Ms (mach number) = 10.0
 */
#define DIVISIONS 20
/* Integrated PDF of the interstellar gas density */
static const double PDF[] = {
	0.0107029097,
	0.0212985946,
	0.0379878324,
	0.0607274862,
	0.0870112819,
	0.1117421531,
	0.1286207158,
	0.1326960245,
	0.1227038406,
	0.1016977334,
	0.0755468486,
	0.0503004598,
	0.0300176109,
	0.0160555948,
	0.0076969705,
	0.0033071430,
	0.0012735674,
	0.0004395643,
	0.0001359721,
	0.0000376963,
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
/* Integrated probability density function of the interstellar gas density */
static const double PDF[] = {1.0,};
/* Density factor: ρ = ρ₀ * F_RHO */
static const double F_RHO[] = {1.0,};

#endif /* #ifdef RHO_PDF */

double *read_ftable(const char *filepath, const int n_rows, const int n_cols);
double rate_of_star_formation(const int index);

#endif /* #ifdef EZ_SFR_H */
