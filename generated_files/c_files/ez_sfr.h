#ifndef EZ_SFR_H
#define EZ_SFR_H
	
/* T [internal_units] * T_MYR = T [Myr] */
#define T_MYR (All.UnitTime_in_s / All.HubbleParam / SEC_PER_YEAR / 1000000)
/* RHO [internal_units] * RHO_COSMO = RHO [cm^(-3)] */
#define RHO_COSMO (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv / PROTONMASS)
/* M [internal_units] * M_COSMO = M [Mₒ] */
#define M_COSMO (All.UnitMass_in_g / SOLAR_MASS)
	
/* Interpolation tables */
#define ETA_NROWS 107  // Number of rows in the η tables
#define ETA_NCOLS 7  // Number of columns in the η tables
#define R_NROWS 7  // Number of rows in the R table
#define R_NCOLS 3  // Number of columns in the R table
/* Paths */
static char *ETA_D_TABLE_PATH = "../code/src/ez_sfr/tables/eta_d.txt";
static char *ETA_I_TABLE_PATH = "../code/src/ez_sfr/tables/eta_i.txt";
static char *R_TABLE_PATH     = "../code/src/ez_sfr/tables/R_Zsn.txt";
	
/* ODE constants */
#define N_EQU 4 /* Number of equations */
#define ODE_CS 2.5735041e+03  /* [Myr * cm^(-3/2)] */
#define ODE_CR 1.2187726e-01  /* [Myr * cm^(-3)] */
#define ODE_CC 5.7491245e+00  /* [Myr * cm^(-3)] */
#define ZEFF 1.2700e-05  /* 1e-3 Zₒ */
#define AW 0.00  /* Weight of the atomic fraction in the computation of the SFR */         
#define MW 1.00  /* Weight of the molecular fraction in the computation of the SFR */

typedef struct DataTable
{
  	double *data;  // Values of the table
  	int n_rows;    // Number of rows in the table
  	int n_cols;    // Number of columns in the table
} 	data_table;

void *read_ftable(const char *file_path, const int n_rows, const int n_cols);
double rate_of_star_formation(const int index);
	
#endif /* #ifdef EZ_SFR_H */
