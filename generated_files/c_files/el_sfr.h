#ifndef EL_SFR_H
#define EL_SFR_H

/* T [internal_units] * T_MYR = T [Myr] */
#define T_MYR (All.UnitTime_in_s / All.HubbleParam / SEC_PER_MEGAYEAR)
/* RHO [internal_units] * RHO_COSMO = RHO [cm^(-3)] */
#define RHO_COSMO (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv / PROTONMASS)
/* M [internal_units] * M_COSMO = M [Mₒ] */
#define M_COSMO (All.UnitMass_in_g / All.HubbleParam / SOLAR_MASS)
/* L [internal_units] * L_CGS = L [cm] */
#define L_CGS (All.UnitLength_in_cm * All.cf_atime / All.HubbleParam)

/* Interpolation tables */
#define ETA_NROWS 107 // Number of rows in the η tables
#define ETA_NCOLS 7   // Number of columns in the η tables
#define R_ZSN_NROWS 7 // Number of rows in the R and Zsn table
#define R_ZSN_NCOLS 2 // Number of columns in the R and Zsn table
#define UVB_NROWS 59  // Number of rows in the UVB table
#define UVB_NCOLS 2   // Number of columns in the UVB table

/* Paths */
static char *ETA_D_TABLE_PATH = "../code/src/el_sfr/tables/eta_d.txt";
static char *ETA_I_TABLE_PATH = "../code/src/el_sfr/tables/eta_i.txt";
static char *R_TABLE_PATH = "../code/src/el_sfr/tables/R.txt";
static char *ZSN_TABLE_PATH = "../code/src/el_sfr/tables/Zsn.txt";
static char *UVB_TABLE_PATH = "../code/src/el_sfr/tables/UVB.txt";

/* ODE constants */

/* Cρ = 100.0 (clumping factor) */

#define N_EQU 6                        /* Number of equations */
#define ODE_CR 8.204976000000000e+00   /* Recombination constant [Myr^(-1) * cm^3 * mp^(-1)] */
#define ODE_CC 1.739395275590551e+01   /* Condensation constant [Myr^(-1) * cm^3 * mp^(-1)] */
#define ODE_CS 1.942876283158012e-02   /* Star formation constant [Myr^(-1) * cm^(3/2) * mp^(-1/2)] */
#define INV_T_DD 4.356568364611260e-04 /* Inverse of the dust loss timescale [Myr^-1] */
#define ODE_CD 5.617360725623994e-04   /* Dust growth constant [Myr^(-1) * mp^(-1) * cm^3] */
#define ODE_CSD -3.149606299212598e-19 /* Dust shielding constant [cm^2 * mp^(-1)] */
#define ODE_CSH2 1.0000e-15            /* Molecular self-shielding constant [cm^2 * mp^(-1)] */
#define ODE_CXD 2.835953313674557e-01  /* Dust initial condition constant [dimensionless] */
#define ZEFF 1.2700e-05                /* Effective metallicity 1e-3 Zₒ */
#define WH2 2.0000e-01                 /* Molecular shielding parameter [dimensionless] */

typedef struct DataTable
{
	double *data;  // Values of the table
	int n_rows;    // Number of rows in the table
	int n_cols;    // Number of columns in the table
} data_table;

void *read_ftable(const char *file_path, const int n_rows, const int n_cols);
double rate_of_star_formation(const int index);

#endif /* #ifdef EL_SFR_H */
	