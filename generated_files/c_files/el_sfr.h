#ifndef EL_SFR_H
#define EL_SFR_H

/* T [internal_units] * T_MYR = T [Myr] */
#define T_MYR (All.UnitTime_in_s / All.HubbleParam / SEC_PER_MEGAYEAR)
/* RHO [internal_units] * RHO_COSMO = RHO [cm^(-3)] */
#define RHO_COSMO (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv / PROTONMASS)
/* M [internal_units] * M_COSMO = M [Mₒ] */
#define M_COSMO (All.UnitMass_in_g / All.HubbleParam / SOLAR_MASS)
/* M [internal_units] * M_CGS = M [mp] */
#define M_CGS (All.UnitMass_in_g / All.HubbleParam / PROTONMASS)
/* L [internal_units] * L_CGS = L [cm] */
#define L_CGS (All.UnitLength_in_cm * All.cf_atime / All.HubbleParam)

/* Interpolation tables */
#define ETA_NROWS 107   // Number of rows in the η tables
#define ETA_NCOLS 101   // Number of columns in the η tables
#define ETA_LOG_AGE_MIN 5.000000000000000e+00 // First age of the η tables, log10(age [yr])
#define R_ZSN_NROWS 100 // Number of rows in the R and Zsn tables
#define R_ZSN_NCOLS 2   // Number of columns in the R and Zsn tables
#define UVB_NROWS 59    // Number of rows in the UVB table
#define UVB_NCOLS 2     // Number of columns in the UVB table

/* Paths */
static char *ETA_D_TABLE_PATH = "../code/src/el_sfr/tables/eta_d.txt";
static char *ETA_I_TABLE_PATH = "../code/src/el_sfr/tables/eta_i.txt";
static char *R_TABLE_PATH = "../code/src/el_sfr/tables/R.txt";
static char *ZSN_TABLE_PATH = "../code/src/el_sfr/tables/Zsn.txt";
static char *UVB_TABLE_PATH = "../code/src/el_sfr/tables/UVB.txt";

/* ODE constants */

/* Cρ = 100.0 (clumping factor) */
/* T_ceil = 40.0 K, μ = 4 / (1 + 3 * 0.76) (Jeans shielding length with temperature ceiling, Safranek-Shrader2017 model L1a) */
/* IMF: Chabrier2003 */
/* Yield model: Portinari1998 */
/* Stellar luminosity model: Millán-Irigoyen2021 */
/* LW background model: Incatasciato2023 */

#define N_EQU 6                         /* Number of equations */
#define ODE_CREC 8.204976000000000e+00  /* Recombination constant [Myr^(-1) * cm^3 * mp^(-1)] */
#define ODE_CCOND 4.149680315962227e+01 /* Condensation constant [Myr^(-1) * cm^3 * mp^(-1)] */
#define ODE_CS 1.942876283158012e-02    /* Star formation constant [Myr^(-1) * cm^(3/2) * mp^(-1/2)] */
#define ODE_CACCSI 4.908642141884553e-01 /* Silicate accretion constant [Myr^(-1) * mp^(-1) * cm^3] */
#define ODE_CACCC 1.820696248189202e-01  /* Graphite accretion constant [Myr^(-1) * mp^(-1) * cm^3] */
#define ODE_CSD 4.907796224212826e-19   /* Dust shielding constant [cm^2 * mp^(-1)] */
#define ODE_CSH2 9.999999999999999e-16             /* Molecular self-shielding constant [cm^2 * mp^(-1)] */
#define ODE_CLSH 2.760278182715268e+20     /* Shielding length constant [mp^(1/2) * cm^(-1/2)] */
#define ODE_CXSI 4.518072289156626e+00     /* Silicate initial condition constant [dimensionless] */
#define ODE_CXC 8.500000000000000e-01      /* Graphite initial condition constant [dimensionless] */
#define ODE_CMAXSI 6.024096385542168e+00   /* Maximum dust mass per unit mass of Si, 1 / f_Si [dimensionless] */
#define ODE_CMAXC 1.000000000000000e+00    /* Maximum dust mass per unit mass of C, 1 / f_C [dimensionless] */
#define Z_FLOOR 1.0e-10                    /* Metallicity below which a cell is considered pristine (numerical noise) [dimensionless] */
#define ODE_CTION 6.3000e-18            /* Photoionization optical depth constant [cm^2 * mp^(-1)] */
#define ODE_CTDISS 1.0500e-19           /* Photodissociation optical depth constant [cm^2 * mp^(-1)] */
#define FD_EFF 8.150297643300318e-06       /* Dust fraction floor for condensation, 1e-3 D_ref [dimensionless] */
#define WH2 2.0000e-01                  /* Molecular shielding parameter [dimensionless] */
#define ABEL97 4.354948800000000e+01              /* LWB dissociation constant [Myr^(-1)] */
#define DELTA_DUST 1.00e-01                /* Dust condensation efficiency [dimensionless] */
#define ODE_CSWMASS 7.443367288063504e+02  /* Nominal value of the mass swept-up by a SNe shock times the efficiency of dust destruction [Msun] */
#define SNII_FRAC_TEST 1.18e-02            /* Characteristic value of the fraction of SNII for testing purposes [Msun^(-1)] */

typedef struct DataTable
{
	double *data;  // Values of the table
	int n_rows;    // Number of rows in the table
	int n_cols;    // Number of columns in the table
} data_table;

void *read_ftable(const char *file_path, const int n_rows, const int n_cols);
double rate_of_star_formation(const int index);

#endif /* #ifdef EL_SFR_H */
	