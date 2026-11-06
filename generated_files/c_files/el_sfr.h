#ifndef EL_SFR_H
#define EL_SFR_H

/* T [internal_units] * T_MYR = T [Myr] */
#define T_MYR (All.UnitTime_in_s / All.HubbleParam / SEC_PER_MEGAYEAR)
/* RHO [internal_units] * RHO_COSMO = RHO [mp * cm^(-3)] */
#define RHO_COSMO (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv / PROTONMASS)
/* M [internal_units] * M_COSMO = M [Msun] */
#define M_COSMO (All.UnitMass_in_g / All.HubbleParam / SOLAR_MASS)
/* M [internal_units] * M_CGS = M [mp] */
#define M_CGS (All.UnitMass_in_g / All.HubbleParam / PROTONMASS)

/* Interpolation tables */
#define ETA_NROWS 285   // Number of rows in the η tables
#define ETA_NCOLS 5     // Number of columns in the η tables
#define ETA_LOG_AGE_MIN 5.000000000000000e+00 // First age of the η tables, log10(age [yr])
#define R_ZSN_NROWS 105 // Number of rows in the R and Zsn tables
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

/* Model choices */
/* IMF: Chabrier2003, 0.15 M⊙ - 100.0 M⊙ (SNII progenitors above 8.0 M⊙) */
/* Yield model: Portinari1998 */
/* Stellar luminosity model: Millán-Irigoyen2021 */
/* UVB model: Haardt2012 */
/* LW background model: Incatasciato2023, J21 constant below z = 6.0 */
/* Clumping factor: Cρ = 100.0 */
/* Dust grains (condensation): a_eff = 0.0158 μm, s = 3.0 g cm^-3 */
/* Dust initial condition: ξ_Si = 0.25, ξ_C = 0.15 (Hirashita2012) */
/* Dust accretion: T_gas = 100.0 K, T_dust = 10.0 K */
/* Dust destruction: ϵ = 0.1, E_SNII = 1.09e51 erg */
/* Shielding length: T = 40.0 K */

#define N_EQU 6                          /* Number of equations */
#define ODE_CREC 8.204976000000000e+00   /* Recombination constant [Myr^(-1) * cm^3 * mp^(-1)] */
#define ODE_CCOND 4.175944115430344e+01  /* Condensation constant [Myr^(-1) * cm^3 * mp^(-1)] */
#define ODE_CS 1.942876283158012e-02     /* Star formation constant [Myr^(-1) * cm^(3/2) * mp^(-1/2)] */
#define ODE_CACCSI 4.939709497212937e-01 /* Silicate accretion constant [Myr^(-1) * cm^3 * mp^(-1)] */
#define ODE_CACCC 1.832219642165083e-01  /* Graphite accretion constant [Myr^(-1) * cm^3 * mp^(-1)] */
#define ODE_CSD 4.907796224212826e-19    /* Dust shielding constant [cm^2 * mp^(-1)] */
#define ODE_CSH2 9.999999999999999e-16   /* Molecular self-shielding constant [cm^2 * mp^(-1)] */
#define ODE_CLSH 2.760278182715268e+20   /* Shielding length constant [mp^(1/2) * cm^(-1/2)] */
#define ODE_CXSI 4.518072289156626e+00   /* Dust mass per unit mass of Si in the initial condition [dimensionless] */
#define ODE_CXC 8.500000000000000e-01    /* Dust mass per unit mass of C in the initial condition [dimensionless] */
#define ODE_CMAXSI 6.024096385542168e+00 /* Maximum dust mass per unit mass of Si [dimensionless] */
#define ODE_CMAXC 1.000000000000000e+00  /* Maximum dust mass per unit mass of C [dimensionless] */
#define Z_FLOOR 1.0e-10                  /* Metallicity below which a cell is considered pristine [dimensionless] */
#define ODE_CTION 6.3000e-18             /* Photoionization optical depth constant [cm^2 * mp^(-1)] */
#define ODE_CTDISS 1.0500e-19            /* Photodissociation optical depth constant [cm^2 * mp^(-1)] */
#define FD_EFF 8.150297643300318e-06     /* Dust fraction floor for condensation [dimensionless] */
#define WH2 2.0000e-01                   /* Molecular self-shielding parameter [dimensionless] */
#define ODE_EXPFAC 8.500000000000000e-04 /* Exponential factor of the molecular self-shielding [dimensionless] */
#define ABEL97 4.354948800000000e+01     /* H2 photodissociation rate per unit J21 [Myr^(-1)] */
#define LWB_A 2.119000000000000e+00      /* J21 fit, constant coefficient [dimensionless] */
#define LWB_B -1.117000000000000e-01     /* J21 fit, linear coefficient [dimensionless] */
#define LWB_C -2.782000000000000e-03     /* J21 fit, quadratic coefficient [dimensionless] */
#define LWB_Z_MIN 6.000000000000000e+00  /* Lowest redshift of the validity range of the J21 fit [dimensionless] */
#define DELTA_DUST 1.00e-01              /* Dust condensation efficiency of the SN ejecta [dimensionless] */
#define ODE_CDD 9.215105644801868e+00    /* Dust destruction constant [dimensionless] */

typedef struct DataTable
{
	double *data;  // Values of the table
	int n_rows;    // Number of rows in the table
	int n_cols;    // Number of columns in the table
} data_table;

void *read_ftable(const char *file_path, const int n_rows, const int n_cols);
double rate_of_star_formation(const int index);

#endif /* #ifndef EL_SFR_H */
	