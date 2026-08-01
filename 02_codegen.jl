### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 65c19790-a9a8-11ee-06de-ab8cb360a743
let
	using Printf, PlutoLinks, Libdl

	using DataFrames, ChaosTools, CSV, DifferentialEquations, Interpolations, LinearAlgebra, Measurements, NaNMath, PlutoUI, QuadGK, SciMLLogging, SpecialFunctions, Symbolics, TikzPictures, Unitful, UnitfulAstro
end

# ╔═╡ 37653018-aaf5-42d1-a938-e05a44f18918
# ╠═╡ skip_as_script = true
#=╠═╡
TableOfContents(title="Code generation", depth=4)
  ╠═╡ =#

# ╔═╡ f028697a-0534-4ccd-8282-1b9a3d8e284b
# ╠═╡ skip_as_script = true
#=╠═╡
md"# Code generation"
  ╠═╡ =#

# ╔═╡ 4b65fd9c-bd69-4e97-a28c-357c6a9c7bda
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Load the model"
  ╠═╡ =#

# ╔═╡ a0ca487a-8a51-48cb-acc0-7b318c793d09
const MODEL = @ingredients("./01_model.jl");

# ╔═╡ 176d2173-b440-4d03-822b-60034921cade
# If we will use an ODE solver from GSL that needs the jacobian
const JACOBIAN = false;

# ╔═╡ 33e83e53-1df4-48a9-85fd-17b50c1b1be6
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Paths"
  ╠═╡ =#

# ╔═╡ 0c6f9d84-30d5-4b4d-8efd-f7414205aa5e
begin
	const GEN_FILES = mkpath("./generated_files")

	# Path to the .c and .h tables
	const C_FILES = mkpath(joinpath(GEN_FILES, "c_files"))

	# Paths to the interpolation tables
	const INTERP_FILES = mkpath(joinpath(GEN_FILES, "interpolation_tables"))
	const ETA_D_TABLE  = joinpath(INTERP_FILES, "eta_d.txt")
	const ETA_I_TABLE  = joinpath(INTERP_FILES, "eta_i.txt")
	const R_TABLE      = joinpath(INTERP_FILES, "R.txt")
	const ZSN_TABLE    = joinpath(INTERP_FILES, "Zsn.txt")
	const UVB_TABLE    = joinpath(INTERP_FILES, "UVB.txt")
	const SDSS_TABLE   = joinpath(INTERP_FILES, "magnitudes.txt")

	# Paths to the C libraries (DLLs or SOs) and GCC
	@static if Sys.iswindows()
	    const LIB_PATH = joinpath(C_FILES, "lib.dll")
		const GCC_PATH = "C:/msys64/mingw64/bin/gcc.exe"
		const GSL_INC = "C:/msys64/mingw64/include" 
		const GSL_LIB = "C:/msys64/mingw64/lib"
	elseif Sys.islinux()
	    const LIB_PATH = joinpath(C_FILES, "lib.so")
		const GCC_PATH = "gcc"
		const GSL_INC = "/home/linuxbrew/.linuxbrew/include"
		const GSL_LIB = "/home/linuxbrew/.linuxbrew/lib"
	end
end;

# ╔═╡ 5eadbd03-da0b-48e1-9027-bc3243389a8d
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Interpolation tables"
  ╠═╡ =#

# ╔═╡ 0f03beab-d57c-4e7c-8dd1-e5eb3425c0c2
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Tables sizes"
  ╠═╡ =#

# ╔═╡ 5769c1c3-d184-41b2-9f23-d490fb1a9df1
begin
	# Stellar ages for the η_diss and η_ion interpolation tables
	η_ages = MODEL.Q_ages

	# Metallicities for the η_diss and η_ion interpolation tables
	η_Zs = collect(range(MODEL.Q_metals[1], MODEL.Q_metals[end], 100))

	# Metallicities for the R and Zsn interpolation tables
	recycling_Zs = collect(range(MODEL.sy_metals[1], MODEL.sy_metals[end], 100))

	# Number of rows in the η tables
	const ETA_NROWS = length(η_ages) + 1   
	
	# Number of columns in the η tables
	const ETA_NCOLS = length(η_Zs) + 1  

	# Number of rows in the R and Zsn tables
	const R_ZSN_NROWS = length(recycling_Zs) 

	# Number of columns in the R and Zsn tables
	const R_ZSN_NCOLS = 2   

	# Number of rows in the UVB table
	const UVB_NROWS = size(MODEL.UVB_TABLE, 1)    

	# Number of columns in the UVB table
	const UVB_NCOLS = 2    
end;

# ╔═╡ 798cafe1-58e6-4c30-8668-59615b1d26ba
# ╠═╡ skip_as_script = true
#=╠═╡
md"### UVB"
  ╠═╡ =#

# ╔═╡ b94ba6f9-f136-4ca9-bb21-207df7a5534b
# ╠═╡ skip_as_script = true
#=╠═╡
# We use the table from Haardt et al. (2012)
cp(
	joinpath(GEN_FILES, "../data/photoionization_rate_hm12.txt"),
	UVB_TABLE;
	force=true,
);
  ╠═╡ =#

# ╔═╡ 65b56c83-f48a-4df5-8e52-326e40c56e78
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Photodissociation and photoionization"
  ╠═╡ =#

# ╔═╡ 0db83126-5520-4438-9f6b-2a6839559787
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Write the interpolation tables for η_diss(age, Z) and η_ion(age, Z)
#
# Each file is named:
# 
#     η_diss(age, Z) -> eta_d.txt
#     η_ion(age, Z)  -> eta_i.txt
# 
# The columns are the tellar metallicity, the rows are the stellar ages
# 
# Arguments
# 
#     path: Where to store the resulting text files
#     ages: List of stellar ages, as log₁₀(age [y])
#     Zs:   List of stellar metallicities
#####################################################################################
function write_η_tables(
	path::String,
	log_ages::Vector{Float64},
	Zs::Vector{Float64},
)::Nothing

	# Allocate memory for η_diss and η_ion matrices
	η_diss = Matrix{Float64}(undef, length(log_ages) + 1, length(Zs) + 1)
	η_ion = Matrix{Float64}(undef, length(log_ages) + 1, length(Zs) + 1)

	# First column
	η_diss[:, 1] .= [0.0, log_ages...]
	η_ion[:, 1]  .= [0.0, log_ages...]
	
	# First row
	η_diss[1, 2:end] .= Zs
	η_ion[1, 2:end]  .= Zs

	for (i, log_age) in pairs(log_ages)

		for (j, Z) in pairs(Zs)

			age = ustrip(u"Myr", exp10(log_age) * u"yr")

			η_diss[i + 1, j + 1], η_ion[i + 1, j + 1] = MODEL.compute_η(age, Z)

		end

	end

	# Write `eta_d.txt`
	CSV.write(
		joinpath(path, "eta_d.txt"), 
		Tables.table(η_diss); 
		delim=' ',
		writeheader=false,
	)

	# Write `eta_i.txt`
	CSV.write(
		joinpath(path, "eta_i.txt"), 
		Tables.table(η_ion); 
		delim=' ',
		writeheader=false,
	)

	return nothing

end;
  ╠═╡ =#

# ╔═╡ 93d71685-815f-4faf-bff8-3bb35945d2bf
# ╠═╡ skip_as_script = true
#=╠═╡
write_η_tables(INTERP_FILES, η_ages, η_Zs)
  ╠═╡ =#

# ╔═╡ 5db32b26-0485-4929-89ec-34c09450555e
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Mass recycling"
  ╠═╡ =#

# ╔═╡ 871a53c5-82aa-4a4e-9627-c759901e7a89
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Write the interpolation tables for R(Z) and Zsn(Z)
# 
# Each file is named:
# 
#     R(Z)   -> R.txt
#     Zsn(Z) -> Zsn.txt
# 
# The columns are:
# 
#     R.txt   -> Z | R
#     Zsn.txt -> Z | Zsn
# 
# Arguments
# 
#     path: Where to store the resulting text files
#     Zs:   List of stellar metallicities
#####################################################################################
function write_recycling_tables(path::String, Zs::Vector{Float64})::Nothing

	# Allocate memory for R and Zsn
	R   = similar(Zs)
	Zsn = similar(Zs)

	for (i, Z) in pairs(Zs)
	    R[i], Zsn[i] = MODEL.compute_recycled_fractions(Z)
	end

	# Write `R.txt`
	CSV.write(
		joinpath(path, "R.txt"), 
		Tables.table([Zs;; R]); 
		delim=' ', 
		writeheader=false,
	)

	# Write `Zsn.txt`
	CSV.write(
		joinpath(path, "Zsn.txt"), 
		Tables.table([Zs;; Zsn]); 
		delim=' ', 
		writeheader=false,
	)

	return nothing

end;
  ╠═╡ =#

# ╔═╡ 76d42a07-ae3f-42f9-8a97-9b2892e55088
# ╠═╡ skip_as_script = true
#=╠═╡
write_recycling_tables(INTERP_FILES, recycling_Zs)
  ╠═╡ =#

# ╔═╡ 95fc4af6-247b-43fd-b271-5f3a588f3dff
# ╠═╡ skip_as_script = true
#=╠═╡
md"### SDSS magnitudes"
  ╠═╡ =#

# ╔═╡ 5fa18c9e-0ae1-402c-aa28-7396a4973105
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Write the interpolation tables for the magnitudes of the g, r, and i SDSS filters
# 
# The columns are:
# 
#     Z | log_age | g | r | i
# 
#       Z:       Metallicty
#       log_age: log₁₀(age [yr])
#       g:       Blue filter
#       r:       Green filter
#       i:       Red filter
# 
# Arguments
# 
#     path: Where to store the resulting text files
#     file: Path to the source file with the magnitudes (HR-pyPopStar 2025)
#     imf:  Target IMF
#####################################################################################
function write_magnitude_table(path::String, file::String, imf::String)::Nothing

	# Load the data
	data = CSV.read(
		file, 
		DataFrame; 
		delim=" ",
		ignorerepeated=true,
		header=[:IMF, :Z, :log_age, :U, :B1, :B2, :V, :R, :u, :g, :r, :i, :z], 
		skipto=2, 
		types=[String, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64]
	)

	# Select the target IMF and the relevant columns
	magnitudes = data[data.IMF .== uppercase(imf), [:Z, :log_age, :g, :r, :i]]

	CSV.write(
		joinpath(path, "magnitudes.txt"), 
		magnitudes; 
		delim=' ', 
		writeheader=false,
	)
	
	return nothing

end;
  ╠═╡ =#

# ╔═╡ 22361146-afa5-40a5-ad79-9b7b91cba8db
# ╠═╡ skip_as_script = true
#=╠═╡
write_magnitude_table(
	INTERP_FILES, 
	MODEL.MAGNITUDE_TABLE, 
	MODEL.IMF_FUNCTIONS[MODEL.IMF][1],
)
  ╠═╡ =#

# ╔═╡ 73503ce5-c8b4-4f1d-a95c-cc7a307fdb8f
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Header file"
  ╠═╡ =#

# ╔═╡ 2017166e-2932-48f3-a9e2-f5423f018742
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Write the header file
# 
# Arguments
# 
#     path: Where to store the header file
# 
# Returns
# 
#     If path === nothing, return the header file as a string
#####################################################################################
function write_c_header(path::Union{String,Nothing})::Union{String,Nothing}

	header = """
#ifndef EL_SFR_H
#define EL_SFR_H

/* T [internal_units] * T_MYR = T [Myr] */
#define T_MYR (All.UnitTime_in_s / All.HubbleParam / SEC_PER_MEGAYEAR)
/* RHO [internal_units] * RHO_COSMO = RHO [cm^(-3)] */
#define RHO_COSMO (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv / PROTONMASS)
/* M [internal_units] * M_COSMO = M [M☉] */
#define M_COSMO (All.UnitMass_in_g / All.HubbleParam / SOLAR_MASS)
/* M [internal_units] * M_CGS = M [mp] */
#define M_CGS (All.UnitMass_in_g / All.HubbleParam / PROTONMASS)
/* L [internal_units] * L_CGS = L [cm] */
#define L_CGS (All.UnitLength_in_cm * All.cf_atime / All.HubbleParam)

/* Interpolation tables */
#define ETA_NROWS $(ETA_NROWS)   // Number of rows in the η tables
#define ETA_NCOLS $(ETA_NCOLS)   // Number of columns in the η tables
#define R_ZSN_NROWS $(R_ZSN_NROWS) // Number of rows in the R and Zsn tables
#define R_ZSN_NCOLS $(R_ZSN_NCOLS)   // Number of columns in the R and Zsn tables
#define UVB_NROWS $(UVB_NROWS)    // Number of rows in the UVB table
#define UVB_NCOLS $(UVB_NCOLS)     // Number of columns in the UVB table

/* Paths */
static char *ETA_D_TABLE_PATH = "../code/src/el_sfr/tables/eta_d.txt";
static char *ETA_I_TABLE_PATH = "../code/src/el_sfr/tables/eta_i.txt";
static char *R_TABLE_PATH = "../code/src/el_sfr/tables/R.txt";
static char *ZSN_TABLE_PATH = "../code/src/el_sfr/tables/Zsn.txt";
static char *UVB_TABLE_PATH = "../code/src/el_sfr/tables/UVB.txt";

/* ODE constants */

/* Cρ = $(@sprintf("%.1f", MODEL.Cρ)) (clumping factor) */
/* R⊙ = $(@sprintf("%.3e", ustrip(u"cm^3 * s^-1", MODEL.Rsun))) cm^3 * s^-1 (formation rate coefficient of H2 on dust grain, at solar metallicity) */
/* IMF: $(MODEL.IMF) */
/* Yield model: $(MODEL.YIELD_MODEL) */
/* Stellar luminosity model: $(MODEL.HR_LUMINOSITY ? "Millán-Irigoyen2021" : "Mollá2009") */
/* LW background model: $(MODEL.LW_SOURCE) */

#define N_EQU $(MODEL.N_EQU)                            /* Number of equations */
#define ODE_CREC $(@sprintf("%.15e", MODEL.c_rec))     /* Recombination constant [Myr^(-1) * cm^3 * mp^(-1)] */
#define ODE_CCOND $(@sprintf("%.15e", MODEL.c_cond))    /* Condensation constant [Myr^(-1) * cm^3 * mp^(-1)] */
#define ODE_CS $(@sprintf("%.15e", MODEL.c_star))       /* Star formation constant [Myr^(-1) * cm^(3/2) * mp^(-1/2)] */
#define ODE_CDG $(@sprintf("%.15e", MODEL.c_dg))      /* Dust growth constant [Myr^(-1) * mp^(-1) * cm^3] */
#define ODE_CSD $(@sprintf("%.15e", MODEL.c_sd))      /* Dust shielding constant [cm^2 * mp^(-1)] */
#define ODE_CSH2 $(@sprintf("%.4e", MODEL.c_sh2))                /* Molecular self-shielding constant [cm^2 * mp^(-1)] */
#define ODE_CXD $(@sprintf("%.15e", MODEL.c_xd))      /* Dust initial condition constant [dimensionless] */
#define ODE_CTION $(@sprintf("%.4e", MODEL.c_τion))               /* Photoionization optical depth constant [cm^2 * mp^(-1)] */
#define ODE_CTDISS $(@sprintf("%.4e", MODEL.c_τdiss))              /* Photodissociation optical depth constant [cm^2 * mp^(-1)] */
#define ZEFF $(@sprintf("%.4e", MODEL.Zeff))                    /* Effective metallicity 1e-3 Zₒ */
#define WH2 $(@sprintf("%.4e", MODEL.ωH2))                     /* Molecular shielding parameter [dimensionless] */
#define ABEL97 $(@sprintf("%.4e", MODEL.abel97))                  /* LWB dissociation constant [Myr^(-1)] */
#define DELTA_DUST $(@sprintf("%.2e", MODEL.δ_d))                /* Dust condensation efficiency [dimensionless] */
#define ODE_CSWMASS $(@sprintf("%.15e", MODEL.c_sw_mass))  /* Nominal value of the mass swept-up by a SNe shock times the efficiency of dust destruction [Msun] */
#define SNII_FRAC_TEST $(@sprintf("%.2e", MODEL.SNII_FRAC))            /* Characteristic value of the fraction of SNII for testing purposes [Msun^(-1)] */

typedef struct DataTable
{
	double *data;  // Values of the table
	int n_rows;    // Number of rows in the table
	int n_cols;    // Number of columns in the table
} data_table;

void *read_ftable(const char *file_path, const int n_rows, const int n_cols);
double rate_of_star_formation(const int index);

#endif /* #ifdef EL_SFR_H */
	"""

	if isnothing(path)
		return header
	else
		mkpath(dirname(path))

		open(path, "w") do file
			write(file, header)
		end

		return nothing
	end

end;
  ╠═╡ =#

# ╔═╡ be89676f-0e5a-48c5-b8e3-61b2d7b892a5
# ╠═╡ skip_as_script = true
#=╠═╡
write_c_header(joinpath(C_FILES, "el_sfr.h"))
  ╠═╡ =#

# ╔═╡ 7e5980f4-2bcc-4c42-a244-5a8e45df9ae1
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Jacobian"
  ╠═╡ =#

# ╔═╡ 85d9b763-ded9-4e95-be86-e519abe904d9
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Write the jacobian as a C function in the `el_sfr.c` file
# 
# Arguments
# 
#     path: Path to the `el_sfr.c` file
#####################################################################################
function write_jacobian(path::String)::Nothing

	# Boilerplate
	head = """
	\n/*! \\brief Evaluate the Jacobian of the systems of equations.
	*
	*  Evaluate the Jacobian matrix of the model, using the following variables:
	*
	*  Ionized gas fraction:    fi(t) = Mi(t) / MC --> y[0]
	*  Atomic gas fraction:     fa(t) = Ma(t) / MC --> y[1]
	*  Molecular gas fraction:  fm(t) = Mm(t) / MC --> y[2]
	*  Stellar fraction:        fs(t) = Ms(t) / MC --> y[3]
	*  Metal fraction:          fZ(t) = MZ(t) / MC --> y[4]
	*  Dust fraction:           fd(t) = Md(t) / MC --> y[5]
	*
	*  where MC = Mi(t) + Ma(t) + Mm(t) + Ms(t) + MZ(t) + Md(t) is the total density of the gas cell,
	*  and each equation has units of Myr^(-1).
	*
	*  \\param[in] t Unused variable to comply with the `gsl_odeiv2_driver_alloc_y_new()` API.
	*  \\param[in] y Values of the variables at which the Jacobian will be evaluated.
	*  \\param[out] dfdy Where the results of evaluating the Jacobian will be stored.
	*  \\param[out] dfdt Where the results of evaluating the time derivatives will be stored.
	*  \\param[in] parameters Parameters for the Jacobian.
	*
	*  \\return Constant `GSL_SUCCESS`, to confirm that the computation was successful.
	*/
	static int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *parameters)
	{
		(void)(t);

		/*
		 * Destructure the parameters
		 *
		 * rho_C: Total cell density           [mp * cm^(-3)]
		 * UVB:   UVB photoionization rate     [Myr^(-1)]
		 * LWB:   LWB Photodissociation  rate  [Myr^(-1)]
		 * eta_d: Photodissociation efficiency [dimensionless]
		 * eta_i: Photoionization efficiency   [dimensionless]
		 * R:     Mass recycling fraction      [dimensionless]
		 * Zsn:   Metals recycling fraction    [dimensionless]
		 * h:     Column height                [cm]
		 */
		double *p    = (double *)parameters;
		double rho_C = p[0];
		double UVB   = p[1];
		double LWB   = p[2];
		double eta_d = p[3];
		double eta_i = p[4];
		double R     = p[5];
		double Zsn   = p[6];
		double h     = p[7];

		gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, $(MODEL.N_EQU), $(MODEL.N_EQU));
		gsl_matrix *m = &dfdy_mat.matrix;

		double aux_var_01 = AUX_VAR_01;
		double aux_var_02 = AUX_VAR_02;
		double aux_var_03 = AUX_VAR_03;
		double aux_var_04 = AUX_VAR_04;
		double aux_var_05 = AUX_VAR_05;
		double aux_var_06 = AUX_VAR_06;
		double aux_var_07 = AUX_VAR_07;
		double aux_var_08 = AUX_VAR_08;
		double aux_var_09 = AUX_VAR_09;
		double aux_var_10 = AUX_VAR_10;
		double aux_var_11 = AUX_VAR_11;
		double aux_var_12 = AUX_VAR_12;
		double aux_var_13 = AUX_VAR_13;

	"""

	tail = """
		dfdt[0] = 0;
		dfdt[1] = 0;
		dfdt[2] = 0;
		dfdt[3] = 0;
		dfdt[4] = 0;
		dfdt[5] = 0;

		return GSL_SUCCESS;
	}
	"""

	c_file_str = read(path, String)

	# Main regex patterns
	jac_block    = r"(?<=// JACOBIAN_START)((?s:.)*?)(?=// JACOBIAN_END)"
	solver_block = r"(?<=gsl_odeiv2_driver_alloc_y_new\(\&sys)((?s:.)*?)(?=;)"
	sys_block    = r"(?<=gsl_odeiv2_system sys = )((?s:.)*?)(?=;)"

	if JACOBIAN

		# Create the C version of the jacobian
		@variables S_t S_ic[1:MODEL.N_EQU] S_parameters[1:MODEL.N_PAR]
	    S_dydt = Vector{Num}(undef, MODEL.N_EQU)
	    MODEL.system!(S_dydt, S_ic, S_parameters, S_t)

		# Compute the Jacobian symbolically
	    jac = Symbolics.jacobian(S_dydt, S_ic)

		# Function block regex pattern
		func_block = r"(?<=\{\n  )(.*?)(?=\n\}\n)"

		matrix = ""

		for i in 1:MODEL.N_EQU

			for j in 1:MODEL.N_EQU

				# Transform the symbolic expresions into C functions
				c_function = build_function(
					jac[i, j],
					S_ic,
					S_parameters,
					S_t;
					target=Symbolics.CTarget(),
				)

				# Replacements for correct formatting
				matrix *= replace(
					 match(func_block, c_function).match,
					"du[0] =" => "\tgsl_matrix_set(m, $(i-1), $(j-1),",
					 ";"      => ");\n",
				)

			end

		end

		aux_var_01 = "sqrt(rho_C)"
		aux_var_02 = "sqrt(1.0 + 9.999999999999999e-16 * y[2] * rho_C * h)"
		aux_var_03 = "exp(-3.149606299212598e-19 * (y[1] + y[2]) * (y[4] + y[5]) * rho_C * h)"
		aux_var_04 = "expm1(-6.2999999999999996e-18 * y[1] * rho_C * h)"
		aux_var_05 = "pow(1.0 + 9.999999999999999e-16 * y[2] * rho_C * h, 2)"
		aux_var_06 = "expm1(-1.05e-19 * y[2] * rho_C * h)"
		aux_var_07 = "exp(-0.00085 * aux_var_02)"
		aux_var_08 = "pow(aux_var_02, 2)"
		aux_var_09 = "pow(1.0 + 9.999999999999999e-16 * y[2] * rho_C * h, 4)"
		aux_var_10 = "rho_C * eta_i * h * aux_var_04 * aux_var_01 * aux_var_03"
		aux_var_11 = "exp(-6.2999999999999996e-18 * y[1] * rho_C * h)"
		aux_var_12 = "pow(y[1] + y[2], 2)"
		aux_var_13 = "exp(-1.05e-19 * y[2] * rho_C * h) * aux_var_03"
	    
		jacobian_string = head * matrix * "\n" * tail

		# Replacements for correct formatting
		jacobian_string = replace(
			jacobian_string,
			"RHS1"    => "y",
			"RHS2[0]" => "rho_C",
			"RHS2[1]" => "UVB",
			"RHS2[2]" => "LWB",
			"RHS2[3]" => "eta_d",
			"RHS2[4]" => "eta_i",
			"RHS2[5]" => "R",
			"RHS2[6]" => "Zsn",
			"RHS2[7]" => "h",
			"+ -"     => "-",
		)

		# Replacements for correct formatting
		jacobian_string = replace(
			jacobian_string,
			"-1 *"       => "-",
			aux_var_01   => "aux_var_01",
			aux_var_02   => "aux_var_02",
			aux_var_03   => "aux_var_03",
			aux_var_04   => "aux_var_04",
			aux_var_05   => "aux_var_05",
			aux_var_06   => "aux_var_06",
			"AUX_VAR_01" => aux_var_01,
			"AUX_VAR_02" => aux_var_02,
			"AUX_VAR_03" => aux_var_03,
			"AUX_VAR_04" => aux_var_04,
			"AUX_VAR_05" => aux_var_05,
			"AUX_VAR_06" => aux_var_06,
		)

		# Replacements for correct formatting
		jacobian_string = replace(
			jacobian_string,
			aux_var_07   => "aux_var_07",
			aux_var_08   => "aux_var_08",
			aux_var_09   => "aux_var_09",
			aux_var_10   => "aux_var_10",
			aux_var_11   => "aux_var_11",
			aux_var_12   => "aux_var_12",
			aux_var_13   => "aux_var_13",
			"AUX_VAR_07" => aux_var_07,
			"AUX_VAR_08" => "aux_var_02 * aux_var_02",
			"AUX_VAR_09" => aux_var_09,
			"AUX_VAR_10" => aux_var_10,
			"AUX_VAR_11" => aux_var_11,
			"AUX_VAR_12" => aux_var_12,
			"AUX_VAR_13" => aux_var_13,
		)

		# Write the jacobian into the main string
		file_with_jacobian = replace(
			c_file_str,
			jac_block    => jacobian_string,
			solver_block => ", gsl_odeiv2_step_msbdf, it * 1e-4, 1e-10, 0.0)",
			sys_block    => "{sf_ode, jacobian, N_EQU, parameters}",
		)

	else

		file_with_jacobian = replace(
			c_file_str,
			jac_block    => "\n",
			solver_block => ", gsl_odeiv2_step_msadams, it * 1e-4, 1e-10, 0.0)",
			sys_block    => "{sf_ode, NULL, N_EQU, parameters}",
		)

	end

	# Write the jacobian into the C file
	write(path, file_with_jacobian)
		
	return nothing
	
end;
  ╠═╡ =#

# ╔═╡ 8d718f6f-1ace-4d5a-8c9d-6b1994cb014d
# ╠═╡ skip_as_script = true
#=╠═╡
write_jacobian(joinpath(C_FILES, "el_sfr.c"))
  ╠═╡ =#

# ╔═╡ 2e4eda53-ef7b-411e-893d-c47de1f4b4df
md"## LWB photodissociation"

# ╔═╡ 8aea8d2e-90ef-4302-875d-d99f952a6ca0
#####################################################################################
# Write the choseen LWB photodissociation into the `el_sfr.c` file
# 
# Arguments
# 
#     path: Path to the `el_sfr.c` file
#####################################################################################
function write_LWB(path::String)::Nothing

	c_file_str = read(path, String)

	# Main regex patterns
	A_block = r"(?<=double LWB_A = )((?s:.)*?)(?=;)"
	B_block = r"(?<=double LWB_B = )((?s:.)*?)(?=;)"
	C_block = r"(?<=double LWB_C = )((?s:.)*?)(?=;)"

	file_with_LWB = replace(
		c_file_str,
		A_block => MODEL.lwb_A,
		B_block => MODEL.lwb_B,
		C_block => MODEL.lwb_C,
	)

	# Write the LWB photodissociation into the C file
	write(path, file_with_LWB)
		
	return nothing
	
end;

# ╔═╡ 8c7ed540-9e7d-428f-bb68-ae7fe91fef42
write_LWB(joinpath(C_FILES, "el_sfr.c"))

# ╔═╡ e36e8f05-138f-4a8f-be3a-23edbac61304
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Dynamic libraries"
  ╠═╡ =#

# ╔═╡ a1e9bf3b-8bf9-4425-a57c-61876cfdaa7a
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Compile `el_sfr.c` as a dynamic C library for Windows (DLL) or Linux (SO)
# 
# Arguments
# 
#     path: Path to the folder with the `el_sfr.c` and `el_sfr.h` files
#####################################################################################
function compile_libraries(path::String)::Nothing

	opt_cmd = `$(GCC_PATH) -Wall -Wno-unused-variable -fpic -shared -O3 -march=native -mtune=native -flto`
	in_out_cmd = `$(path)/el_sfr.c -o $(LIB_PATH)`
	gsl_cmd = `-I$(GSL_INC) -L$(GSL_LIB) -Wl,-rpath,$(GSL_LIB) -lgsl -lgslcblas -lm`

	run(`$(opt_cmd) -D EL_SFR -D TESTING $(in_out_cmd) $(gsl_cmd)`)

	return nothing

end;
  ╠═╡ =#

# ╔═╡ c90f69cf-399e-4831-818d-f6d34036b641
# ╠═╡ skip_as_script = true
#=╠═╡
compile_libraries(C_FILES)
  ╠═╡ =#

# ╔═╡ d4af2366-fa8d-4280-9057-15699d585daa
# ╠═╡ skip_as_script = true
#=╠═╡
md"### GSL ODE solver"
  ╠═╡ =#

# ╔═╡ 6dbbe165-1177-405d-a3cf-661cf96b265e
#####################################################################################
# Integrate the ODEs using the GNU scientific library (GSL)
# 
# Arguments
# 
#     eta_d_table: Path to the η_diss interpolation table
#     eta_i_table: Path to the η_ion interpolation table
#     R_table:     Path to the R interpolation table
#     Zsn_table:   Path to the Zsn interpolation table
#     UVB_table:   Path to the UVB interpolation table
#     library:     Pointer to the dynamic library that has the C functions
#                  void integrate_ode(
#                      double *ic, 
#                      double *parameters, 
#                      double it,
#                  )
#                  double interpolate1D(
#                      double x, 
#                      const char *table_path, 
#                      const int NROWS, 
#                      const int NCOLS,
#                  )
#                  double interpolate2D(
#                      double x, 
#                      double y, 
#                      const char *table_path,
#                  )
#                  double J21(double z)
# 
# Returns
# 
#     An integration function that will use the C routines in `library`
#####################################################################################
function integrate_with_c(
	eta_d_table::String,
	eta_i_table::String,
	R_table::String,
	Zsn_table::String,
	UVB_table::String,
	library::Ptr{Nothing},
)::Function

	# Load the C functions rom `library`
	integrate_ode = Libdl.dlsym(library, :integrate_ode)
	interpolate1D = Libdl.dlsym(library, :interpolate1D)
	interpolate2D = Libdl.dlsym(library, :interpolate2D)
	J21           = Libdl.dlsym(library, :J21)

	# Construct the integration function
	#
	# ICs (`ic`):
	#
	#     fi: Ionized gas fraction   [dimensionless]
	#     fa: Atomic gas fraction    [dimensionless]
	#     fm: Molecular gas fraction [dimensionless]
	#     fs: Stellar fraction       [dimensionless]
	#     fZ: Metal fraction         [dimensionless]
	#     fd: Dust fraction          [dimensionless]
	#
	# Parameters (`base_params`):
	#
	#     ρ_cell: Total cell density [mp * cm^(-3)]
 	#     Z:      Metallicity        [dimensionless]
	#     a:      Scale factor       [dimensionless]
	#     h:      Column height      [cm]
	#
	# Integration time (`it`):
	#
	#     it: Integration time [Myr]
	function integration(
		ic::Vector{Float64},
		base_params::Vector{Float64},
		it::Float64,
	)::Vector{Float64}

		fractions = copy(ic)

		# ODE parameters
		ρ_cell, Z, a, h = base_params

		# Redshift
		z = (1.0 / a) - 1.0

		# Integration time as log10(it [yr]) for the interpolation tables
		log_it = log10(ustrip(u"yr", it * u"Myr"))
		
		# UVB photoionization rate [Myr^(-1)]
		ΓUVB = @ccall $interpolate1D(
			z::Cdouble,
			UVB_table::Cstring,
			UVB_NROWS::Cint,
			UVB_NCOLS::Cint,
		)::Cdouble

		# LWB photodissociation rate [Myr^(-1)]
		ΓLWB = MODEL.abel97 * @ccall $J21(z::Cdouble)::Cdouble

		# Photodissociation efficiency [dimensionless]
		η_diss = @ccall $interpolate2D(
			log_it::Cdouble,
			Z::Cdouble,
			eta_d_table::Cstring,
		)::Cdouble

		# Photoionization efficiency [dimensionless]
		η_ion = @ccall $interpolate2D(
			log_it::Cdouble,
			Z::Cdouble,
			eta_i_table::Cstring,
		)::Cdouble

		# Mass recycling fraction [dimensionless]
		R = @ccall $interpolate1D(
			Z::Cdouble,
			R_table::Cstring,
			R_ZSN_NROWS::Cint,
			R_ZSN_NCOLS::Cint,
		)::Cdouble

		# Metals recycling fraction [dimensionless]
		Zsn = @ccall $interpolate1D(
			Z::Cdouble,
			Zsn_table::Cstring,
			R_ZSN_NROWS::Cint,
			R_ZSN_NCOLS::Cint,
		)::Cdouble

		parameters = [ρ_cell, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h]
		
		@ccall $integrate_ode(
		    fractions::Ptr{Cdouble},
			parameters::Ptr{Cdouble},
		    it::Cdouble,
		)::Cvoid

		return fractions

	end

	return integration

end;

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
ChaosTools = "608a59af-f2a3-5ad4-90b4-758bdf3122a7"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
Libdl = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
PlutoLinks = "0ff47ea0-7a50-410d-8455-4348d5de0420"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
SciMLLogging = "a6db7da4-7206-11f0-1eab-35f2a5dbe1d1"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"

[compat]
CSV = "~0.10.16"
ChaosTools = "~3.5.4"
DataFrames = "~1.8.2"
DifferentialEquations = "~8.0.0"
Interpolations = "~0.16.2"
Measurements = "~2.14.1"
NaNMath = "~1.1.4"
PlutoLinks = "~0.1.8"
PlutoUI = "~0.7.83"
QuadGK = "~2.11.3"
SciMLLogging = "~2.0.0"
SpecialFunctions = "~2.8.0"
Symbolics = "~7.26.0"
TikzPictures = "~3.5.1"
Unitful = "~1.28.0"
UnitfulAstro = "~1.2.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.6"
manifest_format = "2.0"
project_hash = "9d3d3f94d1124376e699c8a76d2639d341d1604e"

[[deps.ADTypes]]
git-tree-sha1 = "bbc22a9a08a0ef6460041086d8a7b27940ed4ffd"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.22.0"
weakdeps = ["ChainRulesCore", "ConstructionBase", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

[[deps.AMD]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "45a1272e3f809d36431e57ab22703c6896b8908f"
uuid = "14f7f29c-3bd6-536c-9a0b-7339e30b5a3e"
version = "0.5.3"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "2eeb2c9bef11013efc6f8f97f32ee59b146b09fb"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.44"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7715e5b2b186c4d9b664d299d2c9e48b9a778c88"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.6.1"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "3d0cabd25fab32390e3bcb82cd67e700aebd9816"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.25.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceAMDGPUExt = "AMDGPU"
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.BracketingNonlinearSolve]]
deps = ["CommonSolve", "ConcreteStructs", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "7ad7171d693ae5552ac43862e7f6b61df4471c2b"
uuid = "70df07ce-3d50-431d-a3e7-ca6ddb60ac1e"
version = "1.12.1"
weakdeps = ["ChainRulesCore", "ForwardDiff"]

    [deps.BracketingNonlinearSolve.extensions]
    BracketingNonlinearSolveChainRulesCoreExt = ["ChainRulesCore", "ForwardDiff"]
    BracketingNonlinearSolveForwardDiffExt = "ForwardDiff"

[[deps.BranchAndPrune]]
deps = ["AbstractTrees"]
git-tree-sha1 = "8bea09a63c60b6b5200e613379947e7c8527992b"
uuid = "d3bc4f2e-91e6-11e9-365e-cd067da536ce"
version = "0.2.2"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "8d8e0b0f350b8e1c91420b5e64e5de774c2f0f4d"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.16"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2ac646d71d0d24b44f3f8c84da8c9f4d70fb67df"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.4+0"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9cb23bbb1127eefb022b022481466c0f1127d430"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "12177ad6b3cad7fd50c8b3825ce24a99ad61c18f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChaosTools]]
deps = ["Combinatorics", "DSP", "DataStructures", "Distances", "Distributions", "DynamicalSystemsBase", "IntervalRootFinding", "LinearAlgebra", "LombScargle", "Neighborhood", "Optim", "ProgressMeter", "Random", "Reexport", "Roots", "SpecialFunctions", "Statistics", "StatsBase"]
git-tree-sha1 = "12f9540406976e23eb2d7396d5bbe6e0e7be7874"
uuid = "608a59af-f2a3-5ad4-90b4-758bdf3122a7"
version = "3.5.4"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "REPL", "UUIDs"]
git-tree-sha1 = "cfb7a2e89e245a9d5016b70323db412b3a7438d5"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "3.0.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "dd91a10d8b8ae06e15706158eaf1a3e87e97b5f5"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.7"
weakdeps = ["EnzymeCore"]

    [deps.CommonSolve.extensions]
    CommonSolveEnzymeCoreExt = "EnzymeCore"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.Compiler]]
git-tree-sha1 = "382d79bfe72a406294faca39ef0c3cef6e6ce1f1"
uuid = "807dbc54-b67e-4c79-8afb-eafe4df6f2e1"
version = "0.1.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "ed1da4eac5ba9b3f6d061c90f3ca6ba190dd6595"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.4"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.CoreMath]]
deps = ["CoreMath_jll"]
git-tree-sha1 = "8c0480f92b1b1796239156a1b9b1bfb1b39499b4"
uuid = "b7a15901-be09-4a0e-87d2-2e66b0e09b5a"
version = "0.1.0"

[[deps.CoreMath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a692a4c1dc59a4b8bc0b6403876eb3250fde2bc3"
uuid = "a38c48d9-6df1-5ac9-9223-b6ada3b5572b"
version = "0.1.0+0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DSP]]
deps = ["Bessels", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "d335b2929e1b6067951a1250df247cc5fab7d40e"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.8.5"
weakdeps = ["OffsetArrays"]

    [deps.DSP.extensions]
    OffsetArraysExt = "OffsetArrays"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "5fab31e2e01e70ad66e3e24c968c264d1cf166d6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.2"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "6fb53a69613a0b2b68a0d12671717d307ab8b24e"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.5"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "BracketingNonlinearSolve", "ConcreteStructs", "DocStringExtensions", "FastBroadcast", "FastClosures", "FastPower", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLLogging", "SciMLOperators", "SciMLStructures", "Setfield", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "c3ab0a3d326ab31f1318d65776fcd9e56bc8916d"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "7.5.5"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDynamicQuantitiesExt = "DynamicQuantities"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseFlexUnitsExt = "FlexUnits"
    DiffEqBaseForwardDiffExt = ["ForwardDiff"]
    DiffEqBaseGTPSAExt = "GTPSA"
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseMooncakeExt = "Mooncake"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseSparseArraysExt = "SparseArrays"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FlexUnits = "76e01b6b-c995-4ce6-8559-91e72a3d4e95"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffEqCallbacks]]
deps = ["ConcreteStructs", "DataStructures", "DiffEqBase", "DifferentiationInterface", "LinearAlgebra", "Markdown", "PrecompileTools", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "450217b2cc2fe3396de5452ca315952573f0d452"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "4.18.0"

    [deps.DiffEqCallbacks.extensions]
    DiffEqCallbacksFunctorsExt = "Functors"

    [deps.DiffEqCallbacks.weakdeps]
    Functors = "d9f16b24-f501-4c13-a1f2-28368ffc5196"

[[deps.DiffEqNoiseProcess]]
deps = ["CommonSolve", "DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "PoissonRandom", "QuadGK", "Random", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "bd078194af24a5674a6bad822669f9489d6b8c10"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.32.0"

    [deps.DiffEqNoiseProcess.extensions]
    DiffEqNoiseProcessOptimExt = "Optim"
    DiffEqNoiseProcessReverseDiffExt = "ReverseDiff"

    [deps.DiffEqNoiseProcess.weakdeps]
    Optim = "429524aa-4258-5aef-a3af-852621145aeb"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "79a2aca180a85c690c58a020d47b426954b590f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.16.0"

[[deps.DifferentialEquations]]
deps = ["OrdinaryDiffEq", "Reexport", "SciMLBase"]
git-tree-sha1 = "076a5e999d38ffbba860d4c1ba49828855b65c14"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "8.0.0"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "2147a95a217cc8a78ec96ee03581adf129468e49"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.18"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = ["GPUArraysCore", "Adapt"]
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceHyperHessiansExt = "HyperHessians"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    HyperHessians = "06b494a0-c8e0-40cc-ad32-d99506a00a6c"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "96f76dcd6cc75cf8eb49109123868499d413f526"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.126"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.DomainSets]]
deps = ["CompositeTypes", "FunctionMaps", "IntervalSets", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "1af39efbaf76fb648432b5efaac0d73af6760407"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.8.0"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"
    DomainSetsRandomExt = "Random"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.DynamicPolynomials]]
deps = ["LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "StarAlgebras", "Test"]
git-tree-sha1 = "5bfabc3827dfdd164359bad6800c115a81280c00"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.6"

[[deps.DynamicalSystemsBase]]
deps = ["ForwardDiff", "LinearAlgebra", "OrdinaryDiffEqTsit5", "Random", "Reexport", "Roots", "SciMLBase", "SparseArrays", "StateSpaceSets", "Statistics", "StochasticDiffEqHighOrder", "SymbolicIndexingInterface"]
git-tree-sha1 = "72103afddf62ddc10aa6f9b4e553e338d45a5dff"
uuid = "6e36e845-645a-534a-86f2-f5d4aa5a06b4"
version = "3.18.1"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.EnzymeCore]]
git-tree-sha1 = "c6ee69ee502060982d12dbaaf3d8fcb4e835a0d1"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.20"
weakdeps = ["Adapt", "ChainRulesCore"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"
    EnzymeCoreChainRulesCoreExt = "ChainRulesCore"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c307cd83373868391f3ac30b41530bc5d5d05d08"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.8.1+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6866aec60ef98e3164cd8d6855225684207e9dff"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.12+0"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra"]
git-tree-sha1 = "7feeed2c9a7fa272189a5561bebf0c4ccaedb6ec"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "1.3.2"

    [deps.FastBroadcast.extensions]
    FastBroadcastPolyesterExt = "Polyester"
    FastBroadcastStaticExt = "Static"

    [deps.FastBroadcast.weakdeps]
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"
    Static = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastPower]]
git-tree-sha1 = "862831f78c7a48681a074ecc9aac09f2de563f71"
uuid = "a4df4552-cc26-4903-aec0-212e50a0e84b"
version = "1.3.1"

    [deps.FastPower.extensions]
    FastPowerEnzymeExt = "Enzyme"
    FastPowerForwardDiffExt = "ForwardDiff"
    FastPowerMeasurementsExt = "Measurements"
    FastPowerMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    FastPowerMooncakeExt = "Mooncake"
    FastPowerReverseDiffExt = "ReverseDiff"
    FastPowerTrackerExt = "Tracker"

    [deps.FastPower.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2f979084d1e13948a3352cf64a25df6bd3b4dca3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.16.0"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "f7017a4f337f8df189fcce98e32b67a1298a2115"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.31.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "2c5d0b0e12088cde2cf84afb2784415b1ea3dfee"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.4.1"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

[[deps.FunctionMaps]]
deps = ["CompositeTypes", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "31bd99a57edf98990d1c21486032963955450e8d"
uuid = "a85aefff-f8ca-4649-a888-c8e5398bc76c"
version = "0.1.2"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers", "PrecompileTools", "TruncatedStacktraces"]
git-tree-sha1 = "b28fca87e487d18ba462317e4a95d7253ae51929"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "1.9.1"

    [deps.FunctionWrappersWrappers.extensions]
    FunctionWrappersWrappersEnzymeExt = ["Enzyme", "EnzymeCore"]
    FunctionWrappersWrappersMooncakeExt = "Mooncake"

    [deps.FunctionWrappersWrappers.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "b0036b392358c80d2d2124746c2bf3d48d457938"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.4+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Inflate", "LinearAlgebra", "Random", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "7eb45fe833a5b7c51cf6d89c5a841d5967e44be3"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.14.0"
weakdeps = ["Distributed", "SharedArrays"]

    [deps.Graphs.extensions]
    GraphsSharedArraysExt = "SharedArrays"

[[deps.HarfBuzz_ICU_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "HarfBuzz_jll", "ICU_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "6ccbc4fdf65c8197738c2d68cc55b74b19c97ac2"
uuid = "655565e8-fb53-5cb3-b0cd-aec1ca0647ea"
version = "2.8.1+0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.ICU_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "20b6765a3016e1fca0c9c93c80d50061b94218b7"
uuid = "a51ab1cf-af8e-5615-a023-bc2c838bba6b"
version = "69.1.0+0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InlineStrings]]
git-tree-sha1 = "8f3d257792a522b4601c24a577954b0a8cd7334d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.5"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"
weakdeps = ["ForwardDiff", "Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "CoreMath", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "921d7e91687e15a2c7c269c226960491fc041832"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.9"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticIrrationalConstantsExt = "IrrationalConstants"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    IrrationalConstants = "92d709cd-6900-40b7-9082-c6be49f344b6"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalRootFinding]]
deps = ["BranchAndPrune", "ForwardDiff", "IntervalArithmetic", "LinearAlgebra", "Reexport", "StaticArrays"]
git-tree-sha1 = "e0bbe14c3e38ddf829f130e698dd53f5d7400cce"
uuid = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
version = "0.6.3"

[[deps.IntervalSets]]
git-tree-sha1 = "79d6bd28c8d9bccc2229784f1bd637689b256377"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.14"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c0c9b76f3520863909825cbecdef58cd63de705a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.5+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "58927c485919bf17ea308d9d82156de1adf4b006"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.10.12"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "PoissonRandom", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "SymbolicIndexingInterface"]
git-tree-sha1 = "0026bbad898bd6b9ce31772be345d4875cfc0e5a"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.29.0"

    [deps.JumpProcesses.extensions]
    JumpProcessesKernelAbstractionsExt = ["Adapt", "KernelAbstractions"]
    JumpProcessesOrdinaryDiffEqCoreExt = "OrdinaryDiffEqCore"

    [deps.JumpProcesses.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    OrdinaryDiffEqCore = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "c4d19f51afc7ba2afbe32031b8f2d21b11c9e26e"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.10.6"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "17b94ecafcfa45e8360a4fc9ca6b583b049e4e37"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.1.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "27ecae93dd25ee0909666e6835051dd684cc035e"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+2"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc3ad4faf30015a3e8094c9b5b7f19e85bdf2386"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d620582b1f0cbe2c72dd1d5bd195a9ce73370ab1"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.0+0"

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "PrecompileTools", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "fd58a77c92e7c8f1db25c9839127d52943a49349"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.9"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Printf"]
git-tree-sha1 = "cef1ba655e8c1f65af9d96c4fffe18bf1a3a3291"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.7.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LinearSolve]]
deps = ["AMD", "ArrayInterface", "ConcreteStructs", "DocStringExtensions", "EnumX", "GPUArraysCore", "InteractiveUtils", "Krylov", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "OpenBLAS_jll", "PrecompileTools", "Preferences", "PureKLU", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLLogging", "SciMLOperators", "Setfield", "SparseColumnPivotedQR", "StaticArraysCore"]
git-tree-sha1 = "70de717540ffc613424866635b99d0c187cf7b59"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "3.85.1"

    [deps.LinearSolve.extensions]
    LinearSolveAMDGPUExt = "AMDGPU"
    LinearSolveAlgebraicMultigridExt = "AlgebraicMultigrid"
    LinearSolveBLISExt = ["blis_jll", "LAPACK_jll"]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = ["cuSOLVER"]
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveCUSOLVERRFExt = ["CUSOLVERRF", "SparseArrays"]
    LinearSolveChainRulesCoreExt = "ChainRulesCore"
    LinearSolveCliqueTreesExt = ["CliqueTrees", "SparseArrays"]
    LinearSolveElementalExt = "Elemental"
    LinearSolveEnzymeExt = ["EnzymeCore", "SparseArrays"]
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveFastLapackInterfaceExt = "FastLapackInterface"
    LinearSolveForwardDiffExt = "ForwardDiff"
    LinearSolveGinkgoExt = ["Ginkgo", "SparseArrays"]
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMUMPSExt = ["MUMPS", "SparseArrays"]
    LinearSolveMetalExt = "Metal"
    LinearSolveMooncakeExt = "Mooncake"
    LinearSolvePETScExt = ["PETSc", "SparseArrays", "SparseMatricesCSR"]
    LinearSolvePETScMPIExt = ["PETSc", "PartitionedArrays", "SparseArrays", "SparseMatricesCSR"]
    LinearSolveParUExt = ["ParU_jll", "SparseArrays"]
    LinearSolvePardisoExt = ["Pardiso", "SparseArrays"]
    LinearSolvePartitionedSolversExt = ["PartitionedArrays", "PartitionedSolvers"]
    LinearSolveRecursiveFactorizationExt = "RecursiveFactorization"
    LinearSolveSTRUMPACKExt = ["SparseArrays", "STRUMPACK_jll"]
    LinearSolveSparseArraysExt = "SparseArrays"
    LinearSolveSparspakExt = ["SparseArrays", "Sparspak"]

    [deps.LinearSolve.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    AlgebraicMultigrid = "2169fc97-5a83-5252-b627-83903c6c433c"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    CUSOLVERRF = "a8cc9031-bad2-4722-94f5-40deabb4245c"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Elemental = "902c3f28-d1ec-5e7e-8399-a24c3845ee38"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    FastLapackInterface = "29a986be-02c6-4525-aec4-84b980013641"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Ginkgo = "4c8bd3c9-ead9-4b5e-a625-08f1338ba0ec"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    LAPACK_jll = "51474c39-65e3-53ba-86ba-03b1b862ec14"
    MUMPS = "55d2b088-9f4e-11e9-26c0-150b02ea6a46"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PETSc = "ace2c81b-2b5f-4b1e-a30d-d662738edfe0"
    ParU_jll = "9e0b026c-e8ce-559c-a2c4-6a3d5c955bc9"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    PartitionedArrays = "5a9dfac6-5c52-46f7-8278-5e2210713be9"
    PartitionedSolvers = "11b65f7f-80ac-401b-9ef2-3db765482d62"
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
    STRUMPACK_jll = "86fbd0b9-476f-557c-b766-62c724b42d8c"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseMatricesCSR = "a0a7dd2c-ebf4-11e9-1f05-cf50bc540ca1"
    Sparspak = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
    blis_jll = "6136c539-28a5-5bf0-87cc-b183200dce32"
    cuSOLVER = "887afef0-6a32-4de5-add4-7827692ba8fc"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll"]
git-tree-sha1 = "38928f7999753af13d4e13966ae15958ff3a917a"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.19.1+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "bba2d9aa057d8f126415de240573e86a8f39d2a1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "1.0.1"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.LombScargle]]
deps = ["FFTW", "LinearAlgebra", "Measurements", "Random", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "d64a0ce7539181136a85fd8fe4f42626387f0f26"
uuid = "fc60dff9-86e7-5f2f-a8a0-edeadbb75bd9"
version = "1.0.3"

[[deps.LoweredCodeUtils]]
deps = ["CodeTracking", "Compiler", "JuliaInterpreter"]
git-tree-sha1 = "0aad96d7b987a5600e260eec50147b254d5ff7e6"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.6.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "54e2fdc38130c05b42be423e90da3bade29b74bd"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.4"
weakdeps = ["SparseArrays"]

    [deps.MaybeInplace.extensions]
    MaybeInplaceSparseArraysExt = "SparseArrays"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf"]
git-tree-sha1 = "cb47f69a1cab9dcec7ff4a5d6e163410d6905866"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.14.1"

    [deps.Measurements.extensions]
    MeasurementsBaseTypeExt = "BaseType"
    MeasurementsJunoExt = "Juno"
    MeasurementsMakieExt = "Makie"
    MeasurementsRecipesBaseExt = "RecipesBase"
    MeasurementsSpecialFunctionsExt = "SpecialFunctions"
    MeasurementsUnitfulExt = "Unitful"

    [deps.Measurements.weakdeps]
    BaseType = "7fbed51b-1ef5-4d67-9085-a4a9b26f478c"
    Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "999dfa2b4f8334c1c23bd307c2b0cb6f97c54591"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.8"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics", "StarAlgebras"]
git-tree-sha1 = "4838893d9b035c2f6967c0d533350e1755b58a70"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.19"
weakdeps = ["ChainRulesCore"]

    [deps.MultivariatePolynomials.extensions]
    MultivariatePolynomialsChainRulesCoreExt = "ChainRulesCore"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "dc5b2c4c111c46bc79ac4405eeb563523b39c004"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.8.0"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "FiniteDiff", "LinearAlgebra"]
git-tree-sha1 = "b3f76b463c7998473062992b246045e6961a074e"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "8.0.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "dbd2e8cd2c1c27f0b584f6661b4309609c5a685e"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.4"

[[deps.NearestNeighbors]]
deps = ["AbstractTrees", "Distances", "StaticArrays"]
git-tree-sha1 = "e2c3bba08dd6dedfe17a17889131b885b8c082f0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.27"

[[deps.Neighborhood]]
deps = ["Distances", "NearestNeighbors", "Random", "Test"]
git-tree-sha1 = "3efba0551882b18a1e97d27fafe8d01fb70b5b1b"
uuid = "645ca80c-8b79-4109-87ea-e1f58159d116"
version = "0.2.5"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "NonlinearSolveBase", "NonlinearSolveFirstOrder", "NonlinearSolveQuasiNewton", "NonlinearSolveSpectralMethods", "PrecompileTools", "Preferences", "Reexport", "SciMLBase", "Setfield", "SimpleNonlinearSolve", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "a6c5719bbb42985c72f4cacbfa49e86bab850d66"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "4.19.1"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = ["NLsolve", "LineSearches"]
    NonlinearSolvePETScExt = ["PETSc", "MPI", "SparseArrays"]
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"
    NonlinearSolveSundialsExt = "Sundials"

    [deps.NonlinearSolve.weakdeps]
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    LineSearches = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    PETSc = "ace2c81b-2b5f-4b1e-a30d-d662738edfe0"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"
    Sundials = "c3572dad-4567-51f8-b174-8c6c989267f4"

[[deps.NonlinearSolveBase]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "Compat", "ConcreteStructs", "DifferentiationInterface", "EnzymeCore", "FastClosures", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "LogExpFunctions", "Markdown", "MaybeInplace", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "RecursiveArrayTools", "SciMLBase", "SciMLJacobianOperators", "SciMLLogging", "SciMLOperators", "SciMLStructures", "Setfield", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "8676d489aac4cae9977d4520db3aa4fbe57d3bf3"
uuid = "be0214bd-f91f-a760-ac4e-3421ce2b2da0"
version = "2.31.0"

    [deps.NonlinearSolveBase.extensions]
    NonlinearSolveBaseBandedMatricesExt = "BandedMatrices"
    NonlinearSolveBaseChainRulesCoreExt = "ChainRulesCore"
    NonlinearSolveBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    NonlinearSolveBaseForwardDiffExt = "ForwardDiff"
    NonlinearSolveBaseLineSearchExt = "LineSearch"
    NonlinearSolveBaseLinearSolveExt = "LinearSolve"
    NonlinearSolveBaseMooncakeExt = "Mooncake"
    NonlinearSolveBaseReverseDiffExt = "ReverseDiff"
    NonlinearSolveBaseSparseArraysExt = "SparseArrays"
    NonlinearSolveBaseSparseMatrixColoringsExt = "SparseMatrixColorings"
    NonlinearSolveBaseTrackerExt = "Tracker"

    [deps.NonlinearSolveBase.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    LineSearch = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
    LinearSolve = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.NonlinearSolveFirstOrder]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConcreteStructs", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLJacobianOperators", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "ce68820a4f421fb5bee7ec4dcf875aff33886bfb"
uuid = "5959db7a-ea39-4486-b5fe-2dd0bf03d60d"
version = "2.1.1"

[[deps.NonlinearSolveQuasiNewton]]
deps = ["ArrayInterface", "CommonSolve", "ConcreteStructs", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLOperators", "StaticArraysCore"]
git-tree-sha1 = "538432ca1aea8bf63db02929bf870501f8a7c64c"
uuid = "9a2c21bd-3a47-402d-9113-8faf9a0ee114"
version = "1.13.1"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveQuasiNewton.extensions]
    NonlinearSolveQuasiNewtonForwardDiffExt = "ForwardDiff"

[[deps.NonlinearSolveSpectralMethods]]
deps = ["CommonSolve", "ConcreteStructs", "LineSearch", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "a3781e12becdf0ce5520bd97ec617e879bf4e9f2"
uuid = "26075421-4e9a-44e1-8bd1-420ed7ad02b2"
version = "1.7.1"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveSpectralMethods.extensions]
    NonlinearSolveSpectralMethodsForwardDiffExt = "ForwardDiff"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3287ec88df50429a934ebc6cf14606215e27b987"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.33+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "libpng_jll"]
git-tree-sha1 = "215a6666fee6d6b3a6e75f2cc22cb767e2dd393a"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.5.5+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["ADTypes", "EnumX", "FillArrays", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "PositiveFactorizations", "Printf", "SparseArrays", "Statistics"]
git-tree-sha1 = "7ecee58b8bd88cc8bbdb90d056d1c6546322eebd"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "2.1.0"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.OrderedCollections]]
git-tree-sha1 = "94ba93778373a53bfd5a0caaf7d809c445292ff4"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.2"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "CommonSolve", "DocStringExtensions", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "SciMLBase"]
git-tree-sha1 = "4c411ca0fa291015bc4c23058799a7280c4223e4"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "7.0.0"

[[deps.OrdinaryDiffEqBDF]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqSDIRK", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "d53667de07ac995298bc7549efcb9dc94d53e66c"
uuid = "6ad6398a-0878-4a85-9266-38940aa047c8"
version = "2.2.0"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "ConcreteStructs", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "FastPower", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "PrecompileTools", "Preferences", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLLogging", "SciMLOperators", "SciMLStructures", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "16bb44e25e6a2be15a47dba831355f6641797ea2"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "4.3.0"

    [deps.OrdinaryDiffEqCore.extensions]
    OrdinaryDiffEqCoreMooncakeExt = "Mooncake"
    OrdinaryDiffEqCorePolyesterExt = "Polyester"
    OrdinaryDiffEqCoreSparseArraysExt = "SparseArrays"

    [deps.OrdinaryDiffEqCore.weakdeps]
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.OrdinaryDiffEqDefault]]
deps = ["ADTypes", "DiffEqBase", "EnumX", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "PrecompileTools", "Preferences", "Reexport", "SciMLBase"]
git-tree-sha1 = "60c06fe39b257d1b222f188902f9648681386652"
uuid = "50262376-6c5a-4cf5-baba-aaf4f84d72d7"
version = "2.2.0"

[[deps.OrdinaryDiffEqDifferentiation]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqCore", "SciMLBase", "SciMLOperators", "SparseMatrixColorings", "StaticArraysCore"]
git-tree-sha1 = "61a383d64ca9632c5cd300ac7a01042dc1f1c7c7"
uuid = "4302a76b-040a-498a-8c04-15b101fed76b"
version = "3.2.0"
weakdeps = ["SparseArrays"]

    [deps.OrdinaryDiffEqDifferentiation.extensions]
    OrdinaryDiffEqDifferentiationSparseArraysExt = "SparseArrays"

[[deps.OrdinaryDiffEqNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "PreallocationTools", "RecursiveArrayTools", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SparseArrays", "StaticArraysCore"]
git-tree-sha1 = "c3e92bdb7bf16385e6f7b505f78601a2708eff28"
uuid = "127b3ac7-2247-4354-8eb6-78cf4e7c58e8"
version = "2.0.0"

[[deps.OrdinaryDiffEqRosenbrock]]
deps = ["ADTypes", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqRosenbrockTableaus", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "00a6f06d131c0778dd1da79cf13ef330021c709a"
uuid = "43230ef6-c299-4910-a778-202eb28ce4ce"
version = "2.3.0"

[[deps.OrdinaryDiffEqRosenbrockTableaus]]
git-tree-sha1 = "968933504b01dc4a8a95ce5fcdd47f03d831883c"
uuid = "b4bd8bb3-f80f-41d2-9b21-73a655b304b9"
version = "2.1.0"

[[deps.OrdinaryDiffEqSDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "ed8aae3b58145b6ea46ac51b792b9348d291e4e3"
uuid = "2d112036-d095-4a1e-ab9a-08536f3ecdbf"
version = "2.7.0"

[[deps.OrdinaryDiffEqTsit5]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "21f58fc5ef0103d7708fc3a2b7b24ef16353d883"
uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a"
version = "2.0.1"

[[deps.OrdinaryDiffEqVerner]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "6061884927ecbd67f44d1ea655379ada4e3fdac9"
uuid = "79d7bb75-1356-48c1-b8c0-6832512096c2"
version = "2.1.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e4cff168707d441cd6bf3ff7e4832bdf34278e4a"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.37"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "468dbe2b510c876dc091b2c74ed52c7c34f48b9b"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.5"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "844a829c8dc9fd0fe62eced22bc2d0dfd66a3f51"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.1.0"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "aea4eede5ab3ee188906d0cf3bbfa36eb543dccc"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.8"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.PoissonRandom]]
deps = ["LogExpFunctions", "PrecompileTools", "Random"]
git-tree-sha1 = "ef5d9046fd90aeb2e0ae271a7550a7d969737157"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.9"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "Setfield", "SparseArrays"]
git-tree-sha1 = "2d99b4c8a7845ab1342921733fa29366dae28b24"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.1"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"
    PolynomialsRecipesBaseExt = "RecipesBase"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.Poppler_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "LibCURL_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "libpng_jll"]
git-tree-sha1 = "7dbfb7f61c3aa5def7b7dad3fa344c1c2858a83b"
uuid = "9c32591e-4766-534b-9725-b71a8799265b"
version = "24.6.0+0"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "PrecompileTools"]
git-tree-sha1 = "e16b73bf892c55d16d53c9c0dbd0fb31cb7e25da"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "1.2.0"

    [deps.PreallocationTools.extensions]
    PreallocationToolsForwardDiffExt = "ForwardDiff"
    PreallocationToolsReverseDiffExt = "ReverseDiff"
    PreallocationToolsSparseConnectivityTracerExt = "SparseConnectivityTracer"

    [deps.PreallocationTools.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "624de6279ab7d94fc9f672f0068107eb6619732c"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.3.2"

    [deps.PrettyTables.extensions]
    PrettyTablesTypstryExt = "Typstry"

    [deps.PrettyTables.weakdeps]
    Typstry = "f0ed7684-a786-439e-b1e3-3b82803b501e"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.PureKLU]]
deps = ["LinearAlgebra", "MuladdMacro", "PrecompileTools", "SparseArrays"]
git-tree-sha1 = "02d7809d74834e7b04ac8858704c4efc34283f8a"
uuid = "0c0d3e7f-3a8b-4f7e-b6f1-9a4d2e7c1f01"
version = "1.0.1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "5e8e8b0ab68215d7a2b14b9921a946fee794749e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.3"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.ReadOnlyArrays]]
git-tree-sha1 = "e6f7ddf48cf141cb312b078ca21cb2d29d0dc11d"
uuid = "988b38a3-91fc-5605-94a2-ee2116b3bd83"
version = "0.2.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "LinearAlgebra", "PrecompileTools", "RecipesBase", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "b2114acf9f93c32953cf2bbc70cb2a9029a986f3"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "4.3.1"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsCUDAExt = "CUDA"
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsFastBroadcastPolyesterExt = ["FastBroadcast", "Polyester"]
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsKernelAbstractionsExt = "KernelAbstractions"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsMooncakeExt = "Mooncake"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStatisticsExt = "Statistics"
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTablesExt = ["Tables"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "31c086583c92ab32d82ebef0d09fbcd6dd2c54a7"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.2.0"

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "InteractiveUtils", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Preferences", "REPL", "UUIDs"]
git-tree-sha1 = "65569cca6282716a14177af1358db91ef79300f0"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.15.0"
weakdeps = ["Distributed"]

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "3eff988b9bd09543783e2e051b0a1eef11f65c2d"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.3.0"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "28154d426e557495aa75097861b18c11f2541ded"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.19"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "Random", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLLogging", "SciMLOperators", "SciMLPublic", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "6f35020fba2d580e324156cbaff4c578afa48940"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "3.20.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseDifferentiationInterfaceExt = "DifferentiationInterface"
    SciMLBaseDistributionsExt = "Distributions"
    SciMLBaseEnzymeExt = "Enzyme"
    SciMLBaseForwardDiffExt = "ForwardDiff"
    SciMLBaseMakieExt = "Makie"
    SciMLBaseMeasurementsExt = "Measurements"
    SciMLBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    SciMLBaseMooncakeExt = "Mooncake"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseReverseDiffExt = "ReverseDiff"
    SciMLBaseStaticArraysExt = "StaticArrays"
    SciMLBaseTrackerExt = "Tracker"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DifferentiationInterface = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLJacobianOperators]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DifferentiationInterface", "FastClosures", "LinearAlgebra", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "7156a5b51cba1bea33a82a036939ead4131f92bc"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.13"

[[deps.SciMLLogging]]
deps = ["Logging", "LoggingExtras", "Preferences"]
git-tree-sha1 = "3f98a53703f925cbd5aac5da4924f82ca34d09ab"
uuid = "a6db7da4-7206-11f0-1eab-35f2a5dbe1d1"
version = "2.0.0"

    [deps.SciMLLogging.extensions]
    SciMLLoggingTracyExt = "Tracy"

    [deps.SciMLLogging.weakdeps]
    Tracy = "e689c965-62c8-4b79-b2c5-8359227902fd"

[[deps.SciMLOperators]]
deps = ["Accessors", "Adapt", "ArrayInterface", "DocStringExtensions", "LinearAlgebra"]
git-tree-sha1 = "3b204078e8574b9de19cac90e0c87c811a87deac"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.22.0"

    [deps.SciMLOperators.extensions]
    SciMLOperatorsLoopVectorizationExt = "LoopVectorization"
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

    [deps.SciMLOperators.weakdeps]
    LoopVectorization = "bdcacae8-1622-11e9-2a5c-532679323890"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"

[[deps.SciMLPublic]]
git-tree-sha1 = "0ba076dbdce87ba230fff48ca9bca62e1f345c9b"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.1"

[[deps.SciMLStructures]]
deps = ["ArrayInterface", "PrecompileTools"]
git-tree-sha1 = "607f6867d0b0553e98fc7f725c9f9f13b4d01a32"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.10.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "084c47c7c5ce5cfecefa0a98dff69eb3646b5a80"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.10"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "08d8bafd57b7ea16dc1f69a8bd57db23a2ceb686"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "2.12.0"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveTrackerExt = "Tracker"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "7ddb0b49c109481b046972c0e4ab02b2127d6a75"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SparseColumnPivotedQR]]
deps = ["LinearAlgebra", "PrecompileTools", "SparseArrays"]
git-tree-sha1 = "aa3872796237441bc82f51b6f094aaf97386a717"
uuid = "a57abbd0-fea5-4d57-96be-5e525945e8e4"
version = "2.1.1"
weakdeps = ["AMD"]

    [deps.SparseColumnPivotedQR.extensions]
    SparseColumnPivotedQRAMDExt = "AMD"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "Random", "SparseArrays"]
git-tree-sha1 = "f63d76c7b7c329cf11badd564fd8ba877b09c3fe"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.27"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsCUDAExt = ["CUDA", "cuSPARSE"]
    SparseMatrixColoringsCliqueTreesExt = "CliqueTrees"
    SparseMatrixColoringsColorsExt = "Colors"
    SparseMatrixColoringsGPUArraysExt = "GPUArrays"
    SparseMatrixColoringsJuMPExt = ["JuMP", "MathOptInterface"]

    [deps.SparseMatrixColorings.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
    GPUArrays = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
    JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
    cuSPARSE = "b26da814-b3bc-49ef-b0ee-c816305aa060"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "6547cbdd8ce32efba0d21c5a40fa96d1a3548f9f"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.8.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StarAlgebras]]
deps = ["LinearAlgebra", "MutableArithmetics", "SparseArrays"]
git-tree-sha1 = "235b1f9d287bbf34083b3d0829343a7942c0ad1c"
uuid = "0c0c59c1-dc5f-42e9-9a8b-b5dc384a6cd1"
version = "0.3.0"

[[deps.StateSpaceSets]]
deps = ["Distances", "LinearAlgebra", "Neighborhood", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "b6cbad6435392da931d2d99fa54b96b09400d2eb"
uuid = "40b095a5-5852-4c12-98c7-d43bf788e795"
version = "2.5.5"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.18"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "c6f18e5a52a176a383f6f6c635e0f81feed1d6d4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.11"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "3f4e1d24289cd974e089c617b1472311a2b1feab"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "2.1.0"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StochasticDiffEqCore]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "FastPower", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LinearAlgebra", "Logging", "MuladdMacro", "OrdinaryDiffEqCore", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLLogging", "SciMLOperators", "SimpleNonlinearSolve", "SparseArrays", "StaticArrays", "StochasticDiffEqLevyArea"]
git-tree-sha1 = "feeb9e364c903555ee3613a2e2bc1a29f6343656"
uuid = "19c5a474-6cd1-4a5f-be79-46dc34e54d7f"
version = "2.0.1"

[[deps.StochasticDiffEqHighOrder]]
deps = ["DiffEqBase", "DiffEqNoiseProcess", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "StochasticDiffEqCore"]
git-tree-sha1 = "43dc5e03e7e5b7457a20f861c9f4b1d4cd6600be"
uuid = "0520c28c-50fd-4d16-9c96-902fc80b3bab"
version = "2.1.0"

[[deps.StochasticDiffEqLevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "083805654fde34d3d32d4dd91bcc20f4c0330db3"
uuid = "90dbc90e-856a-4131-af2c-e8c5aa0f35b5"
version = "2.0.0"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "d05693d339e37d6ab134c5ab53c29fce5ee5d7d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.4"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "64b15330b9e3c91a86bcac92f369c58e382981c6"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.48"
weakdeps = ["PrettyTables"]

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils", "TermInterface"]
git-tree-sha1 = "5085671d2cba1eb02136a3d6661c583e801984c1"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "1.1.0"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "EnumX", "ExproniconLite", "Graphs", "LinearAlgebra", "MacroTools", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "ReadOnlyArrays", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "WeakCacheSets"]
git-tree-sha1 = "db2dd7feba95f6d94d52e7a63185449d57f808c5"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "4.35.0"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsChainRulesCoreExt = "ChainRulesCore"
    SymbolicUtilsDistributionsExt = "Distributions"
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "AbstractPlutoDingetjes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "Preferences", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "a73710e38577a31b32e4d204d362f768bc5badff"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "7.26.0"

    [deps.Symbolics.extensions]
    SymbolicsD3TreesExt = "D3Trees"
    SymbolicsDistributionsExt = "Distributions"
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsHypergeometricFunctionsExt = "HypergeometricFunctions"
    SymbolicsLatexifyExt = ["Latexify", "LaTeXStrings"]
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"
    SymbolicsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Symbolics.weakdeps]
    D3Trees = "e3df1716-f71e-5df9-9e2d-98e193103c45"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    HypergeometricFunctions = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TikzPictures]]
deps = ["LaTeXStrings", "Poppler_jll", "tectonic_jll"]
git-tree-sha1 = "875854f63fbe215b554390efd249bfbef1418c31"
uuid = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
version = "3.5.1"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "57e1b2c9de4bd6f40ecb9de4ac1797b81970d008"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.28.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.UnitfulAngles]]
deps = ["Dates", "Unitful"]
git-tree-sha1 = "79875b1f2e4bf918f0702a5980816955066d9ae2"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.7.2"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "fbe44a0ade62ae5ed0240ad314dfdd5482b90b40"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.2.2"

[[deps.WeakCacheSets]]
git-tree-sha1 = "386050ae4353310d8ff9c228f83b1affca2f7f38"
uuid = "d30d5f5c-d141-4870-aa07-aabb0f5fe7d5"
version = "0.1.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "0716e01c3b40413de5dedbc9c5c69f27cddfddfc"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "80d3930c6347cfce7ccf96bd3bafdf079d9c0390"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.9+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e51150d5ab85cee6fc36726850f0e627ad2e4aba"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.58+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "da8c1f6eee04831f14edcfa5dae611d309807e57"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.3.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.tectonic_jll]]
deps = ["Artifacts", "Fontconfig_jll", "FreeType2_jll", "Graphite2_jll", "HarfBuzz_ICU_jll", "HarfBuzz_jll", "ICU_jll", "JLLWrappers", "Libdl", "OpenSSL_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "b62c5dcf5d80a82e40d58b908b8eca27a54f215b"
uuid = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"
version = "0.15.0+0"
"""

# ╔═╡ Cell order:
# ╠═65c19790-a9a8-11ee-06de-ab8cb360a743
# ╟─37653018-aaf5-42d1-a938-e05a44f18918
# ╟─f028697a-0534-4ccd-8282-1b9a3d8e284b
# ╟─4b65fd9c-bd69-4e97-a28c-357c6a9c7bda
# ╠═a0ca487a-8a51-48cb-acc0-7b318c793d09
# ╠═176d2173-b440-4d03-822b-60034921cade
# ╟─33e83e53-1df4-48a9-85fd-17b50c1b1be6
# ╠═0c6f9d84-30d5-4b4d-8efd-f7414205aa5e
# ╟─5eadbd03-da0b-48e1-9027-bc3243389a8d
# ╟─0f03beab-d57c-4e7c-8dd1-e5eb3425c0c2
# ╠═5769c1c3-d184-41b2-9f23-d490fb1a9df1
# ╟─798cafe1-58e6-4c30-8668-59615b1d26ba
# ╠═b94ba6f9-f136-4ca9-bb21-207df7a5534b
# ╟─65b56c83-f48a-4df5-8e52-326e40c56e78
# ╠═0db83126-5520-4438-9f6b-2a6839559787
# ╠═93d71685-815f-4faf-bff8-3bb35945d2bf
# ╟─5db32b26-0485-4929-89ec-34c09450555e
# ╠═871a53c5-82aa-4a4e-9627-c759901e7a89
# ╠═76d42a07-ae3f-42f9-8a97-9b2892e55088
# ╟─95fc4af6-247b-43fd-b271-5f3a588f3dff
# ╠═5fa18c9e-0ae1-402c-aa28-7396a4973105
# ╠═22361146-afa5-40a5-ad79-9b7b91cba8db
# ╟─73503ce5-c8b4-4f1d-a95c-cc7a307fdb8f
# ╠═2017166e-2932-48f3-a9e2-f5423f018742
# ╠═be89676f-0e5a-48c5-b8e3-61b2d7b892a5
# ╟─7e5980f4-2bcc-4c42-a244-5a8e45df9ae1
# ╠═85d9b763-ded9-4e95-be86-e519abe904d9
# ╠═8d718f6f-1ace-4d5a-8c9d-6b1994cb014d
# ╟─2e4eda53-ef7b-411e-893d-c47de1f4b4df
# ╠═8aea8d2e-90ef-4302-875d-d99f952a6ca0
# ╠═8c7ed540-9e7d-428f-bb68-ae7fe91fef42
# ╟─e36e8f05-138f-4a8f-be3a-23edbac61304
# ╠═a1e9bf3b-8bf9-4425-a57c-61876cfdaa7a
# ╠═c90f69cf-399e-4831-818d-f6d34036b641
# ╟─d4af2366-fa8d-4280-9057-15699d585daa
# ╠═6dbbe165-1177-405d-a3cf-661cf96b265e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
