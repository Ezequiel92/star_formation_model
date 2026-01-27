### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ fed88caa-1520-41f7-adb3-785e5c9529c6
using DataFrames, CSV, DifferentialEquations, Interpolations, LinearAlgebra, Measurements, NaNMath, PlutoUI, QuadGK, SpecialFunctions, Symbolics, TikzPictures, Unitful, UnitfulAstro

# ╔═╡ 734b3b08-061e-4f93-8574-468d824815da
# ╠═╡ skip_as_script = true
#=╠═╡
TableOfContents(title="🌌 SF model", depth=4)
  ╠═╡ =#

# ╔═╡ 800dc762-ce0b-463f-858f-6e8eabbc26b0
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Introduction

## Motivation

In the prevailing cosmological framework, dark matter haloes form from initial density fluctuations. Baryonic gas subsequently collapses within these haloes, and galaxies grow through the accretion of material and the merging of smaller structures ([Rees1977](https://doi.org/10.1093/mnras/179.4.541), [Silk1977](https://doi.org/10.1086/154972), [White1978](https://doi.org/10.1093/mnras/183.3.341), [Blumenthal1984](https://doi.org/10.1038/311517a0)). From the onset of galaxy formation, a range of physical processes -- such as gravitational collapse, gas cooling, mass accretion, star formation (SF), chemical enrichment, and feedback -- operate in tandem to influence the evolution and observable properties of galaxies ([Dalgarno1972](https://doi.org/10.1146/annurev.aa.10.090172.002111), [Efstathiou1992](https://doi.org/10.1093/mnras/256.1.43P), [White1991](https://doi.org/10.1086/170483)). Additionally, interactions and mergers can drive significant changes in morphology and affect the chemical, dynamical, and structural development of galaxies.

While no comprehensive analytical theory yet captures all aspects of galaxy evolution, cosmological simulations have proven essential in advancing our understanding of galaxy formation in a realistic context. Contemporary simulations can now reproduce systems with structural and kinematic features broadly consistent with observations, including those of the Milky Way ([Grand2017](https://doi.org/10.1093/mnras/stx071), [Springel2017](https://doi.org/10.1093/mnras/stx3304), [Libeskind2020](https://doi.org/10.1093/mnras/staa2541), [Vogelsberger2020](https://doi.org/10.1038/s42254-019-0127-2)).

A key challenge in modeling galaxy formation is the treatment of star formation and its associated feedback processes. Stars form primarily in dense gas clouds -- often concentrated in the centers of dark matter haloes and along galactic disks -- once gas cooling enables their collapse. As stars form, they enrich the interstellar medium (ISM) with heavy elements and release energy through supernova explosions, which in turn influence further SF activity. The balance between cooling and heating regulates the star formation rate (SFR), while additional processes such as gas inflow, galactic fountains, internal instabilities, and environmental interactions further affect the ISM in complex ways. One major advantage of cosmological simulations is their ability to capture these processes self-consistently. However, because star formation and stellar feedback occur on scales below the resolution of current simulations, sub-grid models must be introduced. These models allow simulations to include the macroscopic effects of unresolved physics while maintaining computational feasibility on cosmological volumes.

## Previous work

Early models of star formation in simulations relied on a simple empirical relationship between the gas density and the SFR, first introduced by [Katz1992](https://doi.org/10.1086/171366) and [Steinmetz1994](https://doi.org/10.48550/arXiv.astro-ph/9312010). In these models, star particles form from cold, dense gas at a rate proportional to the gas density and inversely proportional to a characteristic timescale. When combined with effective feedback, this prescription successfully reproduces the observed slope and normalization of the SFR-gas surface density relation ([Schmidt1959](https://doi.org/10.1086/146614), [Kennicutt1998](https://doi.org/10.1086/305588)), provided that a free parameter -- the star formation efficiency (SFE) -- is calibrated, and the star formation timescale is assumed to be the gas dynamical or free-fall time ([Scannapieco2006](https://doi.org/10.1111/j.1365-2966.2006.10785.x)).

Observational evidence has shown that the correlation between gas density and SFR becomes stronger when focusing on molecular gas ([Wong2002](https://doi.org/10.1086/339287), [Leroy2008](https://doi.org/10.1088/0004-6256/136/6/2782), [Bolatto2011](https://doi.org/10.1088/0004-637X/741/1/12), [Schruba2011](https://doi.org/10.1088/0004-6256/142/2/37), [Robertson2008](https://doi.org/10.1086/587796), [Thompson2013](https://doi.org/10.1088/0004-637x/780/2/145)). This holds true both at resolved scales ([Baker2021](https://doi.org/10.1093/mnras/stab3672)) and on galaxy-wide scales across different redshifts ([Baker2022](https://doi.org/10.1093/mnras/stac3413)). Additionally, the typical SF timescale can vary depending on the local ISM conditions ([Bigiel2008](https://doi.org/10.1088/0004-6256/136/6/2846), [Bigiel2010](https://doi.org/10.1088/0004-6256/140/5/1194), [Bolatto2011](https://doi.org/10.1088/0004-637X/741/1/12), [Leroy2013](https://doi.org/10.1088/0004-6256/146/2/19)). In response, more sophisticated SF models have been introduced in both simulations and semi-analytic frameworks. These improvements often aim to resolve or better characterize the ISM, enabling models to explicitly follow molecular hydrogen and incorporate $\mathrm{H}_2$-based SF laws ([Murante2010](https://doi.org/10.1111/j.1365-2966.2010.16567.x), [Molla2015](https://doi.org/10.1093/mnras/stv1102), [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635)).

Modeling stellar feedback, in contrast, remains more difficult due to the inherent complexity of the underlying physics and its interactions with the ISM ([Agertz2013](https://doi.org/10.1088/0004-637X/770/1/25)). Different implementations can lead to widely varying predictions of galaxy properties, including stellar mass and morphology, even when starting from the same initial conditions ([Scannapieco2012](https://doi.org/10.1111/j.1365-2966.2012.20993.x)). Nonetheless, considerable progress has been made in improving feedback prescriptions ([Kim2013](https://doi.org/10.1088/0067-0049/210/1/14), [Kim2016](https://doi.org/10.3847/1538-4357/833/2/202)). Since the strength and impact of feedback are closely tied to star formation, the associated numerical parameters are typically model-specific and cannot be directly constrained by observations. Other feedback channels -- such as those from black holes and cosmic rays -- have also been studied, though stellar feedback remains the dominant mechanism for galaxies with masses up to that of the Milky Way ([Naab2017](https://doi.org/10.1146/annurev-astro-081913-040019)).

Crucially, realistic modeling of the ISM is essential to determine the conditions that regulate SF. Several sub-grid models have been proposed to describe the ISM's multiphase structure. For example, [Yepes1997](https://doi.org/10.1093/mnras/284.1.235) and [Hultman1999](https://ui.adsabs.harvard.edu/abs/1999A%26A...347..769H) adapted the classic two-phase ISM model of [McKee1977](https://doi.org/10.1086/155667) in both Eulerian and Lagrangian codes, treating each gas element as a mix of hot and cold phases in pressure equilibrium. Cold clouds form via thermal instability, leading to SF and subsequent feedback that redistributes mass between phases. [Springel2003](https://doi.org/10.1046/j.1365-8711.2003.06206.x) expanded upon these ideas in the \$\texttt{GADGET}\$ code, incorporating starburst-driven galactic winds as an additional feedback mechanism.

More recent implementations have aimed to track molecular hydrogen explicitly. [Pelupessy2006](https://doi.org/10.1086/504366) introduced a sub-grid model for the time-dependent transition from atomic to molecular hydrogen and applied it to isolated dwarf galaxies. [Robertson2008](https://doi.org/10.1086/587796) assumed equilibrium $\mathrm{H}_2$ abundances dependent on metallicity, while [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) used non-equilibrium chemistry in adaptive mesh refinement (AMR) simulations of low-metallicity systems. Their approach was later adapted for SPH codes by [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x), who used it in cosmological simulations of dwarf galaxies up to $z = 0$. [Murante2014](https://doi.org/10.1093/mnras/stu2400) applied a pressure-based $\mathrm{H}_2$ prescription ([Blitz2006](https://doi.org/10.1086/505417)) to a Milky Way-like halo, while [Valentini2022](https://doi.org/10.1093/mnras/stac2110) compared two different molecular fraction models -- one theoretical ([Krumholz2009](https://doi.org/10.1088/0004-637X/699/1/850)) and one phenomenological ([Blitz2006](https://doi.org/10.1086/505417)) -- and found that the former yields a clumpier, more observationally consistent ISM. [Hopkins2014](https://doi.org/10.1093/mnras/stu1738) also used an $\mathrm{H}_2$-regulated SF law, computing molecular fractions from local column density and metallicity following [Krumholz2011](https://doi.org/10.1088/0004-637X/729/1/36)).

A selected list of previous simulations can be seen below

| Reference     | Code |
|:-------------:|:----:|
| [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) and [Gnedin2011](https://doi.org/10.1088/0004-637X/728/2/88) | $\texttt{ART}$ ([Kravtsov1997](https://doi.org/10.1086/313015)) |
| [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x) | $\texttt{GASOLINE}$ ([Wadsley2004](https://doi.org/10.1016/j.newast.2003.08.004)) |
| [Tomassetti2014](https://doi.org/10.1093/mnras/stu2273) | $\texttt{RAMSES}$ ([Teyssier2001](https://doi.org/10.1051/0004-6361:20011817)) |
| [Baczynski2015](https://doi.org/10.1093/mnras/stv1906) | $\texttt{FLASH4}$ ([Fryxell2000](https://doi.org/10.1086/317361)) |
| [Richings2016](https://doi.org/10.1093/mnras/stw327) | $\texttt{GADGET3}$ ([Springel2005](https://doi.org/10.1111/j.1365-2966.2005.09655.x)) |
| [Hu2016](https://doi.org/10.1093/mnras/stw544) | $\texttt{GADGET3}$ ([Springel2005](https://doi.org/10.1111/j.1365-2966.2005.09655.x)) |
| [Katz2017](https://doi.org/10.1093/mnras/stx608) | $\texttt{RAMSES-RT}$ ([Rosdahl2013](https://doi.org/10.1093/mnras/stt1722)) |
| [Pallottini2017](https://doi.org/10.1093/mnras/stx608) | $\texttt{RAMSES}$ ([Teyssier2001](https://doi.org/10.1051/0004-6361:20011817)) |
| [Lupi2017](https://doi.org/10.1093/mnras/stx2874) | $\texttt{GASOLINE2}$ ([Wadsley2017](https://doi.org/10.1093/mnras/stx1643)) |
| [Capelo2018](https://doi.org/10.1093/mnras/stx3355) | $\texttt{GIZMO}$ ([Hopkins2015](https://doi.org/10.1093/mnras/stv195)) |
| [Nickerson2018](https://doi.org/10.1093/mnras/sty1556) and [Nickerson2019](https://doi.org/10.1093/mnras/stz048) | $\texttt{RAMSES-RT}$ ([Rosdahl2013](https://doi.org/10.1093/mnras/stt1722)) |
| [Sillero2021](https://doi.org/10.1093/mnras/stab1015) | $\texttt{GADGET3}$ ([Springel2005](https://doi.org/10.1111/j.1365-2966.2005.09655.x)) |

## $\texttt{AREPO}$

To perform our simulations, we used the $N$-body magnetohydrodynamics (MHD) code $\texttt{AREPO}$ ([Springel2010](https://doi.org/10.1111/j.1365-2966.2009.15715.x)). $\texttt{AREPO}$ is a moving-mesh code designed to evolve MHD and collisionless systems in cosmological volumes. Gravitational forces are calculated using a hybrid TreePM scheme ([Springel2005](https://doi.org/10.1111/j.1365-2966.2005.09655.x)) that combines long-range forces from a particle-mesh (PM) method with short-range forces from an oct-tree algorithm ([Barnes1986](https://doi.org/10.1038/324446a0)), using adaptive time stepping for increased efficiency.

The code employs a dynamic, unstructured mesh based on a Voronoi tessellation of space, enabling a finite-volume solution to the MHD equations. These equations are integrated using a second-order Runge-Kutta scheme and gradient estimates based on a least-squares method for high accuracy ([Pakmor2015](https://doi.org/10.1093/mnras/stv2380)).

A notable feature of $\texttt{AREPO}$ is its ability to perform mesh reconstruction at every time step. This ensures that each cell retains a target mass within a set tolerance, naturally achieving higher resolution in dense regions. The mesh-generating points also move with the gas flow, mitigating the Galilean invariance issues of Eulerian codes and reducing advection errors in supersonic regimes. This quasi-Lagrangian behavior makes $\texttt{AREPO}$ competitive with Smoothed Particle Hydrodynamics (SPH), while overcoming some of SPH's limitations, such as the need for artificial viscosity.

The $\texttt{AREPO}$ framework includes a variety of physical modules. In this work, we adopt the galaxy formation model from the Auriga Project, which includes metal-line and primordial cooling, a uniform UV background, magnetic fields ([Pakmor2014](https://doi.org/10.1088/2041-8205/783/1/L20), [Pakmor2017](https://doi.org/10.1093/mnras/stx1074), [Pakmor2018](https://doi.org/10.1093/mnras/sty2601)), and stellar feedback from both Type II and Type Ia supernovae, as well as mass and metal return from AGB stars ([Vogelsberger2013](https://doi.org/10.1093/mnras/stt1789), [Marinacci2013](https://doi.org/10.1093/mnras/stt2003), [Grand2017](https://doi.org/10.1093/mnras/stx071)).
"""
  ╠═╡ =#

# ╔═╡ b842e98e-34e2-40f2-84b6-c180815c2df3
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Star formation model

## Phases and notation

We modeled the interstellar medium as a multiphase structure made up of six components. Three correspond to different phases of hydrogen: ionized ($\sim 10^4 \, \mathrm{K}$), atomic ($\sim 100 \, \mathrm{K}$), and molecular ($\sim 10 \, \mathrm{K}$) gas. The other components are: metals, stars, and dust.

For each mass-transfer reaction between phases, only the dominant pathway is considered. Although the ISM consists primarily of hydrogen and helium, we model only hydrogen-related reactions, as they govern the key processes involved. These include the photoionization of atoms, the recombination of electrons with ions, the formation of molecular hydrogen from atomic hydrogen, and its destruction via photodissociation by ultraviolet radiation. Additionally, we account for the production of ionized gas and metals by supernovae, the formation and destruction of dust from and into metals, and the formation of stars from molecular gas.

We characterize each phase by its mass fraction with respect to the total mass of the cell,

|               |                                      |
|:-------------:|:------------------------------------:|
| Ionized gas   | $f_i(t) = M_i(t) / M_\mathrm{cell}$ |
| Atomic gas    | $f_a(t) = M_a(t) / M_\mathrm{cell}$ |
| Molecular gas | $f_m(t) = M_m(t) / M_\mathrm{cell}$ |
| Stars         | $f_s(t) = M_s(t) / M_\mathrm{cell}$ |
| Metals        | $f_Z(t) = M_Z(t) / M_\mathrm{cell}$ |
| Dust          | $f_d(t) = M_d(t) / M_\mathrm{cell}$ |

where $M_i$, $M_a$, $M_m$, $M_s$, $M_Z$, and $M_d$ are the corresponding masses and

$\begin{equation}
    M_\mathrm{cell} = M_i(t) + M_a(t) + M_m(t) + M_s(t) + M_Z(t) + M_d(t) \, ,
\end{equation}$

is the total cell mass.

It is important to clarify that these fractions describe the internal composition of each gas cell within our subgrid model. The stellar mass fraction, $f_s$, represents the mass of stars formed within a gas cell that has not yet been converted into a distinct star particle. It is used to calculate the cell's instantaneous SFR, which determines if a new star particle will be generated.

Now, we will show the relation between the number density of the $j$ component, $n_j$, and the dimensionless fractions defined above,

$\begin{equation}
    n_j = \frac{N_j}{V_j} = \frac{M_j}{m_j} \, \frac{1}{V_j} = \frac{M_j}{M_\mathrm{cell}} \,  \frac{M_\mathrm{cell}}{m_j \, V_j} = f_j \, \frac{M_\mathrm{cell}}{m_j \, x_j \, V_\mathrm{cell}} = f_j \, \frac{\rho_\mathrm{cell}}{m_j \, x_j} \, ,
\end{equation}$

where the quantities for the gas cell are

|                      |              |
|:--------------------:|:------------:|
| $V_\mathrm{cell}$    | Volume       |
| $M_\mathrm{cell}$    | Total mass   |
| $\rho_\mathrm{cell}$ | Mass density |

and for the $j$ component are

|       |                                                  |
|:-----:|:------------------------------------------------:|
| $N_j$ | Number of elements (e.g. atoms)                  |
| $V_j$ | Volume                                           |
| $M_j$ | Total mass                                       |
| $m_j$ | Mass of a single element (e.g. atom)             |
| $f_j$ | Mass fraction ($f_j = M_j / M_\mathrm{cell}$)   |
| $x_j$ | Volume fraction ($x_j = V_j / V_\mathrm{cell}$) |

In the context of our model, $V_\mathrm{cell}$, $M_\mathrm{cell}$, $\rho_\mathrm{cell}$, $V_j$, $m_j$, and $x_j$ are constants during the integration of the ODEs.

In the same way, we can write a relation for the mass density of the $j$ component,

$\begin{equation}
    \rho_j = \frac{M_j}{V_j} = \frac{M_j}{M_\mathrm{cell}} \, \frac{M_\mathrm{cell}}{V_j} = f_j \,  \frac{M_\mathrm{cell}}{x_j \, V_\mathrm{cell}} = f_j \, \frac{\rho_\mathrm{cell}}{x_j} \, .
\end{equation}$

So, using these relations we can write any differential equation for the quantities $M_j$, $\rho_j$, and $n_j$ as equations for $f_j$.

So far, we have made only two assumptions: first, that the ISM consists solely of the six components previously described; and second, that $V_\mathrm{cell}$, $M_\mathrm{cell}$, $\rho_\mathrm{cell}$, $V_j$, $m_j$, and $x_j$ are treated as constants.

For simplicity we will adopt $x_i = x_a = x_m = 1$. This implies that the ionized, atomic, and molecular hydrogen phases are well mixed and collectively fill the entire cell volume.

## Volume fractions in the Milky Way

Although we assume $x_i = x_a = x_m = 1.0$, for completeness we show the estimates for the volume fractions in the MW.

Using the values from [Ferrière2001](https://doi.org/10.1103/RevModPhys.73.1031) we have
"""
  ╠═╡ =#

# ╔═╡ 3c0c1fea-1ae6-4c9d-9db1-764937cd5fbc
# ╠═╡ skip_as_script = true
#=╠═╡
let
	# Approximate radius for the MW
	R = (25.0u"kpc" + 30.0u"kpc") / 2.0

	# Approximate height for the MW
	h = (400.0u"pc" + 600.0u"pc") / 2.0

	# Volume (as a disc) of the MW
	V = π * R^2 * h

	# Masses of the hydrogen phases in the MW
	Mm = (1.3e9u"Msun" + 2.5e9u"Msun") / 2.0
	Ma = 6.0e9u"Msun"
	Mi = 1.6e9u"Msun"

	# Number densities of the hydrogen phases in the MW
	nm = exp10((2.0 + 6.0) / 2.0) * u"cm^-3"
	na = (20.0u"cm^-3" + 50.0u"cm^-3") / 2.0
	ni = (0.2u"cm^-3" + 0.5u"cm^-3") / 2.0

	# Reference mass of an 'element' of each hydrogen phase
	mm = 2.0u"mp"
	ma = 1.0u"mp"
	mi = 1.0u"mp"

	# Volume fractions
	xm = Mm / (mm * nm * V)
	xa = Ma / (ma * na * V)
	xi = Mi / (mi * ni * V)

	xt = xm + xa + xi

	# Notice how V cancels when renormalizing, so its specific value does not matter
	xmr = round((xm / xt) * 100, sigdigits=3)
	xar = round((xa / xt) * 100, sigdigits=3)
	xir = round((xi / xt) * 100, sigdigits=3)

	println("Volume fractions in the MW:\n")
	println("\tMolecular hydrogen: $(xmr)%")
	println("\tAtomic hydrogen:    $(xar)%")
	println("\tIonized hydrogen:   $(xir)%")
end;
  ╠═╡ =#

# ╔═╡ 99d0fdc9-e368-462c-a357-86f07624a52e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Physical processes
"""
  ╠═╡ =#

# ╔═╡ 6d6cbebb-bb3b-437d-9835-0a79f36857f2
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### By name
"""
  ╠═╡ =#

# ╔═╡ 71b94af7-76ee-4228-987c-2f22e0951552
# ╠═╡ skip_as_script = true
#=╠═╡
TikzPictures.TikzPicture(
	L"""
		\node[box, white] (stars) {Stars};
		\node[box, white, text width=2em, above=2.5cm of stars] (atom) {HI};
		\node[box, white, text width=2em, right=2cm of atom] (molecule) {\ch{H2}};
		\node[box, white, text width=2em, left=2cm of atom] (ion) {HII};
		\node[box, white, text width=2em, below=0.8cm of stars] (metals) {Z};
		\node[box, white, text width=2em, below=0.8cm of metals] (dust) {Dust};
		\draw[line, white, ->]
		(ion) edge [bend left, "\textcolor{d_pink}{recombination}"] (atom)
		(atom) edge [bend left, "\textcolor{d_orange}{condensation}"] (molecule)
		(molecule) edge [bend left,"\textcolor{d_green}{dissociation}"] (atom)
		(atom) edge [bend left,"\textcolor{d_blue}{ionization}"] (ion)
		(stars) edge [bend left, "\textcolor{d_yellow}{supernova}"] (ion)
		(molecule) edge [bend left, "\textcolor{red}{star formation}"] (stars)
		(stars) edge ["\textcolor{d_yellow}{supernova}"] (metals)
		(metals) edge [bend left, "\textcolor{g_green}{dust growth}"] (dust)
		(dust) edge [bend left, "\textcolor{g_red}{dust destruction}"] (metals);
	""",
	width="75em",
	preamble = """
		\\usepackage{chemformula}
		\\definecolor{d_pink}{HTML}{C721DD}
		\\definecolor{d_orange}{HTML}{D14A00}
		\\definecolor{d_green}{HTML}{008C00}
		\\definecolor{d_blue}{HTML}{007FB1}
		\\definecolor{d_yellow}{HTML}{D1AC00}
		\\definecolor{g_green}{HTML}{88b57b}
		\\definecolor{g_red}{HTML}{b57b7b}
		\\usetikzlibrary{shapes.misc, arrows, positioning, quotes, fit}
		\\tikzset{
    		>=stealth',
    		box/.style={
        		rectangle,
        		rounded corners,
        		draw=black,
        		thick,
        		text width=4em,
        		minimum height=2em,
        		text centered,
    		},
			line/.style = {
				thick,
			},
			every edge quotes/.append style = {
				font=\\small,
				align=center,
				auto,
			},
			myrect/.style={
				rectangle,
				draw,
				inner sep=0pt,
				fit=\\#1,
				thick,
				rounded corners,
			},
		}
	""",
)
  ╠═╡ =#

# ╔═╡ 22c37732-1cf7-4d80-a250-9fd5e4a2f88c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### By equation
"""
  ╠═╡ =#

# ╔═╡ 8da87629-09ff-4677-b5e8-2a3ff5777a8e
# ╠═╡ skip_as_script = true
#=╠═╡
TikzPictures.TikzPicture(
    L"""
        \node[box, white] (stars) at (270:2cm) {Stars};
        \node[box, white, text width=2em] (atom) at (90:2cm) {HI};
        \node[box, white, text width=2em] (molecule) at (0:2cm) {\ch{H2}};
        \node[box, white, text width=2em] (ion) at (180:2cm) {HII};
        \node[box, white, text width=2em] (metals) at (270:4.5cm) {Z};
        \node[box, white, text width=2em] (dust) at (270:7cm) {Dust};
        \draw[line, white, ->]
        (ion) edge [bend left, "$\textcolor{d_pink}{\dfrac{f_i}{\tau_\mathrm{rec}}} - \textcolor{d_blue}{e_\mathrm{ion} \, \eta_\mathrm{ion}^* \, \psi} - \textcolor{d_blue}{e_\mathrm{ion} \, \Gamma_\mathrm{UVB} \, f_a}$"] (atom)
        (atom) edge [bend left, "$\textcolor{d_orange}{\dfrac{f_a}{\tau_\mathrm{cond}}} - \textcolor{d_green}{e_\mathrm{diss} \, \eta_\mathrm{diss}^* \, \psi} - \textcolor{d_green}{e_\mathrm{diss} \, \Gamma_\mathrm{LWB} \, f_m}$"] (molecule)
        (stars) edge [bend left, "$\textcolor{d_yellow}{R \, (1 - Z_\mathrm{SN}) \, \psi(t)}$"] (ion)
        (molecule) edge [bend left, "$\textcolor{red}{\psi}$"] (stars)
        (stars) edge node[midway, xshift=-10mm] {$\textcolor{d_yellow}{R \, Z_\mathrm{SN} \, \psi}$} (metals)
        (metals) edge node[midway, xshift=-25mm] {$\textcolor{g_red}{\dfrac{f_d}{\tau_\mathrm{dd}}} - \textcolor{g_green}{\dfrac{f_Z}{f_Z + f_d}\, \dfrac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m)}$} (dust);
    """,
    width="70em",
    preamble = """
        \\usepackage{chemformula}
        \\definecolor{d_pink}{HTML}{C721DD}
        \\definecolor{d_orange}{HTML}{D14A00}
        \\definecolor{d_green}{HTML}{008C00}
        \\definecolor{d_blue}{HTML}{007FB1}
        \\definecolor{d_yellow}{HTML}{D1AC00}
        \\definecolor{g_green}{HTML}{88b57b}
        \\definecolor{g_red}{HTML}{b57b7b}
        \\usetikzlibrary{shapes.misc, arrows, positioning, quotes, fit}
        \\tikzset{
            >=stealth',
            box/.style={
                rectangle,
                rounded corners,
                draw=black,
                thick,
                text width=4em,
                minimum height=2em,
                text centered,
            },
            line/.style = {
                thick,
            },
            every edge quotes/.append style = {
                font=\\small,
                align=center,
                auto,
            },
            myrect/.style={
                rectangle,
                draw,
                inner sep=0pt,
                fit=\\#1,
                thick,
                rounded corners,
            },
        }
    """,
)
  ╠═╡ =#

# ╔═╡ c5e21675-f120-4555-84be-d99b35592f2d
# ╠═╡ skip_as_script = true
#=╠═╡
md"""

## Equations

The equations presented here describe only the gain terms for each component -- that is, how mass is added through various physical processes. To ensure mass conservation in the model, each gain term is accompanied by a corresponding loss term in the source component. For example, the formation of stars reduces the molecular fraction, while ionization decreases the atomic fraction and increases the ionized fraction. In the full system of differential equations, each such process appears twice: once as a positive term in the equation for the receiving phase, and once with a negative sign in the equation for the source phase. This symmetric formulation guarantees that the total mass within each cell remains conserved under the assumptions of the model.

### Stars

We define $\psi(t)$ as the fractional SFR (SFR per unit of cell mass),

$\begin{equation}
    \left. \frac{\mathrm{d}}{\mathrm{d}t}f_s \right|_{\mathrm{formation}} = \frac{\mathrm{SFR}}{M_\mathrm{cell}} =: \psi \, ,
\end{equation}$

### Ionized gas

The ionized component grows through the ionization of atomic gas and from the remnants of supernova explosions.

The former is produced by the radiation of newborn stars (within the cell) and from the metagalactic ultraviolet background radiation (UVB), so it can be written as

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_i\right|_{\mathrm{ionization}} = e_\mathrm{ion} \, \eta_\mathrm{ion}^* \, \psi + e_\mathrm{ion} \, \Gamma_\mathrm{UVB} \, f_a \, ,
\end{equation}$

where $e_\mathrm{ion} = S_d$ is the ionization shielding factor, $\eta_\mathrm{ion}$ is the ionized mass rate per unit of created stellar mass, and $\Gamma_\mathrm{UVB}$ is the UVB photoionization rate.

The latter, under the instantaneous recycling hypothesis, can be written as

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_i\right|_{\mathrm{recycling}} = R \, (1 - Z_\mathrm{SN}) \, \psi \, ,
\end{equation}$

where $R$ is the total (ionized gas + metals) mass fraction of a stellar population that is returned to the ISM, and $Z_\mathrm{SN}$ is the fraction of the returned mass that is metals.

### Atomic gas

The atomic component grows through the dissociation of hydrogen molecules and the recombination of the ionized gas with free electrons.

The former can be written as

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_a\right|_{\mathrm{dissociation}} = e_\mathrm{diss} \, \eta_\mathrm{diss}^* \, \psi + e_\mathrm{diss} \, \Gamma_\mathrm{LWB} \, f_m \, ,
\end{equation}$

where $e_\mathrm{diss} = S_d \, S_{\mathrm{H}_2}$ is the dissociation shielding factor, $S_{\mathrm{H}_2}$ is the molecular self shielding factor, $\eta_\mathrm{diss}$ is the disassociated mass rate per unit of created stellar mass, and $\Gamma_\mathrm{LWB}$ is the disassociation rate due to Lyman–Werner backgound radiation.

The latter will depend on the mass of ionized gas and on the characteristic time scale of recombination, $\tau_\mathrm{rec}$, so it can be written as

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_a\right|_{\mathrm{recombination}} = \frac{f_i}{\tau_\mathrm{rec}} \, .
\end{equation}$

### Molecular gas

The molecular component grows mainly through the condensation of hydrogen atoms on the surface of dust grains. This process depends on the mass of atomic gas and on the characteristic time scale of condensation, $\tau_\mathrm{cond}$, so it can be written as

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_m\right|_{\mathrm{condensation}} = \frac{f_a}{\tau_\mathrm{cond}} \, .
\end{equation}$

### Metals

The metals in the gas come from the dust and from supernovas.

The former depends on the mass of dust and on a characteristic time scale of dust destruction, $\tau_\mathrm{dd}$, so it can be written as

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_Z\right|_{\mathrm{dust}} = \frac{f_d}{\tau_\mathrm{dd}} \, .
\end{equation}$

While the latter is

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_Z\right|_{\mathrm{recycling}} = R \, Z_\mathrm{SN} \, \psi \, .
\end{equation}$

### Dust

The dust grows from the condensation of metals, and can be written as

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_d\right|_{\mathrm{accretion}} = \frac{f_Z}{f_Z + f_d} \, \frac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m) \, ,
\end{equation}$

where $\tau_\mathrm{dg}$ is the characteristic time scale of dust growth.

From all the above, we can write the system of six ODEs
"""
  ╠═╡ =#

# ╔═╡ d12b4a5f-a6ad-4c91-8795-fe8f93f5c93d
# ╠═╡ skip_as_script = true
#=╠═╡
TikzPictures.TikzPicture(
	L"""
	\node[white] {
  	${\boldmath
	\begin{aligned}
		\dv{}{t}f_i &= \textcolor{d_blue}{S_d \, \eta_\mathrm{ion}^* \, \psi} + \textcolor{d_blue}{S_d \, \Gamma_\mathrm{UVB} \, f_a} + \textcolor{d_yellow}{R \, (1 - Z_\mathrm{SN}) \, \psi} - \textcolor{d_pink}{\frac{f_i}{\tau_\mathrm{rec}}} \, , \\
		\dv{}{t}f_a &= \textcolor{d_pink}{\frac{f_i}{\tau_\mathrm{rec}}} + \textcolor{d_green}{e_\mathrm{diss} \, \Gamma_\mathrm{LWB} \, f_m} + \textcolor{d_green}{e_\mathrm{diss} \, \eta_\mathrm{diss}^* \, \psi} - \textcolor{d_blue}{S_d \, \eta_\mathrm{ion}^* \, \psi} - \textcolor{d_blue}{S_d \, \Gamma_\mathrm{UVB} \, f_a} - \textcolor{d_orange}{\frac{f_a}{\tau_\mathrm{cond}}} \, , \\
		\dv{}{t}f_m &= \textcolor{d_orange}{\frac{f_a}{\tau_\mathrm{cond}}} - \textcolor{d_green}{e_\mathrm{diss}  \, \eta_\mathrm{diss}^* \, \psi} - \textcolor{d_green}{e_\mathrm{diss}  \, \Gamma_\mathrm{LWB} \, f_m} - \textcolor{red}{\psi} \, , \\
		\dv{}{t}f_s &= \textcolor{red}{\psi} - \textcolor{d_yellow}{R \, \psi} \, , \\
		\dv{}{t}f_Z &= \textcolor{d_yellow}{R \, Z_\mathrm{SN} \, \psi} + \textcolor{g_red}{\frac{f_d}{\tau_\mathrm{dd}}} - \textcolor{g_green}{\frac{f_Z}{f_Z + f_d} \, \frac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m)} \, , \\
		\dv{}{t}f_d &= \textcolor{g_green}{\frac{f_Z}{f_Z + f_d} \, \frac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m)} - \textcolor{g_red}{\frac{f_d}{\tau_\mathrm{dd}}} \, .
	\end{aligned}}$
	};
	""",
	width="70em",
	preamble = """
		\\usepackage{chemformula}
		\\usepackage{physics}
		\\setlength{\\jot}{10pt}
		\\definecolor{d_pink}{HTML}{C721DD}
		\\definecolor{d_orange}{HTML}{D14A00}
		\\definecolor{d_green}{HTML}{008C00}
		\\definecolor{d_blue}{HTML}{007FB1}
		\\definecolor{d_yellow}{HTML}{D1AC00}
		\\definecolor{g_green}{HTML}{88b57b}
		\\definecolor{g_red}{HTML}{b57b7b}
	""",
)
  ╠═╡ =#

# ╔═╡ c7eee8bb-b03d-4679-a51e-36fc629294e1
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
And we can explicitly check for mass conservation

$\begin{equation}
    \frac{\mathrm{d}}{\mathrm{d}t}(f_i + f_a + f_m + f_s + f_Z + f_d) = 0 \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 700015f6-aac0-40c5-a455-8d98c1227049
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Units

We have the freedom to choose three independent units in this model. Time, mass, and length. The choice is reflected in the constants $C_\mathrm{star}$, $C_\mathrm{rec}$, $C_\mathrm{cond}$, $C_\mathrm{dg}$, $C_\mathrm{sd}$, and $C_\mathrm{sh2}$.

We will use $\mathrm{[T] = Myr}$, $\mathrm{[M] = mp}$ (proton mass) and $\mathrm{[L] = cm}$.
"""
  ╠═╡ =#

# ╔═╡ 207305eb-4496-4518-a5bc-be85173314a5
begin
	const t_u = u"Myr"
	const m_u = u"mp"
	const l_u = u"cm"
end;

# ╔═╡ a6fdec29-3438-4ebc-af81-7a3baf0175ae
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Physical constants

Below is a list of physical constants that appear explicitly in the equations

*  $\tau_\mathrm{dd}$: Timescale for the destruction of dust into metals.
*  $\eta_\mathrm{diss}$: Rate of molecular gas dissociation by stellar radiation per unit mass of newly formed stars.
*  $\eta_\mathrm{ion}$: Rate of atomic gas ionization by stellar radiation per unit mass of newly formed stars.
*  $R$: Fraction of a stellar population’s mass that is returned to the ISM as ionized gas and metals.
*  $Z_\mathrm{SN}$: Fraction of the returned mass that consists of metals.
*  $\Gamma_\mathrm{UVB}$: Photoionization rate due to the ultraviolet background (UVB).
*  $\Gamma_\mathrm{LWB}$: Photodissociation rate due to the Lyman-Werner background (LWB).

In contrast, the parameters $C_\mathrm{star}$, $C_\mathrm{rec}$, $C_\mathrm{cond}$, $C_\mathrm{dg}$, $C_\mathrm{sd}$, and $C_\mathrm{sh2}$ do not have physical meaning; they are included solely for convenience.
"""
  ╠═╡ =#

# ╔═╡ ee7baf70-c7af-4791-8277-fadc0961e0e7
# ╠═╡ skip_as_script = true
#=╠═╡
md"# Physical processes"
  ╠═╡ =#

# ╔═╡ 72416f96-8a0a-4c13-97d1-2990e49c4cb0
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Star formation

For the formation of stars we will use the notation and definitions commonly used in the field ([McKee2007](https://doi.org/10.1146/annurev.astro.45.051806.110602), [Krumholz2014](https://doi.org/10.1016/j.physrep.2014.02.001), [Krumholz2019](https://doi.org/10.1146/annurev-astro-091918-104430)).

In particular, we will follow [Krumholz2005](https://doi.org/10.1086/431734), but
assuming only a dependence on the mass of molecular gas instead of the total amount of gas, as in [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635)

$\begin{equation}
	\mathrm{SFR} = \frac{\mathrm{d}}{\mathrm{d}t} M_s = \frac{M_m}{\tau_\mathrm{star}} \, ,
\end{equation}$

where $M_m$ is the molecular mass and $\tau_\mathrm{star}$ is the characteristic timescale of star formation (defined by this very relation).

So, the fractional SFR is

$\begin{equation}
	\psi = \frac{\mathrm{SFR}}{M_\mathrm{cell}} = \frac{M_m}{M_\mathrm{cell}} \, \frac{1}{\tau_\mathrm{star}} = \frac{f_m}{\tau_\mathrm{star}} \, .
\end{equation}$

Following [Krumholz2005](https://doi.org/10.1086/431734), we write the characteristic timescale of star formation as

$\begin{equation}
    \tau_\mathrm{star} = \frac{t_\text{ff}}{\epsilon_\text{ff}} \, ,
\end{equation}$

where $\epsilon_\text{ff}$ is the star formation efficiency and $t_\text{ff}$ is the free-fall time, which is the time for a pressure-free spherical cloud to collapse into a point due to its self-gravity,

$\begin{equation}
    t_\text{ff} = \sqrt{\frac{3\pi}{32 \, G \, \rho_\mathrm{cell}}} \, ,
\end{equation}$

where $G$ is the gravitational constant and $\rho_\mathrm{cell}$ is the density of the corresponding gas cell.

In the literature, a value of $\epsilon_\text{ff} \approx 0.01$ is commonly used ([Krumholz2007](https://doi.org/10.1086/509101), [Krumholz2014](https://doi.org/10.1016/j.physrep.2014.02.001)), although there remains significant uncertainty regarding its true value ([Lee2016](https://doi.org/10.3847/1538-4357/833/2/229), [Utomo2018](https://doi.org/10.3847/2041-8213/aacf8f)). In our model, we adopt $\epsilon_\text{ff} = 1.0$ because star formation is self-regulated through the coupled evolution of gas phases, as described by the system of ODEs. This framework naturally determines the effective star formation rate. We note, however, that previous studies have shown this parameter has little impact on the global properties of simulated galaxies ([Li2018](https://doi.org/10.3847/1538-4357/aac9b8), [Brown2022](https://doi.org/10.1093/mnras/stac1164)).

With all the previous definitions, we have

$\begin{equation}
    \tau_\mathrm{star} = \frac{1}{C_\mathrm{star} \, \sqrt{\rho_\mathrm{cell}}} \, ,
\end{equation}$

where

$\begin{equation}
	\begin{split}
	    C_\mathrm{star} &= \sqrt{\frac{32 \, G}{3\pi}} \\
		&= 0.0194 \, \mathrm{mp^{-1/2} \, Myr^{-1} \, cm^{3/2}} \, .
	\end{split}
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 0baae21f-71b5-42fd-be33-31de085928d3
begin
	# Star formation constant
	const C_star = sqrt(32u"G" / 3π)

	# Our choice of units for C_star
	const cs_unit = m_u^(-1/2) * t_u^-1 * l_u^(3/2)

	# Nominal value of C_star in our choice of units
	const c_star = ustrip(Unitful.NoUnits, C_star / cs_unit)

	# Star formation time-scale
	τ_star(ρc) = 1.0 / (c_star * sqrt(ρc))

	# Fractional SFR
	@inline ψ(fm, ρc) = c_star * fm * sqrt(ρc)
end;

# ╔═╡ a7c5650c-9bbc-44a1-8d51-d473864faae7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Recombination

The rate of recombination for hydrogen can be written as ([Osterbrock2006](http://www.worldcat.org/oclc/60611705), [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55), [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x), [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635))

$\begin{equation}
    \frac{\mathrm{d}}{\mathrm{d}t} n_a \biggr\rvert_\mathrm{recombination} = \alpha_H \, n_e \, n_i \, ,
\end{equation}$

where $\alpha_H$ is the recombination coefficient, $n_e$ is the electron number density, and $n_i$ is the ionized hydrogen number density. We note that in our model $n_e = n_i$.

Using the conversion factor already mentioned, we can write

$\begin{equation}
	\begin{split}
	    \frac{\mathrm{d}}{\mathrm{d}t} f_a \biggr\rvert_\mathrm{recombination} &= \frac{m_p}{\rho_\mathrm{cell}} \, \alpha_H \, n_e \, n_i \\
		&= \frac{m_p}{\rho_\mathrm{cell}} \, \alpha_H \, f_i \, \frac{\rho_\mathrm{cell}}{m_p} \, f_i \, \frac{\rho_\mathrm{cell}}{m_p} \\
	    &= \alpha_H \, f_i^{\,2} \, \frac{\rho_\mathrm{cell}}{m_p} \, .
	\end{split}
\end{equation}$

Following the discussion in [Nebrin2023](https://doi.org/10.3847/2515-5172/acd37a), and assuming an optically thick cloud, we will use case B recombination (sum over all hydrogen states except the ground state). Using the values in table 2.1 of [Osterbrock2006](http://www.worldcat.org/oclc/60611705) for $T = 10^4 \, \mathrm{K}$ we have

$\begin{equation}
    \alpha_H = 2.6 \times 10^{-13} \, \mathrm{cm}^3 \, \mathrm{s}^{-1} \, .
\end{equation}$

So, we can write

$\begin{equation}
    \frac{\mathrm{d}}{\mathrm{d}t} f_a \biggr\rvert_\mathrm{recombination} = \alpha_H \, f_i^{\,2} \, \frac{\rho_\mathrm{cell}}{m_p} = \frac{f_i}{\tau_\mathrm{rec}} \, ,
\end{equation}$

where the time scale $\tau_\mathrm{rec}$ is

$\begin{equation}
    \tau_\mathrm{rec} = \frac{m_p}{\alpha_H \, f_i \, \rho_\mathrm{cell}} = \frac{1}{C_\mathrm{rec} \, f_i \, \rho_\mathrm{cell}} \, ,
\end{equation}$

with the constant

$\begin{equation}
	\begin{split}
		C_\mathrm{rec} &= \frac{\alpha_H}{m_p} \\
		&= 8.205 \, \mathrm{mp^{-1} \, Myr^{-1} \, cm^{3}} \, .
	\end{split}
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 49d1a7f7-2bf2-4472-94df-6247b9237ddd
begin
	# Recombination coefficient (case B at 10^4 K)
	const αH = 2.6e-13u"cm^3 * s^-1"

	# Recombination constant
	const C_rec = αH / m_u

	# Our choice of units for C_rec
	const cr_unit = m_u^-1 * t_u^-1 * l_u^3

	# Nominal value of C_rec in our choice of units
	const c_rec = ustrip(Unitful.NoUnits, C_rec / cr_unit)

	# Recombination time-scale
	τ_rec(fi, ρc) = 1.0 / (c_rec * fi * ρc)
end;

# ╔═╡ 4c975b71-a387-49cb-90d9-fc51acefc795
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Case A recombination

For completeness we show how we can compute the case A recombination coefficient, using either fits or analytical formulas from [Seaton1959](https://doi.org/10.1093/mnras/119.2.81), [Black1981](https://doi.org/10.1093/mnras/197.3.553), and [Verner1996](https://doi.org/10.1086/192284).
"""
  ╠═╡ =#

# ╔═╡ d431eb6d-429a-49d6-83b8-b6f21831ec86
# ╠═╡ skip_as_script = true
#=╠═╡
let
	###############################
	# From eq. 36 of Seaton (1959)
	###############################

	D = 5.197e-14u"cm^3*s^-1"

	λ(T) = 157890u"K" / T

	αA_seaton(T) = D * sqrt(λ(T)) * (0.4288 + 0.5 * log(λ(T)) + 0.469 * λ(T)^(-1/3))

	αA_seaton59 = αA_seaton(1e4u"K")

	#############################
	# From eq. 6 of Black (1981)
	#############################

	αA_black(T) = 4.36e-10u"cm^3*s^-1" * T^(-0.7573)

	αA_black81 = αA_black(1e4)

	#####################################
	# From eq. 4 of Verner et al. (1996)
	#####################################

	a  = 7.982e-11u"cm^3*s^-1"
	b  = 0.748
	T0 = 3.148u"K"
	T1 = 7.036e5u"K"

	sqT0(T) = sqrt(T / T0)
	sqT1(T) = sqrt(T / T1)

	αA_verner(T) = a / (sqT0(T) * (1 + sqT0(T))^(1-b) * (1 + sqT1(T))^(1+b))

	αA_verner96 = αA_verner(1e4u"K")

	###########################
	# Mean case A recobination
	###########################

	mean_case_A = (αA_seaton(1e4u"K") + αA_black(1e4) + αA_verner(1e4u"K")) / 3.0
end;
  ╠═╡ =#

# ╔═╡ 98610902-dfff-4fd1-8ffa-8484405aee20
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
$\begin{align}
	\alpha_\text{A}^\text{Seaton59} &= 4.1206 \times 10^{-13} \, \mathrm{cm^3 \, s^{-1}} \, , \\
	\alpha_\text{A}^\text{Black81} &= 4.0765 \times 10^{-13} \, \mathrm{cm^3 \, s^{-1}} \, , \\
	\alpha_\text{A}^\text{Verner96} &= 4.1923 \times 10^{-13} \, \mathrm{cm^3 \, s^{-1}} \, , \\
	\alpha_\text{A}^\text{Mean} &= 4.1298 \times 10^{-13} \, \mathrm{cm^3 \, s^{-1}} \, 
\end{align}$
"""
  ╠═╡ =#

# ╔═╡ 2d6fcef9-df4b-4eec-be5c-a8865c3a1b76
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Condensation

From [Hollenbach1971a](https://doi.org/10.1086/150754) and [Hollenbach1971b](https://doi.org/10.1086/150755) we know that the rate of molecular hydrogen formation due to the condensation of atomic gas in the surface of dust grains is

$\begin{equation}
    \frac{\mathrm{d}}{\mathrm{d}t} n_m \biggr\rvert_\mathrm{condensation} = R_d \, n_\mathrm{H} \, n_a  \, ,
\end{equation}$

where $R_d$ is the formation rate coefficient of $\mathrm{H}_2$ on dust grain, $n_\mathrm{H}$ is the hydrogen nucleus number density, and $n_a$ is the atomic hydrogen number density. $n_\mathrm{H}$ comes from the assumption that the number density of dust grains is proportional to it (see Section II of [Hollenbach1971b](https://doi.org/10.1086/150755)).

In the literature, $n_\mathrm{H}$ is generally taken as $n_\mathrm{H} = n_a + 2 \, n_m$, because most studies consider cold gas clouds ($T \sim 100 \, \mathrm{K}$) dominated by molecular and atomic hydrogen. In contrast, we consider atomic, molecular, and ionized gas within our gas cells; so, we will use $n_\mathrm{H} = n_a + 2 \, n_m + n_i$.

We also note that the expression for $\mathrm{d}n_m / \mathrm{d}t$ is only used in equilibrium equations in most of the early works ([Hollenbach1971a](https://doi.org/10.1086/150754), [Hollenbach1971b](https://doi.org/10.1086/150755), [Jura1974](https://doi.org/10.1086/152975), [Jura1975](https://doi.org/10.1086/153545), [Black1987](https://doi.org/10.1086/165740), [Sternberg1988](https://doi.org/10.1086/166664), [Goldshmidt1995](https://doi.org/10.1086/175168)), while first appearing in an actual differential equation that does not assume equilibrium in [Draine1996](https://doi.org/10.1086/177689).

Using the conversion factor already mentioned we can write

$\begin{align}
    \frac{\mathrm{d}}{\mathrm{d}t} f_m \biggr\rvert_\mathrm{condensation} &= \frac{2 \, m_p}{\rho_\mathrm{cell}} \, R_d \, \left( n_a + 2\, n_m + n_i \right) \, n_a \\
	&= \frac{2 \, m_p}{\rho_\mathrm{cell}} \, R_d \left( f_a \, \frac{\rho_\mathrm{cell}}{m_p} + 2 \, f_m \, \frac{\rho_\mathrm{cell}}{2 \, m_p} + f_i \, \frac{\rho_\mathrm{cell}}{m_p} \right) \, f_a \, \frac{\rho_\mathrm{cell}}{m_p} \\
	&= 2 \, R_d \left( f_a + f_m + f_i \right) \, f_a \, \frac{\rho_\mathrm{cell}}{m_p} \\
	&= \frac{f_a}{\tau_\mathrm{cond}} \, ,
\end{align}$

where the time scale $\tau_\mathrm{cond}$ is

$\begin{equation}
    \tau_\mathrm{cond} = \frac{m_p}{2 \, R_d \, \rho_\mathrm{cell} \, (f_a + f_m + f_i)} \, .
\end{equation}$

A table with several values for $R_d$ is presented below. We note that more than one value for $R_d$ and its dependence on other parameters may be discussed within each reference. In the table, we reflect the fiducial value used by each author

| $R_d \,\, [10^{-17} \, \mathrm{cm^3 \, s^{-1}}]$ | Reference          |
|:-----:|:-------------------------------------------------------------:|
| $1$   | eq. 3 in [Hollenbach1971b](https://doi.org/10.1086/150755)    |
| $5$   | section III in [Jura1974](https://doi.org/10.1086/152975)     |
| $3$   | section V in [Jura1975](https://doi.org/10.1086/153545)       |
| $0.9$ | eq. 6 in [Black1987](https://doi.org/10.1086/165740)          |
| $3$   | section 3 in [Goldshmidt1995](https://doi.org/10.1086/175168) |
| $0.6$ | eq. 18 in [Draine1996](https://doi.org/10.1086/177689)        |
| $3.5$ | section 4.2 in [Wolfire2008](https://doi.org/10.1086/587688)  |

Theoretical and experimental studies have shown that $R_d$ scales with the square root of temperature, $R_d \propto T^{1/2}$ (see eq. 6 in [Black1987](https://doi.org/10.1086/165740)), and is directly proportional to the dust number density, $R_d \propto n_\mathrm{dust}$ ([Hollenbach1971b](https://doi.org/10.1086/150755)). Assuming the simplest dust model, in which $n_\mathrm{dust} \propto Z$, this leads to the scaling relation $R_d \propto T^{1/2} Z$ ([Pelupessy2006](https://doi.org/10.1086/504366), [Wolfire2008](https://doi.org/10.1086/587688)).

Following previous approaches ([Pelupessy2006](https://doi.org/10.1086/504366), [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55), [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x), [Molla2017](https://doi.org/10.1093/mnras/stx419), [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635)), we adopt a simplified model where $R_d$ scales only with metallicity, and we fix the cold neutral gas temperature at approximately $100 \, \mathrm{K}$.

Our reference value for $R_d$ at solar metallicity is taken from [Wolfire2008](https://doi.org/10.1086/587688)

$\begin{equation}
    R_\odot = R_d(Z = Z_\odot) = 3.5 \times 10^{-17} \, \mathrm{cm^3 \, s^{-1}} \, ,
\end{equation}$

yielding the following relation

$\begin{equation}
    R_d = Z \, \frac{R_\odot}{Z_\odot} \, .
\end{equation}$

Variations of $R_\odot$ within one dex are not expected to significantly affect the results. However, we introduce an adjustable global factor -- analogous to the clumping factor $C_\rho$ used in [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) and [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x) -- to account for unresolved sub-cell density fluctuations. We find that a clumping factor of $C_\rho \approx 100$ is required to form a disk, which is consistent with values found under certain physical conditions ([Micic2012](https://doi.org/10.1111/j.1365-2966.2012.20477.x)).

A known issue with this formulation is that $R_d \to 0$ as $Z \to 0$, leading to a divergence in $\tau_\mathrm{cond}$. This results from the oversimplified assumption that molecular hydrogen forms exclusively on dust grains. To address this, we follow [Glover2007](https://doi.org/10.1086/519445) and modify the metallicity dependence as $Z \rightarrow Z + Z_\mathrm{eff}$, where $Z_\mathrm{eff} = 10^{-3} \, Z_\odot$. This correction both eliminates the divergence and accounts for molecular hydrogen formation in environments with metallicities below $10^{-3} \, Z_\odot$.

So, we finally have

$\begin{align}
    \tau_\mathrm{cond} &= \frac{m_p \, Z_\odot}{2 \, R_\odot \, C_\rho \, (Z + Z_\mathrm{eff}) \, \rho_\mathrm{cell} \, (f_a + f_m + f_i)} \\
	&= \frac{1}{C_\mathrm{cond} \, (Z + Z_\mathrm{eff}) \, \rho_\mathrm{cell} \, (f_a + f_m + f_i)}  \\
	&= \frac{1}{C_\mathrm{cond} \, (f_Z + f_d + Z_\mathrm{eff}) \, \rho_\mathrm{cell} \, (f_i + f_a + f_m)} \, ,
\end{align}$

where we use that the total metallicity can be written as $Z = f_Z + f_d$ and

$\begin{equation}
    \begin{split}
        C_\mathrm{cond} &= \frac{2 \, R_\odot \, C_\rho}{m_p \, Z_\odot} \\
        &= 17.394 \, \mathrm{mp^{-1} \, Myr^{-1} \, cm^3} \, .
    \end{split}
\end{equation}$

For the solar metallicity, we find several values in the literature

| $Z_\odot$ | Reference                                                            |
|:--------:|:---------------------------------------------------------------------:|
| $0.0169$ | [Draine1996](https://doi.org/10.1086/177689)                          |
| $0.0153$ | [Grevesse1998](https://doi.org/10.1023/A:1005161325181)               |
| $0.0122$ | [Asplund2006](https://doi.org/10.1553/cia147s76) and [Grevesse2007](https://doi.org/10.1007/s11214-007-9173-7)                                        |
| $0.0134$ | [Asplund2009](https://doi.org/10.1146/annurev.astro.46.060407.145222) |
| $0.0141$ | [Lodders2009](https://doi.org/10.1007/978-3-540-88055-4_34)           |
| $0.0127$ | $\texttt{AREPO}$                                                      |
| $0.0196$ | [Steiger2015](https://doi.org/10.3847/0004-637X/816/1/13)             |

To keep it consistent with $\texttt{AREPO}$, we will use $Z_\odot = 0.0127$, noting that is only $35\%$ off the largest value in the list ([Steiger2015](https://doi.org/10.3847/0004-637X/816/1/13)).
"""
  ╠═╡ =#

# ╔═╡ 568d9fe3-6716-4a9a-ba1d-b9e6fd039150
begin
	# Solar metallicity
    const Zsun = 0.0127

	# Formation rate coefficient of H2 on dust grain, at solar metallicity
    const Rsun = 3.5e-17u"cm^3 * s^-1"

	# Effective metallicity
	const Zeff = 1e-3 * Zsun

	# Clumping factor
	const Cρ = 100.0

	# Condensation constant
	const C_cond = (2 * Rsun * Cρ) / (m_u * Zsun)

	# Our choice of units for C_cond
	const cc_unit = m_u^-1 * t_u^-1 * l_u^3

	# Nominal value of C_cond in our choice of units
	const c_cond = ustrip(Unitful.NoUnits, C_cond / cc_unit)

	# Condensation time-scale
	τ_cond(ρc, Z, fg) = 1.0 / (c_cond * (Z + Zeff) * ρc * fg)
end;

# ╔═╡ 6cab6cb7-a432-40b6-9390-ad0083fe486d
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Dust destruction

Following [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635), we summarize the various physical processes responsible for dust destruction using a single effective timescale, $\tau_\mathrm{dd}$. The main mechanisms contributing to dust destruction include sputtering, collisions with cosmic rays, shocks from supernova (SN) explosions, and radiative torques in the presence of intense, anisotropic radiation fields. In this work, we neglect astration, i.e., the incorporation of dust into stars during star formation.

From [Slavin2015](https://doi.org/10.1088/0004-637X/803/1/7), the dust destruction timescale depends on both the grain composition and the characteristics of the surrounding interstellar medium (ISM). Table 3 from that study provides estimates for carbonaceous and silicate grains under two ISM models: hot ionized medium (HIM) dominated and warm medium (WM) dominated. The HIM scenario assumes dust destruction occurs within warm clouds embedded in a hot ambient medium, while the WM model reflects the fiducial assumptions in their simulation.

| Grain Type          | ISM Model     | SN Mass Interval: Local         | SN Mass Interval: Global        |
|:-------------------:|:-------------:|:-------------------------------:|:-------------------------------:|
| Carbonaceous grains | HIM dominated | $1.6 \pm 0.7 \, \mathrm{Gyr}$   | $1.2 \pm 0.3 \, \mathrm{Gyr}$   |
|                     | WM dominated  | $3.2 \pm 1.4 \, \mathrm{Gyr}$   | $2.6 \pm 0.7 \, \mathrm{Gyr}$   |
| Silicate grains     | HIM dominated | $0.92 \pm 0.39 \, \mathrm{Gyr}$ | $0.72 \pm 0.2 \, \mathrm{Gyr}$  |
|                     | WM dominated  | $2.0 \pm 0.8 \, \mathrm{Gyr}$   | $1.5 \pm 0.4 \, \mathrm{Gyr}$   |

In our model, we will adopt the WM-dominated scenario and the local SN mass interval as our reference, consistent with the fiducial setup in [Slavin2015](https://doi.org/10.1088/0004-637X/803/1/7). This yields the following relevant timescales

| Grain Composition     | $\tau_\mathrm{dd} \, [\mathrm{Gyr}]$ |
|:---------------------:|:------------------------------------:|
| Carbonaceous grains   | $3.2 \pm 1.4$                        |
| Silicate grains       | $2.0 \pm 0.8$                        |

Since we do not distinguish between grain types in our model, we adopt a single dust destruction timescale based on the weighted average (by their uncertanties) of the values above,

$\begin{equation}
	\tau_\mathrm{dd} = 2.3 \, \mathrm{Gyr} \, .
\end{equation}$

This approximation allows us to model dust destruction consistently while minimizing complexity, aligning with the assumptions in previous literature.
"""
  ╠═╡ =#

# ╔═╡ 38e91d0d-1020-44a8-becc-8542fd600104
begin
	# Dust destruction time-scales from Slavin et al. (2015)
	const times = [
		3.2 ± 1.4, # Carbonaceous grains
		2.0 ± 0.8, # Silicate grains
	]

	# Mean dust destruction time-scale
	const Τ_dd = Measurements.value(weightedmean(times)) * u"Gyr"

	# Nominal value of Τ_dd in our choice of units
	const τ_dd = ustrip(t_u, Τ_dd)

	# Nominal value of the inverse of Τ_dd in our choice of units
	const inv_τ_dd = ustrip(t_u^-1, 1.0 / Τ_dd)
end;

# ╔═╡ 4b81f302-9637-4161-b745-8ac39b9e31d3
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Dust creation

Following [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635), we model dust formation as the accretion of metals from the gas phase onto existing dust grains. Based on the formulations in [Dwek1998](https://doi.org/10.1086/305829) (eq. 32) and [Hirashita1999](https://doi.org/10.48550/arXiv.astro-ph/9903259) (eq. 2), the rate of change of the dust mass fraction due to accretion is given by

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_d(t)\right|_{\mathrm{accretion}} = \frac{f_Z}{f_Z + f_d} \, \frac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m) = \frac{f_d}{\tau_\mathrm{dc}} \, ,
\end{equation}$

where $\tau_\mathrm{dg}$ is the dust growth timescale and $\tau_\mathrm{dc}$ is the dust creation time scale, both defined by this expression. The term $f_a + f_m$ represents the neutral gas fraction, and is chosen here as a proxy for the cold gas mass fraction $X_\mathrm{cold}$, consistent with the definition in [Hirashita1999](https://doi.org/10.48550/arXiv.astro-ph/9903259), eq. 29.

The formulation used here differs slightly from that in [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635), where the denominator includes only $f_Z$. By instead using $f_Z + f_d$, we explicitly account for the total mass fraction of metals -- both in gas and dust form. This choice is consistent with interpretations in [Dwek1998](https://doi.org/10.1086/305829), [Hirashita1998](https://doi.org/10.1086/311806), and [Hirashita2018](https://doi.org/10.1093/mnras/sty2838), where the ratio

$\begin{equation}
	\frac{f_Z}{f_Z + f_d} = 1 - \frac{f_d}{f_Z + f_d} \, ,
\end{equation}$

represents the fraction of metals that remain in the gas phase (i.e., not yet incorporated into dust).

The timescale $\tau_\mathrm{dg}$ can be estimated using the expression from [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) (eqs. 19 and 20)

$\begin{equation}
	\tau_\mathrm{dg} = c_a \, A \, \frac{a}{a^0} \, \frac{Z_\odot^d}{Z} \, \frac{n_\mathrm{H}^0}{n_\mathrm{H}} \, \sqrt{\frac{T^0}{T}} \, \frac{S^0}{S} \, ,
\end{equation}$

where

| Parameter        | Description                         |
|:----------------:|:-----------------------------------:|
| $c_a$            | Geometric and integral factor       |
| $A$              | Normalization constant              |
| $a$              | Mean grain radius                   |
| $Z$              | Metallicity                         |
| $n_\mathrm{H}$   | Hydrogen number density             |
| $T$              | Gas temperature                     |
| $S$              | Sticking coefficient                |
| $a^0$            | Fiducial grain radius               |
| $Z_\odot^d$      | Fiducial solar metallicity          |
| $n_\mathrm{H}^0$ | Fiducial hydrogen number density    |
| $T^0$            | Fiducial gas temperature            |
| $S^0$            | Fiducial sticking coefficient       |

Among these parameters, only $Z$ and $n_\mathrm{H}$ vary across simulation cells; all others are treated as fixed constants in our model.
"""
  ╠═╡ =#

# ╔═╡ f525b379-b5e2-48bc-af1a-e57436f6176e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### $c_a$ -- Geometric and integral factor

In contrast to eqs. 19 and 20 from [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x), our formulation includes an additional geometric factor $c_a$. To understand its origin, we need to translate the notation from the Hirashita et al. papers into our own.

From [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x), eq. 30 we can deduce that

$\begin{equation}
	\tau_\mathrm{grow} \equiv \tau_\mathrm{dg} \, ,
\end{equation}$

Then, eq. 33 gives

$\begin{equation}
	\tau_\mathrm{dg} \equiv \tau_\mathrm{grow} \approx \frac{\langle a^3 \rangle_0}{3 \, \langle a^2 \rangle_0 \, a_0} \, \tau \, ,
\end{equation}$

where $\langle a^3 \rangle_0$, $\langle a^2 \rangle_0$, and $a_0$ are the moments of $a$ and can be taken from Table 2, and $\tau$ is a microscopic timescale from eqs. 23 and 24 of the same paper (though we adopt constants from [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x), as explained below).

The $\langle a \rangle$ terms arise because Hirashita distinguishes between

*  $\tau$, a microscopic (local) grain-growth timescale for a single dust grain, and
*  $\tau_\mathrm{grow}$, the macroscopic timescale for growth across an entire cloud.

These are then related by integrals over the grain-size distribution (see eq. 26 in [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x)).

We adopt model C and F models, which are equivalent for our purposes. Their characteristic values are

| Parameter                                       | Value                      |
|:-----------------------------------------------:|:--------------------------:|
| $\langle a^3 \rangle_0 / \langle a^2 \rangle_0$ | $0.0159 \, \mathrm{\mu m}$ |
| $a_0$                                           | $0.1 \, \mathrm{\mu m}$    |

So we define

$\begin{equation}
	c_a = \frac{\langle a^3 \rangle_0}{3 \, \langle a^2 \rangle_0 \, a_0} = 0.053 \, ,
\end{equation}$

The factor of 3 arises when converting from radial growth to mass growth, as shown in [Hirashita2018](https://doi.org/10.1093/mnras/sty2838), eqs. 12 and 13.

Letting

$\begin{align}
	\frac{\mathrm{d} a}{\mathrm{d} t} &= \xi \, \frac{a}{\tau^a} \, , \\
	\frac{\mathrm{d} m}{\mathrm{d} t} &= \xi \, \frac{m}{\tau^m} \, ,
\end{align}$

where $a$ is the grain radius, $m$ its mass, and $\xi$ a dimensionless proportionality factor. Assuming spherical grains with uniform density $s$

$\begin{equation}
	m = \frac{4}{3} \, \pi \, a^3 \, s \, ,
\end{equation}$

we find

$\begin{equation}
	\tau^m = \frac{\tau^a}{3} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ c01cf078-c0cb-4000-bf76-a90cbeb10b89
# Geometric and integral factor
const ca = 0.053;

# ╔═╡ 46c5bd5a-21b1-4f92-8eb1-a11ec2a0c94a
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### $A$ -- Normalization constant

The parameter $A$ depends on the physical properties of the dust grains (composition, size distribution, etc.). In principle, it is determined by microphysics, but in practice its value also depends on the fiducial parameters used in the expression for $\tau_\mathrm{dg}$.

From [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x), the fiducial parameters are

| Parameter        | Fiducial value             |
|:----------------:|:--------------------------:|
| $a^0$            | $0.1 \, \mathrm{\mu m}$    |
| $Z_\odot^d$      | $0.015$                    |
| $n_\mathrm{H}^0$ | $1000 \, \mathrm{cm^{-3}}$ |
| $T^0$            | $50 \, \mathrm{K}$         |
| $S^0$            | $0.3$                      |

Resulting in $A$ values

| Composition | $A$                               |
|:-----------:|:---------------------------------:|
| Silicate    | $6.30 \times 10^7 \, \mathrm{yr}$ |
| Graphite    | $5.59 \times 10^7 \, \mathrm{yr}$ |

Instead, in [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) and [Hirashita2018](https://doi.org/10.1093/mnras/sty2838) we have

| Parameter        | Fiducial value             |
|:----------------:|:--------------------------:|
| $a^0$            | $0.1 \, \mathrm{\mu m}$    |
| $Z_\odot^d$      | $0.02$                     |
| $n_\mathrm{H}^0$ | $1000 \, \mathrm{cm^{-3}}$ |
| $T^0$            | $10 \, \mathrm{K}$         |
| $S^0$            | $0.3$                      |

With corresponding $A$ values

| Composition | $A$                                |
|:-----------:|:----------------------------------:|
| Silicate    | $1.61 \times 10^8 \, \mathrm{yr}$  |
| Graphite    | $0.993 \times 10^8 \, \mathrm{yr}$ |

While one might expect the differences in $A$ to be due solely to the choice of fiducial values, comparing the results reveals discrepancies beyond parameter scaling -- indicating differences in how $A$ itself is computed in each paper.

We will adopt the formulation and fiducial values from [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) and [Hirashita2018](https://doi.org/10.1093/mnras/sty2838), and use the mean of the two compositions

| Parameter        | Fiducial value                   |
|:----------------:|:--------------------------------:|
| $A$              | $1.3 \times 10^8 \, \mathrm{yr}$ |
| $a^0$            | $0.1 \, \mathrm{\mu m}$          |
| $Z_\odot^d$      | $0.02$                           |
| $n_\mathrm{H}^0$ | $1000 \, \mathrm{cm^{-3}}$       |
| $T^0$            | $10 \, \mathrm{K}$               |
| $S^0$            | $0.3$                            |
"""
  ╠═╡ =#

# ╔═╡ 9c5b30d5-08aa-48a8-9ae2-c3b6c432ab89
# Normalization constant
const A = 1.3e8u"yr";

# ╔═╡ ce2383dc-c34f-4b56-8501-42ce0539c95c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### $a$ -- Mean grain radius

The mean grain radius $\langle a \rangle$ is set by the adopted grain-size distribution. Following [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x), and again using models C and F, we have

$\begin{equation}
	\langle a \rangle =  0.00167 \, \mathrm{\mu m} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ a3d1e1bf-c513-4d6b-a43b-3dab0106f1a5
begin
	# Mean grain radius
	const a = 0.00167u"μm"

	# Fiducial mean grain radius
	const a0 = 0.1u"μm"
end;

# ╔═╡ 5680e62b-973b-4a60-bb3f-8785ce07e581
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### $n_\mathrm{H}$ -- Hydrogen number density

The hydrogen number density depends on the phase of the gas surrounding the dust-forming regions (e.g., molecular, atomic, or a mixture of both). For consistency with the previous assumptions, we consider the neutral fraction of the gas

$\begin{equation}
    n_\mathrm{H} = n_a + 2 \, n_m = f_a \, \frac{\rho_\mathrm{cell}}{m_p} + 2 \, f_m \, \frac{\rho_\mathrm{cell}}{2 \, m_p} = f_n \, \frac{\rho_\mathrm{cell}}{m_p} \, ,
\end{equation}$

where $m_p$ is the proton mass, and $f_n = f_a + f_m$ represents the total neutral fraction of hydrogen (atomic + molecular).
"""
  ╠═╡ =#

# ╔═╡ 47bf94da-1368-41bd-ba46-c9c1f75cf44e
# Fiducial hydrogen number density
const nH0 = 1000.0u"cm^-3";

# ╔═╡ 9af6ce74-5678-4b48-8758-37b0d5a6f0e4
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### $T$ -- Gas temperature

As with $n_\mathrm{H}$, the gas temperature depends on the phase. For consistency with the assumption of a neutral medium, we adopt

$\begin{equation}
    T = 100 \, \mathrm{K} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ a9c6a292-086e-4aa0-9856-78d3ab3fbe35
begin
	# Gas temperature
	const T = 100.0u"K"

	# Fiducial gas temperature
	const T0 = 10.0u"K"
end;

# ╔═╡ 5890b699-b7de-47a3-bee7-1e7dd7663fbe
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### $S$ -- Sticking coefficient

For the sticking coefficient we adopt the fit of the results in [Leitch-Devlin1985](https://doi.org/10.1093/mnras/213.2.295) presented by [Grassi2011](https://doi.org/10.1051/0004-6361/200913779), which reads

$\begin{equation}
    S(T_g, T_d) = 0.019 \, T_g \, (0.0017 \, T_d + 0.4) \, \exp(-0.007 \, T_g) \, .
\end{equation}$

where $T_g$ and $T_d$ are the gas and dust temperatures, respectively, in Kelvin.

Assuming the simplest case with $T_g = T_d = 100 \, \mathrm{K}$, we obtain

$\begin{equation}
    S = 0.538 \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ ef79392b-395c-497a-9c0e-dc2cd468f6e1
begin
	# Fit for the sticking coefficient from Grassi et al. (2011)
    Se(Tg, Td) = 0.019 * Tg * (0.0017 * Td + 0.4) * exp(-0.007 * Tg)

	# Sticking coefficient at 100 K
	const S = Se(ustrip(u"K", T), ustrip(u"K", T))

	# Fiducial sticking coefficient
	const S0 = 0.3
end;

# ╔═╡ 8c9ab125-2acb-4732-a9bf-7838e819e4f7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Final form of $\tau_\mathrm{dg}$ and $\tau_\mathrm{dc}$

Using all the components derived above, the dust growth time scale

$\begin{equation}
	\tau_\mathrm{dg} = c_a \, A \, \frac{a}{a^0} \, \frac{Z_\odot^d}{Z} \, \frac{n_\mathrm{H}^0}{n_\mathrm{H}} \, \sqrt{\frac{T^0}{T}} \, \frac{S^0}{S} \, ,
\end{equation}$

can be rearranged into the more useful form

$\begin{equation}
	\tau_\mathrm{dg} = \frac{1}{C_\mathrm{dg} \, Z \, \rho_\mathrm{cell} \, f_n} \, ,
\end{equation}$

where the coefficient $C_\mathrm{dg}$ encapsulates all constants and fiducial parameters

$\begin{align}
	C_\mathrm{dg} &= \frac{a^0 \, \sqrt{T} \, S}{c_a \, A \, a \, Z_\odot^d \, n_\mathrm{H}^0 \, m_p \, \sqrt{T^0} \, S^0} \\
    &= 2.463 \, \mathrm{mp^{-1} \, Myr^{-1} \, cm^3} \, .
\end{align}$

Note that $Z_\odot^d$ here is the dust-phase solar metallicity used in [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) and [Hirashita2018](https://doi.org/10.1093/mnras/sty2838), and not the solar metallicity used elsewhere in our framework.

The dust creation time scale is then

$\begin{align}
	\tau_\mathrm{dc} &= \frac{f_Z + f_d}{f_Z \, f_n} \, \tau_\mathrm{dg} \\
	&= \frac{f_Z + f_d}{f_Z \, f_n} \, \frac{1}{C_\mathrm{dg} \, Z \, \rho_\mathrm{cell} \, f_n} \\
	&= \frac{1}{C_\mathrm{dg} \, f_Z \, f_n^2 \, \rho_\mathrm{cell}} \, ,
\end{align}$

where, following [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x), we took the total metallicity as $Z = f_Z + f_d$.

### Dust accretion term

The dust growth source term becomes

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_d(t)\right|_{\text{accretion}} = \frac{f_d}{\tau_\mathrm{dc}} = C_\mathrm{dg} \, f_Z \, f_d \, f_n^2 \, \rho_\mathrm{cell} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ e74e6af3-7c28-4836-a7a7-db65cca302dc
# Alternative factor for a different way of writing τdg
const κ_d = uconvert(u"Myr", (ca * A * a * sqrt(T0) * S0) / (a0 * sqrt(T) * S))

# ╔═╡ 2a39d6f8-da49-4976-9aa7-889391e55a5d
begin
	# Dust-phase solar metallicity
	const Zdsun = 0.02

	# Dust growth constant
	const C_dg = (a0 * sqrt(T) * S) / (ca * A * a * Zdsun * nH0 * Unitful.mp * sqrt(T0) * S0)

	# Our choice of units for C_dg
	const cdg_unit = m_u^-1 * t_u^-1 * l_u^3

	# Nominal value of C_dg in our choice of units
	const c_dg = ustrip(Unitful.NoUnits, C_dg / cdg_unit)

	# Dust growth time-scale
	τ_dg(ρc, fn, Z) = 1 / (c_dg * Z * ρc * fn)
end;

# ╔═╡ f4399d4f-9e85-4974-9f39-d47acaf1399c
md"""
## Initial Mass Function (IMF)

The initial mass function $\phi(m)$ describes the number of stars formed per unit mass interval. Specifically, $\phi(m) \, \mathrm{d}m$ gives the number of stars with initial masses in the range $[m, m + \mathrm{d}m]$.

To represent a stellar population with total mass $M$, the IMF must satisfy the normalization condition

$\begin{equation}
    M = \int_{m_\mathrm{low}}^{m_\mathrm{high}} m \, \phi(m) \, \mathrm{d}m \, ,
\end{equation}$

One of the simplest and most widely used forms of the IMF is a power-law

$\begin{equation}
    \phi(m) = A \, m^{-\alpha}\, ,
\end{equation}$

where $\alpha$ is the slope and $A$ is a normalization constant. Stellar masses are typically expressed in solar units, i.e., $[m] = M_\odot$.

In what follows we will use the IMF of [Chabrier2003](https://doi.org/10.1086/374879)
"""

# ╔═╡ 6a59ed2e-f040-4e10-bc55-91d2c1dcc97e
#####################################################################################
# Initial mass functions
# 
# The mass range [m_low, m_high] will depend on the application, 
# but a general mass range for every model is 0.1 <= m / M⊙ <= 100.0
# 
# In the following implementations, a specific normalization is not assumed, when
# possible, A = 1 is used for simplicity. To model a population with a specific total
# mass M, the IMF must be rescaled using the appropriate normalization constant
# derived from the integral above
# 
# Sources:
# 
#     Original papers
#     eqs. 6 to 14 in Mollá et al. (2015)
#     eqs. A1 to A4 in Appendix A of Millán-Irigoyen et al. (2021)
#####################################################################################
begin
    # Salpeter (1955) eq. 5
	# https://doi.org/10.1086/145971
    ϕSAL(m::Float64)::Float64 = m^(-2.35)

    # Miller et al. (1979) eq. 30 and Table 7
	# https://doi.org/10.1086/190629
    function ϕMIL(m::Float64)::Float64 
		C1 = 1.09
    	C2 = -1.02
		logm = log10(m)
		return (1 / m) * exp(-C1 * (logm - C2)^2)
	end

    # Ferrini et al. (1992) eq. 2
	# https://doi.org/10.1086/171066
    function ϕFER(m::Float64)::Float64
		logm = log10(m)
        return m^(-0.52) * exp10(-sqrt(2.07 * logm^2 + 1.92 * logm + 0.73))
    end

    # Kroupa et al. (1993) eq. 13
	# https://doi.org/10.1093/mnras/262.3.545
    function ϕKRO_93(m::Float64)::Float64
        if m < 0.5
            return m^(-1.3)
        elseif 0.5 <= m < 1
            return 0.5 * (m^-2.2)
        else
            return 0.5 * (m^-2.7)
        end
    end

	# Kroupa (2002)	eq. 4 and 5 in Table 1
	# https://doi.org/10.1126/science.1067524
	function ϕKRO_02(m::Float64)::Float64
        if m < 0.08
            return m^(-0.3)
		elseif 0.08 <= m < 0.5
			return 0.08 * m^(-1.3)
        elseif 0.5 <= m < 1
            return 0.04 * m^(-2.3)
        else
            return 0.04 * m^(-2.7)
        end
    end

    # Chabrier (2003) eq. 2
	# https://doi.org/10.1086/374879
    function ϕCHA(m::Float64)::Float64
        if m <= 1
			logm = log10(m)
            return (1 / m) * exp(-(logm + 0.6576)^2 / 0.6498)
        else
            return 0.514 * m^(-2.3)
        end
    end

    # Millán-Irigoyen et al. (2020) eq. 14
	# https://doi.org/10.1093/mnras/staa635
    function ϕMILLA(m::Float64)::Float64
        if m < 0.5
            return m^(-1.3)
        else
            return 0.5 * m^(-2.3)
        end
    end

    ϕSAL(m::Quantity)::Float64    = ϕSAL(ustrip(u"Msun", m))
    ϕMIL(m::Quantity)::Float64    = ϕMIL(ustrip(u"Msun", m))
    ϕFER(m::Quantity)::Float64    = ϕFER(ustrip(u"Msun", m))
    ϕKRO_93(m::Quantity)::Float64 = ϕKRO_93(ustrip(u"Msun", m))
    ϕKRO_02(m::Quantity)::Float64 = ϕKRO_02(ustrip(u"Msun", m))
    ϕCHA(m::Quantity)::Float64    = ϕCHA(ustrip(u"Msun", m))
    ϕMILLA(m::Quantity)::Float64  = ϕMILLA(ustrip(u"Msun", m))

    const IMF_FUNCTIONS = Dict(
        "Salpeter1955"        => ["SAL", ϕSAL],
        "Miller1979"          => ["MIL", ϕMIL],
        "Ferrini1990"         => ["FER", ϕFER],
        "Kroupa1993"          => ["KRO_93", ϕKRO_93],
        "Kroupa2002"          => ["KRO_02", ϕKRO_02],
        "Chabrier2003"        => ["CHA", ϕCHA],
        "Millan-Irigoyen2020" => ["MILLA", ϕMILLA],
    )

	# Chosen initial mass function
	const IMF = "Chabrier2003"
	
	# Initial mass functions
	ϕ(m) = IMF_FUNCTIONS[IMF][2](m) * u"Msun^-1"
	mϕ(m) = m * ϕ(m)

	# Lower mass limit for the IMF
	const M_LOW = 0.1u"Msun"    
	
	# Upper mass limit for the IMF
	const M_HIGH = 100.0u"Msun"          

	# IMF normalization
	const IMF_NORM = quadgk(mϕ, M_LOW, M_HIGH)[1]
end;

# ╔═╡ 14ae7f11-1065-46aa-a6ed-eea400c0d2ec
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Photodissociation and photoionization

We define the disassociated mass rate per unit of created stellar mass as

$\begin{equation}
    \eta_\mathrm{diss}^* = \frac{\dot{M}_\mathrm{diss}}{\mathrm{SFR}}  \, ,
\end{equation}$

and the ionized mass rate per unit of created stellar mass as

$\begin{equation}
    \eta_\mathrm{ion}^* = \frac{\dot{M}_\mathrm{ion}}{\mathrm{SFR}}  \, .
\end{equation}$

The mass rates $\dot{M}_\mathrm{diss}$ and $\dot{M}_\mathrm{ion}$ can be computed from the photon production rates as

$\begin{align}
    \dot{M}_\mathrm{diss} &= c_\mathrm{diss} \, \left(1 - e^{-\tau_\mathrm{diss}}\right) \, \dot{N}_\mathrm{diss} \, , \\
    \dot{M}_\mathrm{ion} &= c_\mathrm{ion} \, \left(1 - e^{-\tau_\mathrm{ion}}\right) \, \dot{N}_\mathrm{ion} \, ,
\end{align}$

where $\dot{N}_\mathrm{diss}$ is the number of photodissociating photons produced per unit time (in the Lyman-Werner band, $911.6 \, \mathrm{Å}$ to $1107 \, \mathrm{Å}$), $\dot{N}_\mathrm{ion}$ is the number of ionizing photons produced per unit time (between $0$ and $911.6 \, \mathrm{Å}$), $c_\mathrm{diss}$ and $c_\mathrm{ion}$ are mass conversion factors per photon, and $\tau_\mathrm{diss}$ and $\tau_\mathrm{ion}$ are the optical depths for dissociating and ionizing radiation, respectively.

The exponential attenuation terms, $1 - e^{-\tau}$, represent the probability that a photon is absorbed in the cell and thus can contribute to dissociation or ionization event. In an optically thick regime, nearly all photons interact, whereas in an optically thin medium, most photons escape without interaction.

The optical depth $\tau$ depends on the number density of the absorbing species ($n$), the interaction cross-section ($\sigma$), and the characteristic path length through the gas cell ($h$)

$\begin{equation}
    \tau = n \, \sigma \, h \, .
\end{equation}$

The optical depth for ionizing photons is given by

$\begin{equation}
    \tau_\mathrm{ion} = \left(\frac{f_a \, \rho_\mathrm{cell}}{m_p}\right) \sigma_\mathrm{ion} \, h \, ,
\end{equation}$

where $\sigma_\mathrm{ion}$ is the hydrogen photoionization cross-section at the Lyman limit, given by (Table 2 of [Hirashita2002](https://doi.org/10.1046/j.1365-8711.2002.05968.x))

$\begin{equation}
    \sigma_\mathrm{ion} = 6.3 \times 10^{-18} \, \mathrm{cm^2} \, .
\end{equation}$

Photodissociation by Lyman-Werner photons is approximated using an effective cross-section $\sigma_\mathrm{diss,eff}$ (eq. 22 in [Nickerson2018](https://doi.org/10.1093/mnras/sty1556))

$\begin{equation}
    \sigma_\mathrm{diss,eff} = 2.1 \times 10^{-19} \, \mathrm{cm^2} \, .
\end{equation}$

The corresponding optical depth is

$\begin{equation}
    \tau_\mathrm{diss} = \left(\frac{f_m \, \rho_\mathrm{cell}}{2 \, m_p}\right) \sigma_\mathrm{diss,eff} \, h \, .
\end{equation}$

where the factor of 2 accounts for the two protons in each $\mathrm{H}_2$ molecule.

Photodissociation of molecular hydrogen by Lyman-Werner photons occurs via the two-step Solomon process ([Draine1996](http://doi.org/10.1086/177689), [Incatasciato2023](https://doi.org/10.1093/mnras/stad1008))

1.  **Absorption**: An $\mathrm{H}_2$ molecule absorbs a UV photon in the Lyman-Werner bands (absorption lines in the energy range of $11.2 \, \mathrm{eV} < E < 13.6 \, \mathrm{eV}$). This excites the molecule to a higher electronic state.

$\begin{equation}
	\mathrm{H}_2 + \gamma_{LW} \rightarrow \mathrm{H}_2^* \, .
\end{equation}$

2.  **Decay**: The excited molecule ($\mathrm{H}_2^*$) quickly decays. About $85\%$ of the time, it decays back to a stable vibrational state of the ground electronic level, emitting a photon. However, about $15\%$ of the time, it decays to the vibrational continuum of the ground state, meaning the molecule dissociates into two hydrogen atoms.

$\begin{align}
	&\mathrm{H}_2^* \rightarrow \mathrm{H} + \mathrm{H} \,\, \mathrm{(dissociation, \,\, \sim 15\% \,\, probability)} \\
	&\mathrm{H}_2^* \rightarrow \mathrm{H}_2 + \gamma \,\, \mathrm{(radiative \,\, decay, \,\, \sim 85\% \,\, probability)}
\end{align}$

Finally, the unit conversion coefficients are

$\begin{align}
    c_\mathrm{diss} &= 0.3 \, \mathrm{mp} \, , \\
    c_\mathrm{ion} &= 1.0 \, \mathrm{mp} \, .
\end{align}$

These reflect mass affected per photodissociation or ionization event.
"""
  ╠═╡ =#

# ╔═╡ c8efcb3f-cc6f-4c45-a0f3-55cdb73d5195
begin
	# Cross sections
	const σ_ion = 6.3e-18u"cm^2"
	const σ_diss = 2.1e-19u"cm^2"

	# Optical depth constant
	const C_τion = σ_ion / u"mp"
	const C_τdiss = σ_diss / 2.0u"mp"

	# Our choice of units for C_τion and C_τdiss
	const cτ_unit = l_u^2 * m_u^-1

	# Nominal value of C_τion in our choice of units
	const c_τion = ustrip(Unitful.NoUnits, C_τion / cτ_unit)

	# Nominal value of C_τdiss in our choice of units
	const c_τdiss = ustrip(Unitful.NoUnits, C_τdiss / cτ_unit)

	# Solomon factor 
	const solomon = 0.15
	
	# Mass conversion factors
	const c_diss = solomon * 2.0u"mp"
	const c_ion  = 1.0u"mp"

	# Solar luminosity
	# Normalization factor for the spectral energy distribution (SED)
	const u_sed = 3.82e33u"erg * s^-1 * Å^-1 * Msun^-1"

	# Unit for Q_diss and Q_ion
	const u_q = u"Msun^-1 * s^-1"

	# Wavelength range for the photodissociation of hydrogen molecules
    const λ_diss = (911.6u"Å", 1107.0u"Å")

	# Wavelength range for the photoionization of hydrogen atoms
    const λ_ion = (0.0u"Å", 911.6u"Å")

	# Computed magnitudes table
	# HR-pyPopStar (2025)
    # https://www.fractal-es.com/PopStar/
	const MAGNITUDE_TABLE = "./data/HR-pyPopStar_2025_magnitudes.dat"

	# Spectral energy distribution (SED)
    # https://www.fractal-es.com/PopStar/
	@static if Sys.iswindows()
		const HR_LUMINOSITY = true
	elseif Sys.islinux()
		const HR_LUMINOSITY = false
	end

	 @kwdef struct LuminosityTable
		# Path to the files with the SED data
		sed_path::String
		# Column names in the tables, it has to have :λ and :Ltot
		header::Vector{Symbol}
		# First row with data in the tables 
		skipto::Int64
		# Column separator in the tables 
		delim::String
		# Regex to extract the IMF, Z, and stellar age from the file names
		regex::Regex
		# String to add at the beginning of the second capture (metallicity). 
		str_add::String
	end
	
	if HR_LUMINOSITY
		
		#############################################################################
		# High resolution SED
		# Millán-Irigoyen et al. (2021)
		# https://doi.org/10.1093/mnras/stab1969
		#############################################################################
		
		luminosity_table = LuminosityTable(
			sed_path="./data/luminosity/hr/total/",
			header=[:λ, :Ltot],
			skipto=2,
			delim="\t",
			# SSP-IMF-total_ZXXXX_logtXXXX.dat
			regex=r"SSP-(\w+)-total_Z([0-9.]+)_logt([0-9.]+)\.dat",
			str_add=""
		)

	else

		#############################################################################
		# Low resolution SED
	    # Mollá et al. (2009)
	    # https://doi.org/10.1111/j.1365-2966.2009.15160.x
		#############################################################################

		luminosity_table = LuminosityTable(
			sed_path="./data/luminosity/lr/",
			header=[:λ, :Lstar, :Lneb, :Ltot],
			skipto=1,
			delim=" ",
			# spneb_IMF_xxx_xxx_zXXXX_tXXXX
			regex=r"spneb_([A-Za-z]+)_[\d.]+_[\d.]+_z(\d+)_t(\d+(?:\.\d+)?)",
			str_add="0."
		)
		
	end
end;

# ╔═╡ 4d09ed45-423c-4bd6-802d-59389a966d2e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The number of photons can be computed from $Q(t', Z)$ which is defined as the number of photons produced by a stellar population of one solar mass, of age $t'$, and metallicity $Z$, per unit of time. So, we have the relation

$\begin{equation}
    \dot{N}(t) = \int_0^t \mathrm{SFR}(t - t') \, Q(t', Z(t - t')) \, \mathrm{d}t' \, ,
\end{equation}$

where $\dot{N}(t)$ is the current (time $t$) rate of photon production, $\mathrm{SFR}(t - t')$ is the instantaneous SFR at the moment of birth of the stellar population of current age $t'$ (in units of solar mass per unit of time), and $Z(t - t')$ is defined as the metallicity at that moment. Because most of the contribution to the integral comes from young blue stars that die in the first $10 \ \mathrm{to} \ 100 \, \mathrm{Myr}$, it is possible to approximate

$\begin{align}
    \mathrm{SFR}(t - t') &\approx \mathrm{SFR}(t) \, , \\
	Z(t - t') &\approx Z(t) \, ,
\end{align}$

where we are assuming that the $\mathrm{SFR}$ and $Z$ stayed approximately constant since the younger stars where born (time $t - t'$) until today (time $t$).

So, if we define $\eta$ as 

$\begin{equation}
    \eta = \eta^* \, (1 - e^{-\tau})^{-1} \, .
\end{equation}$

We end up with the general form

$\begin{equation}
    \eta = c \, \frac{\dot{N}}{\mathrm{SFR}} = c \, \int_0^t Q(t', Z) \, \mathrm{d}t' \, ,
\end{equation}$

The value of $Q$ can be calculated using

$\begin{equation}
    Q(t', Z) = \int_{\lambda_1}^{\lambda_2} \frac{L_\lambda(t', Z)}{E_\lambda} \mathrm{d}\lambda = \int_{\lambda_1}^{\lambda_2} \frac{\lambda \, L_\lambda(t', Z)}{h \, c} \, \mathrm{d}\lambda \, ,
\end{equation}$

where $L_\lambda(t', Z)$ is the luminosity per unit of wavelength of a stellar population of one solar mass, age $t'$, and metallicity $Z$; and $E_\lambda = \lambda / (h \, c)$ is the energy of a photon of wavelength $\lambda$. So, integrating between the wavelength of interest, we get the number of photons produced per unit of time for a stellar population of the given characteristics.

The luminosity will not only depend on the age and metallicity of the population, but on the initial mass function (IMF) too, so in principle for each IMF, $t'$, and $Z$ we have a function of luminosity versus $\lambda$.

Using the values from $\texttt{HR-pypopstar}$ ([Millan-Irigoyen2021](https://doi.org/10.1093/mnras/stab1969)) we compute a table of $Q$ for four different IMFs ([Salpeter1955](https://doi.org/10.1086/145971), [Kroupa2001](https://doi.org/10.1046/j.1365-8711.2001.04022.x), [Ferrini1990](https://ui.adsabs.harvard.edu/abs/1990A%26A...231..391F), and [Chabrier2003](https://doi.org/10.1086/374879)), four metallicities ($0.004$, $0.008$, $0.019$, and $0.05$) and more than a hundred ages ranging from $0.1 \, \mathrm{Myr}$ to $13.8 \, \mathrm{Gyr}$.
"""
  ╠═╡ =#

# ╔═╡ b9d7da8b-5f01-4f4f-8929-fd44f9689f3e
##################################################################################
# Compute Q_diss and Q_ion for a given spectral energy distribution (SED) file
#
# Arguments
#
#     file:             Path to the file with the SED data
#     luminosity_table: Metadata of the luminosity table to be used
#
# Returns
#
#     (Q_diss [Msun^(-1) s^(-1)], Q_ion [Msun^(-1) s^(-1)])
##################################################################################
function integrate_SED(
	file::String,
	luminosity_table::LuminosityTable,
)::NTuple{2,Quantity}

	# Load the data
	data = CSV.read(
		file, 
		DataFrame; 
		delim=luminosity_table.delim,
		ignorerepeated=true,
		header=luminosity_table.header, 
		skipto=luminosity_table.skipto, 
		types=Float64,
	)

	# Wavelength
	λ = data[!, :λ] .* u"Å"

	# Fix for repeted wavelengths
	Interpolations.deduplicate_knots!(λ)

	# Total (stellar + nebular) SED per unit wavelength
	Ltot = data[!, :Ltot] .* u_sed

	# SED interpolation function
	interp_Ltot = linear_interpolation(λ, Ltot, extrapolation_bc=Flat())

	# SED integrand
	integrand(λ) = λ * interp_Ltot(λ) / (Unitful.h * Unitful.c)

	Q_diss = quadgk(integrand, λ_diss[1], λ_diss[2])[1] |> u_q
	Q_ion  = quadgk(integrand, λ_ion[1], λ_ion[2])[1] |> u_q

	return Q_diss, Q_ion

end;

# ╔═╡ e527f834-9965-4cdf-bda8-b4ea07e86f27
##################################################################################
# Compute Q_diss and Q_ion for a given IMF
#
# Arguments
#
#     luminosity_table: Metadata of the luminosity table to be used
#     imf:              Target IMF
#
# Returns
#
#     A dataframe with columns:
# 
#       :Z       -> Stellar metallicity
#       :log_age -> Stellar age as log(age [yr])
#       :Q_diss  -> Number of dissociating photons per unit time
#       :Q_ion   -> Number of ionizing photons per unit time
##################################################################################
function compute_Q(luminosity_table::LuminosityTable, imf::String)::DataFrame

	# List the files with the SEDs
    files = readdir(luminosity_table.sed_path; join=true)

	# Parse the filenames
	matches = match.(luminosity_table.regex, basename.(files))

	# Read only the files for the target IMF
	idxs = findall(i->uppercase(matches[i].captures[1]) == imf, eachindex(files))

	target_files   = files[idxs]
	target_matches = matches[idxs]

	# Allocate memory
	Z = Vector{Float64}(undef, length(target_files))
	log_age = Vector{Float64}(undef, length(target_files))
	Q_diss = Vector{Quantity}(undef, length(target_files))
	Q_ion = Vector{Quantity}(undef, length(target_files))

    for (i, (file, m)) in enumerate(zip(target_files, target_matches))

		# Metallicity
		Z[i] = parse(Float64, luminosity_table.str_add * m.captures[2])

		# Stellar age as log(age [yr])
		log_age[i] = parse(Float64, m.captures[3])

		# Number of dissociating and ionizing photons per unit time
		Q_diss[i], Q_ion[i] = integrate_SED(file, luminosity_table)
		
    end

	Q_df = sort!(identity.(DataFrame(; Z, log_age, Q_diss, Q_ion)), [:Z, :log_age])

	return Q_df
	
end;

# ╔═╡ c869e896-cc87-4b12-bba2-84b0cf0964c6
begin
	Q_df = compute_Q(luminosity_table, IMF_FUNCTIONS[IMF][1])

	# List of stellar ages
	const Q_ages = unique(Q_df[!, :log_age])

	# List of metallicities
    const Q_metals = unique(Q_df[!, :Z])

	########################################################
	# Reshape the dataframe into matrices for interpolation
	# Rows    -> :log_age
	# Columns -> :Z
	########################################################
	wide_diss = unstack(Q_df, :log_age, :Z, :Q_diss)
	sort!(wide_diss, :log_age)
	
	wide_ion = unstack(Q_df, :log_age, :Z, :Q_ion)
	sort!(wide_ion, :log_age)

    Qdiss_mat = ustrip.(u_q, Matrix(wide_diss[:, Not(:log_age)]))
	Qion_mat  = ustrip.(u_q, Matrix(wide_ion[:, Not(:log_age)]))

	# Create the interpolation function for Q_diss
	Qdiss_interp = linear_interpolation(
		(Q_ages, Q_metals),
		Qdiss_mat,
		extrapolation_bc=Flat(),
	)

	# Create the interpolation function for Q_ion
	Qion_interp = linear_interpolation(
		(Q_ages, Q_metals),
		Qion_mat,
		extrapolation_bc=Flat(),
	)

	Q_df
end

# ╔═╡ 2c2f97cc-127a-446e-a6cd-e9651df833f0
##################################################################################
# Compute η_diss and η_ion for a given time and metallicity
#
# Arguments
#
#     age: Stellar age [Myr]
#     Z:   Stellar metallicity [dimensionless]
#
# Returns
#
#     (η_diss [dimensionless], η_ion [dimensionless])
##################################################################################

function compute_η(age::Float64, Z::Float64)::NTuple{2,Float64}

	# Integration limits
	t0 = 0.0u"Myr"
	tf = age * u"Myr"

	# Interpolation function for Q_diss at the given Z
	Q_diss(t) = Qdiss_interp(log10(ustrip(u"yr", t)), Z) * u_q

	# Interpolation function for Q_ion at the given Z
	Q_ion(t) = Qion_interp(log10(ustrip(u"yr", t)), Z) * u_q

	η_diss = uconvert(Unitful.NoUnits, c_diss * quadgk(Q_diss, t0, tf)[1])
	η_ion = uconvert(Unitful.NoUnits, c_ion * quadgk(Q_ion, t0, tf)[1])
	
	return η_diss, η_ion

end;

# ╔═╡ 9a756404-f35b-40b9-bc88-4b7e3a7df0b6
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Shielding

Not all the photons discussed in the previous section reach atoms or molecules due to attenuation by dust and molecular gas. This shielding reduces the effective photodissociation and photoionization rates.

Following [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) and [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x) -- based on earlier models by [Draine1996](https://doi.org/10.1086/177689) and [Glover2007](https://doi.org/10.1086/512238) -- the transmitted fraction of the ionizing and dissociating flux can be expressed as

$\begin{align}
    e_\mathrm{ion} &= S_d \, , \\
    e_\mathrm{diss} &= S_d \, S_{\mathrm{H}_2} \, ,
\end{align}$

where $S_d$ accounts for attenuation by dust, and $S_{\mathrm{H}_2}$ accounts for molecular self-shielding. Both terms depend on the local metallicity and gas column density.

From eq. 8 in [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x), the dust shielding term is given by

$\begin{equation}
    S_d = \exp(-\sigma_\mathrm{d, eff} \, Z^* \, (N_\mathrm{HI} + 2 \, N_{\mathrm{H}_2})) \, ,
\end{equation}$

where $\sigma_\mathrm{d, eff} = 4 \times 10^{-21} \, \mathrm{cm^{-2}}$ is the effective dust cross-section per hydrogen atom, $Z^*$ is the metallicity in solar units, and $N_i$ represents the column density of species $i$. Notice that this parameter is based on the empirical fits from [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) (Section 2.2) which is different from the fits in [Draine1996](https://doi.org/10.1086/177689) and [Glover2007](https://doi.org/10.1086/512238).

To estimate the column density, we use a characteristic length scale $h$, which we take as the smoothing length of the gas cell -- i.e., the force softening length in $\texttt{AREPO}$, following [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x).

Using the relation between number density and fractions, we can write

$\begin{equation}
    N_\mathrm{HI} + 2 \, N_{\mathrm{H}_2} = (f_a + f_m) \, \frac{\rho_\mathrm{cell} \, h}{m_p} \, .
\end{equation}$

Substituting, the dust shielding term becomes

$\begin{equation}
    S_d = \exp\left(-C_{sd} \, (f_Z + f_d) \, (f_a + f_m) \, \rho_\mathrm{cell} \, h\right) \, ,
\end{equation}$

where we used $Z = f_Z + f_d$ and

$\begin{align}
    C_{sd} &= \frac{\sigma_\mathrm{d, eff}}{Z_\odot \, m_p} \\
	&= 3.150 \times 10^{-19} \, \mathrm{mp^{-1} \, cm^2} \, .
\end{align}$
"""
  ╠═╡ =#

# ╔═╡ 805ef471-1bd2-4ce2-bf50-eeb9c6b10e4d
begin
	# Effective cross section
	const σdeff = 4.0e-21u"cm^2"

	# Dust shielding constant
	const C_sd = σdeff / (Zsun * u"mp")

	# Unit of C_sd
	const csd_unit = m_u^-1 * l_u^2

	# Nominal value of C_sd in our choice of units
	const c_sd = ustrip(Unitful.NoUnits, C_sd / csd_unit)
end;

# ╔═╡ 6eef3827-1af6-431d-a163-67bb595c337d
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The molecular shielding factor is given in eq. 9 of [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x) as

$\begin{equation}
    S_{\mathrm{H}_2} = \frac{1 - \omega_{\mathrm{H}_2}}{(1 + x)^2} + \frac{\omega_{\mathrm{H}_2}}{(1 + x)^{1/2}} \, \exp(-0.00085 \, (1 + x)^{1/2}) \, ,
\end{equation}$

where $\omega_{\mathrm{H}_2} = 0.2$ and $x = N_{\mathrm{H}_2} \, / \, 5 \times 10^{14} \, \mathrm{cm^{-2}}$.

Notice that in eq. 37 of [Draine1996](https://doi.org/10.1086/177689) and eq. 32 of [Glover2007](https://doi.org/10.1086/512238) the first term is

$\begin{equation}
    \frac{1 - \omega_{\mathrm{H}_2}}{(1 + x / b_5)^2} \, ,
\end{equation}$

where $b_5 = b / 10^5 \, \mathrm{cm \, s^{-1}}$ and $b$ is the Doppler broadening parameter. It is implicit in [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) and [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x) that they assumed $b = 10^5 \, \mathrm{cm \, s^{-1}}$. We will do the same.

We can rewrite $x$ as

$\begin{equation}
    x = \frac{f_m \, \rho_\mathrm{cell} \, h}{2 \, m_p \, x_\mathrm{nf}} = C_{sh2} \, f_m \, \rho_\mathrm{cell} \, h \, ,
\end{equation}$

where $x_\mathrm{nf} = 5 \times 10^{14} \, \mathrm{cm^{-2}}$ and

$\begin{align}
    C_{sh2} &= \frac{1}{2 \, m_p \, x_\mathrm{nf}} \\
	&= 10^{-15} \, \mathrm{mp^{-1} \, cm^2}\, .
\end{align}$
"""
  ╠═╡ =#

# ╔═╡ e14648ff-69f6-47c2-924c-c5ac69e200ed
begin
	# Molecular self-shielding parameters
	const ωH2 = 0.2
	const exp_fac  = 8.5e-4

	# Fiducial column density
	const xnf = 5.0e14u"cm^-2"

	# Molecular self-shielding constant
	const C_sh2 = 1.0 / (2u"mp" * xnf)

	# Unit of C_sh2
	const ch2_unit = m_u^-1 * l_u^2

	# Nominal value of C_sh2 in our choice of units
	const c_sh2 = ustrip(Unitful.NoUnits, C_sh2 / ch2_unit)
end;

# ╔═╡ aaba63dd-8fc3-405d-ad19-8e8368e70019
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### UVB Photoionization

We account for the ionizing effect of the metagalactic ultraviolet background (UVB) using the ionization rates provided in Table 3 of [Haardt2012](https://doi.org/10.1088/0004-637X/746/2/125).

These rates correspond to the optically thin limit, as discussed in Section 9.1 of [Haardt2012](https://doi.org/10.1088/0004-637X/746/2/125). To incorporate attenuation due to dust, we apply the dust shielding factor $S_d$. Additionally, since the tabulated ionization rates $\Gamma_i$ are defined **per ion** (see their eq. 18), the effective rate must be scaled by the neutral atomic hydrogen fraction $f_a$.

The resulting contribution to the time evolution of the atomic hydrogen fraction is

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_a(t)\right|_{\mathrm{UVB}} = -\left. \frac{\mathrm{d}}{\mathrm{d}t}f_i(t)\right|_{\mathrm{UVB}} = -S_d \, \Gamma_\mathrm{UVB} \, f_a \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 03754250-3831-41db-ac11-fe6f5a08f0c1
begin
	# Path to the UVB table from Haardt et al. (2012)
	const UVB_TABLE_PATH = "./data/photoionization_rate_hm12.txt"
	const UVB_TABLE = CSV.read(UVB_TABLE_PATH, DataFrame; header=[:redshift, :ΓUVB])
end;

# ╔═╡ 2d27f581-2c53-409e-afde-812e70bba8bf
##################################################################################
# Compute ΓUVB
#
# Arguments
#
#     a: Scale factor [dimensionless]
#
# Returns
#
#     ΓUVB [Myr^(-1)]
##################################################################################
function compute_UVB(a::Float64)::Float64

	z = (1 / a) - 1

    ΓUVB = linear_interpolation(
		UVB_TABLE[!, :redshift],
		UVB_TABLE[!, :ΓUVB],
		extrapolation_bc=Flat(),
	)

	return ΓUVB(z)

end;

# ╔═╡ 6ddf393c-9f37-43b2-b401-9ab0e7a9e88c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### LWB photodissociation

The overall rate of photodissociation of molecular hydrogen depends on how many LW photons are available and the probability that a molecule will absorb one.

Then, the photodissosiation rate can be modeled as (eq. 7 in [Abel1997](https://doi.org/10.1016/S1384-1076(97)00010-9) and eq. 3 in [Incatasciato2023](https://doi.org/10.1093/mnras/stad1008))

$\begin{equation}
	\frac{\mathrm{d}}{\mathrm{d}t} \rho_m(t) = -\rho_m(t) \int_{LW} \sigma_{pd}(\nu) \, \frac{F_\nu}{h \, \nu} \, \mathrm{d}\nu \, ,
\end{equation}$

where

 *  $F_\nu$ is the energy flux of photons per unit area per unit frequency $[\mathrm{erg \, s^{-1} \, cm^{-2} \, Hz^{-1}}]$,
 *  $h \, \nu$ is the energy of a single photon, $F_\nu / (h \, \nu)$ is therefore the number flux of photons $[\mathrm{photons \, s^{-1} \, cm^{-2} \, Hz^{-1}}]$,
 *  $\sigma_{pd}(\nu)$ is the photodissociation cross-section $[\mathrm{cm^{-2}}]$.

Dividing both sides by the constant $\rho_\text{cell}$ we get the LWB photodissociation term for our equations

$\begin{equation}
	\frac{\mathrm{d}}{\mathrm{d}t} f_m(t) = -f_m(t) \int_{LW} \sigma_{pd}(\nu) \, \frac{F_\nu}{h \, \nu} \, \mathrm{d}\nu =: -f_m(t) \, \Gamma_\mathrm{LWB}(z) \, .
\end{equation}$

Assuming an homogeneous background field we can write the energy flux as

$\begin{equation}
	F_\nu = 4\pi \, J_\nu \, ,
\end{equation}$

where $J_\nu$ is the radiation intensity, wich is generally parametrized as

$\begin{equation}
	J_\nu(z) = J_{21}(z) \times 10^{-21} \, \mathrm{erg \, s^{-1} \, cm^{-2} \, Hz^{-1} \, sr^{-1}} \, .
\end{equation}$

Putting all together we can write the photodissociation rate as

$\begin{equation}
	\Gamma_\mathrm{LWB}(z) = \left(\frac{4\pi}{h} \times 10^{-21}\right) J_{21}(z) \int_{\mathrm{LW \, band}} \frac{\sigma_{pd}(\nu)}{\nu} \mathrm{d}\nu \, .
\end{equation}$

The integral part contains all the complex quantum mechanics of the $\mathrm{H}_2$ molecule's many absorption lines. This integral has been numerically calculated and is treated as a constant value.

From [Abel1997](https://doi.org/10.1016/S1384-1076(97)00010-9) we have (in the optically thin limit)

$\begin{equation}
	\Gamma_\mathrm{LWB}(z) = (1.38 \times 10^{-12} \, \mathrm{s^{-1}}) \, J_{21}(z) \, .
\end{equation}$

Now, If we want to consider the shielding effects of dust and the self-shielding of molecular hydrogen we simply add the terms $S_d \, S_{\mathrm{H}_2}$ previously computed

$\begin{equation}
	\frac{\mathrm{d}}{\mathrm{d}t} f_m(t) = -\frac{\mathrm{d}}{\mathrm{d}t} f_a(t) = -f_m(t) \, S_d \, S_{\mathrm{H}_2} \, \Gamma_\mathrm{LWB} \, .
\end{equation}$

For $J_{21}(z)$ we will use the fit in eq. 9 of [Incatasciato2023](https://doi.org/10.1093/mnras/stad1008)

$\begin{equation}
	\log_{10} J_{21}(z) = A + B \, (1 + z) + C \, (1 + z)^2 \, ,
\end{equation}$

where $A = 2.119$, $B = -1.117 \times 10^{-1}$, and $C = -2.782 \times 10^{-3}$.
"""
  ╠═╡ =#

# ╔═╡ 5ce798b5-4366-4cc1-8d2d-881a827d6dc9
begin
	const abel97 = ustrip(t_u^-1, 1.38e-12u"s^-1")
	const LW_SOURCE = "Incatasciato2023"

	if LW_SOURCE == "Incatasciato2023"
		#############################################################################
		# LW background radiation intensity from Incatasciato et al. (2023)
		# in [10^(-21) erg s^(-1) cm^(-2) Hz^(-1) sr^(-1)]
		#############################################################################
		const lwb_A = 2.119
		const lwb_B = -1.117e-1
		const lwb_C = -2.782e-3
	elseif LW_SOURCE == "Ahn2009"
		#############################################################################
		# LW background radiation intensity from Ahn et al. (2009)
		# in [10^(-21) erg s^(-1) cm^(-2) Hz^(-1) sr^(-1)]
		#############################################################################
		const lwb_A = 3.5368167124911394
		const lwb_B = -1.3476447221880308e-1
		const lwb_C = -9.560105322484195e-3
	elseif LW_SOURCE == "Visbal2014"
		#############################################################################
		# LW background radiation intensity from Visbal et al. (2014)
		# in [10^(-21) erg s^(-1) cm^(-2) Hz^(-1) sr^(-1)]
		#############################################################################
		const lwb_A = 0.7050655910103948
		const lwb_B = 5.670120801889934e-2
		const lwb_C = -3.366590865530786e-3
	else
		throw(ArgumentError("I don't recognize the LW background source $(LW_SOURCE)."))
	end
	
	J21(z) = exp10(lwb_A + lwb_B * (1 + z) + lwb_C * (1 + z)^2)

	#################################################################################
	# Compute ΓLWB
	#
	# Arguments
	#
	#     a: Scale factor [dimensionless]
	#
	# Returns
	#
	#     ΓLWB [Myr^(-1)]
	#################################################################################
	function compute_LWB(a::Float64)::Float64
	
		z = (1 / a) - 1
	
		return abel97 * J21(z)
	
	end
end;

# ╔═╡ 7df63fdf-ef41-44f2-8a55-e5c2c849029c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Mass recycling

We consider two key parameters for mass recycling, $R$, the fraction of stellar mass that is returned to the ISM, under the instantaneous recycling approximation, and $Z_\mathrm{SN}$, the metallicity of the recycled gas, i.e., the fraction of the returned mass composed of metals (the rest is assumed to be ionized hydrogen).

The instantaneous recycling approximation (IRA) assumes that stars above a certain mass die and return material to the ISM instantaneously, while lower-mass stars live indefinitely. This approach was first formalized by [Schmidt1963](https://doi.org/10.1086/147553).

Although more accurate treatments avoid the IRA by considering stellar lifetimes (e.g., using empirical relations as done in [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635), Sections 2.2.2-2.2.3), this introduces a time-dependence in $R$, requiring the evaluation of stellar evolution integrals at every timestep -- greatly increasing computational cost. For simplicity and efficiency, we assume $R$ is constant over the timescales relevant to our ODE integrations.

The recycled mass fraction $R$ can be computed using a stellar yield model, which specifies the mass of each element ejected by stars of different initial masses. Mathematically,

$\begin{equation}
	R = \dfrac{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} (m - m_\mathrm{rem}(m)) \, \phi(m) \, \mathrm{d}m}{\int_{m_\mathrm{low}}^{m_\mathrm{high}} m \, \phi(m) \, \mathrm{d}m} \, ,
\end{equation}$

where $\phi(m)$ is the initial mass function (IMF), $m_\mathrm{rem}(m)$ is the remnant mass for a star of mass $m$, $m_\mathrm{low}$ and $m_\mathrm{high}$ define the full mass range of the IMF, and $m_\mathrm{ir}$ is the mass threshold above which the IRA applies.

The denominator represents the total mass of the stellar population and acts as a normalization factor for the IMF.

Similarly, the metal yield of the recycled gas is given by

$\begin{equation}
	Z_\mathrm{SN} = \dfrac{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} m \, f_Z \, \phi(m) \, \mathrm{d}m}{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} (m - m_\mathrm{rem}(m)) \, \phi(m) \, \mathrm{d}m} \, ,
\end{equation}$

where $f_Z(m)$ is the fraction of the stellar mass returned as metals.

The two previous equations were taken from [_Nucleosynthesis and Chemical Evolution of Galaxies_](https://doi.org/10.1017/CBO9780511812170) by Bernard Pagel (eq. 7.24 and 7.26), and [_Chemical Evolution
of Galaxies_](https://doi.org/10.1007/978-94-010-0967-6) by Francesca Matteucci (eq. 2.74).

A typical choice for the stellar mass limits is:

*  $m_\mathrm{low} = 0.08 \, M_\odot$ (the hydrogen-burning limit),
*  $m_\mathrm{ir} = 8 \, M_\odot$ (approximate threshold for core-collapse supernovae),
*  $m_\mathrm{high} = 100 \, M_\odot$ (upper limit constrained by observations and yield model validity).

[Ascasibar2015](https://doi.org/10.1093/mnras/stv098) found $R \approx 0.18$ and $Z_\mathrm{sn} \approx 0.09$, using the yield model of [Woosley1995](https://doi.org/10.2172/115557) and the IMF of [Kroupa2001](https://doi.org/10.1046/j.1365-8711.2001.04022.x), although the exact mass integration limits used were not specified.

The stellar yields used to compute $R$ and $Z_\mathrm{SN}$ can be selected from the following models

*  [Woosley1995](https://doi.org/10.2172/115557)
*  [Portinari1998](https://doi.org/10.48550/arXiv.astro-ph/9711337)
*  [Chieff2004](https://doi.org/10.1086/392523)
*  [Kobayashi2006](https://doi.org/10.1086/508914)
*  [Heger2010](https://doi.org/10.1088/0004-637X/724/1/341)
*  [Limongi2012](https://doi.org/10.1088/0067-0049/199/2/38)

These datasets are conveniently compiled and discussed in [Mollá2015](https://doi.org/10.1093/mnras/stv1102).

In what follows we will use the stellar yield model of [Portinari1998](https://doi.org/10.48550/arXiv.astro-ph/9711337)
"""
  ╠═╡ =#

# ╔═╡ 64c2a3ee-e7a0-4e03-b9a8-86239e1ca81e
##################################################################################
# Compute the mass fraction of ejected metals for a given stellar yield model
#
# Arguments
#
#     path:     Path to the files with the SED data
#     sy_model: Target stellar yield model
#     header:   Column names of the files
#     metals:   Columns names of the metal yields
#     skipto:   First row with data
#
# Returns
#
#     A dataframe with columns:
# 
#       :Z      -> Metallicity of the stellar population modeled by the IMF
#       :ms     -> Stellar mass
#       :m_rem  -> Remnant mass, after stellar death
#       :fz_rem -> Fraction of the stellar mass ejected as metals to the ISM
##################################################################################
function compute_sy(
	path::String, 
	sy_model::String;
	header::Vector{Symbol}=[
		:Z,    # Metallicity
		:m,    # Stellar mass
		:H,    # Hydrogen
		:D,    # Deuterium
		:He3,  # Helium 3
		:He4,  # Helium 4
		:C12,  # Carbon 12
		:O16,  # Oxygen 16
		:Ne20, # Neon 20
		:Mg24, # Magnesium 24
		:Si28, # Silicon 28
		:S32,  # Sulfur 32
		:Ca40, # Calcium 40
		:Fe56, # Iron 56
		:Mrem, # Remnant mass
		:C13s, # Secondary contribution of carbon 13 
		:N14s, # Secondary contribution of Nitrogen 14
    ],
	metals::Vector{Symbol}=[
		:C12, 
		:O16, 
		:Ne20, 
		:Mg24, 
		:Si28, 
		:S32, 
		:Ca40, 
		:Fe56, 
		:C13s, 
		:N14s,
	],
	skipto::Int64=2,
)::DataFrame

	# List the files with the stellar yields
    files = readdir(path, join=true)

	# Parse the filenames
	filenames = @. uppercase(getindex(split(basename(files), "."), 1))

	# Read only the files for the target IMF
	target_file = files[findfirst(filenames .== sy_model)]

	# Load the data
    df = identity.(CSV.read(target_file, DataFrame; header, skipto))

	# Stellar mass
	ms = df[!, :m] .* u"Msun" 

	# Remnant mass 
    m_rem = df[!, :Mrem] .* u"Msun" 

	# Remnant mass fraction
	fm = m_rem ./ ms

	# Total stellar yield
	q = sum(eachcol(df[!, metals]))

	# Metallicity
	Z = df[!, :Z]

	# Mass fraction of metals ejected by the star 
	# throughout its evolution and death
	# See eqs. 1 and 2 from Mollá et al. (2015)
	# https://doi.org/10.1093/mnras/stv1102
	fz_ejec = @. q + (1 - fm) * Z

	sy_data = sort(identity.(DataFrame(; Z, ms, m_rem, fz_ejec)), [:Z, :ms])

	return sy_data
	
end;

# ╔═╡ 0895d464-a029-410d-9e7d-89cfac2d1615
begin
	const YIELD_MODEL_KEYS = Dict(
        "Woosley1995"   => "WOW",
        "Portinari1998" => "PCB",
        "Chieff2004"    => "CLI",
        "Kobayashi2006" => "KOB",
        "Heger2010"     => "HEG",
        "Limongi2012"   => "LIM",
    )

	# Chosen stellar yield model
    const YIELD_MODEL = "Portinari1998"
	
	# Mass limit for the IR hypothesis
	const M_IR = 8.0u"Msun"              
	
    # Raw stellar yields from Mollá et al. (2015)
	# https://doi.org/10.1093/mnras/stv1102
    const sy_path = "./data/stellar_yields"

	sy_df = compute_sy(sy_path, YIELD_MODEL_KEYS[YIELD_MODEL])

	# Metallicities
    sy_metals = sort!(unique(sy_df[!, :Z]))
	
    # Stellar masses
    sy_masses = sort!(ustrip.(u"Msun", unique(sy_df[!, :ms])))

	########################################################
	# Reshape the dataframe into matrices for interpolation
	# Rows    -> :ms
	# Columns -> :Z
	########################################################

	wide_mrem = unstack(sy_df, :ms, :Z, :m_rem)
	sort!(wide_mrem, :ms)

	wide_fz = unstack(sy_df, :ms, :Z, :fz_ejec)
	sort!(wide_fz, :ms)

    mrem_mat = ustrip.(u"Msun", Matrix(wide_mrem[:, Not(:ms)]))
	fz_mat   = Matrix(wide_fz[:, Not(:ms)])

	# Create the interpolation function for m_rem
	mrem_interp = linear_interpolation(
		(sy_masses, sy_metals),
		mrem_mat,
		extrapolation_bc=Flat(),
	)

	# Create the interpolation function for fz_ejec
	fz_interp = linear_interpolation(
		(sy_masses, sy_metals),
		fz_mat,
		extrapolation_bc=Flat(),
	)

	sy_df
end

# ╔═╡ 07cd9aad-029f-42c6-abe8-ab4a9a2a910c
#####################################################################################
# Compute R and Zsn.
#
# Arguments
#
#     Z: Metallicity [dimensionless]
#
# Returns
#
#     (R [dimensionless], Zsn [dimensionless])
##################################################################################
function compute_recycled_fractions(Z::Float64)
	
	# Remnant mass
	mrem(m) = mrem_interp(ustrip(u"Msun", m), Z) * u"Msun"

	# Fraction of the stellar mass ejected as metals to the ISM
	fz(m) = fz_interp(ustrip(u"Msun", m), Z)

	integrand_R(m)   = (m - mrem(m)) * ϕ(m)
	integrand_Zsn(m) = m * fz(m) * ϕ(m)

	R_int = quadgk(integrand_R, M_IR, M_HIGH)[1]
	Zsn_int = quadgk(integrand_Zsn, M_IR, M_HIGH)[1]

	R = uconvert(Unitful.NoUnits, R_int / IMF_NORM)
	Zsn = uconvert(Unitful.NoUnits, Zsn_int / R_int)

	return R, Zsn

end;

# ╔═╡ 0bdf9dbf-479c-46f6-bd86-50576095cba0
# ╠═╡ skip_as_script = true
#=╠═╡
md"# Implementation"
  ╠═╡ =#

# ╔═╡ 6fe43e3a-2e8f-4708-a3ec-6f5a8088060e
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Constants"
  ╠═╡ =#

# ╔═╡ f863d68f-590e-4b96-8433-dc6b5177539f
begin
	# Number of equations
	const N_EQU = 6  
	
	# Number of parameters
	const N_PAR = 8 

	# Index of each phase in the ODE solution matrix
	const phase_name_to_index = Dict(
		"ionized"   => 1,
		"atomic"    => 2,
		"molecular" => 3,
		"stellar"   => 4,
		"metals"    => 5,
		"dust"      => 6,
	)
end;

# ╔═╡ f5a983bf-ef3a-4d5d-928d-da1097b91ee8
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Equations"
  ╠═╡ =#

# ╔═╡ 9ab0a10b-8165-401a-a2b6-c9726526a906
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Jacobian"
  ╠═╡ =#

# ╔═╡ 1b8af600-56eb-4508-bc52-aa4e138b4c7e
#####################################################################################
# Construct the Jacobian matrix with each term as a Julia function
#
# Arguments
#
#     system: ODE system. This function should follows the API requirements of 
#             DifferentialEquations.jl
#
# Returns
#
#     A matrix of Julia functions
#####################################################################################
function construct_jacobian(system::Function)::Matrix{Function}

	@variables S_t S_ic[1:N_EQU] S_parameters[1:N_PAR]
    S_dydt = Vector{Num}(undef, N_EQU)
    system(S_dydt, S_ic, S_parameters, S_t)

	# Compute the Jacobian symbolically
    jac = Symbolics.jacobian(S_dydt, S_ic)

    # Transform the symbolic expresions into julia functions
    return [
        build_function(jac[i, j], S_ic, S_parameters, S_t, expression=Val{false}) for
        i in 1:N_EQU, j in 1:N_EQU
    ]

end;

# ╔═╡ 4fcf6e4d-b9e7-43d1-badc-d8d8afc004b3
#####################################################################################
# Personalized method of fma() for the computation of the jacobian
#####################################################################################
function Base.fma(x::Symbolics.Num, y::Symbolics.Num, z::Symbolics.Num)
	return x * y + z
end;

# ╔═╡ bd8743d6-8f21-413d-835a-e543926baa09
#####################################################################################
# System of ODEs, where each equation has units of Myr^(-1), and
#
# Ionized gas fraction:    fi(t) = Mi(t) / M_cell --> y[1]
# Atomic gas fraction:     fa(t) = Ma(t) / M_cell --> y[2]
# Molecular gas fraction:  fm(t) = Mm(t) / M_cell --> y[3]
# Stellar fraction:        fs(t) = Ms(t) / M_cell --> y[4]
# Metals fraction:         fZ(t) = MZ(t) / M_cell --> y[5]
# Dust fraction:           fd(t) = Md(t) / M_cell --> y[6]
# 
# This function follows the API requirements of DifferentialEquations.jl
#####################################################################################
function system!(dydt, ic, parameters, t)

	##################################################################
    # Initial conditions
	##################################################################

    fi, fa, fm, fs, fZ, fd = ic

	##################################################################
    # Parameters
	##################################################################
	#
	# ρ_cell: Total cell density           [mp * cm^(-3)]
	# ΓUVB:   UVB photoionization rate     [Myr^(-1)]
	# ΓLWB:   LWB photodissociation rate   [Myr^(-1)]
	# η_diss: Photodissociation efficiency [dimensionless]
	# η_ion:  Photoionization efficiency   [dimensionless]
	# R:      Mass recycling fraction      [dimensionless]
	# Zsn:    Metals recycling fraction    [dimensionless]
	# h:      Column height                [cm]
	##################################################################

    ρ_cell, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h = parameters

	##################################################################
    # Auxiliary equations
	##################################################################

	# Star formation rate [Myr^(-1)]
	ψs = ψ(fm, ρ_cell)

	###################
	# Partial fraction
    ###################

	# Neutral fraction
	fn = fa + fm

	# Gas fraction
	fg = fn + fi

	# Metal fraction
	Zt = fZ + fd

	############
	# Shielding
    ############

	# Dust optical depth 
	τ_d = c_sd * h * ρ_cell * fn * Zt
	
	# Dust shielding 
	sd = exp(-τ_d)

	# Molecular self-shielding optical depth
	xp1     = c_sh2 * h * fm * ρ_cell + 1.0
	xp1_2   = xp1 * xp1
	sq_xp1  = NaNMath.sqrt(xp1)
	τ_sh    = exp_fac * sq_xp1

	# Molecular self-shielding
	exp_xp1 = exp(-τ_sh)
	sh2     = ((1 - ωH2) / xp1_2) + (ωH2 * exp_xp1 / sq_xp1)

	# Combined shielding 
	e_diss = sd * sh2

	#################
	# Net ionization
    #################

	# Recombination [Myr^(-1)]
	recomb = c_rec * fi * fi * ρ_cell

	# Ionization optical depth
	τ_ion = c_τion * fa * ρ_cell * h

	# Stellar ionization [Myr^(-1)]
	ion_prob = -expm1(-τ_ion)
	s_ion    = sd * η_ion * ion_prob * ψs

	# UVB photoionization [Myr^(-1)]
	uvb = sd * ΓUVB * fa

	net_ionization = uvb + s_ion - recomb

	#########################
	# Stellar gas production
    #########################

	# Gas and metals produced at stellar death [Myr^(-1)]
	R_ψs = R * ψs

	s_gas_production = fma(R_ψs, -Zsn, R_ψs)

	###################
	# Net dissociation
    ###################

	# Condensation [Myr^(-1)]
	Zt_eff = Zt + Zeff
	cond   = c_cond * fa * ρ_cell * fg * Zt_eff

	# Dissociation optical depth
	τ_diss = c_τdiss * fm * ρ_cell * h

	# Stellar dissociation [Myr^(-1)]
	diss_prob = -expm1(-τ_diss)
	s_diss    = e_diss * η_diss * diss_prob * ψs

	# LW background dissociation [Myr^(-1)]
	lwb = e_diss * ΓLWB * fm

	net_dissociation = lwb + s_diss - cond

	##################
	# Net dust growth
	##################

	# Dust growth  [Myr^(-1)]
	dg = c_dg * fZ * fd * fn * fn * ρ_cell

	# Dust destruction [Myr^(-1)]
	dd = fd * inv_τ_dd

	net_dust_growth = dg - dd

	##################################################################
    # ODE system
	##################################################################

	dydt[1] = net_ionization + s_gas_production
    dydt[2] = net_dissociation - net_ionization
    dydt[3] = -net_dissociation - ψs
    dydt[4] = fma(R, -ψs, ψs)
	dydt[5] = fma(Zsn, R_ψs, -net_dust_growth)
	dydt[6] = net_dust_growth

end;

# ╔═╡ 80005099-7154-4306-9172-c9a168336e14
const JACOBIAN_FUNCTION = construct_jacobian(system!);

# ╔═╡ c291700e-3a84-49a7-85d6-592cfb3b1a11
#####################################################################################
# Evaluate the Jacobian
#
# Arguments
# 
#     J:          Matrix to save the results, it must have size N_EQU × N_EQU
#     ic:         Initial condition, [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
#     parameters: Parameters for the ODEs, [ρ, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h]
#     t:          Unused variable to comply with the DifferentialEquations.jl API
# 
# Returns
#
#     The Jacobian matrix
#####################################################################################
function jacobian!(
	J::Matrix{Float64},
	ic::Vector{Float64},
	parameters::Vector{Float64},
	t::Float64=0.0,
)

	for i in 1:N_EQU
	    for j in 1:N_EQU
	        J[i, j] = JACOBIAN_FUNCTION[i, j](ic, parameters, t)
	    end
	end

end;

# ╔═╡ 1ec99792-905d-4e1b-a413-ef58143d3c68
# ╠═╡ skip_as_script = true
#=╠═╡
md"## ODE function"
  ╠═╡ =#

# ╔═╡ 4e1140be-3e61-47bd-9f9f-6b5dfbff6ae2
#####################################################################################
# System to be solved by DifferentialEquations.jl
#####################################################################################
ode_function = ODEFunction{true}(
	system!;
	jac=jacobian!,
	tgrad=(dt, ic, p, t) -> nothing,
);

# ╔═╡ 2afd459c-90e9-4105-9121-27e21bb89eeb
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Integration"
  ╠═╡ =#

# ╔═╡ 4173e91c-c01a-4fa2-ae77-0408bf7d9a1b
#####################################################################################
# Solve the system of ODEs
#
# Arguments
# 
#     ic:          Initial conditions, [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
#     base_params: Parameters for the ODEs, [ρ_cell, Z, a, h]
#     tspan:       Integration span, (ti, tf) [Myr]
#     times:       Times at which the solution will be returned [Myr]
#     args:        Positional arguments for the solver of DifferentialEquations.jl
#     kwargs:      Keyword arguments for the solver of DifferentialEquations.jl
# 
# Returns
#
#     A vector of vectors [fi(t), fa(t), fm(t), fs(t), fZ(t), fd(t)] 
#     for each t in `times`
#####################################################################################
function integrate_model(
    ic::Vector{Float64},
	base_params::Vector{Float64},
	tspan::NTuple{2,Float64};
    times::Vector{Float64}=[tspan[2],],
	args::Tuple=(),
    kwargs::NamedTuple=(
		dense=false,
		reltol=1.0e-10,
		verbose=false,
	),
)::Vector{Vector{Float64}}

	ρ_cell, Z, a, h = base_params

	ΓUVB          = compute_UVB(a)
	ΓLWB          = compute_LWB(a)
	η_diss, η_ion = compute_η(tspan[2], Z)
	R, Zsn        = compute_recycled_fractions(Z)

	parameters = [ρ_cell, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h]

    sol = solve(ODEProblem(
		ode_function,
		ic,
		tspan,
		parameters,
	), args...; kwargs...)

    return sol(times).u

end;

# ╔═╡ 24920a1b-355c-43e1-aff2-028e6f25584c
md"# Extras"

# ╔═╡ 8fd489b2-620a-40cd-8283-dc88b7f3584f
# ╠═╡ skip_as_script = true
#=╠═╡
md"## ODE coefficients"
  ╠═╡ =#

# ╔═╡ 59a342e1-930d-40c4-86fc-11ba334d7103
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Stiffness coefficient"
  ╠═╡ =#

# ╔═╡ a2a12511-97d6-4e83-a21d-1c7528f486bc
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Compute the stiffness ratio
#
# Arguments
# 
#     ic:         Initial conditions, [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
#     base_parms: Parameters for the ODEs, [ρ_cell, Z, a, h, it]
#
# Returns
#
#     The stiffness ratio
#####################################################################################
function stiffness_ratio(
	ic::Vector{Float64},
	base_params::Vector{Float64},
)::Float64

	# Parameters for the ODEs
	# ρ_cell: Density [cm^(-3)]
	# Z:      Metallicity [dimensionless]
	# a:      Scale factor [dimensionless]
	# h:      Column height [cm]
	# it:     Integration time [Myr]
	ρ_cell, Z, a, h, it = base_params

	ΓUVB          = compute_UVB(a)
	ΓLWB          = compute_LWB(a)
	η_diss, η_ion = compute_η(it, Z)
	R, Zsn        = compute_recycled_fractions(Z)

	parameters = [ρ_cell, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h]

	# Allocate memory for the Jacobian
	J = Matrix{Float64}(undef, N_EQU, N_EQU)

	# Compute the Jacobian and store it in J
	jacobian!(J, ic, parameters)

	# Get the norm of the real part of the non-zero eigenvalues
	eigen_values = filter(x -> x > eps(typeof(x)), eigvals(J) .|> real .|> abs)

	# Compute the stiffness ratio
	return maximum(eigen_values) / minimum(eigen_values)

end;
  ╠═╡ =#

# ╔═╡ 8ade720c-37b7-4436-8bcd-8e902f52057b
# ╠═╡ skip_as_script = true
#=╠═╡
stiffness_ratio(
	# [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
	[0.3, 0.57, 0.0, 0.0, 0.02, 0.01],
	# [ρ_cell [cm^(-3)], Z [dimensionless], a [dimensionless], h [cm], it [Myr]]
	[1.0, Zsun, 1.0, 3.0e20, 5.0],
)
  ╠═╡ =#

# ╔═╡ b39e416e-a204-4f6b-a6e4-00f6a71d3b27
# ╠═╡ skip_as_script = true
#=╠═╡
md"The model, for a typical set of initial conditions, is very stiff (stiffness ratio >> 1)."
  ╠═╡ =#

# ╔═╡ c414f207-0d4f-4cc4-965c-9b68b5a85400
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Lyapunov spectrum"
  ╠═╡ =#

# ╔═╡ ac57dfbe-f19b-4c5f-a753-5f4162e245d4
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Compute the Lyapunov spectrum
#
# Arguments
# 
#     ic:         Initial conditions, [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
#     base_parms: Parameters for the ODEs, [ρ_cell, Z, a, h, it]
#
# Returns
#
#     The Lyapunov spectrum
#####################################################################################
function lyapunov_spectrum(
	ic::Vector{Float64},
	base_params::Vector{Float64},
)::Vector{Float64}

	# Parameters for the ODEs
	# ρ_cell: Density [cm^(-3)]
	# Z:      Metallicity [dimensionless]
	# a:      Scale factor [dimensionless]
	# h:      Column height [cm]
	# it:     Integration time [Myr]
	ρ_cell, Z, a, h, it = base_params

	ΓUVB          = compute_UVB(a)
	ΓLWB          = compute_LWB(a)
	η_diss, η_ion = compute_η(it, Z)
	R, Zsn        = compute_recycled_fractions(Z)

	parameters = [ρ_cell, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h]
	
	ds = CoupledODEs(system!, ic, parameters)

	# Compute the Lyapunov spectrum
	return lyapunovspectrum(ds, 100; Δt=1e-6)

end;
  ╠═╡ =#

# ╔═╡ b1e67fcd-f70b-4dd2-87e5-74bef1fa201e
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
lyapunov_spectrum(
	# [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
	[0.3, 0.57, 0.0, 0.0, 0.02, 0.01],
	# [ρ_cell [cm^(-3)], Z [dimensionless], a [dimensionless], h [cm], it [Myr]]
	[10.0, Zsun, 1.0, 3.0e20, 10.0],
)
  ╠═╡ =#

# ╔═╡ d1f3d076-b9d6-4b92-b25c-4155b5fc464b
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
lyapunov_spectrum(
	# [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
	[0.3, 0.57, 0.0, 0.0, 0.02, 0.01],
	# [ρ_cell [cm^(-3)], Z [dimensionless], a [dimensionless], h [cm], it [Myr]]
	[1000.0, Zsun, 1.0, 3.0e20, 10.0],
)
  ╠═╡ =#

# ╔═╡ f336b083-99c7-4345-8ada-78b166a98abe
# ╠═╡ skip_as_script = true
#=╠═╡
md"The model, for a typical set of initial conditions, is not chaotic (maximum lyapunov exponent < 0)."
  ╠═╡ =#

# ╔═╡ 74c82c98-f233-4e72-8f74-17bdfeddf884
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Dust initial condition"
  ╠═╡ =#

# ╔═╡ ab1b7695-3cd6-4e67-9cfd-fbaf2f4d1d15
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Schematic diagram for the contribution of each phase to the initial condition.
"""
  ╠═╡ =#

# ╔═╡ 6475b25b-8711-44cd-bc3b-e3d42681bb93
# ╠═╡ skip_as_script = true
#=╠═╡
TikzPictures.TikzPicture(
	L"""
	    % Base line
	    \draw[thick] (0,0) -- (10,0);

	    % Vertical dividers
	    \foreach \x in {0, 5, 10} {
	        \draw[thick] (\x, -0.2) -- (\x, 0.2);
	    }

	    % Top braces
	    \draw[thick] (0,0.5) .. controls (2.5,0.8) and (2.5,0.8) .. (5,0.5);
	    \draw[thick] (5,0.5) .. controls (7.5,0.8) and (7.5,0.8) .. (10,0.5);

	    \node at (2.5,1.0) {metals};
	    \node at (7.5,1.0) {H + He};

	    % Bottom braces
	    \draw[thick] (0,-0.5) .. controls (0.5,-0.8) and (1.5,-0.8) .. (2,-0.5);
	    \node at (1,-1.0) {dust};

	    \draw[thick] (2,-0.5) .. controls (2.5,-0.8) and (4.5,-0.8) .. (5,-0.5);
	    \node at (3.5,-1.0) {$Z$};

	    \draw[thick] (5,-0.5) .. controls (5.5,-0.8) and (7,-0.8) .. (7.5,-0.5);
	    \node at (6.25,-1.0) {$a$};

	    \draw[thick] (7.5,-0.5) .. controls (8,-0.8) and (9.5,-0.8) .. (10,-0.5);
	    \node at (8.75,-1.0) {$i$};

	    % Left label
	    \node[left] at (0,0) {cell};
	""",
	width="75em",
	preamble = """
		\\usepackage{xcolor}
	    \\color{white}
	""",
)
  ╠═╡ =#

# ╔═╡ ea4e58e9-d041-4a6e-b0d8-83e3aef7648b
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Following eq. 18 of [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x), the initial dust mass density is given by

$\begin{equation}
	\rho_d(0) = f_d \, \rho_\mathrm{cell} = \frac{m_\mathrm{X_d}}{f_\mathrm{X_d}}  \, (1 - \xi_0) \, \frac{Z}{Z_\odot} \, \left( \mathrm{\frac{X_d}{H}} \right)_{\!\odot} \, n_H \, ,
\end{equation}$

where $f_\mathrm{X_d}$ is the mass fraction in the dust of element $\mathrm{X_d}$, $m_\mathrm{X_d}$ is the atomic mass of $\mathrm{X_d}$, $\left( \mathrm{X_d / H} \right)_\odot$ is the solar abundance of $\mathrm{X_d}$, $n_H$ is the number density of hydrogen, and $\xi_0$ is the initial gas-phase fraction of $\mathrm{X_d}$, defined by

$\begin{equation}
	\xi(t) = \frac{n_\mathrm{X_d}^\mathrm{gas}}{n_\mathrm{X_d}^\mathrm{tot}} \, ,
\end{equation}$

with $n_\mathrm{X_d}^\mathrm{gas}$ the number density of $\mathrm{X_d}$ in the gas phase, and $n_\mathrm{X_d}^\mathrm{tot}$ the total number density (in gas and dust).

Using the following values

|  Species | $\mathrm{X_d}$ | $f_\mathrm{X_d}$ | $m_\mathrm{X_d}$ | $(\mathrm{X_d / H})_\odot$ |
|:--------:|:--------------:|:----------------:|:----------------:|:--------------------------:|
| Silicate |       Si       |      $0.166$     |      $28.1$      |    $3.55 \times 10^{-5}$   |
| Graphite |        C       |        $1$       |       $12$       |    $3.63 \times 10^{-4}$   |

and adopting the fiducial value $\xi_0 = 0.3$ from [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x), we can express the dust-to-gas mass ratio as

$\begin{equation}
    f_d = C_\mathrm{X_d} \, Z \, f_n \, ,
\end{equation}$

where

$\begin{align}
    n_H &= n_a + 2 \, n_m = f_n \, \frac{\rho_\mathrm{cell}}{m_p} \, , \\
	Z_\odot &= 0.0127 \, .
\end{align}$

The coefficient $C_\mathrm{X_d}$ is defined as

$\begin{equation}
    C_\mathrm{X_d} = \frac{m_\mathrm{X_d}}{f_\mathrm{X_d}}  \, (1 - \xi_0) \, \frac{1}{Z_\odot \, m_p} \, \left( \mathrm{\frac{X_d}{H}} \right)_{\!\odot} \, ,
\end{equation}$

Using the tabulated values, we find

|  Species | $\mathrm{X_d}$ | $C_\mathrm{X_d}$ |
|:--------:|:--------------:|:----------------:|
| Silicate |       Si       |      $0.329$     |
| Graphite |        C       |      $0.238$     |

Thus, our fiducial value for $C_\mathrm{X_d}$ is taken as the average 

$\begin{equation}
	C_\mathrm{X_d} = 0.283 \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ d7cc8a66-220a-4e9b-b42e-a0e085ed3a0f
begin
	# Initial fraction of Si and C in the gas pahse
	const ξ0 = 0.3

	# Initial dust fraction constant
	C_xd(mx, fx, XH) = (mx / fx) * (1.0 - ξ0) * XH / (Zsun * 1.0u"mp")

	# C_xd for silicate
	const silicate = C_xd(28.1u"u", 0.166, 3.55e-5)

	# C_xd for graphite
	const graphite = C_xd(12.0u"u", 1.0, 3.63e-4)

	# Nominal value of C_xd
	const c_xd = ustrip(Unitful.NoUnits, silicate + graphite) / 2.0
end;

# ╔═╡ 477c8f59-97d1-405d-a1c8-64b7e0b9119f
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Efficiency per free-fall time"
  ╠═╡ =#

# ╔═╡ 08e05c65-06a2-4560-92eb-b014dc7c3d70
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
As defined in [Krumholz2011](https://iopscience.iop.org/article/10.1088/0004-637X/745/1/69), the star formation efficiency per free-fall time is given by

$\begin{equation}
	\epsilon_\text{ff} = \text{SFR} \, \frac{t_\text{ff}}{M_g} \, ,
\end{equation}$

where $M_g$ is the gas mass under consideration(in our case, this is the cell mass: $M_\text{cell}$), $\text{SFR}$ the star formation rate (in our case, this refers to the probabilistic star formation rate: $\text{SFR}_p$), and $t_\text{ff}$ is the free-fall time, defined as

$\begin{equation}
	t_\text{ff} = \sqrt{\frac{3 \, \pi}{32 \, G \, \rho}} \, .
\end{equation}$

Substituting into the definition, we can express the efficiency as

$\begin{equation}
	\epsilon_\text{ff} = c \, \frac{\text{SFR}_p}{M_\text{cell} \, \sqrt{\rho_\text{cell}}} \, .
\end{equation}$

where $c = \sqrt{3 \, \pi / 32 \, G}$.

For our model we have

$\begin{equation}
	\text{SFR}_p = \frac{f_s \, M_\text{cell}}{t_i} \, ,
\end{equation}$

where $f_s$ is the stellar fraction from the ODEs and $t_i$ is the integration time. 

So, in our model, we end up with

$\begin{equation}
	\epsilon_\text{ff}^{\,\text{SFM}} = c \, \frac{f_s}{t_i \, \sqrt{\rho_\text{cell}}} \, .
\end{equation}$

In the case of the [Blitz2006](https://doi.org/10.1086/505417) model the $\text{SFR}_p$ is

$\begin{equation}
	\text{SFR}_p = 0.02 \, \frac{P}{P + P_0} \, \frac{M_\text{cell}}{t_\text{ff}} \, ,
\end{equation}$

where $P$ is the gas pressure and $P_0 / k_B = 2.0 \times 10^4 \, \mathrm{K \, cm^{-3}}$. 

So, the star formation efficiency is 

$\begin{equation}
	\epsilon_\text{ff}^{\,\text{BLT}} = 0.02 \, \frac{P}{P + P_0} \, .
\end{equation}$

Finally, in the fiducial model within $\texttt{AREPO}$ the $\text{SFR}_p$ is

$\begin{equation}
	\text{SFR}_p = \frac{x \, M_\text{cell} \, \sqrt{\rho}}{\sqrt{\rho_\text{th}} \, t_0^*} \, ,
\end{equation}$

where $x$ is the cold gas mass fraction, $\rho_\text{th}$ is the threshold density (a constant with a value of $1.03378 \times 10^6$ in internal units of density), and $t_0^*$ is the maximum SFR timescale (a constant with a value of $0.00227$ in internal units of time).

So, the star formation efficiency is 

$\begin{equation}
	\epsilon_\text{ff}^{\,\text{STD}} = c \, x \, ,
\end{equation}$

where $c$ is the constant

$\begin{equation}
	c = \sqrt{\frac{3 \, \pi}{32 \, G \, \rho_\text{th}}} \, \frac{1}{t_0^*} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 2235689d-9c83-4907-aa17-c2624fbeb68d
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Equilibrium fractions"
  ╠═╡ =#

# ╔═╡ df8a9449-851c-4546-97a7-7fd4a270a867
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
To find the per-equation and global equilibrium we take que ODEs

$\begin{align}
	\frac{\mathrm{d}}{\mathrm{d}t} f_i &= e_\mathrm{ion} \, (\eta_\mathrm{ion}^* \, \psi + \Gamma_\mathrm{UVB} \, f_a) + R \, (1 - Z_\mathrm{SN}) \, \psi - \frac{f_i}{\tau_\mathrm{rec}} \, , \\
	\frac{\mathrm{d}}{\mathrm{d}t} f_a &= \, \frac{f_i}{\tau_\mathrm{rec}} + e_\mathrm{diss} \, (\eta_\mathrm{diss}^* \, \psi + \Gamma_\mathrm{LWB} \, f_m) \, - e_\mathrm{ion} \, (\eta_\mathrm{ion}^* \, \psi + \Gamma_\mathrm{UVB} \, f_a) - \frac{f_a}{\tau_\mathrm{cond}} \, , \\
	\frac{\mathrm{d}}{\mathrm{d}t} f_m &= \frac{f_a}{\tau_\mathrm{cond}} - e_\mathrm{diss}  \, \eta_\mathrm{diss}^* \, \psi - e_\mathrm{diss} \, \Gamma_\mathrm{LWB} \, f_m - \psi \, , \\
	\frac{\mathrm{d}}{\mathrm{d}t} f_s &= \psi - R \, \psi \, , \\
	\frac{\mathrm{d}}{\mathrm{d}t} f_Z &= R \, Z_\mathrm{SN} \, \psi + \frac{f_d}{\tau_\mathrm{dd}} - \frac{f_Z}{f_Z + f_d} \, \frac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m) \, , \\
	\frac{\mathrm{d}}{\mathrm{d}t} f_d &= \frac{f_Z}{f_Z + f_d} \, \frac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m) - \frac{f_d}{\tau_\mathrm{dd}} \, ,
\end{align}$

and set the derivatives to $0$.
"""
  ╠═╡ =#

# ╔═╡ e8c51ae9-ded4-450a-80f9-c484c2b991d1
#####################################################################################
# Evaluate the ODEs
# 
# Arguments
# 
#     fractions:  Mass fractions, [fi, fa, fm, fs, fZ, fd]
#     base_parms: Parameters for the ODEs, (ρ_cell, Z, a, h, it)
# 
# Returns
#
#     The ODE evaluated at `fractions` and `base_params`
#####################################################################################
function equilibrium(
    fractions::Vector{Float64},
	base_params::Vector{Float64};
)::Vector{Float64}

	ρ_cell, Z, a, h, it = base_params

	ΓUVB          = compute_UVB(a)
	ΓLWB          = compute_LWB(a)
	η_diss, η_ion = compute_η(it, Z)
	R, Zsn        = compute_recycled_fractions(Z)

	parameters = [ρ_cell, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h]

	dydt = Vector{Float64}(undef, N_EQU)

    system!(dydt, fractions, parameters, 1.0)

    return dydt

end;

# ╔═╡ c96ea0ad-47c3-469e-a2f7-0743e681a6d1
#####################################################################################
# Test how close to 0 are the ODEs evaluated at `fractions` and `base_params`
# 
# Arguments
# 
#     fractions:  Mass fractions, [fi, fa, fm, fs, fZ, fd]
#     base_parms: Parameters for the ODEs, (ρ_cell, Z, a, h, it)
#     limit:      Equilibrium criteria.
# 
# Returns
# 
#     If the system pass the criteria (the derivatives are close enough to 0)
#####################################################################################
function test_equilibrium(
    fractions::Vector{Float64},
	base_params::Vector{Float64};
	limit::Float64=1e-10,
)::Vector{Bool}

	time_scale = base_params[5]

    return (equilibrium(fractions, base_params) ./ (1 / time_scale)) .< limit

end;

# ╔═╡ 0813f6a3-aadc-491c-ad1e-ec54dcbd0d56
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Density PDF"
  ╠═╡ =#

# ╔═╡ de4ba100-5f58-4657-aea0-a3f31eacda65
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Parameters"
  ╠═╡ =#

# ╔═╡ 77ca6eab-c5bf-479c-8982-1483933bbb6e
Base.@kwdef struct PDF_params
    # Density PDF function
    func::Function
    # Power law slope
    α::Float64
    # Dimensionless turbulent forcing parameter
    b::Float64
    # Mach number
    Ms::Float64
    # (min, max) values for s = ln(ρ/ρ₀)
    deviation::NTuple{2,Float64}
    # Number of divisions for the discretization of the density PDF
    divisions::Int64
end;

# ╔═╡ d4a1ec85-e6f5-48ed-9724-202c4dae9ae9
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Mass fractions"
  ╠═╡ =#

# ╔═╡ 7f8ebe64-ea86-4eb8-9939-6682180b9dd2
#####################################################################################
# Compute the density PDF mass fractions
#
# Arguments
# 
#     params:  Parameters for the density PDF
#     log_var: Selects which variable will be used
#              log_var == true:  s = ln(ρ/ρ₀) and logarithmic divisions
#              log_var == false: f = ρ/ρ₀ and linear divisions
# 
# Returns
# 
#     (mass fractions, list of corresponding densities)
#####################################################################################
function mass_fraction(params::PDF_params, log_var::Bool)::NTuple{2,Vector{Float64}}

	if params.divisions == 1
        return [1], [log_var ? 0 : 1.0]
    end

    # Select which variable will be used and how the function will be divided
    # log_var == true: s = ln(ρ/ρ₀) and logarithmic divisions
    # log_var == false: f = ρ/ρ₀ and linear divisions
    dev = log_var ? params.deviation : exp.(params.deviation)

    # Compute the step in the range of the variable s = ln(ρ/ρ₀) or f = ρ/ρ₀
    step = (dev[2] - dev[1]) / params.divisions

    # Compute the range of values of s = ln(ρ/ρ₀) or f = ρ/ρ₀
    points = [dev[1] + step * (i - 0.5) for i in 1:(params.divisions)]

    # Compute the fractions of mass within each division
    mass_f = [
        quadgk(
            x -> params.func(x, params),
            log_var ? point - (step / 2) : log(point - (step / 2)),
            log_var ? point + (step / 2) : log(point + (step / 2)),
            order=10,
            atol=10e-10,
        )[1] for point in points
    ]

    return (mass_f ./ sum(mass_f)), points

end;

# ╔═╡ aec6e4fc-e496-4add-b982-ab60f9f900a0
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Density PDF by Burkhart (2018)"
  ╠═╡ =#

# ╔═╡ 8940cb2b-2c8c-407e-bc6b-426843cf6125
#####################################################################################
# Compute the mass density PDF using the model in Burkhart (2018)
# https://doi.org/10.3847/1538-4357/aad002
#
# Arguments
# 
#     s:      Mass density, as ln(ρ/ρ₀)
#     params: Parameters, as the struct `PDF_params`
# 
# Returns
# 
#     The mass density PDF
#####################################################################################
function pBurkhart2018(s::Float64, params::PDF_params)::Float64

    b = params.b
    Ms = params.Ms
    α = params.α

    σs2 = log(1 + b^2 * Ms^2)
    s0 = -0.5 * σs2
    st = (α - 0.5) * σs2
    C = exp((α - 1) * 0.5 * α * σs2) / sqrt(2π * σs2)

	A_B18 = C * exp(-α * st) / α
	B_B18 = 0.5
	C_B18 = 0.5 * erf((2 * st + σs2) / sqrt(8 * σs2))
    N = 1 / (A_B18 + B_B18 + C_B18)

    if s < st
        return (N / sqrt(2π * σs2)) * exp(-((s - s0)^2) / (2 * σs2))
    else
        return N * C * exp(-α * s)
    end

end;

# ╔═╡ 28955b95-df19-403a-bf79-b68e9be8e1dd
begin
    const PDF_PARAMS = PDF_params(
		func = pBurkhart2018,      # Density PDF function
		α = 2.0,                   # Power law slope
		b = 0.5,                   # Dimensionless turbulent forcing parameter
		Ms = 10.0,                 # Mach number
		deviation = (-6, 6),       # (min, max) values for s = log(ρ/ρ₀)
		divisions = 20,            # Number of divisions for the discretization
								   # of the density PDF
	)

    # Pre computation of the default mass fractions
    # for each division of the density PDF
    const (MASS_FRAC, S_POINTS) = mass_fraction(PDF_PARAMS, true)
    const F_POINTS = exp.(S_POINTS)
end;

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"

[compat]
CSV = "~0.10.15"
DataFrames = "~1.8.1"
DifferentialEquations = "~7.17.0"
Interpolations = "~0.16.2"
Measurements = "~2.14.1"
NaNMath = "~1.1.3"
PlutoUI = "~0.7.79"
QuadGK = "~2.11.2"
SpecialFunctions = "~2.6.1"
Symbolics = "~7.9.1"
TikzPictures = "~3.5.1"
Unitful = "~1.27.0"
UnitfulAstro = "~1.2.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "68de101a0aa2cae5b54affa6ed5bf6c67c48d465"

[[deps.ADTypes]]
git-tree-sha1 = "f7304359109c768cf32dc5fa2d371565bb63b68a"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.21.0"
weakdeps = ["ChainRulesCore", "ConstructionBase", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "856ecd7cebb68e5fc87abecd2326ad59f0f911f3"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.43"

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
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.AlmostBlockDiagonals]]
deps = ["ConcreteStructs"]
git-tree-sha1 = "743abe5e5fe8cff96dad4123f263c0d8eee281c0"
uuid = "a95523ee-d6da-40b5-98cc-27bc505739d5"
version = "0.1.10"

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
git-tree-sha1 = "d81ae5489e13bc03567d4fbbb06c546a5e53c857"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.22.0"

    [deps.ArrayInterface.extensions]
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

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "e0b47732a192dd59b9d079a06d04235e2f833963"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.12.2"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "02fa77c70ba84361b9bc9ff28523bd9d78519265"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.11.0"

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"
    CliqueTreesExt = "CliqueTrees"

    [deps.BandedMatrices.weakdeps]
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.BoundaryValueDiffEq]]
deps = ["ADTypes", "BoundaryValueDiffEqAscher", "BoundaryValueDiffEqCore", "BoundaryValueDiffEqFIRK", "BoundaryValueDiffEqMIRK", "BoundaryValueDiffEqMIRKN", "BoundaryValueDiffEqShooting", "DiffEqBase", "FastClosures", "ForwardDiff", "LinearAlgebra", "Reexport", "SciMLBase"]
git-tree-sha1 = "d6ec33e4516b2e790a64128afdb54f3b536667a7"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "5.18.0"

    [deps.BoundaryValueDiffEq.extensions]
    BoundaryValueDiffEqODEInterfaceExt = "ODEInterface"

    [deps.BoundaryValueDiffEq.weakdeps]
    ODEInterface = "54ca160b-1b9f-5127-a996-1867f4bc2a2c"

[[deps.BoundaryValueDiffEqAscher]]
deps = ["ADTypes", "AlmostBlockDiagonals", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield"]
git-tree-sha1 = "47c833c459738a3f27c5b458ecf7832a4731ef4d"
uuid = "7227322d-7511-4e07-9247-ad6ff830280e"
version = "1.8.0"

[[deps.BoundaryValueDiffEqCore]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "ForwardDiff", "LineSearch", "LinearAlgebra", "Logging", "NonlinearSolveFirstOrder", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseConnectivityTracer", "SparseMatrixColorings"]
git-tree-sha1 = "b7b4d8cc80f116eab2eb6124dba58ea7aef31b85"
uuid = "56b672f2-a5fe-4263-ab2d-da677488eb3a"
version = "1.11.1"

[[deps.BoundaryValueDiffEqFIRK]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "325e6981a414cfa5181218936c23f0e16dee8f08"
uuid = "85d9eb09-370e-4000-bb32-543851f73618"
version = "1.9.0"

[[deps.BoundaryValueDiffEqMIRK]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "da6ae5e564ad06ced4d7504929c58130558007dd"
uuid = "1a22d4ce-7765-49ea-b6f2-13c8438986a6"
version = "1.9.0"

[[deps.BoundaryValueDiffEqMIRKN]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "609c2d03ea024df0d475fee483b93cf0e87c29d6"
uuid = "9255f1d6-53bf-473e-b6bd-23f1ff009da4"
version = "1.8.0"

[[deps.BoundaryValueDiffEqShooting]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "ba9bd1f31b58bfd5e48a56da0a426bcbd3462546"
uuid = "ed55bfe0-3725-4db6-871e-a1dc9f42a757"
version = "1.9.0"

[[deps.BracketingNonlinearSolve]]
deps = ["CommonSolve", "ConcreteStructs", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "750f782fcc7e09283be7d8a7aa687a95e4911b60"
uuid = "70df07ce-3d50-431d-a3e7-ca6ddb60ac1e"
version = "1.6.2"
weakdeps = ["ChainRulesCore", "ForwardDiff"]

    [deps.BracketingNonlinearSolve.extensions]
    BracketingNonlinearSolveChainRulesCoreExt = ["ChainRulesCore", "ForwardDiff"]
    BracketingNonlinearSolveForwardDiffExt = "ForwardDiff"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Preferences", "Static"]
git-tree-sha1 = "f3a21d7fc84ba618a779d1ed2fcca2e682865bab"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.7"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

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
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

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
git-tree-sha1 = "78ea4ddbcf9c241827e7035c3a03e2e456711470"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.6"

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
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d8928e9169ff76c6281f39a659f9bca3a573f24c"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.1"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DiffEqBase", "FastBroadcast", "ForwardDiff", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqRosenbrock", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLLogging", "SimpleNonlinearSolve", "SymbolicIndexingInterface"]
git-tree-sha1 = "c4b5c6d91acd09165f8cbabc4eefdc10350bce68"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.67.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "BracketingNonlinearSolve", "ConcreteStructs", "DocStringExtensions", "EnzymeCore", "FastBroadcast", "FastClosures", "FastPower", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "0514bf55835444420ce81f8e32c5e2369d4af456"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.199.0"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
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
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
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
git-tree-sha1 = "f17b863c2d5d496363fe36c8d8535cc6a33c9952"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "4.12.0"

    [deps.DiffEqCallbacks.extensions]
    DiffEqCallbacksFunctorsExt = "Functors"

    [deps.DiffEqCallbacks.weakdeps]
    Functors = "d9f16b24-f501-4c13-a1f2-28368ffc5196"

[[deps.DiffEqNoiseProcess]]
deps = ["CommonSolve", "DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "PoissonRandom", "QuadGK", "Random", "Random123", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "76fbb6985dda2aba32c97148540ad2c043e741e3"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.26.0"

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
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "JumpProcesses", "LinearAlgebra", "LinearSolve", "NonlinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "1df783c534cd0c4a865a397b1c4801771b5cbb07"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.17.0"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "5e6897d988addbfe7d9ad2ee467cc0c91001aae4"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.15"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = "GPUArraysCore"
    DifferentiationInterfaceGTPSAExt = "GTPSA"
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

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "fbcc7610f6d8348428f722ecbe0e6cfe22e672c6"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.123"

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
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "c249d86e97a7e8398ce2068dce4c078a1c3464de"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.16"

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
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "3f50fa86c968fc1a9e006c07b6bc40ccbb1b704d"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.4"

[[deps.EnumX]]
git-tree-sha1 = "7bebc8aad6ee6217c78c5ddcf7ed289d65d0263e"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.6"

[[deps.EnzymeCore]]
git-tree-sha1 = "990991b8aa76d17693a98e3a915ac7aa49f08d1a"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.18"
weakdeps = ["Adapt", "ChainRulesCore"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"
    EnzymeCoreChainRulesCoreExt = "ChainRulesCore"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "cc294ead6a85e975a8519dd4a0a6cb294eeb18d1"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.30.0"
weakdeps = ["StaticArrays"]

    [deps.ExponentialUtilities.extensions]
    ExponentialUtilitiesStaticArraysExt = "StaticArrays"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.FastAlmostBandedMatrices]]
deps = ["ArrayInterface", "ArrayLayouts", "BandedMatrices", "ConcreteStructs", "LazyArrays", "LinearAlgebra", "MatrixFactorizations", "PrecompileTools", "Reexport"]
git-tree-sha1 = "3733dfb413e2a87c790cdf34f32f2c6a6f7fb95e"
uuid = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
version = "0.1.6"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "ab1b34570bcdf272899062e1a56285a53ecaae08"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.3.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "0044e9f5e49a57e88205e8f30ab73928b05fe5b6"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "1.1.0"

[[deps.FastPower]]
git-tree-sha1 = "1dd291358e0e3a31f77dbbb76ce7abbcca38d410"
uuid = "a4df4552-cc26-4903-aec0-212e50a0e84b"
version = "1.3.0"

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
git-tree-sha1 = "9340ca07ca27093ff68418b7558ca37b05f8aeb1"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.29.0"

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
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "b2977f86ed76484de6f29d5b36f2fa686f085487"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.3.1"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "a694e2a57394e409f7a11ee0977362a9fafcb8c7"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.6"

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
git-tree-sha1 = "031d63d09bd3e6e319df66bb466f5c3e8d147bee"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.4"
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

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

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

[[deps.IntervalSets]]
git-tree-sha1 = "d966f85b3b7a8e49d034d27a189e9a4874b4391a"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.13"
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

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "PoissonRandom", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "SymbolicIndexingInterface"]
git-tree-sha1 = "0593944281222db3b3f3bc1cfe2c464a17f9eea1"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.21.0"

    [deps.JumpProcesses.extensions]
    JumpProcessesKernelAbstractionsExt = ["Adapt", "KernelAbstractions"]

    [deps.JumpProcesses.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "125d65fe5042faf078383312dd060adf11d90802"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.10.5"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

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

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "41d433e5854d7a67e8ab2b04962a713fcbcffcf1"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.9.5"

    [deps.LazyArrays.extensions]
    LazyArraysBandedMatricesExt = "BandedMatrices"
    LazyArraysBlockArraysExt = "BlockArrays"
    LazyArraysBlockBandedMatricesExt = "BlockBandedMatrices"
    LazyArraysStaticArraysExt = "StaticArrays"

    [deps.LazyArrays.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "56513a09b8e0ae6485f34401ea9e2f31357958ec"
uuid = "2d8b4e74-eb68-11e8-0fb9-d5eb67b50637"
version = "1.0.0"

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
git-tree-sha1 = "3acf07f130a76f87c041cfb2ff7d7284ca67b072"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.2+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2a7a12fc0a4e7fb773450d17975322aa77142106"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.2+0"

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "PrecompileTools", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "9f7253c0574b4b585c8909232adb890930da980a"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.6"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Printf"]
git-tree-sha1 = "738bdcacfef25b3a9e4a39c28613717a6b23751e"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.6.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "GPUArraysCore", "InteractiveUtils", "Krylov", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "OpenBLAS_jll", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLLogging", "SciMLOperators", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "1fdeafa25801ec2e0caa88fbe0afd7c41d3d087a"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "3.57.0"

    [deps.LinearSolve.extensions]
    LinearSolveAMDGPUExt = "AMDGPU"
    LinearSolveBLISExt = ["blis_jll", "LAPACK_jll"]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveCUSOLVERRFExt = ["CUSOLVERRF", "SparseArrays"]
    LinearSolveCliqueTreesExt = ["CliqueTrees", "SparseArrays"]
    LinearSolveEnzymeExt = ["EnzymeCore", "SparseArrays"]
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveFastLapackInterfaceExt = "FastLapackInterface"
    LinearSolveForwardDiffExt = "ForwardDiff"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolveMooncakeExt = "Mooncake"
    LinearSolvePardisoExt = ["Pardiso", "SparseArrays"]
    LinearSolveRecursiveFactorizationExt = "RecursiveFactorization"
    LinearSolveSparseArraysExt = "SparseArrays"
    LinearSolveSparspakExt = ["SparseArrays", "Sparspak"]

    [deps.LinearSolve.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    CUSOLVERRF = "a8cc9031-bad2-4722-94f5-40deabb4245c"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    FastLapackInterface = "29a986be-02c6-4525-aec4-84b980013641"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    LAPACK_jll = "51474c39-65e3-53ba-86ba-03b1b862ec14"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Sparspak = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
    blis_jll = "6136c539-28a5-5bf0-87cc-b183200dce32"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll"]
git-tree-sha1 = "8e6a74641caf3b84800f2ccd55dc7ab83893c10b"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.17.0+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

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

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "3bb3cf4685f1c90f22883f4c4bb6d203fa882b79"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "3.1.3"
weakdeps = ["BandedMatrices"]

    [deps.MatrixFactorizations.extensions]
    MatrixFactorizationsBandedMatricesExt = "BandedMatrices"

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
git-tree-sha1 = "53f817d3e84537d84545e0ad749e483412dd6b2a"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "d38b8653b1cdfac5a7da3b819c0a8d6024f9a18c"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.13"
weakdeps = ["ChainRulesCore"]

    [deps.MultivariatePolynomials.extensions]
    MultivariatePolynomialsChainRulesCoreExt = "ChainRulesCore"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "22df8573f8e7c593ac205455ca088989d0a2c7a0"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.7"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "FiniteDiff", "LinearAlgebra"]
git-tree-sha1 = "b3f76b463c7998473062992b246045e6961a074e"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "8.0.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "NonlinearSolveBase", "NonlinearSolveFirstOrder", "NonlinearSolveQuasiNewton", "NonlinearSolveSpectralMethods", "PrecompileTools", "Preferences", "Reexport", "SciMLBase", "SciMLLogging", "SimpleNonlinearSolve", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "5db62e7c18af752a7a5260ed7576e7429ca87be4"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "4.14.0"

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
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "Compat", "ConcreteStructs", "DifferentiationInterface", "EnzymeCore", "FastClosures", "LinearAlgebra", "Markdown", "MaybeInplace", "Preferences", "Printf", "RecursiveArrayTools", "SciMLBase", "SciMLJacobianOperators", "SciMLLogging", "SciMLOperators", "SciMLStructures", "Setfield", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "6613e7839ba583c63e9e13f8a310027086f66af6"
uuid = "be0214bd-f91f-a760-ac4e-3421ce2b2da0"
version = "2.11.1"

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
git-tree-sha1 = "df31d105d8e7254447256a44606f2a7e98b61aba"
uuid = "5959db7a-ea39-4486-b5fe-2dd0bf03d60d"
version = "1.11.1"

[[deps.NonlinearSolveQuasiNewton]]
deps = ["ArrayInterface", "CommonSolve", "ConcreteStructs", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLOperators", "StaticArraysCore"]
git-tree-sha1 = "ade27e8e9566b6cec63ee62f6a6650a11cf9a2eb"
uuid = "9a2c21bd-3a47-402d-9113-8faf9a0ee114"
version = "1.12.0"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveQuasiNewton.extensions]
    NonlinearSolveQuasiNewtonForwardDiffExt = "ForwardDiff"

[[deps.NonlinearSolveSpectralMethods]]
deps = ["CommonSolve", "ConcreteStructs", "LineSearch", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "eafd027b5cd768f19bb5de76c0e908a9065ddd36"
uuid = "26075421-4e9a-44e1-8bd1-420ed7ad02b2"
version = "1.6.0"
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

[[deps.OpenBLAS32_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "46cce8b42186882811da4ce1f4c7208b02deb716"
uuid = "656ef2d0-ae68-5445-9ca0-591084a874a2"
version = "0.3.30+0"

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

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqAdamsBashforthMoulton", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqExplicitRK", "OrdinaryDiffEqExponentialRK", "OrdinaryDiffEqExtrapolation", "OrdinaryDiffEqFIRK", "OrdinaryDiffEqFeagin", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqHighOrderRK", "OrdinaryDiffEqIMEXMultistep", "OrdinaryDiffEqLinear", "OrdinaryDiffEqLowOrderRK", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqNordsieck", "OrdinaryDiffEqPDIRK", "OrdinaryDiffEqPRK", "OrdinaryDiffEqQPRK", "OrdinaryDiffEqRKN", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqSSPRK", "OrdinaryDiffEqStabilizedIRK", "OrdinaryDiffEqStabilizedRK", "OrdinaryDiffEqSymplecticRK", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SparseArrays", "Static", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "dc70b19f140e9cbd02864c96ccfc833cbc6b8dd7"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.106.0"

[[deps.OrdinaryDiffEqAdamsBashforthMoulton]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqLowOrderRK", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "8307937159c3aeec5f19f4b661d82d96d25a3ff1"
uuid = "89bda076-bce5-4f1c-845f-551c83cdda9a"
version = "1.9.0"

[[deps.OrdinaryDiffEqBDF]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqSDIRK", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "156f2623ac97e7cf340848ba606f1226998980af"
uuid = "6ad6398a-0878-4a85-9266-38940aa047c8"
version = "1.14.0"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "ConcreteStructs", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "FastBroadcast", "FastClosures", "FastPower", "FillArrays", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLLogging", "SciMLOperators", "SciMLStructures", "Static", "StaticArrayInterface", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "862b950dd2a4d406b6e9cbb71f6e3b2043647665"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "3.2.0"

    [deps.OrdinaryDiffEqCore.extensions]
    OrdinaryDiffEqCoreEnzymeCoreExt = "EnzymeCore"
    OrdinaryDiffEqCoreMooncakeExt = "Mooncake"

    [deps.OrdinaryDiffEqCore.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"

[[deps.OrdinaryDiffEqDefault]]
deps = ["ADTypes", "DiffEqBase", "EnumX", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "PrecompileTools", "Preferences", "Reexport", "SciMLBase"]
git-tree-sha1 = "eef7a901d2a38462c13c12a39da7876c7b9b4a25"
uuid = "50262376-6c5a-4cf5-baba-aaf4f84d72d7"
version = "1.12.0"

[[deps.OrdinaryDiffEqDifferentiation]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqCore", "SciMLBase", "SciMLOperators", "SparseMatrixColorings", "StaticArrayInterface", "StaticArrays"]
git-tree-sha1 = "c3706545346a550a2669d8bcfe6db683af04a21c"
uuid = "4302a76b-040a-498a-8c04-15b101fed76b"
version = "1.22.0"
weakdeps = ["SparseArrays"]

    [deps.OrdinaryDiffEqDifferentiation.extensions]
    OrdinaryDiffEqDifferentiationSparseArraysExt = "SparseArrays"

[[deps.OrdinaryDiffEqExplicitRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "9d9b6bc309574d95acbf52e0f98a163f670e8dee"
uuid = "9286f039-9fbf-40e8-bf65-aa933bdc4db0"
version = "1.8.0"

[[deps.OrdinaryDiffEqExponentialRK]]
deps = ["ADTypes", "DiffEqBase", "ExponentialUtilities", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "65f2e40d7e9b1415c41838ec762777a4c36e4804"
uuid = "e0540318-69ee-4070-8777-9e2de6de23de"
version = "1.12.0"

[[deps.OrdinaryDiffEqExtrapolation]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastPower", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "e2f3ebd6cd7ed9c8d551fb10192644e8f6dd3cbb"
uuid = "becaefa8-8ca2-5cf9-886d-c06f3d2bd2c4"
version = "1.13.0"

[[deps.OrdinaryDiffEqFIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastGaussQuadrature", "FastPower", "LinearAlgebra", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "cbb6a36f09f1357a526c55a0a6805b60121eafb8"
uuid = "5960d6e9-dd7a-4743-88e7-cf307b64f125"
version = "1.20.0"

[[deps.OrdinaryDiffEqFeagin]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "b123f64a8635a712ceb037a7d2ffe2a1875325d3"
uuid = "101fe9f7-ebb6-4678-b671-3a81e7194747"
version = "1.8.0"

[[deps.OrdinaryDiffEqFunctionMap]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "cbd291508808caf10cf455f974c2025e886ed2a3"
uuid = "d3585ca7-f5d3-4ba6-8057-292ed1abd90f"
version = "1.9.0"

[[deps.OrdinaryDiffEqHighOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "9584dcc90cf10216de7aa0f2a1edc0f54d254cf6"
uuid = "d28bc4f8-55e1-4f49-af69-84c1a99f0f58"
version = "1.9.0"

[[deps.OrdinaryDiffEqIMEXMultistep]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "23602428114124a3e3df85fcbc5b461c79fb91bf"
uuid = "9f002381-b378-40b7-97a6-27a27c83f129"
version = "1.11.0"

[[deps.OrdinaryDiffEqLinear]]
deps = ["DiffEqBase", "ExponentialUtilities", "LinearAlgebra", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "c92913fa5942ed9bc748f3e79a5c693c8ec0c3d7"
uuid = "521117fe-8c41-49f8-b3b6-30780b3f0fb5"
version = "1.10.0"

[[deps.OrdinaryDiffEqLowOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "78223e34d4988070443465cd3f2bdc38d6bd14b0"
uuid = "1344f307-1e59-4825-a18e-ace9aa3fa4c6"
version = "1.10.0"

[[deps.OrdinaryDiffEqLowStorageRK]]
deps = ["Adapt", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "StaticArrays"]
git-tree-sha1 = "708c362418bd4503fd158f4f4e53151fbe57b46a"
uuid = "b0944070-b475-4768-8dec-fb6eb410534d"
version = "1.11.0"

[[deps.OrdinaryDiffEqNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "PreallocationTools", "RecursiveArrayTools", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "9f0be4bd586829a28a04c8f923598497f56ac226"
uuid = "127b3ac7-2247-4354-8eb6-78cf4e7c58e8"
version = "1.19.0"

[[deps.OrdinaryDiffEqNordsieck]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "05f3319c3bf1440897dc613194eb3db4d2d3e692"
uuid = "c9986a66-5c92-4813-8696-a7ec84c806c8"
version = "1.8.0"

[[deps.OrdinaryDiffEqPDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "7d63467f59f6504672ba93226f156f99c6095f60"
uuid = "5dd0a6cf-3d4b-4314-aa06-06d4e299bc89"
version = "1.10.0"

[[deps.OrdinaryDiffEqPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "Reexport", "SciMLBase"]
git-tree-sha1 = "baa77b7f874cda1f58f8c793fc7a9778e78a91c5"
uuid = "5b33eab2-c0f1-4480-b2c3-94bc1e80bda1"
version = "1.8.0"

[[deps.OrdinaryDiffEqQPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "9e351a8f923c843adb48945318437e051f6ee139"
uuid = "04162be5-8125-4266-98ed-640baecc6514"
version = "1.8.0"

[[deps.OrdinaryDiffEqRKN]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "98e22bd36729281743f77dd87a6036a9c3611370"
uuid = "af6ede74-add8-4cfd-b1df-9a4dbb109d7a"
version = "1.9.0"

[[deps.OrdinaryDiffEqRosenbrock]]
deps = ["ADTypes", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "e4605c3930703b5d38083ce1a998ee824dd67266"
uuid = "43230ef6-c299-4910-a778-202eb28ce4ce"
version = "1.22.0"

[[deps.OrdinaryDiffEqSDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "5d0a230f4e431e53af19502eaea8778f8f15edd4"
uuid = "2d112036-d095-4a1e-ab9a-08536f3ecdbf"
version = "1.11.0"

[[deps.OrdinaryDiffEqSSPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "StaticArrays"]
git-tree-sha1 = "8abc61382a0c6469aa9c3bff2d61c9925a088320"
uuid = "669c94d9-1f4b-4b64-b377-1aa079aa2388"
version = "1.11.0"

[[deps.OrdinaryDiffEqStabilizedIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqStabilizedRK", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "1719060baf014a3c1a6506113bc09d82a0903f0e"
uuid = "e3e12d00-db14-5390-b879-ac3dd2ef6296"
version = "1.10.0"

[[deps.OrdinaryDiffEqStabilizedRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "d156a972fa7bc37bf8377d33a7d51d152e354d4c"
uuid = "358294b1-0aab-51c3-aafe-ad5ab194a2ad"
version = "1.8.0"

[[deps.OrdinaryDiffEqSymplecticRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "9b783806fe2dc778649231cb3932cb71b63222d9"
uuid = "fa646aed-7ef9-47eb-84c4-9443fc8cbfa8"
version = "1.11.0"

[[deps.OrdinaryDiffEqTsit5]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "8be4cba85586cd2efa6c76d1792c548758610901"
uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a"
version = "1.9.0"

[[deps.OrdinaryDiffEqVerner]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "7c50b87bc8f00a38ef7acb8457828585f9ba4300"
uuid = "79d7bb75-1356-48c1-b8c0-6832512096c2"
version = "1.10.0"

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
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

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

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3ac7038a98ef6977d44adeadc73cc6f596c08109"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.79"

[[deps.PoissonRandom]]
deps = ["LogExpFunctions", "Random"]
git-tree-sha1 = "67afbcbe9e184d6729a92a022147ed4cf972ca7b"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.7"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "6f7cd22a802094d239824c57d94c8e2d0f7cfc7d"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.18"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

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

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "PrecompileTools"]
git-tree-sha1 = "c05b4c6325262152483a1ecb6c69846d2e01727b"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.34"

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
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "522f093a29b31a93e34eaea17ba055d850edea28"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.1"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "c5a07210bd060d6a8491b0ccdee2fa0235fc00bf"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.1.2"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

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

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "dbe5fd0b334694e905cb9fda73cd8554333c46e2"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.1"

[[deps.RandomNumbers]]
deps = ["Random"]
git-tree-sha1 = "c6ec94d2aaba1ab2ff983052cf6a606ca5985902"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.6.0"

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
git-tree-sha1 = "da082125fe8ca1b5d347c107c4a0ebb874734525"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.46.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsKernelAbstractionsExt = "KernelAbstractions"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStatisticsExt = "Statistics"
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTablesExt = ["Tables"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
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

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "2f609ec2295c452685d3142bc4df202686e555d2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.16"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLLogging", "SciMLOperators", "SciMLPublic", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "ba4214bef5dc42bb7dcba837300a558fdcd51ba0"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.135.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseDifferentiationInterfaceExt = "DifferentiationInterface"
    SciMLBaseDistributionsExt = "Distributions"
    SciMLBaseEnzymeExt = "Enzyme"
    SciMLBaseForwardDiffExt = "ForwardDiff"
    SciMLBaseMLStyleExt = "MLStyle"
    SciMLBaseMakieExt = "Makie"
    SciMLBaseMeasurementsExt = "Measurements"
    SciMLBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    SciMLBaseMooncakeExt = "Mooncake"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseReverseDiffExt = "ReverseDiff"
    SciMLBaseTrackerExt = "Tracker"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DifferentiationInterface = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    MLStyle = "d8e11817-5142-5d16-987a-aa16d5891078"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLJacobianOperators]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DifferentiationInterface", "FastClosures", "LinearAlgebra", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "e96d5e96debf7f80a50d0b976a13dea556ccfd3a"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.12"

[[deps.SciMLLogging]]
deps = ["Logging", "LoggingExtras", "Preferences"]
git-tree-sha1 = "7eebb9985e35b123e12025a3a2ad020cd6059f71"
uuid = "a6db7da4-7206-11f0-1eab-35f2a5dbe1d1"
version = "1.8.0"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "d1d14b15bbebf48dc80e8a7cfe640e2d835e22ea"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.14.1"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

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
git-tree-sha1 = "ebe7e59b37c400f694f52b58c93d26201387da70"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.9"

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
git-tree-sha1 = "315da09948861edbc6d18e066c08903487bb580d"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "2.10.0"

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
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

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

[[deps.SparseConnectivityTracer]]
deps = ["ADTypes", "DocStringExtensions", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "322365aa23098275562cbad6a1c2539ee40d9618"
uuid = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
version = "1.1.3"

    [deps.SparseConnectivityTracer.extensions]
    SparseConnectivityTracerChainRulesCoreExt = "ChainRulesCore"
    SparseConnectivityTracerLogExpFunctionsExt = "LogExpFunctions"
    SparseConnectivityTracerNNlibExt = "NNlib"
    SparseConnectivityTracerNaNMathExt = "NaNMath"
    SparseConnectivityTracerSpecialFunctionsExt = "SpecialFunctions"

    [deps.SparseConnectivityTracer.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "Random", "SparseArrays"]
git-tree-sha1 = "6ed48d9a3b22417c765dc273ae3e1e4de035e7c8"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.23"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsCUDAExt = "CUDA"
    SparseMatrixColoringsCliqueTreesExt = "CliqueTrees"
    SparseMatrixColoringsColorsExt = "Colors"
    SparseMatrixColoringsJuMPExt = ["JuMP", "MathOptInterface"]

    [deps.SparseMatrixColorings.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
    JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "49440414711eddc7227724ae6e570c7d5559a086"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.3.1"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eee1b9ad8b29ef0d936e3ec9838c7ec089620308"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.16"
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
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.SteadyStateDiffEq]]
deps = ["ConcreteStructs", "DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NonlinearSolveBase", "Reexport", "SciMLBase"]
git-tree-sha1 = "7b32737ebda77355ee61cfa9e59d376de3604629"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "2.9.0"

[[deps.StochasticDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FastPower", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLLogging", "SciMLOperators", "SimpleNonlinearSolve", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8158222b220f99334f57f98338003d02b2360990"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.92.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "83151ba8065a73f53ca2ae98bc7274d817aa30f2"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.8"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a3c1536470bf8c5e02096ad4853606d7c8f62721"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.2"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse32_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "libblastrampoline_jll"]
git-tree-sha1 = "1d43a4874b879f381b8a3a978f0ebe837cfd0922"
uuid = "ca45d3f4-326b-53b0-9957-23b75aacb3f2"
version = "7.12.1+0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.Sundials]]
deps = ["Accessors", "ArrayInterface", "CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "LinearSolve", "Logging", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll", "SymbolicIndexingInterface"]
git-tree-sha1 = "2d27edb89b7c555a57b8f22bfde92d6828d11cee"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "5.1.0"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS32_jll", "SuiteSparse32_jll"]
git-tree-sha1 = "a872f379c836e9cb5734485ca0681b192a59b98b"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "7.5.0+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "94c58884e013efff548002e8dc2fdd1cb74dfce5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.46"
weakdeps = ["PrettyTables"]

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils", "TermInterface"]
git-tree-sha1 = "5085671d2cba1eb02136a3d6661c583e801984c1"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "1.1.0"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "EnumX", "ExproniconLite", "LinearAlgebra", "MacroTools", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "ReadOnlyArrays", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "WeakCacheSets"]
git-tree-sha1 = "669af5a0ead11fe4f310f2e3778a24769dae20ec"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "4.14.0"

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
git-tree-sha1 = "53d84d58e5fd1dfab5446d2b31eb1153bcc9f554"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "7.9.1"

    [deps.Symbolics.extensions]
    SymbolicsD3TreesExt = "D3Trees"
    SymbolicsDistributionsExt = "Distributions"
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLatexifyExt = ["Latexify", "LaTeXStrings"]
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"
    SymbolicsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Symbolics.weakdeps]
    D3Trees = "e3df1716-f71e-5df9-9e2d-98e193103c45"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
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

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "d969183d3d244b6c33796b5ed01ab97328f2db85"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.5"

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
git-tree-sha1 = "c25751629f5baaa27fef307f96536db62e1d754e"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.27.0"

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
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

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
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

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
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

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
git-tree-sha1 = "6ab498eaf50e0495f89e7a5b582816e2efb95f64"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.54+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

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
# ╠═fed88caa-1520-41f7-adb3-785e5c9529c6
# ╟─734b3b08-061e-4f93-8574-468d824815da
# ╟─800dc762-ce0b-463f-858f-6e8eabbc26b0
# ╟─b842e98e-34e2-40f2-84b6-c180815c2df3
# ╟─3c0c1fea-1ae6-4c9d-9db1-764937cd5fbc
# ╟─99d0fdc9-e368-462c-a357-86f07624a52e
# ╟─6d6cbebb-bb3b-437d-9835-0a79f36857f2
# ╟─71b94af7-76ee-4228-987c-2f22e0951552
# ╟─22c37732-1cf7-4d80-a250-9fd5e4a2f88c
# ╟─8da87629-09ff-4677-b5e8-2a3ff5777a8e
# ╟─c5e21675-f120-4555-84be-d99b35592f2d
# ╟─d12b4a5f-a6ad-4c91-8795-fe8f93f5c93d
# ╟─c7eee8bb-b03d-4679-a51e-36fc629294e1
# ╟─700015f6-aac0-40c5-a455-8d98c1227049
# ╠═207305eb-4496-4518-a5bc-be85173314a5
# ╟─a6fdec29-3438-4ebc-af81-7a3baf0175ae
# ╟─ee7baf70-c7af-4791-8277-fadc0961e0e7
# ╟─72416f96-8a0a-4c13-97d1-2990e49c4cb0
# ╠═0baae21f-71b5-42fd-be33-31de085928d3
# ╟─a7c5650c-9bbc-44a1-8d51-d473864faae7
# ╠═49d1a7f7-2bf2-4472-94df-6247b9237ddd
# ╟─4c975b71-a387-49cb-90d9-fc51acefc795
# ╟─d431eb6d-429a-49d6-83b8-b6f21831ec86
# ╟─98610902-dfff-4fd1-8ffa-8484405aee20
# ╟─2d6fcef9-df4b-4eec-be5c-a8865c3a1b76
# ╠═568d9fe3-6716-4a9a-ba1d-b9e6fd039150
# ╟─6cab6cb7-a432-40b6-9390-ad0083fe486d
# ╠═38e91d0d-1020-44a8-becc-8542fd600104
# ╟─4b81f302-9637-4161-b745-8ac39b9e31d3
# ╟─f525b379-b5e2-48bc-af1a-e57436f6176e
# ╠═c01cf078-c0cb-4000-bf76-a90cbeb10b89
# ╟─46c5bd5a-21b1-4f92-8eb1-a11ec2a0c94a
# ╠═9c5b30d5-08aa-48a8-9ae2-c3b6c432ab89
# ╟─ce2383dc-c34f-4b56-8501-42ce0539c95c
# ╠═a3d1e1bf-c513-4d6b-a43b-3dab0106f1a5
# ╟─5680e62b-973b-4a60-bb3f-8785ce07e581
# ╠═47bf94da-1368-41bd-ba46-c9c1f75cf44e
# ╟─9af6ce74-5678-4b48-8758-37b0d5a6f0e4
# ╠═a9c6a292-086e-4aa0-9856-78d3ab3fbe35
# ╟─5890b699-b7de-47a3-bee7-1e7dd7663fbe
# ╠═ef79392b-395c-497a-9c0e-dc2cd468f6e1
# ╟─8c9ab125-2acb-4732-a9bf-7838e819e4f7
# ╠═e74e6af3-7c28-4836-a7a7-db65cca302dc
# ╠═2a39d6f8-da49-4976-9aa7-889391e55a5d
# ╟─f4399d4f-9e85-4974-9f39-d47acaf1399c
# ╠═6a59ed2e-f040-4e10-bc55-91d2c1dcc97e
# ╟─14ae7f11-1065-46aa-a6ed-eea400c0d2ec
# ╠═c8efcb3f-cc6f-4c45-a0f3-55cdb73d5195
# ╟─4d09ed45-423c-4bd6-802d-59389a966d2e
# ╠═b9d7da8b-5f01-4f4f-8929-fd44f9689f3e
# ╠═e527f834-9965-4cdf-bda8-b4ea07e86f27
# ╠═c869e896-cc87-4b12-bba2-84b0cf0964c6
# ╠═2c2f97cc-127a-446e-a6cd-e9651df833f0
# ╟─9a756404-f35b-40b9-bc88-4b7e3a7df0b6
# ╠═805ef471-1bd2-4ce2-bf50-eeb9c6b10e4d
# ╟─6eef3827-1af6-431d-a163-67bb595c337d
# ╠═e14648ff-69f6-47c2-924c-c5ac69e200ed
# ╟─aaba63dd-8fc3-405d-ad19-8e8368e70019
# ╠═03754250-3831-41db-ac11-fe6f5a08f0c1
# ╟─2d27f581-2c53-409e-afde-812e70bba8bf
# ╟─6ddf393c-9f37-43b2-b401-9ab0e7a9e88c
# ╠═5ce798b5-4366-4cc1-8d2d-881a827d6dc9
# ╟─7df63fdf-ef41-44f2-8a55-e5c2c849029c
# ╟─64c2a3ee-e7a0-4e03-b9a8-86239e1ca81e
# ╟─0895d464-a029-410d-9e7d-89cfac2d1615
# ╟─07cd9aad-029f-42c6-abe8-ab4a9a2a910c
# ╟─0bdf9dbf-479c-46f6-bd86-50576095cba0
# ╟─6fe43e3a-2e8f-4708-a3ec-6f5a8088060e
# ╠═f863d68f-590e-4b96-8433-dc6b5177539f
# ╟─f5a983bf-ef3a-4d5d-928d-da1097b91ee8
# ╠═bd8743d6-8f21-413d-835a-e543926baa09
# ╟─9ab0a10b-8165-401a-a2b6-c9726526a906
# ╠═1b8af600-56eb-4508-bc52-aa4e138b4c7e
# ╠═4fcf6e4d-b9e7-43d1-badc-d8d8afc004b3
# ╠═80005099-7154-4306-9172-c9a168336e14
# ╠═c291700e-3a84-49a7-85d6-592cfb3b1a11
# ╟─1ec99792-905d-4e1b-a413-ef58143d3c68
# ╠═4e1140be-3e61-47bd-9f9f-6b5dfbff6ae2
# ╟─2afd459c-90e9-4105-9121-27e21bb89eeb
# ╠═4173e91c-c01a-4fa2-ae77-0408bf7d9a1b
# ╟─24920a1b-355c-43e1-aff2-028e6f25584c
# ╟─8fd489b2-620a-40cd-8283-dc88b7f3584f
# ╟─59a342e1-930d-40c4-86fc-11ba334d7103
# ╠═a2a12511-97d6-4e83-a21d-1c7528f486bc
# ╠═8ade720c-37b7-4436-8bcd-8e902f52057b
# ╟─b39e416e-a204-4f6b-a6e4-00f6a71d3b27
# ╟─c414f207-0d4f-4cc4-965c-9b68b5a85400
# ╠═ac57dfbe-f19b-4c5f-a753-5f4162e245d4
# ╠═b1e67fcd-f70b-4dd2-87e5-74bef1fa201e
# ╠═d1f3d076-b9d6-4b92-b25c-4155b5fc464b
# ╟─f336b083-99c7-4345-8ada-78b166a98abe
# ╟─74c82c98-f233-4e72-8f74-17bdfeddf884
# ╟─ab1b7695-3cd6-4e67-9cfd-fbaf2f4d1d15
# ╟─6475b25b-8711-44cd-bc3b-e3d42681bb93
# ╟─ea4e58e9-d041-4a6e-b0d8-83e3aef7648b
# ╠═d7cc8a66-220a-4e9b-b42e-a0e085ed3a0f
# ╟─477c8f59-97d1-405d-a1c8-64b7e0b9119f
# ╟─08e05c65-06a2-4560-92eb-b014dc7c3d70
# ╟─2235689d-9c83-4907-aa17-c2624fbeb68d
# ╟─df8a9449-851c-4546-97a7-7fd4a270a867
# ╠═e8c51ae9-ded4-450a-80f9-c484c2b991d1
# ╠═c96ea0ad-47c3-469e-a2f7-0743e681a6d1
# ╟─0813f6a3-aadc-491c-ad1e-ec54dcbd0d56
# ╟─de4ba100-5f58-4657-aea0-a3f31eacda65
# ╠═77ca6eab-c5bf-479c-8982-1483933bbb6e
# ╟─d4a1ec85-e6f5-48ed-9724-202c4dae9ae9
# ╠═7f8ebe64-ea86-4eb8-9939-6682180b9dd2
# ╟─aec6e4fc-e496-4add-b982-ab60f9f900a0
# ╠═8940cb2b-2c8c-407e-bc6b-426843cf6125
# ╠═28955b95-df19-403a-bf79-b68e9be8e1dd
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
