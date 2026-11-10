### A Pluto.jl notebook ###
# v1.1.0

using Markdown
using InteractiveUtils

# ╔═╡ fed88caa-1520-41f7-adb3-785e5c9529c6
using DataFrames, ChaosTools, CSV, DifferentialEquations, Interpolations, LinearAlgebra, OrdinaryDiffEqBDF, PlutoUI, QuadGK, Symbolics, TikzPictures, Unitful, UnitfulAstro

# ╔═╡ 734b3b08-061e-4f93-8574-468d824815da
# ╠═╡ skip_as_script = true
#=╠═╡
TableOfContents(title="🌌 SF model", depth=3)
  ╠═╡ =#

# ╔═╡ 313c5c60-1c11-4fa3-a71f-16df0d2d0ebc
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Introduction

In the standard cosmological framework, galaxies form as baryonic gas cools and collapses within dark matter haloes, growing through accretion and mergers ([Rees1977](https://doi.org/10.1093/mnras/179.4.541), [White1978](https://doi.org/10.1093/mnras/183.3.341), [Blumenthal1984](https://doi.org/10.1038/311517a0)). Their evolution is shaped by the interplay of gravitational collapse, gas cooling, star formation (SF), chemical enrichment, and feedback ([Efstathiou1992](https://doi.org/10.1093/mnras/256.1.43P), [White1991](https://doi.org/10.1086/170483)). Cosmological simulations capture these processes self-consistently and now reproduce systems broadly consistent with observations, including Milky Way analogues ([Grand2017](https://doi.org/10.1093/mnras/stx071), [Vogelsberger2020](https://doi.org/10.1038/s42254-019-0127-2)). However, SF and stellar feedback occur below the resolution of these simulations and must be described through sub-grid models.

Early SF prescriptions formed stars from cold, dense gas at a rate proportional to the gas density divided by a characteristic timescale ([Katz1992](https://doi.org/10.1086/171366), [Steinmetz1994](https://doi.org/10.48550/arXiv.astro-ph/9312010)). With effective feedback and a calibrated star formation efficiency, this reproduces the observed Kennicutt-Schmidt relation ([Schmidt1959](https://doi.org/10.1086/146614), [Kennicutt1998](https://doi.org/10.1086/305588), [Scannapieco2006](https://doi.org/10.1111/j.1365-2966.2006.10785.x)). Observations, however, show that the SFR correlates more tightly with molecular than with total gas, on both resolved and galaxy-wide scales ([Leroy2008](https://doi.org/10.1088/0004-6256/136/6/2782), [Bigiel2008](https://doi.org/10.1088/0004-6256/136/6/2846), [Baker2021](https://doi.org/10.1093/mnras/stab3672), [Baker2022](https://doi.org/10.1093/mnras/stac3413)), and that SF timescales depend on local ISM conditions ([Leroy2013](https://doi.org/10.1088/0004-6256/146/2/19)). Stellar feedback is even harder to model. Different implementations yield widely different galaxies from identical initial conditions ([Scannapieco2012](https://doi.org/10.1111/j.1365-2966.2012.20993.x)), and its parameters are model-specific and not directly constrained by observations, although it remains the dominant feedback channel up to Milky Way masses ([Naab2017](https://doi.org/10.1146/annurev-astro-081913-040019)).

Because the ISM sets the conditions for SF, considerable effort has gone into sub-grid descriptions of its multiphase structure, from two-phase models of hot and cold gas in pressure equilibrium ([McKee1977](https://doi.org/10.1086/155667), [Yepes1997](https://doi.org/10.1093/mnras/284.1.235), [Springel2003](https://doi.org/10.1046/j.1365-8711.2003.06206.x)) to models that explicitly track molecular hydrogen and adopt $\mathrm{H}_2$-based SF laws. The latter compute $\mathrm{H}_2$ fractions using time-dependent or equilibrium chemistry ([Pelupessy2006](https://doi.org/10.1086/504366), [Robertson2008](https://doi.org/10.1086/587796), [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55)), pressure-based prescriptions ([Blitz2006](https://doi.org/10.1086/505417), [Murante2014](https://doi.org/10.1093/mnras/stu2400)), or theoretical models based on column density and metallicity ([Krumholz2009](https://doi.org/10.1088/0004-637X/699/1/850), [Hopkins2014](https://doi.org/10.1093/mnras/stu1738)), with the latter producing a clumpier, more realistic ISM ([Valentini2022](https://doi.org/10.1093/mnras/stac2110)). A selection of such simulations is listed below.

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

Our simulations use the moving-mesh magnetohydrodynamics code $\texttt{AREPO}$ ([Springel2010](https://doi.org/10.1111/j.1365-2966.2009.15715.x)), which solves the MHD equations on a Voronoi mesh that moves with the flow and is regenerated every time step, combining the adaptive resolution of Lagrangian methods with the accuracy of finite-volume schemes ([Pakmor2015](https://doi.org/10.1093/mnras/stv2380)). Gravity is computed with a TreePM scheme ([Springel2005](https://doi.org/10.1111/j.1365-2966.2005.09655.x)). We adopt the Auriga galaxy formation model, which includes primordial and metal-line cooling, a uniform UV background, magnetic fields, Type II and Type Ia supernova feedback, and mass and metal return from AGB stars ([Vogelsberger2013](https://doi.org/10.1093/mnras/stt1789), [Pakmor2017](https://doi.org/10.1093/mnras/stx1074), [Grand2017](https://doi.org/10.1093/mnras/stx071)).
"""
  ╠═╡ =#

# ╔═╡ b842e98e-34e2-40f2-84b6-c180815c2df3
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Star formation model

## Phases and notation

We modeled the interstellar medium as a multiphase structure made up of six components. Three correspond to different phases of hydrogen: ionized ($\sim 10^4 \, \mathrm{K}$), atomic ($\sim 100 \, \mathrm{K}$), and molecular ($\sim 10 \, \mathrm{K}$) gas. The other components are: metals, stars, and dust.

For each mass-transfer reaction between phases, only the dominant pathway is considered. Although the ISM consists primarily of hydrogen and helium, we model only hydrogen-related reactions, as they govern the key processes involved. These include the photoionization of atoms, the recombination of electrons with ions, the formation of molecular hydrogen from atomic hydrogen, and its destruction via photodissociation by ultraviolet radiation. Additionally, we account for the production of ionized gas, metals, and dust by supernovae, the formation and destruction of dust from and into metals, and the formation of stars from molecular gas.

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

It is important to clarify that these fractions describe the internal composition of each gas cell within our subgrid model. The stellar mass fraction, $f_s$, represents the mass of stars formed within a gas cell that has not yet been converted into a distinct star particle. It is used to compute the cell's instantaneous SFR, which determines if a new star particle will be generated.

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

In the context of our model, $V_\mathrm{cell}$, $M_\mathrm{cell}$, $\rho_\mathrm{cell}$, $V_j$, $m_j$, and $x_j$ are constants for a given gas cell.

In the same way, we can write a relation for the mass density of the $j$ component,

$\begin{equation}
    \rho_j = \frac{M_j}{V_j} = \frac{M_j}{M_\mathrm{cell}} \, \frac{M_\mathrm{cell}}{V_j} = f_j \,  \frac{M_\mathrm{cell}}{x_j \, V_\mathrm{cell}} = f_j \, \frac{\rho_\mathrm{cell}}{x_j} \, .
\end{equation}$

So, using these relations we can write any differential equation for the quantities $M_j$, $\rho_j$, and $n_j$ as equations for $f_j$.

So far, we have made only two assumptions: first, that the ISM is made solely of the six components previously described; and second, that $V_\mathrm{cell}$, $M_\mathrm{cell}$, $\rho_\mathrm{cell}$, $V_j$, $m_j$, and $x_j$ are constants.

For simplicity we will also adopt $x_i = x_a = x_m = 1$. This implies that the ionized, atomic, and molecular hydrogen phases are well mixed and collectively fill the entire cell volume.
"""
  ╠═╡ =#

# ╔═╡ 99d0fdc9-e368-462c-a357-86f07624a52e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Overview of the physical processes

The diagrams below summarize the model. Each color identifies a physical process, and it is the same in the three figures. Each arrow is a physical process, which moves mass from one component to another.
"""
  ╠═╡ =#

# ╔═╡ 8febbbdb-b31c-4232-8a8b-15ef6a87557b
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	#################################################################################
	# LaTeX preamble shared by the TikZ diagrams of the model
	#
	# TIKZ_COLORS: One color per physical process, chosen for a dark background
	# TIKZ_STYLES: Styles of the diagrams of the physical processes
	#################################################################################

	const TIKZ_COLORS = raw"""
		\usepackage{amsmath}
		% Colors for a dark background, one per physical process
		% (lighter versions of the colorblind-safe palette of Okabe & Ito 2008)
		\definecolor{c_rec}{HTML}{E58FC9}   % Recombination
		\definecolor{c_ion}{HTML}{56B4E9}   % Photoionization
		\definecolor{c_cond}{HTML}{E69F00}  % Condensation
		\definecolor{c_diss}{HTML}{3DD6A3}  % Photodissociation
		\definecolor{c_sf}{HTML}{FF7B72}    % Star formation
		\definecolor{c_star}{HTML}{F0E442}  % Stellar yields
		\definecolor{c_dg}{HTML}{B5E37A}    % Dust growth
		\definecolor{c_dd}{HTML}{B4A7F5}    % Dust destruction
	"""

	const TIKZ_STYLES = raw"""
		\usetikzlibrary{arrows.meta, positioning}
		\tikzset{
		    >={Stealth[length=2.4mm]},
		    box/.style={
		        rectangle,
		        rounded corners=3pt,
		        draw=white,
		        thick,
		        text=white,
		        minimum width=3.4em,
		        minimum height=2.2em,
		        align=center,
		    },
		    flow/.style={
		        ->,
		        thick,
		    },
		}
	"""
end;
  ╠═╡ =#

# ╔═╡ 71b94af7-76ee-4228-987c-2f22e0951552
# ╠═╡ skip_as_script = true
#=╠═╡
TikzPictures.TikzPicture(
	L"""
		% Phases (the same positions in both diagrams)
		\node[box] (atom) {H\,\textsc{i}};
		\node[box, left=2.6cm of atom] (ion) {H\,\textsc{ii}};
		\node[box, right=2.6cm of atom] (molecule) {H$_2$};
		\node[box, below=2.6cm of atom] (stars) {Stars};
		\node[box, below left=1.1cm and 0.6cm of stars] (metals) {Metals};
		\node[box, below right=1.1cm and 0.6cm of stars] (dust) {Dust};

		% Processes, one arrow per direction
		\draw[flow, c_rec]  (ion) to[bend left=30] node[above] {recombination} (atom);
		\draw[flow, c_ion]  (atom) to[bend left=30] node[below] {photoionization} (ion);
		\draw[flow, c_cond] (atom) to[bend left=30] node[above] {condensation} (molecule);
		\draw[flow, c_diss] (molecule) to[bend left=30] node[below] {photodissociation} (atom);
		\draw[flow, c_sf]   (molecule) to[bend left=25] node[below right] {star formation} (stars);
		\draw[flow, c_star] (stars) to[bend left=25] node[below left] {recycling} (ion);
		\draw[flow, c_star] (stars) to[bend right=20] node[above left] {enrichment} (metals);
		\draw[flow, c_star] (stars) to[bend left=20] node[above right] {dust production} (dust);
		\draw[flow, c_dg]   (metals) to[bend left=20] node[above] {dust growth} (dust);
		\draw[flow, c_dd]   (dust) to[bend left=20] node[below] {dust destruction} (metals);
	""",
	width="75em",
	preamble=TIKZ_COLORS * TIKZ_STYLES,
)
  ╠═╡ =#

# ╔═╡ d12b4a5f-a6ad-4c91-8795-fe8f93f5c93d
# ╠═╡ skip_as_script = true
#=╠═╡
TikzPictures.TikzPicture(
	L"""
		\node[text=white] {
		${\boldmath
		\begin{aligned}
		    \dv{}{t}f_i &= \textcolor{c_ion}{e_\mathrm{ion} \, \eta_\mathrm{ion}^* \, \psi} + \textcolor{c_ion}{e_\mathrm{ion} \, \Gamma_\mathrm{UVB} \, f_a} + \textcolor{c_star}{R \, (1 - Z_\mathrm{SN}) \, \psi} - \textcolor{c_rec}{\frac{f_i}{\tau_\mathrm{rec}}} \, , \\
		    \dv{}{t}f_a &= \textcolor{c_rec}{\frac{f_i}{\tau_\mathrm{rec}}} + \textcolor{c_diss}{e_\mathrm{diss} \, \eta_\mathrm{diss}^* \, \psi} + \textcolor{c_diss}{e_\mathrm{diss} \, \Gamma_\mathrm{LWB} \, f_m} - \textcolor{c_ion}{e_\mathrm{ion} \, \eta_\mathrm{ion}^* \, \psi} - \textcolor{c_ion}{e_\mathrm{ion} \, \Gamma_\mathrm{UVB} \, f_a} - \textcolor{c_cond}{\frac{f_a}{\tau_\mathrm{cond}}} \, , \\
		    \dv{}{t}f_m &= \textcolor{c_cond}{\frac{f_a}{\tau_\mathrm{cond}}} - \textcolor{c_diss}{e_\mathrm{diss} \, \eta_\mathrm{diss}^* \, \psi} - \textcolor{c_diss}{e_\mathrm{diss} \, \Gamma_\mathrm{LWB} \, f_m} - \textcolor{c_sf}{\psi} \, , \\
		    \dv{}{t}f_s &= \textcolor{c_sf}{\psi} - \textcolor{c_star}{R \, \psi} \, , \\
		    \dv{}{t}f_Z &= \textcolor{c_star}{R \, Z_\mathrm{SN} \, (1 - \delta_\mathrm{d}) \, \psi} + \textcolor{c_dd}{C_\mathrm{dd} \, f_d \, \psi} - \textcolor{c_dg}{\xi \, \frac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m)} \, , \\
		    \dv{}{t}f_d &= \textcolor{c_star}{R \, Z_\mathrm{SN} \, \delta_\mathrm{d} \, \psi} + \textcolor{c_dg}{\xi \, \frac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m)} - \textcolor{c_dd}{C_\mathrm{dd} \, f_d \, \psi} \, .
		\end{aligned}}$
		};
	""",
	width="70em",
	preamble=TIKZ_COLORS * raw"""
		\usepackage{physics}
		\setlength{\jot}{10pt}
	""",
)
  ╠═╡ =#

# ╔═╡ c7eee8bb-b03d-4679-a51e-36fc629294e1
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
We can explicitly check for mass conservation

$\begin{equation}
    \frac{\mathrm{d}}{\mathrm{d}t}(f_i + f_a + f_m + f_s + f_Z + f_d) = 0 \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 8da87629-09ff-4677-b5e8-2a3ff5777a8e
# ╠═╡ skip_as_script = true
#=╠═╡
TikzPictures.TikzPicture(
	L"""
		% Phases (the same positions in both diagrams)
		\node[box] (atom) {H\,\textsc{i}};
		\node[box, left=4cm of atom] (ion) {H\,\textsc{ii}};
		\node[box, right=4cm of atom] (molecule) {H$_2$};
		\node[box, below=2.9cm of atom] (stars) {Stars};
		\node[box, below left=1.5cm and 1.3cm of stars] (metals) {Metals};
		\node[box, below right=1.5cm and 1.3cm of stars] (dust) {Dust};

		% Net flows: positive in the direction of the arrow
		\draw[flow, white] (ion) -- node[above=20pt, align=center] {
		    $\textcolor{c_rec}{\dfrac{f_i}{\tau_\mathrm{rec}}}
		    - \textcolor{c_ion}{e_\mathrm{ion} \left(\eta_\mathrm{ion}^* \, \psi + \Gamma_\mathrm{UVB} \, f_a\right)}$
		} (atom);
		\draw[flow, white] (atom) -- node[above=20pt, align=center] {
		    $\textcolor{c_cond}{\dfrac{f_a}{\tau_\mathrm{cond}}}
		    - \textcolor{c_diss}{e_\mathrm{diss} \left(\eta_\mathrm{diss}^* \, \psi + \Gamma_\mathrm{LWB} \, f_m\right)}$
		} (molecule);
		\draw[flow, white] (molecule) to[bend left=25] node[below right] {$\textcolor{c_sf}{\psi}$} (stars);
		\draw[flow, white] (stars) to[bend left=25] node[below left] {$\textcolor{c_star}{R \, (1 - Z_\mathrm{SN}) \, \psi}$} (ion);
		\draw[flow, white] (stars) to[bend right=20] node[above left] {$\textcolor{c_star}{R \, Z_\mathrm{SN} \, (1 - \delta_\mathrm{d}) \, \psi}$} (metals);
		\draw[flow, white] (stars) to[bend left=20] node[above right] {$\textcolor{c_star}{R \, Z_\mathrm{SN} \, \delta_\mathrm{d} \, \psi}$} (dust);
		\draw[flow, white] (metals) -- node[below=8pt] {
		    $\textcolor{c_dg}{\xi \, \dfrac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m)}
		    - \textcolor{c_dd}{C_\mathrm{dd} \, f_d \, \psi}$
		} (dust);
	""",
	width="75em",
	preamble=TIKZ_COLORS * TIKZ_STYLES,
)
  ╠═╡ =#

# ╔═╡ 700015f6-aac0-40c5-a455-8d98c1227049
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Units

We have the freedom to choose three independent units in this model. Time, mass, and length. The choice is reflected in the constants $C_\mathrm{star}$, $C_\mathrm{rec}$, $C_\mathrm{cond}$, $C_\mathrm{acc,Si}$, $C_\mathrm{acc,C}$, $C_\mathrm{sd}$, $C_\mathrm{sh2}$, $C_\mathrm{lsh}$, $C_\mathrm{\tau ion}$, and $C_\mathrm{\tau diss}$.

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

*  $\eta_\mathrm{diss}^*$: Rate of molecular gas dissociation by stellar radiation per unit mass of newly formed stars.
*  $\eta_\mathrm{ion}^*$: Rate of atomic gas ionization by stellar radiation per unit mass of newly formed stars.
*  $R$: Fraction of a stellar population’s mass that is returned to the ISM as ionized gas and metals.
*  $Z_\mathrm{SN}$: Fraction of the returned mass that consists of metals.
*  $\Gamma_\mathrm{UVB}$: Photoionization rate due to the ultraviolet background (UVB).
*  $\Gamma_\mathrm{LWB}$: Photodissociation rate due to the Lyman-Werner background (LWB).

In contrast, the parameters $C_\mathrm{star}$, $C_\mathrm{rec}$, $C_\mathrm{cond}$, $C_\mathrm{acc,Si}$, $C_\mathrm{acc,C}$, $C_\mathrm{sd}$, $C_\mathrm{sh2}$, $C_\mathrm{lsh}$, $C_\mathrm{\tau ion}$, $C_\mathrm{\tau diss}$, and $C_\mathrm{dd}$ do not have physical meaning; they are included solely for convenience.
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
assuming only a dependence on the mass of molecular gas, as in [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635), instead of the total amount of gas

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

In the literature, a value of $\epsilon_\text{ff} \approx 0.01$ is commonly used ([Krumholz2007](https://doi.org/10.1086/509101), [Krumholz2014](https://doi.org/10.1016/j.physrep.2014.02.001)), although there is still significant uncertainty regarding its true value ([Lee2016](https://doi.org/10.3847/1538-4357/833/2/229), [Utomo2018](https://doi.org/10.3847/2041-8213/aacf8f)). In our model, we adopt $\epsilon_\text{ff} = 1.0$ because star formation is self-regulated through the coupled evolution of gas phases, as described by the system of ODEs. This framework naturally determines the effective star formation rate. We note, however, that previous studies have shown this parameter has little impact on the global properties of simulated galaxies ([Li2018](https://doi.org/10.3847/1538-4357/aac9b8), [Brown2022](https://doi.org/10.1093/mnras/stac1164)).

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

Using $n_a = f_a \, \rho_\mathrm{cell} / m_p$, we can write

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

# ╔═╡ 80bae1ac-9f61-4512-84c9-c562ad7a61ba
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Condensation

The rate of molecular hydrogen formation on the surface of dust grains is ([Hollenbach1971a](https://doi.org/10.1086/150754), [Hollenbach1971b](https://doi.org/10.1086/150755))

$\begin{equation}
    \frac{\mathrm{d}}{\mathrm{d}t} n_m \biggr\rvert_\mathrm{condensation} = \frac{1}{2} \, \gamma \, \langle v_a \rangle \, n_d \, n_a \, \langle \sigma_d \rangle \, ,
\end{equation}$

where $\langle v_a \rangle$ is the average thermal velocity of hydrogen atoms, $\langle \sigma_d \rangle$ is the mean geometric cross-section of a dust grain, $n_d$ and $n_a$ are the number densities of dust grains and atomic hydrogen respectively, and $\gamma$ is the recombination coefficient, defined as the fraction of atoms striking a grain that eventually form a molecule. As far as we know, most of the early works only use this expression in equilibrium equations ([Hollenbach1971a](https://doi.org/10.1086/150754), [Hollenbach1971b](https://doi.org/10.1086/150755), [Jura1974](https://doi.org/10.1086/152975), [Jura1975](https://doi.org/10.1086/153545), [Black1987](https://doi.org/10.1086/165740), [Sternberg1988](https://doi.org/10.1086/166664), [Goldshmidt1995](https://doi.org/10.1086/175168)), and it first appears in a differential equation that does not assume equilibrium in [Draine1996](https://doi.org/10.1086/177689).

Following [Hollenbach1971b](https://doi.org/10.1086/150755) (see Section II), we use

$\begin{equation}
    \gamma \, \langle v_a \rangle = 5 \times 10^4 \, \mathrm{cm \, s^{-1}} \, ,
\end{equation}$

which they find to be valid within a factor of 2 for gas temperatures between $4$ and $200 \, \mathrm{K}$.

To relate the cross-section of the grains to the dust content of the cell, we assume spherical grains of material density $s$, following a fixed size distribution with moments $\langle a^2 \rangle$ and $\langle a^3 \rangle$. Then, the total grain cross-section and the dust mass density are

$\begin{align}
    n_d \, \langle \sigma_d \rangle &= \pi \, \langle a^2 \rangle \, n_d \, , \\
    f_d \, \rho_\mathrm{cell} &= \frac{4}{3} \, \pi \, s \, \langle a^3 \rangle \, n_d \, ,
\end{align}$

so

$\begin{equation}
    n_d \, \langle \sigma_d \rangle = \frac{3}{4} \, \frac{f_d \, \rho_\mathrm{cell}}{s \, a_\mathrm{eff}} \, ,
\end{equation}$

where $a_\mathrm{eff} = \langle a^3 \rangle / \langle a^2 \rangle$ is the effective grain radius. We use the values

$\begin{align}
    a_\mathrm{eff} &= \frac{\langle a^3 \rangle_0}{\langle a^2 \rangle_0} = 0.0158 \, \mathrm{μm} \, , \\
    s &= 3 \, \mathrm{g \, cm^{-3}} \, ,
\end{align}$

where $a_\mathrm{eff}$ comes from models C and F of [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x) (Table 2), the same size distribution we use for dust growth, and $s$ is taken from footnote 5 of [Ragone-Figueroa2024](https://doi.org/10.1051/0004-6361/202451344), which lies between the values for silicate ($3.3 \, \mathrm{g \, cm^{-3}}$) and graphite ($2.26 \, \mathrm{g \, cm^{-3}}$) in Table 1 of [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x). Notice that the number-weighted mean radius of the same distribution, $\langle a \rangle = 0.00167 \, \mathrm{μm}$, is not the relevant quantity here, as it is dominated by the smallest grains.

Using $f_m = 2 \, m_p \, n_m / \rho_\mathrm{cell}$ and $n_a = f_a \, \rho_\mathrm{cell} / m_p$, the condensation rate becomes

$\begin{align}
    \frac{\mathrm{d}}{\mathrm{d}t} f_m \biggr\rvert_\mathrm{condensation} &= \frac{2 \, m_p}{\rho_\mathrm{cell}} \, \frac{1}{2} \, \gamma \, \langle v_a \rangle \, n_d \, n_a \, \langle \sigma_d \rangle \\
	&= \frac{m_p}{\rho_\mathrm{cell}} \, \gamma \, \langle v_a \rangle \, f_a \, \frac{\rho_\mathrm{cell}}{m_p} \, \frac{3}{4} \, \frac{f_d \, \rho_\mathrm{cell}}{s \, a_\mathrm{eff}} \\
	&= \frac{3}{4} \, \frac{\gamma \, \langle v_a \rangle \, \rho_\mathrm{cell}}{s \, a_\mathrm{eff}} \, f_d \, f_a \\
	&= \frac{f_a}{\tau_\mathrm{cond}} \, ,
\end{align}$

with the time scale

$\begin{equation}
    \tau_\mathrm{cond} = \frac{4 \, s \, a_\mathrm{eff}}{3 \, \gamma \, \langle v_a \rangle \, \rho_\mathrm{cell} \, f_d} \, .
\end{equation}$

This time scale has two limitations, which we correct as follows.

First, the condensation rate goes to $0$ as $f_d \to 0$, so $\tau_\mathrm{cond}$ diverges. This results from the oversimplified assumption that molecular hydrogen forms exclusively on dust grains. To address it, we follow [Glover2007](https://doi.org/10.1086/519445), who modify the metallicity dependence as $Z \rightarrow Z + Z_\mathrm{eff}$, with $Z_\mathrm{eff} = 10^{-3} \, Z_\odot$. As our rate depends on the dust instead of the metallicity, we use the equivalent floor $f_d \rightarrow f_d + f_d^\mathrm{eff}$, with

$\begin{equation}
    f_d^\mathrm{eff} = 10^{-3} \, \mathcal{D}_\mathrm{ref} \, ,
\end{equation}$

i.e. $10^{-3}$ of the Milky Way dust abundance. This correction both eliminates the divergence and accounts for the formation of molecular hydrogen in environments with very low dust content.

Second, the rate is proportional to $n_a \, n_d$, i.e. to the square of the gas density, so computing it with the mean density of the cell underestimates it when the gas is clumpy on scales below the resolution. Following [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) and [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x), we account for these unresolved density fluctuations with a clumping factor, $C_\rho = \langle \rho^2 \rangle / \langle \rho \rangle^2$, which multiplies the rate (for a detailed justification see Section 3.3 of [Micic2012](https://doi.org/10.1111/j.1365-2966.2012.20477.x)). The clumping factor is also used for the recombination of ionized gas during the epoch of reionization, but we only apply it to the condensation rate. 

The values used in the literature range between $1$ and $100$:

| Source | Context | Clumping factor | Simulation |
| :----: | :----: | :----: | :----: |
| [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) | Section 2.2 | $10$ (fiducial) and tested values of $1$, $10$, and $100$ | Section 2.1 |
| [Gnedin2011](https://doi.org/10.1088/0004-637X/728/2/88) | Section A.7 | $10 - 30$ | Section 2.1 |
| [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x) | Section 2.2.1 | $10$ (taken from [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55)) | Section 2.4 |
| [Micic2012](https://doi.org/10.1111/j.1365-2966.2012.20477.x) | Section 3.3 | $\sim 1-100$ (See $C_\mathrm{\rho, \, eff}$ in Fig. 7 and 8) | -- |
| [Davé2016](https://doi.org/10.1093/mnras/stw1862) | Section 2.3 | $30$ | Section 2.10 |
| [Trayford2025](https://doi.org/10.1093/mnras/staf2040) | Section 3.8 | $1$ to $100$ (density-dependent) | Section 2.5 |

We adopt $C_\rho = 100$, as we find that $C_\rho \gtrsim 100$ is required to form a disk in our cosmological simulations that uses the level 4 initial conditions of the Auriga project ([Grand2017](https://doi.org/10.1093/mnras/stx071)).

Including both corrections, the condensation time scale is

$\begin{align}
    \tau_\mathrm{cond} &= \frac{4 \, s \, a_\mathrm{eff}}{3 \, C_\rho \, \gamma \, \langle v_a \rangle \, \rho_\mathrm{cell} \, (f_d + f_d^\mathrm{eff})} \\
	&= \frac{1}{C_\mathrm{cond} \, \rho_\mathrm{cell} \, (f_d + f_d^\mathrm{eff})} \, ,
\end{align}$

where

$\begin{equation}
    \begin{split}
        C_\mathrm{cond} &= \frac{3 \, C_\rho \, \gamma \, \langle v_a \rangle}{4 \, s \, a_\mathrm{eff}} \\
        &= 41.759 \, \mathrm{mp^{-1} \, Myr^{-1} \, cm^3} \, .
    \end{split}
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ dc45a37a-b939-499f-b9d2-4140aee92282
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Effective grain radius, ā = ⟨a^3⟩ / ⟨a^2⟩, of the size distribution of models
# C and F of Hirashita & Kuo (2011): n(a) ∝ a^(-3.5) between 0.001 and 0.25 μm, 
# with the moments weighted by the number of grains (their eq. 27)
#####################################################################################
let
	# Power-law index, and minimum and maximum grain radii [μm]
	r    = 3.5
	amin = 0.001
	amax = 0.25

	# ℓ-th moment of a power-law size distribution (eq. 27 of Hirashita & Kuo 2011)
	moment(ℓ) = (1 - r) / (ℓ + 1 - r) * (amax^(ℓ + 1 - r) - amin^(ℓ + 1 - r)) / (amax^(1 - r) - amin^(1 - r))

	println(
		"ā = ⟨a^3⟩ / ⟨a^2⟩ = $(round(moment(3) / moment(2); sigdigits=3) * u"μm")"
	)
end
  ╠═╡ =#

# ╔═╡ 0a6364c3-851e-4f0b-b140-f4f0e633fb36
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Initial Mass Function (IMF)

The initial mass function $\phi(m)$ describes the number of stars formed per unit mass interval, i.e. $\phi(m) \, \mathrm{d}m$ is the number of stars with initial masses in the range $[m, m + \mathrm{d}m]$.

We use the IMF of single stars of [Chabrier2003](https://doi.org/10.1086/374879), as given in eq. A4 of [Millan-Irigoyen2021](https://doi.org/10.1093/mnras/stab1969) (Appendix A)

$\begin{equation}
    \phi(m) =
    \begin{cases}
        \dfrac{0.068}{m} \, \exp\left[ -\dfrac{(\log_{10} m - \log_{10} m_c)^2}{2 \, \sigma^2} \right] & \mathrm{for} \ m_\mathrm{low} < m \le 1 \, \mathrm{M_\odot} \, , \\[2ex]
        0.0193 \, m^{-2.3} & \mathrm{for} \ 1 \, \mathrm{M_\odot} < m \le m_\mathrm{high} \, ,
    \end{cases}
\end{equation}$

with

$\begin{align}
    \sigma &= 0.69 \, , \\
    m_c &= 0.079 \, \mathrm{M_\odot} \, , \\
    m_\mathrm{low} &= 0.15 \, \mathrm{M_\odot} \, , \\
    m_\mathrm{high} &= 100 \, \mathrm{M_\odot} \, .
\end{align}$

With these amplitudes, $\phi(m)$ has a small jump at $1 \, \mathrm{M_\odot}$, where the log-normal part gives $0.0190$ and the power law $0.0193$ (a difference below $2\%$). Notice also that eq. A4 of [Millan-Irigoyen2021](https://doi.org/10.1093/mnras/stab1969) is printed with $\log m + \log 0.079$, which is a typo, as it would put the peak of the log-normal part at $12.7 \, \mathrm{M_\odot}$.

The overall amplitude is irrelevant, because every quantity we derive from $\phi(m)$ is divided by the total mass of the stellar population

$\begin{equation}
    M = \int_{m_\mathrm{low}}^{m_\mathrm{high}} m \, \phi(m) \, \mathrm{d}m \, ,
\end{equation}$

as in [Millan-Irigoyen2021](https://doi.org/10.1093/mnras/stab1969) (eq. A5), where the IMF is normalized to a total mass of $1 \, \mathrm{M_\odot}$.
"""
  ╠═╡ =#

# ╔═╡ f9bcdb15-ce25-49a4-9837-53d042065151
#####################################################################################
# Initial mass function
#
# Chabrier (2003) IMF for single stars, as used by HR-pyPopStar
# Millán-Irigoyen et al. (2021), Appendix A, eqs. A4 and A5
# https://doi.org/10.1093/mnras/stab1969
#
# Eq. A4 is printed with log(m) + log(0.079), which is a typo.
# Chabrier (2003) has a minus sign
#
# The overall amplitude is irrelevant, because every use of ϕ is divided by the
# total mass IMF_NORM
#####################################################################################
begin
	# Chosen initial mass function
	const IMF = "Chabrier2003"

	# Key of the chosen IMF in the HR-pyPopStar file names
	const IMF_KEY = "CHA"

	# Lower mass limit for the IMF
	const IMF_M_LOW = 0.15u"Msun"

	# Upper mass limit for the IMF
	const IMF_M_HIGH = 100.0u"Msun"

	# Transition mass between the log-normal and the power-law parts [M⊙]
	const IMF_M_BREAK = 1.0

	# Number of stars per unit mass interval [M⊙^(-1)]
	function ϕ(m::Quantity)
		
		x = ustrip(u"Msun", m)
		
		if x > IMF_M_BREAK
			ϕ_x = 0.0193 * x^(-2.3)
		else
			ϕ_x = 0.068 * exp(-(log10(x) - log10(0.079))^2 / (2.0 * 0.69^2)) / x
		end
		
		return ϕ_x * u"Msun^-1"
		
	end

	# IMF normalization (total mass of the stellar population) [M⊙]
	# The integral is split at the break of ϕ, where it is discontinuous
	const IMF_NORM = quadgk(
		m -> m * ϕ(m), 
		IMF_M_LOW, 
		IMF_M_BREAK * u"Msun", 
		IMF_M_HIGH,
	)[1]
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

Then, the optical depth for ionizing photons is given by

$\begin{equation}
    \tau_\mathrm{ion} = \left(\frac{f_a \, \rho_\mathrm{cell}}{m_p}\right) \sigma_\mathrm{ion} \, h \, ,
\end{equation}$

where $\sigma_\mathrm{ion}$ is the hydrogen photoionization cross-section at the Lyman limit, given by (Table 2 of [Hirashita2002](https://doi.org/10.1046/j.1365-8711.2002.05968.x))

$\begin{equation}
    \sigma_\mathrm{ion} = 6.3 \times 10^{-18} \, \mathrm{cm^2} \, .
\end{equation}$

So, $\tau_\mathrm{ion} = C_\mathrm{\tau ion} \, f_a \, \rho_\mathrm{cell} \, h$, with

$\begin{equation}
    C_\mathrm{\tau ion} = \frac{\sigma_\mathrm{ion}}{m_p} = 6.3 \times 10^{-18} \, \mathrm{mp^{-1} \, cm^2} \, .
\end{equation}$

The corresponding optical depth for dissociation is

$\begin{equation}
    \tau_\mathrm{diss} = \left(\frac{f_m \, \rho_\mathrm{cell}}{2 \, m_p}\right) \sigma_\mathrm{diss,eff} \, h \, ,
\end{equation}$

where the factor of 2 accounts for the two protons in each $\mathrm{H}_2$ molecule, and we use the effective cross-section $\sigma_\mathrm{diss,eff}$ (eq. 22 in [Nickerson2018](https://doi.org/10.1093/mnras/sty1556))

$\begin{equation}
    \sigma_\mathrm{diss,eff} = 2.1 \times 10^{-19} \, \mathrm{cm^2} \, .
\end{equation}$

So, $\tau_\mathrm{diss} = C_\mathrm{\tau diss} \, f_m \, \rho_\mathrm{cell} \, h$, with

$\begin{equation}
    C_\mathrm{\tau diss} = \frac{\sigma_\mathrm{diss,eff}}{2 \, m_p} = 1.05 \times 10^{-19} \, \mathrm{mp^{-1} \, cm^2} \, .
\end{equation}$

Photodissociation of molecular hydrogen by Lyman-Werner photons occurs via the two-step Solomon process ([Draine1996](http://doi.org/10.1086/177689), [Incatasciato2023](https://doi.org/10.1093/mnras/stad1008))

1.  Absorption: An $\mathrm{H}_2$ molecule absorbs a UV photon in the Lyman-Werner bands (absorption lines in the energy range of $11.2 \, \mathrm{eV} < E < 13.6 \, \mathrm{eV}$). This excites the molecule to a higher electronic state.

$\begin{equation}
	\mathrm{H}_2 + \gamma_{LW} \rightarrow \mathrm{H}_2^* \, .
\end{equation}$

2.  Decay: The excited molecule ($\mathrm{H}_2^*$) quickly decays. About $85\%$ of the time, it decays back to a stable vibrational state of the ground electronic level, emitting a photon. However, about $15\%$ of the time, it decays to the vibrational continuum of the ground state, meaning the molecule dissociates into two hydrogen atoms.

$\begin{align}
	&\mathrm{H}_2^* \rightarrow \mathrm{H} + \mathrm{H} \,\, \mathrm{(dissociation, \,\, \sim 15\% \,\, probability)} \\
	&\mathrm{H}_2^* \rightarrow \mathrm{H}_2 + \gamma \,\, \mathrm{(radiative \,\, decay, \,\, \sim 85\% \,\, probability)}
\end{align}$

Using that probability, the unit conversion coefficients are

$\begin{align}
    c_\mathrm{diss} &= 0.3 \, \mathrm{mp} \, , \\
    c_\mathrm{ion} &= 1.0 \, \mathrm{mp} \, .
\end{align}$

These reflect the mass affected per photodissociation or ionization event.
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
end;

# ╔═╡ d59793b2-3414-4051-ac2b-33b6cb8063e8
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The number of photons can be computed from $Q(t', Z)$, defined as the number of photons produced per unit of time by a stellar population of one solar mass, age $t'$, and metallicity $Z$. The current (time $t$) rate of photon production is then

$\begin{equation}
    \dot{N}(t) = \int_0^t \mathrm{SFR}(t - t') \, Q(t', Z(t - t')) \, \mathrm{d}t' \, ,
\end{equation}$

where $\mathrm{SFR}(t - t')$ is the instantaneous SFR (in units of solar mass per unit of time) at the moment of birth of the stellar population of current age $t'$, and $Z(t - t')$ is the metallicity at that moment. The figure below illustrates this integral.
"""
  ╠═╡ =#

# ╔═╡ 225887e5-e214-48a0-8a7c-7814b4e66c83
# ╠═╡ skip_as_script = true
#=╠═╡
TikzPictures.TikzPicture(
	L"""
		\def\xzero{0}
		\def\xyoung{4.0}
		\def\xend{11}
		\def\ybase{0}
		\def\ytop{3.0}

		% Axes
		\draw[->] (\xzero, \ybase) -- (\xend + 0.8, \ybase) node[right] {$t'$};
		\draw[->] (\xzero, \ybase) -- (\xzero, \ytop + 0.3) node[above] {$Q(t', Z)$};

		\draw (\xzero, \ybase) node[below left] {$0$};
		\draw (\xend, \ybase) node[below] {$t$};
		\draw[dashed] (\xend, \ybase) -- (\xend, \ytop + 0.3);

		% Photon production rate per unit stellar mass
		\draw[curve, c_star, domain=0:11, samples=120, variable=\x]
		    plot ({\x}, {\ybase + 2.7 * exp(-1.6 * \x) + 0.10 * exp(-0.05 * \x)});

		% Stellar populations of age t', with mass dM = SFR(t - t') dt'
		\foreach \x in {0.3, 0.9, 1.5, 2.1, 2.7, 3.3, 4, 5, 6, 7, 8, 9, 9.8, 10.6}{
		    \pgfmathsetmacro{\h}{2.7 * exp(-1.6 * \x) + 0.10 * exp(-0.05 * \x)}
		    \draw[dM] (\x - 0.2, \ybase) rectangle (\x + 0.2, \h);
		}

		\draw[->, thin] (4.1, 2.35) to[out=180, in=60] (0.95, 0.80);
		\node[align=left, anchor=west] at (4.2, 2.35) {
		    A population born at\\
		    $t - t'$, with mass\\
		    $\mathrm{d}M = \mathrm{SFR}(t - t') \, \mathrm{d}t'$,\\
		    has a current age $t'$\\
		    and emits photons at\\
		    a rate $Q(t', Z) \, \mathrm{d}M$
		};

		% Young populations, which dominate the integral
		\draw[decorate, decoration={calligraphic brace, amplitude=6pt, raise=2pt}, thick, c_ion, pen colour=c_ion]
		    (\xyoung - 0.5, \ybase - 0.15) -- (\xzero, \ybase - 0.15);
		\node[below=10pt, align=center, c_ion] at ({(\xzero + \xyoung - 0.5) / 2}, \ybase - 0.15) {
		    $t' \lesssim 10$--$100 \, \mathrm{Myr}$\\[1pt]
		    \footnotesize (young stars dominate $Q$)
		};
	""",
	options=raw"x=0.9cm, y=1.5cm, color=white, font=\small, >=Latex",
	height="40em",
	preamble=TIKZ_COLORS * raw"""
		\usepackage{amssymb}
		\usetikzlibrary{arrows.meta, decorations.pathreplacing, calligraphy}
		\tikzset{
		    curve/.style={very thick},
		    dM/.style={fill=c_ion, fill opacity=0.35, draw=c_ion},
		}
	""",
)
  ╠═╡ =#

# ╔═╡ a33249d1-4b85-494c-82c9-0b5a164e32d1
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Each bar is a stellar population of current age $t'$ (the axes are not to scale), formed at $t - t'$ with a mass $\mathrm{d}M = \mathrm{SFR}(t - t') \, \mathrm{d}t'$, which emits photons at a rate $Q(t', Z) \, \mathrm{d}M$. As $Q(t', Z)$ (yellow curve) decays quickly with age, most of the contribution to the integral comes from young blue stars, which die in the first $10$ to $100 \, \mathrm{Myr}$. So, we can approximate

$\begin{align}
    \mathrm{SFR}(t - t') &\approx \mathrm{SFR}(t) \, , \\
	Z(t - t') &\approx Z(t) \, ,
\end{align}$

i.e. we assume that the $\mathrm{SFR}$ and $Z$ stayed approximately constant since the young stars were born (time $t - t'$) until today (time $t$).

Defining $\eta$ as $\eta^*$ without the factor $1 - e^{-\tau}$ of the photons absorbed within the cell,

$\begin{equation}
    \eta = \eta^* \, (1 - e^{-\tau})^{-1} \, ,
\end{equation}$

we end up with the general form

$\begin{equation}
    \eta = c \, \frac{\dot{N}}{\mathrm{SFR}} = c \, \int_0^t Q(t', Z) \, \mathrm{d}t' \, .
\end{equation}$

The value of $Q$ can be calculated using

$\begin{equation}
    Q(t', Z) = \int_{\lambda_1}^{\lambda_2} \frac{L_\lambda(t', Z)}{E_\lambda} \mathrm{d}\lambda = \int_{\lambda_1}^{\lambda_2} \frac{\lambda \, L_\lambda(t', Z)}{h \, c} \, \mathrm{d}\lambda \, ,
\end{equation}$

where $L_\lambda(t', Z)$ is the luminosity per unit of wavelength of a stellar population of one solar mass, age $t'$, and metallicity $Z$, and $E_\lambda = h \, c / \lambda$ is the energy of a photon of wavelength $\lambda$. So, integrating between the wavelengths of interest, we get the number of photons produced per unit of time by a stellar population of the given characteristics.

The luminosity depends not only on the age and metallicity of the population, but also on the initial mass function (IMF), so in principle for each IMF, $t'$, and $Z$ we have a function of luminosity versus $\lambda$.

Using the values from $\texttt{HR-pypopstar}$ ([Millan-Irigoyen2021](https://doi.org/10.1093/mnras/stab1969)) we compute a table of $Q$ for the IMF of [Chabrier2003](https://doi.org/10.1086/374879) (see the section Initial Mass Function), four metallicities ($0.004$, $0.008$, $0.02$, and $0.05$) and more than a hundred ages ranging from $0.1 \, \mathrm{Myr}$ to $15 \, \mathrm{Gyr}$. We use the total (stellar + nebular) SEDs. The nebular continuum contributes less than $10^{-5}$ of the photons in both bands, so the total and the stellar SEDs give the same $Q$.
"""
  ╠═╡ =#

# ╔═╡ 82d801d9-49f4-4d74-93c0-388b09dc76f1
begin
	# Solar luminosity
	# Normalization factor for the spectral energy distribution (SED)
	const u_sed = 3.82e33u"erg * s^-1 * Å^-1 * Msun^-1"

	# Unit for Q_diss and Q_ion
	const u_q = u"Msun^-1 * s^-1"

	# Wavelength range for the photodissociation of hydrogen molecules
	# Lyman-Werner band
    const λ_diss = (911.6u"Å", 1107.0u"Å")

	# Wavelength range for the photoionization of hydrogen atoms
    const λ_ion = (0.0u"Å", 911.6u"Å")

	# Computed magnitudes table
	# HR-pyPopStar (2025)
    # https://www.fractal-es.com/PopStar/
	const MAGNITUDE_TABLE = "./data/HR-pyPopStar_2025_magnitudes.dat"

	# Spectral energy distribution (SED)
    # https://www.fractal-es.com/PopStar/
	# If inside Github actions use the low resolution luminosities
	const HR_LUMINOSITY = get(ENV, "GITHUB_ACTIONS", "false") != "true"

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

	#################################################################################
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
	#################################################################################
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
	
	end

	#################################################################################
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
	#################################################################################
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
		Z       = Vector{Float64}(undef, length(target_files))
		log_age = Vector{Float64}(undef, length(target_files))
		Q_diss  = Vector{Quantity}(undef, length(target_files))
		Q_ion   = Vector{Quantity}(undef, length(target_files))
	
	    for (i, (file, m)) in enumerate(zip(target_files, target_matches))
	
			# Metallicity
			Z[i] = parse(Float64, luminosity_table.str_add * m.captures[2])
	
			# Stellar age as log(age [yr])
			log_age[i] = parse(Float64, m.captures[3])
	
			# Number of dissociating and ionizing photons per unit time
			Q_diss[i], Q_ion[i] = integrate_SED(file, luminosity_table)
			
	    end
	
		Q_df = sort!(
			identity.(DataFrame(; Z, log_age, Q_diss, Q_ion)), 
			[:Z, :log_age],
		)
	
		return Q_df
		
	end
end;

# ╔═╡ c869e896-cc87-4b12-bba2-84b0cf0964c6
begin
	Q_df = compute_Q(luminosity_table, IMF_KEY)

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

	# Integration limits, with the ages of the SED tables below `age` as breakpoints
	#
	# Q is only known at the ages of the SED tables, and it is interpolated linearly
	# in log10(age) between them, so the integrand is smooth within each interval,
	# but has a kink at each age of the table. Also, Q is large only for young
	# populations (the first few Myr), and drops by orders of magnitude after.
	#
	# quadgk places its nodes for a smooth integrand over the whole interval, so for
	# large ages a single interval [0, age] puts few or no nodes in the first few
	# Myr. It can then miss the peak, where most of the integral is, and its error
	# estimate does not notice it, so it returns a wrong value.
	#
	# Splitting the integral at the table ages gives quadgk smooth intervals, and
	# forces it to resolve the young ages, independently of `age`
	knots  = filter(t -> t < age, exp10.(Q_ages .- 6.0))   # [Myr]
	limits = [0.0; knots; age] .* u"Myr"

	# Interpolation function for Q_diss at the given Z
	Q_diss(t) = Qdiss_interp(log10(ustrip(u"yr", t)), Z) * u_q

	# Interpolation function for Q_ion at the given Z
	Q_ion(t) = Qion_interp(log10(ustrip(u"yr", t)), Z) * u_q

	η_diss = uconvert(Unitful.NoUnits, c_diss * quadgk(Q_diss, limits...)[1])
	η_ion  = uconvert(Unitful.NoUnits, c_ion * quadgk(Q_ion, limits...)[1])

	return η_diss, η_ion

end;

# ╔═╡ 9a756404-f35b-40b9-bc88-4b7e3a7df0b6
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Shielding

Not all the photons discussed in the previous section reach atoms or molecules, due to attenuation by dust and molecular gas. This shielding reduces the effective photodissociation and photoionization rates.

Following [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) (eqs. 2a and 2b) and [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x) (eqs. 1 and 2), the transmitted fractions of the ionizing and dissociating fluxes are

$\begin{align}
    e_\mathrm{ion} &= S_d \, , \\
    e_\mathrm{diss} &= S_d \, S_{\mathrm{H}_2} \, ,
\end{align}$

where $S_d$ accounts for the attenuation by dust, and $S_{\mathrm{H}_2}$ for the molecular self-shielding. $S_d$ depends on the dust column density, and $S_{\mathrm{H}_2}$ on the molecular column density. The functional forms of both factors follow [Draine1996](https://doi.org/10.1086/177689) and [Glover2007](https://doi.org/10.1086/512238).

From eq. 8 of [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x), the dust shielding factor is

$\begin{equation}
    S_d = \exp\left[ -\sigma_\mathrm{d, eff} \, \frac{Z}{Z_\odot} \, (N_\mathrm{HI} + 2 \, N_{\mathrm{H}_2}) \right] \, ,
\end{equation}$

where $N_i$ is the column density of species $i$, and $\sigma_\mathrm{d, eff} = 4 \times 10^{-21} \, \mathrm{cm^2}$ is an effective dust cross-section per hydrogen nucleus. This value comes from [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) (Section 2.2), who treat $\sigma_\mathrm{d, eff}$ as an adjustable parameter. They fit it by comparing the atomic-to-molecular transition in their simulations with the molecular fractions measured by FUSE in the Milky Way, the LMC, and the SMC (their Figure 6), and their value is slightly different from the one in [Glover2007](https://doi.org/10.1086/512238).

Both works assume that the dust-to-gas ratio is proportional to the gas metallicity ([Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55), Section 2.2, and [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x), Section 2.2.1). So $\sigma_\mathrm{d, eff}$ is the cross-section per hydrogen nucleus of gas with solar metallicity, i.e. with the dust abundance of the Milky Way, and $Z / Z_\odot$ is only a proxy for the dust abundance relative to that of the Milky Way.

Since we follow the dust explicitly, we undo that substitution

$\begin{equation}
    \frac{Z}{Z_\odot} \rightarrow \frac{\mathcal{D}}{\mathcal{D}_\mathrm{ref}} \, ,
\end{equation}$

where $\mathcal{D} = \rho_d / (m_p \, n_\mathrm{H})$ is the dust-to-hydrogen mass ratio, and $\mathcal{D}_\mathrm{ref}$ its value in the Milky Way. With $N_\mathrm{H} = N_\mathrm{HI} + 2 \, N_{\mathrm{H}_2}$, the dust optical depth is

$\begin{equation}
    \tau_d = \sigma_\mathrm{d, eff} \, \frac{\mathcal{D}}{\mathcal{D}_\mathrm{ref}} \, N_\mathrm{H} = \frac{\sigma_\mathrm{d, eff}}{\mathcal{D}_\mathrm{ref} \, m_p} \, \Sigma_d \, ,
\end{equation}$

where $\Sigma_d = \mathcal{D} \, m_p \, N_\mathrm{H}$ is the dust mass column density, and $\sigma_\mathrm{d, eff} / (\mathcal{D}_\mathrm{ref} \, m_p)$ is the effective opacity per unit of dust mass.

To estimate the column density, we use the shielding length $L_\mathrm{sh}$ (see the section Shielding length), assuming that the shielding material has the density and composition of the cell, so

$\begin{equation}
    \Sigma_d = f_d \, \rho_\mathrm{cell} \, L_\mathrm{sh} \, .
\end{equation}$

Substituting, the dust shielding factor becomes

$\begin{equation}
    S_d = \exp\left(-C_\mathrm{sd} \, f_d \, \rho_\mathrm{cell} \, L_\mathrm{sh}\right) \, ,
\end{equation}$

where

$\begin{equation}
    C_\mathrm{sd} = \frac{\sigma_\mathrm{d, eff}}{\mathcal{D}_\mathrm{ref} \, m_p} \, .
\end{equation}$

[Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) do not give the dust abundance of the gas their fit refers to, so we need an independent estimate of $\mathcal{D}_\mathrm{ref}$. We take it from [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x), who set the abundance of silicate and graphite grains to reproduce the mean extinction curve of the Milky Way per hydrogen nucleus (Section 4.1). For each grain species, with key element $\mathrm{X}$ (Si for silicate and C for graphite), eq. 18 of [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) at solar metallicity gives the dust mass density

$\begin{equation}
    \rho_{d, \mathrm{X}} = \frac{m_\mathrm{X}}{f_\mathrm{X}} \, (1 - \xi_\mathrm{X}) \left( \mathrm{\frac{X}{H}} \right)_{\!\odot} n_\mathrm{H} \, ,
\end{equation}$

where $(\mathrm{X/H})_\odot$ is the solar abundance of $\mathrm{X}$ (by number), $m_\mathrm{X}$ its atomic mass, $f_\mathrm{X}$ the mass fraction of $\mathrm{X}$ in the grains, and $\xi_\mathrm{X}$ the fraction of $\mathrm{X}$ that remains in the gas phase, so $1 - \xi_\mathrm{X}$ is the fraction locked in dust. The values are ([Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x), Table 1, with the abundances and masses from [Cox2002](https://doi.org/10.1007/978-1-4612-1186-0))

| Species | $\mathrm{X}$ | $f_\mathrm{X}$ | $m_\mathrm{X} \, [\mathrm{u}]$ | $(\mathrm{X/H})_\odot$ | $\xi_\mathrm{X}$ |
| :----: | :----: | :----: | :----: | :----: | :----: |
| Silicate | Si | $0.166$ | $28.1$ | $3.55 \times 10^{-5}$ | $0.25$ |
| Graphite | C | $1$ | $12$ | $3.63 \times 10^{-4}$ | $0.15$ |

Summing both species, and defining $c_\mathrm{X} = (1 - \xi_\mathrm{X}) / f_\mathrm{X}$ (the dust mass per unit mass of $\mathrm{X}$, so $c_\mathrm{Si} = 4.52$ and $c_\mathrm{C} = 0.85$)

$\begin{equation}
    \mathcal{D}_\mathrm{ref} = \frac{\rho_d}{m_p \, n_\mathrm{H}} = c_\mathrm{Si} \, \frac{m_\mathrm{Si}}{m_p} \left( \mathrm{\frac{Si}{H}} \right)_{\!\odot} + c_\mathrm{C} \, \frac{m_\mathrm{C}}{m_p} \left( \mathrm{\frac{C}{H}} \right)_{\!\odot} = 8.15 \times 10^{-3} \, ,
\end{equation}$

i.e. $1.363 \times 10^{-26} \, \mathrm{g}$ of dust per hydrogen nucleus. Identifying this value with the dust abundance behind $\sigma_\mathrm{d, eff}$ is our assumption, which links the two works. The solar-metallicity gas of [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) has the dust of the Milky Way, and the grains of [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) reproduce its extinction. 

The corresponding opacity is

$\begin{align}
    C_\mathrm{sd} &= \frac{\sigma_\mathrm{d, eff}}{\mathcal{D}_\mathrm{ref} \, m_p} \\
    &= 4.908 \times 10^{-19} \, \mathrm{mp^{-1} \, cm^2} \, .
\end{align}$
"""
  ╠═╡ =#

# ╔═╡ 6eef3827-1af6-431d-a163-67bb595c337d
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The molecular shielding factor is given by eq. 9 of [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x) as

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
    x = \frac{f_m \, \rho_\mathrm{cell} \, L_\mathrm{sh}}{2 \, m_p \, x_\mathrm{nf}} = C_\mathrm{sh2} \, f_m \, \rho_\mathrm{cell} \, L_\mathrm{sh} \, ,
\end{equation}$

where $x_\mathrm{nf} = 5 \times 10^{14} \, \mathrm{cm^{-2}}$ and

$\begin{align}
    C_\mathrm{sh2} &= \frac{1}{2 \, m_p \, x_\mathrm{nf}} \\
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

# ╔═╡ 79fa6ccd-e137-4178-9376-c1222dd0a8f7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
#### Shielding length

The shielding factors $S_d$ and $S_{\mathrm{H}_2}$ describe the attenuation, by the surrounding material, of radiation coming from outside the cell. So, the column densities must be computed with a length that represents the physical structure of the gas.

Following [Safranek-Shrader2017](https://doi.org/10.1093/mnras/stw2647) (Section 2.5), we approximate the column densities as the product of the local number density and a shielding length, $N = n \, L_\mathrm{sh}$, and use the Jeans length (their eq. 27)

$\begin{equation}
    L_\mathrm{sh} = \left( \frac{\pi \, c_s^2}{G \, \rho_\mathrm{cell}} \right)^{\!1/2} \, ,
\end{equation}$

where $c_s$ is the sound speed. The physical motivation is that the gas able to become molecular sits in dense clumps undergoing gravitational instability, which are shielded by a column of gas of size $\sim L_\mathrm{sh}$. [Safranek-Shrader2017](https://doi.org/10.1093/mnras/stw2647) compare several local approximations for $L_\mathrm{sh}$ with full ray tracing, and find that the Jeans length computed with a temperature ceiling of $T_\mathrm{ceil} = 40 \, \mathrm{K}$ (their model L1a) is among the best performing. Without the ceiling (model L1), warm diffuse gas has an unphysically large shielding length, as $L_\mathrm{sh} \propto T^{1/2}$.

In star-forming cells $\texttt{AREPO}$ uses the two-phase model of [Springel2003](https://doi.org/10.1046/j.1365-8711.2003.06206.x). The internal energy of these cells represents the effective pressure of the two-phase medium, and not the temperature of the shielded gas, so we cannot use it. The cold clouds have a temperature $T_c \simeq 1000 \, \mathrm{K}$, but [Springel2003](https://doi.org/10.1046/j.1365-8711.2003.06206.x) note that their results do not depend on this choice as long as $T_c \ll 10^4 \, \mathrm{K}$, i.e. it is not a physical estimate of the temperature of the molecular gas. Applying the ceiling of model L1a to the cold phase, $T = \min(T_c, T_\mathrm{ceil}) = T_\mathrm{ceil}$, so we will use the isothermal sound speed

$\begin{equation}
    c_s = \left( \frac{k_B \, T_\mathrm{ceil}}{\mu \, m_p} \right)^{\!1/2} = 0.52 \, \mathrm{km \, s^{-1}} \, ,
\end{equation}$

with $\mu = 4 / (1 + 3 \, X_\mathrm{H}) = 1.22$, the mean molecular weight of neutral gas that $\texttt{AREPO}$ uses for the cold phase.

[Safranek-Shrader2017](https://doi.org/10.1093/mnras/stw2647) do not specify their convention for the sound speed (isothermal or adiabatic, and the mean molecular weight), so the normalization of $L_\mathrm{sh}$ is uncertain at the $\sim 30\%$ level.

Then

$\begin{equation}
    L_\mathrm{sh} = C_\mathrm{lsh} \, \rho_\mathrm{cell}^{-1/2} \, ,
\end{equation}$

where

$\begin{align}
    C_\mathrm{lsh} &= \left( \frac{\pi \, c_s^2}{G} \right)^{\!1/2} \\
	&= 2.760 \times 10^{20} \, \mathrm{mp^{1/2} \, cm^{-1/2}} \, ,
\end{align}$

i.e. $L_\mathrm{sh} = 89.5 \, \mathrm{pc}$ for $\rho_\mathrm{cell} = 1 \, \mathrm{mp \, cm^{-3}}$.

With a constant sound speed the column densities scale as $N \propto \rho_\mathrm{cell}^{1/2}$, close to the scaling $N \propto n_\mathrm{H}^{0.42}$ that [Safranek-Shrader2017](https://doi.org/10.1093/mnras/stw2647) find with ray tracing (their eq. 29).

The optical depths $\tau_\mathrm{ion}$ and $\tau_\mathrm{diss}$ used for $\eta_\mathrm{ion}^*$ and $\eta_\mathrm{diss}^*$ (see the section Photodissociation and photoionization) answer a different question than the shielding factors. They set the fraction of the photons emitted by the stars formed in the cell that is absorbed before leaving it, $1 - e^{-\tau}$. As the cell is the control volume of our equations, the relevant path length is set by its geometry, so for these optical depths we keep the cell size $h = 2 \, (3 \, V_\mathrm{cell} / 4 \pi)^{1/3}$ as the length scale.
"""
  ╠═╡ =#

# ╔═╡ e7f0842c-954a-40e1-a112-ec03631e5270
begin
	# Ceiling temperature for the Jeans length
	const T_ceil = 40.0u"K"

	# Hydrogen mass fraction (AREPO HYDROGEN_MASSFRAC)
	const X_H = 0.76

	# Mean molecular weight of the neutral gas (AREPO convention for the cold phase)
	const μ_sh = 4.0 / (1.0 + 3.0 * X_H)

	# Isothermal sound speed at T_ceil
	const cs_sh = sqrt(u"k" * T_ceil / (μ_sh * u"mp"))

	# Shielding length constant
	const C_Lsh = sqrt(π * cs_sh^2 / u"G")

	# Our choice of units for C_Lsh
	const clsh_unit = m_u^(1/2) * l_u^(-1/2)

	# Nominal value of C_Lsh in our choice of units
	const c_lsh = ustrip(Unitful.NoUnits, C_Lsh / clsh_unit)

	# Shielding length [cm]
	L_sh(ρc) = c_lsh / sqrt(ρc)
end;

# ╔═╡ aaba63dd-8fc3-405d-ad19-8e8368e70019
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### UVB Photoionization

We account for the ionizing effect of the metagalactic ultraviolet background (UVB) using the ionization rates provided in Table 3 of [Haardt2012](https://doi.org/10.1088/0004-637X/746/2/125).

These rates correspond to the optically thin limit, as discussed in Section 9.1 of [Haardt2012](https://doi.org/10.1088/0004-637X/746/2/125). To incorporate attenuation due to dust, we apply the dust shielding factor $S_d$. Additionally, since the tabulated photoionization rates $\Gamma_i$ are defined per atom of the species being ionized (see their eq. 18), i.e. per neutral hydrogen atom for hydrogen, the effective rate must be scaled by the atomic hydrogen fraction $f_a$.

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

	#################################################################################
	# Compute ΓUVB
	#
	# Arguments
	#
	#     a: Scale factor [dimensionless]
	#
	# Returns
	#
	#     ΓUVB [Myr^(-1)]
	#################################################################################
	function compute_UVB(a::Float64)::Float64
	
		z = (1 / a) - 1
	
	    ΓUVB = linear_interpolation(
			UVB_TABLE[!, :redshift],
			UVB_TABLE[!, :ΓUVB],
			extrapolation_bc=Flat(),
		)
	
		return ΓUVB(z)
	
	end
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
	\left. \frac{\mathrm{d}}{\mathrm{d}t} f_m(t)\right|_{\mathrm{LWB}} = -f_m(t) \int_{LW} \sigma_{pd}(\nu) \, \frac{F_\nu}{h \, \nu} \, \mathrm{d}\nu =: -f_m(t) \, \Gamma_\mathrm{LWB}(z) \, .
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
	\left. \frac{\mathrm{d}}{\mathrm{d}t} f_m(t)\right|_{\mathrm{LWB}} = \left. -\frac{\mathrm{d}}{\mathrm{d}t} f_a(t)\right|_{\mathrm{LWB}} = -f_m(t) \, S_d \, S_{\mathrm{H}_2} \, \Gamma_\mathrm{LWB} \, .
\end{equation}$

For $J_{21}(z)$ we will use the fit (eq. 9) of [Incatasciato2023](https://doi.org/10.1093/mnras/stad1008)

$\begin{equation}
	\log_{10} J_{21}(z) = A + B \, (1 + z) + C \, (1 + z)^2 \, ,
\end{equation}$

where $A = 2.119$, $B = -1.117 \times 10^{-1}$, and $C = -2.782 \times 10^{-3}$.

We restrict the redshift to the validity range of this fit, $z \geq 6$, so below $z = 6$ we keep the background at its value at $z = 6$

$\begin{equation}
	\Gamma_\mathrm{LWB}(z) = \Gamma_\mathrm{LWB}(\max(z, 6)) \, .
\end{equation}$

We tested other choices for $z < 6$, and the results did not change significantly.
"""
  ╠═╡ =#

# ╔═╡ 5ce798b5-4366-4cc1-8d2d-881a827d6dc9
begin
	const abel97 = ustrip(t_u^-1, 1.38e-12u"s^-1")
	const LW_SOURCE = "Incatasciato2023"

	#################################################################################
	# LW background radiation intensity from Incatasciato et al. (2023), 
	# eq. 9 in [10^(-21) erg s^(-1) cm^(-2) Hz^(-1) sr^(-1)]
	#################################################################################
	const lwb_A = 2.119
	const lwb_B = -1.117e-1
	const lwb_C = -2.782e-3

	J21(z) = exp10(lwb_A + lwb_B * (1 + z) + lwb_C * (1 + z)^2)

	# Lowest redshift of the validity range of the J21 fit
	const LWB_Z_MIN = 6.0

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
	
		return abel97 * J21(max(z, LWB_Z_MIN))
	
	end
end;

# ╔═╡ 7df63fdf-ef41-44f2-8a55-e5c2c849029c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Mass recycling

We consider two key parameters for mass recycling, $R$, the fraction of stellar mass that is returned to the ISM, under the instantaneous recycling approximation, and $Z_\mathrm{SN}$, the metallicity of the recycled gas, i.e., the fraction of the returned mass composed of metals (the rest is assumed to be ionized hydrogen).

The instantaneous recycling approximation (IRA) assumes that stars above a certain mass die and return material to the ISM instantaneously, while lower-mass stars live indefinitely. This approach was first formalized by [Schmidt1963](https://doi.org/10.1086/147553).

Although more accurate treatments avoid the IRA by considering stellar lifetimes (e.g., using empirical relations as done in [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635), Sections 2.2.2-2.2.3), this introduces a time-dependence in $R$, requiring the evaluation of stellar evolution integrals at every timestep, greatly increasing computational cost. For simplicity and efficiency, we assume $R$ is constant over the timescales relevant to our ODE integrations.

The recycled mass fraction $R$ can be computed using a stellar yield model, which specifies the mass of each element ejected by stars of different initial masses. Mathematically,

$\begin{equation}
	R = \dfrac{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} (m - m_\mathrm{rem}(m)) \, \phi(m) \, \mathrm{d}m}{\int_{m_\mathrm{low}}^{m_\mathrm{high}} m \, \phi(m) \, \mathrm{d}m} \, ,
\end{equation}$

where $\phi(m)$ is the IMF, $m_\mathrm{rem}(m)$ is the remnant mass for a star of mass $m$, $m_\mathrm{low}$ and $m_\mathrm{high}$ define the full mass range of the IMF, and $m_\mathrm{ir}$ is the mass threshold above which the IRA applies.

The denominator represents the total mass of the stellar population and acts as a normalization factor for the IMF.

Similarly, the metal yield of the recycled gas is given by

$\begin{equation}
	Z_\mathrm{SN} = \dfrac{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} m \, f_Z \, \phi(m) \, \mathrm{d}m}{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} (m - m_\mathrm{rem}(m)) \, \phi(m) \, \mathrm{d}m} \, ,
\end{equation}$

where $f_Z(m)$ is the fraction of the stellar mass returned as metals.

The two previous equations were taken from [_Nucleosynthesis and Chemical Evolution of Galaxies_](https://doi.org/10.1017/CBO9780511812170) by Bernard Pagel (eq. 7.24 and 7.26), and [_Chemical Evolution
of Galaxies_](https://doi.org/10.1007/978-94-010-0967-6) by Francesca Matteucci (eq. 2.74).

We use the mass limits:

*  $m_\mathrm{low} = 0.15 \, M_\odot$ (lower limit of the IMF of $\texttt{HR-pyPopStar}$),
*  $m_\mathrm{ir} = 8 \, M_\odot$ (approximate threshold for core-collapse supernovae, $\texttt{SNII\_MinMass\_Msun}$ in $\texttt{AREPO}$),
*  $m_\mathrm{high} = 100 \, M_\odot$ (upper limit of the IMF, the same in $\texttt{HR-pyPopStar}$).

[Ascasibar2015](https://doi.org/10.1093/mnras/stv098) found $R \approx 0.18$ and $Z_\mathrm{sn} \approx 0.09$, using the yield model of [Woosley1995](https://doi.org/10.2172/115557) and the IMF of [Kroupa2001](https://doi.org/10.1046/j.1365-8711.2001.04022.x), although the exact mass integration limits used were not specified.

To compute $R$ and $Z_\mathrm{SN}$ we use the stellar yields of [Portinari1998](https://doi.org/10.48550/arXiv.astro-ph/9711337), as compiled and discussed in [Mollá2015](https://doi.org/10.1093/mnras/stv1102).
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
	# Chosen stellar yield model
    const YIELD_MODEL = "Portinari1998"

	# Key of the chosen yield model in the files of Mollá et al. (2015)
	const YIELD_MODEL_KEY = "PCB"
	
	# Mass limit for the IR hypothesis
	const M_IR = 8.0u"Msun"              
	
    # Raw stellar yields from Mollá et al. (2015)
	# https://doi.org/10.1093/mnras/stv1102
    const sy_path = "./data/stellar_yields"

	sy_df = compute_sy(sy_path, YIELD_MODEL_KEY)

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

	R_int   = quadgk(integrand_R, M_IR, IMF_M_HIGH)[1]
	Zsn_int = quadgk(integrand_Zsn, M_IR, IMF_M_HIGH)[1]

	R   = uconvert(Unitful.NoUnits, R_int / IMF_NORM)
	Zsn = uconvert(Unitful.NoUnits, Zsn_int / R_int)

	return R, Zsn

end;

# ╔═╡ cf02b69e-f6e5-4312-9a85-d43b61fac435
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Metallicity of the SN ejecta, Z_SN, at the lowest and highest metallicities of
# the yield tables of Portinari et al. (1998)
#####################################################################################
let
	R_Zsn_lowest  = compute_recycled_fractions(first(sy_metals))
	R_Zsn_highest = compute_recycled_fractions(last(sy_metals))
	
	R_lowest  = R_Zsn_lowest[1]
	R_highest = R_Zsn_highest[1]
	
	Zsn_lowest  = R_Zsn_lowest[2]
	Zsn_highest = R_Zsn_highest[2]

	println("R_lowest    = $(round(R_lowest; sigdigits=3))")
	println("R_highest   = $(round(R_highest; sigdigits=3))")
	println()
	println("Zsn_lowest  = $(round(Zsn_lowest; sigdigits=3))")
	println("Zsn_highest = $(round(Zsn_highest; sigdigits=3))")
end;
  ╠═╡ =#

# ╔═╡ 6cab6cb7-a432-40b6-9390-ad0083fe486d
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Dust destruction

The main mechanisms contributing to dust destruction are sputtering, collisions with cosmic rays, shocks from supernova (SN) explosions, and radiative torques in the presence of intense, anisotropic radiation fields. In this work, we will neglect astration, i.e., the incorporation of dust into stars during star formation.

[Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635) summarized the various physical processes responsible for dust destruction using a single effective timescale, $\tau_\mathrm{dd} \sim 0.8 - 1.8 \, \mathrm{Gyr}$.

We will follow [McKinnon2016](https://doi.org/10.1093/mnras/stw253) and compute $\tau_\mathrm{dd}$ from the Taylor-von Neumann-Sedov blast wave model. Assuming that the main mechanism of dust destruction is SNII shocks. $\tau_\mathrm{dd}$ can be written as (eq. 7 in [McKinnon2016](https://doi.org/10.1093/mnras/stw253)),

$\begin{equation}
	\tau_\mathrm{dd} = \frac{M_\mathrm{cell}}{\epsilon \, \gamma \, M_s(100)} \, ,
\end{equation}$
where $M_\mathrm{cell}$ is the mass of the gas cell, $\epsilon$ is the efficiency with which grains are destroyed in SN shocks, $\gamma$ is the local Type II SN rate, and $M_s(100)$ is the mass of gas shocked to at least $100 \, \mathrm{km \, s^{−1}}$. [McKinnon2016](https://doi.org/10.1093/mnras/stw253) take $\epsilon = 0.3$, but they note that there is a range of physically plausible values. We adopt $\epsilon = 0.1$, the empirical value estimated by [McKee1989](https://doi.org/10.1017/S0074180900125434) and used by [Aoyama2017](https://doi.org/10.1093/mnras/stw3061), [Hirashita2018](https://doi.org/10.1093/mnras/sty2838) (Section 2.2), and [Granato2021](https://doi.org/10.1093/mnras/stab362) (Section 2.3.4). These are the same works on which we will base our dust condensation efficiency for SNe, $\delta_\text{d}$ (see the section Dust production by stars), so dust production and destruction by SNe come from a consistent set of parameters.

The swept-up mass, $M_s(100)$, can be obtained from the Taylor-von Neumann-Sedov theory:

The radius of a blast wave with energy $E$ is

$\begin{equation}
	R(t) = \xi \left(\frac{E \, t^2}{\rho}\right)^{\!1/5} \, ,
\end{equation}$
where $\xi \simeq 1.15167$ is a dimensionless similarity constant and $\rho$ is the mean density of the medium.

So the shock velocity is

$\begin{equation}
	v_s = \frac{dR}{dt} = \frac{2}{5} \, \xi \left(\frac{E}{\rho}\right)^{\!1/5} t^{-3/5} \, .
\end{equation}$

Solving for $t$ and replacing in $R(t)$ gives

$\begin{equation}
	R = \xi^{5/3} \left(\frac{2}{5}\right)^{\!2/3} \left(\frac{E}{\rho}\right)^{\!1/3} v_s^{-2/3} \, .
\end{equation}$

The swept-up mass is

$\begin{equation}
	M_s = \frac{4\pi}{3} \, \rho \, R^3 \, ,
\end{equation}$

so replacing in $R^3$ gives

$\begin{equation}
	M_s = \frac{4\pi}{3} \, \xi^{5} \, \frac{4}{25} \, \frac{E}{v_s^2} \, .
\end{equation}$

Using $E = 1.09 \times 10^{51} \, \mathrm{erg}$ (Table 2 in [McKinnon2016](https://doi.org/10.1093/mnras/stw253)) and the fiducial velocity $v_s = 100 \, \mathrm{km \, s^{−1}}$ we get

$\begin{equation}
	M_s = 7443.3 \, \mathrm{M_\odot} \, ,
\end{equation}$
i.e. $6829 \, \mathrm{M_\odot}$ per $10^{51} \, \mathrm{erg}$, which reproduces the value 

$\begin{equation}
	M_s(100) = E / (\sigma \, v_s^2) = 6800 \, \mathrm{M_\odot} \, (E / 10^{51} \, \mathrm{erg}) \, ,
\end{equation}$

of [McKee1989](https://doi.org/10.1017/S0074180900125434) (eq. 6), used in [McKinnon2016](https://doi.org/10.1093/mnras/stw253) (eq. 8), [Hirashita2018](https://doi.org/10.1093/mnras/sty2838), and [Granato2021](https://doi.org/10.1093/mnras/stab362).

The SNII rate is

$\begin{equation}
	\gamma = \mathrm{SFR} \times f_\mathrm{SNII} \, ,
\end{equation}$

where $\mathrm{SFR} = \psi \, M_\mathrm{cell}$ (see the section Star formation), and $f_\mathrm{SNII}$ is the number of SNII progenitors per unit stellar mass formed. It depends only on the IMF, so it is a constant of the model, and we compute it with our IMF (see the section Initial Mass Function), with SNII progenitors between $m_\mathrm{ir} = 8 \, \mathrm{M_\odot}$ and $m_\mathrm{high} = 100 \, \mathrm{M_\odot}$ (the same limits we use for $R$ and $Z_\mathrm{SN}$)

$\begin{equation}
	f_\mathrm{SNII} = \frac{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} \phi(m) \, \mathrm{d}m}{\int_{m_\mathrm{low}}^{m_\mathrm{high}} m \, \phi(m) \, \mathrm{d}m} = 0.0124 \, \mathrm{M_\odot^{-1}} \, .
\end{equation}$

So, $\tau_\mathrm{dd}$ is

$\begin{equation}
	\tau_\mathrm{dd} = \frac{1}{\epsilon \, M_s \, f_\mathrm{SNII} \, \psi} = \frac{1}{C_\mathrm{dd} \, \psi} \, ,
\end{equation}$

where we have defined the constant

$\begin{align}
	C_\mathrm{dd} &= \epsilon \, M_s \, f_\mathrm{SNII} \\ 
	&= 9.2 \, .
\end{align}$
"""
  ╠═╡ =#

# ╔═╡ 38e91d0d-1020-44a8-becc-8542fd600104
begin
	###########################################################
	# Dust destruction time-scales from McKinnon et al. (2016)
	###########################################################

	# Sedov-Taylor dimensionless similarity constant
	# for an adibatic index of γ = 5/3
	const ξ_SNII = 1.15167

	# SNII energy
	const E_SNII = 1.09e51u"erg"

	# Shock velocity
	const vs = 100.0u"km * s^-1"

	# Numerical constant
	const C_SNII = (4.0/3.0) * π * (4.0 / 25.0)

	# Swept-up mass
	const Ms = C_SNII * ξ_SNII^5 * E_SNII * vs^-2

	# Efficiency with which grains are destroyed in SN shocks
	# Value at a shock velocity of 100 km/s in McKee (1989), as used in Aoyama2017, Hirashita2018, and Granato2021
	# https://doi.org/10.1093/mnras/stw3061
	# https://doi.org/10.1093/mnras/sty2838
	# https://doi.org/10.1093/mnras/stab362
    const ϵ_SNII = 0.1

	# Number of SNII progenitors per unit stellar mass formed
	const f_SNII = ustrip(u"Msun^-1", quadgk(ϕ, M_IR, IMF_M_HIGH)[1] / IMF_NORM)

	# Dust destruction constant [dimensionless]
	const c_dd = ustrip(u"Msun", ϵ_SNII * Ms) * f_SNII
end;

# ╔═╡ 9bf08113-14a4-468a-8766-a6a9b1ed0577
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# M_100 per 10^51 erg of explosion energy, to compare with the 6800 M⊙ per 10^51 erg
# of McKee (1989)
#####################################################################################
println("McKee (1989): $(round(u"Msun", Ms * (1.0e51u"erg" / E_SNII); sigdigits=4))")
  ╠═╡ =#

# ╔═╡ b5ede1bb-8361-4380-b943-c9780d4a6253
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Dust creation

### Dust production by stars

After supernova explosions, the ISM is enriched with both metals and dust. To model dust production from SNe, one can either adopt tabulated dust yields for different SN types (e.g. [Millan-Irigoyen2021](https://doi.org/10.1093/mnras/stab1969)), or assume that a fixed fraction of the metal yields condenses into dust via an effective condensation efficiency (e.g. [Aoyama2017](https://doi.org/10.1093/mnras/stw3061), [Hirashita2019](https://doi.org/10.1093/mnras/sty2838), [Granato2021](https://doi.org/10.1093/mnras/stab362), [Choban2026](
https://doi.org/10.48550/arXiv.2603.08504)). The former approach attempts to capture the complex microphysics of dust formation in stellar ejecta, while the latter provides a simplified, phenomenological prescription. Given the large theoretical and observational uncertainties in SN dust yields (arising from processes such as reverse-shock destruction and poorly constrained condensation conditions) we adopt the simpler approach, encapsulating the relevant physics in a single parameter, the dust condensation efficiency $\delta_\text{d}$.

This contribution to the overall dust budget is expected to be most relevant at early times (i.e. high redshift and low metallicity), whereas at later stages (low redshift and high metallicity) dust growth by accretion of gas-phase metals onto pre-existing grains dominates, rendering the results largely insensitive to the assumed SN condensation efficiency ([Inoue2011](https://doi.org/10.5047/eps.2011.02.013), [Valiante2011](https://doi.org/10.1111/j.1365-2966.2011.19168.x)).

The value of the dust condensation efficiency, $\delta_\text{d}$, depends on both the supernova (SN) type and the resulting dust composition. However, since we do not explicitly distinguish between different SN channels or grain species, we adopt a single effective value of $\delta_\text{d} = 0.1$, consistent with previous works (e.g. [Inoue2011](https://doi.org/10.5047/eps.2011.02.013) (Table 1), [Aoyama2017](https://doi.org/10.1093/mnras/stw3061) (Section 4.4.1), and [Hirashita2018](https://doi.org/10.1093/mnras/sty2838) (Section 2.1)). This choice is further motivated by the fact that, although theoretical models predict that a large fraction of metals ($\sim 0.5 – 0.8$) may initially condense into dust in SN ejecta, a significant portion is subsequently destroyed by reverse shocks before being injected into the ISM. Therefore, the effective efficiency used in simulations should be interpreted as the surviving fraction of dust that successfully reaches the ISM, rather than the initial condensation efficiency ([Zhukovska2008](https://doi.org/10.1051/0004-6361:20077789)).
"""
  ╠═╡ =#

# ╔═╡ 027dbe5e-4993-41c2-983d-289101c118fb
# Dust condensation efficiency
const δ_d = 0.1;

# ╔═╡ 4b81f302-9637-4161-b745-8ac39b9e31d3
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Dust accretion

Following [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635), dust grows by the accretion of gas-phase metals onto pre-existing grains. Here we derive the accretion term, based on [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) (Section 2.1), [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x) (Section 2), and [Hirashita1999](https://doi.org/10.48550/arXiv.astro-ph/9903259). We use $\mathrm{X}$ for the key element of each grain species, $Z_\mathrm{X}$ for its mass fraction in the cell, $n_\mathrm{X}$ for its gas-phase number density, and $\mathcal{N}_\mathrm{X}(a)$ for the grain size distribution.

#### Microscopic relations

We consider two grain species, silicate and graphite, each characterized by a key element $\mathrm{X}$ (Si for silicate, C for graphite). A grain of species $\mathrm{X}$ is a sphere of radius $a$ and constant material density $s_\mathrm{X}$, so

$\begin{align}
	m &= \frac{4}{3} \, \pi \, a^3 \, s_\mathrm{X} \, , \\
	\sigma_d &= \pi \, a^2 \, ,
\end{align}$

where $m$ is the mass of the grain and $\sigma_d$ its cross-section.

The atoms of $\mathrm{X}$ (mass $m_\mathrm{X}$, gas-phase number density $n_\mathrm{X}$) follow a Maxwellian distribution at temperature $T$, so the number of collisions per unit time on a grain is

$\begin{equation}
	\dot{N}_\mathrm{col} = \pi \, a^2 \, n_\mathrm{X} \, \bar{v} = \pi \, a^2 \, n_\mathrm{X} \left( \frac{8 \, k \, T}{\pi \, m_\mathrm{X}} \right)^{\!1/2} \, ,
\end{equation}$

where we used the mean thermal speed $\bar{v} = (8 \, k \, T / \pi \, m_\mathrm{X})^{1/2}$ (eq. 3 of [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x)).

Only a fraction $S$ of the collisions (the sticking coefficient) incorporate the atom permanently into the grain, and the growth is limited by the arrival of the key element $\mathrm{X}$. Each atom of $\mathrm{X}$ that sticks is accompanied by the other elements of the grain in stoichiometric proportion, and adds a total mass $m_\mathrm{X} / f_\mathrm{X}$ to the grain, where $f_\mathrm{X}$ is the mass fraction of $\mathrm{X}$ in the grain material (Section 2.2 of [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x), Section 2.1 of [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x)). Then, the mass of a grain grows as (eq. 4 of [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x))

$\begin{equation}
	\dot{m} = \frac{m_\mathrm{X}}{f_\mathrm{X}} \, S \, \dot{N}_\mathrm{col} = \frac{m_\mathrm{X}}{f_\mathrm{X}} \, S \, \pi \, a^2 \, n_\mathrm{X} \left( \frac{8 \, k \, T}{\pi \, m_\mathrm{X}} \right)^{\!1/2} \, .
\end{equation}$

As the grains are spheres of constant density,

$\begin{equation}
	\dot{m} = \dot{V} \, s_\mathrm{X} = 4 \, \pi \, s_\mathrm{X} \, a^2 \, \dot{a} \, ,
\end{equation}$

and equating both expressions we get (eq. 5 of [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x))

$\begin{equation}
	\dot{a} = \frac{S \, m_\mathrm{X} \, n_\mathrm{X}}{f_\mathrm{X} \, s_\mathrm{X}} \left( \frac{k \, T}{2 \, \pi \, m_\mathrm{X}} \right)^{\!1/2} \, ,
\end{equation}$

which does not depend on $a$.

Let $n_\mathrm{X,tot}$ be the number density of $\mathrm{X}$ in gas and dust, and

$\begin{equation}
	\xi_\mathrm{X} := \frac{n_\mathrm{X}}{n_\mathrm{X,tot}}
\end{equation}$

the fraction of $\mathrm{X}$ in the gas phase. We define the time scale $\tau_\mathrm{X}(a)$ as

$\begin{equation}
	\tau_\mathrm{X}(a) := \xi_\mathrm{X} \, \frac{a}{\dot{a}} \, ,
\end{equation}$

i.e. $a$ divided by the growth rate the grain would have if all of $\mathrm{X}$ were in the gas phase, $\dot{a} / \xi_\mathrm{X}$. 

Using the expression for $\dot{a}$ with $n_\mathrm{X} = \xi_\mathrm{X} \, n_\mathrm{X,tot}$

$\begin{equation}
	\frac{\dot{a}}{\xi_\mathrm{X}} = \frac{S \, m_\mathrm{X} \, n_\mathrm{X,tot}}{f_\mathrm{X} \, s_\mathrm{X}} \left( \frac{k \, T}{2 \, \pi \, m_\mathrm{X}} \right)^{\!1/2} \, ,
\end{equation}$

and replacing in the definition of $\tau_\mathrm{X}(a)$ we get

$\begin{equation}
	\tau_\mathrm{X}(a) = \frac{a \, f_\mathrm{X} \, s_\mathrm{X}}{n_\mathrm{X,tot} \, m_\mathrm{X} \, S} \left( \frac{2 \, \pi \, m_\mathrm{X}}{k \, T} \right)^{\!1/2} \, ,
\end{equation}$

which is eq. 19 of [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x) with $a = a_0$.

Let $\rho_\mathrm{cell}$ be the density of the cell and $Z_\mathrm{X}$ the mass fraction of element $\mathrm{X}$ in the cell, counting gas and dust ($\texttt{AREPO}$ tracks it for Si and C). By definition

$\begin{equation}
	n_\mathrm{X,tot} = \frac{Z_\mathrm{X} \, \rho_\mathrm{cell}}{m_\mathrm{X}} \quad \Longrightarrow \quad n_\mathrm{X,tot} \, m_\mathrm{X} = Z_\mathrm{X} \, \rho_\mathrm{cell} \, ,
\end{equation}$

and replacing in the expression for $\tau_\mathrm{X}(a)$

$\begin{equation}
	\tau_\mathrm{X}(a) = \frac{a \, f_\mathrm{X} \, s_\mathrm{X}}{Z_\mathrm{X} \, \rho_\mathrm{cell} \, S} \left( \frac{2 \, \pi \, m_\mathrm{X}}{k \, T} \right)^{\!1/2} \, .
\end{equation}$

#### Macroscopic relations

Let $\mathcal{N}_\mathrm{X}(a) \, \mathrm{d}a$ be the number density of grains of species $\mathrm{X}$ with radii in $[a, a + \mathrm{d}a]$, and

$\begin{equation}
	\langle a^k \rangle := \frac{1}{N_\mathrm{X}} \int_0^\infty a^k \, \mathcal{N}_\mathrm{X}(a) \, \mathrm{d}a
\end{equation}$

the moments weighted by the number of grains, where

$\begin{equation}
	N_\mathrm{X} = \int_0^\infty \mathcal{N}_\mathrm{X}(a) \, \mathrm{d}a
\end{equation}$

is the total number density of grains. We assume that both species have the same size distribution, so the moments $\langle a^k \rangle$ do not depend on $\mathrm{X}$.

The macroscopic dust density of species $\mathrm{X}$ is

$\begin{align}
	\rho_{d,\mathrm{X}} &= \int_0^\infty m(a) \, \mathcal{N}_\mathrm{X}(a) \, \mathrm{d}a \\
	&= \frac{4 \, \pi}{3} \, s_\mathrm{X} \int_0^\infty a^3 \, \mathcal{N}_\mathrm{X}(a) \, \mathrm{d}a \\
	&= \frac{4 \, \pi}{3} \, s_\mathrm{X} \, N_\mathrm{X} \langle a^3 \rangle \, ,
\end{align}$

and its time derivative is

$\begin{align}
	\frac{\mathrm{d}\rho_{d,\mathrm{X}}}{\mathrm{d}t} &= \int_0^\infty \dot{m}(a) \, \mathcal{N}_\mathrm{X}(a) \, \mathrm{d}a \\
	&= \int_0^\infty 4 \, \pi \, a^2 \, s_\mathrm{X} \, \dot{a} \, \mathcal{N}_\mathrm{X}(a) \, \mathrm{d}a \\
	&= 4 \, \pi \, s_\mathrm{X} \, \dot{a} \int_0^\infty a^2 \, \mathcal{N}_\mathrm{X}(a) \, \mathrm{d}a \\
	&= 4 \, \pi \, s_\mathrm{X} \, \dot{a} \, N_\mathrm{X} \langle a^2 \rangle \, ,
\end{align}$

where we used that $\dot{a}$ does not depend on $a$. So

$\begin{align}
	\frac{\dot{\rho}_{d,\mathrm{X}}}{\rho_{d,\mathrm{X}}} &= \dfrac{4 \, \pi \, s_\mathrm{X} \, \dot{a} \, N_\mathrm{X} \langle a^2 \rangle}{\dfrac{4 \, \pi}{3} \, s_\mathrm{X} \, N_\mathrm{X} \langle a^3 \rangle} \\
	&= 3 \, \dot{a} \, \frac{\langle a^2 \rangle}{\langle a^3 \rangle} \, ,
\end{align}$

and defining $\bar{a} := \langle a^3 \rangle / \langle a^2 \rangle$

$\begin{equation}
	\dot{\rho}_{d,\mathrm{X}} = \frac{3 \, \dot{a}}{\bar{a}} \, \rho_{d,\mathrm{X}} \, .
\end{equation}$

As $\tau_\mathrm{X}(a) \propto a$, $a / \tau_\mathrm{X}(a)$ does not depend on $a$, so $\dot{a} = \xi_\mathrm{X} \, a / \tau_\mathrm{X}(a)$ holds for any $a$, in particular for $\bar{a}$. Then

$\begin{equation}
	\dot{\rho}_{d,\mathrm{X}} = \frac{3 \, \xi_\mathrm{X}}{\tau_\mathrm{X}(\bar{a})} \, \rho_{d,\mathrm{X}} = \frac{\xi_\mathrm{X} \, \rho_{d,\mathrm{X}}}{\tau_{\mathrm{dg},\mathrm{X}}} \, ,
\end{equation}$

where $\tau_{\mathrm{dg},\mathrm{X}} = \tau_\mathrm{X}(\bar{a}) / 3$ is the dust growth time scale of species $\mathrm{X}$. Using the expression for $\tau_\mathrm{X}(a)$ evaluated at $a = \bar{a}$

$\begin{equation}
	\tau_{\mathrm{dg},\mathrm{X}} = \frac{\tau_\mathrm{X}(\bar{a})}{3} = \frac{\bar{a} \, f_\mathrm{X} \, s_\mathrm{X}}{3 \, Z_\mathrm{X} \, \rho_\mathrm{cell} \, S} \left( \frac{2 \, \pi \, m_\mathrm{X}}{k \, T} \right)^{\!1/2} \, .
\end{equation}$

We adopt models C and F of [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x) (Table 2), which are equivalent for our purposes. Both use a power-law size distribution, $\mathcal{N}(a) \propto a^{-3.5}$, between $a_\mathrm{min} = 0.001 \, \mathrm{\mu m}$ and $a_\mathrm{max} = 0.25 \, \mathrm{\mu m}$. With the moments of eq. 27 of [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x), the ratio reduces exactly to

$\begin{equation}
	\bar{a} = \frac{\langle a^3 \rangle_0}{\langle a^2 \rangle_0} = \sqrt{a_\mathrm{min} \, a_\mathrm{max}} = 0.0158 \, \mathrm{\mu m} \, ,
\end{equation}$
the same value $a_\mathrm{eff}$ we use for the condensation of molecular hydrogen.

#### From the grain species to a single dust component

In our model we follow a single dust fraction $f_d$, which adds silicate and graphite. The mass of $\mathrm{X}$ in dust is $(1 - \xi_\mathrm{X}) \, Z_\mathrm{X} \, M_\mathrm{cell}$, and the mass of the grains built on $\mathrm{X}$ is that mass divided by $f_\mathrm{X}$. Dividing by the mass of the cell, the dust fraction of species $\mathrm{X}$ is

$\begin{equation}
	f_{d,\mathrm{X}} = (1 - \xi_\mathrm{X}) \, \frac{Z_\mathrm{X}}{f_\mathrm{X}} \, ,
\end{equation}$

and the total dust fraction is

$\begin{equation}
	f_d = (1 - \xi_\mathrm{Si}) \, \frac{Z_\mathrm{Si}}{f_\mathrm{Si}} + (1 - \xi_\mathrm{C}) \, \frac{Z_\mathrm{C}}{f_\mathrm{C}} \, .
\end{equation}$

This only restates the definitions, so it holds at any time. As we follow a single dust fraction, we can only follow one depletion, so we assume $\xi_\mathrm{Si} = \xi_\mathrm{C} = \xi$. Then

$\begin{equation}
	f_d = (1 - \xi) \, f_d^\mathrm{max} \, , \qquad f_d^\mathrm{max} = \frac{Z_\mathrm{Si}}{f_\mathrm{Si}} + \frac{Z_\mathrm{C}}{f_\mathrm{C}} \, ,
\end{equation}$

where $f_d^\mathrm{max}$ is the dust fraction for $\xi = 0$, i.e. with all the Si and C in grains (see the section Dust initial condition). So

$\begin{equation}
	\xi = 1 - \frac{f_d}{f_d^\mathrm{max}} \, ,
\end{equation}$

and accretion stops when $f_d$ reaches $f_d^\mathrm{max}$.

Notice that the fraction of all the metals in the gas phase, $f_Z / (f_Z + f_d) = 1 - f_d / Z$, is the same expression for $\xi$ with $f_d^\mathrm{max} = Z$, i.e. assuming that all the metals can end up in dust, including the oxygen not bound in silicates and elements that do not condense into these grains (e.g. N and Ne).

Summing $\dot{\rho}_{d,\mathrm{X}}$ over both species, with $\xi_\mathrm{X} = \xi$ and $w_\mathrm{X} := \rho_{d,\mathrm{X}} / \rho_d$,

$\begin{equation}
	\dot{\rho}_d = \sum_\mathrm{X} \frac{\xi \, w_\mathrm{X} \, \rho_d}{\tau_{\mathrm{dg},\mathrm{X}}} = \frac{\xi \, \rho_d}{\tau_\mathrm{dg}} \, ,
\end{equation}$

where the dust growth time scale is given by the mean of the rates of each species, weighted by the composition of the dust

$\begin{equation}
	\frac{1}{\tau_\mathrm{dg}} = \sum_\mathrm{X} \frac{w_\mathrm{X}}{\tau_{\mathrm{dg},\mathrm{X}}} \, .
\end{equation}$

Using the expression for $\tau_{\mathrm{dg},\mathrm{X}}$ we finally get

$\begin{equation}
	\frac{1}{\tau_\mathrm{dg}} = \frac{3 \, S \, \rho_\mathrm{cell}}{\bar{a}} \left( \frac{k \, T}{2 \, \pi} \right)^{\!1/2} \sum_\mathrm{X} \frac{w_\mathrm{X} \, Z_\mathrm{X}}{f_\mathrm{X} \, s_\mathrm{X} \, m_\mathrm{X}^{1/2}} \, .
\end{equation}$

The weights $w_\mathrm{X}$ can be written in terms of the known quantities $Z_\mathrm{X}$ and $f_\mathrm{X}$ as

$\begin{align}
	w_\mathrm{X} &:= \frac{\rho_{d,\mathrm{X}}}{\rho_d} = \frac{f_{d,\mathrm{X}}}{f_d} \\
	&= \frac{(1 - \xi_\mathrm{X}) \, \dfrac{Z_\mathrm{X}}{f_\mathrm{X}}}{\sum_\mathrm{Y} (1 - \xi_\mathrm{Y}) \, \dfrac{Z_\mathrm{Y}}{f_\mathrm{Y}}} \\
	&= \frac{\dfrac{Z_\mathrm{X}}{f_\mathrm{X}}}{\dfrac{Z_\mathrm{Si}}{f_\mathrm{Si}} + \dfrac{Z_\mathrm{C}}{f_\mathrm{C}}} \, ,
\end{align}$

where we used $\xi_\mathrm{Si} = \xi_\mathrm{C} = \xi$. So

$\begin{align}
	w_\mathrm{Si} &= \frac{Z_\mathrm{Si} \, f_\mathrm{C}}{Z_\mathrm{Si} \, f_\mathrm{C} + Z_\mathrm{C} \, f_\mathrm{Si}} \, , \\
	w_\mathrm{C} &= \frac{Z_\mathrm{C} \, f_\mathrm{Si}}{Z_\mathrm{Si} \, f_\mathrm{C} + Z_\mathrm{C} \, f_\mathrm{Si}} \, .
\end{align}$

#### From the grain population to the cell

Dust and depletion are uniform across phases (Section 2 of [Hirashita1999](https://doi.org/10.48550/arXiv.astro-ph/9903259), citing [McKee1989](https://doi.org/10.1017/S0074180900125434)), so the mass of dust exposed to accretion is $M_d \, X_\mathrm{cold}$, where $X_\mathrm{cold}$ is the mass fraction of cold gas. For this reason, [Hirashita1999](https://doi.org/10.48550/arXiv.astro-ph/9903259) (eq. 2) writes the accretion term as

$\begin{equation}
	\left[ \frac{\mathrm{d}M_d}{\mathrm{d}t} \right]_\mathrm{acc} = \frac{M_d \, X_\mathrm{cold} \, \xi}{\tau_\mathrm{dg}} \, .
\end{equation}$

With $X_\mathrm{cold} \simeq f_a + f_m$, and the expression for $\xi$, we get

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_d(t)\right|_{\mathrm{accretion}} = \left( 1 - \frac{f_d}{f_d^\mathrm{max}} \right) (f_a + f_m) \, \frac{f_d}{\tau_\mathrm{dg}} = \frac{f_d}{\tau_\mathrm{da}} \, ,
\end{equation}$

where $\tau_\mathrm{da}$ is the dust accretion time scale, defined by this expression.

#### Grain properties

The properties of each species are taken from Table 1 of [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x)

| Species  | $\mathrm{X}$ | $f_\mathrm{X}$ | $s_\mathrm{X}$                 | $m_\mathrm{X}$         |
|:--------:|:------------:|:--------------:|:------------------------------:|:----------------------:|
| Silicate | Si           | $0.166$        | $3.3 \, \mathrm{g \, cm^{-3}}$  | $28.1 \, \mathrm{u}$   |
| Graphite | C            | $1$            | $2.26 \, \mathrm{g \, cm^{-3}}$ | $12 \, \mathrm{u}$     |

As a consistency check, evaluating $\tau_\mathrm{X}(a)$ with the fiducial values of [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) ($a = 0.1 \, \mathrm{\mu m}$, $T = 10 \, \mathrm{K}$, $S = 0.3$, and $n_\mathrm{H} = 1000 \, \mathrm{cm^{-3}}$), and their solar abundances ($n_\mathrm{X,tot} = (\mathrm{X/H})_\odot \, n_\mathrm{H}$, with $(\mathrm{Si/H})_\odot = 3.55 \times 10^{-5}$ and $(\mathrm{C/H})_\odot = 3.63 \times 10^{-4}$ from their Table 1), gives

$\begin{align}
	\tau_\mathrm{Si} &= 1.61 \times 10^8 \, \mathrm{yr} \, , \\ 
	\tau_\mathrm{C} &= 0.994 \times 10^8 \, \mathrm{yr} \, ,
\end{align}$ 

the normalization constants $A_\mathrm{X}$ of eqs. 19 and 20 of [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x).
"""
  ╠═╡ =#

# ╔═╡ 866df30a-7a45-4821-a4ef-fa886040d9c0
# ╠═╡ skip_as_script = true
#=╠═╡
let
	###############################################################################
	# Consistency check: τ_X(a) with the fiducial values and solar abundances of
	# Hirashita (2012) reproduces the normalization constants A_X of their eqs. 19-20
	###############################################################################

	# Fiducial values of Hirashita (2012)
	a  = 0.1u"μm"
	T  = 10.0u"K"
	S  = 0.3
	nH = 1000.0u"cm^-3"

	# Silicate
	f_Si  = 0.166
	m_Si  = 28.1u"u"
	s_Si  = 3.3u"g * cm^-3"

	# Graphite
	f_C  = 1.0
	m_C  = 12.0u"u"
	s_C  = 2.26u"g * cm^-3"
	
	# Growth time scale, with n_X,tot = (X/H) * nH
	τ(; fx, mx, sx, XH) = (a * fx * sx) / (XH * nH * mx * S) * sqrt((2 * π * mx) / (1.0u"k" * T))

	τ_Si = ustrip(u"yr", τ(; fx=f_Si, mx=m_Si, sx=s_Si, XH=3.55e-5))
	τ_C  = ustrip(u"yr", τ(; fx=f_C, mx=m_C, sx=s_C, XH=3.63e-4))

	println("Hirashita (2012):\n")
	println("\tA_Si = $(round(τ_Si; sigdigits=3)) yr")
	println("\tA_C = $(round(τ_C; sigdigits=3)) yr")
end;
  ╠═╡ =#

# ╔═╡ 9af6ce74-5678-4b48-8758-37b0d5a6f0e4
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
#### $T$ -- Gas temperature

The gas temperature depends on the phase. For consistency with the assumption that accretion happens in the cold neutral gas, we adopt

$\begin{equation}
    T = 100 \, \mathrm{K} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ a9c6a292-086e-4aa0-9856-78d3ab3fbe35
begin
	# Gas temperature
	const T = 100.0u"K"
end;

# ╔═╡ 5890b699-b7de-47a3-bee7-1e7dd7663fbe
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
#### $S$ -- Sticking coefficient

For the sticking coefficient we adopt the [Grassi2011](https://doi.org/10.1051/0004-6361/200913779) fit of the [Leitch-Devlin1985](https://doi.org/10.1093/mnras/213.2.295) results, which reads

$\begin{equation}
    S(T_g, T_d) = 0.019 \, T_g \, (0.0017 \, T_d + 0.4) \, \exp(-0.007 \, T_g) \, .
\end{equation}$

where $T_g$ and $T_d$ are the gas and dust temperatures, respectively, in Kelvin.

We use the gas temperature $T_g = T = 100 \, \mathrm{K}$, and a dust temperature $T_d = 10 \, \mathrm{K}$ (as assumed by [Safranek-Shrader2017](https://doi.org/10.1093/mnras/stw2647), Section 2.4), since grains in neutral gas are much colder than the gas. We obtain

$\begin{equation}
    S = 0.393 \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ ef79392b-395c-497a-9c0e-dc2cd468f6e1
begin
	# Fit for the sticking coefficient from Grassi et al. (2011)
    Se(Tg, Td) = 0.019 * Tg * (0.0017 * Td + 0.4) * exp(-0.007 * Tg)

	# Dust temperature
	# Safranek-Shrader et al. (2017), Section 2.4
	# https://doi.org/10.1093/mnras/stw2647
	const Td = 10.0u"K"

	# Sticking coefficient for a gas at T and dust at Td
	const S = Se(ustrip(u"K", T), ustrip(u"K", Td))
end;

# ╔═╡ 8c9ab125-2acb-4732-a9bf-7838e819e4f7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
#### Final form of $\tau_\mathrm{dg}$ and $\tau_\mathrm{da}$

From the expression for $1 / \tau_\mathrm{dg}$ derived above, the dust growth time scale can be written as

$\begin{equation}
	\tau_\mathrm{dg} = \frac{1}{C_\mathrm{dg} \, \rho_\mathrm{cell}} \, ,
\end{equation}$

where

$\begin{equation}
	C_\mathrm{dg} = \sum_\mathrm{X} w_\mathrm{X} \, C_\mathrm{acc,X} \, Z_\mathrm{X} \, , \qquad C_\mathrm{acc,X} = \frac{3 \, S}{\bar{a} \, f_\mathrm{X} \, s_\mathrm{X}} \left( \frac{k \, T}{2 \, \pi \, m_\mathrm{X}} \right)^{\!1/2} \, ,
\end{equation}$

with 

$\begin{align}
	C_\mathrm{acc,Si} &= 0.4940 \, \mathrm{mp^{-1} \, Myr^{-1} \, cm^3} \, , \\
	C_\mathrm{acc,C} &= 0.1832 \, \mathrm{mp^{-1} \, Myr^{-1} \, cm^3} \, .
\end{align}$

Unlike the other constants of the model, $C_\mathrm{dg}$ depends on the cell, through the abundances $Z_\mathrm{Si}$ and $Z_\mathrm{C}$ (which also set the weights $w_\mathrm{X}$), so it is computed at the start of each integration.

The accretion term can be written as $f_d / \tau_\mathrm{da}$, which defines the dust accretion time scale

$\begin{equation}
	\tau_\mathrm{da} = \frac{\tau_\mathrm{dg}}{\xi \, f_n} = \frac{1}{C_\mathrm{dg} \, (1 - f_d / f_d^\mathrm{max}) \, f_n \, \rho_\mathrm{cell}} \, ,
\end{equation}$

where $f_n = f_a + f_m$.

#### Dust accretion term

The dust growth source term finally becomes

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_d(t)\right|_{\text{accretion}} = C_\mathrm{dg} \left( 1 - \frac{f_d}{f_d^\mathrm{max}} \right) f_d \, f_n \, \rho_\mathrm{cell} \, .
\end{equation}$

Cells without Si and C have no dust and no accretion, so we set $C_\mathrm{dg} = 0$, and $f_d^\mathrm{max} = 1$ in the ODEs to avoid dividing by zero.

As $f_d^\mathrm{max}$ is kept constant during the integration, dust produced by SNe can exceed it in cells with very low metallicity. Accretion can only add metals to the grains, so in the ODEs we use $\xi = \max(0, 1 - f_d / f_d^\mathrm{max})$. Which means that above $f_d^\mathrm{max}$ there is nothing left to accrete, and the term is $0$ instead of turning dust back into metals. For $0 \leq f_d \leq f_d^\mathrm{max}$ this is the expression above.

The metallicity and the Si and C abundances computed by $\texttt{AREPO}$ can have unphysical values (e.g. slightly negative, or inconsistent with each other) when they are close to $0$. So, before computing the initial conditions and the accretion parameters, they are clamped at $0$, cells with $Z < 10^{-10}$ (about $10^{-8} \, Z_\odot$) are considered pristine ($Z = Z_\mathrm{Si} = Z_\mathrm{C} = 0$), and $Z_\mathrm{Si}$ and $Z_\mathrm{C}$ are scaled down if their sum exceeds $Z$.
"""
  ╠═╡ =#

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
	const N_PAR = 10 

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

# ╔═╡ 1ec99792-905d-4e1b-a413-ef58143d3c68
# ╠═╡ skip_as_script = true
#=╠═╡
md"## ODE function"
  ╠═╡ =#

# ╔═╡ 2afd459c-90e9-4105-9121-27e21bb89eeb
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Integration"
  ╠═╡ =#

# ╔═╡ ab1b7695-3cd6-4e67-9cfd-fbaf2f4d1d15
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Initial conditions

Schematic diagram for the contribution of each $\texttt{AREPO}$ component to the initial condition:
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

# ╔═╡ 14a6478e-939f-4e9c-b905-4cec39917d47
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Gas initial conditions

Each time a cell enters the model, the gas phases and the stars are set from the current state of the cell in $\texttt{AREPO}$. Star-forming cells follow the two-phase model of [Springel2003](https://doi.org/10.1046/j.1365-8711.2003.06206.x): a mass fraction $x$ of cold clouds embedded in a hot phase. We take $x$ from the effective equation of state of $\texttt{AREPO}$, and the abundance of ionized hydrogen in the hot phase, $y_\mathrm{HII} = n_\mathrm{HII} / n_\mathrm{H}$, from its cooling module. Assuming that the cold clouds are fully neutral, the ionized and neutral fractions of the gas are

$\begin{align}
    f_i^\mathrm{gas} &= (1 - x) \, y_\mathrm{HII} \, , \\
    f_n^\mathrm{gas} &= x + (1 - x) \, (1 - y_\mathrm{HII}) \, ,
\end{align}$

with $x$ and $y_\mathrm{HII}$ clamped to $[0, 1]$. As we do not follow helium (see the section Phases and notation), all the gas that is not metals is distributed among the gas phases, so with the metallicity $Z$ of the cell (gas + dust) the initial conditions are

$\begin{align}
    f_i(0) &= (1 - Z) \, f_i^\mathrm{gas} \, , \\
    f_a(0) &= (1 - Z) \, f_n^\mathrm{gas} \, , \\
    f_m(0) &= 0 \, , \\
    f_s(0) &= 0 \, , \\
    f_Z(0) &= Z - f_d(0) \, ,
\end{align}$

where $f_d(0)$ is given in the next section. All the neutral gas starts as atomic, and the molecular gas forms during the integration. As $f_i^\mathrm{gas} + f_n^\mathrm{gas} = 1$, the fractions add up to one.

### Dust initial condition

Following eqs. 17 and 18 of [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x), the initial mass density of a grain species with key element $\mathrm{X}$ (Si for silicate, C for graphite) is

$\begin{equation}
	\rho_{d,\mathrm{X}}(0) = \frac{m_\mathrm{X}}{f_\mathrm{X}} \, (1 - \xi_\mathrm{X}) \, n_\mathrm{X}^\mathrm{tot} \, ,
\end{equation}$

where $f_\mathrm{X}$ is the mass fraction of element $\mathrm{X}$ in the grain material, $m_\mathrm{X}$ is the atomic mass of $\mathrm{X}$, $n_\mathrm{X}^\mathrm{tot}$ is the number density of $\mathrm{X}$ in both gas and dust, and $\xi_\mathrm{X}$ is the fraction of $\mathrm{X}$ in the gas phase

$\begin{equation}
	\xi_\mathrm{X} = \frac{n_\mathrm{X}^\mathrm{gas}}{n_\mathrm{X}^\mathrm{tot}} \, .
\end{equation}$

[Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) closes this expression with $n_\mathrm{X}^\mathrm{tot} = (Z / Z_\odot) \, (\mathrm{X/H})_\odot \, n_\mathrm{H}$ (eq. 17), which assumes a solar abundance pattern. Instead, $\texttt{AREPO}$ tracks the mass of nine elements in each cell (H, He, C, N, O, Ne, Mg, Si, and Fe; [McKinnon2017](https://doi.org/10.1093/mnras/stx467)), so we can use the cell abundances directly

$\begin{equation}
	m_\mathrm{X} \, n_\mathrm{X}^\mathrm{tot} = Z_\mathrm{X} \, \rho_\mathrm{cell} \, ,
\end{equation}$

where $Z_\mathrm{X} = M_\mathrm{X} / M_\mathrm{cell}$ is the mass fraction of element $\mathrm{X}$ (gas + dust). Silicate and graphite are separate mass reservoirs, so their contributions are summed, and the initial dust fraction is

$\begin{equation}
	f_d(0) = \sum_\mathrm{X} \frac{\rho_{d,\mathrm{X}}(0)}{\rho_\mathrm{cell}} = c_\mathrm{Si} \, Z_\mathrm{Si} + c_\mathrm{C} \, Z_\mathrm{C} \, ,
\end{equation}$

with

$\begin{equation}
	c_\mathrm{X} = \frac{1 - \xi_\mathrm{X}}{f_\mathrm{X}} \, .
\end{equation}$

We adopt the following values

|  Species | $\mathrm{X}$ | $f_\mathrm{X}$ | $\xi_\mathrm{X}$ | $c_\mathrm{X}$ |
|:--------:|:------------:|:--------------:|:----------------:|:--------------:|
| Silicate |      Si      |     $0.166$    |      $0.25$      |     $4.52$     |
| Graphite |       C      |       $1$      |      $0.15$      |     $0.85$     |

where $f_\mathrm{Si} = 0.166$ corresponds to the composition $\mathrm{Mg_{1.1} Fe_{0.9} Si O_4}$ (Table 1 of [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x)), and $\xi_\mathrm{Si}$ and $\xi_\mathrm{C}$ are the values that [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) (Section 4.1) uses so the initial grain abundance reproduces the mean extinction curve of the Milky Way. They are close to that paper's fiducial $\xi(0) = 0.3$, which roughly matches the depletion in the diffuse medium (see Section 2.3 of [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x), where $\xi(0) = 0.1 - 0.7$ is considered a reasonable range). Thus, $f_d(0)$ represents the dust content of the diffuse ISM before any growth in the dense gas, which is then followed by the ODEs.

#### Maximum dust fraction

Dust can grow until all the key elements are in grains ($\xi_\mathrm{X} = 0$), so the maximum dust fraction is

$\begin{equation}
	f_d^\mathrm{max} = \min\left( \frac{Z_\mathrm{Si}}{f_\mathrm{Si}} + \frac{Z_\mathrm{C}}{f_\mathrm{C}}, Z \right) \, ,
\end{equation}$

where the cap at $Z$ only matters for unrealistic abundance patterns, and guarantees $f_Z \geq 0$. For the present-day solar photosphere of [Asplund2009](https://doi.org/10.1146/annurev.astro.46.060407.145222), $f_d^\mathrm{max} / Z = 0.476$, which means that dust can lock at most about half of the metals, because most of the oxygen and elements such as N and Ne do not condense into silicate or graphite grains. For cells without Si and C, $f_d^\mathrm{max} = 0$, so there is no dust and no accretion.

#### Cells that enter the model more than once

After an integration, the stochastic star formation of $\texttt{AREPO}$ may keep a cell as gas, and the cell can enter the model again. The gas phases, the stars, and the metallicity are always reset from the current state of the cell in $\texttt{AREPO}$, so the model remains a closed box that, given the state of a cell, returns its SFR and the fraction of each component. The exception is the dust, because its time scales (accretion and destruction, see the section Dust creation) are much longer than the integration time, so resetting it at every entry would erase its evolution. At the end of each integration, we store the final fractions of dust and metals of the cell, $f_d$ and $f_Z$, and the next time the cell enters the model we compute from them the dust-to-metal ratio of the previous integration,

$\begin{equation}
	\mathrm{DZ} = \frac{f_d}{f_Z + f_d} \, ,
\end{equation}$

and combine it with the current metallicity of the cell

$\begin{equation}
	f_d(0) = \min(\mathrm{DZ} \, Z, f_d^\mathrm{max}) \, .
\end{equation}$

The dust-to-metal ratio is an intensive quantity, so it can be combined with the metallicity of the cell after $\texttt{AREPO}$ has changed its metal content (stellar enrichment, or mass exchange with the neighboring cells). 

The initial condition of [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) given above is used only the first time a cell enters the model, and when the stored $f_Z + f_d$ is not positive, as the dust-to-metal ratio is then undefined. We do not reset the dust-to-metal ratio of cells that leave the star-forming regime and enter it again later, as the dust processes outside star-forming cells (e.g. thermal sputtering in hot gas) are not included in the model.
"""
  ╠═╡ =#

# ╔═╡ d7cc8a66-220a-4e9b-b42e-a0e085ed3a0f
begin
	# Fraction of the key element in the gas phase 
	# Hirashita 2012, Sect. 4.1
	const ξ_Si = 0.25   # Silicate (key element Si)
	const ξ_C  = 0.15   # Graphite (key element C)

	# Mass fraction of the key element in the grain material 
	# Hirashita 2012, Table 1
	const f_Si = 0.166  # Silicate (Mg1.1 Fe0.9 Si O4)
	const f_C  = 1.0    # Graphite

	# Initial dust fraction coefficients
	const c_Si = (1.0 - ξ_Si) / f_Si
	const c_C  = (1.0 - ξ_C) / f_C

	#################################################################################
	# Initial dust fraction
	#
	# Arguments
	#
	#     Z_Si: Mass fraction of Si in the cell (gas + dust)
	#     Z_C:  Mass fraction of C in the cell (gas + dust)
	#     Z:    Total metallicity of the cell (gas + dust)
	#
	# Returns
	#
	#     fd(0), capped at Z so that fZ(0) = Z - fd(0) >= 0
	#################################################################################
	initial_dust_fraction(Z_Si, Z_C, Z) = min(c_Si * Z_Si + c_C * Z_C, Z)

	# Maximum dust mass per unit mass of the key element (all of it in grains)
	const cmax_Si = 1.0 / f_Si  # Silicate
	const cmax_C  = 1.0 / f_C   # Graphite

	#################################################################################
	# Maximum dust fraction (all Si and C in grains)
	#
	# Arguments
	#
	#     Z_Si: Mass fraction of Si in the cell (gas + dust)
	#     Z_C:  Mass fraction of C in the cell (gas + dust)
	#     Z:    Total metallicity of the cell (gas + dust)
	#
	# Returns
	#
	#     fd_max, capped at Z
	#################################################################################
	max_dust_fraction(Z_Si, Z_C, Z) = min(cmax_Si * Z_Si + cmax_C * Z_C, Z)

	# Present-day solar photosphere 
	# Asplund 2009, Tables 1 and 4
	const X_H_A09  = 0.7381 # Hydrogen mass fraction
	const Z_A09    = 0.0134 # Metallicity
	const Z_C_A09  = X_H_A09 * 10^(8.43 - 12) * 12.011 / 1.008 # Carbon mass fraction
	const Z_Si_A09 = X_H_A09 * 10^(7.51 - 12) * 28.086 / 1.008 # Silicon mass fraction

	# Atomic mass of the key element 
	# Hirashita 2012, Table 1
	const m_Si = 28.1u"u"  # Silicon
	const m_C  = 12.0u"u"  # Carbon

	# Solar abundance of the key element, by number relative to H 
	# Hirashita 2012, Table 1
	const SiH_sun = 3.55e-5  # Silicon
	const CH_sun  = 3.63e-4  # Carbon

	# Reference (Milky Way) dust-to-hydrogen mass ratio 
	# Hirashita 2012, eq. 18 and Sect. 4.1
	const D_ref = ustrip(
		Unitful.NoUnits,
		(c_Si * m_Si * SiH_sun + c_C * m_C * CH_sun) / 1.0u"mp",
	)
end;

# ╔═╡ 568d9fe3-6716-4a9a-ba1d-b9e6fd039150
begin
	# Dust grains internal density
    const s_d = 3.0u"g * cm^-3"

	# Effective grain radius ⟨a³⟩/⟨a²⟩
    const a_eff = 0.0158u"μm"
	
    # Recombination coefficient (γ) times 
	# the average thermal velocity of hydrogen atoms (⟨va⟩)
	const γva = 5.0e4u"cm * s^-1"

	# Dust fraction floor (1e-3 of the Milky Way dust-to-hydrogen mass ratio)
	const fd_eff = 1e-3 * D_ref

	# Clumping factor
	const Cρ = 100.0

	# Condensation constant
	const C_cond = (3 * Cρ * γva) / (4 * s_d * a_eff)

	# Our choice of units for C_cond
	const cc_unit = m_u^-1 * t_u^-1 * l_u^3

	# Nominal value of C_cond in our choice of units
	const c_cond = ustrip(Unitful.NoUnits, C_cond / cc_unit)

	# Condensation time-scale
	τ_cond(ρc, fd) = 1.0 / (c_cond * ρc * (fd + fd_eff))
end;

# ╔═╡ 805ef471-1bd2-4ce2-bf50-eeb9c6b10e4d
begin
	# Effective cross section
	const σdeff = 4.0e-21u"cm^2"

	# Dust shielding constant
	const C_sd = σdeff / (D_ref * u"mp")

	# Unit of C_sd
	const csd_unit = m_u^-1 * l_u^2

	# Nominal value of C_sd in our choice of units
	const c_sd = ustrip(Unitful.NoUnits, C_sd / csd_unit)
end;

# ╔═╡ bd8743d6-8f21-413d-835a-e543926baa09
#####################################################################################
# System of ODEs, where each equation has units of Myr^(-1), and
#
# Ionized gas fraction:   fi(t) = Mi(t) / M_cell --> y[1]
# Atomic gas fraction:    fa(t) = Ma(t) / M_cell --> y[2]
# Molecular gas fraction: fm(t) = Mm(t) / M_cell --> y[3]
# Stellar fraction:       fs(t) = Ms(t) / M_cell --> y[4]
# Metals fraction:        fZ(t) = MZ(t) / M_cell --> y[5]
# Dust fraction:          fd(t) = Md(t) / M_cell --> y[6]
# 
# This function follows the API requirements of DifferentialEquations.jl
#####################################################################################
function system!(dydt, ic, parameters, t)

	##################################################################
    # Initial conditions
	##################################################################

    fi, fa, fm, _, _, fd = ic

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
	# h:      Cell size                    [cm]
	# fd_max: Maximum dust fraction        [dimensionless]
	# c_dg:   Dust growth constant         [Myr^(-1) * mp^(-1) * cm^3]
	##################################################################

    ρ_cell, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h, fd_max, c_dg = parameters

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

	############
	# Shielding
    ############

	# Shielding length [cm]
	Lsh = L_sh(ρ_cell)

	# Dust optical depth 
	τ_d = c_sd * fd * ρ_cell * Lsh
	
	# Dust shielding 
	sd = exp(-τ_d)

	# Molecular self-shielding optical depth
	xp1     = c_sh2 * max(fm, 0.0) * ρ_cell * Lsh + 1.0
	xp1_2   = xp1 * xp1
	sq_xp1  = sqrt(xp1)
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
	recomb = c_rec * fi * max(fi, 0.0) * ρ_cell

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
	cond = c_cond * fa * ρ_cell * (fd + fd_eff)

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

	# Dust yield by SNe [Myr^(-1)]
	dg_SNe = δ_d * Zsn * R_ψs

	# Gas-phase fraction of the key elements
	ξ = max(0.0, 1.0 - fd / fd_max)

	# Dust accretion [Myr^(-1)]
	dg_acc = c_dg * ξ * fd * fn * ρ_cell

	# Dust destruction [Myr^(-1)]
	dd = c_dd * fd * ψs

	net_dust_growth = dg_SNe + dg_acc - dd

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
#     parameters: Parameters for the ODEs, 
#                 [ρ, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h, fd_max, c_dg]
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

# ╔═╡ 4e1140be-3e61-47bd-9f9f-6b5dfbff6ae2
#####################################################################################
# System to be solved by DifferentialEquations.jl
#####################################################################################
ode_function = ODEFunction{true}(
	system!;
	jac=jacobian!,
	tgrad=(dt, ic, p, t) -> nothing,
);

# ╔═╡ 2a39d6f8-da49-4976-9aa7-889391e55a5d
begin
	# Material density of the grains (Hirashita 2012, Table 1)
	const s_Si = 3.3u"g * cm^-3"   # Silicate
	const s_C  = 2.26u"g * cm^-3"  # Graphite

	# Accretion constant of each species, 
	# c_X = 3 * S / (ā * f_X * s_X) * sqrt(k * T / (2π * m_X))
	# ā = ⟨a³⟩/⟨a²⟩ = a_eff (models C and F of Hirashita 2011, Table 2)
	const C_acc_Si = (3 * S) / (a_eff * f_Si * s_Si) * sqrt(u"k" * T / (2 * π * m_Si))
	const C_acc_C  = (3 * S) / (a_eff * f_C * s_C) * sqrt(u"k" * T / (2 * π * m_C))

	# Our choice of units for C_acc
	const cacc_unit = m_u^-1 * t_u^-1 * l_u^3

	# Nominal value of C_acc in our choice of units
	const c_acc_Si = ustrip(Unitful.NoUnits, C_acc_Si / cacc_unit)
	const c_acc_C  = ustrip(Unitful.NoUnits, C_acc_C / cacc_unit)

	#################################################################################
	# Accretion parameters of a cell
	#
	# Arguments
	#
	#     Z_Si: Mass fraction of Si in the cell (gas + dust)
	#     Z_C:  Mass fraction of C in the cell (gas + dust)
	#     Z:    Total metallicity of the cell (gas + dust)
	#
	# Returns
	#
	#     fd_max: Maximum dust fraction (all Si and C in grains), capped at Z.
	#             For cells without Si nor C (no accretion) it is set to 1
	#             to avoid dividing by 0 in the ODEs
	#     c_dg:   Dust growth constant [Myr^(-1) * mp^(-1) * cm^3],
	#             It is 0 for cells without Si nor C
	#################################################################################
	function accretion_parameters(Z_Si, Z_C, Z)

		# Maximum silicate and graphite fractions
		y_Si = cmax_Si * Z_Si
		y_C  = cmax_C * Z_C

		# Maximum dust fraction
		fd_max = min(y_Si + y_C, Z)

		if fd_max <= 0.0
			return 1.0, 0.0
		end

		# Dust growth constant
		c_dg = (c_acc_Si * Z_Si * y_Si + c_acc_C * Z_C * y_C) / (y_Si + y_C)

		return fd_max, c_dg

	end
end;

# ╔═╡ 4173e91c-c01a-4fa2-ae77-0408bf7d9a1b
#####################################################################################
# Solve the system of ODEs
#
# Arguments
# 
#     ic:          Initial conditions, [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
#     base_params: Parameters for the ODEs, [ρ_cell, Z, a, h, Z_Si, Z_C]
#     tspan:       Integration span, (ti, tf) [Myr]
#     times:       Times at which the solution will be returned [Myr]
#     args:        Positional arguments for the solver of DifferentialEquations.jl
#     kwargs:      Keyword arguments for the solver of DifferentialEquations.jl
# 
# Returns
#
#     A vector of vectors [fi(t), fa(t), fm(t), fs(t), fZ(t), fd(t)] 
#     for each t in `times`. It throws an error if the solver did not succeed
#####################################################################################
function integrate_model(
    ic::Vector{Float64},
	base_params::Vector{Float64},
	tspan::NTuple{2,Float64};
    times::Vector{Float64}=[tspan[2],],
	args::Tuple=(QNDF(),),
    kwargs::NamedTuple=(
		dense=false,
		reltol=1.0e-10,
		abstol=1.0e-12,
		verbose=SciMLLogging.None(),
	),
)::Vector{Vector{Float64}}

	ρ_cell, Z, a, h, Z_Si, Z_C = base_params

	ΓUVB          = compute_UVB(a)
	ΓLWB          = compute_LWB(a)
	η_diss, η_ion = compute_η(tspan[2], Z)
	R, Zsn        = compute_recycled_fractions(Z)
	fd_max, c_dg  = accretion_parameters(Z_Si, Z_C, Z)

	parameters = [ρ_cell, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h, fd_max, c_dg]

    sol = solve(ODEProblem(
		ode_function,
		ic,
		tspan,
		parameters,
	), args...; saveat=times, kwargs...)

	# Check the return code of the solver
	if !occursin("Success", string(sol.retcode))
		error("integrate_model: the solver failed with return code $(sol.retcode)")
	end

    return sol(times).u

end;

# ╔═╡ 0fc36dec-9f45-4f49-a44f-7884ec33d7be
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Star formation rate passed to $\texttt{AREPO}$

After the integration, the model gives $\texttt{AREPO}$ the SFR of the cell, which $\texttt{AREPO}$ uses to compute the probability of forming a star particle of mass $m_*$ during a time step $\Delta t$

$\begin{equation}
    p = \frac{M_\mathrm{cell}}{m_*} \left[ 1 - \exp\left( -\frac{\mathrm{SFR} \, \Delta t}{M_\mathrm{cell}} \right) \right] \, .
\end{equation}$

The expected stellar mass formed in a time step is $p \, m_* \simeq \mathrm{SFR} \, \Delta t$, so this SFR sets the mass that $\texttt{AREPO}$ turns into star particles.

In the ODEs, the stellar fraction follows

$\begin{equation}
    \frac{\mathrm{d}}{\mathrm{d}t} f_s = \psi - R \, \psi \, ,
\end{equation}$

i.e. $f_s$ is the mass that remains in stars after the instantaneous return of the mass of the massive stars, which goes back to the gas phases ($R \, (1 - Z_\mathrm{SN}) \, \psi$ to $f_i$, and $R \, Z_\mathrm{SN} \, \psi$ to $f_Z$ and $f_d$).

The code $\texttt{AREPO}$ uses one of two mechanisms to account for the mass returned by stars, never both:

*  Without stellar evolution, the star particles do not return mass, and the SFR is reduced by the mass fraction returned under the IRA, $\beta$ ($\texttt{FactorSN}$): $\mathrm{SFR} = (1 - \beta) \, M_\mathrm{cold} / t_*$ ([Springel2003](https://doi.org/10.1046/j.1365-8711.2003.06206.x)).
*  With stellar evolution ($\texttt{GFM\_STELLAR\_EVOLUTION = 0}$, the configuration of our simulations), each star particle is a stellar population that returns mass and metals to the gas (SNe II, AGB stars, and SNe Ia) following its evolution, and the SFR is the gross rate, $\mathrm{SFR} = M_\mathrm{cold} / t_*$, without the factor $1 - \beta$.

So, if we passed the net rate $f_s \, M_\mathrm{cell} / t_i$, the star particles would get only a fraction $(1 - R)$ of the stellar mass formed in the model, and then they would return a similar fraction again through the stellar evolution of $\texttt{AREPO}$, which means that the mass of the massive stars would be returned twice.

To avoid this, we pass the gross rate, i.e. the mass formed by $\psi$ during the integration

$\begin{equation}
    \mathrm{SFR} = \frac{M_\mathrm{cell}}{t_i} \int_0^{t_i} \psi \, \mathrm{d}t \, ,
\end{equation}$

where $t_i$ is the integration time. As $R$ is constant during the integration, integrating the equation for $f_s$ with $f_s(0) = 0$ gives

$\begin{equation}
    f_s(t_i) = (1 - R) \int_0^{t_i} \psi \, \mathrm{d}t \, ,
\end{equation}$

so

$\begin{equation}
    \mathrm{SFR} = \frac{f_s(t_i) \, M_\mathrm{cell}}{(1 - R) \, t_i} \, .
\end{equation}$

The return of mass within the ODEs is still needed, as it sets the gas phases and the metals during the integration, and closes the mass budget of the six fractions. The factor $1 / (1 - R)$ only concerns the SFR given to $\texttt{AREPO}$, and it makes our SFR directly comparable with the one of the fiducial model in $\texttt{AREPO}$ (see the section Efficiency per free-fall time).
"""
  ╠═╡ =#

# ╔═╡ 24920a1b-355c-43e1-aff2-028e6f25584c
# ╠═╡ skip_as_script = true
#=╠═╡
md"# Extras"
  ╠═╡ =#

# ╔═╡ ad67d733-4365-4b37-ab45-1b936c80d617
begin
	#################################################################################
	# Example cell for the tests below
	#
	# Solar metallicity (Zsun) with the abundance pattern of Asplund (2009), 
	# a two-phase split typical of star-forming gas (90% of cold neutral clouds, 
	# and a fully ionized hot phase), and the gas initial conditions of the model
	#################################################################################

	# Solar metallicity used in AREPO
    const Zsun = 0.0127

	# Metallicity and abundances of Si and C
	const EXAMPLE_Z    = Zsun
	const EXAMPLE_X_Si = Z_Si_A09 * Zsun / Z_A09
	const EXAMPLE_X_C  = Z_C_A09 * Zsun / Z_A09

	# Ionized and neutral fractions of the gas
	const EXAMPLE_FI_GAS = 0.1
	const EXAMPLE_FN_GAS = 1.0 - EXAMPLE_FI_GAS

	# Initial dust fraction (first entry)
	const EXAMPLE_FD = initial_dust_fraction(EXAMPLE_X_Si, EXAMPLE_X_C, EXAMPLE_Z)

	# [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
	const EXAMPLE_IC = [
		(1.0 - EXAMPLE_Z) * EXAMPLE_FI_GAS,
		(1.0 - EXAMPLE_Z) * EXAMPLE_FN_GAS,
		0.0,
		0.0,
		EXAMPLE_Z - EXAMPLE_FD,
		EXAMPLE_FD,
	]

	# Parameters
	example_base_params(ρ_cell) = [
		ρ_cell,       # [cm^(-3)]
		EXAMPLE_Z,    # [dimensionless]
		1.0,          # [dimensionless]
		3.0e20,       # [cm]
		EXAMPLE_X_Si, # [dimensionless]
		EXAMPLE_X_C,  # [dimensionless]
	]
end;

# ╔═╡ eb211ddb-1ba6-4ef0-8712-4aa7c4320c26
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Condensation rate coefficient compared with LOZ25, for a cell with Z = Zsun and 
# the solar abundance pattern (Asplund et al. 2009, scaled to Zsun).
#####################################################################################
let
	# LOZ25 (Wolfire et al. 2008 normalisation and Z_eff of Glover & Jappsen 2007)
	# Eq. 21: 2 R⊙ Cρ (Z + Z_eff) (f_a + f_m + f_i) / (m_p Z⊙)
	# f_a + f_m + f_i = 1.0
	R_sun = 3.5e-17u"cm^3 * s^-1"
	Z_eff = 1e-3 * Zsun
	c_old = ustrip(
		Unitful.NoUnits, 
		2 * R_sun * Cρ * (Zsun + Z_eff) / (u"mp" * Zsun) / cc_unit
	)

	# This work, with the initial dust fraction of a solar cell
	# c_cond (f_d + f_d^eff), with f_d = f_d(0)
	fd0   = initial_dust_fraction(
		Z_Si_A09 * Zsun / Z_A09, Z_C_A09 * Zsun / Z_A09, 
		Zsun,
	)
	c_new = c_cond * (fd0 + fd_eff)

	# Coefficients in [mp^-1 Myr^-1 cm^3], and relative change of the rate
	println("LOZ25           = $(round(c_old; sigdigits=3))")
	println("this_work       = $(round(c_new; sigdigits=3))")
	println("relative_change = $(round(c_new / c_old - 1.0; sigdigits=2))")
	
end
  ╠═╡ =#

# ╔═╡ 3a2fce00-7fee-4c04-9ed0-8e7cdc5bcc0a
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Dust growth time scale, τ_dg = 1 / (C_dg ρ_cell), for Z = Zsun 
# with the solar abundance pattern (Asplund et al. 2009, scaled to Zsun), 
# at ρ_cell = 10^3 mp cm^(-3)
#####################################################################################
let
	# Si and C mass fractions for the solar abundance pattern, scaled to Zsun
	Z_Si = Z_Si_A09 * Zsun / Z_A09
	Z_C  = Z_C_A09 * Zsun / Z_A09

	# Dust growth constant [Myr^(-1) * mp^(-1) * cm^3]
	_, c_dg = accretion_parameters(Z_Si, Z_C, Zsun)

	# Cell density [mp * cm^(-3)]
	ρ_cell = 1.0e3

	println("τ_dg(1.0e3) = $(round(1.0 / (c_dg * ρ_cell); sigdigits=3) * u"Myr")")
end
  ╠═╡ =#

# ╔═╡ 3871738a-a126-4996-ad95-9674564a4931
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Dust-to-metal ratio of the initial condition, f_d(0) / Z, for Z = Zsun with the
# solar abundance pattern (Asplund et al. 2009, scaled to Zsun)
#####################################################################################
let
	# Si and C mass fractions for the solar abundance pattern, scaled to Zsun
	Z_Si = Z_Si_A09 * Zsun / Z_A09
	Z_C  = Z_C_A09 * Zsun / Z_A09

	println(
		"DZ(⊙) = $(round(initial_dust_fraction(Z_Si, Z_C, Zsun) / Zsun; sigdigits=2))"
	)
end
  ╠═╡ =#

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
#     base_params: Parameters for the ODEs, [ρ_cell, Z, a, h, Z_Si, Z_C]
#     it:          Integration time [Myr]
#
# Returns
#
#     The stiffness ratio
#####################################################################################
function stiffness_ratio(
	ic::Vector{Float64},
	base_params::Vector{Float64},
	it::Float64,
)::Float64

	# Parameters for the ODEs
	# ρ_cell: Density [cm^(-3)]
	# Z:      Metallicity [dimensionless]
	# a:      Scale factor [dimensionless]
	# h:      Cell size [cm]
	# Z_Si:   Mass fraction of Si [dimensionless]
	# Z_C:    Mass fraction of C [dimensionless]
	ρ_cell, Z, a, h, Z_Si, Z_C = base_params

	ΓUVB          = compute_UVB(a)
	ΓLWB          = compute_LWB(a)
	η_diss, η_ion = compute_η(it, Z)
	R, Zsn        = compute_recycled_fractions(Z)
	fd_max, c_dg  = accretion_parameters(Z_Si, Z_C, Z)

	parameters = [ρ_cell, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h, fd_max, c_dg]

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
stiffness_ratio(EXAMPLE_IC, example_base_params(1.0), 5.0)
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
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Compute the Lyapunov spectrum
#
# Arguments
# 
#     ic:         Initial conditions, [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
#     base_params: Parameters for the ODEs, [ρ_cell, Z, a, h, Z_Si, Z_C]
#     it:          Integration time [Myr]
#
# Returns
#
#     The Lyapunov spectrum
#####################################################################################
function lyapunov_spectrum(
	ic::Vector{Float64},
	base_params::Vector{Float64},
	it::Float64,
)::Vector{Float64}

	# Parameters for the ODEs
	# ρ_cell: Density [cm^(-3)]
	# Z:      Metallicity [dimensionless]
	# a:      Scale factor [dimensionless]
	# h:      Cell size [cm]
	# Z_Si:   Mass fraction of Si [dimensionless]
	# Z_C:    Mass fraction of C [dimensionless]
	ρ_cell, Z, a, h, Z_Si, Z_C = base_params

	ΓUVB          = compute_UVB(a)
	ΓLWB          = compute_LWB(a)
	η_diss, η_ion = compute_η(it, Z)
	R, Zsn        = compute_recycled_fractions(Z)
	fd_max, c_dg  = accretion_parameters(Z_Si, Z_C, Z)

	parameters = [ρ_cell, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h, fd_max, c_dg]
	
	ds = CoupledODEs(system!, ic, parameters)

	# Compute the Lyapunov spectrum
	return lyapunovspectrum(ds, 100; Δt=1e-6)

end;
  ╠═╡ =#

# ╔═╡ b1e67fcd-f70b-4dd2-87e5-74bef1fa201e
# ╠═╡ skip_as_script = true
#=╠═╡
lyapunov_spectrum(EXAMPLE_IC, example_base_params(10.0), 10.0)
  ╠═╡ =#

# ╔═╡ d1f3d076-b9d6-4b92-b25c-4155b5fc464b
# ╠═╡ skip_as_script = true
#=╠═╡
lyapunov_spectrum( EXAMPLE_IC, example_base_params(1000.0), 10.0,)
  ╠═╡ =#

# ╔═╡ f336b083-99c7-4345-8ada-78b166a98abe
# ╠═╡ skip_as_script = true
#=╠═╡
md"The model, for a typical set of initial conditions, is chaotic (maximum lyapunov exponent > 0)."
  ╠═╡ =#

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

where $M_g$ is the gas mass under consideration(in our case, this is the cell mass $M_\text{cell}$), $\text{SFR}$ the star formation rate (in our case, this refers to the probabilistic star formation rate $\text{SFR}_p$), and $t_\text{ff}$ is the free-fall time, defined as

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
	\text{SFR}_p = \frac{f_s \, M_\text{cell}}{(1 - R) \, t_i} \, ,
\end{equation}$

where $f_s$ is the stellar fraction from the ODEs and $t_i$ is the integration time (see the section Star formation rate passed to $\texttt{AREPO}$ for the factor $1 - R$).

So, in our model, we end up with

$\begin{equation}
	\epsilon_\text{ff}^{\,\text{SFM}} = c \, \frac{f_s}{(1 - R) \, t_i \, \sqrt{\rho_\text{cell}}} \, .
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
To find the per-equation and global equilibrium we take the ODEs and set the derivatives to $0$.
"""
  ╠═╡ =#

# ╔═╡ e8c51ae9-ded4-450a-80f9-c484c2b991d1
#####################################################################################
# Evaluate the ODEs
# 
# Arguments
# 
#     fractions:  Mass fractions, [fi, fa, fm, fs, fZ, fd]
#     base_params: Parameters for the ODEs, [ρ_cell, Z, a, h, Z_Si, Z_C]
#     it:          Integration time [Myr]
# 
# Returns
#
#     The ODE evaluated at `fractions` and `base_params`
#####################################################################################
function equilibrium(
    fractions::Vector{Float64},
	base_params::Vector{Float64},
	it::Float64,
)::Vector{Float64}

	ρ_cell, Z, a, h, Z_Si, Z_C = base_params

	ΓUVB          = compute_UVB(a)
	ΓLWB          = compute_LWB(a)
	η_diss, η_ion = compute_η(it, Z)
	R, Zsn        = compute_recycled_fractions(Z)
	fd_max, c_dg  = accretion_parameters(Z_Si, Z_C, Z)

	parameters = [ρ_cell, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h, fd_max, c_dg]

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
#     fractions:   Mass fractions, [fi, fa, fm, fs, fZ, fd]
#     base_params: Parameters for the ODEs, [ρ_cell, Z, a, h, Z_Si, Z_C]
#     it:          Integration time [Myr]
#     limit:       Equilibrium criteria.
# 
# Returns
# 
#     If the system pass the criteria (the derivatives are close enough to 0)
#####################################################################################
function test_equilibrium(
    fractions::Vector{Float64},
	base_params::Vector{Float64},
	it::Float64;
	limit::Float64=1e-10,
)::Vector{Bool}

	# Change of each fraction during the integration time, relative to `limit`
    return abs.(equilibrium(fractions, base_params, it) .* it) .< limit

end;

# ╔═╡ cd0e57c7-3fa0-4875-9a38-249b3a0bb134
# ╠═╡ skip_as_script = true
#=╠═╡
let
	# Example cell at ρ_cell = 1 mp cm^(-3)
	base_params = example_base_params(1.0)

	# Integration time [Myr]
	it = 10.0

	# Fractions at the end of the integration
	fractions = integrate_model(EXAMPLE_IC, base_params, (0.0, it))[end]

	# Phases whose change over `it`, at the current rates, is below 10^(-3)
	test_equilibrium(fractions, base_params, it; limit=1e-3)
end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
ChaosTools = "608a59af-f2a3-5ad4-90b4-758bdf3122a7"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
OrdinaryDiffEqBDF = "6ad6398a-0878-4a85-9266-38940aa047c8"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"

[compat]
CSV = "~1.3.0"
ChaosTools = "~3.5.4"
DataFrames = "~1.8.3"
DifferentialEquations = "~8.1.1"
Interpolations = "~0.16.3"
OrdinaryDiffEqBDF = "~2.4.14"
PlutoUI = "~0.7.84"
QuadGK = "~2.11.3"
Symbolics = "~7.48.0"
TikzPictures = "~3.5.2"
Unitful = "~1.29.0"
UnitfulAstro = "~1.2.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.3"
manifest_format = "2.1"
project_hash = "4c36dc326a9fcc9c390f4432247a21dbb6dc1d8d"

[[deps.ADTypes]]
deps = ["PrecompileTools"]
git-tree-sha1 = "629de23e1c16911b439dabd2303c08af9575b226"
registries = "General"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.24.0"
weakdeps = ["ChainRulesCore", "ConstructionBase", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

[[deps.AMD]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "344f3d221a6bee75925b5c97a30cfd870f12a138"
registries = "General"
uuid = "14f7f29c-3bd6-536c-9a0b-7339e30b5a3e"
version = "0.5.4"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
registries = "General"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "e71ee7b4aa06b045259a7d6101e1cb45ad140bce"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.1"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
registries = "General"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "8c0bc12397dabf5d2dfc4b5763ec463854cc44d3"
registries = "General"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.46"

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
git-tree-sha1 = "eec8a0ad4e4b673afba213931504c360fffb8f97"
registries = "General"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.8.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
registries = "General"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
registries = "General"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "1aec1ff0dcaa83a484e7f72097562fa3102c6722"
registries = "General"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.30.2"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceAMDGPUExt = "AMDGPU"
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceFillArraysExt = "FillArrays"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceGPUArraysCoreTrackerExt = ["GPUArraysCore", "Tracker"]
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
    FillArrays = "1a297f60-69ca-5386-bcde-b61e274b549b"
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
registries = "General"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
registries = "General"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
registries = "General"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.BinaryHeaps]]
deps = ["PrecompileTools"]
git-tree-sha1 = "0adb6d8d0a4bae93a07d79c8d59eaf617f4c59da"
registries = "General"
uuid = "b2a6c25c-c996-4615-94ab-6e519ffb8690"
version = "1.1.0"

[[deps.BracketingNonlinearSolve]]
deps = ["CommonSolve", "ConcreteStructs", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLLogging"]
git-tree-sha1 = "b1b04c236ff7ed92e69d64b620448b3da00933d4"
registries = "General"
uuid = "70df07ce-3d50-431d-a3e7-ca6ddb60ac1e"
version = "1.12.8"
weakdeps = ["ChainRulesCore", "ForwardDiff"]

    [deps.BracketingNonlinearSolve.extensions]
    BracketingNonlinearSolveChainRulesCoreExt = ["ChainRulesCore", "ForwardDiff"]
    BracketingNonlinearSolveForwardDiffExt = "ForwardDiff"

[[deps.BranchAndPrune]]
deps = ["AbstractTrees"]
git-tree-sha1 = "8bea09a63c60b6b5200e613379947e7c8527992b"
registries = "General"
uuid = "d3bc4f2e-91e6-11e9-365e-cd067da536ce"
version = "0.2.2"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
registries = "General"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "0db159493e83ecb2d122f9c40138d89b70d43963"
registries = "General"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.1.1"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
registries = "General"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.CSV]]
deps = ["CodecZlib", "DataStrings", "Dates", "Downloads", "Durations", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "Printf", "Tables", "Unicode"]
git-tree-sha1 = "feb15aa037b0004ed9a988a05bb6f5833a670ac7"
registries = "General"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "1.3.0"

    [deps.CSV.extensions]
    CSVDataDecimalsExt = "DataDecimals"
    CSVFilePathsBaseExt = "FilePathsBase"
    CSVInlineStringsExt = "InlineStrings"

    [deps.CSV.weakdeps]
    DataDecimals = "3e2245cb-6932-498d-a7cf-39d39de8bde3"
    FilePathsBase = "48062228-2e41-5def-b9a4-89aafe57970f"
    InlineStrings = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2ac646d71d0d24b44f3f8c84da8c9f4d70fb67df"
registries = "General"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.4+0"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9cb23bbb1127eefb022b022481466c0f1127d430"
registries = "General"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.2"

[[deps.ChainRulesCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "295fa9b898db109bc56df62c0fc20bd105cec00e"
registries = "General"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.8"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChaosTools]]
deps = ["Combinatorics", "DSP", "DataStructures", "Distances", "Distributions", "DynamicalSystemsBase", "IntervalRootFinding", "LinearAlgebra", "LombScargle", "Neighborhood", "Optim", "ProgressMeter", "Random", "Reexport", "Roots", "SpecialFunctions", "Statistics", "StatsBase"]
git-tree-sha1 = "12f9540406976e23eb2d7396d5bbe6e0e7be7874"
registries = "General"
uuid = "608a59af-f2a3-5ad4-90b4-758bdf3122a7"
version = "3.5.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "970758a3d591a2a5c2a907c53f2e2f8c1b1d3537"
registries = "General"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.9"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "af523e59de94de05b9b367a2ca4d63233aacd222"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.5"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
registries = "General"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
deps = ["PrecompileTools"]
git-tree-sha1 = "6c389fa857f6ca5a95474b52a52023fd77f24cb7"
registries = "General"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.14"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
registries = "General"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
deps = ["PrecompileTools"]
git-tree-sha1 = "7c4ea0adfb7238bf4678aa36325863fab6163f6f"
registries = "General"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.2.2"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
registries = "General"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
registries = "General"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
registries = "General"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a72c3b5ce6d2a477f55b5c9b8756e91e695c67f6"
registries = "General"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.8"

[[deps.ConstructionBase]]
git-tree-sha1 = "b1bc10667320dd5367aff450a992c951ea38fbbc"
registries = "General"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.1"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.CoreMath]]
deps = ["CoreMath_jll"]
git-tree-sha1 = "8c0480f92b1b1796239156a1b9b1bfb1b39499b4"
registries = "General"
uuid = "b7a15901-be09-4a0e-87d2-2e66b0e09b5a"
version = "0.1.0"

[[deps.CoreMath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a692a4c1dc59a4b8bc0b6403876eb3250fde2bc3"
registries = "General"
uuid = "a38c48d9-6df1-5ac9-9223-b6ada3b5572b"
version = "0.1.0+0"

[[deps.Crayons]]
git-tree-sha1 = "54b76cbb40d9a0f5368c880725b2f141da77c94f"
registries = "General"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.2.0"

[[deps.DSP]]
deps = ["Bessels", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "dd73b345b24e52d7b57550b69113261ef2d79c92"
registries = "General"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.8.7"
weakdeps = ["OffsetArrays"]

    [deps.DSP.extensions]
    OffsetArraysExt = "OffsetArrays"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
registries = "General"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fd798b61b03b15dda0715aac807885397e6ded88"
registries = "General"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.3"

[[deps.DataStrings]]
git-tree-sha1 = "b8032263d364177a0b192ea15f277df4e6ab07e2"
registries = "General"
uuid = "48204bd6-5611-42ea-b167-fc1d713f9be2"
version = "1.1.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
registries = "General"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
registries = "General"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
registries = "General"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "BracketingNonlinearSolve", "DocStringExtensions", "FastBroadcast", "FastClosures", "FastPower", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "RespecializeParams", "SciMLBase", "SciMLLogging", "SciMLOperators", "SciMLStructures", "Setfield", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "47230845a4944d292c07079838bdfee25bffebf8"
registries = "General"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "7.21.6"

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

[[deps.DiffEqNoiseProcess]]
deps = ["CommonSolve", "DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "PoissonRandom", "PrecompileTools", "QuadGK", "Random", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "891caab6046d57486f779b529bfbfebcdeada50a"
registries = "General"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.36.5"

    [deps.DiffEqNoiseProcess.extensions]
    DiffEqNoiseProcessOptimExt = "Optim"
    DiffEqNoiseProcessReverseDiffExt = "ReverseDiff"

    [deps.DiffEqNoiseProcess.weakdeps]
    Optim = "429524aa-4258-5aef-a3af-852621145aeb"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
registries = "General"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "79a2aca180a85c690c58a020d47b426954b590f8"
registries = "General"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.16.0"

[[deps.DifferentialEquations]]
deps = ["OrdinaryDiffEq", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "d167917674ca92bbf51b56a3e00d019792695110"
registries = "General"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "8.1.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "4b6690f3770daf596a0a092203a2e6130d34c53b"
registries = "General"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.22"

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
registries = "General"
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
git-tree-sha1 = "3c8a0a9a6d4a10bdfb6b751bd2b6051ed3e25fd4"
registries = "General"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.127"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
registries = "General"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.DomainSets]]
deps = ["CompositeTypes", "FunctionMaps", "IntervalSets", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "6ac879a64381881dc24f96a15e13e80204439906"
registries = "General"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.8.3"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"
    DomainSetsRandomExt = "Random"
    DomainSetsRecipesBaseExt = "RecipesBase"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.Durations]]
deps = ["Dates"]
git-tree-sha1 = "407b209d76168dcc4355f9dcc4d3b3aaa60dfb3e"
registries = "General"
uuid = "41ac09f5-0a95-4810-aa18-c88bf48ee130"
version = "1.4.1"

    [deps.Durations.extensions]
    DurationsArrowExt = "Arrow"
    DurationsTimeZonesExt = "TimeZones"

    [deps.Durations.weakdeps]
    Arrow = "69666777-d1a9-59fb-9406-91d4454c9d45"
    TimeZones = "f269a46b-ccf7-5d73-abea-4c690281aa53"

[[deps.DynamicPolynomials]]
deps = ["LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "StarAlgebras"]
git-tree-sha1 = "b95c4325301d86a6a1b59790b71f7269a8702dd5"
registries = "General"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.8"

[[deps.DynamicalSystemsBase]]
deps = ["DiffEqNoiseProcess", "ForwardDiff", "LinearAlgebra", "OrdinaryDiffEqTsit5", "Random", "RecursiveArrayTools", "Reexport", "Roots", "SciMLBase", "SparseArrays", "StateSpaceSets", "Statistics", "StochasticDiffEqHighOrder", "SymbolicIndexingInterface"]
git-tree-sha1 = "de7d92a7f648aff91688262c7be6d1cb73c7a389"
registries = "General"
uuid = "6e36e845-645a-534a-86f2-f5d4aa5a06b4"
version = "3.19.3"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
registries = "General"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.EnzymeCore]]
git-tree-sha1 = "d88c803e26799193db42e47f6e1e28eeb0322e4d"
registries = "General"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.27"
weakdeps = ["Adapt", "ChainRulesCore"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"
    EnzymeCoreChainRulesCoreExt = "ChainRulesCore"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f2659a735f9360fedfa29925812380f322679e93"
registries = "General"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.10.0+0"

[[deps.ExprTools]]
git-tree-sha1 = "d2e49e7efd29719d6f28b891b0e0e159daa9d2b4"
registries = "General"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.11"

[[deps.ExpressionExplorer]]
git-tree-sha1 = "fdcc9d4a253e0041d37c95487f9ae027c3ec002d"
registries = "General"
uuid = "21656369-7473-754a-2065-74616d696c43"
version = "1.1.6"

[[deps.ExproniconLite]]
git-tree-sha1 = "250aa1562e840b38c72b44a9ebef0f39f6b6e950"
registries = "General"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.16"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
registries = "General"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6866aec60ef98e3164cd8d6855225684207e9dff"
registries = "General"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.12+0"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "093bf5e9ea63b36032240033ec2f149088fc0926"
registries = "General"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "1.5.1"

    [deps.FastBroadcast.extensions]
    FastBroadcastPolyesterExt = "Polyester"
    FastBroadcastStaticExt = "Static"
    FixedSizeArraysExt = "FixedSizeArrays"

    [deps.FastBroadcast.weakdeps]
    FixedSizeArrays = "3821ddf9-e5b5-40d5-8e25-6813ab96b5e2"
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"
    Static = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
registries = "General"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastPower]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a6bcbd090fc6b6f022fc0ad96f77ba31b9577282"
registries = "General"
uuid = "a4df4552-cc26-4903-aec0-212e50a0e84b"
version = "1.5.1"

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

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "831f86fde78ab7f6470d6ba4d1161a11fb6b4ec0"
registries = "General"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.18.2"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FindFirstFunctions]]
deps = ["PrecompileTools"]
git-tree-sha1 = "93eb9afa3a53b29d8a49714b83878225402b3281"
registries = "General"
uuid = "64ca27bc-2ba2-4a57-88aa-44e436879224"
version = "3.4.1"

    [deps.FindFirstFunctions.extensions]
    FindFirstFunctionsEnzymeCoreExt = "EnzymeCore"
    FindFirstFunctionsMooncakeExt = "Mooncake"

    [deps.FindFirstFunctions.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "5031f23e040bf17082e5b52422d77b5e844eefb1"
registries = "General"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.33.0"

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
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "e1e9ba3a8f5d71e6ae698f138a1362aeabb42e9f"
registries = "General"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.18.3+0"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a55669aa77ae8bbdc7310c4027112500ff5ec874"
registries = "General"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.4.7"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d95a7a3a7e0f2fffc6291bd33035830aab01310b"
registries = "General"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.15.0+0"

[[deps.FunctionMaps]]
deps = ["CompositeTypes", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "31bd99a57edf98990d1c21486032963955450e8d"
registries = "General"
uuid = "a85aefff-f8ca-4649-a888-c8e5398bc76c"
version = "0.1.2"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
registries = "General"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "2bcce3ad6f6977d617928d7707fdc86ac83cce03"
registries = "General"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "1.13.0"

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
git-tree-sha1 = "f60436514e6d73f8ee069ae4ee0ff038fe96b69d"
registries = "General"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.2"

[[deps.Gamma]]
deps = ["LogExpFunctions"]
git-tree-sha1 = "462ab5e6578928370f632705e96af82f6c1c09a1"
registries = "General"
uuid = "a0844989-3bd2-4988-8bea-c9407ab0941b"
version = "1.2.1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
registries = "General"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "b0036b392358c80d2d2124746c2bf3d48d457938"
registries = "General"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.4+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "69ffb934a5c5b7e086a0b4fee3427db2556fba6e"
registries = "General"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.16+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Inflate", "LinearAlgebra", "Random", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "dbb2b976f605b220dfe139d6db9f3859680da09e"
registries = "General"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.15.0"
weakdeps = ["Distributed", "SharedArrays"]

    [deps.Graphs.extensions]
    GraphsSharedArraysExt = "SharedArrays"

[[deps.HarfBuzz_ICU_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "HarfBuzz_jll", "ICU_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "6ccbc4fdf65c8197738c2d68cc55b74b19c97ac2"
registries = "General"
uuid = "655565e8-fb53-5cb3-b0cd-aec1ca0647ea"
version = "2.8.1+0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
registries = "General"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["Gamma", "LinearAlgebra"]
git-tree-sha1 = "c28dceeadc1c9c699656a3ac9d63f79cfaa17472"
registries = "General"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.31"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
registries = "General"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
registries = "General"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.ICU_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "20b6765a3016e1fca0c9c93c80d50061b94218b7"
registries = "General"
uuid = "a51ab1cf-af8e-5615-a023-bc2c838bba6b"
version = "69.1.0+0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
registries = "General"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
registries = "General"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InlineStrings]]
git-tree-sha1 = "6e0a59a838f58259f49474e913e39c23e8ca1366"
registries = "General"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "2.0.1"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"
    SerializationExt = "Serialization"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
    Serialization = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "c72458f1962faeb003bf23cbdb75164fe6280906"
registries = "General"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.4"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
registries = "General"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "48922d06068130f87e43edef52382e6a94305ae6"
registries = "General"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.3"
weakdeps = ["ForwardDiff", "Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "CoreMath", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "4f051d490760061dcc650a3c6c641bb0c06bb182"
registries = "General"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.13"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticIrrationalConstantsExt = "IrrationalConstants"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticMakieExt = "Makie"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    IrrationalConstants = "92d709cd-6900-40b7-9082-c6be49f344b6"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalRootFinding]]
deps = ["BranchAndPrune", "ForwardDiff", "IntervalArithmetic", "LinearAlgebra", "Reexport", "StaticArrays"]
git-tree-sha1 = "e0bbe14c3e38ddf829f130e698dd53f5d7400cce"
registries = "General"
uuid = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
version = "0.6.3"

[[deps.IntervalSets]]
git-tree-sha1 = "0db6aea5b64caa1c47e9fdd27394f0296f05e8bf"
registries = "General"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.15"

    [deps.IntervalSets.extensions]
    IntervalSetsMakieExt = "Makie"
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

    [deps.IntervalSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.InverseFunctions]]
git-tree-sha1 = "5cc6ad3e2bfc75a11eed92632ee341d37f4a4536"
registries = "General"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.18"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
registries = "General"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "d22e8bdb3f03c300a103ff23c0dd5a654c6701de"
registries = "General"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.7"

[[deps.IterTools]]
git-tree-sha1 = "123311cf9834c6ee6a8347bba9af68f255615391"
registries = "General"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.11.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
registries = "General"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
registries = "General"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
registries = "General"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e06b451acfb1c6fa5ea753b5042207226271f455"
registries = "General"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.2.1+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.JumpProcesses]]
deps = ["ADTypes", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "PoissonRandom", "PrecompileTools", "Random", "RecursiveArrayTools", "SciMLBase", "SimpleNonlinearSolve", "StaticArrays", "SymbolicIndexingInterface"]
git-tree-sha1 = "8bab411bfa528f97ac0f8961288b37c567ed5e80"
registries = "General"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "10.1.1"

    [deps.JumpProcesses.extensions]
    JumpProcessesForwardDiffExt = "ForwardDiff"
    JumpProcessesKernelAbstractionsExt = ["Adapt", "KernelAbstractions"]
    JumpProcessesOrdinaryDiffEqCoreExt = "OrdinaryDiffEqCore"

    [deps.JumpProcesses.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    OrdinaryDiffEqCore = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "b13076e4fff6ea0d0a5ddd9b939a11ca2027755a"
registries = "General"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.10.12"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d7c65411dda4bfa0beb7490476122b4f805cf31"
registries = "General"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.3.0+0"

[[deps.LHLFactorization]]
deps = ["LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "51d99d9892d0f95b4c9d4d7310785971f2d48298"
registries = "General"
uuid = "2faa5264-e118-4071-8864-e10559f68c7c"
version = "2.2.2"

    [deps.LHLFactorization.extensions]
    LHLFactorizationPolyesterExt = "Polyester"
    LHLFactorizationSparseExt = ["PureKLU", "SparseArrays"]

    [deps.LHLFactorization.weakdeps]
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"
    PureKLU = "0c0d3e7f-3a8b-4f7e-b6f1-9a4d2e7c1f01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e5b100780d4d30d63b4618d7930d48af409c1772"
registries = "General"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "23.1.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
registries = "General"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f88f3ccef05a6a72a0cf0ed417c8fd68530f4ab2"
registries = "General"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "1.0.0"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.22.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.7+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "OpenSSL_jll", "Zlib_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.104+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "27ecae93dd25ee0909666e6835051dd684cc035e"
registries = "General"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+2"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
registries = "General"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b8a91037379e1fd661e9b5ff3d7a770424b1f9f2"
registries = "General"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.4+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "aebd334d06cee9f24cea70bd19a39749daf73881"
registries = "General"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.3+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9e887184f306b28dbfcff88cefbb0ec80f1aa9e0"
registries = "General"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.4+0"

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "PrecompileTools", "ReactantCore", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "53cf132c453dcfcf10f1cb9a701a07c99b7e842e"
registries = "General"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.20"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Printf"]
git-tree-sha1 = "a09a85e94e681def2446752cddaf9401a4d82b6d"
registries = "General"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.8.2"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.13.0"

[[deps.LinearSolve]]
deps = ["AMD", "ArrayInterface", "ConcreteStructs", "DocStringExtensions", "EnumX", "GPUArraysCore", "InteractiveUtils", "Krylov", "LHLFactorization", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "OpenBLAS_jll", "PrecompileTools", "Preferences", "PureKLU", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLLogging", "SciMLOperators", "SciMLStructures", "Setfield", "SparseArrays", "SparseColumnPivotedQR", "StaticArraysCore"]
git-tree-sha1 = "e3e84c5c330ea65568a85aef00031e2b870b0e99"
registries = "General"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "5.18.3"

    [deps.LinearSolve.extensions]
    LinearSolveAMDGPUExt = "AMDGPU"
    LinearSolveAMGXExt = ["AMGX", "CUDA"]
    LinearSolveAlgebraicMultigridExt = "AlgebraicMultigrid"
    LinearSolveArnoldiMethodExt = "ArnoldiMethod"
    LinearSolveArpackExt = "Arpack"
    LinearSolveBLISExt = ["blis_jll", "LAPACK_jll"]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCKTSOExt = "CKTSO"
    LinearSolveCUDAExt = ["cuSOLVER"]
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveCUSOLVERRFExt = "CUSOLVERRF"
    LinearSolveChainRulesCoreExt = "ChainRulesCore"
    LinearSolveCliqueTreesExt = "CliqueTrees"
    LinearSolveConjugateGradientsExt = "ConjugateGradients"
    LinearSolveElementalExt = "Elemental"
    LinearSolveEnzymeExt = "EnzymeCore"
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveFastLapackInterfaceExt = "FastLapackInterface"
    LinearSolveForwardDiffExt = "ForwardDiff"
    LinearSolveGinkgoExt = "Ginkgo"
    LinearSolveHSLExt = "HSL"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveJacobiDavidsonExt = "JacobiDavidson"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMUMPSExt = "MUMPS"
    LinearSolveMetalExt = "Metal"
    LinearSolveMooncakeExt = "Mooncake"
    LinearSolvePETScExt = ["PETSc", "SparseMatricesCSR"]
    LinearSolvePETScMPIExt = ["PETSc", "PartitionedArrays", "SparseMatricesCSR"]
    LinearSolveParUExt = "ParU_jll"
    LinearSolvePardisoExt = "Pardiso"
    LinearSolvePartitionedSolversExt = ["PartitionedArrays", "PartitionedSolvers"]
    LinearSolvePureUMFPACKExt = "PureUMFPACK"
    LinearSolveReactantExt = "Reactant"
    LinearSolveRecursiveFactorizationExt = ["RecursiveFactorization", "TriangularSolve"]
    LinearSolveSLATEExt = "SLATE_jll"
    LinearSolveSTRUMPACKExt = "STRUMPACK_jll"
    LinearSolveSparspakExt = "Sparspak"
    LinearSolveSpecializingFactorizationsExt = "SpecializingFactorizations"
    LinearSolveSuperLUDISTExt = "SuperLUDIST"

    [deps.LinearSolve.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    AMGX = "c963dde9-0319-47f5-bf0c-b07d3c80ffa6"
    AlgebraicMultigrid = "2169fc97-5a83-5252-b627-83903c6c433c"
    ArnoldiMethod = "ec485272-7323-5ecc-a04f-4719b315124d"
    Arpack = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CKTSO = "8e543fc2-2ef1-4e2d-a99e-39beaa046845"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    CUSOLVERRF = "a8cc9031-bad2-4722-94f5-40deabb4245c"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    ConjugateGradients = "f59de78d-195d-4e7b-a078-2e47da4c3ad6"
    Elemental = "902c3f28-d1ec-5e7e-8399-a24c3845ee38"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    FastLapackInterface = "29a986be-02c6-4525-aec4-84b980013641"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Ginkgo = "4c8bd3c9-ead9-4b5e-a625-08f1338ba0ec"
    HSL = "34c5aeac-e683-54a6-a0e9-6e0fdc586c50"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    JLArrays = "27aeb0d3-9eb9-45fb-866b-73c2ecf80fcb"
    JacobiDavidson = "11c68b98-9c9b-11e8-267b-bbb95576cead"
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
    PureUMFPACK = "b7e1f0a2-3c4d-4e5f-9a0b-1c2d3e4f5a6b"
    Reactant = "3c362404-f566-11ee-1572-e11a4b42c853"
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
    SLATE_jll = "cb01ff19-587c-516c-bcf9-532d49da298d"
    STRUMPACK_jll = "86fbd0b9-476f-557c-b766-62c724b42d8c"
    SparseMatricesCSR = "a0a7dd2c-ebf4-11e9-1f05-cf50bc540ca1"
    Sparspak = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
    SpecializingFactorizations = "fa08b7a1-13d3-4faf-875d-5cbc1520e3f3"
    SuperLUDIST = "4cd002a6-0da4-410d-a012-232df062f478"
    TriangularSolve = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
    blis_jll = "6136c539-28a5-5bf0-87cc-b183200dce32"
    cuSOLVER = "887afef0-6a32-4de5-add4-7827692ba8fc"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll"]
git-tree-sha1 = "38928f7999753af13d4e13966ae15958ff3a917a"
registries = "General"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.19.1+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "28d620bcdb37a3fb02de9f521928ee453690f414"
registries = "General"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "1.0.3"

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
registries = "General"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.LombScargle]]
deps = ["FFTW", "LinearAlgebra", "Measurements", "Random", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "d64a0ce7539181136a85fd8fe4f42626387f0f26"
registries = "General"
uuid = "fc60dff9-86e7-5f2f-a8a0-edeadbb75bd9"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
registries = "General"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
registries = "General"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
registries = "General"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools", "PrecompileTools"]
git-tree-sha1 = "f2cde0ae772162f20287b803e623966d0d94dea9"
registries = "General"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.8"
weakdeps = ["SparseArrays"]

    [deps.MaybeInplace.extensions]
    MaybeInplaceSparseArraysExt = "SparseArrays"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf"]
git-tree-sha1 = "cb47f69a1cab9dcec7ff4a5d6e163410d6905866"
registries = "General"
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
git-tree-sha1 = "d9d20d6eb2a3290ba5c70a1c013b45305cc39c9e"
registries = "General"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "3f234cf937eb1d58f613440ecb022da784ad34d4"
registries = "General"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.13"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.10.20"

[[deps.MuladdMacro]]
deps = ["PrecompileTools"]
git-tree-sha1 = "283bf85d4a767481dd924dff0eee1735e95f449e"
registries = "General"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.7"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics", "StarAlgebras"]
git-tree-sha1 = "65f1160de342181533efa7061f86539d420ab81d"
registries = "General"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.20"

    [deps.MultivariatePolynomials.extensions]
    LatexifyExt = "Latexify"
    MultivariatePolynomialsChainRulesCoreExt = "ChainRulesCore"

    [deps.MultivariatePolynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "ffd1161361638e5f0bb796a095b0fdc230d6f019"
registries = "General"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.9.3"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "FiniteDiff", "LinearAlgebra"]
git-tree-sha1 = "f96d38936d92d610318dec4d3b5ef37b0373c20f"
registries = "General"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "8.0.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "676c0f1c1c34224d6869675047a54b42fdc7b8c5"
registries = "General"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.5"

[[deps.NearestNeighbors]]
deps = ["AbstractTrees", "Distances", "StaticArrays"]
git-tree-sha1 = "576eb4656529c12e77a46b17c23103dfba9fa570"
registries = "General"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.29"

[[deps.Neighborhood]]
deps = ["Distances", "NearestNeighbors", "Random", "Test"]
git-tree-sha1 = "3efba0551882b18a1e97d27fafe8d01fb70b5b1b"
registries = "General"
uuid = "645ca80c-8b79-4109-87ea-e1f58159d116"
version = "0.2.5"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "NonlinearSolveBase", "NonlinearSolveFirstOrder", "NonlinearSolveQuasiNewton", "NonlinearSolveSpectralMethods", "PrecompileTools", "Preferences", "Reexport", "SciMLBase", "Setfield", "SimpleNonlinearSolve", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "90afae5b2b33dcd1098e4c1a7b6066f18741bf6f"
registries = "General"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "4.32.0"

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
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "Compat", "ConcreteStructs", "DifferentiationInterface", "EnumX", "EnzymeCore", "FastClosures", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "LogExpFunctions", "Markdown", "MaybeInplace", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "RecursiveArrayTools", "RespecializeParams", "SciMLBase", "SciMLJacobianOperators", "SciMLLogging", "SciMLOperators", "SciMLStructures", "Setfield", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "f83b1fb97c02b99390346e13a71a027a59dc2459"
registries = "General"
uuid = "be0214bd-f91f-a760-ac4e-3421ce2b2da0"
version = "2.55.1"

    [deps.NonlinearSolveBase.extensions]
    NonlinearSolveBaseBandedMatricesExt = "BandedMatrices"
    NonlinearSolveBaseChainRulesCoreExt = "ChainRulesCore"
    NonlinearSolveBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    NonlinearSolveBaseForwardDiffExt = "ForwardDiff"
    NonlinearSolveBaseLineSearchExt = "LineSearch"
    NonlinearSolveBaseLinearSolveExt = "LinearSolve"
    NonlinearSolveBaseMooncakeExt = ["ChainRulesCore", "Mooncake"]
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
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConcreteStructs", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLJacobianOperators", "SciMLLogging", "SciMLOperators", "Setfield", "Sobol", "StaticArraysCore"]
git-tree-sha1 = "95a105b647b656dfbd7985ede8eb0eaf4b836505"
registries = "General"
uuid = "5959db7a-ea39-4486-b5fe-2dd0bf03d60d"
version = "2.10.0"

[[deps.NonlinearSolveQuasiNewton]]
deps = ["ArrayInterface", "CommonSolve", "ConcreteStructs", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLLogging", "SciMLOperators", "StaticArraysCore"]
git-tree-sha1 = "db0bca65e9c849333023a0d5827813947d003f09"
registries = "General"
uuid = "9a2c21bd-3a47-402d-9113-8faf9a0ee114"
version = "1.15.3"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveQuasiNewton.extensions]
    NonlinearSolveQuasiNewtonForwardDiffExt = "ForwardDiff"

[[deps.NonlinearSolveSpectralMethods]]
deps = ["CommonSolve", "ConcreteStructs", "LineSearch", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLLogging"]
git-tree-sha1 = "c3e5a641a867b7f37155a78086cf3f22134c088f"
registries = "General"
uuid = "26075421-4e9a-44e1-8bd1-420ed7ad02b2"
version = "1.8.3"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveSpectralMethods.extensions]
    NonlinearSolveSpectralMethodsForwardDiffExt = "ForwardDiff"

[[deps.OffsetArrays]]
git-tree-sha1 = "750d0f2b5aac238d50dfee5f7a44fae3175a3305"
registries = "General"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.18.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "38a93f17e431141c6470bb67a88952a7c4f0e928"
registries = "General"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.34+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "libpng_jll"]
git-tree-sha1 = "215a6666fee6d6b3a6e75f2cc22cb767e2dd393a"
registries = "General"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.5.5+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.8+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
registries = "General"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["ADTypes", "EnumX", "FillArrays", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "PositiveFactorizations", "Printf", "SparseArrays", "Statistics"]
git-tree-sha1 = "81f338ad984b25bc82471db4d235e7b8bf21b8d6"
registries = "General"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "2.3.2"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.OrderedCollections]]
git-tree-sha1 = "b079a6affebde167dad8ffe0d670b53917fbcb68"
registries = "General"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "2.0.3"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "CommonSolve", "DiffEqBase", "DocStringExtensions", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "SciMLBase", "SciMLLogging"]
git-tree-sha1 = "bae9b8e5f1d00e51c554648e5094e215164a19f8"
registries = "General"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "7.8.1"

[[deps.OrdinaryDiffEqBDF]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqSDIRK", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "06dddbb32c023b4fc42f620bc7ff98f72d3d1f5b"
registries = "General"
uuid = "6ad6398a-0878-4a85-9266-38940aa047c8"
version = "2.4.14"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "BinaryHeaps", "CommonSolve", "ConstructionBase", "DiffEqBase", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "FastPower", "FindFirstFunctions", "FunctionWrappers", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "PrecompileTools", "Preferences", "Printf", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLLogging", "SciMLOperators", "SciMLStructures", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "5765dd5605c88f1db2b3630df6f291751a96cd16"
registries = "General"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "4.19.1"

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
git-tree-sha1 = "1ce6ba1f52c94059089f835f173fddda90233959"
registries = "General"
uuid = "50262376-6c5a-4cf5-baba-aaf4f84d72d7"
version = "2.6.2"

[[deps.OrdinaryDiffEqDifferentiation]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqCore", "SciMLBase", "SciMLOperators", "SparseMatrixColorings", "StaticArraysCore"]
git-tree-sha1 = "7bab95310a175e7b24b419cdc49d234c35a13a22"
registries = "General"
uuid = "4302a76b-040a-498a-8c04-15b101fed76b"
version = "3.12.5"
weakdeps = ["SparseArrays"]

    [deps.OrdinaryDiffEqDifferentiation.extensions]
    OrdinaryDiffEqDifferentiationSparseArraysExt = "SparseArrays"

[[deps.OrdinaryDiffEqNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConstructionBase", "DiffEqBase", "FastBroadcast", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MuladdMacro", "NonlinearSolve", "NonlinearSolveBase", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "PreallocationTools", "RecursiveArrayTools", "SciMLBase", "SciMLOperators", "SciMLPublic", "SimpleNonlinearSolve", "SparseArrays", "StaticArraysCore"]
git-tree-sha1 = "19dcb7f30de00a97d4c4e39094f398eeca31080a"
registries = "General"
uuid = "127b3ac7-2247-4354-8eb6-78cf4e7c58e8"
version = "2.9.10"

[[deps.OrdinaryDiffEqRosenbrock]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqRosenbrockTableaus", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "4d6354525db068b9d1578389274b7e6cabc663b4"
registries = "General"
uuid = "43230ef6-c299-4910-a778-202eb28ce4ce"
version = "2.7.5"

[[deps.OrdinaryDiffEqRosenbrockTableaus]]
git-tree-sha1 = "bab0d4ccee6fcb89a786b89207774314523e2e5d"
registries = "General"
uuid = "b4bd8bb3-f80f-41d2-9b21-73a655b304b9"
version = "2.4.2"

[[deps.OrdinaryDiffEqSDIRK]]
deps = ["ADTypes", "CommonSolve", "ConstructionBase", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "9e2bd5473cf2fb7c813f7c3f93a1f5a61b5a699a"
registries = "General"
uuid = "2d112036-d095-4a1e-ab9a-08536f3ecdbf"
version = "2.9.8"

[[deps.OrdinaryDiffEqTsit5]]
deps = ["CommonSolve", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "371de3eb5985febd9d7babb0937520653f2b15ca"
registries = "General"
uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a"
version = "2.1.5"

[[deps.OrdinaryDiffEqVerner]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "5a66a18325e678f260e849a62840cb4903ece511"
registries = "General"
uuid = "79d7bb75-1356-48c1-b8c0-6832512096c2"
version = "2.4.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.47.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "123266c25174ef6c8d4718920abc206452cf8de6"
registries = "General"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.41"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools"]
git-tree-sha1 = "e5ae52d1f4f19a9591d615cb0f1ef8d415b40933"
registries = "General"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "3.0.6"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
registries = "General"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "9869376035f47756b8b8e711fd0d06a06d51da6f"
registries = "General"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.84"

[[deps.PoissonRandom]]
deps = ["LogExpFunctions", "PrecompileTools", "Random"]
git-tree-sha1 = "7b7becc6ada13786479af8c5d08e7f7fe7d68301"
registries = "General"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.14"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "Setfield", "SparseArrays"]
git-tree-sha1 = "d86ea3a5364b306b8547e02741830e69c2d04448"
registries = "General"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.6"

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
registries = "General"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.Poppler_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "LibCURL_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "libpng_jll"]
git-tree-sha1 = "7dbfb7f61c3aa5def7b7dad3fa344c1c2858a83b"
registries = "General"
uuid = "9c32591e-4766-534b-9725-b71a8799265b"
version = "24.6.0+0"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
registries = "General"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "cd00e28071a75c98664663f87c2ae447d25c0503"
registries = "General"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "1.7.1"

    [deps.PreallocationTools.extensions]
    PreallocationToolsEnzymeCoreExt = "EnzymeCore"
    PreallocationToolsForwardDiffExt = "ForwardDiff"
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
registries = "General"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "aef1987d90dcc2e2cbd919ad10b2d25bcce4f418"
registries = "General"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.6.1"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "StyledStrings", "Tables"]
git-tree-sha1 = "eeac4c9337f876c3c32b8da732141953567f57c9"
registries = "General"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.5.1"

    [deps.PrettyTables.extensions]
    PrettyTablesExcelExt = "XLSX"
    PrettyTablesTypstryExt = "Typstry"
    PrettyTablesWriteDocxExt = "WriteDocx"

    [deps.PrettyTables.weakdeps]
    Typstry = "f0ed7684-a786-439e-b1e3-3b82803b501e"
    WriteDocx = "d049ceea-54ee-41d7-a26f-ba29db3b6599"
    XLSX = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
registries = "General"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Printf"]
git-tree-sha1 = "ce5ba3edfdeba107a721748d53b5f5887f415853"
registries = "General"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.13.0"

    [deps.ProgressMeter.extensions]
    ProgressMeterDistributedExt = "Distributed"
    ProgressMeterIJuliaExt = "IJulia"

    [deps.ProgressMeter.weakdeps]
    Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
registries = "General"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.PureKLU]]
deps = ["LinearAlgebra", "PrecompileTools", "SparseArrays"]
git-tree-sha1 = "ad0e84647fe7b3c2ae02aaeace1522139aaf6a27"
registries = "General"
uuid = "0c0d3e7f-3a8b-4f7e-b6f1-9a4d2e7c1f01"
version = "1.6.0"
weakdeps = ["ForwardDiff"]

    [deps.PureKLU.extensions]
    PureKLUForwardDiffExt = "ForwardDiff"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "5e8e8b0ab68215d7a2b14b9921a946fee794749e"
registries = "General"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.3"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["Base64", "Dates", "FileWatching", "InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
registries = "General"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.ReactantCore]]
deps = ["ExpressionExplorer", "MacroTools"]
git-tree-sha1 = "1fe4e616752a0c4818ae88a504fc59b0fb016659"
registries = "General"
uuid = "a3311ec8-5e00-46d5-b541-4f83e724a433"
version = "0.1.26"

[[deps.ReadOnlyArrays]]
git-tree-sha1 = "e6f7ddf48cf141cb312b078ca21cb2d29d0dc11d"
registries = "General"
uuid = "988b38a3-91fc-5605-94a2-ee2116b3bd83"
version = "0.2.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "0c0dcffc464cda2ff33fc8ffc1dd1bdd53cd98be"
registries = "General"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.4.0"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "LinearAlgebra", "PrecompileTools", "RecipesBase", "SciMLPublic", "SciMLStructures", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "47fe0f2a25b5e98412c65227f8a649d26bb7f1b9"
registries = "General"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "4.5.6"

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
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
registries = "General"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.ResettableStacks]]
deps = ["PrecompileTools", "StaticArrays"]
git-tree-sha1 = "ed16b7ff60555aa33a8013b9a165c404da3685b3"
registries = "General"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.4.0"

[[deps.RespecializeParams]]
deps = ["FunctionWrappersWrappers", "PrecompileTools"]
git-tree-sha1 = "140c074bb63518a0e1ec2c8aba2a8a5719bcf21c"
registries = "General"
uuid = "9fe22ead-9e00-4db2-8b46-706a60d40f5e"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
registries = "General"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d40b2fe70437b01397d2a4d5b020008da4e7019"
registries = "General"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.2+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "3eff988b9bd09543783e2e051b0a1eef11f65c2d"
registries = "General"
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
registries = "General"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "PrecompileTools", "SHA", "Serialization"]
git-tree-sha1 = "9bb075002b40849582905394bb3031fb77aa102d"
registries = "General"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.27"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FindFirstFunctions", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "LoggingExtras", "Markdown", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "Random", "RecipesBase", "RecursiveArrayTools", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLPublic", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "85f85c3c5b949e7d691047e758e3fdb77c21310c"
registries = "General"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "3.58.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseDifferentiationInterfaceExt = "DifferentiationInterface"
    SciMLBaseDistributionsExt = "Distributions"
    SciMLBaseEnzymeExt = "Enzyme"
    SciMLBaseForwardDiffExt = "ForwardDiff"
    SciMLBaseFunctionPropertiesExt = "FunctionProperties"
    SciMLBaseMakieExt = ["Makie", "GeometryBasics"]
    SciMLBaseMeasurementsExt = "Measurements"
    SciMLBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    SciMLBaseMooncakeExt = "Mooncake"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseReactantCoreExt = "ReactantCore"
    SciMLBaseReactantExt = "Reactant"
    SciMLBaseReverseDiffExt = ["ReverseDiff", "ForwardDiff"]
    SciMLBaseStaticArraysExt = "StaticArrays"
    SciMLBaseTrackerExt = "Tracker"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore", "ZygoteRules", "FillArrays"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DifferentiationInterface = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FillArrays = "1a297f60-69ca-5386-bcde-b61e274b549b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    FunctionProperties = "f62d2435-5019-4c03-9749-2d4c77af0cbc"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Reactant = "3c362404-f566-11ee-1572-e11a4b42c853"
    ReactantCore = "a3311ec8-5e00-46d5-b541-4f83e724a433"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"
    ZygoteRules = "700de1a5-db45-46bc-99cf-38207098b444"

[[deps.SciMLJacobianOperators]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DifferentiationInterface", "FastClosures", "LinearAlgebra", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "781d5413c2991e581da86979b0d9f86f9a034a1c"
registries = "General"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.20"

[[deps.SciMLLogging]]
deps = ["Logging", "LoggingExtras", "PrecompileTools", "Preferences"]
git-tree-sha1 = "3aafee7edc2ec7e71f50dcc6816cb031fac390a0"
registries = "General"
uuid = "a6db7da4-7206-11f0-1eab-35f2a5dbe1d1"
version = "2.1.0"

    [deps.SciMLLogging.extensions]
    SciMLLoggingTracyExt = "Tracy"

    [deps.SciMLLogging.weakdeps]
    Tracy = "e689c965-62c8-4b79-b2c5-8359227902fd"

[[deps.SciMLOperators]]
deps = ["Accessors", "Adapt", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "3feb8a7bce88bb3ef051a4eda8106004ad040b6f"
registries = "General"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.30.4"

    [deps.SciMLOperators.extensions]
    SciMLOperatorsLoopVectorizationExt = "LoopVectorization"
    SciMLOperatorsSparseArraysExt = "SparseArrays"

    [deps.SciMLOperators.weakdeps]
    LoopVectorization = "bdcacae8-1622-11e9-2a5c-532679323890"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SciMLPublic]]
deps = ["PrecompileTools"]
git-tree-sha1 = "74685afb51732a464fbce79a72708f7c4203ceb7"
registries = "General"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.3.0"

[[deps.SciMLStructures]]
deps = ["ArrayInterface", "PrecompileTools"]
git-tree-sha1 = "4b8fc0ffe22b4ecba7a02f1fec0b3757db014d0c"
registries = "General"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.10.6"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "4910c2c2f334d74735b7a8fd02f837ead352576a"
registries = "General"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.11"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "89880a2cca5cabc68c9023f58e544beda1d2fc1f"
registries = "General"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.2.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLLogging", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "8d9ec24ec7bd126abf97b44ba5bdd6d280572095"
registries = "General"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "2.14.6"

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
registries = "General"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

[[deps.Sobol]]
deps = ["DelimitedFiles", "Random"]
git-tree-sha1 = "5a74ac22a9daef23705f010f72c81d6925b19df8"
registries = "General"
uuid = "ed01d8cd-4d21-5b2a-85b4-cc3bdc58bad4"
version = "1.5.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
registries = "General"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.13.0"

[[deps.SparseColumnPivotedQR]]
deps = ["LinearAlgebra", "PrecompileTools", "SparseArrays"]
git-tree-sha1 = "8cf1a59c78dd6f900feb366e5a574f99ed278ef1"
registries = "General"
uuid = "a57abbd0-fea5-4d57-96be-5e525945e8e4"
version = "2.1.8"
weakdeps = ["AMD"]

    [deps.SparseColumnPivotedQR.extensions]
    SparseColumnPivotedQRAMDExt = "AMD"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "Random", "SparseArrays"]
git-tree-sha1 = "0730f14b07c68fd3b3a6fae89a829a5124a7ae7e"
registries = "General"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.29"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsAMDGPUExt = "AMDGPU"
    SparseMatrixColoringsCUDAExt = ["CUDA", "cuSPARSE"]
    SparseMatrixColoringsCliqueTreesExt = "CliqueTrees"
    SparseMatrixColoringsColorsExt = "Colors"
    SparseMatrixColoringsGPUArraysExt = "GPUArrays"
    SparseMatrixColoringsJuMPExt = ["JuMP", "MathOptInterface"]

    [deps.SparseMatrixColorings.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
    GPUArrays = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
    JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
    cuSPARSE = "b26da814-b3bc-49ef-b0ee-c816305aa060"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "88ef0d98204515bc53dd4aaea5f7788882d1b66c"
registries = "General"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.9.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StarAlgebras]]
deps = ["LinearAlgebra", "MutableArithmetics", "SparseArrays"]
git-tree-sha1 = "95c5ea2eef117ba1310ee14aa1737d3ea2d262ae"
registries = "General"
uuid = "0c0c59c1-dc5f-42e9-9a8b-b5dc384a6cd1"
version = "0.3.1"

    [deps.StarAlgebras.extensions]
    StarAlgebrasJuMPExt = "JuMP"
    StarAlgebrasMathOptInterfaceExt = "MathOptInterface"

    [deps.StarAlgebras.weakdeps]
    JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.StateSpaceSets]]
deps = ["Distances", "LinearAlgebra", "Neighborhood", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "b6cbad6435392da931d2d99fa54b96b09400d2eb"
registries = "General"
uuid = "40b095a5-5852-4c12-98c7-d43bf788e795"
version = "2.5.5"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "8a8a52f887f38d44014e24ee46b8a7311c625f6d"
registries = "General"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.24"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "bf7e2d0d9ed9ac2a652fd598c0f784910cc759f2"
registries = "General"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.5"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e2b53ce13a53367e96601081e33d34746b571bad"
registries = "General"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.5"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
registries = "General"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "adb9da019510162e67a4493fc235c23203d8b09e"
registries = "General"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.13"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91a5737baed20ee31f3faea0e51f57461f6a689e"
registries = "General"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "2.2.1"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StochasticDiffEqCore]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "DiffEqBase", "DiffEqNoiseProcess", "FastPower", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LinearAlgebra", "Logging", "MuladdMacro", "OrdinaryDiffEqCore", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLLogging", "SciMLOperators", "SimpleNonlinearSolve", "SparseArrays", "StaticArrays", "StochasticDiffEqLevyArea"]
git-tree-sha1 = "0c59edcc1d619c58a5fdba92b14d47d44590d713"
registries = "General"
uuid = "19c5a474-6cd1-4a5f-be79-46dc34e54d7f"
version = "2.2.5"

[[deps.StochasticDiffEqHighOrder]]
deps = ["DiffEqBase", "DiffEqNoiseProcess", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "StochasticDiffEqCore"]
git-tree-sha1 = "6e416fdb64d6c26538f39cd3259173a31318d6ad"
registries = "General"
uuid = "0520c28c-50fd-4d16-9c96-902fc80b3bab"
version = "2.2.1"

[[deps.StochasticDiffEqLevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "0ebc2c0670353f529e6b407a8a737692a622e5fd"
registries = "General"
uuid = "90dbc90e-856a-4131-af2c-e8c5aa0f35b5"
version = "2.1.1"

[[deps.StringManipulation]]
deps = ["PrecompileTools", "StyledStrings"]
git-tree-sha1 = "9e3a02d73a8f0a9be04e677dbef758aab40f8bea"
registries = "General"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.6.1"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.10.1+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "PrecompileTools", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "42dd4c193a1bacea5b6adf78b6224de907b06351"
registries = "General"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.56"
weakdeps = ["PrettyTables"]

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

[[deps.SymbolicLimits]]
deps = ["PrecompileTools", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "cdf3f5760ff2ece5abe16f26a92a3f7da0eb8276"
registries = "General"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "1.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "EnumX", "ExproniconLite", "Graphs", "LinearAlgebra", "MacroTools", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "OrderedCollections", "PrecompileTools", "ReadOnlyArrays", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "WeakCacheSets"]
git-tree-sha1 = "61a82af727479af3d96988079fe3e05cdd57d695"
registries = "General"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "4.55.0"

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
deps = ["ADTypes", "AbstractPlutoDingetjes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "IntervalSets", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "Preferences", "Primes", "RecipesBase", "RuntimeGeneratedFunctions", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "49a65a35d6fb99d152ad7266f8acd278da7e8d83"
registries = "General"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "7.48.0"

    [deps.Symbolics.extensions]
    SymbolicsD3TreesExt = "D3Trees"
    SymbolicsDistributionsExt = "Distributions"
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsHypergeometricFunctionsExt = "HypergeometricFunctions"
    SymbolicsLatexifyExt = ["Latexify", "LaTeXStrings"]
    SymbolicsNemoExt = "Nemo"
    SymbolicsStaticArraysExt = "StaticArrays"
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
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
registries = "General"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "a94d9bdda1b7bed0046cea645639ab3f62196fac"
registries = "General"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.14.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
registries = "General"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

[[deps.TermInterface]]
deps = ["PrecompileTools"]
git-tree-sha1 = "d6984451701c6462ce685d0e2d4f2c47cd966773"
registries = "General"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TikzPictures]]
deps = ["LaTeXStrings", "Poppler_jll", "tectonic_jll"]
git-tree-sha1 = "75e5927a4d4707d2f4bd8332320f3deef6f16d7c"
registries = "General"
uuid = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
version = "3.5.2"

[[deps.TimerOutputs]]
deps = ["ExprTools", "PrettyTables", "Printf", "Tables"]
git-tree-sha1 = "2244393f04554e24b54896b430b2ab3850a2232d"
registries = "General"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "1.2.2"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
registries = "General"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
registries = "General"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
registries = "General"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
deps = ["Serialization"]
git-tree-sha1 = "eb3b85ca904a1d1647b15acc091b113b8874b648"
registries = "General"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.7.3"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "1f0f9f401753701a7e4113b5056ca38d33875b55"
registries = "General"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.29.0"

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
registries = "General"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.7.2"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "fbe44a0ade62ae5ed0240ad314dfdd5482b90b40"
registries = "General"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.2.2"

[[deps.WeakCacheSets]]
git-tree-sha1 = "386050ae4353310d8ff9c228f83b1affca2f7f38"
registries = "General"
uuid = "d30d5f5c-d141-4870-aa07-aabb0f5fe7d5"
version = "0.1.0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
registries = "General"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "80d3930c6347cfce7ccf96bd3bafdf079d9c0390"
registries = "General"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.9+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e52eca002a11c30a858185efdfb15311e1c7a6bf"
registries = "General"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "ce91b54a5f90c4f27df6950706234b34e49cf79e"
registries = "General"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.14+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
registries = "General"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
registries = "General"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
registries = "General"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a54cf54677c9f6eac0d97e9d757c233b56bc21b9"
registries = "General"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.13+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
registries = "General"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
registries = "General"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.2+0"

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "1fcd6cf05d23adab0b5ffc496d020537215751dc"
registries = "General"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.59+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.70.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ed2d8d4ad96c03d3abe3382cdbbdea9489b63636"
registries = "General"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2023.1.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.8.3+0"

[[deps.tectonic_jll]]
deps = ["Artifacts", "Fontconfig_jll", "FreeType2_jll", "Graphite2_jll", "HarfBuzz_ICU_jll", "HarfBuzz_jll", "ICU_jll", "JLLWrappers", "Libdl", "OpenSSL_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "b62c5dcf5d80a82e40d58b908b8eca27a54f215b"
registries = "General"
uuid = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"
version = "0.15.0+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

# ╔═╡ Cell order:
# ╠═fed88caa-1520-41f7-adb3-785e5c9529c6
# ╟─734b3b08-061e-4f93-8574-468d824815da
# ╟─313c5c60-1c11-4fa3-a71f-16df0d2d0ebc
# ╟─b842e98e-34e2-40f2-84b6-c180815c2df3
# ╟─99d0fdc9-e368-462c-a357-86f07624a52e
# ╟─8febbbdb-b31c-4232-8a8b-15ef6a87557b
# ╟─71b94af7-76ee-4228-987c-2f22e0951552
# ╟─d12b4a5f-a6ad-4c91-8795-fe8f93f5c93d
# ╟─c7eee8bb-b03d-4679-a51e-36fc629294e1
# ╟─8da87629-09ff-4677-b5e8-2a3ff5777a8e
# ╟─700015f6-aac0-40c5-a455-8d98c1227049
# ╠═207305eb-4496-4518-a5bc-be85173314a5
# ╟─a6fdec29-3438-4ebc-af81-7a3baf0175ae
# ╟─ee7baf70-c7af-4791-8277-fadc0961e0e7
# ╟─72416f96-8a0a-4c13-97d1-2990e49c4cb0
# ╠═0baae21f-71b5-42fd-be33-31de085928d3
# ╟─a7c5650c-9bbc-44a1-8d51-d473864faae7
# ╠═49d1a7f7-2bf2-4472-94df-6247b9237ddd
# ╟─80bae1ac-9f61-4512-84c9-c562ad7a61ba
# ╠═568d9fe3-6716-4a9a-ba1d-b9e6fd039150
# ╠═dc45a37a-b939-499f-b9d2-4140aee92282
# ╠═eb211ddb-1ba6-4ef0-8712-4aa7c4320c26
# ╟─0a6364c3-851e-4f0b-b140-f4f0e633fb36
# ╠═f9bcdb15-ce25-49a4-9837-53d042065151
# ╟─14ae7f11-1065-46aa-a6ed-eea400c0d2ec
# ╠═c8efcb3f-cc6f-4c45-a0f3-55cdb73d5195
# ╟─d59793b2-3414-4051-ac2b-33b6cb8063e8
# ╟─225887e5-e214-48a0-8a7c-7814b4e66c83
# ╟─a33249d1-4b85-494c-82c9-0b5a164e32d1
# ╠═82d801d9-49f4-4d74-93c0-388b09dc76f1
# ╠═c869e896-cc87-4b12-bba2-84b0cf0964c6
# ╠═2c2f97cc-127a-446e-a6cd-e9651df833f0
# ╟─9a756404-f35b-40b9-bc88-4b7e3a7df0b6
# ╠═805ef471-1bd2-4ce2-bf50-eeb9c6b10e4d
# ╟─6eef3827-1af6-431d-a163-67bb595c337d
# ╠═e14648ff-69f6-47c2-924c-c5ac69e200ed
# ╟─79fa6ccd-e137-4178-9376-c1222dd0a8f7
# ╠═e7f0842c-954a-40e1-a112-ec03631e5270
# ╟─aaba63dd-8fc3-405d-ad19-8e8368e70019
# ╠═03754250-3831-41db-ac11-fe6f5a08f0c1
# ╟─6ddf393c-9f37-43b2-b401-9ab0e7a9e88c
# ╠═5ce798b5-4366-4cc1-8d2d-881a827d6dc9
# ╟─7df63fdf-ef41-44f2-8a55-e5c2c849029c
# ╠═64c2a3ee-e7a0-4e03-b9a8-86239e1ca81e
# ╠═0895d464-a029-410d-9e7d-89cfac2d1615
# ╠═07cd9aad-029f-42c6-abe8-ab4a9a2a910c
# ╠═cf02b69e-f6e5-4312-9a85-d43b61fac435
# ╟─6cab6cb7-a432-40b6-9390-ad0083fe486d
# ╠═38e91d0d-1020-44a8-becc-8542fd600104
# ╠═9bf08113-14a4-468a-8766-a6a9b1ed0577
# ╟─b5ede1bb-8361-4380-b943-c9780d4a6253
# ╠═027dbe5e-4993-41c2-983d-289101c118fb
# ╟─4b81f302-9637-4161-b745-8ac39b9e31d3
# ╠═866df30a-7a45-4821-a4ef-fa886040d9c0
# ╟─9af6ce74-5678-4b48-8758-37b0d5a6f0e4
# ╠═a9c6a292-086e-4aa0-9856-78d3ab3fbe35
# ╟─5890b699-b7de-47a3-bee7-1e7dd7663fbe
# ╠═ef79392b-395c-497a-9c0e-dc2cd468f6e1
# ╟─8c9ab125-2acb-4732-a9bf-7838e819e4f7
# ╠═2a39d6f8-da49-4976-9aa7-889391e55a5d
# ╠═3a2fce00-7fee-4c04-9ed0-8e7cdc5bcc0a
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
# ╟─ab1b7695-3cd6-4e67-9cfd-fbaf2f4d1d15
# ╟─6475b25b-8711-44cd-bc3b-e3d42681bb93
# ╟─14a6478e-939f-4e9c-b905-4cec39917d47
# ╠═d7cc8a66-220a-4e9b-b42e-a0e085ed3a0f
# ╠═3871738a-a126-4996-ad95-9674564a4931
# ╟─0fc36dec-9f45-4f49-a44f-7884ec33d7be
# ╟─24920a1b-355c-43e1-aff2-028e6f25584c
# ╠═ad67d733-4365-4b37-ab45-1b936c80d617
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
# ╟─477c8f59-97d1-405d-a1c8-64b7e0b9119f
# ╟─08e05c65-06a2-4560-92eb-b014dc7c3d70
# ╟─2235689d-9c83-4907-aa17-c2624fbeb68d
# ╟─df8a9449-851c-4546-97a7-7fd4a270a867
# ╠═e8c51ae9-ded4-450a-80f9-c484c2b991d1
# ╠═c96ea0ad-47c3-469e-a2f7-0743e681a6d1
# ╠═cd0e57c7-3fa0-4875-9a38-249b3a0bb134
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
