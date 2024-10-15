### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ fed88caa-1520-41f7-adb3-785e5c9529c6
using ChaosTools, DataFrames, DataFramesMeta, DelimitedFiles, DifferentialEquations, Interpolations, LinearAlgebra, PlutoUI, QuadGK, SpecialFunctions, Symbolics, TikzPictures, Trapz, Unitful, UnitfulAstro

# ╔═╡ 734b3b08-061e-4f93-8574-468d824815da
# ╠═╡ skip_as_script = true
#=╠═╡
TableOfContents(title="🌌 SF model", depth=4)
  ╠═╡ =#

# ╔═╡ 53ad27ee-c5be-4d29-9ed5-21b6b64de42b
# ╠═╡ skip_as_script = true
#=╠═╡
md"# Star formation model"
  ╠═╡ =#

# ╔═╡ a8f2765b-444f-4f3a-8a19-fbcfa6e79d2c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Motivation

The star formation rate (SFR) is a key characteristic of galaxies. In the context of the standard cosmological model ($\Lambda$CDM), the SFR is given by a combination of processes that take place over the course of a galaxy's lifetime, such as gas cooling, star formation, chemical enrichment, and feedback from supernovae and galactic nuclei. These processes are influenced by factors like mergers, interactions, and mass accretion, which affect the amount and properties of the star forming gas. The density of a gas cloud is believed to be the most important factor in determining its star formation rate, although the details of this process are not yet fully understood. Observationally, the total gas density is found to be correlated to the star formation rate ([Kennicutt1998](https://doi.org/10.1086/305588)), and this correlation is even stronger with molecular gas ([Wong2002](https://doi.org/10.1086/339287), [Bigiel2008](https://doi.org/10.1088/0004-6256/136/6/2846)). The underlying reason is the intrinsic relation between the molecular gas mass and the SFR, which can be found at resolved scales ([Baker2021](https://doi.org/10.1093/mnras/stab3672)) and at integrated (i.e. galaxy-wide) scales across redshifts ([Baker2022](https://doi.org/10.1093/mnras/stac3413)).

Because the formation of dark matter halos and galaxies is highly non-linear, numerical simulations have become the preferred tool to investigate how galaxies form and evolve from early times up to the present. This type of simulation naturally includes mergers/interactions and continuous gas accretion. However, there are still significant uncertainties in the modeling of the baryonic component, since the physical processes that affect baryons – such as star formation, feedback, and chemical enrichment – take place at scales that are too small to be resolved directly. As a result, these processes are introduced using sub-grid physics involving several adjustable parameters that are not always independent of one another or constrained by observations. This can lead to inconsistencies in the predictions of different models ([Scannapieco2012](https://doi.org/10.1111/j.1365-2966.2012.20993.x), [Zhu2016](https://doi.org/10.3847/0004-637X/831/1/52), [Naab2017](https://doi.org/10.1146/annurev-astro-081913-040019)).

Because of its importance in galaxy formation, it is critical for simulations to accurately describe the star formation process at the scales that can be resolved and the associated feedback effects.
"""
  ╠═╡ =#

# ╔═╡ c23a574d-4150-47a8-bcaf-45ae5a8c7c66
# ╠═╡ skip_as_script = true
#=╠═╡
md"""

## Previous work

As the precision of $\mathrm{H_2}$ measurements increases, so too must the sophistication of simulations that model its formation and evolution. Many simulations use an equilibrium model to define the $\mathrm{H_2}$ content, which assumes that the chemistry in each volume element is in a state of equilibrium given only by local variables ([Krumholz2008](https://doi.org/10.1086/592490), [McKee2010](https://doi.org/10.1088/0004-637X/709/1/308), [Krumholz2011](https://doi.org/10.1088/0004-637X/729/1/36), and [Krumholz2013](https://doi.org/10.1093/mnras/stt1780)). However, this assumption is not always valid; as the formation and destruction of $\mathrm{H_2}$ can be influenced by a variety of factors, including temperature, density, and the presence of ultraviolet (UV) radiation.

In [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) it is described, for the first time, a non-equilibrium chemical network for molecular hydrogen. This network uses rate equations to track the formation and destruction of $\mathrm{H_2}$ in each volume element, considering a variety of factors that can influence its abundance. In the last decade, several non-equilibrium models have been developed and implemented in hydrodynamical simulations. In general, they use radiative transfer to model the radiation field, which, in conjunction with a chemical network (a set of ODEs), evolves the abundance of molecular hydrogen and other species of hydrogen and helium. A summary table of previous work can be seen below

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

As an alternative to the computationally expensive radiative transfer process, semi-analytical models (SAMs) have been developed for the multiphase structure of the interstellar medium (MP ISM). Broadly, there are two ways to model the MP ISM, one considers the physical properties of the gas (hot and cold phases), and the other its chemical composition (hydrogen phases).

The former was pioneer by [Field1969](https://doi.org/10.1086/180324) (see [Cowie1977](https://doi.org/10.1086/154911), [McKee1977a](https://doi.org/10.1086/155350), and for a review [Cox2005](https://doi.org/10.1146/annurev.astro.43.072103.150615)), within the context of pure SAMs. The model developed by [McKee1977b](https://doi.org/10.1086/155667) was first incorporated into numerical simulation of galaxy formation by [Yepes1997](https://doi.org/10.1093/mnras/284.1.235) (Eulerian simulations) and [Hultman1999](https://ui.adsabs.harvard.edu/abs/1999A%26A...347..769H) (Lagrangian simulations). These works were later extended by [Springel2003](https://doi.org/10.1046/j.1365-8711.2003.06206.x), adding galactic winds driven by star formation as feedback.

[Monaco2004](https://doi.org/10.1111/j.1365-2966.2004.07916.x) developed a SAM in a similar vein to [Springel2003](https://doi.org/10.1046/j.1365-8711.2003.06206.x), providing the theoretical foundation for MUPPI (MUlti-Phase Particle Integrator) ([Murante2010](https://doi.org/10.1111/j.1365-2966.2010.16567.x)), a sub-resolution MP ISM model that adds stellar feedback to $\texttt{GADGET2}$ ([Springel2005](https://doi.org/10.1111/j.1365-2966.2005.09655.x)). MUPPI separates a gas particle in a hot and cold phase if a set of conditions for its density and temperature are met. It then evolves those components, plus a stellar phase and an energy term (energy of the hot gas), using a set of four ODEs.

Based on [Ferrini1992](https://doi.org/10.1086/171066) and later work, [Mollá2015](https://doi.org/10.1111/j.1365-2966.2005.08782.x) developed a SAM to follow the metallicity in galaxies. These chemical evolution models (CEMs) were subsequently improved and extended in [Mollá2015](https://doi.org/10.1093/mnras/stv1102), [Mollá2016](https://doi.org/10.1093/mnras/stw1723), [Mollá2017](https://doi.org/10.1093/mnras/stx419) and [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635). The latter tracks five components: the molecular, atomic, and ionized phases of hydrogen, dust, and the stars.

We will use a sub-resolution SAM, following closely the one developed by [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635) but implemented within the hydrodynamical code $\texttt{AREPO}$ ([Springel2010](https://doi.org/10.1111/j.1365-2966.2009.15715.x)), like the way MUPPI is integrated with $\texttt{GADGET3}$.
"""
  ╠═╡ =#

# ╔═╡ 417aad38-6928-4b13-9286-3f11efcffb99
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Arepo

We use the N-body, magnetohydrodynamics (MHD) code $\texttt{AREPO}$ ([Springel2010](https://doi.org/10.1111/j.1365-2966.2009.15715.x)) to perform the simulations, which we briefly describe here.

The code $\texttt{AREPO}$ is a moving-mesh code that tracks collisionless dynamics in a cosmological setting. Gravitational forces are computed using a conventional TreePM method ([Springel2005](https://doi.org/10.1111/j.1365-2966.2005.09655.x)), which employs a fast Fourier Transform technique for long-range forces and a hierarchical oct-tree algorithm ([Barnes1986](https://doi.org/10.1038/324446a0)) for short-range forces, in conjunction with adaptive time-stepping.

The code $\texttt{AREPO}$ employs a dynamic unstructured mesh, constructed from a Voronoi tessellation of a set of mesh-generating points, to track the MHD. This allows for a finite-volume discretization of the MHD equations. The MHD equations are solved using a second-order Runge–Kutta integration scheme with high-accuracy least-squares spatial gradient estimators of primitive variables ([Pakmor2015](https://doi.org/10.1093/mnras/stv2380)), which enhance the estimators in the original version of $\texttt{AREPO}$.

A distinctive feature of $\texttt{AREPO}$ is the ability to transform the mesh at any time-step through a mesh reconstruction, a capability not found in standard grid-based methods. The mesh construction ensures that each cell contains a specific target mass (within a certain tolerance), meaning that areas of high density are resolved with more cells than areas of low density. Additionally, the mesh generating points can move with the fluid flow, so each cell of the newly constructed mesh moves approximately with the fluid. This allows $\texttt{AREPO}$ to overcome the Galilean non-invariance problem that standard Eulerian mesh codes have and significantly reduce the advection errors that appear in complex supersonic flows. The quasi-Lagrangian nature of the method makes it comparable to other Lagrangian methods such as SPH, although it eliminates several limitations of the SPH method, such as the absence of artificial viscosity, and provides higher accuracy in the hydrodynamics of underdense regions.
"""
  ╠═╡ =#

# ╔═╡ b842e98e-34e2-40f2-84b6-c180815c2df3
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Phases and notation

We will model the interstellar medium as a multi-phase structure made up of four components. Three are the different states of hydrogen: ionized gas with a temperature of $\sim \! 10^4 \, \mathrm{K}$, atomic gas with a temperature of $\sim \! 100 \, \mathrm{K}$, and molecular gas with a temperature of $\sim \! 10 \, \mathrm{K}$. The final component is the stars.

For every reaction that transfers mass between the phases, we will only use the dominant channel. So, even though the gas is mostly made up of hydrogen and helium, only the hydrogen reactions are incorporated into the model. The processes involved are the photoionization of atoms, the recombination of electrons with ions, the conversion of atomic hydrogen into molecular hydrogen, and the destruction of the latter through photodissociation caused by UV light. In addition, we consider the formation of ionized gas by supernovas and the influence of the molecular gas on the SFR.

We characterized each phase by its mass fraction with respect to the total mass of the cell,

|||
|:-------------:|:--------------------:|
| Ionized gas   | $f_i(t) := M_i / M_\mathrm{cell}$ |
| Atomic gas    | $f_a(t) := M_a / M_\mathrm{cell}$ |
| Molecular gas | $f_m(t) := M_m / M_\mathrm{cell}$ |
| Stars         | $f_s(t) := M_s / M_\mathrm{cell}$ |

where $M_i$, $M_a$, $M_m$, $M_s$ are the corresponding masses and

$\begin{equation}
    M_\mathrm{cell} := M_i(t) + M_a(t) + M_m(t) + M_s(t) \, ,
\end{equation}$

is the total cell mass.

Now, we will compute the relation between the number density of the $j$ component ($n_j$) and the dimensionless fractions defined above ($f_j := M_j / M_\mathrm{cell}$).

$\begin{equation}
    n_j := \frac{N_j}{V_j} = \frac{M_j}{m_j} \, \frac{1}{V_j} = \frac{M_j}{M_\mathrm{cell}} \,  \frac{M_\mathrm{cell}}{m_j \, V_j} = f_j \, \frac{M_\mathrm{cell}}{m_j \, x_j \, V_\mathrm{cell}} = f_j \, \frac{\rho_\mathrm{cell}}{m_j \, x_j} \, ,
\end{equation}$

where the quantities are

|||
|:--------:|:----------------------------------------------------------:|
| $V_\mathrm{cell}$    | Volume of the cell                                         |
| $M_\mathrm{cell}$    | Total mass of the cell                                     |
| $\rho_\mathrm{cell}$ | Mass density of the cell                                   |
| $N_j$    | Number of elements (e.g. atoms) of the $j$ component       |
| $V_j$    | Total volume of the $j$ component                          |
| $M_j$    | Total mass of the $j$ component                            |
| $m_j$    | Mass of a single element (e.g. atoms) of the $j$ component |
| $f_j$    | $(:= M_j / M_\mathrm{cell})$ Mass fraction of the $j$ component        |
| $x_j$    | $(:= V_j / V_\mathrm{cell})$ Volume fraction of the $j$ component      |

In the context of our model, $V_\mathrm{cell}$, $\rho_\mathrm{cell}$, $M_\mathrm{cell}$, $m_j$, $V_j$, and $x_j$ are constants.

In the same way, we can write a relation for the mass density of the $j$ component,

$\begin{equation}
    \rho_j := \frac{M_j}{V_j} = \frac{M_j}{M_\mathrm{cell}} \, \frac{M_\mathrm{cell}}{V_j} = f_j \,  \frac{M_\mathrm{cell}}{x_j \, V_\mathrm{cell}} = f_j \, \frac{\rho_\mathrm{cell}}{x_j} \, ,
\end{equation}$

So, using these relations we can write any differential equation for the quantities $M_j$, $\rho_j$, and $n_j$ as an equation for $f_j$.

In our model, we have done only two hypotheses until now. First, that the ISM is only made up of the four components already mentioned, and second, that $V_\mathrm{cell}$, $\rho_\mathrm{cell}$, $M_\mathrm{cell}$, $m_j$, $V_j$, and $x_j$ are constants.

For simplicity we will adopt $x_i = x_a = x_m = 1.0$, which is like assuming that the three hydrogen phases occupy the whole cell.
"""
  ╠═╡ =#

# ╔═╡ 99d0fdc9-e368-462c-a357-86f07624a52e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Physical relationships
"""
  ╠═╡ =#

# ╔═╡ 6d6cbebb-bb3b-437d-9835-0a79f36857f2
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Physical processes by name
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
		\draw[line, white, ->]
		(ion) edge [bend left, "\textcolor{d_pink}{recombination}"] (atom)
		(atom) edge [bend left, "\textcolor{d_orange}{condensation}"] (molecule)
		(molecule) edge [bend left,"\textcolor{d_green}{dissociation}"] (atom)
		(atom) edge [bend left,"\textcolor{d_blue}{ionization}"] (ion)
		(stars) edge [bend left, "\textcolor{d_yellow}{supernova}"] (ion)
		(molecule) edge [bend left, "\textcolor{red}{star formation}"] (stars);
	""",
	width="75em",
	preamble = """
		\\usepackage{chemformula}
		\\definecolor{d_pink}{HTML}{C721DD}
		\\definecolor{d_orange}{HTML}{D14A00}
		\\definecolor{d_green}{HTML}{008C00}
		\\definecolor{d_blue}{HTML}{007FB1}
		\\definecolor{d_yellow}{HTML}{D1AC00}
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
### Physical processes as mathematical expressions
"""
  ╠═╡ =#

# ╔═╡ 2acad39a-23b7-4004-bf0e-59e73b914f01
# ╠═╡ skip_as_script = true
#=╠═╡
TikzPictures.TikzPicture(
	L"""
		\node[box, white] (stars) at (180:2cm) {Stars};
		\node[box, white, text width=2em] (atom) at (0:2cm) {HI};
		\node[box, white, text width=2em] (molecule) at (270:2cm) {\ch{H2}};
		\node[box, white, text width=2em] (ion) at (90:2cm) {HII};
		\draw[line, white, ->]
		(ion) edge [bend left, "$\textcolor{d_pink}{\frac{f_i(t)}{\tau_\mathrm{rec}(t)}} - \textcolor{d_blue}{\eta_\text{ion} \, \psi(t)}$"] (atom)
		(atom) edge [bend left, "$\textcolor{d_orange}{\frac{f_a(t)}{\tau_\mathrm{cond}(t)}} - \textcolor{d_green}{\eta_\text{diss} \, \psi(t)}$"] (molecule)
		(stars) edge [bend left, "$\textcolor{d_yellow}{R \, \psi(t)}$"] (ion)
		(molecule) edge [bend left, "$\textcolor{red}{\psi(t)}$"] (stars);
	""",
	width="75em",
	preamble = """
		\\usepackage{chemformula}
		\\definecolor{d_pink}{HTML}{C721DD}
		\\definecolor{d_orange}{HTML}{D14A00}
		\\definecolor{d_green}{HTML}{008C00}
		\\definecolor{d_blue}{HTML}{007FB1}
		\\definecolor{d_yellow}{HTML}{D1AC00}
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

# ╔═╡ 82991902-223e-42bb-80f4-9d6260f8a040
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Equations"
  ╠═╡ =#

# ╔═╡ c5e21675-f120-4555-84be-d99b35592f2d
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Stars

We define $\psi(t)$ as the fractional star formation rate (SFR per unit of cell mass),

$\begin{equation}
    \left. \frac{\mathrm{d}}{\mathrm{d}t}f_s(t) \right|_{\text{SFR}} =: \psi(t) \, ,
\end{equation}$

### Ionized gas

The ionized component grows through the ionization of atomic gas and from the remnants of supernova explosions.

The former is produced mainly by the radiation of newborn stars, so it is given by

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_i(t)\right|_{\text{ion.}} = \eta_\text{ion} \, \psi(t) \, ,
\end{equation}$

where $\eta_\text{ion}$ is the ionized mass rate per unit of created stellar mass. All the physics of the ionization process are summarized in this parameter, which is assumed constant during the integration of the differential equations.

The latter, under the instantaneous recycling hypothesis, can be written as

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_i(t)\right|_{\text{recyc.}} = R \, \psi(t) \, ,
\end{equation}$

where $R$ is the mass fraction of a stellar population that is returned to the ISM, where we assumed that all the returned mass is in the form of ionized gas. Ignoring the metal enrichment is a good approximation that does not alter significantly the results. This parameter it is also assumed constant while integrating the differential equations.

### Atomic gas

The atomic component grows through the dissociation of hydrogen molecules and the recombination of the ionized gas with free electrons.

The former, as with the ionized gas, is given by

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_a(t)\right|_{\text{diss.}} = \eta_\text{diss} \, \psi(t) \, ,
\end{equation}$

where $\eta_\text{diss}$ is the disassociated mass rate per unit of created stellar mass, also assumed constant during the ODEs integration.

The latter will depend on the mass of ionized gas and the time scale of recombination ($\tau_\mathrm{rec}$), so it is given by

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_a(t)\right|_{\text{recon.}} = \frac{f_i(t)}{\tau_\mathrm{rec}(t)} \, .
\end{equation}$

### Molecular gas

The molecular component gains mass mainly by the condensation of hydrogen atoms on the surface of dust grains. This process depends on the mass of atomic gas and the characteristic time scale of condensation ($\tau_\mathrm{cond}$). We are putting all the dust physics into this time parameter,

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_m(t)\right|_{\text{cond.}} = \frac{f_a(t)}{\tau_\mathrm{cond}(t)} \, ,
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 8dcc45e1-192a-4007-9d6e-6e4149b563d7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
From all the above, we can write the system of four ODEs,
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
		\dv{}{t}f_i(t) &= - \textcolor{d_pink}{\frac{f_i(t)}{\tau_\mathrm{rec}(t)}} + \textcolor{d_blue}{\eta_\text{ion} \, \psi(t)} + \textcolor{d_yellow}{R \, \psi(t)} \, , \\
		\dv{}{t}f_a(t) &= \textcolor{d_pink}{\frac{f_i(t)}{\tau_\mathrm{rec}(t)}} - \textcolor{d_blue}{\eta_\text{ion} \, \psi(t)} - \textcolor{d_orange}{\frac{f_a(t)}{\tau_\mathrm{cond}(t)}} + \textcolor{d_green}{\eta_\text{diss} \, \psi(t)} \, , \\
		\dv{}{t}f_m(t) &= \textcolor{d_orange}{\frac{f_a(t)}{\tau_\mathrm{cond}(t)}} - \textcolor{d_green}{\eta_\text{diss} \, \psi(t)} - \textcolor{red}{\psi(t)} \, , \\
		\dv{}{t}f_s(t) &= \textcolor{red}{\psi(t)} - \textcolor{d_yellow}{R \, \psi(t)} \, ,
	\end{aligned}}$
	};
	""",
	width="75em",
	preamble = """
		\\usepackage{chemformula}
		\\usepackage{physics}
		\\setlength{\\jot}{10pt}
		\\definecolor{d_pink}{HTML}{C721DD}
		\\definecolor{d_orange}{HTML}{D14A00}
		\\definecolor{d_green}{HTML}{008C00}
		\\definecolor{d_blue}{HTML}{007FB1}
		\\definecolor{d_yellow}{HTML}{D1AC00}
	""",
)
  ╠═╡ =#

# ╔═╡ c7eee8bb-b03d-4679-a51e-36fc629294e1
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
And we can explicitly check mass conservation,

$\begin{equation}
    \frac{\mathrm{d}}{\mathrm{d}t}(f_i + f_a + f_m + f_s) = 0 \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 7112c2fc-0cd6-407f-afaf-56a49e02fa57
# ╠═╡ skip_as_script = true
#=╠═╡
md"## ODE coefficients"
  ╠═╡ =#

# ╔═╡ 6b135719-c120-4b69-aed3-28647fc601bf
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Stiffness coefficient"
  ╠═╡ =#

# ╔═╡ f036dffd-8aab-4337-9f1a-a1870da284df
# ╠═╡ skip_as_script = true
#=╠═╡
md"The model, for a typical set of initial conditions, is very stiff (stiffness ratio >> 1)."
  ╠═╡ =#

# ╔═╡ fd7e6a35-72a0-4133-b43c-95ee73da5fc4
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Lyapunov spectrum"
  ╠═╡ =#

# ╔═╡ 4d09ff8b-fd43-4b4a-ba9d-030e9d2b6cec
# ╠═╡ skip_as_script = true
#=╠═╡
md"The model, for a typical set of initial conditions, is not chaotic (maximum lyapunov exponent < 0). But, we note that with a higher density the system can turn chaotic (maximum lyapunov exponent > 0)."
  ╠═╡ =#

# ╔═╡ 700015f6-aac0-40c5-a455-8d98c1227049
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Units

We have the freedom to choose three independent units in this model. Time, mass, and length. The choice is reflected in the constants $C_\mathrm{star}$, $C_\mathrm{rec}$, and $C_\mathrm{cond}$.

Following the standard in astronomy and astrophysics, we will use $\mathrm{[T] = Myr}$, $\mathrm{[M] = mp}$ and $\mathrm{[L] = cm}$, where $\mathrm{mp}$ is the proton mass.
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
# Parameters

*  $\psi(t)$: Star formation rate.

*  $\tau_\mathrm{rec}$: Time scale of atomic gas formation, from ionized gas, generally called recombination time.

*  $\tau_\mathrm{cond}$: Time scale of molecular gas formation, from atomic gas, generally called condensation (or cloud formation) time.

*  $\eta_\mathrm{diss}$: Rate of molecular gas dissociation by stars per unit of created stellar mass.

*  $\eta_\mathrm{ion}$: Rate of atomic gas ionization by stars per unit of created stellar mass.

*  $R$: Mass of ionized gas produced per unit of created stellar mass.
"""
  ╠═╡ =#

# ╔═╡ ee7baf70-c7af-4791-8277-fadc0961e0e7
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Physical processes"
  ╠═╡ =#

# ╔═╡ 72416f96-8a0a-4c13-97d1-2990e49c4cb0
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Star formation

For the formation of stars we will use the notation and definitions commonly used in the field (for a review see [McKee2007](https://doi.org/10.1146/annurev.astro.45.051806.110602), [Krumholz2014](https://doi.org/10.1016/j.physrep.2014.02.001), and [Krumholz2019](https://doi.org/10.1146/annurev-astro-091918-104430)).

In particular, we will follow [Krumholz2005](https://doi.org/10.1086/431734), but taking into account the correlation between molecular hydrogen and star formation ([Bigiel2008](https://doi.org/10.1088/0004-6256/136/6/2846), [Bigiel2010](https://doi.org/10.1088/0004-6256/140/5/1194), [Wong2002](https://doi.org/10.1086/339287), [Robertson2008](https://doi.org/10.1086/587796), [Halle2013](https://doi.org/10.1051/0004-6361/201220952), [Thompson2013](https://doi.org/10.1088/0004-637x/780/2/145)).

So, our model is

$\begin{equation}
	\mathrm{SFR} := \frac{\mathrm{d}}{\mathrm{d}t} M_s = \frac{M_m}{\tau_\mathrm{star}}
\end{equation}$
where $M_m$ is the molecular mass and $\tau_\mathrm{star}$ is the characteristic timescale of star formation (defined by this very relation).

Given that we want equations for the dimensionless fraction, we will use

$\begin{equation}
	\psi := \frac{\mathrm{SFR}}{M_\mathrm{cell}} = \frac{M_m}{\tau_\mathrm{star}} \, \frac{1}{M_\mathrm{cell}} = \frac{f_m}{\tau_\mathrm{star}} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 099a7feb-019a-4e37-9332-592c9b1dc0ed
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Following [Krumholz2005](https://doi.org/10.1086/431734), we write the characteristic timescale of star formation as

$\begin{equation}
    \tau_\mathrm{star} =  \frac{t_\text{ff}}{\epsilon_\text{ff}} \, ,
\end{equation}$

where $\epsilon_\text{ff}$ is the star formation efficiency (in the literature is often used $\epsilon_\text{ff} \approx 0.01$ [Krumholz2007](https://doi.org/10.1086/509101), [Krumholz2014](https://doi.org/10.1016/j.physrep.2014.02.001)) and $t_\text{ff}$ is the free-fall time, which is the time for a pressure-free spherical cloud to collapse into a point due to its self-gravity.

The free-fall time can be written as

$\begin{equation}
    t_\text{ff} = \sqrt{\frac{3\pi}{32 \, G\, \rho_g}} \, ,
\end{equation}$

where $\rho_g$ is the density of the gas cloud.

There is a lot of uncertainty for the parameter $\epsilon_\text{ff}$ ([Lee2016](https://doi.org/10.3847/1538-4357/833/2/229) and [Utomo2018](https://doi.org/10.3847/2041-8213/aacf8f)). We will use $\epsilon_\text{ff} = 1.0$ because we already consider the low efficiency of star formation when only molecular hydrogen is used to form stars. We note though, that this parameter has been shown to have little influence on the global properties of simulated galaxies ([Li2018](https://doi.org/10.3847/1538-4357/aac9b8) and [Brown2022](https://doi.org/10.1093/mnras/stac1164)).

With all the previous definitions, we have

$\begin{equation}
    \tau_\mathrm{star} = \frac{C_\mathrm{star}}{\sqrt{\rho_g}} \, ,
\end{equation}$
where

$\begin{equation}
    C_\mathrm{star} = \sqrt{\frac{3\pi}{32 \, G}} \, .
\end{equation}$

Given that $\rho_g$ is the density of an individual cold gas cloud, which is unresolved within a cell, we will simply use $\rho_g = \rho_\mathrm{cell}$.
"""
  ╠═╡ =#

# ╔═╡ 0baae21f-71b5-42fd-be33-31de085928d3
begin
    const ϵff     = 1.0
	const C_star  = sqrt(3π / 32u"G") / ϵff
	const c_star  = ustrip(t_u * l_u^-(3/2), C_star / sqrt(m_u))

	τ_star(ρ_cell) = c_star / sqrt(ρ_cell)
	ψ(fm, ρ_cell) = fm / τ_star(ρ_cell)
end;

# ╔═╡ a7c5650c-9bbc-44a1-8d51-d473864faae7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Recombination

The rate of recombination for hydrogen atoms can be written as ([Osterbrock2006](http://www.worldcat.org/oclc/60611705), [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55), [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x), and [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635))

$\begin{equation}
    \frac{d}{dt} n_a \biggr\rvert_\mathrm{recomb.} = \alpha_H(T) \, n_e \, n_i \, ,
\end{equation}$

where $\alpha_H(T)$ is the recombination coefficient, $n_e$ is the electron number density, and $n_i$ is the ionized hydrogen number density. We note that, by mass balance and within our model, $n_e = n_i$.

Using the conversion factor already mentioned, we can write

$\begin{align}
    \frac{d}{dt} f_a \biggr\rvert_\mathrm{recomb.} &= \frac{x_a \, m_p}{\rho_\mathrm{cell}} \, \alpha_H(T) \, n_e \, n_i \\
	&= \frac{x_a \, m_p}{\rho_\mathrm{cell}} \, \alpha_H(T) \, f_i \, \frac{\rho_\mathrm{cell}}{x_i \, m_p} \, f_i \, \frac{\rho_\mathrm{cell}}{x_i \, m_p} \\
    &= \alpha_H(T) \, f_i^{\,2} \, \frac{\rho_\mathrm{cell}}{m_p} \, \frac{x_a}{x_i^{\,2}}  \, .
\end{align}$

We can readily find fits for $\alpha_H(T)$ in the literature ([Seaton1959](https://doi.org/10.1093/mnras/119.2.81), [Black1981](https://doi.org/10.1093/mnras/197.3.553), [Verner1996](https://doi.org/10.1086/192284), and [Osterbrock2006](http://www.worldcat.org/oclc/60611705)). In particular, if we take $T = 10^4 \, \mathrm{K}$ for the ionized phase and use case B recombination assuming an optically thick cloud ([Nebrin2023](https://doi.org/10.3847/2515-5172/acd37a)), we get

$\begin{align}
    \alpha_H(10^4 \, \mathrm{K}) = \alpha_H = 2.6 \times 10^{-13} \, \mathrm{cm}^3 \, \mathrm{s}^{-1} \, .
\end{align}$

So, we can write

$\begin{equation}
    \frac{d}{dt} f_a \biggr\rvert_\mathrm{recomb.} = \alpha_H \, f_i^{\,2} \, \frac{\rho_\mathrm{cell}}{m_p} \, \frac{x_a}{x_i^{\,2}} = \frac{f_i}{\tau_\mathrm{rec}} \, ,
\end{equation}$

where the time scale $\tau_\mathrm{rec}$ is

$\begin{equation}
    \tau_\mathrm{rec} = \frac{m_p \, x_i^{\,2}}{\alpha_H \, f_i \, \rho_\mathrm{cell} \, x_a} = \frac{C_\mathrm{rec}}{f_i \, \rho_\mathrm{cell}} \, ,
\end{equation}$

with the constant

$\begin{equation}
	C_\mathrm{rec} = \frac{m_p}{\alpha_H} \, \frac{x_i^{\,2}}{x_a} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 49d1a7f7-2bf2-4472-94df-6247b9237ddd
begin
	const αH    = 2.6e-13u"cm^3 * s^-1"
	const C_rec = m_u / αH
	const c_rec = ustrip(t_u * l_u^-3, C_rec / m_u)

	τ_rec(fi, ρ_cell) = c_rec / (fi * ρ_cell)
end;

# ╔═╡ 2d6fcef9-df4b-4eec-be5c-a8865c3a1b76
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Condensation

From at least [Hollenbach1971a](https://doi.org/10.1086/150754) and [Hollenbach1971b](https://doi.org/10.1086/150755) we know that the rate of molecular hydrogen formation due to the condensation of atomic gas in the surface of dust grains is

$\begin{equation}
    \frac{d}{dt} n_m \biggr\rvert_\mathrm{cond.} = R_d \, n_H \, n_a  \, ,
\end{equation}$

where $R_d$ is the formation rate coefficient of $H_2$ on dust grain, $n_H$ is the hydrogen nucleus number density, and $n_a$ is the atomic hydrogen number density. $n_H$ comes from the assumption that the number density of dust grains is proportional to it ([Hollenbach1971b](https://doi.org/10.1086/150755) Section II).

In the literature, $n_H$ is generally taken as $n_H = n_a + 2 \, n_m$, because most studies consider cold gas clouds ($T \sim 100 \, \mathrm{K}$) dominated by molecular and atomic hydrogen gas. In contrast, we consider atomic, molecular, and ionized gas within our gas cell; so, we will use $n_H = n_a + 2\, n_m + n_i$.

We also note that the expression for $\frac{d}{dt} n_m$ is only used in equilibrium equations in most of the early works ([Hollenbach1971a](https://doi.org/10.1086/150754), [Hollenbach1971b](https://doi.org/10.1086/150755), [Jura1974](https://doi.org/10.1086/152975), [Jura1975](https://doi.org/10.1086/153545), [Black1987](https://doi.org/10.1086/165740), [Sternberg1988](https://doi.org/10.1086/166664), and [Goldshmidt1995](https://doi.org/10.1086/175168)), while first appearing in an actual differential equation that does not assume equilibrium in [Draine1996](https://doi.org/10.1086/177689).

Using the conversion factor already mentioned we can write

$\begin{align}
    \frac{d}{dt} f_m \biggr\rvert_\mathrm{cond.} &= \frac{2 \, m_p \, x_m}{\rho_\mathrm{cell}} \, R_d \, \left( n_a + 2\, n_m + n_i \right) \, n_a \\
	&= \frac{2 \, m_p \, x_m}{\rho_\mathrm{cell}} \, R_d \left( f_a \, \frac{\rho_\mathrm{cell}}{m_p \, x_a} + 2 \, f_m \, \frac{\rho_\mathrm{cell}}{2 \, m_p \, x_m} + f_i \, \frac{\rho_\mathrm{cell}}{m_p \, x_i} \right) \, f_a \, \frac{\rho_\mathrm{cell}}{m_p \, x_a} \\
	&= 2 \, R_d \left( \frac{f_a}{x_a} + \frac{f_m}{x_m} + \frac{f_i}{x_i} \right) \, f_a \, \frac{\rho_\mathrm{cell}}{m_p} \, \frac{x_m}{x_a} \\
	&= \frac{f_a}{\tau_\mathrm{cond}} \, ,
\end{align}$

where the time scale $\tau_\mathrm{cond}$ is

$\begin{equation}
    \tau_\mathrm{cond} := \frac{m_p \, x_a}{2 \, R_d \, \rho_\mathrm{cell} \, x_m \left( \dfrac{f_a}{x_a} + \dfrac{f_m}{x_m} + \dfrac{f_i}{x_i} \right)} \, .
\end{equation}$

A table with several values for $R_d$ is presented below. We note that more than one value for $R_d$ and its dependence on other parameters may be discussed within each reference. In the table, we reflect the fiducial value used by each author

| $R_d \,\, [10^{-17} \, \mathrm{cm^3 \, s^{-1}}]$ | Reference |
|:-----:|:--------------------------------------------------:|
| $1$   | [Hollenbach1971b](https://doi.org/10.1086/150755)  |
| $5$   | [Jura1974](https://doi.org/10.1086/152975)         |
| $3$   | [Jura1975](https://doi.org/10.1086/153545)         |
| $9$   | [Black1987](https://doi.org/10.1086/165740)        |
| $3$   | [Goldshmidt1995](https://doi.org/10.1086/175168)   |
| $6$   | [Draine1996](https://doi.org/10.1086/177689)       |
| $3.5$ | [Wolfire2008](https://doi.org/10.1086/587688)      |

Theoretical and experimental work has shown that $R_d \propto T^{1/2}$ ([Black1987](https://doi.org/10.1086/165740)) and $R_d \propto n_\mathrm{dust}$ ([Hollenbach1971b](https://doi.org/10.1086/150755)). Assuming the simplest dust model $n_\mathrm{dust} \propto Z$, we have $R_d \propto T^{1/2} \, Z$ ([Pelupessy2006](https://doi.org/10.1086/504366) and [Wolfire2008](https://doi.org/10.1086/587688)).

Following previous prescriptions ([Pelupessy2006](https://doi.org/10.1086/504366), [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55), [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x), [Mollá2017](https://doi.org/10.1093/mnras/stx419), [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635)) we will only scale $R_d$ with the metallicity and use $\sim 100\,\mathrm{K}$ for the temperature of the cold neutral gas.

Our reference value for $R_d$ at $Z = Z_\odot$ is ([Wolfire2008](https://doi.org/10.1086/587688))

$\begin{equation}
    R_\odot := R_d(Z = Z_\odot) = 3.5 \times 10^{-17} \, \mathrm{cm^3 \, s^{-1}} \, ,
\end{equation}$

so, we have

$\begin{equation}
    R_d = Z \, \frac{R_\odot}{Z_\odot} \, .
\end{equation}$

The value of $R_\odot$ within $1 \, \mathrm{dex}$ should not affect significantly the results, but we will use an adjustable global factor (like the clumping factor $C_\rho$ in [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) and [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x)) to account for all the uncertainties.

The fact that the final expression for $R_d$ would give $0$ when $Z = 0$ is problematic because $\tau_\mathrm{cond}$ would diverge. The problem originates in the incorrect assumption that the only conversion channel for $\mathrm{HI} \rightarrow \mathrm{H_2}$ is the condensation on the surface of dust grains. To solve this we do the replacement $Z \rightarrow Z + Z_\mathrm{eff}$ with $Z_\mathrm{eff} = 10^{-3} \, Z_\odot$, which eliminates the divergence and takes into account the molecular formation that occurs below $10^{-3} \, Z_\odot$ ([Glover2007](https://doi.org/10.1086/519445)).

So, we finally have

$\begin{align}
    \tau_\mathrm{cond} &= \frac{m_p \, x_a \, Z_\odot}{2 \, R_\odot \, C_\rho \, (Z + Z_\mathrm{eff}) \, \rho_\mathrm{cell} \, x_m \left( \dfrac{f_a}{x_a} + \dfrac{f_m}{x_m} + \dfrac{f_i}{x_i} \right)} \\
	&= \frac{C_\mathrm{cond}}{(Z + Z_\mathrm{eff}) \, \rho_\mathrm{cell} \left( \dfrac{f_a}{x_a} + \dfrac{f_m}{x_m} + \dfrac{f_i}{x_i} \right)} \, ,
\end{align}$

where

$\begin{equation}
	C_\mathrm{cond} = \frac{m_p \, x_a \, Z_\odot}{2 \, R_\odot \, C_\rho \, x_m} \, .
\end{equation}$

For the solar metallicity, we find several values in the literature

| $Z_\odot$ | Reference                                                            |
|:--------:|:---------------------------------------------------------------------:|
| $0.0169$ | [Draine1996](https://doi.org/10.1086/177689)                          |
| $0.0153$ | [Grevesse1998](https://doi.org/10.1023/A:1005161325181)               |
| $0.0122$ | [Asplund2006](https://doi.org/10.1553/cia147s76) and [Grevesse2007](https://doi.org/10.1007/s11214-007-9173-7)                                        |
| $0.0134$ | [Asplund2009](https://doi.org/10.1146/annurev.astro.46.060407.145222) |
| $0.0141$ | [Lodders2009](https://doi.org/10.1007/978-3-540-88055-4_34)           |
| $0.0127$ | $\texttt{AREPO}$                                                                 |
| $0.0196$ | [Steiger2015](https://doi.org/10.3847/0004-637X/816/1/13)             |

To keep it consistent with the $\texttt{AREPO}$ codebase, we will use $Z_\odot = 0.0127$, noting that is only $35\%$ off the largest value in the list ([Steiger2015](https://doi.org/10.3847/0004-637X/816/1/13)).

We found that we need a $C_\rho \approx 100$ to form a disk, which is a reasonable value for the clumping factor under certain circumstances ([Micic2012](https://doi.org/10.1111/j.1365-2966.2012.20477.x)). 
"""
  ╠═╡ =#

# ╔═╡ 568d9fe3-6716-4a9a-ba1d-b9e6fd039150
begin
    const Rsun     = 3.5e-17u"cm^3 * s^-1"
    const Zsun     = 0.0127
	const Zeff     = 1e-3 * Zsun
	const Cρ       = 100.0
	const C_cond   = (m_u * Zsun) / (2 * Rsun * Cρ)
	const c_cond   = ustrip(t_u * l_u^-3, C_cond / m_u)

	τ_cond(fs, ρ_cell, Z) = c_cond / (ρ_cell * (Z + Zeff) * (1.0 - fs))
end;

# ╔═╡ 43ee281f-1a16-445d-894d-23e0319b1fd0
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Photodissociation efficiency

We define the disassociated mass rate per unit of created stellar mass as

$\begin{equation}
    \eta_\mathrm{diss} = \frac{\dot{M}_\mathrm{diss}}{\mathrm{SFR}}  \, ,
\end{equation}$

and the ionized mass rate per unit of created stellar mass as

$\begin{equation}
    \eta_\mathrm{ion} = \frac{\dot{M}_\mathrm{ion}}{\mathrm{SFR}}  \, .
\end{equation}$

The mass rates, $\dot{M}_\mathrm{diss}$ and $\dot{M}_\mathrm{ion}$, can be computed from the photon production rate

$\begin{align}
    \dot{M}_\mathrm{diss} &= \dot{N}_\mathrm{diss} \, c_\mathrm{diss} \, , \\
    \dot{M}_\mathrm{ion} &= \dot{N}_\mathrm{ion} \, c_\mathrm{ion} \, ,
\end{align}$

where $\dot{N}_\mathrm{diss}$ is the number of photodissociating photons produced per unit time (in the Lyman–Werner band, $912\,\mathrm{Å}$ to $1107\,\mathrm{Å}$), $\dot{N}_\mathrm{ion}$ the number of ionizing photons produced per unit time (between $0$ and $912\,Å$), and $c_\mathrm{diss}$ and $c_\mathrm{ion}$ are the unit conversion factors (proton mass into solar mass). The factors $c_\mathrm{diss}$ and $c_\mathrm{ion}$ allow us to consider that the reaction may not be $100\%$ efficient too.

For the ionization reaction, each photon will produce one proton and we assume $100\%$ efficiency.
"""
  ╠═╡ =#

# ╔═╡ d90371c3-3359-41cb-9bc6-27484beddb3c
const c_ion = ustrip(1.0u"mp" |> u"Msun")

# ╔═╡ 1f59c4ca-37c4-4492-8670-4bfdda78bd65
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
For the molecular dissociation reaction, we will consider the numerical factor given by [Draine1996](https://doi.org/10.1086/177689), where it is shown that dust grains may absorb up to $\sim \! 60$ percent of the photons capable of dissociating hydrogen molecules and a large fraction of the remaining photons excite different rotational and vibrational states, reducing their dissociation probability to $\sim \! 15$ percent. So, we end up with an efficiency factor of $0.4 \times 0.15 = 0.06$.
$c_\mathrm{diss}$ has an extra factor of two because each photon contributes with two protons (from the dissociated molecule) to the atomic gas.
"""
  ╠═╡ =#

# ╔═╡ 1657d7c2-87af-4079-8cf7-7a1164a95cc8
const c_diss = 0.4 * 0.15 * 2.0 * ustrip(1.0u"mp" |> u"Msun")

# ╔═╡ 4d09ed45-423c-4bd6-802d-59389a966d2e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The number of photons can be computed from $Q(t', Z)$ which is defined as the number of photons produced by a stellar population of one solar mass, of age $t'$, and metallicity $Z$, per unit of time. So, we have the relation

$\begin{equation}
    \dot{N}(t) = \int_0^t \mathrm{SFR}(t - t') \, Q(t', Z(t - t')) \, \mathrm{d}t' \, ,
\end{equation}$

where $\mathrm{SFR}(t - t')$ is the instantaneous SFR at the moment of birth of the stellar population of current age $t'$ (in units of solar mass per unit of time), and $Z = Z(t - t')$ is defined as the metallicity at that moment. Because most of the contribution to the integral comes from young blue stars that die in the first $10 \ \mathrm{to} \ 100 \, \mathrm{Myr}$, it is possible to approximate

$\begin{align}
    \mathrm{SFR}(t - t') &\approx \mathrm{SFR}(t) \, , \\
	Z(t - t') &\approx Z(t) \, .
\end{align}$

So, we end up with the general form

$\begin{equation}
    \eta = c \, \frac{\dot{N}}{\mathrm{SFR}} = c \, \int_0^t Q(t', Z) \, \mathrm{d}t'\, ,
\end{equation}$

The value of $Q$ can be calculated using

$\begin{equation}
    Q(t', Z) = \int_{\lambda_1}^{\lambda_2} \frac{L_\lambda(t', Z)}{E_\lambda} \mathrm{d}\lambda = \int_{\lambda_1}^{\lambda_2} \frac{\lambda \, L_\lambda(t', Z)}{h \, c} \, \mathrm{d}\lambda \, ,
\end{equation}$

where $L_\lambda(t', Z)$ is the luminosity per unit of wavelength of a stellar population of one solar mass, age $t'$, and metallicity $Z$; and $E_\lambda = \lambda / (h \, c)$ is the energy of a photon of wavelength $\lambda$. So, integrating between the wavelength of interest, we get the number of photons produced per unit of time for a stellar population of the given characteristics.

The luminosity will not only depend on the age and metallicity of the population, but on the initial mass function (IMF)) too, so in principle for each IMF, $t'$, and $Z$ we have a function of luminosity versus $\lambda$.

Using the values from PopStar by [Mollá2009](https://doi.org/10.1111/j.1365-2966.2009.15160.x) we compute a table of $Q$ for six IMFs ([Salpeter1955](https://doi.org/10.1086/145971) in two mass ranges: $0.85\,\mathrm{M}_\odot$ to $120\,\mathrm{M}_\odot$ and $0.15\,\mathrm{M}_\odot$ to $100\,\mathrm{M}_\odot$, [Kroupa2001](https://doi.org/10.1046/j.1365-8711.2001.04022.x), [Ferrini1990](https://ui.adsabs.harvard.edu/abs/1990A%26A...231..391F), and [Chabrier2003](https://doi.org/10.1086/374879)), six metallicities (0.0001, 0.0004, 0.004, 0.008, 0.02, and 0.05) and for ages between $0.1\,\mathrm{Myr}$, and $15\,\mathrm{Gyr}$.
"""
  ╠═╡ =#

# ╔═╡ 35083c71-0a7c-4b97-94ef-4d06ecbd2ee8
begin

    # Raw luminosity data from
    # https://www.fractal-es.com/PopStar/#download (PopStar2009)
    data_folders = readdir("./data/luminosity", join=true)

    # Regex patterns to parse the filenames, IMF_mlow_mup_zXXXX_tXXXX
    patterns = [
        r"(?<=spneb_)(.*?)(?=_z)",  # IMF_mlow_mup
        r".+?(?=_)",                # IMF
        r"(?<=_)(.*?)(?=_)",        # mlow
        r"[^_]+$",                  # mup
        r"(?<=z)(.*?)(?=_)",        # metallicity
        r"(?<=t)(.*)",              # log(age)
    ]

    # Wavelength range for the photodissociation of Hydrogen molecules
    λ_range_d = (912.0u"Å", 1107.0u"Å")
	# Wavelength range for the photoionization of Hydrogen atoms
    λ_range_i = (0.0u"Å", 912.0u"Å")

	# Unit conversion factor for the spectral energy distributions
	u_sed = 3.82e33u"erg * s^-1 * Å^-1"

    data_in_files = Vector{DataFrame}(undef, length(data_folders))

    @inbounds for (i, dir) in pairs(data_folders)
        files = readdir(dir, join=true)

        IMF_mlow_mup = getfield.(match.(patterns[1], basename.(files)), :match)
        IMF          = getfield.(match.(patterns[2], IMF_mlow_mup), :match)
        mlow         = getfield.(match.(patterns[3], IMF_mlow_mup), :match)
        mup          = getfield.(match.(patterns[4], IMF_mlow_mup), :match)
        Zmet         = getfield.(match.(patterns[5], basename.(files)), :match)
        ages         = getfield.(match.(patterns[6], basename.(files)), :match)

        Q_diss = Vector{Quantity}(undef, length(files))
		Q_ion  = Vector{Quantity}(undef, length(files))

        @inbounds for (i, file) in pairs(files)

			# Load data
            data = readdlm(file)
            df = identity.(DataFrame(data, [:λ, :Ls, :Lneb, :Ltot]))

			# Set units
            # Wavelength
            df[!, 1] = df[!, 1] .* u"Å"
            # Stellar spectral energy distributions per unit wavelength
            df[!, 2] = df[!, 2] .* u_sed
            # Nebular spectral energy distributions per unit wavelength
            df[!, 3] = df[!, 3] .* u_sed
            # Total spectral energy distributions per unit wavelength
            df[!, 4] = df[!, 4] .* u_sed

            # Spectral energy distribution integration
			let
                λ         = @subset(df, λ_range_d[1] .< :λ .< λ_range_d[2])
                integrand = λ[!, 1] .* λ[!, 2] ./ (Unitful.h * Unitful.c)
                Q_diss[i] = trapz(λ[!, 1], integrand) |> u"s^-1"
            end
            let
                λ         = @subset(df, λ_range_i[1] .< :λ .< λ_range_i[2])
                integrand = λ[!, 1] .* λ[!, 2] ./ (Unitful.h * Unitful.c)
                Q_ion[i]  = trapz(λ[!, 1], integrand) |> u"s^-1"
            end


        end

        data_in_files[i] = identity.(
            DataFrame(
                :IMF     => uppercase.(IMF),        # Initial mass function
                :mlow    => parse.(Float64, mlow),  # Min. mass of the IMF
                :mup     => parse.(Float64, mup),   # Max. mass of the IMF
                :Zmet    => parse.(Float64, "0." .* Zmet),  # Metallicities
                :log_age => parse.(Float64, ages),          # Stellar ages
                :Q_diss  => Q_diss,  # Number of dissociating photons per unit time
				:Q_ion   => Q_ion,  # Number of dissociating photons per unit time
            )
        )
    end

    Q_data = sort(vcat(data_in_files...), [:IMF, :mlow, :Zmet, :log_age])

	# Separate the computed Q values by IMF
	Salpeter1955A = @select(
        @subset(Q_data, :IMF .== "SAL", :mlow .== 0.85),
        $(Not([:IMF, :mlow, :mup])),
    )
    Salpeter1955B = @select(
        @subset(Q_data, :IMF .== "SAL", :mlow .== 0.15),
        $(Not([:IMF, :mlow, :mup])),
    )
    Ferrini1990 = @select(
        @subset(Q_data, :IMF .== "FER"),
        $(Not([:IMF, :mlow, :mup])),
    )
    Kroupa2001 = @select(
        @subset(Q_data, :IMF .== "KRO"),
        $(Not([:IMF, :mlow, :mup])),
    )
    Chabrier2003 = @select(
        @subset(Q_data, :IMF .== "CHA"),
        $(Not([:IMF, :mlow, :mup])),
    )

    Q_by_imf = Dict(
		"Salpeter1955A" => Salpeter1955A,
		"Salpeter1955B" => Salpeter1955B,
		"Ferrini1990"   => Ferrini1990,
		"Kroupa2001"    => Kroupa2001,
		"Chabrier2003"  => Chabrier2003,
	)

	Q_data
end

# ╔═╡ 76a10b7c-1135-49f1-a298-36c59bb94b37
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
In what follows we will use the IMF by [Chabrier2003](https://doi.org/10.1086/374879)
"""
  ╠═╡ =#

# ╔═╡ 72794eb3-7f09-4445-906a-af9756dceebd
begin
	const IMF      = "Chabrier2003"             # Initial mass function
	const Q_df     = Q_by_imf[IMF]
	const Q_ages   = unique(Q_df[!, :log_age])  # List of stellar ages
    const Q_metals = unique(Q_df[!, :Zmet])     # List of metallicities
end;

# ╔═╡ 1f190b07-164b-4ba9-86f8-edcf6e7ab3a9
##################################################################################
# Compute η_diss and η_ion
#
# age: Age [Myr]
# Z:   Metallicity [dimensionless]
##################################################################################

function photodissociation_efficiency(age::Float64, Z::Float64)::NTuple{2,Float64}

	# Allocate memory for η_diss and η_ion
	η_diss = Matrix{Float64}(undef, length(Q_ages), length(Q_metals))
	η_ion  = Matrix{Float64}(undef, length(Q_ages), length(Q_metals))

	@inbounds for (i, log_age) in pairs(Q_ages)

		@inbounds for (j, Zmet) in pairs(Q_metals)

	        sub_df = @subset(Q_df, :Zmet .== Zmet, :log_age .<= log_age)

	        # Set the values of the axes, with an extra point,
			# to integrate from age 0 onwards
	        ages = [0.0, exp10.(sub_df[!, :log_age])...] .* u"yr"

			q_diss       = sub_df[!, :Q_diss]
			Q_diss       = [q_diss[1], q_diss...]
			η_diss[i, j] = uconvert(Unitful.NoUnits, trapz(ages, Q_diss) * c_diss)

			q_ion       = sub_df[!, :Q_ion]
			Q_ion       = [q_ion[1], q_ion...]
	        η_ion[i, j] = uconvert(Unitful.NoUnits, trapz(ages, Q_ion) * c_ion)

	    end

	end

	max_age = log10(age * 10^6)

	ifunc_η_diss = linear_interpolation(
		(Q_ages, Q_metals),
		η_diss,
		extrapolation_bc=Flat(),
	)
	ifunc_η_ion = linear_interpolation(
		(Q_ages, Q_metals),
		η_ion,
		extrapolation_bc=Flat(),
	)

	return ifunc_η_diss(max_age, Z), ifunc_η_ion(max_age, Z)

end;

# ╔═╡ 7df63fdf-ef41-44f2-8a55-e5c2c849029c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Mass recycling

There are two mass recycling parameters: $R$ which is defined as the mass fraction of a stellar population that is returned to the ISM under the instantaneous recycling hypothesis (stars under a certain mass live forever and stars above that mass die instantly), first explicitly used by [Schmidt1963](https://doi.org/10.1086/147553), and $Z_\mathrm{sn}$ which is the fraction of the returned gas that is composed of metals (the rest is assumed to be ionized gas). We only need $R$ but for completeness we will compute $Z_\mathrm{sn}$ too.

Notice that the instantaneous recycling hypothesis can be avoided by considering the lifetimes of the stars (using empirical relations) as it is done in [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635) (sections 2.2.2 and 2.2.3). This would effectively make $R$ time-dependent (the integrals below would have to be computed at each evaluation of the equations) increasing significantly the computational cost. To avoid this, we will simply assume $R$ constant within the ODEs integration time scales.

A stellar yield model gives the amount of each element (as a fraction of the stellar mass) that is returned to the ISM by stars with masses between $m$ and $m + \mathrm{d}m$, so $R$ can be computed as

$\begin{equation}
	R = \dfrac{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} (m - m_\mathrm{rem}(m)) \, \phi(m) \, \mathrm{d}m}{\int_{m_\mathrm{low}}^{m_\mathrm{high}} m \, \phi(m) \, \mathrm{d}m} \, ,
\end{equation}$

where $\phi(m)$ is the ISM, $m_\mathrm{low}$ and $m_\mathrm{high}$ are the extremes in the mass range of the ISM, $m_\mathrm{ir}$ is the mass limit for the instantaneous recycling hypothesis, and $m_\mathrm{rem}(m)$ is the remnant stellar mass given by the yield model ([Pipino2014](https://doi.org/10.1093/mnras/stu579) and [Ascasibar2015](https://doi.org/10.1093/mnras/stv098)).

Notice that the denominator in the expression for $R$ is the total mass of the stellar population modeled by $\phi(m)$, so it is just a normalization. This is needed because the IMF is in general defined except for a global constant.

Using the same notation, we can calculate $Z_\mathrm{sn}$ as

$\begin{equation}
	Z_\mathrm{sn} = \dfrac{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} m \, f_Z \, \phi(m) \, \mathrm{d}m}{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} (m - m_\mathrm{rem}(m)) \, \phi(m) \, \mathrm{d}m} \, ,
\end{equation}$

where $f_Z$ is the fraction of the stellar mass that is returned to the ISM as metals.

Some traditional choices for the masses are $m_\mathrm{ir} = 8 \, M_\odot$, $m_\mathrm{low} = 0.08 \, M_\odot$ (the limit for hydrogen fusion), and $m_\mathrm{high} = 100 \, M_\odot$ (the order of magnitude for the upper limit of validity of the yield models given by observational limitations).

[Ascasibar2015](https://doi.org/10.1093/mnras/stv098) got $R \approx 0.18$ and $Z_\mathrm{sn} \approx 0.09$, using the yield model of [Woosley1995](https://doi.org/10.2172/115557) and the IMF of [Kroupa2001](https://doi.org/10.1046/j.1365-8711.2001.04022.x), even though it is not clear which mass limits were used.

The two previous equations were taken from [_Nucleosynthesis and Chemical Evolution of Galaxies_](https://doi.org/10.1017/CBO9780511812170) by Bernard Pagel (eq. 7.24 and 7.26), and [_Chemical Evolution
of Galaxies_](https://doi.org/10.1007/978-94-010-0967-6) by Francesca Matteucci (eq. 2.74).

We will consider the following stellar yield models,

  -  [Woosley et al. 1995](https://doi.org/10.2172/115557)
  -  [Portinari et al. 1998](https://doi.org/10.48550/arXiv.astro-ph/9711337)
  -  [Chieff et al. 2004](https://doi.org/10.1086/392523)
  -  [Kobayashi et al. 2006](https://doi.org/10.1086/508914)
  -  [Heger et al. 2010](https://doi.org/10.1088/0004-637X/724/1/341)
  -  [Limongi et al. 2012](https://doi.org/10.1088/0067-0049/199/2/38)

compiled by [Mollá2015](https://doi.org/10.1093/mnras/stv1102), which are summarized in the following table, where
  - `model`: Stellar yield model.
  - `s_Z`: Metallicity of the stellar population modeled by the IMF.
  - `s_m`: Stellar mass.
  - `m_rem`: Remnant mass after stellar death.
  - `zf_rem`: Fraction of the stellar mass ejected as metals to the ISM.
"""
  ╠═╡ =#

# ╔═╡ 0895d464-a029-410d-9e7d-89cfac2d1615
begin
    # Raw stellar yields from
    # Mollá et al. 2015 (https://doi.org/10.1093/mnras/stv1102)
    sy_files = readdir("./data/stellar_yields", join=true)

	columns = [
		:s_Z, :s_m, :H, :D, :He3, :He4, :C12, :O16, :Ne20, :Mg24,
		:Si28, :S32, :Ca40, :Fe56, :m_rem, :C13s, :N14s,
    ]
	regex = r"^(.+?)(\.[^.]*$|$)"
    sy_data = DataFrame[]

    for file in sy_files
        data = readdlm(file, skipstart=1)
        df   = identity.(DataFrame(data, columns))
        name = uppercase(
            getindex(getproperty(match(regex, basename(file)), :captures), 1),
        )

        df[!, 2]  = df[!, 2] .* u"Msun"  # Stellar mass
        df[!, 15] = df[!, 15] .* u"Msun" # Remnant mass

        insertcols!(df, 1, :model => fill(name, length(df[!, 1])))

        # See eq. 1 and 2 from
		# Mollá et al. 2015 (https://doi.org/10.1093/mnras/stv1102)
        @transform! df begin
            :zf_rem = (
                :C12 .+ :O16 .+ :Ne20 .+ :Mg24 .+ :Si28 .+ :S32 .+ :Ca40 .+ :Fe56
				.+ :C13s .+ :N14s .+ (1 .- :m_rem ./ :s_m) .* :s_Z
            )
        end

        select!(df, [:model, :s_Z, :s_m, :m_rem, :zf_rem])
        push!(sy_data, df)
    end

    sy_data = sort(vcat(sy_data...), [:model, :s_Z, :s_m])

    # Metallicities
    sy_metals = unique(sy_data[!, :s_Z])
    # Masses
    sy_masses = unique(sy_data[!, :s_m])

    # Columns:
    #   model:  Stellar yield model
    #   s_Z:    Metallicity of the stellar population modeled by the IMF
    #   s_m:    Stellar mass
    #   m_rem:  Remnant mass, after stellar death
    #   zf_rem: Fraction of the stellar mass ejected as metals to the ISM
    sy_data
end

# ╔═╡ ada0a462-e4cf-438d-b2a7-35d7280e955a
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Initial mass functions (IMF)

The initial mass function $\phi(m)$ gives the number of stars between masses $m$ and $m + \mathrm{d}m$.

For a given population of total mass $M$, we have

$\begin{equation}
    M = \int_{m_\mathrm{low}}^{m_\mathrm{high}} m \, \phi(m) \, \mathrm{d}m \, ,
\end{equation}$

which allows to normalize $\phi(m)$ for a population of mass $M$, within the range $[m_\mathrm{low}, m_\mathrm{high}]$.

There are many models for $\phi(m)$, but one of the simplest is the power law

$\begin{equation}
    \phi(m) = A \, m^{-\alpha}\, ,
\end{equation}$

where it is assumed that $[m] = M_\odot$, and $A$ is the normalization constant.

The following implementations don't have a specific choice for normalization (when possible $A = 1$), so they must be multiplied by the correct constant if one wants a given value of $M$.
"""
  ╠═╡ =#

# ╔═╡ 6a59ed2e-f040-4e10-bc55-91d2c1dcc97e
#####################################################################################
# Implementation of several IMF
#####################################################################################

begin
    # Salpeter 1955 (https://doi.org/10.1086/145971)
    # This model is valid for 0.4 <= m / M⊙ <= 10
    ϕSAL(m::Float64)::Float64 = m^(-2.35)

    # Miller et al. 1979 (https://doi.org/10.1086/190629)
    # This model is valid for 0.1 <= m / M⊙ <= 62
    const C1_MIL = 1.09
    const C2_MIL = -1.02
    ϕMIL(m::Float64)::Float64 = m^(-1) * exp(-(log10(m) - C2_MIL)^2 / (1 / C1_MIL))

    # Ferrini et al. 1990 (https://ui.adsabs.harvard.edu/abs/1990A%26A...231..391F)
    # Ferrini et al. 1992 (https://doi.org/10.1086/171066)
    # From the papers it is not clear the range of validity for this model, but it is
    # generally accepted that no model is valid outside 0.072 <= m / M⊙ <= 100
    function ϕFER(m::Real)::Real
        return m^(-0.52) * exp10(-sqrt(0.73 + log10(m) * (1.92 + 2.07 * log10(m))))
    end

    # Kroupa 1993 (https://doi.org/10.1093/mnras/262.3.545)
    # This model is valid for m / M⊙ >= 0.072
    function ϕKRO_93(m::Real)::Real
        if m < 0.5
            return m^(-1.2)
        elseif 0.5 <= m < 1
            return 0.5 * (m^-2.2)
        else
            return 0.5 * (m^-2.7)
        end
    end

    # Kroupa 2001 (https://doi.org/10.1046/j.1365-8711.2001.04022.x)
    # This model is valid for m / M⊙ >= 0.072
    function ϕKRO_01(m::Real)::Real
        if m < 0.08
            return m^-0.3
        elseif 0.08 <= m < 0.5
            return 0.08 * (m^-1.3)
        else
            return 0.0386 * (m^-2.35)
        end
    end

    # Chabrier 2003 (https://doi.org/10.1086/374879)
    # This model is valid for m / M⊙ <= 10
    # (above m = 1 M⊙ uses Salpeter (1955) results)
    function ϕCHA(m::Real)::Real
        if m <= 1
            return m^(-1) * exp(-(log10(m) - log10(0.22))^2 / (2 * 0.57^2))
        else
            return 0.514 * m^(-2.3)
        end
    end

    # Weidner 2005 (https://doi.org/10.1086/429867)
    # This model is valid for m / M⊙ >= 0.072
    function ϕWEI(m::Real)::Real
        if m < 0.08
            return m^(-0.3)
        elseif 0.08 <= m < 0.5
            return 0.08 * m^(-1.3)
        elseif 0.5 <= m < 1
            return 0.0386 * m^(-2.35)
        else
            return 0.0386 * m^(-2.7)
        end
    end

    # Millán-Irigoyen et al. 2020 (https://doi.org/10.1093/mnras/staa635)
    # This model is valid for 0.1 <= m / M⊙ <= 50
    function ϕMILLA(m::Real)::Real
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
    ϕKRO_01(m::Quantity)::Float64 = ϕKRO_01(ustrip(u"Msun", m))
    ϕCHA(m::Quantity)::Float64    = ϕCHA(ustrip(u"Msun", m))
    ϕWEI(m::Quantity)::Float64    = ϕWEI(ustrip(u"Msun", m))
    ϕMILLA(m::Quantity)::Float64  = ϕMILLA(ustrip(u"Msun", m))

    const imf_funcs = Dict(
        "Salpeter1955"        => ["SAL", ϕSAL],
        "Miller1979"          => ["MIL", ϕMIL],
        "Ferrini1990"         => ["FER", ϕFER],
        "Kroupa1993"          => ["KRO93", ϕKRO_93],
        "Kroupa2001"          => ["KRO01", ϕKRO_01],
        "Chabrier2003"        => ["CHA", ϕCHA],
        "Weidner2005"         => ["WEI", ϕWEI],
        "Millán-Irigoyen2020" => ["MILLA", ϕMILLA],
    )

	const yield_model_keys = Dict(
        "Woosley1995"   => "WOW",
        "Portinari1998" => "PCB",
        "Chieff2004"    => "CLI",
        "Kobayashi2006" => "KOB",
        "Heger2010"     => "HEG",
        "Limongi2012"   => "LIM",
    )
end;

# ╔═╡ 137c43b7-a541-4ca0-a6d0-0ad3579b345a
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
In what follows we will use the stellar yield model by [Portinari1998](https://doi.org/10.48550/arXiv.astro-ph/9711337)
"""
  ╠═╡ =#

# ╔═╡ c195ae82-0e44-4181-b863-bdd365059f4b
begin
	const Y_MODEL = "Portinari1998"  # Stellar yield model
	const M_LOW   = 0.1u"Msun"       # Lower mass limit for the ISM
	const M_HIGH  = 100.0u"Msun"     # Upper mass limit for the ISM
	const M_IR    = 8.0u"Msun"       # Mass limit for the IR hypothesis
end;

# ╔═╡ 14e326a5-d2bf-4873-9cf7-b57be9416da2
#####################################################################################
# Compute R and Zsn for a given stellar metallicity
#
# Z: Metallicity [dimensionless]
#####################################################################################

function recycled_fractions(Z::Float64)::NTuple{2,Float64}

	m_low    = ustrip(u"Msun", M_LOW)
	m_high   = ustrip(u"Msun", M_HIGH)
	m_ir     = ustrip(u"Msun", M_IR)
	imf_func = imf_funcs[IMF][2]
	Rs       = similar(sy_metals)
	Zsns     = similar(sy_metals)

	norm   = quadgk(m -> m * imf_func(m), m_low, m_high)[1] * u"Msun^2"
    sub_df = @subset(
		sy_data,
		:model .== yield_model_keys[Y_MODEL],
		m_ir .* u"Msun" .< :s_m .< m_high .* u"Msun",
	)

	@inbounds for (i, Z) in pairs(sy_metals)
		data = @subset(sub_df, :s_Z .== Z)

		stellar_mass  = data[:, :s_m]
		remnant_mass  = data[!, :m_rem]
		remnant_metal = data[!, :zf_rem]

	    R_int = trapz(
			stellar_mass,
			(stellar_mass .- remnant_mass) .* imf_func.(stellar_mass),
		)
	    Zsn_int = trapz(
			stellar_mass,
			stellar_mass .* remnant_metal .* imf_func.(stellar_mass),
		)

	    Rs[i]   = R_int / norm
		Zsns[i] = Zsn_int / R_int
	end

	R_vs_Z   = linear_interpolation(sy_metals, Rs, extrapolation_bc=Flat())
	Zsn_vs_Z = linear_interpolation(sy_metals, Zsns, extrapolation_bc=Flat())

	return R_vs_Z(Z), Zsn_vs_Z(Z)

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
	const N_EQU = 4  # Number of equations
	const N_PAR = 5  # Number of parameters

	# Index of each phase in the ODE solution  matrix
	const phase_name_to_index = Dict(
		"ionized"   => 1,
		"atomic"    => 2,
		"molecular" => 3,
		"stellar"   => 4,
	)
end;

# ╔═╡ f5a983bf-ef3a-4d5d-928d-da1097b91ee8
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Equations"
  ╠═╡ =#

# ╔═╡ bd8743d6-8f21-413d-835a-e543926baa09
#####################################################################################
# System of ODEs, where each equation has units of time⁻¹, and
#
# Ionized gas fraction:    fᵢ(t) = Mᵢ(t) / M_cell --> y[1]
# Atomic gas fraction:     fₐ(t) = Mₐ(t) / M_cell --> y[2]
# Molecular gas fraction:  fₘ(t) = Mₘ(t) / M_cell --> y[3]
# Stellar fraction:        fₛ(t) = Mₛ(t) / M_cell --> y[4]
#####################################################################################

function system!(dydt, ic, parameters, t)

    # Initial conditions
    fi, fa, fm, fs = ic

    # Parameters
	#
	# ρ_cell: Total cell density                                 [mp * cm⁻³]
	# Z:      Metallicity                                        [dimensionless]
	# η_diss: Photodissociation efficiency of hydrogen molecules [dimensionless]
	# η_ion:  Photoionization efficiency of hydrogen atoms       [dimensionless]
	# R:      Mass recycling fraction                            [dimensionless]
    ρ_cell, Z, η_diss, η_ion, R = parameters

    # Auxiliary equations
	recombination   = fi / τ_rec(fi, ρ_cell)
    cloud_formation = fa / τ_cond(fs, ρ_cell, Z)
	sfr             = ψ(fm, ρ_cell)

    # ODE system
	dydt[1] = -recombination + (η_ion + R) * sfr
    dydt[2] = -cloud_formation + recombination + (η_diss - η_ion) * sfr
    dydt[3] = cloud_formation - (1 + η_diss) * sfr
    dydt[4] = (1 - R) * sfr

end;

# ╔═╡ bae7dbde-0453-472d-bf1f-5171ecbb7967
#####################################################################################
# Compute the Lyapunov spectrum
#
# ic:         Initial conditions, [fi(0), fa(0), fm(0), fs(0)]
# base_parms: Parameters, (ρ [cm⁻³], Z [dimensionless], it [Myr])
#####################################################################################

function lyapunov_spectrum(
	ic::Vector{Float64},
	base_params::NTuple{3,Float64},
)::Vector{Float64}

	# Construct the parameters for the ODEs
	ρ             = base_params[1]  # Density [cm⁻³]
	Z             = base_params[2]  # Metallicity [dimensionless]
	it            = base_params[3]  # Integration time [Myr]
	η_diss, η_ion = photodissociation_efficiency(it, Z)
	R, _          = recycled_fractions(Z)

	parameters = [ρ, Z, η_diss, η_ion, R]
	ds = CoupledODEs(system!, ic, parameters)

	# Compute the Lyapunov spectrum
	return lyapunovspectrum(ds, 100)

end;

# ╔═╡ ae2c3c1c-aa3f-4be6-b209-38b9b5b77244
# ╠═╡ skip_as_script = true
#=╠═╡
lyapunov_spectrum(
	[0.8, 0.2, 0.0, 0.0],     # Initial conditions, [fi(0), fa(0), fm(0), fs(0)]
	(100.0, 0.5*Zsun, 10.0),  # Parameters, (ρ [cm⁻³], Z [dimensionless], it [Myr])
)
  ╠═╡ =#

# ╔═╡ 0d4877b6-0659-441e-81a9-43c59635a14f
# ╠═╡ skip_as_script = true
#=╠═╡
lyapunov_spectrum(
	[0.8, 0.2, 0.0, 0.0],      # Initial conditions, [fi(0), fa(0), fm(0), fs(0)]
	(1000.0, 0.5*Zsun, 10.0),  # Parameters, (ρ [cm⁻³], Z [dimensionless], it [Myr])
)
  ╠═╡ =#

# ╔═╡ 9ab0a10b-8165-401a-a2b6-c9726526a906
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Jacobian"
  ╠═╡ =#

# ╔═╡ 1b8af600-56eb-4508-bc52-aa4e138b4c7e
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

# ╔═╡ 80005099-7154-4306-9172-c9a168336e14
const JACOBIAN_FUNCTION = construct_jacobian(system!);

# ╔═╡ c291700e-3a84-49a7-85d6-592cfb3b1a11
#####################################################################################
# Evaluate the Jacobian
#
# J:          Matrix to save the results, it must have size N_EQU × N_EQU
# ic:         Initial condition, [fi(0), fa(0), fm(0), fs(0)]
# parameters: Parameters for the ODEs, [ρ, Z, η_diss, η_ion, R]
# t:          Unused variable to comply with with the DifferentialEquations.jl API
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

# ╔═╡ 40da0211-1132-412c-a769-05bd4b22b69b
#####################################################################################
# Compute the stiffness ratio
#
# ic:         Initial conditions, [fi(0), fa(0), fm(0), fs(0)]
# base_parms: Parameters, (ρ [cm⁻³], Z [dimensionless], it [Myr])
#####################################################################################

function stiffness_ratio(
	ic::Vector{Float64},
	base_params::NTuple{3,Float64},
)::Float64

	# Construct the parameters for the ODEs
	ρ             = base_params[1]  # Density [cm⁻³]
	Z             = base_params[2]  # Metallicity [dimensionless]
	it            = base_params[3]  # Integration time [Myr]
	η_diss, η_ion = photodissociation_efficiency(it, Z)
	R, _          = recycled_fractions(Z)

	parameters = [ρ, Z, η_diss, η_ion, R]
	J = Matrix{Float64}(undef, N_EQU, N_EQU)

	# Compute the Jacobian and store it in J
	jacobian!(J, ic, parameters)

	# Get the norm of the real part of the non-zero eigenvalues
	eigen_values = filter(x -> x > eps(typeof(x)), eigvals(J) .|> real .|> abs)

	# Compute the stiffness ratio
	return maximum(eigen_values) / minimum(eigen_values)

end;

# ╔═╡ 61db358f-a57d-485c-9bb3-cfd1297c0ca3
# ╠═╡ skip_as_script = true
#=╠═╡
stiffness_ratio(
	[0.8, 0.2, 0.0, 0.0],     # Initial conditions, [fi(0), fa(0), fm(0), fs(0)]
	(100.0, 0.5*Zsun, 10.0),  # Parameters, (ρ [cm⁻³], Z [dimensionless], it [Myr])
)
  ╠═╡ =#

# ╔═╡ 1ec99792-905d-4e1b-a413-ef58143d3c68
# ╠═╡ skip_as_script = true
#=╠═╡
md"## ODE Function"
  ╠═╡ =#

# ╔═╡ 4e1140be-3e61-47bd-9f9f-6b5dfbff6ae2
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
# ic:          Initial conditions, [fi(0), fa(0), fm(0), fs(0)]
# base_params: Parameters for the ODEs, [ρ_cell, Z]
# tspan:       Integration span, (ti, tf) [Myr]
# times:       Times at which the solution will be returned [Myr]
# args:        Positional arguments for the solver of DifferentialEquations.jl
# kwargs:      Keyword arguments for the solver of DifferentialEquations.jl
#####################################################################################

function integrate_model(
    ic::Vector{Float64},
	base_params::Vector{Float64},
	tspan::NTuple{2,Float64};
    times::Vector{Float64}=[tspan[2],],
    args::Tuple=(Rodas4P(),),
    kwargs::NamedTuple=(
		dense=false,
		force_dtmin=true,
		dtmin=(tspan[2] - tspan[1]) / 1.0e8,
		reltol=1.0e-10,
		maxiters=1.0e8,
		verbose=false,
	),
)::Vector{Vector{Float64}}

	ρ_cell        = base_params[1]
	Z             = base_params[2]
	η_diss, η_ion = photodissociation_efficiency(exp10(Q_ages[end] - 6), Z)
	R, _          = recycled_fractions(Z)

	parameters = [ρ_cell, Z, η_diss, η_ion, R]

	# println("With Julia:")
	# @show [Z, η_diss, η_ion]
	# println()

    sol = solve(ODEProblem(
		ode_function,
		ic,
		tspan,
		parameters,
	), args...; kwargs...)

    return sol(times).u

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
# params:  Parameters for the density PDF
# log_var: Selects which variable will be used and how the function will be divided
#   log_var == true:  s = ln(ρ/ρ₀) and logarithmic divisions
#   log_var == false: f = ρ/ρ₀ and linear divisions
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
# Density PDF acording to Burkhart (2018)
# https://doi.org/10.3847/1538-4357/aad002
#####################################################################################

function pBurkhart2018(s::Float64, params::PDF_params)::Float64

    b = params.b
    Ms = params.Ms
    α = params.α

    σs2 = log(1 + b^2 * Ms^2)
    s0 = -0.5 * σs2
    st = (α - 0.5) * σs2
    C = exp((α - 1) * 0.5 * α * σs2) / sqrt(2π * σs2)
    N = 1 / ((C * exp(-α * st)) / α + 0.5 + 0.5 * erf((2 * st + σs2) / sqrt(8 * σs2)))

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
ChaosTools = "608a59af-f2a3-5ad4-90b4-758bdf3122a7"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataFramesMeta = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
Trapz = "592b5752-818d-11e9-1e9a-2b8ca4a44cd1"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"

[compat]
ChaosTools = "~3.2.1"
DataFrames = "~1.7.0"
DataFramesMeta = "~0.15.3"
DelimitedFiles = "~1.9.1"
DifferentialEquations = "~7.14.0"
Interpolations = "~0.15.1"
PlutoUI = "~0.7.60"
QuadGK = "~2.11.1"
SpecialFunctions = "~2.4.0"
Symbolics = "~6.11.0"
TikzPictures = "~3.5.0"
Trapz = "~2.0.3"
Unitful = "~1.21.0"
UnitfulAstro = "~1.2.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.0"
manifest_format = "2.0"
project_hash = "3aacf9f0e91687db763bb4ab60f68b1ac2e03e64"

[[deps.ADTypes]]
git-tree-sha1 = "eea5d80188827b35333801ef97a40c2ed653b081"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.9.0"
weakdeps = ["ChainRulesCore", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesEnzymeCoreExt = "EnzymeCore"

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
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown"]
git-tree-sha1 = "b392ede862e506d451fc1616e79aa6f4c673dab8"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.38"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsDatesExt = "Dates"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsTestExt = "Test"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
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
git-tree-sha1 = "3640d077b6dafd64ceb8fd5c1ec76f7ca53bcf76"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.16.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
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
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "0dd7edaff278e346eb0ca07a7e75c9438408a3ce"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.10.3"
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
git-tree-sha1 = "a2c85f53ddcb15b4099da59867868bd40f005579"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.7.5"
weakdeps = ["SparseArrays"]

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bijections]]
git-tree-sha1 = "d8b0439d2be438a5f2cd68ec158fe08a7b2595b7"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.BoundaryValueDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "BandedMatrices", "ConcreteStructs", "DiffEqBase", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "Logging", "NonlinearSolve", "OrdinaryDiffEq", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseDiffTools"]
git-tree-sha1 = "22f36cc7199c3b0f34979c6925cf8ec05281409f"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "5.10.0"

    [deps.BoundaryValueDiffEq.extensions]
    BoundaryValueDiffEqODEInterfaceExt = "ODEInterface"

    [deps.BoundaryValueDiffEq.weakdeps]
    ODEInterface = "54ca160b-1b9f-5127-a996-1867f4bc2a2c"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "5a97e67919535d6841172016c9530fd69494e5ec"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.6"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "32abd86e3c2025db5172aa182b982debed519834"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.1"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.Chain]]
git-tree-sha1 = "9ae9be75ad8ad9d26395bf625dea9beac6d519f1"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.6.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "3e4b134270b372f2ed4d4d0e936aabaefc1802bc"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChaosTools]]
deps = ["Combinatorics", "DSP", "DataStructures", "Distances", "Distributions", "DynamicalSystemsBase", "IntervalRootFinding", "LinearAlgebra", "LombScargle", "Neighborhood", "Optim", "ProgressMeter", "Random", "Reexport", "Roots", "SpecialFunctions", "Statistics", "StatsBase"]
git-tree-sha1 = "c73aa5f86db02f08e538861ca73bfc5bedf94909"
uuid = "608a59af-f2a3-5ad4-90b4-758bdf3122a7"
version = "3.2.1"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

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
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

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
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"
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

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "0df00546373af8eee1598fb4b2ba480b1ebe895c"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.10"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataFramesMeta]]
deps = ["Chain", "DataFrames", "MacroTools", "OrderedCollections", "Reexport", "TableMetadataTools"]
git-tree-sha1 = "7042a6ad5910dc9edeae814e1110209752a1c996"
uuid = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
version = "0.15.3"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack"]
git-tree-sha1 = "066f60231c1b0ae2905ffd2651e207accd91f627"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.48.1"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "Statistics", "TruncatedStacktraces"]
git-tree-sha1 = "9c36ac68cf6c59a1c2569aadd7a62c47ed5c8eb5"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.158.1"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseSparseArraysExt = "SparseArrays"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "Functors", "LinearAlgebra", "Markdown", "NonlinearSolve", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "19dbd44d18bbfdfcf5e56c99cea9b0ed23df350a"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "3.9.1"
weakdeps = ["OrdinaryDiffEq", "OrdinaryDiffEqCore", "Sundials"]

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "ab1e6515ce15f01316a9825b02729fefa51726bd"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.23.0"

    [deps.DiffEqNoiseProcess.extensions]
    DiffEqNoiseProcessReverseDiffExt = "ReverseDiff"

    [deps.DiffEqNoiseProcess.weakdeps]
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
git-tree-sha1 = "d851f2ca05f3cec9988f081b047a778a58b48aaf"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.14.0"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "e61c544812d735a53ffbf406a80a73b7c2425532"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.6.9"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = "Enzyme"
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = "ForwardDiff"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = "PolyesterForwardDiff"
    DifferentiationInterfaceReverseDiffExt = "ReverseDiff"
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "66c4c81f259586e8f002eacebc177e1fb06363b0"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.11"
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
git-tree-sha1 = "d7477ecdafb813ddee2ae727afa94e9dcb5f3fb0"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.112"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "490392af2c7d63183bfa2c8aaa6ab981c5ba7561"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.14"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "bbf1ace0781d9744cb697fb856bd2c3f6568dadb"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.0"

[[deps.DynamicalSystemsBase]]
deps = ["ForwardDiff", "LinearAlgebra", "OrdinaryDiffEqTsit5", "Reexport", "Roots", "SciMLBase", "SparseArrays", "StateSpaceSets", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "5da2aa65f991b7c6d2052231d8dfb018d6c6af85"
uuid = "6e36e845-645a-534a-86f2-f5d4aa5a06b4"
version = "3.11.2"
weakdeps = ["StochasticDiffEq"]

    [deps.DynamicalSystemsBase.extensions]
    StochasticSystemsBase = "StochasticDiffEq"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EnzymeCore]]
git-tree-sha1 = "9c3a42611e525352e9ad5e4134ddca5c692ff209"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.4"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

[[deps.ErrorfreeArithmetic]]
git-tree-sha1 = "d6863c556f1142a061532e79f611aa46be201686"
uuid = "90fa49ef-747e-5e6f-a989-263ba693cf1a"
version = "0.5.2"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "8e18940a5ba7f4ddb41fe2b79b6acaac50880a86"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.26.1"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.Expronicon]]
deps = ["MLStyle", "Pkg", "TOML"]
git-tree-sha1 = "fc3951d4d398b5515f91d7fe5d45fc31dccb3c9b"
uuid = "6b7a57c9-7cc1-4fdf-b7f5-e857abae3636"
version = "0.8.5"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4d81ed14783ec49ce9f2e168208a12ce1815aa25"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+1"

[[deps.FastAlmostBandedMatrices]]
deps = ["ArrayInterface", "ArrayLayouts", "BandedMatrices", "ConcreteStructs", "LazyArrays", "LinearAlgebra", "MatrixFactorizations", "PrecompileTools", "Reexport"]
git-tree-sha1 = "3f03d94c71126b6cfe20d3cbcc41c5cd27e1c419"
uuid = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
version = "0.1.4"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "ab1b34570bcdf272899062e1a56285a53ecaae08"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.3.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "cbf5edddb61a43669710cbc2241bc08b36d9e660"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "2.0.4"

[[deps.FastRounding]]
deps = ["ErrorfreeArithmetic", "LinearAlgebra"]
git-tree-sha1 = "6344aa18f654196be82e62816935225b3b9abe44"
uuid = "fa42c844-2597-5d31-933b-ebd51ab2693f"
version = "0.3.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "b10bdafd1647f57ace3885143936749d61638c3b"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.26.0"

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
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "64d8e93700c7a3f28f717d265382d52fac9fa1c1"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.12"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "af49a0851f8113fcfae2ef5027c6d49d0acec39b"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.4"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "674ff0db93fffcd11a3573986e550d66cd4fd71f"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.5+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "1dc470db8b1131cfc7fb4c115de89fe391b9e780"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.12.0"

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

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8e070b599339d622e9a081d17230d74a5c473293"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.17"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "7c4195be1649ae622304031ed46a2f4df989f1eb"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.24"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.ICU_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "20b6765a3016e1fca0c9c93c80d50061b94218b7"
uuid = "a51ab1cf-af8e-5615-a023-bc2c838bba6b"
version = "69.1.0+0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InlineStrings]]
git-tree-sha1 = "45521d31238e87ee9f9732561bfee12d4eebd52d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.2"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "10bd689145d2c3b2a9844005d01087cc1194e79e"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "FastRounding", "LinearAlgebra", "Markdown", "Random", "RecipesBase", "RoundingEmulator", "SetRounding", "StaticArrays"]
git-tree-sha1 = "5ab7744289be503d76a944784bac3f2df7b809af"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.20.9"

[[deps.IntervalRootFinding]]
deps = ["ForwardDiff", "IntervalArithmetic", "LinearAlgebra", "Polynomials", "Reexport", "StaticArrays"]
git-tree-sha1 = "b92e9e2b356146918c4f3f3845571abcf0501594"
uuid = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
version = "0.5.11"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
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
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

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
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "25ee0be4d43d0269027024d75a24c24d6c6e590c"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.4+0"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "6e0fb267f2c869df59b6bf9f8cb0f2794f2d85e9"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.13.7"
weakdeps = ["FastBroadcast"]

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "07649c499349dad9f08dde4243a4c597064663e9"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.6.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "267dad6b4b7b5d529c76d40ff48d33f7e94cb834"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.6"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "854a9c268c43b77b0a27f22d7fab8d33cdb3a731"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "360f6039babd6e4d6364eff0d4fc9120834a2d9a"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.2.1"

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
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "dc82b3c3640a4362f09e4d7594b4387a820857e4"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.3"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "e4c3be53733db1051cc15ecf573b1042b3a712a1"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.3.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "FastLapackInterface", "GPUArraysCore", "InteractiveUtils", "KLU", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "PrecompileTools", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Sparspak", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "8941ad4bdd83768359801982e143008349b1a827"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "2.35.0"

    [deps.LinearSolve.extensions]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveEnzymeExt = "EnzymeCore"
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = "Pardiso"
    LinearSolveRecursiveArrayToolsExt = "RecursiveArrayTools"

    [deps.LinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

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

[[deps.LombScargle]]
deps = ["FFTW", "LinearAlgebra", "Measurements", "Random", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "d64a0ce7539181136a85fd8fe4f42626387f0f26"
uuid = "fc60dff9-86e7-5f2f-a8a0-edeadbb75bd9"
version = "1.0.3"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "8084c25a250e00ae427a379a5b607e7aed96a2dd"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.171"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "16a726dba99685d9e94c8d0a8f655383121fc608"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "3.0.1"
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

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "Requires"]
git-tree-sha1 = "bdcde8ec04ca84aef5b124a17684bf3b302de00e"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.11.0"

    [deps.Measurements.extensions]
    MeasurementsBaseTypeExt = "BaseType"
    MeasurementsJunoExt = "Juno"
    MeasurementsRecipesBaseExt = "RecipesBase"
    MeasurementsSpecialFunctionsExt = "SpecialFunctions"
    MeasurementsUnitfulExt = "Unitful"

    [deps.Measurements.weakdeps]
    BaseType = "7fbed51b-1ef5-4d67-9085-a4a9b26f478c"
    Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
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

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "8d39779e29f80aa6c071e7ac17101c6e31f075d7"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3eba928678787843e504c153a9b8e80d7d73ab17"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.5.0"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "3cebfc94a0754cc329ebc3bab1e6c89621e791ad"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.20"

[[deps.Neighborhood]]
deps = ["Distances", "NearestNeighbors", "Random", "Test"]
git-tree-sha1 = "fdea60ca30d724e76cc3b3d90d7f9d29d3d5cab5"
uuid = "645ca80c-8b79-4109-87ea-e1f58159d116"
version = "0.2.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "LazyArrays", "LineSearch", "LineSearches", "LinearAlgebra", "LinearSolve", "MaybeInplace", "PrecompileTools", "Preferences", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLJacobianOperators", "SciMLOperators", "Setfield", "SimpleNonlinearSolve", "SparseArrays", "SparseConnectivityTracer", "SparseMatrixColorings", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "4d8944f32db2b07a2bdf8477e878bcb9c9ea2308"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "3.15.1"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveBandedMatricesExt = "BandedMatrices"
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = "NLsolve"
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"

    [deps.NonlinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"

[[deps.OffsetArrays]]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ad31332567b189f508a3ea8957a2640b1147ab00"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "d9b79c4eed437421ac4285148fcadf42e0700e89"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.9.4"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqAdamsBashforthMoulton", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqExplicitRK", "OrdinaryDiffEqExponentialRK", "OrdinaryDiffEqExtrapolation", "OrdinaryDiffEqFIRK", "OrdinaryDiffEqFeagin", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqHighOrderRK", "OrdinaryDiffEqIMEXMultistep", "OrdinaryDiffEqLinear", "OrdinaryDiffEqLowOrderRK", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqNordsieck", "OrdinaryDiffEqPDIRK", "OrdinaryDiffEqPRK", "OrdinaryDiffEqQPRK", "OrdinaryDiffEqRKN", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqSSPRK", "OrdinaryDiffEqStabilizedIRK", "OrdinaryDiffEqStabilizedRK", "OrdinaryDiffEqSymplecticRK", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "Static", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "cd892f12371c287dc50d6ad3af075b088b6f2d48"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.89.0"

[[deps.OrdinaryDiffEqAdamsBashforthMoulton]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqLowOrderRK", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "8e3c5978d0531a961f70d2f2730d1d16ed3bbd12"
uuid = "89bda076-bce5-4f1c-845f-551c83cdda9a"
version = "1.1.0"

[[deps.OrdinaryDiffEqBDF]]
deps = ["ArrayInterface", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqSDIRK", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "b4498d40bf35da0b6d22652ff2e9d8820590b3c6"
uuid = "6ad6398a-0878-4a85-9266-38940aa047c8"
version = "1.1.2"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "FastBroadcast", "FastClosures", "FillArrays", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleUnPack", "Static", "StaticArrayInterface", "StaticArraysCore", "TruncatedStacktraces"]
git-tree-sha1 = "f4b7b11e2e4c0f4a7fe15d2edcec1e1ce2917d67"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "1.7.1"
weakdeps = ["EnzymeCore"]

    [deps.OrdinaryDiffEqCore.extensions]
    OrdinaryDiffEqCoreEnzymeCoreExt = "EnzymeCore"

[[deps.OrdinaryDiffEqDefault]]
deps = ["DiffEqBase", "EnumX", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "PrecompileTools", "Preferences", "Reexport"]
git-tree-sha1 = "c8223e487d58bef28a3535b33ddf8ffdb44f46fb"
uuid = "50262376-6c5a-4cf5-baba-aaf4f84d72d7"
version = "1.1.0"

[[deps.OrdinaryDiffEqDifferentiation]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqCore", "SciMLBase", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays"]
git-tree-sha1 = "e63ec633b1efa99e3caa2e26a01faaa88ba6cef9"
uuid = "4302a76b-040a-498a-8c04-15b101fed76b"
version = "1.1.0"

[[deps.OrdinaryDiffEqExplicitRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "TruncatedStacktraces"]
git-tree-sha1 = "4dbce3f9e6974567082ce5176e21aab0224a69e9"
uuid = "9286f039-9fbf-40e8-bf65-aa933bdc4db0"
version = "1.1.0"

[[deps.OrdinaryDiffEqExponentialRK]]
deps = ["DiffEqBase", "ExponentialUtilities", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqVerner", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "f63938b8e9e5d3a05815defb3ebdbdcf61ec0a74"
uuid = "e0540318-69ee-4070-8777-9e2de6de23de"
version = "1.1.0"

[[deps.OrdinaryDiffEqExtrapolation]]
deps = ["DiffEqBase", "FastBroadcast", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "fea595528a160ed5cade9eee217a9691b1d97714"
uuid = "becaefa8-8ca2-5cf9-886d-c06f3d2bd2c4"
version = "1.1.0"

[[deps.OrdinaryDiffEqFIRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLOperators"]
git-tree-sha1 = "795221c662698851328cb7787965ab4a180d9468"
uuid = "5960d6e9-dd7a-4743-88e7-cf307b64f125"
version = "1.1.1"

[[deps.OrdinaryDiffEqFeagin]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "a7cc74d3433db98e59dc3d58bc28174c6c290adf"
uuid = "101fe9f7-ebb6-4678-b671-3a81e7194747"
version = "1.1.0"

[[deps.OrdinaryDiffEqFunctionMap]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "925a91583d1ab84f1f0fea121be1abf1179c5926"
uuid = "d3585ca7-f5d3-4ba6-8057-292ed1abd90f"
version = "1.1.1"

[[deps.OrdinaryDiffEqHighOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "103e017ff186ac39d731904045781c9bacfca2b0"
uuid = "d28bc4f8-55e1-4f49-af69-84c1a99f0f58"
version = "1.1.0"

[[deps.OrdinaryDiffEqIMEXMultistep]]
deps = ["DiffEqBase", "FastBroadcast", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Reexport"]
git-tree-sha1 = "9f8f52aad2399d7714b400ff9d203254b0a89c4a"
uuid = "9f002381-b378-40b7-97a6-27a27c83f129"
version = "1.1.0"

[[deps.OrdinaryDiffEqLinear]]
deps = ["DiffEqBase", "ExponentialUtilities", "LinearAlgebra", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "0f81a77ede3da0dc714ea61e81c76b25db4ab87a"
uuid = "521117fe-8c41-49f8-b3b6-30780b3f0fb5"
version = "1.1.0"

[[deps.OrdinaryDiffEqLowOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "d4bb32e09d6b68ce2eb45fb81001eab46f60717a"
uuid = "1344f307-1e59-4825-a18e-ace9aa3fa4c6"
version = "1.2.0"

[[deps.OrdinaryDiffEqLowStorageRK]]
deps = ["Adapt", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "StaticArrays"]
git-tree-sha1 = "590561f3af623d5485d070b4d7044f8854535f5a"
uuid = "b0944070-b475-4768-8dec-fb6eb410534d"
version = "1.2.1"

[[deps.OrdinaryDiffEqNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "PreallocationTools", "RecursiveArrayTools", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "StaticArrays"]
git-tree-sha1 = "a2a4119f3e35f7982f78e17beea7b12485d179e9"
uuid = "127b3ac7-2247-4354-8eb6-78cf4e7c58e8"
version = "1.2.1"

[[deps.OrdinaryDiffEqNordsieck]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "ef44754f10e0dfb9bb55ded382afed44cd94ab57"
uuid = "c9986a66-5c92-4813-8696-a7ec84c806c8"
version = "1.1.0"

[[deps.OrdinaryDiffEqPDIRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "Reexport", "StaticArrays"]
git-tree-sha1 = "a8b7f8107c477e07c6a6c00d1d66cac68b801bbc"
uuid = "5dd0a6cf-3d4b-4314-aa06-06d4e299bc89"
version = "1.1.0"

[[deps.OrdinaryDiffEqPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "Reexport"]
git-tree-sha1 = "da525d277962a1b76102c79f30cb0c31e13fe5b9"
uuid = "5b33eab2-c0f1-4480-b2c3-94bc1e80bda1"
version = "1.1.0"

[[deps.OrdinaryDiffEqQPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "332f9d17d0229218f66a73492162267359ba85e9"
uuid = "04162be5-8125-4266-98ed-640baecc6514"
version = "1.1.0"

[[deps.OrdinaryDiffEqRKN]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "41c09d9c20877546490f907d8dffdd52690dd65f"
uuid = "af6ede74-add8-4cfd-b1df-9a4dbb109d7a"
version = "1.1.0"

[[deps.OrdinaryDiffEqRosenbrock]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "96b47cdd12cb4ce8f70d701b49f855271a462bd4"
uuid = "43230ef6-c299-4910-a778-202eb28ce4ce"
version = "1.2.0"

[[deps.OrdinaryDiffEqSDIRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "f6683803a58de600ab7a26d2f49411c9923e9721"
uuid = "2d112036-d095-4a1e-ab9a-08536f3ecdbf"
version = "1.1.0"

[[deps.OrdinaryDiffEqSSPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "StaticArrays"]
git-tree-sha1 = "7dbe4ac56f930df5e9abd003cedb54e25cbbea86"
uuid = "669c94d9-1f4b-4b64-b377-1aa079aa2388"
version = "1.2.0"

[[deps.OrdinaryDiffEqStabilizedIRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "348fd6def9a88518715425025eadd58517017325"
uuid = "e3e12d00-db14-5390-b879-ac3dd2ef6296"
version = "1.1.0"

[[deps.OrdinaryDiffEqStabilizedRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "1b0d894c880e25f7d0b022d7257638cf8ce5b311"
uuid = "358294b1-0aab-51c3-aafe-ad5ab194a2ad"
version = "1.1.0"

[[deps.OrdinaryDiffEqSymplecticRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "4e8b8c8b81df3df17e2eb4603115db3b30a88235"
uuid = "fa646aed-7ef9-47eb-84c4-9443fc8cbfa8"
version = "1.1.0"

[[deps.OrdinaryDiffEqTsit5]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "96552f7d4619fabab4038a29ed37dd55e9eb513a"
uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a"
version = "1.1.0"

[[deps.OrdinaryDiffEqVerner]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "81d7841e73e385b9925d5c8e4427f2adcdda55db"
uuid = "79d7bb75-1356-48c1-b8c0-6832512096c2"
version = "1.1.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "6d38fea02d983051776a856b7df75b30cf9a3c1f"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.16"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "3aa2bb4982e575acd7583f01531f241af077b163"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.2.13"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.Poppler_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "02148a0cb2532f22c0589ceb75c110e168fb3d1f"
uuid = "9c32591e-4766-534b-9725-b71a8799265b"
version = "21.9.0+0"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff"]
git-tree-sha1 = "6c62ce45f268f3f958821a1e5192cf91c75ae89c"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.24"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "cb420f77dc474d23ee47ca8d14c90810cafe69e7"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.6"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "8f6bc219586aef8baf0ff9a5fe16ee9c70cb65e4"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.2"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "cda3b045cf9ef07a08ad46731f5a3165e56cf3da"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.1"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "4743b43e5a9c4a2ede372de7061eed81795b12e7"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.0"

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

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "b034171b93aebc81b3e1890a036d13a9c4a9e3e0"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.27.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "PrecompileTools", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "6db1a75507051bc18bfa131fbc7c3f169cc4b2f6"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.23"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "3a7c7e5c3f015415637f5debdf8a674aa2c979c4"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.1"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "04c968137612c4a5629fa531334bb81ad5680f00"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.13"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "456f610ca2fbd1c14f5fcf31c6bfadc55e7d66e0"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.43"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "Expronicon", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "50ed64cd5ad79b0bef71fdb6a11d10c3448bfef0"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.56.1"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLJacobianOperators]]
deps = ["ADTypes", "ConcreteStructs", "ConstructionBase", "DifferentiationInterface", "FastClosures", "LinearAlgebra", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "991d2a8900e687e2c693d587daa739c8fda01366"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.0"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "e39c5f217f9aca640c8e27ab21acf557a3967db5"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.10"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "25514a6f200219cd1073e4ff23a6324e4a7efe64"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.5.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ff11acffdb082493657550959d4feb4b6149e73a"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.5"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SetRounding]]
git-tree-sha1 = "d7a25e439d07a17b7cdf97eecee504c50fedf5f6"
uuid = "3cc68bcd-71a2-5612-b932-767ffbe40ab0"
version = "0.2.1"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "DiffResults", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "MaybeInplace", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "44021f3efc023be3871195d8ad98b865001a2fa1"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "1.12.3"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveTrackerExt = "Tracker"
    SimpleNonlinearSolveZygoteExt = "Zygote"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SparseConnectivityTracer]]
deps = ["ADTypes", "DocStringExtensions", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "35d346f6aa78814859f1f69cd4f41c24108afa00"
uuid = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
version = "0.6.7"

    [deps.SparseConnectivityTracer.extensions]
    SparseConnectivityTracerDataInterpolationsExt = "DataInterpolations"
    SparseConnectivityTracerLogExpFunctionsExt = "LogExpFunctions"
    SparseConnectivityTracerNNlibExt = "NNlib"
    SparseConnectivityTracerNaNMathExt = "NaNMath"
    SparseConnectivityTracerSpecialFunctionsExt = "SpecialFunctions"

    [deps.SparseConnectivityTracer.weakdeps]
    DataInterpolations = "82cc6244-b520-54b8-b5a6-8a565e85f1d0"
    LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Random", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "b906758c107b049b6b71599b9f928d9b14e5554a"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.23.0"

    [deps.SparseDiffTools.extensions]
    SparseDiffToolsEnzymeExt = "Enzyme"
    SparseDiffToolsPolyesterExt = "Polyester"
    SparseDiffToolsPolyesterForwardDiffExt = "PolyesterForwardDiff"
    SparseDiffToolsSymbolicsExt = "Symbolics"
    SparseDiffToolsZygoteExt = "Zygote"

    [deps.SparseDiffTools.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DataStructures", "DocStringExtensions", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "ccc32032d8f6790ad198c99fb8ef20d8100a0de4"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.7"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsColorsExt = "Colors"

    [deps.SparseMatrixColorings.weakdeps]
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StateSpaceSets]]
deps = ["Distances", "LinearAlgebra", "Neighborhood", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "d14f6fdb5ff0a3432831fbc69c17f066740f049c"
uuid = "40b095a5-5852-4c12-98c7-d43bf788e795"
version = "2.3.0"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools"]
git-tree-sha1 = "87d51a3ee9a4b0d2fe054bdd3fc2436258db2603"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.1.1"

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
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

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
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.SteadyStateDiffEq]]
deps = ["ConcreteStructs", "DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "Reexport", "SciMLBase"]
git-tree-sha1 = "920acf6ae36c86f23969fea1d317e040dbfccf53"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "2.4.1"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "b49030a24847e8c3e85eed3359825297ffb70dfa"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.69.1"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "f35f6ab602df8413a50c4a25ca14de821e8605fb"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.7"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a6b1675a536c5ad1a60e5a5153e1fee12eb146e3"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "82304990120934137261aa085e0d05a412fe0825"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.25.0"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "SuiteSparse_jll", "libblastrampoline_jll"]
git-tree-sha1 = "91db7ed92c66f81435fe880947171f1212936b14"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.3+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "b708e5c01d4aad69c148241bda3346fda1c550d2"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.32"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "fabf4650afe966a2ba646cabd924c3fd43577fc3"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TermInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "3927e02dc7648a45ec6aa592bcd8374094a44740"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "3.7.1"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "PrecompileTools", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "2226d810512c678d2ec9c2a9b2e227c2ebc43573"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.11.0"

    [deps.Symbolics.extensions]
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableMetadataTools]]
deps = ["DataAPI", "Dates", "TOML", "Tables", "Unitful"]
git-tree-sha1 = "c0405d3f8189bb9a9755e429c6ea2138fca7e31f"
uuid = "9ce81f87-eacc-4366-bf80-b621a3098ee2"
version = "0.1.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

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
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.TikzPictures]]
deps = ["LaTeXStrings", "Poppler_jll", "Requires", "tectonic_jll"]
git-tree-sha1 = "79e2d29b216ef24a0f4f905532b900dcf529aa06"
uuid = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
version = "3.5.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3a6f063d690135f5c1ba351412c82bae4d1402bf"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.25"

[[deps.Trapz]]
git-tree-sha1 = "79eb0ed763084a3e7de81fe1838379ac6a23b6a0"
uuid = "592b5752-818d-11e9-1e9a-2b8ca4a44cd1"
version = "2.0.3"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "be986ad9dac14888ba338c2554dcfec6939e1393"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.2.1"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulAngles]]
deps = ["Dates", "Unitful"]
git-tree-sha1 = "79875b1f2e4bf918f0702a5980816955066d9ae2"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.7.2"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "da7577e6a726959b14f7451674d00b78d10ca30f"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.2.1"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "e7f5b81c65eb858bed630fe006837b935518aca5"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.70"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "1165b0443d0eca63ac1e32b8c0eb69ed2f4f8127"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "555d1076590a6cc2fdee2ef1469451f872d8b41b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "b70c870239dc3d7bc094eb2d6be9b73d27bef280"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.44+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.tectonic_jll]]
deps = ["Artifacts", "Fontconfig_jll", "FreeType2_jll", "Graphite2_jll", "HarfBuzz_ICU_jll", "HarfBuzz_jll", "ICU_jll", "JLLWrappers", "Libdl", "OpenSSL_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "54867b00af20c70b52a1f9c00043864d8b926a21"
uuid = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"
version = "0.13.1+0"
"""

# ╔═╡ Cell order:
# ╠═fed88caa-1520-41f7-adb3-785e5c9529c6
# ╟─734b3b08-061e-4f93-8574-468d824815da
# ╟─53ad27ee-c5be-4d29-9ed5-21b6b64de42b
# ╟─a8f2765b-444f-4f3a-8a19-fbcfa6e79d2c
# ╟─c23a574d-4150-47a8-bcaf-45ae5a8c7c66
# ╟─417aad38-6928-4b13-9286-3f11efcffb99
# ╟─b842e98e-34e2-40f2-84b6-c180815c2df3
# ╟─99d0fdc9-e368-462c-a357-86f07624a52e
# ╟─6d6cbebb-bb3b-437d-9835-0a79f36857f2
# ╠═71b94af7-76ee-4228-987c-2f22e0951552
# ╟─22c37732-1cf7-4d80-a250-9fd5e4a2f88c
# ╟─2acad39a-23b7-4004-bf0e-59e73b914f01
# ╟─82991902-223e-42bb-80f4-9d6260f8a040
# ╟─c5e21675-f120-4555-84be-d99b35592f2d
# ╟─8dcc45e1-192a-4007-9d6e-6e4149b563d7
# ╟─d12b4a5f-a6ad-4c91-8795-fe8f93f5c93d
# ╟─c7eee8bb-b03d-4679-a51e-36fc629294e1
# ╟─7112c2fc-0cd6-407f-afaf-56a49e02fa57
# ╟─6b135719-c120-4b69-aed3-28647fc601bf
# ╠═40da0211-1132-412c-a769-05bd4b22b69b
# ╠═61db358f-a57d-485c-9bb3-cfd1297c0ca3
# ╟─f036dffd-8aab-4337-9f1a-a1870da284df
# ╟─fd7e6a35-72a0-4133-b43c-95ee73da5fc4
# ╠═bae7dbde-0453-472d-bf1f-5171ecbb7967
# ╠═ae2c3c1c-aa3f-4be6-b209-38b9b5b77244
# ╠═0d4877b6-0659-441e-81a9-43c59635a14f
# ╟─4d09ff8b-fd43-4b4a-ba9d-030e9d2b6cec
# ╟─700015f6-aac0-40c5-a455-8d98c1227049
# ╠═207305eb-4496-4518-a5bc-be85173314a5
# ╟─a6fdec29-3438-4ebc-af81-7a3baf0175ae
# ╟─ee7baf70-c7af-4791-8277-fadc0961e0e7
# ╟─72416f96-8a0a-4c13-97d1-2990e49c4cb0
# ╟─099a7feb-019a-4e37-9332-592c9b1dc0ed
# ╠═0baae21f-71b5-42fd-be33-31de085928d3
# ╟─a7c5650c-9bbc-44a1-8d51-d473864faae7
# ╠═49d1a7f7-2bf2-4472-94df-6247b9237ddd
# ╟─2d6fcef9-df4b-4eec-be5c-a8865c3a1b76
# ╠═568d9fe3-6716-4a9a-ba1d-b9e6fd039150
# ╟─43ee281f-1a16-445d-894d-23e0319b1fd0
# ╟─d90371c3-3359-41cb-9bc6-27484beddb3c
# ╟─1f59c4ca-37c4-4492-8670-4bfdda78bd65
# ╟─1657d7c2-87af-4079-8cf7-7a1164a95cc8
# ╟─4d09ed45-423c-4bd6-802d-59389a966d2e
# ╠═35083c71-0a7c-4b97-94ef-4d06ecbd2ee8
# ╟─76a10b7c-1135-49f1-a298-36c59bb94b37
# ╠═72794eb3-7f09-4445-906a-af9756dceebd
# ╠═1f190b07-164b-4ba9-86f8-edcf6e7ab3a9
# ╟─7df63fdf-ef41-44f2-8a55-e5c2c849029c
# ╟─0895d464-a029-410d-9e7d-89cfac2d1615
# ╟─ada0a462-e4cf-438d-b2a7-35d7280e955a
# ╠═6a59ed2e-f040-4e10-bc55-91d2c1dcc97e
# ╟─137c43b7-a541-4ca0-a6d0-0ad3579b345a
# ╠═c195ae82-0e44-4181-b863-bdd365059f4b
# ╠═14e326a5-d2bf-4873-9cf7-b57be9416da2
# ╟─0bdf9dbf-479c-46f6-bd86-50576095cba0
# ╟─6fe43e3a-2e8f-4708-a3ec-6f5a8088060e
# ╠═f863d68f-590e-4b96-8433-dc6b5177539f
# ╟─f5a983bf-ef3a-4d5d-928d-da1097b91ee8
# ╠═bd8743d6-8f21-413d-835a-e543926baa09
# ╟─9ab0a10b-8165-401a-a2b6-c9726526a906
# ╠═1b8af600-56eb-4508-bc52-aa4e138b4c7e
# ╠═80005099-7154-4306-9172-c9a168336e14
# ╠═c291700e-3a84-49a7-85d6-592cfb3b1a11
# ╟─1ec99792-905d-4e1b-a413-ef58143d3c68
# ╠═4e1140be-3e61-47bd-9f9f-6b5dfbff6ae2
# ╟─2afd459c-90e9-4105-9121-27e21bb89eeb
# ╠═4173e91c-c01a-4fa2-ae77-0408bf7d9a1b
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
