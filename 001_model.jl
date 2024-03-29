### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ b03ee99c-27f4-47df-bba5-2ea3dabdb45d
using CairoMakie, ChaosTools, DataFrames, DataFramesMeta, DelimitedFiles, DifferentialEquations, Interpolations, LinearAlgebra, PGFPlotsX, PlutoUI, QuadGK, SpecialFunctions, Symbolics, TikzPictures, Trapz, Unitful, UnitfulAstro

# ╔═╡ 08df960b-fd82-43ba-a9dc-bf5e83af587e
# ╠═╡ skip_as_script = true
#=╠═╡
TableOfContents(title="🌌 SF model", depth=4)
  ╠═╡ =#

# ╔═╡ cbd51460-8ef0-49eb-8219-14986d8421e4
# ╠═╡ skip_as_script = true
#=╠═╡
md"# Star formation model"
  ╠═╡ =#

# ╔═╡ 5814a7b3-8420-4a57-a2a2-d8c59db29a99
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Motivation

The star formation rate (SFR) is a key characteristic of galaxies. In the context of the standard cosmological model ($\Lambda$CDM), the SFR is given by a combination of processes that take place over the course of a galaxy's lifetime, such as gas cooling, star formation, chemical enrichment, and feedback from supernovae and galactic nuclei. These processes are influenced by factors like mergers, interactions, and mass accretion, which affect the amount and properties of the gas from which stars form. The density of a gas cloud is believed to be the most important factor in determining its star formation rate, although the details of this process are not yet fully understood. Observationally, the total gas density is found to be correlated to the star formation rate ([Kennicutt1998](https://doi.org/10.1086/305588)), and this correlation is even stronger with molecular gas ([Wong2002](https://doi.org/10.1086/339287), [Bigiel2008](https://doi.org/10.1088/0004-6256/136/6/2846)). The underlying reason is the intrinsic relation between molecular gas mass and SFR, which can be found at resolved scales ([Baker2021](https://doi.org/10.1093/mnras/stab3672)) and at integrated (i.e. galaxy-wide) scales across redshifts ([Baker2022](https://doi.org/10.1093/mnras/stac3413)).

Because the formation of dark matter halos and galaxies is highly non-linear, numerical simulations have become the preferred tool to investigate how galaxies form and evolve from early times up to the present. This type of simulation naturally includes mergers/interactions and continuous gas accretion. However, there are still significant uncertainties in the modeling of the baryonic component, since the physical processes that affect baryons – such as star formation, feedback, and chemical enrichment – take place at scales that are too small to be resolved directly. As a result, these processes are introduced using sub-grid physics, which involves several adjustable parameters that are not always independent of one another or constrained by observations. This can lead to inconsistencies in the predictions of different models ([Scannapieco2012](https://doi.org/10.1111/j.1365-2966.2012.20993.x), [Zhu2016](https://doi.org/10.3847/0004-637X/831/1/52), [Naab2017](https://doi.org/10.1146/annurev-astro-081913-040019)).

Because of its importance in galaxy formation, it is critical for simulations to accurately describe the star formation process at the scales that can be resolved and the associated feedback effects.
"""
  ╠═╡ =#

# ╔═╡ 8eb6540d-f5b0-45e6-883c-0cc213e67e45
# ╠═╡ skip_as_script = true
#=╠═╡
md"""

## Previous work

As the precision of $\mathrm{H_2}$ measurements increases, so too must the sophistication of simulations that model its formation and evolution. Many simulations use an equilibrium model to define the $\mathrm{H_2}$ content, which assumes that the chemistry in each volume element is in a state of equilibrium given only by local variables ([Krumholz2008](https://doi.org/10.1086/592490), [McKee2010](https://doi.org/10.1088/0004-637X/709/1/308), [Krumholz2011](https://doi.org/10.1088/0004-637X/729/1/36), and [Krumholz2013](https://doi.org/10.1093/mnras/stt1780)). However, this assumption is not always valid; as the formation and destruction of $\mathrm{H_2}$ can be influenced by a variety of factors, including temperature, density, and the presence of UV radiation.

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

As an alternative to the computationally expensive radiative transfer process, semi-analytical models (SAM) have been developed for the multiphase structure of the interstellar medium (MP ISM). Broadly, there are two ways to model the MP ISM, one considers on the physical properties of the gas (hot and cold phases), and the other on its chemical composition (hydrogen phases).

The former was pioneer by [Field1969](https://doi.org/10.1086/180324) (see [Cowie1977](https://doi.org/10.1086/154911), [McKee1977a](https://doi.org/10.1086/155350), and for a review [Cox2005](https://doi.org/10.1146/annurev.astro.43.072103.150615)), within the context of pure SAMs. The model developed by [McKee1977b](https://doi.org/10.1086/155667) was first incorporated into numerical simulation of galaxy formation by [Yepes1997](https://doi.org/10.1093/mnras/284.1.235) (Eulerian simulations) and [Hultman1999](https://ui.adsabs.harvard.edu/abs/1999A%26A...347..769H) (Lagrangian simulations). These works were later extended by [Springel2003](https://doi.org/10.1046/j.1365-8711.2003.06206.x), adding galactic winds driven by star formation as a form of feedback.

[Monaco2004](https://doi.org/10.1111/j.1365-2966.2004.07916.x) developed a SAM in a similar vein to [Springel2003](https://doi.org/10.1046/j.1365-8711.2003.06206.x), providing the theoretical foundation for MUPPI (MUlti-Phase Particle Integrator) ([Murante2010](https://doi.org/10.1111/j.1365-2966.2010.16567.x)), a sub-resolution MP ISM model that adds stellar feedback to $\texttt{GADGET2}$ ([Springel2005](https://doi.org/10.1111/j.1365-2966.2005.09655.x)). MUPPI separates a gas particle in a hot and cold phase if a set of conditions for its density and temperature are met. It then evolves those components, plus a stellar phase and an energy term (energy of the hot gas), using a set of four ODEs.

Based on [Ferrini1992](https://doi.org/10.1086/171066) and later work, [Mollá2015](https://doi.org/10.1111/j.1365-2966.2005.08782.x) developed a SAM to follow the metallicity in galaxies. These chemical evolution models (CEMs) were subsequently improved and extended in [Mollá2015](https://doi.org/10.1093/mnras/stv1102), [Mollá2016](https://doi.org/10.1093/mnras/stw1723), [Mollá2017](https://doi.org/10.1093/mnras/stx419) and [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635). The latter tracks five components: the molecular, atomic, and ionized phases of hydrogen; plus dust and stars.

We will use a sub-resolution SAM, following closely the one developed by [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635) but implemented within the hydrodynamical code $\texttt{Arepo}$ ([Springel2010](https://doi.org/10.1111/j.1365-2966.2009.15715.x)), like the way MUPPI is integrated with $\texttt{GADGET3}$ ([Springel2005](https://doi.org/10.1111/j.1365-2966.2005.09655.x)).
"""
  ╠═╡ =#

# ╔═╡ a842b24e-8d26-41ab-9de3-91632aede893
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Phases and notation

We will model the interstellar medium as a multi-phase structure made up of four components. Three are the different states of hydrogen: ionized gas with a temperature of $\sim \! 10^4 \, \mathrm{K}$, atomic gas with a temperature of $\sim \! 100 \, \mathrm{K}$, and molecular gas with a temperature of $\sim \! 10 \, \mathrm{K}$. The final component is the stars.

For every reaction that transfers mass between the phases, we will only use the dominant channel. So, even though the gas is mostly made up of hydrogen and helium, only the hydrogen reactions are incorporated into the model. The processes involved are the photoionization of atoms, the recombination of electrons with ions, the conversion of atomic hydrogen into molecular hydrogen, and the destruction of the latter through photodissociation caused by ultraviolet (UV) light. In addition, we consider the formation of ionized gas by supernovas and the influence of the molecular gas on the SFR.

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

Now, we will compute the relation between number densities ($n_j$) and the dimensionless fractions defined above ($f_j := M_j / M_\mathrm{cell}$).

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

In our model, we have done only two hypotheses until now; first that the ISM is only made up of the four components already mentioned, and second, that $V_j$ and $x_j$ are constants.

Using the values in Table 1 of [Ferrière2001](https://doi.org/10.1103/RevModPhys.73.1031) for the Milky Way,

| Component          | $n \, [\mathrm{cm^{-3}}]$ | $M \, [10^9 \, \mathrm{M_\odot}]$ |
|:------------------:|:-------------------------:|:---------------------------------:|
| Molecular gas      | $10^2 - 10^6$             | $\sim 1.3 - 2.5$                  |
| Cold atomic gas    | $20 - 50$                 | $\gtrsim 6.0$                     |
| Warm ionized gas   | $0.2 - 0.5$               | $\gtrsim 1.6$                     |

and the relation from before

$\begin{equation}
    x_j = \frac{M_j}{m_j \, n_j \, V_j} \, .
\end{equation}$

We can estimate the volume fractions by order of magnitude as

|||
|:-----:|:--------------:|
| $x_i$ | $\sim 1.0$     |
| $x_a$ | $\sim 10^{-2}$ |
| $x_m$ | $\sim 10^{-5}$ |


We will use the $x_j$ as adjustable parameters.
"""
  ╠═╡ =#

# ╔═╡ dadc3e3a-ebe7-4f13-a03b-ab988094321a
begin
	const xi = 1.0
	const xa = 1.0
	const xm = 1.0
end;

# ╔═╡ 64787011-b5b8-42be-b6e4-37ebc5138b3e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Physical relationships
"""
  ╠═╡ =#

# ╔═╡ 14c7f574-0623-4254-b8f7-97984d32351c
# ╠═╡ skip_as_script = true
#=╠═╡
TikzPictures.TikzPicture(
	L"""
		\node[box, white] (stars) {Stars};
		\node[box, white, text width=2em, above=2cm of stars] (atom) {HI};
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
	width="55em",
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

# ╔═╡ b2b23d48-9c3d-44d3-9106-745eecc9b561
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Physical processes of the model by name.
"""
  ╠═╡ =#

# ╔═╡ 047bbd39-9cf9-4bd7-b38e-16aa505b0b08
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Equations"
  ╠═╡ =#

# ╔═╡ 35e194f5-20dc-4391-b761-3696fe0bc117
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

The former is assumed to come mainly from the radiation of newborn stars, so it is given by

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_i(t)\right|_{\text{ion.}} = \eta_\text{ion} \, \psi(t) \, ,
\end{equation}$

where $\eta_\text{ion}$ is the ionized mass rate per unit of created stellar mass. All the physics of the ionization process are summarized in this parameter.

The latter, under the instantaneous recycling hypothesis, can be written as

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_i(t)\right|_{\text{recyc.}} = R \, \psi(t) \, ,
\end{equation}$

where $R$ is the mass fraction of a stellar population that is returned to the ISM, where we assumed that all the returned mass is in the form of ionized gas. Ignoring the metal enrichment is a good approximation that does not alter significantly the results.

### Atomic gas

The atomic component grows through the dissociation of hydrogen molecules and the recombination of the ionized gas with free electrons.

The former, as with the ionized gas, is given by

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_a(t)\right|_{\text{diss.}} = \eta_\text{diss} \, \psi(t) \, ,
\end{equation}$

where $\eta_\text{diss}$ is the disassociated mass rate per unit of created stellar mass.

The latter will depend on the mass of ionized gas present and the time scale of recombination ($\tau_\mathrm{rec}$), so it is given by

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_a(t)\right|_{\text{recon.}} = \frac{f_i(t)}{\tau_\mathrm{rec}(t)} \, .
\end{equation}$

### Molecular gas

The molecular component gains mass mainly by the condensation of hydrogen atoms on the surface of dust grains. This process depends on the mass of atomic gas and the characteristic time scale of condensation ($\tau_\mathrm{cond}$). We are condensing all the dust physics into this time parameter,

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_m(t)\right|_{\text{cond.}} = \frac{f_a(t)}{\tau_\mathrm{cond}(t)} \, ,
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 43eafb0f-08a8-4e53-9017-50b97ac48a52
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
	width="55em",
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

# ╔═╡ c843a9e7-c0e1-42c1-bace-c866f777232f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Physical processes of the model as mathematical expressions.
"""
  ╠═╡ =#

# ╔═╡ 70078b44-4d66-49b9-930e-74261df8be78
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
From all the above, we can write the system of four ODEs,
"""
  ╠═╡ =#

# ╔═╡ 2fe0dc4c-da44-4fc8-bef8-1fa615a0fe4a
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
	width="45em",
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

# ╔═╡ 744a9591-c7f1-496e-9bb4-47df2c8937dd
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
And from the equations we can explicitly check mass conservation,

$\begin{equation}
    \frac{\mathrm{d}}{\mathrm{d}t}(f_i + f_a + f_m + f_s) = 0 \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 34b04cf3-dabe-4364-8124-c5f3f351edb2
# ╠═╡ skip_as_script = true
#=╠═╡
md"## ODE coefficients"
  ╠═╡ =#

# ╔═╡ af69ab25-0f06-4837-ac35-acbe38a4ffb1
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Stiffness coefficient"
  ╠═╡ =#

# ╔═╡ 57ade87f-f93c-4b43-a737-a1b44f8af4fc
# ╠═╡ skip_as_script = true
#=╠═╡
md"The model, for a typical set of initial conditions, is very stiff (stiffness ratio >> 1)."
  ╠═╡ =#

# ╔═╡ 7e3183ae-b409-457a-8e35-4dc87894c580
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Lyapunov spectrum"
  ╠═╡ =#

# ╔═╡ 408721dd-bd85-47b4-a75c-6e3c188f5d64
# ╠═╡ skip_as_script = true
#=╠═╡
md"The model, for a typical set of initial conditions, is not chaotic (maximum lyapunov exponent < 0). But, we note that with a higher density the system can turn chaotic (maximum lyapunov exponent > 0)."
  ╠═╡ =#

# ╔═╡ 534a1049-8de5-4b07-abec-c5a3456627c0
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Units

We have the freedom to choose three independent units in this model. Time, mass, and length. The choice is reflected in the constants $C_\mathrm{star}$, $C_\mathrm{rec}$, and $C_\mathrm{cond}$.

Following the standard in astronomy and astrophysics, we will use $\mathrm{[T] = Myr}$, $\mathrm{[M] = mp}$ and $\mathrm{[L] = cm}$, where $\mathrm{mp}$ is the proton mass.
"""
  ╠═╡ =#

# ╔═╡ 86b692f1-0268-40f3-b4a2-d54c9828346d
begin
	const t_u = u"Myr"
	const m_u = u"mp"
	const l_u = u"cm"
end;

# ╔═╡ eaf272c7-4162-4a9a-92e3-9835c6158394
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Parameters

*  $\psi(t)$: Star formation rate.

*  $\tau_\mathrm{rec}$: Time scale of atomic gas formation from ionized gas, generally called recombination time.

*  $\tau_\mathrm{cond}$: Time scale of molecular gas formation from atomic gas, generally called condensation (or cloud formation) time.

*  $\eta_\mathrm{diss}$: Rate of molecular gas dissociation by stars per unit of created stellar mass.

*  $\eta_\mathrm{ion}$: Rate of atomic gas ionization by stars per unit of created stellar mass.

*  $R$: Mass of ionized gas produced per unit of created stellar mass.

*  $Z_\mathrm{sn}$: Mass of metals produced per unit of ionized gas created during the stellar life cycle.
"""
  ╠═╡ =#

# ╔═╡ dc6fd12b-c821-4e20-a896-25c8aab9df94
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Physical processes"
  ╠═╡ =#

# ╔═╡ ac553b12-4857-4cc1-8ea2-fe9e8863b429
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Star formation

For the star formation we will use the notation and definition commonly used in the field (for a review see [McKee2007](https://doi.org/10.1146/annurev.astro.45.051806.110602), [Krumholz2014](https://doi.org/10.1016/j.physrep.2014.02.001), and [Krumholz2019](https://doi.org/10.1146/annurev-astro-091918-104430)).

In particular, we will follow [Krumholz2005](https://doi.org/10.1086/431734), but taking into account the strong correlation between molecular hydrogen and star formation ([Bigiel2008](https://doi.org/10.1088/0004-6256/136/6/2846), [Bigiel2010](https://doi.org/10.1088/0004-6256/140/5/1194), [Wong2002](https://doi.org/10.1086/339287), [Robertson2008](https://doi.org/10.1086/587796), [Halle2013](https://doi.org/10.1051/0004-6361/201220952), [Thompson2013](https://doi.org/10.1088/0004-637x/780/2/145)). So, we end up with the simple model

$\begin{equation}
	\mathrm{SFR} := \frac{\mathrm{d}}{\mathrm{d}t} M_s = \frac{M_m}{\tau_\mathrm{star}}
\end{equation}$
where $M_m$ is the molecular mass and $\tau_\mathrm{star}$ is the characteristic timescale for star formation (defined by this very relation).

Given that we want the equations for the dimensionless fraction; the SFR enters as follows,

$\begin{equation}
	\psi := \frac{\mathrm{SFR}}{M_\mathrm{cell}} = \frac{f_m}{\tau_\mathrm{star}} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ cf488d0e-3294-45b2-b40c-fa18062c97d2
begin
	const AW = 0.0
	const MW = 1.0
end;

# ╔═╡ bc67cf22-4caa-497d-aae9-e5d1191468e2
ψ(fa, fm, τ_star) = (AW * fa + MW * fm) / τ_star;

# ╔═╡ 1d27ec35-65ca-4c94-9e8d-54d1c11e759f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Following [Krumholz2005](https://doi.org/10.1086/431734), we write the characteristic timescale for star formation as

$\begin{equation}
    \tau_\mathrm{star} =  \frac{t_\text{ff}}{\epsilon_\text{ff}} \, ,
\end{equation}$

where $\epsilon_\text{ff}$ is the star formation efficiency, (in the literature is often used $\epsilon_\text{ff} \approx 0.01$ [Krumholz2007](https://doi.org/10.1086/509101), [Krumholz2014](https://doi.org/10.1016/j.physrep.2014.02.001)), and $t_\text{ff}$ is the free-fall time which is the time for a pressure-free spherical cloud to collapse into a point due to its self-gravity.

The free-fall time can be written as

$\begin{equation}
    t_\text{ff} = \sqrt{\frac{3\pi}{32 \, G\, \rho_g}} \, ,
\end{equation}$

where $\rho_g$ is the density of the gas cloud.

There is a lot of uncertainty for the parameter $\epsilon_\text{ff}$ ([Lee2016](https://doi.org/10.3847/1538-4357/833/2/229) and [Utomo2018](https://doi.org/10.3847/2041-8213/aacf8f)). We will use $\epsilon_\text{ff} = 1.0$ because we already consider the efficiency of star formation as only molecular hydrogen is used to form stars. We note though, that this parameter has been shown to have little influence on the global properties of simulated galaxies ([Li2018](https://doi.org/10.3847/1538-4357/aac9b8) and [Brown2022](https://doi.org/10.1093/mnras/stac1164)).

With all the previous definitions, we have

$\begin{equation}
    \tau_\mathrm{star} = \frac{C_\mathrm{star}}{\sqrt{\rho_g}} \, ,
\end{equation}$
where

$\begin{equation}
    C_\mathrm{star} = \sqrt{\frac{3\pi}{32 \, G}} \, .
\end{equation}$

Given that $\rho_g$ is the density of an individual cold gas cloud, which is unresolved within a cell, we use the fact that in general $\rho_g \gg \rho_\mathrm{cell}$, to write

$\begin{align}
    \rho_g &= \rho_m + \rho_a = \rho_\mathrm{cell} \, \left( \frac{f_m}{x_m} + \frac{f_a}{x_a} \right) .
\end{align}$
"""
  ╠═╡ =#

# ╔═╡ 68732d91-805a-4663-9166-f8483213a8d2
begin
    const ϵff     = 1.0
	const C_star  = sqrt(3π / 32u"G") / ϵff
	const c_star  = ustrip(t_u * l_u^-(3/2), C_star / sqrt(m_u))
end

# ╔═╡ 27281e53-e519-4ad0-af5d-59fb0e208534
τ_star(fa, fm, ρ_cell) = c_star / sqrt(ρ_cell * (fm / xm + fa / xa));

# ╔═╡ 327fd38a-5ff6-4ac4-8d29-694272d9d46f
# ╠═╡ skip_as_script = true
#=╠═╡
let
	set_theme!(theme_black())

	f = Figure()
	ax = CairoMakie.Axis(
		f[1,1],
		xlabel=L"\rho_\mathrm{cell} \, / \, \mathrm{%$l_u^{-3}}",
		ylabel=L"\tau_\mathrm{star} \, / \, \mathrm{%$t_u}",
		title=L"f_a + f_m = 0.7",
		xlabelsize=32,
		ylabelsize=32,
		titlesize=35,
		xticklabelsize=25,
		yticklabelsize=25,
		xscale=log10,
		yscale=log10,
	)

	ρ_cell = exp10.(range(-1, 5, 30))

	for fa in range(0.1, 0.5, 3)
		label = L"f_a = %$(fa)"
		lines!(ax, ρ_cell, ρ -> τ_star(fa, 0.7 - fa, ρ); linewidth=3, label)
	end

	axislegend(; position=:rt, labelsize=25)

	f
end
  ╠═╡ =#

# ╔═╡ 9a24d3ac-a238-4eef-afc0-00fa7ef51475
# ╠═╡ skip_as_script = true
#=╠═╡
let
	ρ_cell = exp10.(range(-1, 5, 30))
	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"\rho_\mathrm{cell} \, / \, \mathrm{%$l_u^{-3}}",
	        ylabel=L"\tau_\mathrm{star} \, / \, \mathrm{%$t_u}",
	        xmode="log",
	        ymode="log",
			cycle_list_name="color list",
			"no marks",
			"thick",
	    },
	)

	for fa in range(0.1, 0.5, 3)
		push!(ax, PGFPlotsX.Plot(Coordinates(ρ_cell, τ_star.(fa, 0.7 - fa, ρ_cell))))
		push!(ax, PGFPlotsX.LegendEntry(L"f_a = %$(fa)"))
	end

	mkpath("generated_files/plots")

	pgfsave("generated_files/plots/tau_star-vs-density.pdf", ax)
end
  ╠═╡ =#

# ╔═╡ 897909e2-dcad-4ef6-9161-fd3654160dba
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

We can readily find fits for $\alpha_H(T)$ in the literature ([Seaton1959](https://doi.org/10.1093/mnras/119.2.81), [Black1981](https://doi.org/10.1093/mnras/197.3.553), [Verner1996](https://doi.org/10.1086/192284), and [Osterbrock2006](http://www.worldcat.org/oclc/60611705)). In particular, if we take $T = 10^4 \, \mathrm{K}$ for the ionized phase and use case B recombination (assuming an optically thick cloud [Nebrin2023](https://doi.org/10.3847/2515-5172/acd37a)), we get

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

# ╔═╡ 00030fd8-a9db-4903-b2ed-21a64db30588
begin
	const αH    = 2.6e-13u"cm^3 * s^-1"
	const C_rec = (m_u / αH) * (xi^2 / xa)
	const c_rec = ustrip(t_u * l_u^-3, C_rec / m_u)
end

# ╔═╡ d4f91aa3-183a-4abf-8f7a-7a05d4333e3a
τ_rec(fi, ρ_cell) = c_rec / (fi * ρ_cell);

# ╔═╡ 7e824ce1-1f82-48cc-a3c4-1acfba0e2100
# ╠═╡ skip_as_script = true
#=╠═╡
let
	set_theme!(theme_black())

	f = Figure()
	ax = CairoMakie.Axis(
		f[1,1],
		xlabel=L"\rho_\mathrm{cell} \, / \, \mathrm{%$l_u^{-3}}",
		ylabel=L"\tau_\mathrm{rec} \, / \, \mathrm{%$t_u}",
		xlabelsize=32,
		ylabelsize=32,
		titlesize=35,
		xticklabelsize=25,
		yticklabelsize=25,
		xscale=log10,
		yscale=log10,
	)

	ρ_cell = exp10.(range(-1, 5, 30))

	for fi in range(0.1, 0.9, 5)
		label = L"f_i = %$(fi)"
		lines!(ax, ρ_cell, ρ -> τ_rec(fi, ρ); linewidth=3, label)
	end

	axislegend(; position=:rt, labelsize=28)

	f
end
  ╠═╡ =#

# ╔═╡ f2cc8e6f-737c-42e3-bb81-42c50d62cf78
# ╠═╡ skip_as_script = true
#=╠═╡
let
	ρ_cell = exp10.(range(-1, 5, 30))
	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"\rho_\mathrm{cell} \, / \, \mathrm{%$l_u^{-3}}",
	        ylabel=L"\tau_\mathrm{rec} \, / \, \mathrm{%$t_u}",
	        xmode="log",
	        ymode="log",
			cycle_list_name="color list",
			"no marks",
			"thick",
	    },
	)

	for fi in range(0.1, 0.9, 5)
		push!(ax, PGFPlotsX.Plot(Coordinates(ρ_cell, τ_rec.(fi, ρ_cell))))
		push!(ax, PGFPlotsX.LegendEntry(L"f_i = %$(fi)"))
	end

	mkpath("generated_files/plots")

	pgfsave("generated_files/plots/tau_rec-vs-density.pdf", ax)
end
  ╠═╡ =#

# ╔═╡ 4a7eb24b-0874-49a3-9b08-4ffb6a7f0ce7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Condensation

From at least [Hollenbach1971a](https://doi.org/10.1086/150754) and [Hollenbach1971b](https://doi.org/10.1086/150755), we know that the rate of molecular hydrogen formation due to the condensation of atomic gas in the surface of dust grains is

$\begin{equation}
    \frac{d}{dt} n_m \biggr\rvert_\mathrm{cond.} = R_d \, n_H \, n_a  \, ,
\end{equation}$

where $R_d$ is the formation rate coefficient of $H_2$ on dust grain, $n_H$ is the hydrogen nucleus number density, and $n_a$ is the atomic hydrogen number density. $n_H$ comes from the assumption that the number density of dust grains is proportional to $n_H$ ([Hollenbach1971b](https://doi.org/10.1086/150755) Section II). 

In the literature, $n_H$ is generally taken as $n_H = n_a + 2 \, n_m$, because most studies consider cold gas clouds ($T \sim 100 \, \mathrm{K}$) dominated by molecular and atomic hydrogen gas. In contrast, we consider atomic, molecular, and ionized gas within our gas cell; so, it appears that we should be using $n_H = n_a + 2\, n_m + n_i$. But, considering that the difference should be small (the star forming cells are cold, so $f_i \ll 1$) and that the experimental values are measured with $n_H = n_a + 2\, n_m$ in mind, we will use that definition.  

We also note that the expression for $\frac{d}{dt} n_m$ is only used in equilibrium equations in most of the early works ([Hollenbach1971a](https://doi.org/10.1086/150754), [Hollenbach1971b](https://doi.org/10.1086/150755), [Jura1974](https://doi.org/10.1086/152975), [Jura1975](https://doi.org/10.1086/153545), [Black1987](https://doi.org/10.1086/165740), [Sternberg1988](https://doi.org/10.1086/166664), and [Goldshmidt1995](https://doi.org/10.1086/175168)), while first appearing in an actual differential equation that does not assume equilibrium in [Draine1996](https://doi.org/10.1086/177689).

Using the conversion factor already mentioned we can write

$\begin{align}
    \frac{d}{dt} f_m \biggr\rvert_\mathrm{cond.} &= \frac{2 \, m_p \, x_m}{\rho_\mathrm{cell}} \, R_d \, \left(n_a + 2\, n_m\right) \, n_a \\
	&= \frac{2 \, m_p \, x_m}{\rho_\mathrm{cell}} \, R_d \left( f_a \, \frac{\rho_\mathrm{cell}}{m_p \, x_a} + 2 \, f_m \, \frac{\rho_\mathrm{cell}}{2 \, m_p \, x_m} \right) \, f_a \, \frac{\rho_\mathrm{cell}}{m_p \, x_a} \\
	&= 2 \, R_d \left( \frac{f_a}{x_a} + \frac{f_m}{x_m} \right) \, f_a \, \frac{\rho_\mathrm{cell}}{m_p} \, \frac{x_m}{x_a} \\
	&= \frac{f_a}{\tau_\mathrm{cond}} \, ,
\end{align}$

where the time scale $\tau_\mathrm{cond}$ is

$\begin{equation}
    \tau_\mathrm{cond} := \frac{m_p \, x_a}{2 \, R_d \, \rho_\mathrm{cell} \, x_m \left( \dfrac{f_a}{x_a} + \dfrac{f_m}{x_m} \right)} \, .
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

Theoretical and experimental work has shown that $R_d \propto T^{1/2}$ ([Black1987](https://doi.org/10.1086/165740)) and $R_d \propto n_\mathrm{dust}$ ([Hollenbach1971b](https://doi.org/10.1086/150755)). Assuming the simplest dust model $n_\mathrm{dust} \propto Z$; we have $R_d \propto T^{1/2} \, Z$ ([Pelupessy2006](https://doi.org/10.1086/504366) and [Wolfire2008](https://doi.org/10.1086/587688)).

Following previous prescriptions ([Pelupessy2006](https://doi.org/10.1086/504366), [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55), [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x), [Mollá2017](https://doi.org/10.1093/mnras/stx419), [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635)) we will only scale $R_d$ with the metallicity; using $\sim 100\,\mathrm{K}$ for the temperature of the cold neutral gas.

Our reference value for $Z_\odot$ is ([Wolfire2008](https://doi.org/10.1086/587688))

$\begin{equation}
    R_\odot := R_d(Z = Z_\odot) = 3.5 \times 10^{-17}\mathrm{cm^3 \, s^{-1}} \, ,
\end{equation}$

so, we have

$\begin{equation}
    R_d = Z \, \frac{R_\odot}{Z_\odot} \, .
\end{equation}$

The exact value of $R_\odot$ within $1 \, \mathrm{dex}$ does not affect significantly the results. We will use an adjustable global factor, as was done with the clumping factor $C_\rho$ in [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) and [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x), to account for all the uncertainties.

The fact that the final expression for $R_d$ would give $0$ when $Z = 0$ is problematic because $\tau_\mathrm{cond}$ would diverge. The problem originates in the incorrect assumption that the only conversion channel for $\mathrm{HI} \rightarrow \mathrm{H_2}$ is the condensation on the surface of dust grains. To solve this we do the replacement $Z \rightarrow Z + Z_\mathrm{eff}$ with $Z_\mathrm{eff} = 10^{-3} \, Z_\odot$; which eliminates the divergence and takes into account the molecular formation that occurs below $10^{-3} \, Z_\odot$ ([Glover2007](https://doi.org/10.1086/519445)).

So, we finally have

$\begin{align}
    \tau_\mathrm{cond} &= \frac{m_p \, x_a \, Z_\odot}{2 \, R_\odot \, C_\rho \, (Z + Z_\mathrm{eff}) \, \rho_\mathrm{cell} \, x_m \left( \dfrac{f_a}{x_a} + \dfrac{f_m}{x_m} \right)} \\
	&= \frac{C_\mathrm{cond}}{(Z + Z_\mathrm{eff}) \, \rho_\mathrm{cell} \left( \dfrac{f_a}{x_a} + \dfrac{f_m}{x_m} \right)} \, ,
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
| $0.0127$ | Arepo                                                                 |
| $0.0196$ | [Steiger2015](https://doi.org/10.3847/0004-637X/816/1/13)             |

To keep it consistent with the Arepo codebase, we will use $Z_\odot = 0.0127$, noting that is only $35\%$ off the largest value in the list ([Steiger2015](https://doi.org/10.3847/0004-637X/816/1/13)).
"""
  ╠═╡ =#

# ╔═╡ f2a6676f-457a-476a-9ce7-c336aa9bf47f
begin
    const Rsun     = 3.5e-17u"cm^3 * s^-1"
    const Zsun     = 0.0127
	const Zeff     = 1e-3 * Zsun
	const Cρ       = 100.0
	const C_cond   = (m_u * xa * Zsun) / (2 * Rsun * Cρ * xm)
	const c_cond   = ustrip(t_u * l_u^-3, C_cond / m_u)
end

# ╔═╡ 1734df7f-1309-4ebd-a021-5f75f0bb78b2
τ_cond(fa, fm, ρ_cell, Z) = c_cond / (ρ_cell * (Z + Zeff) * (fa / xa + fm / xm));

# ╔═╡ 4f7de8a3-7f59-4a7b-8980-53390e52e0d1
# ╠═╡ skip_as_script = true
#=╠═╡
let
	set_theme!(theme_black())

	f = Figure()
	ax = CairoMakie.Axis(
		f[1,1],
		xlabel=L"\rho_\mathrm{cell} \, / \, \mathrm{%$l_u^{-3}}",
		ylabel=L"\tau_\mathrm{cond} \, / \, \mathrm{%$t_u}",
		title=L"f_a = 0.6 \,\, \mathrm{and} \,\, f_m = 0.1",
		xlabelsize=32,
		ylabelsize=32,
		titlesize=35,
		xticklabelsize=25,
		yticklabelsize=25,
		xscale=log10,
		yscale=log10,
	)

	ρ_cell = exp10.(range(-1, 5, 30))

	for Zs in range(0, 1.5, 4)
		label = L"Z \, / \, Z_\odot = %$(Zs)"
		lines!(ax, ρ_cell, ρ -> τ_cond(0.6, 0.1, ρ, Zs * Zsun); linewidth=3, label)
	end

	axislegend(; position=:rt, labelsize=28)

	f
end
  ╠═╡ =#

# ╔═╡ 7ad23ea4-9887-4de5-8a5c-c37ebef736b8
# ╠═╡ skip_as_script = true
#=╠═╡
let
	ρ_cell = exp10.(range(-1, 5, 30))
	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"\rho_\mathrm{cell} \, / \, \mathrm{%$l_u^{-3}}",
	        ylabel=L"\tau_\mathrm{cond} \, / \, \mathrm{%$t_u}",
	        xmode="log",
	        ymode="log",
			cycle_list_name="color list",
			"no marks",
			"thick",
	    },
	)

	for Zs in range(0, 1.5, 4)
		push!(
			ax,
			PGFPlotsX.Plot(Coordinates(ρ_cell, τ_cond.(0.6, 0.1, ρ_cell, Zs * Zsun))),
		)
		push!(ax, PGFPlotsX.LegendEntry(L"Z \, / \, Z_\odot = %$(Zs)"))
	end

	mkpath("generated_files/plots")

	pgfsave("generated_files/plots/tau_cond-vs-density.pdf", ax)
end
  ╠═╡ =#

# ╔═╡ 3767c7f9-a0bc-467a-a20a-5e5a266111c7
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
    \eta_\mathrm{ion} = \frac{\dot{M}_\mathrm{ion}}{SFR}  \, .
\end{equation}$

The mass rates, $\dot{M}_\mathrm{diss}$ and $\dot{M}_\mathrm{ion}$, can be computed from the photon production rate in the corresponding energy range,

$\begin{align}
    \dot{M}_\mathrm{diss} &= \dot{N}_\mathrm{ion} \, f_\mathrm{diss} \, , \\
    \dot{M}_\mathrm{ion} &= \dot{N}_\mathrm{diss} \, f_\mathrm{ion} \, ,
\end{align}$

where $\dot{N}_\mathrm{diss}$ is the number of photodissociating photons produced per unit time (in the Lyman–Werner band, $912\,\mathrm{Å}$ to $1107\,\mathrm{Å}$), $\dot{N}_\mathrm{ion}$ the number of ionizing photons produced per unit time (between $0$ and $912\,Å$), and $f_\mathrm{diss}$ and $f_\mathrm{ion}$ are the unit conversion factors (proton mass into solar mass). The factors $f_\mathrm{diss}$ and $f_\mathrm{ion}$ allow us to consider that the reaction may not be $100\%$ efficient too.

For the ionization reaction, each photon will produce one proton and we assume $100\%$ efficiency.
"""
  ╠═╡ =#

# ╔═╡ f65d84cd-ab5f-4270-98ba-568792d1fec1
const f_ion = ustrip(1.0u"mp" |> u"Msun")

# ╔═╡ 34faac11-85a2-44dc-bd8d-1a71656fccf4
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
For the molecular dissociation reaction, we will consider the numerical factor given by [Draine1996](https://doi.org/10.1086/177689). It is shown that dust grains may absorb up to $\sim \! 60$ percent of the photons capable of dissociating hydrogen molecules and a large fraction of the remaining photons excite different rotational and vibrational states, reducing their dissociation probability to $\sim \! 15$ percent. So, we end up with an efficiency factor of $0.4 \times 0.15 = 0.06$.
$f_\mathrm{d}$ has an extra factor of two because each photon contributes with two protons (from the dissociated molecule) to the atomic gas.
"""
  ╠═╡ =#

# ╔═╡ f8b02d00-ff30-480e-b5eb-e150e4678c95
const f_diss = 0.4 * 0.15 * 2.0 * ustrip(1.0u"mp" |> u"Msun")

# ╔═╡ 44c88ad8-a8c3-45e3-9a56-be3ce5bf66fa
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The number of photons can be computed from $Q(t', Z)$ which is defined as the number of photons produced by a stellar population of one solar mass, of age $t'$, and metallicity $Z$, per unit of time. So, we have the relation

$\begin{equation}
    \dot{N}(t) = \int_0^t \mathrm{SFR}(t - t') \, Q(t', Z(t - t')) \mathrm{d}t' \, ,
\end{equation}$

where $\mathrm{SFR}(t - t')$ is the instantaneous SFR at the moment of birth of the stellar population of current age $t'$ (in units of solar mass per unit of time), and $Z = Z(t - t')$ is defined as the metallicity at that moment. Because most of the contribution to the integral comes from young blue stars that die in the first $10 \ \mathrm{to} \ 100 \, \mathrm{Myr}$, it is possible to approximate

$\begin{align}
    \mathrm{SFR}(t - t') &\approx \mathrm{SFR}(t) \, , \\
	Z(t - t') &\approx Z(t) \, .
\end{align}$

So, we end up with

$\begin{equation}
    \eta = f \, \frac{\dot{N}}{\mathrm{SFR}} = f \, \int_0^t Q(t', Z) \mathrm{d}t'\, ,
\end{equation}$

The value of $Q$ can be calculated using

$\begin{equation}
    Q(t', Z) = \int_{\lambda_1}^{\lambda_2} \frac{L_\lambda(t', Z)}{E_\lambda} \mathrm{d}\lambda = \int_{\lambda_1}^{\lambda_2} \frac{\lambda \, L_\lambda(t', Z)}{h \, c} \mathrm{d}\lambda \, ,
\end{equation}$

where $L_\lambda(t', Z)$ is the luminosity per unit of wavelength of a stellar population of one solar mass, age $t'$, and metallicity $Z$; and $E_\lambda = \lambda / (h \, c)$ is the energy of a photon of wavelength $\lambda$. So, integrating between the wavelength of interest, we get the number of photons produced per unit of time for a stellar population of the given characteristics.

The luminosity will not only depend on the age and metallicity of the population, but on the IMF (initial mass function) too, so in principle for each IMF, $t'$, and $Z$ we have a function of luminosity versus $\lambda$.

Using the values from PopStar by [Mollá2009](https://doi.org/10.1111/j.1365-2966.2009.15160.x) we compute a table of $Q$ for six IMFs ([Salpeter1955](https://doi.org/10.1086/145971) in two mass ranges: $0.85\,\mathrm{M}_\odot$ to $120\,\mathrm{M}_\odot$ and $0.15\,\mathrm{M}_\odot$ to $100\,\mathrm{M}_\odot$, [Kroupa2001](https://doi.org/10.1046/j.1365-8711.2001.04022.x), [Ferrini1990](https://ui.adsabs.harvard.edu/abs/1990A%26A...231..391F), and [Chabrier2003](https://doi.org/10.1086/374879)), six metallicities (0.0001, 0.0004, 0.004, 0.008, 0.02, and 0.05) and for ages between $0.1\,\mathrm{Myr}$, and $15\,\mathrm{Gyr}$.
"""
  ╠═╡ =#

# ╔═╡ 448e1dee-4628-4c14-9d6f-dc165b2e826e
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

# ╔═╡ 5311b7cc-7199-45c8-b5e7-20d3ceb191b7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
In what follows we will use the IMF by [Chabrier2003](https://doi.org/10.1086/374879)
"""
  ╠═╡ =#

# ╔═╡ 3e637368-6bdb-4d22-9a4a-df23c6682c2f
begin
	const IMF      = "Chabrier2003"             # Initial mass function
	const Q_df     = Q_by_imf[IMF]
	const Q_ages   = unique(Q_df[!, :log_age])  # List of stellar ages
    const Q_metals = unique(Q_df[!, :Zmet])     # List of metallicities
end;

# ╔═╡ ef65a096-cc2a-4ce6-a06b-8671c99ca777
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
			η_diss[i, j] = uconvert(Unitful.NoUnits, trapz(ages, Q_diss) * f_diss)

			q_ion       = sub_df[!, :Q_ion]
			Q_ion       = [q_ion[1], q_ion...]
	        η_ion[i, j] = uconvert(Unitful.NoUnits, trapz(ages, Q_ion) * f_ion)

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

# ╔═╡ a0294888-90cf-4e5b-a4b8-ce2c63bdae7a
# ╠═╡ skip_as_script = true
#=╠═╡
let
	set_theme!(theme_black())

	f = Figure()
	ax = CairoMakie.Axis(
		f[1,1],
		xlabel=L"\mathrm{stellar \,\, age \, / \, Myr}",
		ylabel=L"\eta_\mathrm{diss}",
		title=L"\mathrm{IMF}: \,\, \mathrm{%$IMF}",
		xlabelsize=32,
		ylabelsize=32,
		xticklabelsize=25,
		yticklabelsize=25,
		titlesize=30,
		xscale=log10,
	)

	ages = exp10.(range(-1, 3, 100))
	ηd(age, Zs) = photodissociation_efficiency(age, Zs * Zsun)[1]

	for Zs in range(0, 1.5, 4)
		label = L"Z \, / \, Z_\odot = %$(Zs)"
		lines!(ax, ages, x -> ηd(x, Zs); linewidth=3, label)
	end

	axislegend(ax; position=:rb, labelsize=25)

	f
end
  ╠═╡ =#

# ╔═╡ 303da7a3-e574-4f7e-9acc-83ed83a09f69
# ╠═╡ skip_as_script = true
#=╠═╡
let
	ages            = exp10.(range(-1, 3, 100))
	η_diss(age, Zs) = photodissociation_efficiency(age, Zs * Zsun)[1]

	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"\mathrm{stellar \,\, age \, / \, Myr}",
	        ylabel=L"\eta_\mathrm{diss}",
	        xmode="log",
			cycle_list_name="color list",
			"no marks",
			"thick",
			legend_style={at = Coordinate(0.75, 0.4), anchor = "north"},
	    },
	)

	for Zs in range(0, 1.5, 4)
		push!(ax, PGFPlotsX.Plot(Coordinates(ages, η_diss.(ages, Zs))))
		push!(ax, PGFPlotsX.LegendEntry(L"Z \, / \, Z_\odot = %$(Zs)"))
	end

	mkpath("generated_files/plots")

	pgfsave("generated_files/plots/eta_diss-vs-stellar_age.pdf", ax)
end
  ╠═╡ =#

# ╔═╡ 994f97fb-1c30-4825-9b29-35fe4ade8fb3
# ╠═╡ skip_as_script = true
#=╠═╡
let
	set_theme!(theme_black())

	f = Figure()
	ax = CairoMakie.Axis(
		f[1,1],
		xlabel=L"\mathrm{stellar \,\, age \, / \, Myr}",
		ylabel=L"\eta_\mathrm{ion}",
		title=L"\mathrm{IMF}: \,\, \mathrm{%$IMF}",
		xlabelsize=32,
		ylabelsize=32,
		xticklabelsize=25,
		yticklabelsize=25,
		titlesize=30,
		xscale=log10,
	)

	ages = exp10.(range(-1, 3, 100))
	η_ion(age, Zs) = photodissociation_efficiency(age, Zs * Zsun)[2]

	for Zs in range(0, 1.5, 4)
		label = L"Z \, / \, Z_\odot = %$(Zs)"
		lines!(ax, ages, x -> η_ion(x, Zs); linewidth=3, label)
	end

	axislegend(ax; position=:rb, labelsize=25)

	f
end
  ╠═╡ =#

# ╔═╡ 493e649f-58ce-4b90-9b61-f9f8a6146ed9
# ╠═╡ skip_as_script = true
#=╠═╡
let
	ages           = exp10.(range(-1, 3, 100))
	η_ion(age, Zs) = photodissociation_efficiency(age, Zs * Zsun)[2]

	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"\mathrm{stellar \,\, age \, / \, Myr}",
	        ylabel=L"\eta_\mathrm{ion}",
	        xmode="log",
			cycle_list_name="color list",
			"no marks",
			"thick",
			legend_style={at = Coordinate(0.75, 0.4), anchor = "north"},
	    },
	)

	for Zs in range(0, 1.5, 4)
		push!(ax, PGFPlotsX.Plot(Coordinates(ages, η_ion.(ages, Zs))))
		push!(ax, PGFPlotsX.LegendEntry(L"Z \, / \, Z_\odot = %$(Zs)"))
	end

	mkpath("generated_files/plots")

	pgfsave("generated_files/plots/eta_ion-vs-stellar_age.pdf", ax)
end
  ╠═╡ =#

# ╔═╡ 533b3cd0-c1f6-4ecd-b196-4ed35bf77135
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Mass recycling

There are two mass recycling parameters: $R$ which is defined as the mass fraction of a stellar population that is returned to the ISM under the instantaneous recycling hypothesis (stars under a certain mass live forever and stars above that mass die instantly), and $Z_\mathrm{sn}$ which is the fraction of the returned gas that is composed of metals (the rest is assumed to be ionized gas).

Notice that the instantaneous recycling hypothesis can be avoided by considering the lifetimes of the stars (using empirical relations) as it is done in [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635) (sections 2.2.2 and 2.2.3). This would effectively make $R$ time-dependent (the integrals below would have to be computed at each evaluation of the equations) increasing significantly the computational cost, so we avoid this and assume a constant $R$ during the ODEs integration time scales.

A stellar yield model gives the amount (as a fraction of the stellar mass) of each element that is returned to the ISM by stars with masses between $m$ and $m + \mathrm{d}m$, so $R$ can be computed as

$\begin{equation}
	R = \dfrac{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} (m - m_\mathrm{rem}(m)) \, \phi(m) \, \mathrm{d}m}{\int_{m_\mathrm{low}}^{m_\mathrm{high}} m \, \phi(m) \, \mathrm{d}m} \, ,
\end{equation}$

where $\phi(m)$ is the ISM, $m_\mathrm{low}$ and $m_\mathrm{high}$ are the extremes in the mass range of the ISM, $m_\mathrm{ir}$ is the mass limit for the instantaneous recycling hypothesis, and $m_\mathrm{rem}(m)$ is the remnant stellar mass given by the yield model ([Pipino2014](https://doi.org/10.1093/mnras/stu579) and [Ascasibar2015](https://doi.org/10.1093/mnras/stv098)).

Notice that the denominator in the expression for $R$ is the total mass of the stellar population modeled by $\phi(m)$, so it is just a normalization, which is needed because the IMF is in general defined except for a global constant.

Using the same notation, we can calculate $Z_\mathrm{sn}$ as

$\begin{equation}
	Z_\mathrm{sn} = \dfrac{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} m \, f_Z \, \phi(m) \, \mathrm{d}m}{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} (m - m_\mathrm{rem}(m)) \, \phi(m) \, \mathrm{d}m} \, ,
\end{equation}$

where $f_Z$ is the fraction of the stellar mass that is returned to the ISM as metals.

Some traditional choices for the masses are $m_\mathrm{ir} = 8 \, M_\odot$, $m_\mathrm{low} = 0.08 \, M_\odot$ (the limit for hydrogen fusion), and $m_\mathrm{high} = 100 \, M_\odot$ (the order of magnitude for the upper limit of validity of the yield models given by the observational limitations).

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

# ╔═╡ be85ba3b-5439-4cf3-bb14-d24d61a283c3
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

# ╔═╡ b3a260b6-eb31-43a0-9fd6-60a507984319
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Initial mass functions (IMF)

The initial mass function $\phi(m)$ gives the number of stars between masses $m$ and $m + \mathrm{d}m$. For a given population of mass $M$, we have

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

# ╔═╡ 5ba3a0c1-6107-45a1-9b1d-5c323b9a7145
######################################################################################
# Implementation of several IMF
######################################################################################

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

# ╔═╡ 7cbf5573-032e-4ddd-9575-f387a577c93e
begin
	const Y_MODEL = "Portinari1998"  # Stellar yield model
	const M_LOW   = 0.1u"Msun"       # Lower mass limit for the ISM
	const M_HIGH  = 100.0u"Msun"     # Upper mass limit for the ISM
	const M_IR    = 8.0u"Msun"       # Mass limit for the IR hypothesis
end;

# ╔═╡ 1b044783-0f5f-4321-abda-35e5b7ae67c4
######################################################################################
# Compute R and Zsn for a given stellar metallicity
#
# Z: Metallicity [dimensionless]
######################################################################################

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

# ╔═╡ 1d0a66d0-5791-4dc9-a5d4-7882b0e91767
# ╠═╡ skip_as_script = true
#=╠═╡
let
	set_theme!(theme_black())

	f = Figure()
	ax = CairoMakie.Axis(
		f[1,1],
		xlabel=L"Z \, / \, Z_\odot",
		ylabel=L"R",
		xlabelsize=32,
		ylabelsize=32,
		xticklabelsize=25,
		yticklabelsize=25,
	)

	metalicities = exp10.(range(-1, 0.2, 100))
	R(Zs) = recycled_fractions(Zs * Zsun)[1]

	lines!(ax, metalicities, R; linewidth=3)

	f
end
  ╠═╡ =#

# ╔═╡ bb4ee2dc-c069-415e-bac5-0130950f3941
# ╠═╡ skip_as_script = true
#=╠═╡
let
	metalicities = exp10.(range(-1, 0.2, 100))
	R(Zs) = recycled_fractions(Zs * Zsun)[1]

	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"Z \, / \, Z_\odot",
	        ylabel=L"R",
			"no marks",
			"thick",
	    },
	)

	push!(ax, PGFPlotsX.Plot(Coordinates(metalicities, R.(metalicities))))

	mkpath("generated_files/plots")

	pgfsave("generated_files/plots/R-vs-metalicity.pdf", ax)
end
  ╠═╡ =#

# ╔═╡ 041916ac-4cb0-4630-a227-043fae52264d
# ╠═╡ skip_as_script = true
#=╠═╡
let
	set_theme!(theme_black())

	f = Figure()
	ax = CairoMakie.Axis(
		f[1,1],
		xlabel=L"Z \, / \, Z_\odot",
		ylabel=L"Z_\mathrm{sn}",
		xlabelsize=32,
		ylabelsize=32,
		xticklabelsize=25,
		yticklabelsize=25,
	)

	metalicities = exp10.(range(-1, 0.2, 100))
    Zsn(Zs) = recycled_fractions(Zs * Zsun)[2]

	lines!(ax, metalicities, Zsn; linewidth=3)

	f
end
  ╠═╡ =#

# ╔═╡ 57ea5e31-d156-4df5-bb77-0bc01b3559af
# ╠═╡ skip_as_script = true
#=╠═╡
let
	metalicities = exp10.(range(-1, 0.2, 100))
    Zsn(Zs) = recycled_fractions(Zs * Zsun)[2]

	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"Z \, / \, Z_\odot",
	        ylabel=L"Z_\mathrm{sn}",
			"no marks",
			"thick",
	    },
	)

	push!(ax, PGFPlotsX.Plot(Coordinates(metalicities, Zsn.(metalicities))))

	mkpath("generated_files/plots")

	pgfsave("generated_files/plots/Zsn-vs-metalicity.pdf", ax)
end
  ╠═╡ =#

# ╔═╡ 9666bdc8-cbc0-4757-9bd8-a76477c252eb
# ╠═╡ skip_as_script = true
#=╠═╡
md"# Implementation"
  ╠═╡ =#

# ╔═╡ ca9a233b-d3ca-4a76-a3d8-f29884ac9484
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Constants"
  ╠═╡ =#

# ╔═╡ d8bee772-3979-42cd-9e38-8df0925b4e6b
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

# ╔═╡ e2e4ae4f-dcdc-4999-88f2-853378be859a
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Equations"
  ╠═╡ =#

# ╔═╡ 177f8253-6c35-495b-9119-ce5e8e15cba8
######################################################################################
# System of ODEs, where each equation has units of time⁻¹, and
#
# Ionized gas fraction:    fᵢ(t) = Mᵢ(t) / M_cell --> y[1]
# Atomic gas fraction:     fₐ(t) = Mₐ(t) / M_cell --> y[2]
# Molecular gas fraction:  fₘ(t) = Mₘ(t) / M_cell --> y[3]
# Stellar fraction:        fₛ(t) = Mₛ(t) / M_cell --> y[4]
######################################################################################

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
    cloud_formation = fa / τ_cond(fa, fm, ρ_cell, Z)
	sfr             = ψ(fa, fm, τ_star(fa, fm, ρ_cell))

    # ODE system
	dydt[1] = -recombination + (η_ion + R) * sfr
    dydt[2] = -cloud_formation + recombination + (η_diss - η_ion) * sfr
    dydt[3] = cloud_formation - (1 + η_diss) * sfr
    dydt[4] = (1 - R) * sfr

end;

# ╔═╡ d5446d26-4b59-46fa-a1f0-8374d5e05194
######################################################################################
# Compute the Lyapunov spectrum
#
# ic:         Initial conditions, [fi(0), fa(0), fm(0), fs(0)]
# base_parms: Parameters, (ρ [cm⁻³], Z [dimensionless], it [Myr])
######################################################################################

function lyapunov(
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

# ╔═╡ 7a744b2e-56ee-4162-a23d-dca9d0657608
# ╠═╡ skip_as_script = true
#=╠═╡
lyapunov(
	[0.8, 0.2, 0.0, 0.0],     # Initial conditions, [fi(0), fa(0), fm(0), fs(0)]
	(100.0, 0.5*Zsun, 10.0),  # Parameters, (ρ [cm⁻³], Z [dimensionless], it [Myr])
)
  ╠═╡ =#

# ╔═╡ d0ecada4-13d1-4977-b66e-d992a1adae15
# ╠═╡ skip_as_script = true
#=╠═╡
lyapunov(
	[0.8, 0.2, 0.0, 0.0],      # Initial conditions, [fi(0), fa(0), fm(0), fs(0)]
	(1000.0, 0.5*Zsun, 10.0),  # Parameters, (ρ [cm⁻³], Z [dimensionless], it [Myr])
)
  ╠═╡ =#

# ╔═╡ b3969810-ab25-4e91-ad5a-80560b80977e
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Jacobian"
  ╠═╡ =#

# ╔═╡ 2620d8a6-030d-4a6f-911c-6552072ff7a1
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

# ╔═╡ b76d4669-26dc-48cb-930f-5e40dd40a9f1
const JACOBIAN_FUNCTION = construct_jacobian(system!);

# ╔═╡ 69b8d934-c031-413b-9c86-3fbd64be5a4a
######################################################################################
# Evaluate the Jacobian
#
# J:          Matrix to save the results, it must have size N_EQU × N_EQU
# ic:         Initial condition, [fi(0), fa(0), fm(0), fs(0)]
# parameters: Parameters for the ODEs, [ρ, Z, η_diss, η_ion, R]
# t:          Unused variable to comply with with the DifferentialEquations.jl API
######################################################################################

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

# ╔═╡ 806db782-1734-4112-ab7b-84e03f4c342d
######################################################################################
# Compute the stiffness ratio
#
# ic:         Initial conditions, [fi(0), fa(0), fm(0), fs(0)]
# base_parms: Parameters, (ρ [cm⁻³], Z [dimensionless], it [Myr])
######################################################################################

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

# ╔═╡ 1743b6dc-0a4e-4a02-90fe-3bc47833421a
# ╠═╡ skip_as_script = true
#=╠═╡
stiffness_ratio(
	[0.8, 0.2, 0.0, 0.0],     # Initial conditions, [fi(0), fa(0), fm(0), fs(0)]
	(100.0, 0.5*Zsun, 10.0),  # Parameters, (ρ [cm⁻³], Z [dimensionless], it [Myr])
)
  ╠═╡ =#

# ╔═╡ 35ac9289-ba53-453f-9d9e-ef3499949a98
# ╠═╡ skip_as_script = true
#=╠═╡
md"## ODE Function"
  ╠═╡ =#

# ╔═╡ 64e7e6aa-4265-4de2-a9c6-474d125b45cc
ode_function = ODEFunction{true}(
	system!;
	jac=jacobian!,
	tgrad=(dt, ic, p, t) -> nothing,
);

# ╔═╡ ff37c5e7-bb7e-4e3a-a16c-b8b7c4072528
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Density PDF"
  ╠═╡ =#

# ╔═╡ ccfafe47-841a-4c8a-b609-de34b453f7ee
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Parameters"
  ╠═╡ =#

# ╔═╡ 5e402cd6-ebcb-4641-8ce8-2f7c6c37a5ea
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

# ╔═╡ e3a6d4ab-fb1c-4fa7-b4ca-50552b720ab5
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Mass fractions"
  ╠═╡ =#

# ╔═╡ 0c4fb275-240f-4689-857c-fb30120af9a5
######################################################################################
# Compute the density PDF mass fractions
#
# params:  Parameters for the density PDF
# log_var: Selects which variable will be used and how the function will be divided
#   log_var == true:  s = ln(ρ/ρ₀) and logarithmic divisions
#   log_var == false: f = ρ/ρ₀ and linear divisions
######################################################################################

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

# ╔═╡ e050e2e9-30c1-4dc1-a57b-4f86821e3965
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Density PDF by Burkhart (2018)"
  ╠═╡ =#

# ╔═╡ d8be4069-e589-495f-a6d8-0b1fcc048dd2
######################################################################################
# Density PDF acording to Burkhart (2018)
# https://doi.org/10.3847/1538-4357/aad002
######################################################################################

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

# ╔═╡ 16b56ce2-6b02-4d53-a7a7-7cd7f96f26c5
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

# ╔═╡ 4607856c-7472-4131-a2ee-29f7150f5cb4
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Integration"
  ╠═╡ =#

# ╔═╡ bbb7263a-91e4-4a23-9e5f-416b6b7fcf6e
######################################################################################
# Solve the system of ODEs
#
# ic:          Initial conditions, [fi(0), fa(0), fm(0), fs(0)]
# base_params: Parameters for the ODEs, [ρ_cell, Z]
# tspan:       Integration span, (ti, tf) [Myr]
# times:       Times at which the solution will be returned [Myr]
# args:        Positional arguments for the solver of DifferentialEquations.jl
# kwargs:      Keyword arguments for the solver of DifferentialEquations.jl
######################################################################################

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
	η_diss, η_ion = photodissociation_efficiency(tspan[2], Z)
	R, _          = recycled_fractions(Z)

	parameters = [ρ_cell, Z, η_diss, η_ion, R]

    sol = solve(ODEProblem(
		ode_function,
		ic,
		tspan,
		parameters,
	), args...; kwargs...)

    return sol(times).u

end;

# ╔═╡ 0cac7dbf-31bb-459b-9291-8092cfb49354
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Plots"
  ╠═╡ =#

# ╔═╡ e65aef30-bf57-4d7d-b878-276ea6c0ef4a
# ╠═╡ skip_as_script = true
#=╠═╡
let
	mkpath("generated_files/plots")

	ic          = [0.15, 0.85, 0.0, 0.0]
	ρ_cell_list = exp10.(range(-1, 5, 100))
	base_params = [[ρ_cell, 0.01 * Zsun] for ρ_cell in ρ_cell_list]
	it          = 10.0

	# Star fraction
	idx = phase_name_to_index["stellar"]
	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"\rho_\mathrm{cell} \, / \, \mathrm{cm^{-3}}",
	        ylabel=L"\mathrm{M_s / M_\mathrm{cell}}",
	        xmode="log",
			ymode="log",
			cycle_list_name="color list",
			"no marks",
			"thick",
	    },
	)

	s_f = [
		integrate_model(ic, base_param, (0.0, it))[end][idx] for
		base_param in base_params
	]
	push!(ax, PGFPlotsX.Plot(Coordinates(ρ_cell_list, s_f)))

	pgfsave("generated_files/plots/star_fraction-vs-density.pdf", ax)

	# Molecular fraction
	idx         = phase_name_to_index["molecular"]

	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"\rho_\mathrm{cell} \, / \, \mathrm{cm^{-3}}",
	        ylabel=L"\mathrm{M_m / M_\mathrm{cell}}",
	        xmode="log",
			ymode="log",
			cycle_list_name="color list",
			"no marks",
			"thick",
	    },
	)

	m_f = [
		integrate_model(ic, base_param, (0.0, it))[end][idx] for
		base_param in base_params
	]
	push!(ax, PGFPlotsX.Plot(Coordinates(ρ_cell_list, m_f)))

	pgfsave("generated_files/plots/molecular_fraction-vs-density.pdf", ax)

	# Atomic fraction
	idx = phase_name_to_index["atomic"]

	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"\rho_\mathrm{cell} \, / \, \mathrm{cm^{-3}}",
	        ylabel=L"\mathrm{M_a / M_\mathrm{cell}}",
	        xmode="log",
			ymode="log",
			cycle_list_name="color list",
			"no marks",
			"thick",
	    },
	)

	a_f = [
		integrate_model(ic, base_param, (0.0, it))[end][idx] for
		base_param in base_params
	]
	push!(ax, PGFPlotsX.Plot(Coordinates(ρ_cell_list, a_f)))

	pgfsave("generated_files/plots/atomic_fraction-vs-density.pdf", ax)

	# Ionized fraction
	idx = phase_name_to_index["ionized"]

	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"\rho_\mathrm{cell} \, / \, \mathrm{cm^{-3}}",
	        ylabel=L"\mathrm{M_i / M_\mathrm{cell}}",
	        xmode="log",
			ymode="log",
			cycle_list_name="color list",
			"no marks",
			"thick",
	    },
	)

	i_f = [
		integrate_model(ic, base_param, (0.0, it))[end][idx] for
		base_param in base_params
	]
	push!(ax, PGFPlotsX.Plot(Coordinates(ρ_cell_list, i_f)))

	pgfsave("generated_files/plots/ionized_fraction-vs-density.pdf", ax)
end
  ╠═╡ =#

# ╔═╡ 729bb0e3-3d60-4ac6-803e-2568a6970d44
# ╠═╡ skip_as_script = true
#=╠═╡
let
	mkpath("generated_files/plots")

	ic          = [0.5, 0.5, 0.0, 0.0]
	Z_list      = exp10.(range(-2, 0.0, 100))
	base_params = [[1.0e2, Z * Zsun] for Z in Z_list]
	it          = 10.0

    # Star fraction
	idx = phase_name_to_index["stellar"]
	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"Z \, / \, Z_\odot",
	        ylabel=L"\mathrm{M_s / M_\mathrm{cell}}",
	        xmode="log",
			ymode="log",
			cycle_list_name="color list",
			"no marks",
			"thick",
	    },
	)

	s_f = [
		integrate_model(ic, base_param, (0.0, it))[end][idx] for
		base_param in base_params
	]
	push!(ax, PGFPlotsX.Plot(Coordinates(Z_list, s_f)))

	pgfsave("generated_files/plots/star_fraction-vs-metalicity.pdf", ax)

	# Molecular fraction
	idx = phase_name_to_index["molecular"]
	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"Z \, / \, Z_\odot",
	        ylabel=L"\mathrm{M_m / M_\mathrm{cell}}",
	        xmode="log",
			ymode="log",
			cycle_list_name="color list",
			"no marks",
			"thick",
	    },
	)

	m_f = [
		integrate_model(ic, base_param, (0.0, it))[end][idx] for
		base_param in base_params
	]
	push!(ax, PGFPlotsX.Plot(Coordinates(Z_list, m_f)))

	pgfsave("generated_files/plots/molecular_fraction-vs-metalicity.pdf", ax)

	# Atomic fraction
	idx = phase_name_to_index["atomic"]
	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"Z \, / \, Z_\odot",
	        ylabel=L"\mathrm{M_a / M_\mathrm{cell}}",
	        xmode="log",
			ymode="log",
			cycle_list_name="color list",
			"no marks",
			"thick",
	    },
	)

	a_f = [
		integrate_model(ic, base_param, (0.0, it))[end][idx] for
		base_param in base_params
	]
	push!(ax, PGFPlotsX.Plot(Coordinates(Z_list, a_f)))

	pgfsave("generated_files/plots/atomic_fraction-vs-metalicity.pdf", ax)

	# Ionized fraction
	idx = phase_name_to_index["ionized"]
	@pgf ax = PGFPlotsX.Axis(
	    {
	        xlabel=L"Z \, / \, Z_\odot",
	        ylabel=L"\mathrm{M_i / M_\mathrm{cell}}",
	        xmode="log",
			ymode="log",
			cycle_list_name="color list",
			"no marks",
			"thick",
	    },
	)

	i_f = [
		integrate_model(ic, base_param, (0.0, it))[end][idx] for
		base_param in base_params
	]
	push!(ax, PGFPlotsX.Plot(Coordinates(Z_list, i_f)))

	pgfsave("generated_files/plots/ionized_fraction-vs-metalicity.pdf", ax)

end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
ChaosTools = "608a59af-f2a3-5ad4-90b4-758bdf3122a7"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataFramesMeta = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PGFPlotsX = "8314cec4-20b6-5062-9cdb-752b83310925"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
Trapz = "592b5752-818d-11e9-1e9a-2b8ca4a44cd1"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"

[compat]
CairoMakie = "~0.10.6"
ChaosTools = "~3.1.2"
DataFrames = "~1.6.0"
DataFramesMeta = "~0.14.1"
DelimitedFiles = "~1.9.1"
DifferentialEquations = "~7.10.0"
Interpolations = "~0.15.1"
PGFPlotsX = "~1.6.0"
PlutoUI = "~0.7.51"
QuadGK = "~2.8.1"
SpecialFunctions = "~2.3.1"
Symbolics = "~5.5.0"
TikzPictures = "~3.5.0"
Trapz = "~2.0.3"
Unitful = "~1.12.3"
UnitfulAstro = "~1.2.0"

[extras]
CPUSummary = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "86fb73f929a6aea624c57bd72c41e0e07a694036"

[[deps.ADTypes]]
git-tree-sha1 = "41c37aa88889c171f1300ceac1313c06e891d245"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.6"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Preferences", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "c3c29bf6363b3ac3e421dc8b2ba8e33bdacbd245"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.32.5"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractLattices]]
git-tree-sha1 = "222ee9e50b98f51b5d78feb93dd928880df35f06"
uuid = "398f06c4-4d28-53ec-89ca-5b2656b7603d"
version = "0.3.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "793501dcd3fa7ce8d375a2c878dca2296232686e"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cde29ddf7e5726c9fb511f340244ea3481267608"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.7.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "bbec08a37f8722786d87bedf84eae19c020c4efa"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.7.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "b08a4043e1c14096ef8efe4dd97e07de5cacf240"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.4.5"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "0da671c730d79b8f9a88a391556ec695ea921040"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.0.2"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "0b816941273b5b162be122a6c94d706e3b3125ca"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.17.38"
weakdeps = ["SparseArrays"]

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "c9b163bd832e023571e86d0b90d9de92a9879088"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.6"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.BoundaryValueDiffEq]]
deps = ["ArrayInterface", "BandedMatrices", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "NonlinearSolve", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "TruncatedStacktraces", "UnPack"]
git-tree-sha1 = "f7392ce20e6dafa8fee406142b1764de7d7cd911"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "4.0.1"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "601f7e7b3d36f18790e2caf83a882d88e9b71ff1"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.4"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

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

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools", "SHA"]
git-tree-sha1 = "5e21a254d82c64b1a4ed9dbdc7e87c5d9cf4a686"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.10.12"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.Chain]]
git-tree-sha1 = "8c4920235f6c561e401dfe569beb8b924adad003"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "2118cb2765f8197b08e5958cdd17c165427425ee"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.19.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChaosTools]]
deps = ["Combinatorics", "DSP", "Distances", "Distributions", "DynamicalSystemsBase", "IntervalRootFinding", "LinearAlgebra", "LombScargle", "Neighborhood", "Optim", "ProgressMeter", "Random", "Reexport", "Roots", "SpecialFunctions", "Statistics", "StatsBase"]
git-tree-sha1 = "a7d881f0ec1dbe089e4edd6315894e3293f8ae6f"
uuid = "608a59af-f2a3-5ad4-90b4-758bdf3122a7"
version = "3.1.2"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "70232f82ffaab9dc52585e0dd043b5e0c6b714f1"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.12"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "886826d76ea9e72b35fcd000e535588f7b60f21d"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

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
git-tree-sha1 = "f7f4319567fe769debfcf7f8c03d8da1dd4e2fb0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.9"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataFramesMeta]]
deps = ["Chain", "DataFrames", "MacroTools", "OrderedCollections", "Reexport"]
git-tree-sha1 = "6970958074cd09727b9200685b8631b034c0eb16"
uuid = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
version = "0.14.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefaultApplication]]
deps = ["InteractiveUtils"]
git-tree-sha1 = "c0dfa5a35710a193d83f03124356eef3386688fc"
uuid = "3f0dd361-4fe0-5fc6-8523-80b14ec94d85"
version = "1.1.0"

[[deps.DelaunayTriangulation]]
deps = ["DataStructures", "EnumX", "ExactPredicates", "Random", "SimpleGraphs"]
git-tree-sha1 = "26eb8e2331b55735c3d305d949aabd7363f07ba7"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "0.8.11"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack"]
git-tree-sha1 = "e40378efd2af7658d0a0579aa9e15b17137368f4"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.44.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DataStructures", "DocStringExtensions", "EnumX", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "Requires", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Static", "StaticArraysCore", "Statistics", "Tricks", "TruncatedStacktraces", "ZygoteRules"]
git-tree-sha1 = "0d9982e8dee851d519145857e79a17ee33ede154"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.130.0"

    [deps.DiffEqBase.extensions]
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"
    DiffEqBaseZygoteExt = "Zygote"

    [deps.DiffEqBase.weakdeps]
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "Functors", "LinearAlgebra", "Markdown", "NLsolve", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "d0b94b3694d55e7eedeee918e7daee9e3b873399"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.35.0"
weakdeps = ["OrdinaryDiffEq", "Sundials"]

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "319377c927a4aa1f491228b2ac23f3554a3497c6"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.20.0"

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
git-tree-sha1 = "96a19f498504e4a3b39524196b73eb60ccef30e9"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.10.0"

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

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "9242eec9b7e2e14f9952e8ea1c7e31a50501d587"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.104"

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
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "51b4b84d33ec5e0955b55ff4b748b99ce2c3faa9"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.6.7"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "fea68c84ba262b121754539e6ea0546146515d4f"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.5.3"

[[deps.DynamicalSystemsBase]]
deps = ["ForwardDiff", "LinearAlgebra", "OrdinaryDiffEq", "Reexport", "Roots", "SciMLBase", "SparseArrays", "StateSpaceSets", "Statistics"]
git-tree-sha1 = "e375b22ff119e960d4f0262b60bf48e4fba1ad88"
uuid = "6e36e845-645a-534a-86f2-f5d4aa5a06b4"
version = "3.4.3"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EnzymeCore]]
deps = ["Adapt"]
git-tree-sha1 = "2efe862de93cd87f620ad6ac9c9e3f83f1b2841b"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.6.4"

[[deps.ErrorfreeArithmetic]]
git-tree-sha1 = "d6863c556f1142a061532e79f611aa46be201686"
uuid = "90fa49ef-747e-5e6f-a989-263ba693cf1a"
version = "0.5.2"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArraysCore", "Test"]
git-tree-sha1 = "276e83bc8b21589b79303b9985c321024ffdf59c"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.5"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "602e4585bcbd5a25bc06f514724593d13ff9e862"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.25.0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.Extents]]
git-tree-sha1 = "2140cd04483da90b2da7f99b2add0750504fc39c"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "ec22cbbcd01cba8f41eecd7d44aac1f23ee985e3"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.7.2"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "a6e756a880fc419c8b41592010aebe6a5ce09136"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.2.8"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "b12f05108e405dadcc2aff0008db7f831374e051"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "2.0.0"

[[deps.FastRounding]]
deps = ["ErrorfreeArithmetic", "LinearAlgebra"]
git-tree-sha1 = "6344aa18f654196be82e62816935225b3b9abe44"
uuid = "fa42c844-2597-5d31-933b-ebd51ab2693f"
version = "0.3.1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "299dc33549f68299137e51e6d49a13b5b1da9673"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "5b93957f6dcd33fc343044af3d48c215be2562f1"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.9.3"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "73d1214fec245096717847c62d389a5d2ac86504"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.22.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "055626e1a35f6771fe99060e835b72ca61a52621"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.1"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

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
git-tree-sha1 = "9a68d75d466ccc1218d0552a8e1631151c569545"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.5"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "d4f85701f569584f2cff7ba67a137d03f0cfb7d0"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.3"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "424a5a6ce7c5d97cca7bcc4eac551b97294c54af"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.9"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "899050ace26649433ef1af25bc17a815b3db52b7"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.9.0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "f57a64794b336d4990d90f80b147474b869b1bc4"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.2"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "ExprTools", "Logging", "MultivariatePolynomials", "Primes", "Random", "SIMD", "SnoopPrecompile"]
git-tree-sha1 = "44f595de4f6485ab5ba71fe257b5eadaa3cf161e"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.4.4"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6df9cd6ee79fc59feab33f63a1b3c9e95e2461d5"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.2"

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
git-tree-sha1 = "eb8fed28f4994600e29beef49744639d985a04b2"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.16"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

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
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "fc5d1d3443a124fde6e92d0260cd9e064eba69f8"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.1"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "bca20b2f5d00c4fbc192c3212da8fa79f4688009"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.7"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3d09a9f60edf77f8a4d99f9e015e8fbf9989605d"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.7+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "ea8031dea4aff6bd41f1df8f2fdfb25b33626381"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.4"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5fdf2fe6724d8caabf43b557b84ce53f3b7e2f6b"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.0.2+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

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
deps = ["Dates", "Random"]
git-tree-sha1 = "3d8866c029dd6b16e69e0d4a939c4dfcb98fac47"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.8"
weakdeps = ["Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

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
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "fa6d0bcff8583bac20f1ffa708c3913ca605c611"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.5"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60b1194df0a3298f460063de985eae7b01bc011a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.1+0"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "UnPack"]
git-tree-sha1 = "c451feb97251965a9fe40bacd62551a72cc5902c"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.10.1"
weakdeps = ["FastBroadcast"]

    [deps.JumpProcesses.extensions]
    JumpProcessFastBroadcastExt = "FastBroadcast"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "884c2968c2e8e7e6bf5956af88cb46aa745c854b"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.4.1"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "fee018a29b60733876eb557804b5b109dd3dd8a7"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.8"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "8a6837ec02fe5fb3def1abc907bb802ef11a0729"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.5"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "f12f2225c999886b69273f84713d1b9cb66faace"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.15.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "62edfee3211981241b57ff1cedf4d74d79519277"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.15"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

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
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "3a994404d3f6709610701c7dabfc03fed87a81f8"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.1"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearAlgebraX]]
deps = ["LinearAlgebra", "Mods", "Permutations", "Primes", "SimplePolynomials"]
git-tree-sha1 = "89ed93300377e0742ae8a7423f7543c8f5eb73a4"
uuid = "9b3f67b0-2d00-526e-9884-9e4938f8fb88"
version = "0.2.5"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ConcreteStructs", "DocStringExtensions", "EnumX", "EnzymeCore", "FastLapackInterface", "GPUArraysCore", "InteractiveUtils", "KLU", "Krylov", "Libdl", "LinearAlgebra", "MKL_jll", "PrecompileTools", "Preferences", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Sparspak", "SuiteSparse", "UnPack"]
git-tree-sha1 = "9f27ba34f5821a0495efb09ea3a465c31326495a"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "2.10.0"

    [deps.LinearSolve.extensions]
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveEnzymeExt = "Enzyme"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = "Pardiso"

    [deps.LinearSolve.weakdeps]
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

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

[[deps.LombScargle]]
deps = ["FFTW", "LinearAlgebra", "Measurements", "Random", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "d64a0ce7539181136a85fd8fe4f42626387f0f26"
uuid = "fc60dff9-86e7-5f2f-a8a0-edeadbb75bd9"
version = "1.0.3"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "0f5648fbae0d015e3abe5867bca2b362f67a5894"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.166"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "72dc3cf284559eb8f53aa593fe62cb33f83ed0c0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.0.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "b211c553c199c111d998ecdaf7623d1b89b69f93"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.12"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Setfield", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "StableHashTraits", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "35fa3c150cd96fd77417a23965b7037b90d6ffc9"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.19.12"

[[deps.MakieCore]]
deps = ["Observables", "REPL"]
git-tree-sha1 = "9b11acd07f21c4d035bd4156e789532e8ee2cc70"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.9"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "96ca8a313eb6437db5ffe946c457a401bbb8ce1d"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

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
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mods]]
git-tree-sha1 = "9d292c7fb23e9a756094f8617a0f10e3b9582f47"
uuid = "7475f97c-0381-53b1-977b-4c60186c8d62"
version = "2.2.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.Multisets]]
git-tree-sha1 = "8d852646862c96e226367ad10c8af56099b4047e"
uuid = "3b2b4ff1-bcff-5658-a3ee-dbcf1ce5ac09"
version = "0.4.4"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "6ffb234d6d7c866a75c1879d2099049d3a35a83a"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.3"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "806eea990fb41f9b36f1253e5697aa645bf6a9f8"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.4.0"

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
git-tree-sha1 = "ded64ff6d4fdd1cb68dfcbb818c69e144a5b2e4c"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.16"

[[deps.Neighborhood]]
deps = ["Distances", "NearestNeighbors", "Random", "Test"]
git-tree-sha1 = "fdea60ca30d724e76cc3b3d90d7f9d29d3d5cab5"
uuid = "645ca80c-8b79-4109-87ea-e1f58159d116"
version = "0.2.4"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ArrayInterface", "DiffEqBase", "EnumX", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "PrecompileTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SparseArrays", "SparseDiffTools", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "e10debcea868cd6e51249e8eeaf191c25f68a640"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "1.10.1"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "6a731f2b5c03157418a20c12195eb4b74c8f8621"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.13.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "a4ca623df1ae99d09bc9868b008262d0c0ac1e4f"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.4+0"

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
git-tree-sha1 = "a12e56c72edee3ce6b96667745e6cbbe5498f200"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "01f85d9269b13fedc61e63cc72ee2213565f7a72"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.8"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "IfElse", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLNLSolve", "SciMLOperators", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "f0f43037c0ba045e96f32d65858eb825a211b817"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.58.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PGFPlotsX]]
deps = ["ArgCheck", "Dates", "DefaultApplication", "DocStringExtensions", "MacroTools", "OrderedCollections", "Parameters", "Requires", "Tables"]
git-tree-sha1 = "3e7a0345b9f37da2cd770a5d47bb5cb6e62c7a81"
uuid = "8314cec4-20b6-5062-9cdb-752b83310925"
version = "1.6.0"
weakdeps = ["Colors", "Contour", "Measurements", "StatsBase"]

    [deps.PGFPlotsX.extensions]
    ColorsExt = "Colors"
    ContourExt = "Contour"
    MeasurementsExt = "Measurements"
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "67186a2bc9a90f9f85ff3cc8277868961fb57cbd"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.3"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4745216e94f71cb768d58330b059c9b76f32cb66"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.14+0"

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

[[deps.Permutations]]
deps = ["Combinatorics", "LinearAlgebra", "Random"]
git-tree-sha1 = "c7745750b8a829bc6039b7f1f0981bcda526a946"
uuid = "2ae35dd2-176d-5d53-8349-f30d82d94d4f"
version = "0.4.19"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "862942baf5663da528f66d24996eb6da85218e76"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "bd7c69c7f7173097e7b5e1be07cee2b8b7447f51"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.54"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "fca25670784a1ae44546bcb17288218310af2778"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.9"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "3aa2bb4982e575acd7583f01531f241af077b163"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.2.13"
weakdeps = ["ChainRulesCore", "MakieCore", "MutableArithmetics"]

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

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
git-tree-sha1 = "13078a8afc9f88ede295c17c3c673eb52c05b0b4"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.16"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "88b895d13d53b5577fd53379d913b9ab9ac82660"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.1"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "1d05623b5952aed1307bf8b43bec8b8d1ef94b6e"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.5"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "00099623ffee15972c16111bcf84c58a0051257c"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.9.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "c860e84651f58ce240dd79e5d9e055d55234c35a"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.6.2"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "b8a399e95663485820000f26b6a43c794e166a49"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.4"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

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
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "d7087c013e8a496ff396bae843b1e16d9a30ede8"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.10"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "PrecompileTools", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "8bc86c78c7d8e2a5fe559e3721c0f9c9e303b2ed"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.21"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

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

[[deps.RingLists]]
deps = ["Random"]
git-tree-sha1 = "f39da63aa6d2d88e0c1bd20ed6a3ff9ea7171ada"
uuid = "286e9d63-9694-5540-9e3c-4e6708fa07b2"
version = "0.2.8"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.Roots]]
deps = ["ChainRulesCore", "CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "0f1d92463a020321983d04c110f476c274bafe2e"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.22"

    [deps.Roots.extensions]
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
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
git-tree-sha1 = "6aacc5eefe8415f47b3e34214c1d79d2674a0ba2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.12"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "d8911cc125da009051fb35322415641d02d9e37f"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.6"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "3aac6d68c5e57449f5b9b865c9ba50ac2970c4cf"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.42"

[[deps.SciMLBase]]
deps = ["ADTypes", "ArrayInterface", "ChainRulesCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FillArrays", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces", "ZygoteRules"]
git-tree-sha1 = "916b8a94c0d61fa5f7c5295649d3746afb866aff"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.98.1"

    [deps.SciMLBase.extensions]
    ZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLNLSolve]]
deps = ["DiffEqBase", "LineSearches", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "765b788339abd7d983618c09cfc0192e2b6b15fd"
uuid = "e9a6253c-8580-4d32-9898-8661bb511710"
version = "0.1.9"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "51ae235ff058a64815e0a2c34b1db7578a06813d"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.7"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "0e7508ff27ba32f26cd459474ca2ede1bc10991f"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SetRounding]]
git-tree-sha1 = "d7a25e439d07a17b7cdf97eecee504c50fedf5f6"
uuid = "3cc68bcd-71a2-5612-b932-767ffbe40ab0"
version = "0.2.1"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "db0219befe4507878b1a90e07820fed3e62c289d"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.4.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleGraphs]]
deps = ["AbstractLattices", "Combinatorics", "DataStructures", "IterTools", "LightXML", "LinearAlgebra", "LinearAlgebraX", "Optim", "Primes", "Random", "RingLists", "SimplePartitions", "SimplePolynomials", "SimpleRandom", "SparseArrays", "Statistics"]
git-tree-sha1 = "f65caa24a622f985cc341de81d3f9744435d0d0f"
uuid = "55797a34-41de-5266-9ec1-32ac4eb504d3"
version = "0.8.6"

[[deps.SimpleNonlinearSolve]]
deps = ["ArrayInterface", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "PackageExtensionCompat", "PrecompileTools", "Reexport", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "15ff97fa4881133caa324dacafe28b5ac47ad8a2"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "0.1.23"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveNNlibExt = "NNlib"

    [deps.SimpleNonlinearSolve.weakdeps]
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"

[[deps.SimplePartitions]]
deps = ["AbstractLattices", "DataStructures", "Permutations"]
git-tree-sha1 = "e9330391d04241eafdc358713b48396619c83bcb"
uuid = "ec83eff0-a5b5-5643-ae32-5cbf6eedec9d"
version = "0.3.1"

[[deps.SimplePolynomials]]
deps = ["Mods", "Multisets", "Polynomials", "Primes"]
git-tree-sha1 = "7063828369cafa93f3187b3d0159f05582011405"
uuid = "cc47b68c-3164-5771-a705-2bc0097375a0"
version = "0.2.17"

[[deps.SimpleRandom]]
deps = ["Distributions", "LinearAlgebra", "Random"]
git-tree-sha1 = "3a6fb395e37afab81aeea85bae48a4db5cd7244a"
uuid = "a6525b86-64cd-54fa-8f65-62fc48bdc0e8"
version = "0.3.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Random", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "Tricks", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "c281e11db4eacb36a292a054bac83c5a0aca2a26"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.15.0"

    [deps.SparseDiffTools.extensions]
    SparseDiffToolsEnzymeExt = "Enzyme"
    SparseDiffToolsSymbolicsExt = "Symbolics"
    SparseDiffToolsZygoteExt = "Zygote"

    [deps.SparseDiffTools.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableHashTraits]]
deps = ["Compat", "SHA", "Tables", "TupleTools"]
git-tree-sha1 = "5a26dfe46e2cb5f5eca78114c7d49548b9597e71"
uuid = "c5dd0088-6c3f-4803-b00e-f31a60c170fa"
version = "1.1.3"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StateSpaceSets]]
deps = ["Distances", "LinearAlgebra", "Neighborhood", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "fcfc8b9ce0f43bcb39e6bd8664855dc501a63886"
uuid = "40b095a5-5852-4c12-98c7-d43bf788e795"
version = "1.4.5"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "f295e0a1da4ca425659c57441bcb59abb035a4bc"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.8"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Requires", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "5d66818a39bb04bf328e92bc933ec5b4ee88e436"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.5.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "42d5373c10272d14ef49cc68ffc22df3b93c549a"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.8.2"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.SteadyStateDiffEq]]
deps = ["DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "2ca69f4be3294e4cd987d83d6019037d420d9fc1"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.16.1"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "b341540a647b39728b6d64eaeda82178e848f76e"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.62.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "d6415f66f3d89c615929af907fdc6a3e17af0d8c"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.2"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.StructArrays]]
deps = ["Adapt", "ConstructionBase", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "0a3db38e4cce3c54fe7a71f831cd7b6194a54213"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.16"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "71dc65a2d7decdde5500299c9b04309e0138d1b4"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.20.1"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "SuiteSparse_jll", "libblastrampoline_jll"]
git-tree-sha1 = "ba4d38faeb62de7ef47155ed321dce40a549c305"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.2+0"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "f8ab052bfcbdb9b48fad2c80c873aa0d0344dfe5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "2f3fa844bcd33e40d8c29de5ee8dded7a0a70422"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.4.0"

[[deps.Symbolics]]
deps = ["ArrayInterface", "Bijections", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "ac7f8825d029b568f82dbf2cb49da9cebcadaffb"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.5.3"

    [deps.Symbolics.extensions]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

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
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "34cc045dd0aaa59b8bbe86c644679bc57f1d5bd0"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.8"

[[deps.TikzPictures]]
deps = ["LaTeXStrings", "Poppler_jll", "Requires", "tectonic_jll"]
git-tree-sha1 = "79e2d29b216ef24a0f4f905532b900dcf529aa06"
uuid = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
version = "3.5.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f548a9e9c490030e545f72074a41edfd0e5bcdd7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.23"

[[deps.TranscodingStreams]]
git-tree-sha1 = "1fbeaaca45801b4ba17c251dd8603ef24801dd84"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.2"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Trapz]]
git-tree-sha1 = "79eb0ed763084a3e7de81fe1838379ac6a23b6a0"
uuid = "592b5752-818d-11e9-1e9a-2b8ca4a44cd1"
version = "2.0.3"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "fadebab77bf3ae041f77346dd1c290173da5a443"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.20"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.TupleTools]]
git-tree-sha1 = "155515ed4c4236db30049ac1495e2969cc06be9d"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.4.3"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "bb37ed24f338bc59b83e3fc9f32dd388e5396c53"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.12.4"

[[deps.UnitfulAngles]]
deps = ["Dates", "Unitful"]
git-tree-sha1 = "d6cfdb6ddeb388af1aea38d2b9905fa014d92d98"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.6.2"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "05adf5e3a3bd1038dd50ff6760cddd42380a7260"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.2.0"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "7209df901e6ed7489fe9b7aa3e46fb788e15db85"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.65"

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
git-tree-sha1 = "801cbe47eae69adc50f36c3caec4758d2650741b"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.2+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

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
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "9d749cd449fb448aeca4feee9a2f4186dbb5d184"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.4"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "93284c28274d9e75218a416c65ec49d0e0fcdf3d"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.40+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.tectonic_jll]]
deps = ["Artifacts", "Fontconfig_jll", "FreeType2_jll", "Graphite2_jll", "HarfBuzz_ICU_jll", "HarfBuzz_jll", "ICU_jll", "JLLWrappers", "Libdl", "OpenSSL_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "54867b00af20c70b52a1f9c00043864d8b926a21"
uuid = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"
version = "0.13.1+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╠═b03ee99c-27f4-47df-bba5-2ea3dabdb45d
# ╟─08df960b-fd82-43ba-a9dc-bf5e83af587e
# ╟─cbd51460-8ef0-49eb-8219-14986d8421e4
# ╟─5814a7b3-8420-4a57-a2a2-d8c59db29a99
# ╟─8eb6540d-f5b0-45e6-883c-0cc213e67e45
# ╟─a842b24e-8d26-41ab-9de3-91632aede893
# ╠═dadc3e3a-ebe7-4f13-a03b-ab988094321a
# ╟─64787011-b5b8-42be-b6e4-37ebc5138b3e
# ╟─14c7f574-0623-4254-b8f7-97984d32351c
# ╟─b2b23d48-9c3d-44d3-9106-745eecc9b561
# ╠═047bbd39-9cf9-4bd7-b38e-16aa505b0b08
# ╟─35e194f5-20dc-4391-b761-3696fe0bc117
# ╟─43eafb0f-08a8-4e53-9017-50b97ac48a52
# ╟─c843a9e7-c0e1-42c1-bace-c866f777232f
# ╟─70078b44-4d66-49b9-930e-74261df8be78
# ╟─2fe0dc4c-da44-4fc8-bef8-1fa615a0fe4a
# ╟─744a9591-c7f1-496e-9bb4-47df2c8937dd
# ╠═34b04cf3-dabe-4364-8124-c5f3f351edb2
# ╠═af69ab25-0f06-4837-ac35-acbe38a4ffb1
# ╠═806db782-1734-4112-ab7b-84e03f4c342d
# ╠═1743b6dc-0a4e-4a02-90fe-3bc47833421a
# ╟─57ade87f-f93c-4b43-a737-a1b44f8af4fc
# ╟─7e3183ae-b409-457a-8e35-4dc87894c580
# ╠═d5446d26-4b59-46fa-a1f0-8374d5e05194
# ╠═7a744b2e-56ee-4162-a23d-dca9d0657608
# ╠═d0ecada4-13d1-4977-b66e-d992a1adae15
# ╟─408721dd-bd85-47b4-a75c-6e3c188f5d64
# ╠═534a1049-8de5-4b07-abec-c5a3456627c0
# ╠═86b692f1-0268-40f3-b4a2-d54c9828346d
# ╟─eaf272c7-4162-4a9a-92e3-9835c6158394
# ╟─dc6fd12b-c821-4e20-a896-25c8aab9df94
# ╟─ac553b12-4857-4cc1-8ea2-fe9e8863b429
# ╠═cf488d0e-3294-45b2-b40c-fa18062c97d2
# ╠═bc67cf22-4caa-497d-aae9-e5d1191468e2
# ╟─1d27ec35-65ca-4c94-9e8d-54d1c11e759f
# ╠═68732d91-805a-4663-9166-f8483213a8d2
# ╠═27281e53-e519-4ad0-af5d-59fb0e208534
# ╟─327fd38a-5ff6-4ac4-8d29-694272d9d46f
# ╠═9a24d3ac-a238-4eef-afc0-00fa7ef51475
# ╠═897909e2-dcad-4ef6-9161-fd3654160dba
# ╠═00030fd8-a9db-4903-b2ed-21a64db30588
# ╠═d4f91aa3-183a-4abf-8f7a-7a05d4333e3a
# ╠═7e824ce1-1f82-48cc-a3c4-1acfba0e2100
# ╠═f2cc8e6f-737c-42e3-bb81-42c50d62cf78
# ╟─4a7eb24b-0874-49a3-9b08-4ffb6a7f0ce7
# ╠═f2a6676f-457a-476a-9ce7-c336aa9bf47f
# ╠═1734df7f-1309-4ebd-a021-5f75f0bb78b2
# ╠═4f7de8a3-7f59-4a7b-8980-53390e52e0d1
# ╠═7ad23ea4-9887-4de5-8a5c-c37ebef736b8
# ╟─3767c7f9-a0bc-467a-a20a-5e5a266111c7
# ╟─f65d84cd-ab5f-4270-98ba-568792d1fec1
# ╟─34faac11-85a2-44dc-bd8d-1a71656fccf4
# ╟─f8b02d00-ff30-480e-b5eb-e150e4678c95
# ╟─44c88ad8-a8c3-45e3-9a56-be3ce5bf66fa
# ╟─448e1dee-4628-4c14-9d6f-dc165b2e826e
# ╠═5311b7cc-7199-45c8-b5e7-20d3ceb191b7
# ╠═3e637368-6bdb-4d22-9a4a-df23c6682c2f
# ╠═ef65a096-cc2a-4ce6-a06b-8671c99ca777
# ╟─a0294888-90cf-4e5b-a4b8-ce2c63bdae7a
# ╠═303da7a3-e574-4f7e-9acc-83ed83a09f69
# ╠═994f97fb-1c30-4825-9b29-35fe4ade8fb3
# ╠═493e649f-58ce-4b90-9b61-f9f8a6146ed9
# ╟─533b3cd0-c1f6-4ecd-b196-4ed35bf77135
# ╟─be85ba3b-5439-4cf3-bb14-d24d61a283c3
# ╟─b3a260b6-eb31-43a0-9fd6-60a507984319
# ╠═5ba3a0c1-6107-45a1-9b1d-5c323b9a7145
# ╠═7cbf5573-032e-4ddd-9575-f387a577c93e
# ╠═1b044783-0f5f-4321-abda-35e5b7ae67c4
# ╠═1d0a66d0-5791-4dc9-a5d4-7882b0e91767
# ╠═bb4ee2dc-c069-415e-bac5-0130950f3941
# ╟─041916ac-4cb0-4630-a227-043fae52264d
# ╠═57ea5e31-d156-4df5-bb77-0bc01b3559af
# ╟─9666bdc8-cbc0-4757-9bd8-a76477c252eb
# ╟─ca9a233b-d3ca-4a76-a3d8-f29884ac9484
# ╠═d8bee772-3979-42cd-9e38-8df0925b4e6b
# ╟─e2e4ae4f-dcdc-4999-88f2-853378be859a
# ╠═177f8253-6c35-495b-9119-ce5e8e15cba8
# ╠═b3969810-ab25-4e91-ad5a-80560b80977e
# ╠═2620d8a6-030d-4a6f-911c-6552072ff7a1
# ╠═b76d4669-26dc-48cb-930f-5e40dd40a9f1
# ╠═69b8d934-c031-413b-9c86-3fbd64be5a4a
# ╟─35ac9289-ba53-453f-9d9e-ef3499949a98
# ╠═64e7e6aa-4265-4de2-a9c6-474d125b45cc
# ╟─ff37c5e7-bb7e-4e3a-a16c-b8b7c4072528
# ╠═ccfafe47-841a-4c8a-b609-de34b453f7ee
# ╠═5e402cd6-ebcb-4641-8ce8-2f7c6c37a5ea
# ╠═e3a6d4ab-fb1c-4fa7-b4ca-50552b720ab5
# ╠═0c4fb275-240f-4689-857c-fb30120af9a5
# ╟─e050e2e9-30c1-4dc1-a57b-4f86821e3965
# ╠═d8be4069-e589-495f-a6d8-0b1fcc048dd2
# ╠═16b56ce2-6b02-4d53-a7a7-7cd7f96f26c5
# ╟─4607856c-7472-4131-a2ee-29f7150f5cb4
# ╠═bbb7263a-91e4-4a23-9e5f-416b6b7fcf6e
# ╟─0cac7dbf-31bb-459b-9291-8092cfb49354
# ╠═e65aef30-bf57-4d7d-b878-276ea6c0ef4a
# ╠═729bb0e3-3d60-4ac6-803e-2568a6970d44
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
