### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ b03ee99c-27f4-47df-bba5-2ea3dabdb45d
using CairoMakie, DataFrames, DataFramesMeta, DelimitedFiles, DifferentialEquations, Interpolations, LinearAlgebra, PlutoUI, QuadGK, Symbolics, TikzPictures, Trapz, Unitful, UnitfulAstro

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

The star formation rate (SFR) is a key characteristic of galaxies. In the context of the standard cosmological model, the SFR is determined by a combination of various processes that take place over the course of a galaxy's lifetime, such as gas cooling, star formation, chemical enrichment, and feedback from supernovae and galactic nuclei. These processes are influenced by factors like mergers, interactions, and mass accretion, which affect the amount and properties of the gas from which stars form. The density of a gas cloud is believed to be the most important factor in determining its star formation rate, although the details of this process are not yet fully understood. Observationally, the total gas density is found to be correlated to the star formation rate ([Kennicutt1998](https://doi.org/10.1086/305588)), and this correlation is even stronger when considering the molecular gas ([Wong2002](https://doi.org/10.1086/339287), [Bigiel2008](https://doi.org/10.1088/0004-6256/136/6/2846)). The underlying reason is the intrinsic relation between molecular gas mass and SFR, which can be found at resolved scales ([Baker2021](https://doi.org/10.1093/mnras/stab3672)) and at integrated (i.e. galaxy wide) scales across redshifts ([Baker2022](https://doi.org/10.1093/mnras/stac3413)).

As the formation of dark matter halos and galaxies is highly non-linear, numerical simulations have become the preferred tool to investigate how galaxies form and evolve from early times up to the present.
This type of simulations naturally include mergers/interactions and continuous gas accretion, processes that may induce changes in the SFR.  However, there are still significant uncertainties in the modelling of the evolution of the baryonic component, since the physical processes that affect baryons – such as star formation, feedback, and chemical enrichment – take place at scales that are too small to be resolved directly. As a result, these processes are introduced using sub-grid physics, which involves several adjustable parameters that are not always independent of one another or constrained observationally. This can lead to inconsistencies in the predictions of different models ([Scannapieco2012](https://doi.org/10.1111/j.1365-2966.2012.20993.x)).

Because of its importance in galaxy formation, it is critical for simulations to accurately describe the star formation process at the scales that can be resolved, as well as the associated feedback effects.
"""
  ╠═╡ =#

# ╔═╡ 8eb6540d-f5b0-45e6-883c-0cc213e67e45
md"""
## Previous work

Broadly, there has been two ways to model the multiphase structure of the interstellar medium (MP ISM), one is based on the physical properties of the gas (hot and cold phases), and the other on its chemical composition (molecular and atomic phases).

The former was pioneer by [Field1969](https://doi.org/10.1086/180324) (see [Cowie1977](https://doi.org/10.1086/154911), [McKee1977a](https://doi.org/10.1086/155350), and for a review [Cox2005](https://doi.org/10.1146/annurev.astro.43.072103.150615)), within the context of semi analytic models (SAMs). The model develop by [McKee1977b](https://doi.org/10.1086/155667) was first incorporated into numerical simulation of galaxy formation by [Yepes1997](https://doi.org/10.1093/mnras/284.1.235) (Eulerian) and [Hultman1999](https://ui.adsabs.harvard.edu/abs/1999A%26A...347..769H) (Lagrangian). These works were later extended by [Springel2003](https://doi.org/10.1046/j.1365-8711.2003.06206.x), adding galactic winds driven by star formation as a form of feedback. 

[Monaco2004](https://doi.org/10.1111/j.1365-2966.2004.07916.x) developed a semi-analytic model in a similar vein to [Springel2003](https://doi.org/10.1046/j.1365-8711.2003.06206.x), providing the theoretical foundation for MUPPI (MUlti-Phase Particle Integrator) ([Murante2010](https://doi.org/10.1111/j.1365-2966.2010.16567.x)), a sub-resolution MP ISM model that adds stellar feedback to $\texttt{GADGET-2}$ ([Springel2005](https://doi.org/10.1111/j.1365-2966.2005.09655.x)). MUPPI separates a gas particle in a hot and cold phase if a set of conditions for its density and temperature are match. It then evolves those components, plus a stellar phase and an energy term (energy of the hot gas) using a set of four ODEs. In each SPH integration step the equations uses the gas state (pressure, density, entropy) as parameters to compute the ICs (in the case that the particle is entering to the multiphase state for the first time), and to evolve the equations. The state of the internal variables (hot mass, cold mass, stellar mass and hot phase energy) are kept to be used in the next SPH step.

The idea to model the formation and evolution of the molecular gas, within an MP ISM, was first implemented in [Pelupessy2006](https://doi.org/10.1086/504366). This was later expanded as a set of ODEs that follow the interaction between atomic gas, molecular gas and dust in [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55), and further develop by [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x) to integrate it into SPH simulations.

A SAM with molecular, atomic an ionized gas was first implemented by [Berry2014](https://doi.org/10.1093/mnras/stu613) and later [Somerville2015](https://doi.org/10.1093/mnras/stv1877), but as far as we are aware there are no MP ISM models with the three gas phases crafted into hydrodynamical codes.

Based on [Ferrini1992](https://doi.org/10.1086/171066) and later work, [Mollá2015](https://doi.org/10.1111/j.1365-2966.2005.08782.x) develop a SAM to follow the metal component in galaxies. These chemical evolution models (CEMs) where subsequently improve and extended in [Mollá2015](https://doi.org/10.1093/mnras/stv1102), [Mollá2016](https://doi.org/10.1093/mnras/stw1723), [Mollá2017](https://doi.org/10.1093/mnras/stx419) and [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635). The latter being a full SAM that models the MP ISM, considering the molecular, atomic and ionized phases of Hydrogen, plus dust and the stellar component. Our model follows closely the one developed by [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635), but implemented within the hydrodynamical code $\texttt{Arepo}$, the same way MUPPI is into $\texttt{GADGET}$.
"""

# ╔═╡ a842b24e-8d26-41ab-9de3-91632aede893
# ╠═╡ skip_as_script = true
#=╠═╡
md"""

## Phases and notation

Given that the interstellar medium (ISM) is not homogeneous, we will model it as a multi-phase structure made up of four components. Three of which are the phases of Hydrogen, the ionized phase with a temperature of $\sim \! 10^4 \, \mathrm{K}$, the atomic phase with a temperature of $\sim \! 100 \, \mathrm{K}$, and the molecular phase with a temperature of $\sim \! 10 \, \mathrm{K}$. The final component is the stars.

For every reaction that transfers mass between the phases, we will consider only the dominant channel. So, even though the gas is made up of Hydrogen and Helium, only the Hydrogen reactions are incorporated into the model. The processes involved are the photoionization of atoms; the recombination of electrons with ions; the conversion of atomic hydrogen into molecular hydrogen; and the destruction of the latter owing to the photodissociation caused by ultraviolet (UV) light from the young stellar population. In addition, we consider the formation by supernovas of ionized gas, and the influence of the neutral gas on the star formation rate.

We characterized the mass of the phases by their density fraction, with respect to the total cell density,

*  Ionized gas, $i_f(t) ≔ \rho_i / \rho_C \, ,$
*  Atomic gas, $a_f(t) ≔ \rho_a / \rho_C \, ,$
*  Molecular gas, $m_f(t) ≔ \rho_m / \rho_C \, ,$
*  Stars, $s_f(t) ≔ \rho_s / \rho_C \, ,$

where $\rho_i$, $\rho_a$, $\rho_m$, $\rho_s$ are the corresponding volume densities, and

$\begin{equation}
    \rho_C(t) ≔ \rho_i(t) + \rho_a(t) + \rho_m(t) + \rho_s(t) \, ,
\end{equation}$

is the total cell density.
"""
  ╠═╡ =#

# ╔═╡ 64787011-b5b8-42be-b6e4-37ebc5138b3e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Physical relationships

The following diagram shows the reactions between the different phases,
"""
  ╠═╡ =#

# ╔═╡ 14c7f574-0623-4254-b8f7-97984d32351c
# ╠═╡ skip_as_script = true
#=╠═╡
TikzPicture(
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

# ╔═╡ 43eafb0f-08a8-4e53-9017-50b97ac48a52
TikzPicture(
	L"""
		\node[box, white] (stars) at (180:2cm) {Stars};
		\node[box, white, text width=2em] (atom) at (0:2cm) {HI};
		\node[box, white, text width=2em] (molecule) at (270:2cm) {\ch{H2}};
		\node[box, white, text width=2em] (ion) at (90:2cm) {HII};
		\draw[line, white, ->]
		(ion) edge [bend left, "$\textcolor{d_pink}{\frac{i_f(t)}{\tau_R(t)}} - \textcolor{d_blue}{\eta_\text{ion} \, \psi(t)}$"] (atom)
		(atom) edge [bend left, "$\textcolor{d_orange}{\frac{a_f(t)}{\tau_C(t)}} - \textcolor{d_green}{\eta_\text{diss} \, \psi(t)}$"] (molecule)
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

# ╔═╡ 047bbd39-9cf9-4bd7-b38e-16aa505b0b08
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Equations"
  ╠═╡ =#

# ╔═╡ 35e194f5-20dc-4391-b761-3696fe0bc117
md"""
### Ionized gas

The ionized component growths through the ionization of atomic gas, and from the remnants of supernova explosions. 

The former is assumed to come mainly from the radiation of newborn stars, so it is given by

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}i_f(t)\right|_{\text{ion}} = \eta_\text{ion} \, \psi(t) \, ,
\end{equation}$

where $\psi(t)$ is the star formation rate 

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}s_f(t)\right|_{\text{SFR}} = \psi(t) \, ,
\end{equation}$

and $\eta_\text{ion}$ is the ionized mass rate per unit of created stellar mass. All the physics of the ionization process are summarized in the parameter $\eta_\text{ion}$.

The latter, under the instantaneous recycling hypothesis, can be written as 

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}i_f(t)\right|_{\text{recyc}} = R \, \psi(t) \, ,
\end{equation}$

where $R$ is the mass fraction of a stellar population that is returned to the ISM, where we assumed that all the returned mass is in the form of ionized gas. Ignoring the metal enrichment is a good approximation that does not alter the results.

### Atomic gas

The atomic component accumulates mass through the dissociation of a Hydrogen molecules, and the recombination of the ionized gas with free electrons.

The former, as with the ionized gas, is given by

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}a_f(t)\right|_{\text{diss}} = \eta_\text{diss} \, \psi(t) \, ,
\end{equation}$

where $\eta_\text{diss}$ is the disassociated mass rate per unit of created stellar mass.

The latter will depend on the mass of ionized gas present, and the time scale of recombination ($\tau_R$), so is given by 

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}a_f(t)\right|_{\text{recon}} = \frac{i_f(t)}{\tau_R(t)} \, .
\end{equation}$

### Molecular gas

The molecular component gains mass mainly by the condensation of Hydrogen atoms in the surface of dust grains. This process depends on the mass of atomic gas, and the characteristic time scale of condensation ($\tau_C$). We are ignoring all the dust physics, and condensing that into the single time parameter,

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}m_f(t)\right|_{\text{cond}} = \frac{a_f(t)}{\tau_C(t)} \, ,
\end{equation}$
"""

# ╔═╡ 70078b44-4d66-49b9-930e-74261df8be78
md"""
From all the above, we can write the following system of four ODEs,
"""

# ╔═╡ 2fe0dc4c-da44-4fc8-bef8-1fa615a0fe4a
# ╠═╡ skip_as_script = true
#=╠═╡
TikzPicture(
	L"""
	\node[white] {
  	${\boldmath
	\begin{aligned}
		\dv{}{t}i_f(t) &= - \textcolor{d_pink}{\frac{i_f(t)}{\tau_R(t)}} + \textcolor{d_blue}{\eta_\text{ion} \, \psi(t)} + \textcolor{d_yellow}{R \, \psi(t)} \, , \\
		\dv{}{t}a_f(t) &= \textcolor{d_pink}{\frac{i_f(t)}{\tau_R(t)}} - \textcolor{d_blue}{\eta_\text{ion} \, \psi(t)} - \textcolor{d_orange}{\frac{a_f(t)}{\tau_C(t)}} + \textcolor{d_green}{\eta_\text{diss} \, \psi(t)} \, , \\
		\dv{}{t}m_f(t) &= \textcolor{d_orange}{\frac{a_f(t)}{\tau_C(t)}} - \textcolor{d_green}{\eta_\text{diss} \, \psi(t)} - \textcolor{red}{\psi} \, , \\
		\dv{}{t}s_f(t) &= \textcolor{red}{\psi(t)} - \textcolor{d_yellow}{R \, \psi(t)} \, ,
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

From the equations we can explicitly check mass conservation,

$\begin{equation}
    \frac{\mathrm{d}}{\mathrm{d}t}(i_f(t) + a_f(t) + m_f(t) + s_f(t)) = 0 \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ af69ab25-0f06-4837-ac35-acbe38a4ffb1
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Stiffness coefficient"
  ╠═╡ =#

# ╔═╡ 57ade87f-f93c-4b43-a737-a1b44f8af4fc
# ╠═╡ skip_as_script = true
#=╠═╡
md"The model, for a typical set of initial conditions, is very stiff (stiffness coefficient >> 1)."
  ╠═╡ =#

# ╔═╡ 534a1049-8de5-4b07-abec-c5a3456627c0
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Units

We have the freedom to chose three independent units for this model. Time, mass and length. This choice is reflected in the the constants $C_S$, $C_R$, and $C_C$. 

Following the standard in astronomy and astrophysics, we choose $\mathrm{[T] = Myr}$, $\mathrm{[M] = mp}$ and $\mathrm{[L] = cm}$, where $\mathrm{mp}$ is the proton mass.
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

*  $\tau_R$: Time scale of atomic gas formation from ionized gas, generally called recombination time.

*  $\tau_C$: Time scale of molecular gas formation from atomic gas, generally called condensation (or cloud formation) time.

*  $\eta_d$: Rate of molecular gas dissociation by stars, per unit of created stellar mass.

*  $\eta_i$: Rate of atomic gas ionization by stars, per unit of created stellar mass.

*  $R$: Mass of ionized gas produced per unit of created stellar mass.

*  $Z_{SN}$: Mass of metals produced per unit of ionized gas created during the stellar life cycle.
"""
  ╠═╡ =#

# ╔═╡ ac553b12-4857-4cc1-8ea2-fe9e8863b429
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Star formation rate

For the star formation rate we take into account the strong correlation between molecular Hydrogen and star formation ([Bigiel2008](https://doi.org/10.1088/0004-6256/136/6/2846), [Bigiel2010](https://doi.org/10.1088/0004-6256/140/5/1194), [Wong2002](https://doi.org/10.1086/339287), [Robertson2008](https://doi.org/10.1086/587796), [Halle2013](https://doi.org/10.1051/0004-6361/201220952), [Thompson2013](https://doi.org/10.1088/0004-637x/780/2/145)). In particular we will follow [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635) with

$\begin{equation}
	\mathrm{SFR} = \frac{A_w \, m_a + M_w \, m_m}{\tau_S}
\end{equation}$
where $A_w$ y $M_w$ are dimensionless free parameters that weight the contribution of the atomic mass, $m_a$, and the molecular mass, $m_m$, respectively. $\tau_S$ is the characteristic timescale for star formation.

Given that the equations are written per unit of volume and cell density, $\rho_C$, the SFR enters the equations as

$\begin{equation}
	\psi = \mathrm{SFR} \, / \, V / \rho_C = \frac{A_w \, a_f + M_w \, m_f}{\tau_S} \, ,
\end{equation}$
where $a_f$ and $m_f$ are the mass fractions of atomic and molecular gas, respectively.
"""
  ╠═╡ =#

# ╔═╡ cf488d0e-3294-45b2-b40c-fa18062c97d2
begin
	const AW = 0.0
	const MW = 1.0
end;

# ╔═╡ bc67cf22-4caa-497d-aae9-e5d1191468e2
ψ(a, m, τS) = (AW * a + MW * m) / τS;

# ╔═╡ dc6fd12b-c821-4e20-a896-25c8aab9df94
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Time scales"
  ╠═╡ =#

# ╔═╡ 1d27ec35-65ca-4c94-9e8d-54d1c11e759f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Star formation time

Following [Krumholz2019](https://doi.org/10.1146/annurev-astro-091918-104430), we define $\tau_S$ as the characteristic timescale for star formation; i.e., all the gas that can be converted into stars has done so in a timescale $\tau_S$. In particular, we have

$\begin{equation}
    \tau_S = \frac{\epsilon_\star}{\epsilon_\text{ff}}\,t_\text{ff} \, ,
\end{equation}$
where $\epsilon_\star$ is the star formation efficiency (mass fraction of gas which will be converted into stars, in the literature is of $\mathcal{O}(\epsilon) \approx 1$), $t_\text{ff}$ is the free-fall time, and $\epsilon_\text{ff}$ is the star-formation efficiency per free fall time (the fraction of a cloud's mass that is transformed into stars per cloud free-fall time, in the literature is of $\mathcal{O}(\epsilon_\text{ff}) \approx 0.01$).

We have that $t_\text{ff}$ and $\epsilon_\star$ can be written as

$\begin{align}
    t_\text{ff} &= \sqrt{\frac{3\pi}{32 \, G\, \rho_g}} \, , \\
    \epsilon_\star &= s_f(t \rightarrow \infty) \, ,
\end{align}$

where $\rho_g = \rho_i + \rho_a + \rho_m$ is the density of the gas, and $s_f$ the stellar cell mass fraction.

There is a lot of uncertainty for the parameter $\epsilon_\star$ ([Lee2016](https://doi.org/10.3847/1538-4357/833/2/229) and [Utomo2018](https://doi.org/10.3847/2041-8213/aacf8f)), so we will follow [Matzner2000](https://doi.org/10.1086/317785) using $\epsilon_\star = 0.5$. This results in a net efficiency of star formation of 

$\begin{equation}
	f_\star = \frac{\epsilon_\text{ff}}{\epsilon_\star} = 0.02 \, ,
\end{equation}$

in relative agreement with [Kennicutt1998](https://doi.org/10.1086/305588) and [Evans2014](https://doi.org/10.1088/0004-637X/782/2/114). We note though that this parameter has been shown to have little influence on the global properties of simulated galaxies ([Li2018](https://doi.org/10.3847/1538-4357/aac9b8) and [Brown2022](https://doi.org/10.1093/mnras/stac1164)).

With all the previous definitions, we have 

$\begin{equation}
    \tau_S = \frac{C_S}{\sqrt{\rho_g}} = \frac{C_S}{\sqrt{(1 - s_f) \, \rho_C}} \, ,
\end{equation}$
where

$\begin{equation}
    C_S = \frac{\epsilon_\star}{\epsilon_\text{ff}} \, \sqrt{\frac{3\pi}{32 \, G}} \, .
\end{equation}$

and where we have used

$\begin{align}
    \rho_g &= \rho_C - \rho_s \\
	&= \left( 1 - \frac{\rho_s}{\rho_C} \right) \, \rho_C \\
	&= \left( 1 - s_f \right) \, \rho_C \, .
\end{align}$
"""
  ╠═╡ =#

# ╔═╡ 68732d91-805a-4663-9166-f8483213a8d2
begin
    const ϵff = 0.01
    const ϵₛ  = 0.5
	const CS  = (ϵₛ / ϵff) * sqrt(3π / 32u"G")
	const cs  = ustrip(t_u * l_u^-(3/2), CS / sqrt(m_u))
end

# ╔═╡ 27281e53-e519-4ad0-af5d-59fb0e208534
τS(s_f, ρC) = cs / sqrt((1 - s_f) * ρC);

# ╔═╡ 327fd38a-5ff6-4ac4-8d29-694272d9d46f
# ╠═╡ skip_as_script = true
#=╠═╡
let
	set_theme!(theme_black())
	
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel=L"ρ_C \, / \, \mathrm{cm^{-3}}", 
		ylabel=L"\tau_S \, / \, \mathrm{Myr}",
		xlabelsize=32,
		ylabelsize=32,
		xticklabelsize=25,
		yticklabelsize=25,
		xscale=log10,
		yscale=log10,
	)

	ρC = exp10.(range(-1, 5, 30))

	for s_f in range(0.1, 0.5, 3)
		label = L"s_f = %$(s_f)"
		lines!(ax, ρC, ρ -> τS(s_f, ρ); linewidth=3, label)
	end

	axislegend(; position=:rt, labelsize=25)

	f
end
  ╠═╡ =#

# ╔═╡ 897909e2-dcad-4ef6-9161-fd3654160dba
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Recombination time

From dimensional analysis, we have that the characteristic time for the recombination reaction

$\begin{equation}
    e^- + p^+ \longleftrightarrow H + \gamma \, ,
\end{equation}$

is given by

$\begin{equation}
    \tau_R = \frac{1}{n_e \langle\sigma \, v\rangle} \, , 
\end{equation}$

where $n_e$ is the number fraction of electrons, and $\langle\sigma\,v\rangle$ is the recombination rate.

From [Osterbrock2006](http://www.worldcat.org/oclc/60611705) (pg. 22) we have that

$\begin{equation}
    \langle\sigma\,v\rangle \approx \alpha_B(10^4 \, K) = 2.59 \times 10^{-13} \, \mathrm{cm}^3 \, \mathrm{s}^{-1} \, ,
\end{equation}$

where we only used case B recombination, which considers all but one transition channel. A transition directly to the ground state does not result in net recombination, because the emitted photon has enough energy to ionize another atom. Thus it is a better approximation to exclude that case.

The electron density is $n_e = n_i$, where $n_i$ in the number density of ionized Hydrogen atoms, so $n_e$ is essentially the same quantity as $\rho_i$, where the only difference is a conversion factor for the different units, $n_e = n_i = \rho_i / m_p$, where $m_p$ is the proton mass.

We have then

$\begin{equation}
    \tau_R = \frac{C_R}{\rho_i} = \frac{C_R}{i_f \, \rho_C} \, , 
\end{equation}$
where

$\begin{equation}
	C_R = \frac{m_p}{\alpha_B(10^4 \, K)} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 00030fd8-a9db-4903-b2ed-21a64db30588
begin
	const αB = 2.59e-13u"cm^3 * s^-1"
	const CR = m_u / αB
	const cr = ustrip(t_u * l_u^-3, CR / m_u)
end

# ╔═╡ d4f91aa3-183a-4abf-8f7a-7a05d4333e3a
τR(i_f, ρC) = cr / (i_f * ρC);

# ╔═╡ 7e824ce1-1f82-48cc-a3c4-1acfba0e2100
# ╠═╡ skip_as_script = true
#=╠═╡
let
	set_theme!(theme_black())
	
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel=L"ρ_C \, / \, \mathrm{cm^{-3}}", 
		ylabel=L"\tau_R \, / \, \mathrm{Myr}",
		xlabelsize=32,
		ylabelsize=32,
		xticklabelsize=25,
		yticklabelsize=25,
		xscale=log10,
		yscale=log10,
	)

	ρC = exp10.(range(-1, 5, 30))

	for i_f in range(0.1, 0.9, 5)
		label = L"i_f = %$(i_f)"
		lines!(ax, ρC, ρ -> τR.(i_f, ρ); linewidth=3, label)
	end

	axislegend(; position=:rt, labelsize=25)

	f
end
  ╠═╡ =#

# ╔═╡ 4a7eb24b-0874-49a3-9b08-4ffb6a7f0ce7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Condensation time

We will only consider Hydrogen reactions as a first approximation. There are several channels for the formation of molecular Hydrogen, but the most efficient ones involve the interaction of $H$ atoms on the surface of dust grains, so we have ([Mollá2017](https://doi.org/10.1093/mnras/stx419))

$\begin{equation}
    \tau_C = \frac{1}{2\,n_\mathrm{dust} \, \langle\sigma v\rangle_\mathrm{dust}}  \, , 
\end{equation}$

where $n_\mathrm{dust}$ is the number density of dust grains in the ISM, and $\langle\sigma v\rangle_\mathrm{dust}$ is the thermally averaged cross-section for the formation of molecular Hydrogen (see [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635) and reference therein).

We can write the density and cross-section as

$\begin{equation}
    n_\mathrm{dust} \, \langle\sigma v\rangle_\mathrm{dust} \approx (Z + Z_\mathrm{eff}) \, n_{ng} \, \frac{\langle\sigma v\rangle_\odot}{Z_\odot} \, , 
\end{equation}$

where $n_g$ denotes the neutral gas number density, $Z$ is the metallicity, $Z_\odot$ the solar metallicity, $\langle\sigma v\rangle_\odot = 6 \times 10^{-17} \, \mathrm{cm}^3 \, \mathrm{s}^{-1}$ (see note below), and $Z_\mathrm{eff} \approx 10^{-3} \, Z_\odot$ ([Glover2007](https://doi.org/10.1086/519445) pg. 10) is an initial value of metallicity needed to kickstart the star formation process. We need this because the initial abundance of metals and dust grains is zero, and stars only form from molecular clouds. This initial value accounts for all other channels of molecular formation, a detailed study of which would have a minimal impact on the results.

In regards to notation we note that previous works call the term $Z \, \langle\sigma v\rangle_\mathrm{dust}$, the grain-surface $\mathrm{H_2}$ formation rate coefficient (units of $\mathrm{cm^3 \, s^{-1}}$): $R$ in [Goldshmidt1995](https://doi.org/10.1086/175168) and [Draine1996](https://doi.org/10.1086/177689), $R_f$ in [Pelupessy2006](https://doi.org/10.1086/504366), and $R_d$ in [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) and [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x). For $100 \, \mathrm{K}$ and solar metallicity [Draine1996](https://doi.org/10.1086/177689) gives $R = 6 \times 10^{-17} \, \mathrm{cm}^3 \, \mathrm{s}^{-1}$, from which we got the expresion $\langle\sigma v\rangle_d = 6 \times 10^{-17} \, \mathrm{cm}^3 \, \mathrm{s}^{-1} \, Z_\odot^{-1} = \langle\sigma v\rangle_\odot \, Z_\odot^{-1}$.

A global factor multiplying $R$ can be added ([Pelupessy2006](https://doi.org/10.1086/504366) (pg.1025) and [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x) (pg. 3061) use the concentration fractor), to account for all the uncertanties.

We have $n_{ng} = \rho_{ng} / m_p$, where we used that number density $n$ is essentially the same quantity as $\rho$, the only difference being the proton mass working as a conversion factor for the different units. So, the characteristic time is given by

$\begin{equation}
    \tau_C = \frac{C_C}{(a_f + m_f) \, \rho_C \, (Z + Z_\mathrm{eff})} \, ,
\end{equation}$

where

$\begin{equation}
	C_C = \frac{Z_\odot \, m_p}{2 \, \langle\sigma v\rangle_\odot} \, .
\end{equation}$

We find several values of the solar metallicity in the literature

 * [Asplund2006](https://doi.org/10.1553/cia147s76) and [Grevesse2007](https://doi.org/10.1007/s11214-007-9173-7): $Z_\odot = 0.0122$ 
 * Used in Arepo: $Z_\odot = 0.0127$
 * [Asplund2009](https://doi.org/10.1146/annurev.astro.46.060407.145222): $Z_\odot = 0.0134$
 * [Lodders2009](https://doi.org/10.1007/978-3-540-88055-4_34): $Z_\odot = 0.0141$
 * [Caffau2010](https://doi.org/10.1007/s11207-010-9541-4): $Z_\odot = 0.0153$ 
 * [Grevesse1998](https://doi.org/10.1023/A:1005161325181): $Z_\odot = 0.0169$
 * [Steiger2016](https://doi.org/10.3847/0004-637X/816/1/13): $Z_\odot = 0.0196$

For consitency with the codebase, we will use $Z_\odot = 0.0127$, noting that is only 35% off the largest value in the list ([Steiger2016](https://doi.org/10.3847/0004-637X/816/1/13)), which is not significant for such a simple model with many other uncertanties.
"""
  ╠═╡ =#

# ╔═╡ 1335861d-d539-4986-b9c1-d883a4fe8405
md"""
 * [Hollenbach1971](https://doi.org/10.1086/150755) (expresion explicita para R)
 * [Jura1974](https://doi.org/10.1086/152975) (dn/dt = R*n*nH)
 * [Jura1975](https://doi.org/10.1086/153545)
 * [Black1987](https://doi.org/10.1086/165740)
 * [Sternberg1988](https://doi.org/10.1086/166664)
 * [Goldshmidt1995](https://doi.org/10.1086/175168)
 * [Draine1996](https://doi.org/10.1086/177689)
 * [Cazaux2002](https://doi.org/10.1086/342607)
 * [Cazaux2004](https://doi.org/10.1086/381775)
 * [Pelupessy2006](https://doi.org/10.1086/504366) (R como funcion de T y Z)
 * [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55)
 * [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x)
 * [Mollá2017](https://doi.org/10.1093/mnras/stx419)
 * [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635)
"""

# ╔═╡ f2a6676f-457a-476a-9ce7-c336aa9bf47f
begin
    const σv   = 6e-17u"cm^3 * s^-1"
    const Zsun = 0.0127
	const Zeff = 1e-3 * Zsun
	const CC   = (Zsun * m_u) / (2 * σv)
	const cc   = ustrip(t_u * l_u^-3, CC / m_u)
end

# ╔═╡ 1734df7f-1309-4ebd-a021-5f75f0bb78b2
τC(a_f, m_f, ρC, Z) = cc / ((a_f + m_f) * ρC * (Z + Zeff));

# ╔═╡ 4f7de8a3-7f59-4a7b-8980-53390e52e0d1
# ╠═╡ skip_as_script = true
#=╠═╡
let
	set_theme!(theme_black())
	
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel=L"ρ_C \, / \, \mathrm{cm^{-3}}", 
		ylabel=L"\tau_C \, / \, \mathrm{Myr}",
		xlabelsize=32,
		ylabelsize=32,
		titlesize=35,
		xticklabelsize=25,
		yticklabelsize=25,
		xscale=log10,
		yscale=log10,
		title=L"a_f + m_f = 0.5"
	)

	ρC = exp10.(range(-1, 5, 30))

	for Zs in range(0, 2, 5)
		label = L"Z \, / \, Z_\odot = %$(Zs)"
		lines!(ax, ρC, ρ -> τC.(0.0, 0.5, ρ, Zs * Zsun); linewidth=3, label)
	end

	axislegend(; position=:rt, labelsize=25)

	f
end
  ╠═╡ =#

# ╔═╡ 3767c7f9-a0bc-467a-a20a-5e5a266111c7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Photodissociation efficiency

We define the disassociated mass rate per unit of created stellar mass, as 

$\begin{equation}
    \eta_d = \frac{\dot{M}_d}{\mathrm{SFR}}  \, , 
\end{equation}$

and the ionized mass rate per unit of created stellar mass as

$\begin{equation}
    \eta_i = \frac{\dot{M}_i}{SFR}  \, .
\end{equation}$

The mass rates, $\dot{M}_d$ and $\dot{M}_i$, can be computed from the photon production rate, in the corresponding energy range,

$\begin{align}
    \dot{M}_d &= \dot{N}_i \, f_d \, , \\
    \dot{M}_i &= \dot{N}_d \, f_i \, ,
\end{align}$

where $\dot{N}_d$ is the number of photodissociating photons produced per unit time (in the Lyman–Werner band, $912\,\mathrm{Å}$ to $1107\,\mathrm{Å}$), $\dot{N}_i$ the number of ionizing photons produced per unit time (between $0$ and $912\,Å$), and $f_d$ and $f_i$ are the unit conversion factors (proton mass into solar mass). The $f$ factors allow us to take into account that the reaction may not be $100\%$ efficient too.

For the ionization reaction, each photon will produce one proton, and we assume $100\%$ efficiency.
"""
  ╠═╡ =#

# ╔═╡ f65d84cd-ab5f-4270-98ba-568792d1fec1
const f_i = ustrip(1.0u"mp" |> u"Msun")

# ╔═╡ 34faac11-85a2-44dc-bd8d-1a71656fccf4
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
For the molecular dissociation reaction, we have to consider the numerical factor given by [Draine1996](https://doi.org/10.1086/177689), where it is shown that dust grains may absorb up to $\sim \! 60$ percent of the photons capable of dissociating Hydrogen molecules, and a large fraction of the remaining photons excite different rotational and vibrational states, reducing their dissociation probability to $\sim \! 15$ percent, so we end up with an efficiency factor of $0.4 \times 0.15 = 0.06$.
$f_\mathrm{d}$ has an extra factor of two because each photon contributes with two protons (from the dissociated molecule) to the atomic gas.
"""
  ╠═╡ =#

# ╔═╡ f8b02d00-ff30-480e-b5eb-e150e4678c95
const f_d = 0.4 * 0.15 * 2.0 * ustrip(1.0u"mp" |> u"Msun")

# ╔═╡ 44c88ad8-a8c3-45e3-9a56-be3ce5bf66fa
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The number of photons can be computed from $Q(t', Z)$, which is defined as the number of photons produced by a stellar population of one solar mass, of age $t'$ and metallicity $Z$, per unit time. So, we have the relation

$\begin{equation}
    \dot{N}(t) = \int_0^t \mathrm{SFR}(t - t') \, Q(t', Z(t - t')) \mathrm{d}t' \, ,
\end{equation}$

where $\mathrm{SFR}(t - t')$ is the instantaneous SFR at the moment of birth of the stellar population of current age $t'$ (in units of solar mass per unit of time), and $Z = Z(t - t')$ is defined as the metallicity at that moment. Because most of the contribution to the integral comes from young blue stars that die in the first $10 \ \mathrm{to} \ 100 \, \mathrm{Myr}$, it is possible to approximate

$\begin{align}
    \mathrm{SFR}(t - t') &\approx \mathrm{SFR}(t) \, , \\
	Z(t - t') &\approx Z(t) \, .
\end{align}$

So we end up with 

$\begin{equation}
    \eta = f \, \frac{\dot{N}}{\mathrm{SFR}} = f \, \int_0^t Q(t', Z) \mathrm{d}t'\, ,
\end{equation}$

The value of $Q$ can be calculated using

$\begin{equation}
    Q(t', Z) = \int_{\lambda_1}^{\lambda_2} \frac{L_\lambda(t', Z)}{E_\lambda} \mathrm{d}\lambda = \int_{\lambda_1}^{\lambda_2} \frac{\lambda \, L_\lambda(t', Z)}{h \, c} \mathrm{d}\lambda \, ,
\end{equation}$

where $L_\lambda(t', Z)$ is the luminosity per unit of wavelength of a stellar population of one solar mass, age $t'$ and metallicity $Z$, and $E_\lambda = \lambda / (h \, c)$ is the energy of a photon of wavelength $\lambda$. So integrating between the wavelength of interest, we get the number of photons produced per unit time for a stellar population of the given characteristics.

The luminosity will not only depend on the age and metallicity of the population, but on the IMF (initial mass function) too, so in principle for each IMF, $t'$ and $Z$ we have a function of luminosity versus $\lambda$.

Using the values from PopStar by [Mollá2009](https://doi.org/10.1111/j.1365-2966.2009.15160.x) we compute a table of $Q$ for six IMFs ([Salpeter1955](https://doi.org/10.1086/145971) in two mass ranges, $0.85\,\mathrm{M}_\odot$ to $120\,\mathrm{M}_\odot$ and $0.15\,\mathrm{M}_\odot$ to $100\,\mathrm{M}_\odot$, [Kroupa2001](https://doi.org/10.1046/j.1365-8711.2001.04022.x), [Ferrini1990](https://ui.adsabs.harvard.edu/abs/1990A%26A...231..391F) and [Chabrier2003](https://doi.org/10.1086/374879)), six metallicities (0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05) and for ages between $0.1\,\mathrm{Myr}$ and $15\,\mathrm{Gyr}$.
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

        Qd = Vector{Quantity}(undef, length(files))
		Qi = Vector{Quantity}(undef, length(files))

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
                λ = @subset(df, λ_range_d[1] .< :λ .< λ_range_d[2])
                integrand = λ[!, 1] .* λ[!, 2] ./ (Unitful.h * Unitful.c)
                Qd[i] = trapz(λ[!, 1], integrand) |> u"s^-1"
            end
            let
                λ = @subset(df, λ_range_i[1] .< :λ .< λ_range_i[2])
                integrand = λ[!, 1] .* λ[!, 2] ./ (Unitful.h * Unitful.c)
                Qi[i] = trapz(λ[!, 1], integrand) |> u"s^-1"
            end
            
			
        end

        data_in_files[i] = identity.(
            DataFrame(
                :IMF     => uppercase.(IMF),        # Initial mass function
                :mlow    => parse.(Float64, mlow),  # Min. mass of the IMF
                :mup     => parse.(Float64, mup),   # Max. mass of the IMF
                :Zmet    => parse.(Float64, "0." .* Zmet),  # Metallicities
                :log_age => parse.(Float64, ages),          # Stellar ages
                :Qd      => Qd,  # Number of dissociating photons per unit time
				:Qi      => Qi,  # Number of dissociating photons per unit time
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
In what follows we will use the IMF of [Chabrier2003](https://doi.org/10.1086/374879)
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
# Compute ηd and ηi
#
# age: Age [Myr]
# Z:   Metallicity [dimensionless]
##################################################################################
	
function photodissociation_efficiency(age::Float64, Z::Float64)::NTuple{2,Float64}

	# Allocate memory for the ηs
	ηd = Matrix{Float64}(undef, length(Q_ages), length(Q_metals))
	ηi = Matrix{Float64}(undef, length(Q_ages), length(Q_metals))
		
	@inbounds for (i, log_age) in pairs(Q_ages)
			
		@inbounds for (j, Zmet) in pairs(Q_metals)
				
	        sub_df = @subset(Q_df, :Zmet .== Zmet, :log_age .<= log_age)
				
	        # Set the values of the axes, with an extra point, 
			# to integrate from age 0 onwards
	        ages = [0.0, exp10.(sub_df[!, :log_age])...] .* u"yr"
				
			qd = sub_df[!, :Qd]
			Qd = [qd[1], qd...]
			ηd[i, j] = uconvert(Unitful.NoUnits, trapz(ages, Qd) * f_d)

			qi = sub_df[!, :Qi]
			Qi = [qi[1], qi...]
	        ηi[i, j] = uconvert(Unitful.NoUnits, trapz(ages, Qi) * f_i)
				
	    end
					
	end
		
	max_age = log10(age * 10^6)
		
	ifunc_ηd = linear_interpolation((Q_ages, Q_metals), ηd, extrapolation_bc=Flat())
	ifunc_ηi = linear_interpolation((Q_ages, Q_metals), ηi, extrapolation_bc=Flat())

	return ifunc_ηd(max_age, Z), ifunc_ηi(max_age, Z)
		
end;

# ╔═╡ a0294888-90cf-4e5b-a4b8-ce2c63bdae7a
# ╠═╡ skip_as_script = true
#=╠═╡
let
	set_theme!(theme_black())
	
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel=L"\mathrm{stellar \,\, age \, / \, Myr}", 
		ylabel=L"\eta_d",
		xlabelsize=32,
		ylabelsize=32,
		xticklabelsize=25,
		yticklabelsize=25,
		xscale=log10,
	)

	ages = exp10.(range(-1, 3, 30))
	ηd(age, Zs) = photodissociation_efficiency(age, Zs * Zsun)[1]

	for Zs in range(0, 2, 5)
		label = L"Z \, / \, Z_\odot = %$(Zs)"
		lines!(ax, ages, x -> ηd(x, Zs); linewidth=3, label)
	end

	axislegend(ax; position=:rb, labelsize=25)

	f
end
  ╠═╡ =#

# ╔═╡ 994f97fb-1c30-4825-9b29-35fe4ade8fb3
# ╠═╡ skip_as_script = true
#=╠═╡
let
	set_theme!(theme_black())
		
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel=L"\mathrm{stellar \,\, age \, / \, Myr}", 
		ylabel=L"\eta_i",
		xlabelsize=32,
		ylabelsize=32,
		xticklabelsize=25,
		yticklabelsize=25,
		xscale=log10,
	)

	ages = exp10.(range(-1, 3, 30))
	ηi(age, Zs) = photodissociation_efficiency(age, Zs * Zsun)[2] 

	for Zs in range(0, 2, 5)
		label = L"Z \, / \, Z_\odot = %$(Zs)"
		lines!(ax, ages, x -> ηi(x, Zs); linewidth=3, label)
	end

	axislegend(ax; position=:rb, labelsize=25)

	f
end
  ╠═╡ =#

# ╔═╡ 533b3cd0-c1f6-4ecd-b196-4ed35bf77135
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Mass recycling

There are two mass recycling parameters, $R$ which is defined as the mass fraction of a stellar population that is returned to the ISM under the instantaneous recycling hypothesis (stars under certain mass live forever, and stars above that mass die instantly), and $Z_\mathrm{SN}$ which is the fraction of the returned gas that is composed of metals (the rest is assumed to be ionized gas). 

Notice that the instantaneous recycling hypothesis can be avoided by considering the lifetimes of the stars (using empirical relations) as it is done in [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635) (sections 2.2.2 and 2.2.3). This would effectively make $R$ time dependant (the integrals below would have to be computed at each evaluation of the equations) increasing significantly the computational cost, so we avoid this and assume $R$ constant during the time scales of the problem.

A stellar yield model gives the amount (as a fraction of the stellar mass) of each element that is returned to the ISM by stars with masses between $m$ and $m + \mathrm{d}m$, so they can be used to compute

$\begin{equation}
	R = \dfrac{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} (m - m_\mathrm{rem}(m)) \, \phi(m) \mathrm{d}m}{\int_{m_\mathrm{low}}^{m_\mathrm{high}} m \, \phi(m) \mathrm{d}m} \, ,
\end{equation}$

where $\phi(m)$ is the ISM, $m_\mathrm{low}$ and $m_\mathrm{high}$ are the extremes in the mass range of the ISM, $m_\mathrm{ir}$ is the mass limit for the instantaneous recycling hypothesis, and $m_\mathrm{rem}(m)$ is the remnant stellar mass given by the yield model, [Pipino2014](https://doi.org/10.1093/mnras/stu579) ans [Ascasibar2015](https://doi.org/10.1093/mnras/stv098).

Notice that the denominator in the expression for $R$ is the total mass of the stellar population modeled by $\phi(m)$, so it is just a normalization, given that the IMF is generally defined except for a global constant.

Using the same notation we can calculate $Z_\mathrm{SN}$ as

$\begin{equation}
	Z_\mathrm{SN} = \dfrac{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} m \, f_Z \, \phi(m) \mathrm{d}m}{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} (m - m_\mathrm{rem}(m)) \, \phi(m) \mathrm{d}m} \, ,
\end{equation}$

where $f_Z$ is the fraction of the stellar mass that is returned to the ISM as metals.

Some traditional choices for the masses are $m_\mathrm{ir} = 8 \, M_\odot$, $m_\mathrm{low} = 0.08 \, M_\odot$ (the limit for hydrogen fusion), and $m_\mathrm{high} = 100 \, M_\odot$ (the order of magnitude for the upper limit of validity of the yield models, given the observational limitations).

[Ascasibar2015](https://doi.org/10.1093/mnras/stv098) got $R \approx 0.18$ and $Z_{SN} \approx 0.09$, using the yield model of [Woosley1995](https://doi.org/10.2172/115557) and the IMF of [Kroupa2001](https://doi.org/10.1046/j.1365-8711.2001.04022.x), even though it is not clear which mass limits were used.

The two previous equations were taken from [_Nucleosynthesis and Chemical Evolution of Galaxies_](https://doi.org/10.1017/CBO9780511812170) by Bernard Pagel (eq. 7.24 and 7.26), and [_Chemical Evolution
of Galaxies_](https://doi.org/10.1007/978-94-010-0967-6) by Francesca Matteucci (eq. 2.74). 

We will consider the following stellar yields models,

  -  [Woosley et al. 1995](https://doi.org/10.2172/115557)
  -  [Portinari et al. 1998](https://doi.org/10.48550/arXiv.astro-ph/9711337)
  -  [Chieff et al. 2004](https://doi.org/10.1086/392523)
  -  [Kobayashi et al. 2006](https://doi.org/10.1086/508914)
  -  [Heger et al. 2010](https://doi.org/10.1088/0004-637X/724/1/341)
  -  [Limongi et al. 2012](https://doi.org/10.1088/0067-0049/199/2/38)

compiled by [Mollá2015](https://doi.org/10.1093/mnras/stv1102), which are summarized in the following table, where
  - `model`: Stellar yield model
  - `s_Z`: Metallicity of the stellar population modeled by the IMF
  - `s_m`: Stellar mass
  - `m_rem`: Remnant mass, after stellar death
  - `zf_rem`: Fraction of the stellar mass ejected as metals to the ISM
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

The initial mass function $\phi(m)$ gives the number of stars between masses $m$ and $m + \mathrm{d}m$, for a given population of total mass $M$, following the relation

$\begin{equation}            
    M = \int_{m_\mathrm{low}}^{m_\mathrm{high}} m \, \phi(m) \, \mathrm{d}m \, ,
\end{equation}$

which allows to normalize $\phi(m)$ for the population of mass $M$, within the range $[m_\mathrm{low}, m_\mathrm{high}]$.

There are many models for $\phi(m)$, but one of the simplest is the power law

$\begin{equation}            
    \phi(m) = A \, m^{-\alpha}\, ,
\end{equation}$

where it is assumed that $[m] = M_\odot$, and $A$ is the normalization constant.

The following implementations don't have a specific choice for normalization (when posible $A = 1$), so they have to be multiplied by the correct constant if one wants a given value of $M$.
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

	R_vs_Z = linear_interpolation(sy_metals, Rs, extrapolation_bc=Flat())
	Zsn_vs_Z = linear_interpolation(sy_metals, Zsns, extrapolation_bc=Flat())

	return R_vs_Z(Z), Zsn_vs_Z(Z)
	
end;

# ╔═╡ 1d0a66d0-5791-4dc9-a5d4-7882b0e91767
# ╠═╡ skip_as_script = true
#=╠═╡
let
	set_theme!(theme_black())
	
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel=L"Z \, / \, Z_\odot", 
		ylabel=L"R",
		xlabelsize=32,
		ylabelsize=32,
		xticklabelsize=25,
		yticklabelsize=25,
	)

	metalicities = exp10.(range(-1, 0.5, 50))
	R(Zs) = recycled_fractions(Zs * Zsun)[1]

	lines!(ax, metalicities, R; linewidth=3)

	f
end
  ╠═╡ =#

# ╔═╡ 041916ac-4cb0-4630-a227-043fae52264d
# ╠═╡ skip_as_script = true
#=╠═╡
let
	set_theme!(theme_black())
	
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel=L"Z \, / \, Z_\odot", 
		ylabel=L"Zsn",
		xlabelsize=32,
		ylabelsize=32,
		xticklabelsize=25,
		yticklabelsize=25,
	)

	metalicities = exp10.(range(-1, 0.5, 50))
    Zsn(Zs) = recycled_fractions(Zs * Zsun)[2]

	lines!(ax, metalicities, Zsn; linewidth=3)

	f
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
	
	# Index in the ODE solution for each phase
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
# Ionized gas fraction:    i(t) = ρᵢ(t) / ρC --> y[1]
# Atomic gas fraction:     a(t) = ρₐ(t) / ρC --> y[2]
# Molecular gas fraction:  m(t) = ρₘ(t) / ρC --> y[3]
# Stellar fraction:        s(t) = ρₛ(t) / ρC --> y[4]
######################################################################################

function system!(dydt, ic, parameters, t)

    # Initial conditions
    i, a, m, s = ic

    # Parameters
	#
	# ρC: Total cell density                                 [mp * cm⁻³]
	# Z:  Metallicity                                        [dimensionless]
	# ηd: Photodissociation efficiency of Hydrogen molecules [dimensionless]
	# ηi: Photoionization efficiency of Hydrogen atoms       [dimensionless]
	# R:  Mass recycling fraction                            [dimensionless]
    ρC, Z, ηd, ηi, R = parameters
  
    # Auxiliary equations
	recombination   = i / τR(i, ρC)
    cloud_formation = a / τC(a, m, ρC, Z)
	sfr             = ψ(a, m, τS(s, ρC))

    # ODE system
	dydt[1] = -recombination + (ηi + R) * sfr
    dydt[2] = -cloud_formation + recombination + (ηd - ηi) * sfr
    dydt[3] = cloud_formation - (1 + ηd) * sfr
    dydt[4] = (1 - R) * sfr
	
end;

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
# J:          Matrix to save the results, it has to have size N_EQU × N_EQU
# ic:         Initial condition, [i(0), a(0), m(0), s(0)]
# parameters: Parameters for the ODEs, [ρ, Z, ηd, ηi, R]
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
# Compute the stiffness coefficient
# 
# ic:         Initial conditions, [i(0), a(0), m(0), s(0)]
# base_parms: Parameters, (ρ, Z, it)
######################################################################################

function stiffness_coefficient(
	ic::Vector{Float64},
	base_parms::NTuple{3,Float64},
)::Float64

	# Construct the parameters for the ODEs
	ρ      = base_parms[1]  # Density [cm⁻³]
	Z      = base_parms[2]  # Metallicity [dimensionless]
	it     = base_parms[3]  # Integration time [Myr]
	ηd, ηi = photodissociation_efficiency(it, Z)
	R, _   = recycled_fractions(Z)
		
	parameters = [ρ, Z, ηd, ηi, R]
	J = Matrix{Float64}(undef, N_EQU, N_EQU)
	
	# Compute the Jacobian and save it in J
	jacobian!(J, ic, parameters)
	
	# Get the norm of the real part of the non-zero eigenvalues
	eigen_values = filter(x -> x > eps(typeof(x)), eigvals(J) .|> real .|> abs)
	
	# Stiffness coefficient
	return maximum(eigen_values) / minimum(eigen_values)
	
end;

# ╔═╡ 1743b6dc-0a4e-4a02-90fe-3bc47833421a
# ╠═╡ skip_as_script = true
#=╠═╡
stiffness_coefficient(
	[0.8, 0.2, 0.0, 0.0], # Initial conditions, [i(0), a(0), m(0), s(0)]
	(1000.0, 0.01, 6.0),  # Parameters for the ODEs, (ρ, Z, it) 
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

# ╔═╡ 4607856c-7472-4131-a2ee-29f7150f5cb4
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Integration"
  ╠═╡ =#

# ╔═╡ bbb7263a-91e4-4a23-9e5f-416b6b7fcf6e
######################################################################################
# Solve the system of ODEs
#
# ic:          Initial conditions, [i(0), a(0), m(0), s(0)]
# base_params: Parameters for the ODEs, [ρC, Z]
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
    args::Tuple=(),
    kwargs::NamedTuple=(
		dense=false, 
		force_dtmin=true, 
		dtmin=(tspan[2] - tspan[1]) / 1.0e8,
		reltol=1.0e-10, 
		maxiters=1.0e8,
		verbose=false,
	),
)::Vector{Vector{Float64}}

	ρC     = base_params[1]
	Z      = base_params[2]
	ηd, ηi = photodissociation_efficiency(tspan[2], Z)
	R, _   = recycled_fractions(Z)
		
	parameters = [ρC, Z, ηd, ηi, R]

    sol = solve(ODEProblem(
		ode_function, 
		ic, 
		tspan, 
		parameters,
	), args...; kwargs...)
		
    return sol(times).u
	
end;

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataFramesMeta = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
Trapz = "592b5752-818d-11e9-1e9a-2b8ca4a44cd1"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"

[compat]
CairoMakie = "~0.10.4"
DataFrames = "~1.5.0"
DataFramesMeta = "~0.14.0"
DifferentialEquations = "~7.7.0"
Interpolations = "~0.14.7"
PlutoUI = "~0.7.50"
QuadGK = "~2.8.1"
Symbolics = "~5.3.1"
TikzPictures = "~3.4.2"
Trapz = "~2.0.3"
Unitful = "~1.12.3"
UnitfulAstro = "~1.2.0"

[extras]
CPUSummary = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "9fa7ab136d1758479978b4d76a7c780b32e0fd2d"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "3ee5c58774f4487a5bf2bb05e39d91ff5022b4cc"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.29.4"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "16b6dbc4cf7caee4e1e75c49485ec67b667098a0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.3.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cc37d689f599e8df4f464b2fa3870ff7db7492ef"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.1"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "38911c7737e123b28182d89027f4216cfc8a9da7"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.3"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4aff5fa660eb95c2e0deb6bcdabe4d9a96bc4667"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.8.18"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["Printf", "ScanByte", "TranscodingStreams"]
git-tree-sha1 = "d50976f217489ce799e366d9561d56a98a30d7fe"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.2"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "SnoopPrecompile", "SparseArrays"]
git-tree-sha1 = "6ef8fc1d77b60f41041d59ce61ef9eb41ed97a83"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.17.18"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.BoundaryValueDiffEq]]
deps = ["BandedMatrices", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase", "SparseArrays"]
git-tree-sha1 = "ed8e837bfb3d1e3157022c9636ec1c722b637318"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "2.11.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "Static"]
git-tree-sha1 = "2c144ddb46b552f72d7eafe7cc2f50746e41ea21"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.2"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA", "SnoopPrecompile"]
git-tree-sha1 = "2aba202861fd2b7603beb80496b6566491229855"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.10.4"

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
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "485193efd2176b88e6622a39a246f8c5b600e74e"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.6"

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
git-tree-sha1 = "be6ab11021cd29f0344d5c4357b163af05a48cba"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.21.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

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
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "89a9db8d28102b094992472d333674bd1a83ce2a"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.1"

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

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "aa51303df86f8626a962fccb878430cdb0a97eee"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.5.0"

[[deps.DataFramesMeta]]
deps = ["Chain", "DataFrames", "MacroTools", "OrderedCollections", "Reexport"]
git-tree-sha1 = "7f13b2f9fa5fc843a06596f1cc917ed1a3d6740b"
uuid = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
version = "0.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack"]
git-tree-sha1 = "89f3fbfe78f9d116d1ed0721d65b0b2cf9b36169"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.42.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DataStructures", "Distributions", "DocStringExtensions", "EnumX", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "Static", "StaticArraysCore", "Statistics", "Tricks", "TruncatedStacktraces", "ZygoteRules"]
git-tree-sha1 = "988bbd7283aaee5c34cd3cc09e78e7c45a931c5b"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.123.0"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "Markdown", "NLsolve", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "63b6be7b396ad395825f3cc48c56b53bfaf7e69d"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.26.1"

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArrays", "Statistics"]
git-tree-sha1 = "2c4ed3eedb87579bfe9f20ecc2440de06b9f3b89"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.16.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "a4ad7ef19d2cdc2eff57abbbe68032b1cd0bd8f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.13.0"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "JumpProcesses", "LinearAlgebra", "LinearSolve", "NonlinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "ac145e3d718157c679fc4febf2fcef73ec77b067"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.7.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "49eba9ad9f7ead780bfb7ee319f962c811c6d3b2"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.8"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "13027f188d26206b9e7b863036f87d2f2e7d013a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.87"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "698124109da77b6914f64edd696be8dccf90229e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.6.6"

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
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "8b84876e31fa39479050e2d3395c4b3b210db8b0"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "Printf", "SnoopPrecompile", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "fb7dbef7d2631e2d02c49e2750f7447648b0ec9b"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.24.0"

[[deps.ExprTools]]
git-tree-sha1 = "c1d06d129da9f55715c6c212866f5b1bddc5fa00"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.9"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "f9818144ce7c8c41edf5c4c179c684d92aa4d9fe"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.6.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "d1248fceea0b26493fd33e8e9e8c553270da03bd"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.2.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c1293a93193f0ae94be7cf338d33e162c39d8788"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "1.2.9"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "7072f1e3e5a8be51d525d64f63d3ec1287ff2790"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.11"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "03fcb1c42ec905d15b305359603888ec3e65f886"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.19.0"

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
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "cabd77ab6a6fdff49bfd24af2ebe76e6e018a2b4"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.0.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "38a92e40157100e796690421e34a11c107205c86"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.0"

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

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "1cd7f0af1aa58abc02ea1d872953a97359cb87fa"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.4"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "0eb6de0b312688f852f347171aba888658e29f20"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "303202358e38d2b01ba46844b92e48a3c238fd9e"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.6"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

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
git-tree-sha1 = "1cf1d7dcb4bc32d7b4a5add4232db3750c27ecb4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.8.0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "678d136003ed5bceaab05cf64519e3f956ffa4ba"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random", "SnoopPrecompile"]
git-tree-sha1 = "b6c3e9e1eb8dcc6fd9bc68fe08dcc7ab22710de6"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.3.4"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "734fd90dd2f920a2f1921d5388dcebe805b262dc"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.14"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "432b5b03176f8182bd6841fbfc42c718506a2d5f"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.15"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "c54b581a83008dc7f292e205f4c409ab5caa0f04"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.10"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "342f789fd041a55166764c351da1710db97ce0e0"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.6"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "36cbaebed194b292590cba2593da27b34763804a"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.8"

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
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0cb9352ef2e01574eeebdb102948a58740dcaf83"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.1.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

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
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "106b6aa272f294ba47e96bd3acbabdc0407b5c60"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.2"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "50bd271af7f6cc23be7d24c8c4804809bb5d05ae"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.6.3"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "764164ed65c30738750965d55652db9c94c59bfe"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.4.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "4a9513ad756e712177bd342ba6c022b515ed8d76"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.6"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "dd90aacbfb622f898a97c2a4411ac49101ebab8a"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.0"

[[deps.KrylovKit]]
deps = ["ChainRulesCore", "GPUArraysCore", "LinearAlgebra", "Printf"]
git-tree-sha1 = "1a5e1d9941c783b0119897d29f2eb665d876ecf3"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.6.0"

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

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "cd04158424635efd05ff38d5f55843397b7416a9"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.14.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "099e356f267354f46ba65087981a77da23a279b7"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.0"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "88b8f66b604da079a627b6fb2860d3704a6729a1"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.14"

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
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

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

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "DocStringExtensions", "EnumX", "FastLapackInterface", "GPUArraysCore", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SnoopPrecompile", "SparseArrays", "Sparspak", "SuiteSparse", "UnPack"]
git-tree-sha1 = "4a4f8cc7a59fadbb02d1852d1e0cef5dca3a9460"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.42.0"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDTypes", "SLEEFPirates", "SnoopPrecompile", "SpecialFunctions", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "defbfba8ddbccdc8ca3edb4a96a6d6fd3cd33ebd"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.157"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "MiniQhull", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Setfield", "Showoff", "SignedDistanceFields", "SnoopPrecompile", "SparseArrays", "StableHashTraits", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "74657542dc85c3b72b8a5a9392d57713d8b7a999"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.19.4"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "9926529455a331ed73c19ff06d16906737a876ed"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.3"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test", "UnicodeFun"]
git-tree-sha1 = "8f52dbaa1351ce4cb847d95568cb29e62a307d93"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.6"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.MiniQhull]]
deps = ["QhullMiniWrapper_jll"]
git-tree-sha1 = "9dc837d180ee49eeb7c8b77bb1c860452634b0d1"
uuid = "978d7f02-9e05-4691-894f-ae31a51d76ca"
version = "0.4.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "eaa98afe2033ffc0629f9d0d83961d66a021dfcc"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3295d296288ab1a0a2528feb424b854418acff57"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.2.3"

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

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "5ae7ca23e13855b3aba94550f26146c01d259267"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ArrayInterface", "DiffEqBase", "EnumX", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "a6000c813371cd3cd9cbbdf8a356fc3a97138d92"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "1.6.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "82d7c9e310fe55aa54996e6f7f94674e2a38fcb4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.9"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

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
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "a89b11f0f354f06099e4001c151dffad7ebab015"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.5"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "IfElse", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLNLSolve", "SimpleNonlinearSolve", "SimpleUnPack", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "b639e192c0422c849aeda7240382375961d0cb4b"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.50.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "f809158b27eba0c18c269cf2a2be6ed751d3e81d"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.17"

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
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "84a314e3926ba9ec66ac097e3635e270986b0f10"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.9+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "0fe4e7c4d8ff4c70bfa507f0dd96fa161b115777"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.3"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Poppler_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "e11443687ac151ac6ef6699eb75f964bed8e1faa"
uuid = "9c32591e-4766-534b-9725-b71a8799265b"
version = "0.87.0+2"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "Requires"]
git-tree-sha1 = "f739b1b3cc7b9949af3b35089931f2b58c289163"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.12"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "2e47054ffe7d0a8872e977c0d09eb4b3d162ebde"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.0.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "548793c7859e28ef026dba514752275ee871169f"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.3"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "311a2aa90a64076ea0fac2ad7492e914e6feeb81"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QhullMiniWrapper_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Qhull_jll"]
git-tree-sha1 = "607cf73c03f8a9f83b36db0b86a3a9c14179621f"
uuid = "460c41e3-6112-5d7f-b78c-b6823adb3f2d"
version = "1.0.0+1"

[[deps.Qhull_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "238dd7e2cc577281976b9681702174850f8d4cbc"
uuid = "784f63db-0788-585a-bace-daefebcd302b"
version = "8.0.1001+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "552f30e847641591ba3f39fd1bed559b9deb0ef3"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.6.1"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

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
git-tree-sha1 = "6d7bb727e76147ba18eed998700998e17b8e4911"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.4"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "140cddd2c457e4ebb0cdc7c2fd14a7fbfbdf206e"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.3"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "SnoopPrecompile", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "9088515ad915c99026beb5436d0a09cd8c18163e"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.18"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

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
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "f139e81a81e6c29c40f1971c9e5309b09c03f2c3"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.6"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "8b20084a97b004588125caebf418d8cab9e393d1"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.4"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "cda0aece8080e992f6370491b08ef3909d1c04e7"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.38"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "2436b15f376005e8790e318329560dcc67188e84"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.3"

[[deps.SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SnoopPrecompile", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces"]
git-tree-sha1 = "392d3e28b05984496af37100ded94dc46fa6c8de"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.91.7"

[[deps.SciMLNLSolve]]
deps = ["DiffEqBase", "LineSearches", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "2e1606c282fae6bd9aed4f159695774a44b9c75f"
uuid = "e9a6253c-8580-4d32-9898-8661bb511710"
version = "0.1.4"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "e61e48ef909375203092a6e83508c8416df55a83"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.2.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "77d3c4726515dca71f6d80fbb5e251088defe305"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.18"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

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

[[deps.SimpleNonlinearSolve]]
deps = ["ArrayInterface", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Reexport", "Requires", "SciMLBase", "SnoopPrecompile", "StaticArraysCore"]
git-tree-sha1 = "54c78ac3cc0343a16785adabe5bbf4063c737967"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "0.1.14"

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
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "e19ac47477c9a8fcca06dab5e5471417d5d9d723"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.31.0"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

[[deps.StableHashTraits]]
deps = ["CRC32c", "Compat", "Dates", "SHA", "Tables", "TupleTools", "UUIDs"]
git-tree-sha1 = "0b8b801b8f03a329a4e86b44c5e8a7d7f4fe10a3"
uuid = "c5dd0088-6c3f-4803-b00e-f31a60c170fa"
version = "0.3.1"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "08be5ee09a7632c32695d954a602df96a877bf0d"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.6"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "33040351d2403b84afce74dae2e22d3f5b18edcb"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.4.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "fd9a77cfd87116a27b2121c1988045f428b35a36"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.22"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

[[deps.SteadyStateDiffEq]]
deps = ["DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "04a7d0bb1c824857ba0bb0c17bc5950dccbfdd5d"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.14.0"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "073da86200349ddf4ef8bc3e3f3acd62e1d554f7"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.60.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "b3e9c174a9df77ed7b66fc0aa605def3351a0653"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.4.13"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "521a0e828e98bb69042fec1809c1b5a680eb7389"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.15"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+0"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "Reexport", "SciMLBase", "SnoopPrecompile", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "a4e8491c163d74fb92905c6443e59558f5e609a4"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.16.0"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg", "SuiteSparse_jll"]
git-tree-sha1 = "04777432d74ec5bc91ca047c9e0e0fd7f81acdb6"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.1+0"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "f8ab052bfcbdb9b48fad2c80c873aa0d0344dfe5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "5cb1f963f82e7b81305102dd69472fcd3e0e1483"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.0.5"

[[deps.Symbolics]]
deps = ["ArrayInterface", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "e23ec62c083ca8f15a4b7174331b3b8d1c511e47"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.3.1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Tectonic]]
deps = ["Pkg"]
git-tree-sha1 = "0b3881685ddb3ab066159b2ce294dc54fcf3b9ee"
uuid = "9ac5f52a-99c6-489f-af81-462ef484790f"
version = "0.8.0"

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
git-tree-sha1 = "c97f60dd4f2331e1a495527f80d242501d2f9865"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.1"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "8621f5c499a8aa4aa970b1ae381aae0ef1576966"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.4"

[[deps.TikzPictures]]
deps = ["LaTeXStrings", "Poppler_jll", "Requires", "Tectonic"]
git-tree-sha1 = "4e75374d207fefb21105074100034236fceed7cb"
uuid = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
version = "3.4.2"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f548a9e9c490030e545f72074a41edfd0e5bcdd7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.23"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "0b829474fed270a4b0ab07117dce9b9a2fa7581a"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.12"

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
git-tree-sha1 = "31eedbc0b6d07c08a700e26d31298ac27ef330eb"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.19"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "7bc1632a4eafbe9bd94cf1a784a9a4eb5e040a91"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.3.0"

[[deps.TupleTools]]
git-tree-sha1 = "3c712976c47707ff893cf6ba4354aa14db1d8938"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.3.0"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

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
git-tree-sha1 = "d5f4ec8c22db63bd3ccb239f640e895cfde145aa"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.2"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "b182207d4af54ac64cbc71797765068fdeff475d"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.64"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "977aed5d006b840e2e40c0b48984f7463109046d"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.3"

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
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

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
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

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
# ╟─64787011-b5b8-42be-b6e4-37ebc5138b3e
# ╟─14c7f574-0623-4254-b8f7-97984d32351c
# ╟─43eafb0f-08a8-4e53-9017-50b97ac48a52
# ╟─047bbd39-9cf9-4bd7-b38e-16aa505b0b08
# ╟─35e194f5-20dc-4391-b761-3696fe0bc117
# ╟─70078b44-4d66-49b9-930e-74261df8be78
# ╟─2fe0dc4c-da44-4fc8-bef8-1fa615a0fe4a
# ╟─744a9591-c7f1-496e-9bb4-47df2c8937dd
# ╟─af69ab25-0f06-4837-ac35-acbe38a4ffb1
# ╠═806db782-1734-4112-ab7b-84e03f4c342d
# ╠═1743b6dc-0a4e-4a02-90fe-3bc47833421a
# ╟─57ade87f-f93c-4b43-a737-a1b44f8af4fc
# ╟─534a1049-8de5-4b07-abec-c5a3456627c0
# ╠═86b692f1-0268-40f3-b4a2-d54c9828346d
# ╟─eaf272c7-4162-4a9a-92e3-9835c6158394
# ╟─ac553b12-4857-4cc1-8ea2-fe9e8863b429
# ╠═cf488d0e-3294-45b2-b40c-fa18062c97d2
# ╠═bc67cf22-4caa-497d-aae9-e5d1191468e2
# ╟─dc6fd12b-c821-4e20-a896-25c8aab9df94
# ╟─1d27ec35-65ca-4c94-9e8d-54d1c11e759f
# ╠═68732d91-805a-4663-9166-f8483213a8d2
# ╠═27281e53-e519-4ad0-af5d-59fb0e208534
# ╟─327fd38a-5ff6-4ac4-8d29-694272d9d46f
# ╟─897909e2-dcad-4ef6-9161-fd3654160dba
# ╠═00030fd8-a9db-4903-b2ed-21a64db30588
# ╠═d4f91aa3-183a-4abf-8f7a-7a05d4333e3a
# ╟─7e824ce1-1f82-48cc-a3c4-1acfba0e2100
# ╟─4a7eb24b-0874-49a3-9b08-4ffb6a7f0ce7
# ╠═1335861d-d539-4986-b9c1-d883a4fe8405
# ╠═f2a6676f-457a-476a-9ce7-c336aa9bf47f
# ╠═1734df7f-1309-4ebd-a021-5f75f0bb78b2
# ╟─4f7de8a3-7f59-4a7b-8980-53390e52e0d1
# ╟─3767c7f9-a0bc-467a-a20a-5e5a266111c7
# ╟─f65d84cd-ab5f-4270-98ba-568792d1fec1
# ╟─34faac11-85a2-44dc-bd8d-1a71656fccf4
# ╟─f8b02d00-ff30-480e-b5eb-e150e4678c95
# ╟─44c88ad8-a8c3-45e3-9a56-be3ce5bf66fa
# ╟─448e1dee-4628-4c14-9d6f-dc165b2e826e
# ╟─5311b7cc-7199-45c8-b5e7-20d3ceb191b7
# ╠═3e637368-6bdb-4d22-9a4a-df23c6682c2f
# ╠═ef65a096-cc2a-4ce6-a06b-8671c99ca777
# ╟─a0294888-90cf-4e5b-a4b8-ce2c63bdae7a
# ╟─994f97fb-1c30-4825-9b29-35fe4ade8fb3
# ╟─533b3cd0-c1f6-4ecd-b196-4ed35bf77135
# ╟─be85ba3b-5439-4cf3-bb14-d24d61a283c3
# ╟─b3a260b6-eb31-43a0-9fd6-60a507984319
# ╠═5ba3a0c1-6107-45a1-9b1d-5c323b9a7145
# ╠═7cbf5573-032e-4ddd-9575-f387a577c93e
# ╠═1b044783-0f5f-4321-abda-35e5b7ae67c4
# ╟─1d0a66d0-5791-4dc9-a5d4-7882b0e91767
# ╟─041916ac-4cb0-4630-a227-043fae52264d
# ╟─9666bdc8-cbc0-4757-9bd8-a76477c252eb
# ╟─ca9a233b-d3ca-4a76-a3d8-f29884ac9484
# ╠═d8bee772-3979-42cd-9e38-8df0925b4e6b
# ╟─e2e4ae4f-dcdc-4999-88f2-853378be859a
# ╠═177f8253-6c35-495b-9119-ce5e8e15cba8
# ╟─b3969810-ab25-4e91-ad5a-80560b80977e
# ╠═2620d8a6-030d-4a6f-911c-6552072ff7a1
# ╠═b76d4669-26dc-48cb-930f-5e40dd40a9f1
# ╠═69b8d934-c031-413b-9c86-3fbd64be5a4a
# ╟─35ac9289-ba53-453f-9d9e-ef3499949a98
# ╠═64e7e6aa-4265-4de2-a9c6-474d125b45cc
# ╟─4607856c-7472-4131-a2ee-29f7150f5cb4
# ╠═bbb7263a-91e4-4a23-9e5f-416b6b7fcf6e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
