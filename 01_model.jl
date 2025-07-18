### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ fed88caa-1520-41f7-adb3-785e5c9529c6
using ChaosTools, DataFrames, DataFramesMeta, DelimitedFiles, DifferentialEquations, Interpolations, LinearAlgebra, Measurements, NaNMath, PlutoUI, QuadGK, SpecialFunctions, Symbolics, TikzPictures, Trapz, Unitful, UnitfulAstro

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

# ╔═╡ 800dc762-ce0b-463f-858f-6e8eabbc26b0
# ╠═╡ skip_as_script = true
#=╠═╡
md"""

## Motivation

In the standard cosmological model, dark matter haloes collapse from small density perturbations. Gas then condenses within these haloes, and galaxies grow further via the accretion of matter and the merger of smaller substructures ([Rees1977](https://doi.org/10.1093/mnras/179.4.541), [Silk1977](https://doi.org/10.1086/154972), [White1978](https://doi.org/10.1093/mnras/183.3.341), [Blumenthal1984](https://doi.org/10.1038/311517a0)). Several physical processes act together from the early stages of galaxy formation. These processes shape galaxy properties from the earliest stages onward. Gravitational collapse, gas cooling, mass accretion, star formation (SF), chemical contamination, and various forms of feedback act together, in an interconnected way, as galaxies evolve, determining their properties ([Dalgarno1972](https://doi.org/10.1146/annurev.aa.10.090172.002111), [Efstathiou1992](https://doi.org/10.1093/mnras/256.1.43P), [White1991](https://doi.org/10.1086/170483)). The occurrence of mergers, interactions, and mass accretion might also trigger morphological transformations and variations in the dynamical, chemical, and structural properties of galactic systems.

Although we still lack a full analytical understanding of how each process shapes galaxy evolution, cosmological simulations over the past few decades have been instrumental in understanding how galaxies form and evolve within a cosmological context. Current simulations are able to produce systems that resemble the overall properties of real galaxies, including our galaxy ([Grand2017](https://doi.org/10.1093/mnras/stx071), [Springel2018](https://doi.org/10.1093/mnras/stx3304), [Libeskind2020](https://doi.org/10.1093/mnras/staa2541)).

One of the crucial ingredients in simulations of the growth and evolution of galaxies is the proper modeling of SF and the corresponding feedback. Dense gas clouds collapse and fragment—particularly in the centers of dark matter haloes and along the disks of spiral galaxies—leading to star formation. Initially triggered by gas cooling -- which determines the amount and distribution of dense gas -- SF serves as the primary source of heavy elements (metals) and provides energetic feedback to the interstellar medium (ISM). The balance between cooling and the heating resulting from supernova explosions is key in determining the SF rates (SFR) over time and space. Additionally, factors such as the circulation of gas due to feedback, accretion of fresh gas from the intergalactic medium and galactic fountains, internal instabilities, interactions and collisions influence the ISM properties in non-trivial ways. One of the most important advantages of cosmological simulations is that all these processes participating in the evolution of galaxies are naturally accounted for. However, neither SF nor stellar feedback can be included from first principles, as the scales at which they occur are not resolved. For this reason, sub-grid recipes are included to account for the impact of unresolved physics on larger scales, thereby allowing galaxy formation simulations to be run on cosmological volumes at high resolutions.
"""
  ╠═╡ =#

# ╔═╡ 03943fac-642a-4063-b043-a6b05b417c7d
# ╠═╡ skip_as_script = true
#=╠═╡
md"""

## Previous work

In the case of SF, models assuming a simple correlation between the SFR density and the gas density have been used since the pioneering works of [Katz1992](https://doi.org/10.1086/171366) and [Steinmetz1994](https://doi.org/10.48550/arXiv.astro-ph/9312010). In such models, cold dense gas transforms into star particles at a rate that is proportional to the gas density and inversely proportional to a typical timescale. Combined with effective feedback, this approach reproduces the observed slope of the SFR–gas density relation ([Schmidt1959](https://doi.org/10.1086/146614), [Kennicutt1998](https://doi.org/10.1086/305588)), as well as the SFR levels, as long as a tunable proportionality factor is introduced, known as the SF efficiency (SFE), and the typical SF timescale is the dynamical or free-fall time of the gas ([Scannapieco2006](https://doi.org/10.1111/j.1365-2966.2006.10785.x)).

During the last decade, more recent observations made it evident that the correlation between the gas density and the SFR is even stronger when molecular gas is considered ([Wong2002](https://doi.org/10.1086/339287), [Leroy2008](https://doi.org/10.1088/0004-6256/136/6/2782), [Bolatto2011](https://doi.org/10.1088/0004-637X/741/1/12), [Schruba2011](https://doi.org/10.1088/0004-6256/142/2/37), [Robertson2008](https://doi.org/10.1086/587796), [Thompson2013](https://doi.org/10.1088/0004-637x/780/2/145)), both at resolved scales ([Baker2021](https://doi.org/10.1093/mnras/stab3672)) and at integrated (i.e. galaxy-wide) scales across different redshifts ([Baker2022](https://doi.org/10.1093/mnras/stac3413)). Additionally, the typical timescale of SF may vary depending on the properties of the ISM where SF occurs ([Bigiel2008](https://doi.org/10.1088/0004-6256/136/6/2846), [Bigiel2010](https://doi.org/10.1088/0004-6256/140/5/1194), [Bolatto2011](https://doi.org/10.1088/0004-637X/741/1/12), [Leroy2013](https://doi.org/10.1088/0004-6256/146/2/19))). Motivated by these results, a number of works have developed more sophisticated SF models within simulations and semi-analytic models. These models rely on a better characterization of the ISM (either increasing numerical resolution or developing better sub-grid models), opening up the possibility to follow the content of molecular hydrogen as a function of time and to introduce $\mathrm{H_2}$-based SF laws ([Murante2010](https://doi.org/10.1111/j.1365-2966.2010.16567.x), [Molla2015](https://doi.org/10.1093/mnras/stv1102), [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635)).

In contrast to SF, modeling the effects of feedback in simulations has proven extremely complex. In part, the reason for this is that a detailed understanding of stellar feedback and its interaction with the surrounding medium has not been fully achieved, and there is also still discussion on numerical aspects of the modeling (see, e.g. [Agertz2013](https://doi.org/10.1088/0004-637X/770/1/25)). Moreover, different ways of implementing stellar feedback show significant variations in predictions for the stellar masses, morphologies, and other galaxy properties, even for a given initial condition [Scannapieco2012](https://doi.org/10.1111/j.1365-2966.2012.20993.x), although significant advances have been made in this direction ([Kim2013](https://doi.org/10.1088/0067-0049/210/1/14), [Kim2016](https://doi.org/10.3847/1538-4357/833/2/202)). The modeling of feedback is naturally linked to the modeling of SF, and thus the various numerical parameters that are included to regulate SF and feedback levels are specific to a given model and do not correspond to physically measurable parameters. Finally, additional forms of feedback including black holes and cosmic rays have also been explored as they also participate in the process of galaxy formation, but may not be as important as stellar feedback for galaxies up to the Milky Way (MW) mass scale [Naab2017](https://doi.org/10.1146/annurev-astro-081913-040019).

While feedback is what ultimately regulates the SF activity in galaxies, it is clear that an adequate description of the ISM is key to determine the conditions for SF. Various sub-resolution models which describe the multiphase structure of the ISM have been proposed, with the aim of providing a better modeling of SF in the context of simulations. [Yepes1997](https://doi.org/10.1093/mnras/284.1.235) and [Hultman1999](https://ui.adsabs.harvard.edu/abs/1999A%26A...347..769H) incorporated in Eulerian and Lagrangian codes the [McKee1977](https://doi.org/10.1086/155667) model, assuming that each gas element was composed of a hot and cold phase in pressure equilibrium. Cold clouds form in the hot medium due to thermal instability and produce stars which transfer mass back to the hot phase and can evaporate the cold clouds. [Springel2003](https://doi.org/10.1046/j.1365-8711.2003.06206.x) improved and extended these models in the $\texttt{GADGET}$ code, also adding galactic winds driven by SF as a form of feedback.

As for numerical implementations including molecular hydrogen, in [Pelupessy2006](https://doi.org/10.1086/504366) a sub-grid model for the formation of molecular gas was implemented, using a time-dependent model for the $\mathrm{HI} \rightarrow \mathrm{H_2}$ transition. The model was used for an integrated study of SF and $\mathrm{H}_2$ gas on a dwarf galaxy, in the setting of non-cosmological galaxy models. In the isolated galaxy formation simulations of [Robertson2008](https://doi.org/10.1086/587796), $\mathrm{H_2}$ was tracked assuming equilibrium abundances based on the gas metallicity. [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) implemented a model for tracking the non-equilibrium abundances of $\mathrm{H_2}$ in their adaptive mesh refinement (AMR) code, and ran cosmological simulations of low-metallicity galaxies, exploring the interactions between atomic gas, molecular gas, and dust. [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x) adapted this model to be incorporated into their Smoothed Particle Hydrodynamics (SPH) code and ran simulations of isolated dwarf galaxies in a cosmological context up to $z = 0$. [Murante2014](https://doi.org/10.1093/mnras/stu2400) implemented a molecular-based SF law in cosmological simulations of a MW-type halo. In their model, SF is linked to the molecular abundance within gas particles, that is computed using the phenomenological relation of [Blitz2006](https://doi.org/10.1086/505417) which only depends on gas pressure. [Valentini2022](https://doi.org/10.1093/mnras/stac2110) used cosmological simulations of a MW-size halo, focusing on the low-metallicity regime, to compare the predictions of two different models for computing the molecular fraction: the theoretical model of [Krumholz2009](https://doi.org/10.1088/0004-637X/699/1/850) and the phenomenological approach of [Blitz2006](https://doi.org/10.1086/505417). They found that the former model produces a clumpier ISM and a more complex $\mathrm{H_2}$ distribution, in better agreement with observations of nearby disc galaxies. [Hopkins2014](https://doi.org/10.1093/mnras/stu1738) also assumed an $\mathrm{H_2}$-regulated SF prescription for their galaxy formation simulation, computing the molecular hydrogen density with the local column density and metallicity as in [Krumholz2011](https://doi.org/10.1088/0004-637X/729/1/36).

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
"""
  ╠═╡ =#

# ╔═╡ 417aad38-6928-4b13-9286-3f11efcffb99
# ╠═╡ skip_as_script = true
#=╠═╡
md"""

## $\texttt{AREPO}$

We use the $N$-body, magnetohydrodynamics (MHD) code $\texttt{AREPO}$ to perform the simulations. For a more comprehensive explanation, see [Springel2010](https://doi.org/10.1111/j.1365-2966.2009.15715.x). $\texttt{AREPO}$ is a moving-mesh code that tracks MHD and collisionless dynamics in a cosmological setting. Gravitational forces are computed using a conventional TreePM method [Springel2005](https://doi.org/10.1111/j.1365-2966.2005.09655.x), which employs a Fast Fourier Transform technique for long-range forces and a hierarchical oct-tree algorithm [Barnes1986](https://doi.org/10.1038/324446a0) for short-range forces, in conjunction with adaptive time-stepping. $\texttt{AREPO}$ employs a dynamic unstructured mesh, constructed from a Voronoi tessellation of the simulation box. This allows for a finite-volume discretization of the MHD equations. The MHD equations are solved using a second-order Runge–Kutta integration scheme with high-accuracy least-squares spatial gradient estimators of primitive variables [Pakmor2015](https://doi.org/10.1093/mnras/stv2380). 

A distinctive feature of $\texttt{AREPO}$ is the ability to transform the mesh at any timestep through a mesh reconstruction, a capability not found in standard grid-based methods. The mesh reconstruction ensures that each cell contains a specific target mass (within a certain tolerance), meaning that areas of high density are resolved with more cells than areas of low density. Additionally, the mesh generating points move with the fluid flow, allowing $\texttt{AREPO}$ to overcome the Galilean non-invariance problem that standard Eulerian mesh codes have and significantly reduce the advection errors that appear in complex supersonic flows. The quasi-Lagrangian nature of the method makes it comparable to other Lagrangian methods such as SPH, although it eliminates several of its limitations, such as the need for artificial viscosity.

Several modules to model different physical phenomena exist within $\texttt{AREPO}$. In particular, we apply the galaxy formation model used in the Auriga Project which includes primordial and metal-line cooling, a uniform ultraviolet background field for reionization, magnetic fields ([Pakmor2014](https://doi.org/10.1088/2041-8205/783/1/L20), [Pakmor2017](https://doi.org/10.1093/mnras/stx1074), [Pakmor2018](https://doi.org/10.1093/mnras/sty2601)), energetic and chemical feedback from Type II supernovae, and mass loss/metal return owing to Type Ia supernovae and asymptotic giant branch stars ([Vogelsberger2013](https://doi.org/10.1093/mnras/stt1789), [Marinacci2013](https://doi.org/10.1093/mnras/stt2003), [Grand2017](https://doi.org/10.1093/mnras/stx071)). We do not include the active galactic nuclei feedback module, given that its effects are subdominant for a MW-mass galaxy.
"""
  ╠═╡ =#

# ╔═╡ b842e98e-34e2-40f2-84b6-c180815c2df3
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Phases and notation

We will model the interstellar medium as a multiphase structure made up of six components. Three correspond to different phases of hydrogen: ionized ($\sim 10^4 \, \mathrm{K}$), atomic ($\sim 100 \, \mathrm{K}$), and molecular ($\sim 10 \, \mathrm{K}$) gas. The other components are: metals, stars, and dust.

For every reaction that transfers mass between the phases, we consider only the dominant reaction pathway. So, although the ISM is primarily composed of hydrogen and helium, only hydrogen-related reactions are modeled, as these dominate the relevant processes. The processes involved are the photoionization of atoms, the recombination of electrons with ions, the conversion of atomic hydrogen into molecular hydrogen, and the destruction of the latter through photodissociation caused by UV light. In addition, we consider the formation of ionized gas and metals by supernovas, the formation and destruction of dust from and into metals, and the influence of the molecular gas on the SFR.

We characterize each phase by its mass fraction with respect to the total mass of the cell,

|||
|:-------------:|:------------------------------------:|
| Ionized gas   | $f_i(t) := M_i(t) / M_\mathrm{cell}$ |
| Atomic gas    | $f_a(t) := M_a(t) / M_\mathrm{cell}$ |
| Molecular gas | $f_m(t) := M_m(t) / M_\mathrm{cell}$ |
| Stars         | $f_s(t) := M_s(t) / M_\mathrm{cell}$ |
| Metals        | $f_Z(t) := M_Z(t) / M_\mathrm{cell}$ |
| Dust          | $f_d(t) := M_d(t) / M_\mathrm{cell}$ |

where $M_i$, $M_a$, $M_m$, $M_s$, $M_Z$, and $M_d$ are the corresponding masses and

$\begin{equation}
    M_\mathrm{cell} := M_i(t) + M_a(t) + M_m(t) + M_s(t) + M_Z(t) + M_d(t) \, ,
\end{equation}$

is the total cell mass.

It is important to clarify that these fractions describe the internal composition of each gas cell within our subgrid model and do not directly correspond to physical particle masses in the simulation. The stellar mass fraction, $f_s$, represents the mass of stars formed within a gas cell that has not yet been converted into a distinct star particle in the simulation. It is used to calculate the cell's instantaneous star formation rate, which determines when a new star particle is generated.

Now, we will compute the relation between the number density of the $j$ component ($n_j$) and the dimensionless fractions defined above ($f_j := M_j / M_\mathrm{cell}$).

$\begin{equation}
    n_j := \frac{N_j}{V_j} = \frac{M_j}{m_j} \, \frac{1}{V_j} = \frac{M_j}{M_\mathrm{cell}} \,  \frac{M_\mathrm{cell}}{m_j \, V_j} = f_j \, \frac{M_\mathrm{cell}}{m_j \, x_j \, V_\mathrm{cell}} = f_j \, \frac{\rho_\mathrm{cell}}{m_j \, x_j} \, ,
\end{equation}$

where the quantities are

|||
|:--------:|:---------------------------------------------------------------------:|
| $V_\mathrm{cell}$    | Volume of the cell                                        |
| $M_\mathrm{cell}$    | Mass of the cell                                          |
| $\rho_\mathrm{cell}$ | Mass density of the cell                                  |
| $N_j$    | Number of elements (e.g. atoms) of the $j$ component                  |
| $V_j$    | Volume of the $j$ component                                           |
| $M_j$    | Mass of the $j$ component                                             |
| $m_j$    | Mass of a single element (e.g. atom) of the $j$ component             |
| $f_j$    | Mass fraction of the $j$ component $(f_j := M_j / M_\mathrm{cell})$   |
| $x_j$    | Volume fraction of the $j$ component $(x_j := V_j / V_\mathrm{cell})$ |

In the context of our model, $V_\mathrm{cell}$, $M_\mathrm{cell}$, $\rho_\mathrm{cell}$, $V_j$, $m_j$, and $x_j$ are constants during the integration of the ODEs.

In the same way, we can write a relation for the mass density of the $j$ component,

$\begin{equation}
    \rho_j := \frac{M_j}{V_j} = \frac{M_j}{M_\mathrm{cell}} \, \frac{M_\mathrm{cell}}{V_j} = f_j \,  \frac{M_\mathrm{cell}}{x_j \, V_\mathrm{cell}} = f_j \, \frac{\rho_\mathrm{cell}}{x_j} \, ,
\end{equation}$

So, using these relations we can write any differential equation for the quantities $M_j$, $\rho_j$, and $n_j$ as an equation for $f_j$.

In our model, we have done only two hypotheses until now. First, that the ISM is only made up of the six components already mentioned, and second, that  $V_\mathrm{cell}$, $M_\mathrm{cell}$, $\rho_\mathrm{cell}$, $V_j$, $m_j$, and $x_j$ are constants.

For simplicity we will adopt $x_i = x_a = x_m = 1$. This implies that the ionized, atomic, and molecular hydrogen phases are assumed to be well mixed and to collectively fill the entire cell volume.
"""
  ╠═╡ =#

# ╔═╡ 4fc9bda3-1f0f-41c2-8005-5b8613271d4a
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Volume fractions in the Milky Way

Although we assume $x_i = x_a = x_m = 1.0$, for completeness we show the stimates for the volume fractions of the hydrogen phases in the MW.

Using the values from [Ferrière2001](https://doi.org/10.1103/RevModPhys.73.1031) we have
"""
  ╠═╡ =#

# ╔═╡ 3c0c1fea-1ae6-4c9d-9db1-764937cd5fbc
# ╠═╡ skip_as_script = true
#=╠═╡
let
	# Approximate MW radius
	R = (25.0u"kpc" + 30.0u"kpc") / 2.0

	# Approximate MW height
	h = (400.0u"pc" + 600.0u"pc") / 2.0

	# MW volume (as a disc)
	V = π * R^2 * h

	# Masses of the hydrogen phases in the MW
	Mm = (1.3e9u"Msun" + 2.5e9u"Msun") / 2.0
	Ma = 6.0e9u"Msun"
	Mi = 1.6e9u"Msun"

	# Number densities of the hydrogen phases in the MW
	nm = exp10((2.0 + 6.0) / 2.0) * u"cm^-3"
	na = (20.0u"cm^-3" + 50.0u"cm^-3") / 2.0
	ni = (0.2u"cm^-3" + 0.5u"cm^-3") / 2.0

	# Reference mass of an element of each hydrogen phase
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

	println("Volume fractions in the MW (Ferrière 2001):\n")
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
### As mathematical expressions
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
        (ion) edge [bend left, "$\textcolor{d_pink}{\dfrac{f_i}{\tau_\mathrm{rec}}} - \textcolor{d_blue}{S_d \, \eta_\mathrm{ion} \, \psi} - \textcolor{d_blue}{S_d \, \Gamma_\mathrm{UVB} \, f_a}$"] (atom)
        (atom) edge [bend left, "$\textcolor{d_orange}{\dfrac{f_a}{\tau_\mathrm{cond}}} - \textcolor{d_green}{e_\mathrm{diss} \, \eta_\mathrm{diss} \, \psi} - \textcolor{d_green}{e_\mathrm{diss} \, \Gamma_\mathrm{LWB} \, f_m}$"] (molecule)
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

The equations presented here describe only the gain terms for each component -- that is, how mass is added through various physical processes. To ensure mass conservation in the model, each gain term is accompanied by a corresponding loss term in the source component. For example, the formation of stars reduces the gas fractions, while ionization decreases the atomic gas fraction and increases the ionized fraction. In the full system of differential equations, each such process appears twice: once as a positive term in the equation for the receiving phase, and once with a negative sign in the equation for the source phase. This symmetric formulation guarantees that the total mass within each cell remains conserved under the assumptions of the model.

### Stars

We define $\psi(t)$ as the fractional SFR (SFR per unit of cell mass),

$\begin{equation}
    \left. \frac{\mathrm{d}}{\mathrm{d}t}f_s \right|_{\mathrm{formation}} =: \psi \, ,
\end{equation}$

### Ionized gas

The ionized component grows through the ionization of atomic gas and from the remnants of supernova explosions.

The former is produced by the radiation of newborn stars (within the cell) and from the metagalactic ultraviolet background radiation (UVB), so it can be written as

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_i\right|_{\mathrm{ionization}} = S_d \, \eta_\mathrm{ion} \, \psi + S_d \, \Gamma_\mathrm{UVB} \, f_a \, ,
\end{equation}$

where $S_d$ is the dust shielding factor, $\eta_\mathrm{ion}$ is the ionized mass rate per unit of created stellar mass, and $\Gamma_\mathrm{UVB}$ is the UVB photoionization rate (which already includes dust shielding).

The latter, under the instantaneous recycling hypothesis, can be written as

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_i\right|_{\mathrm{recycling}} = R \, (1 - Z_\mathrm{SN}) \, \psi \, ,
\end{equation}$

where $R$ is the total (ionized gas + metals) mass fraction of a stellar population that is returned to the ISM, and $Z_\mathrm{SN}$ is the fraction of the returned mass that is metals.

### Atomic gas

The atomic component grows through the dissociation of hydrogen molecules and the recombination of the ionized gas with free electrons.

The former can be written as

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_a\right|_{\mathrm{dissociation}} = e_\mathrm{diss} \, \eta_\mathrm{diss} \, \psi + e_\mathrm{diss} \, \Gamma_\mathrm{LWB} \, f_m \, ,
\end{equation}$

where $e_\mathrm{diss} = S_d \, S_\mathrm{H_2}$ is the dissociation shielding factor, $S_\mathrm{H_2}$ is the molecular self shielding factor, $\eta_\mathrm{diss}$ is the disassociated mass rate per unit of created stellar mass, and $\Gamma_\mathrm{LWB}$ is the disassociation rate due to Lyman–Werner backgound radiation.

The latter will depend on the mass of ionized gas and on the characteristic time scale of recombination ($\tau_\mathrm{rec}$), so it can be written as

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_a\right|_{\mathrm{recombination}} = \frac{f_i}{\tau_\mathrm{rec}} \, .
\end{equation}$

### Molecular gas

The molecular component gains mass mainly by the condensation of hydrogen atoms on the surface of dust grains. This process depends on the mass of atomic gas and on the characteristic time scale of condensation ($\tau_\mathrm{cond}$).

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_m\right|_{\mathrm{condensation}} = \frac{f_a}{\tau_\mathrm{cond}} \, .
\end{equation}$

### Metals

The metals grow from the destruction of dust and from the supernovas. 

The former depends on the mass of dust and on the characteristic time scale ($\tau_\mathrm{dd}$), so it can be written as

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
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_d\right|_{\mathrm{accretion}} = \frac{f_Z}{f_Z + f_d} \, \frac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m) \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 8dcc45e1-192a-4007-9d6e-6e4149b563d7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
From all the above, we can write the system of six ODEs,
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
		\dv{}{t}f_i &= - \textcolor{d_pink}{\frac{f_i}{\tau_\mathrm{rec}}} + \textcolor{d_blue}{S_d \, \eta_\mathrm{ion} \, \psi} + \textcolor{d_blue}{S_d \, \Gamma_\mathrm{UVB} \, f_a} + \textcolor{d_yellow}{R \, (1 - Z_\mathrm{SN}) \, \psi} \, , \\
		\dv{}{t}f_a &= \textcolor{d_pink}{\frac{f_i}{\tau_\mathrm{rec}}} - \textcolor{d_blue}{S_d \, \eta_\mathrm{ion} \, \psi} - \textcolor{d_blue}{S_d \, \Gamma_\mathrm{UVB} \, f_a} - \textcolor{d_orange}{\frac{f_a}{\tau_\mathrm{cond}}} + \textcolor{d_green}{e_\mathrm{diss} \, \eta_\mathrm{diss} \, \psi} + \textcolor{d_green}{e_\mathrm{diss} \, \Gamma_\mathrm{LWB} \, f_m} \, , \\
		\dv{}{t}f_m &= \textcolor{d_orange}{\frac{f_a}{\tau_\mathrm{cond}}} - \textcolor{d_green}{e_\mathrm{diss}  \, \eta_\mathrm{diss} \, \psi} - \textcolor{d_green}{e_\mathrm{diss}  \, \Gamma_\mathrm{LWB} \, f_m} - \textcolor{red}{\psi} \, , \\
		\dv{}{t}f_s &= \textcolor{red}{\psi} - \textcolor{d_yellow}{R \, \psi} \, , \\
		\dv{}{t}f_Z &= \textcolor{d_yellow}{R \, Z_\mathrm{SN} \, \psi} + \textcolor{g_red}{\frac{f_d}{\tau_\mathrm{dd}}} - \textcolor{g_green}{\frac{f_Z}{f_Z + f_d} \, \frac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m)} \, , \\
		\dv{}{t}f_d &= \textcolor{g_green}{\frac{f_Z}{f_Z + f_d} \, \frac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m)} - \textcolor{g_red}{\frac{f_d}{\tau_\mathrm{dd}}} \, ,
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
And we can explicitly check for mass conservation,

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

*  $\tau_\mathrm{dd}$: Time scale of dust destruction into metals.

*  $\eta_\mathrm{diss}$: Rate of molecular gas dissociation by stellar radiation per unit of created stellar mass.

*  $\eta_\mathrm{ion}$: Rate of atomic gas ionization by stellar radiation per unit of created stellar mass.

*  $R$: Mass fraction of a stellar population that is returned to the ISM as ionzed gas and metals.

*  $Z_\mathrm{SN}$: Fraction of the returned mass that is metals.

*  $\Gamma_\mathrm{UVB}$: UVB photoionization rate.

*  $\Gamma_\mathrm{LWB}$: LWB photodissociation rate.
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

For the formation of stars we will use the notation and definitions commonly used in the field (see [McKee2007](https://doi.org/10.1146/annurev.astro.45.051806.110602), [Krumholz2014](https://doi.org/10.1016/j.physrep.2014.02.001), and [Krumholz2019](https://doi.org/10.1146/annurev-astro-091918-104430)).

In particular, we will follow [Krumholz2005](https://doi.org/10.1086/431734), but 
assuming a only dependance on the mass of molecular gas instead of the total amount of gas, as in [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635)

$\begin{equation}
	\mathrm{SFR} := \frac{\mathrm{d}}{\mathrm{d}t} M_s = \frac{M_m}{\tau_\mathrm{star}}
\end{equation}$
where $M_m$ is the molecular mass and $\tau_\mathrm{star}$ is the characteristic timescale of star formation (defined by this very relation).

So, the fractional SFR is

$\begin{equation}
	\psi := \frac{\mathrm{SFR}}{M_\mathrm{cell}} = \frac{M_m}{\tau_\mathrm{star}} \, \frac{1}{M_\mathrm{cell}} = \frac{f_m}{\tau_\mathrm{star}} \, .
\end{equation}$

Following [Krumholz2005](https://doi.org/10.1086/431734), we write the characteristic timescale of star formation as

$\begin{equation}
    \tau_\mathrm{star} = \frac{t_\text{ff}}{\epsilon_\text{ff}} \, ,
\end{equation}$

where $\epsilon_\text{ff}$ is the star formation efficiency (in the literature is often used $\epsilon_\text{ff} \approx 0.01$ [Krumholz2007](https://doi.org/10.1086/509101), [Krumholz2014](https://doi.org/10.1016/j.physrep.2014.02.001)) and $t_\text{ff}$ is the free-fall time (which is the time for a pressure-free spherical cloud to collapse into a point due to its self-gravity)

$\begin{equation}
    t_\text{ff} = \sqrt{\frac{3\pi}{32 \, G \, \rho_\mathrm{g}}} \, ,
\end{equation}$

where $G$ is the gravitational constant and $\rho_\mathrm{g}$, in the context of our model, is the density of the corresponding cell, $\rho_\mathrm{cell}$.

There is a lot of uncertainty for the parameter $\epsilon_\text{ff}$ ([Lee2016](https://doi.org/10.3847/1538-4357/833/2/229) and [Utomo2018](https://doi.org/10.3847/2041-8213/aacf8f)). We will adopt $\epsilon_\text{ff} = 1.0$ because star formation is self-regulated in our model via the coupled evolution of the gas phases described by the ODEs, which naturally determines the effective star formation rate. We note though, that this parameter has been shown to have little influence on the global properties of simulated galaxies ([Li2018](https://doi.org/10.3847/1538-4357/aac9b8) and [Brown2022](https://doi.org/10.1093/mnras/stac1164)).

With all the previous definitions, we have

$\begin{equation}
    \tau_\mathrm{star} = \frac{1}{C_\mathrm{star} \, \sqrt{\rho_\mathrm{g}}} \, ,
\end{equation}$
where

$\begin{equation}
    C_\mathrm{star} = \sqrt{\frac{32 \, G}{3\pi}} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 0baae21f-71b5-42fd-be33-31de085928d3
begin
	const C_star  = sqrt(32u"G" / 3π)
	const cs_unit = m_u^(-1/2) * t_u^-1 * l_u^(3/2)
	const c_star  = ustrip(Unitful.NoUnits, C_star / cs_unit)

	τ_star(ρc) = 1.0 / (c_star * sqrt(ρc))
	@inline ψ(fm, ρc)  = c_star * fm * sqrt(ρc)
end;

# ╔═╡ a7c5650c-9bbc-44a1-8d51-d473864faae7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Recombination

The rate of recombination for hydrogen can be written as ([Osterbrock2006](http://www.worldcat.org/oclc/60611705), [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55), [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x), and [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635))

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

We can readily find fits for the case A recombination of $\alpha_H(T)$ (sum over all hydrogen states) in the literature ([Seaton1959](https://doi.org/10.1093/mnras/119.2.81), [Black1981](https://doi.org/10.1093/mnras/197.3.553), and [Verner1996](https://doi.org/10.1086/192284)).

Following the discussion in [Nebrin2023](https://doi.org/10.3847/2515-5172/acd37a), and assuming an optically thick cloud, we will use case B recombination (sum over all hydrogen states except the ground state). Using the values in table 2.1 of [Osterbrock2006](http://www.worldcat.org/oclc/60611705) for $T = 10^4 \, \mathrm{K}$ we have

$\begin{align}
    \alpha_H(10^4 \, \mathrm{K}) = \alpha_H = 2.6 \times 10^{-13} \, \mathrm{cm}^3 \, \mathrm{s}^{-1} \, .
\end{align}$

So, we can write

$\begin{equation}
    \frac{d}{dt} f_a \biggr\rvert_\mathrm{recomb.} = \alpha_H \, f_i^{\,2} \, \frac{\rho_\mathrm{cell}}{m_p} \, \frac{x_a}{x_i^{\,2}} = \frac{f_i}{\tau_\mathrm{rec}} \, ,
\end{equation}$

where the time scale $\tau_\mathrm{rec}$ is

$\begin{equation}
    \tau_\mathrm{rec} = \frac{m_p \, x_i^{\,2}}{\alpha_H \, f_i \, \rho_\mathrm{cell} \, x_a} = \frac{1}{C_\mathrm{rec} \, f_i \, \rho_\mathrm{cell}} \, ,
\end{equation}$

with the constant

$\begin{equation}
	C_\mathrm{rec} = \frac{\alpha_H}{m_p} \, ,
\end{equation}$

where we already took $x_i = x_a = 1$.
"""
  ╠═╡ =#

# ╔═╡ 49d1a7f7-2bf2-4472-94df-6247b9237ddd
begin
	const αH      = 2.6e-13u"cm^3 * s^-1"
	const C_rec   = αH / m_u 
	const cr_unit = m_u^-1 * t_u^-1 * l_u^3
	const c_rec   = ustrip(Unitful.NoUnits, C_rec / cr_unit)

	τ_rec(fi, ρ_cell) = 1.0 / (c_rec * fi * ρ_cell)
end;

# ╔═╡ 4c975b71-a387-49cb-90d9-fc51acefc795
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Recombination fits

For completeness we show how we can compute the case A recombination coefficient, using either a fit or an analytical formula from [Seaton1959](https://doi.org/10.1093/mnras/119.2.81), [Black1981](https://doi.org/10.1093/mnras/197.3.553), and [Verner1996](https://doi.org/10.1086/192284).
"""
  ╠═╡ =#

# ╔═╡ cea20293-a562-4b30-94c1-aabad784fbfc
# ╠═╡ skip_as_script = true
#=╠═╡
let
	#################################
	# From eq. 36 of Seaton (1959)
	#################################

	D = 5.197e-14u"cm^3*s^-1"
	λ(T) = 157890u"K" / T

	αA_seaton(T) = D * sqrt(λ(T)) * (0.4288 + 0.5 * log(λ(T)) + 0.469 * λ(T)^(-1/3))

	αA_seaton(1e4u"K")
end
  ╠═╡ =#

# ╔═╡ 5b08171f-5ee3-431c-8fc8-f666a4f6ee5d
# ╠═╡ skip_as_script = true
#=╠═╡
let
	##############################
	# From eq. 6 of Black (1981)
	##############################

	αA_black(T) = 4.36e-10u"cm^3*s^-1" * T^(-0.7573)

	αA_black(1e4)
end
  ╠═╡ =#

# ╔═╡ 6c533a69-6d71-4fea-8aad-224c0bf7d53c
# ╠═╡ skip_as_script = true
#=╠═╡
let
	###############################
	# From eq. 4 of Verner (1996)
	###############################

	a  = 7.982e-11u"cm^3*s^-1"
	b  = 0.748
	T0 = 3.148u"K"
	T1 = 7.036e5u"K"

	sqT0(T) = sqrt(T / T0)
	sqT1(T) = sqrt(T / T1)

	αA_verner(T) = a / (sqT0(T) * (1 + sqT0(T))^(1-b) * (1 + sqT1(T))^(1+b))

	αA_verner(1e4u"K")
end
  ╠═╡ =#

# ╔═╡ 2d6fcef9-df4b-4eec-be5c-a8865c3a1b76
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Condensation

From at least [Hollenbach1971a](https://doi.org/10.1086/150754) and [Hollenbach1971b](https://doi.org/10.1086/150755) we know that the rate of molecular hydrogen formation due to the condensation of atomic gas in the surface of dust grains is

$\begin{equation}
    \frac{d}{dt} n_m \biggr\rvert_\mathrm{cond.} = R_d \, n_H \, n_a  \, ,
\end{equation}$

where $R_d$ is the formation rate coefficient of $H_2$ on dust grain, $n_H$ is the hydrogen nucleus number density, and $n_a$ is the atomic hydrogen number density. $n_H$ comes from the assumption that the number density of dust grains is proportional to $n_H$ ([Hollenbach1971b](https://doi.org/10.1086/150755) Section II).

In the literature, $n_H$ is generally taken as $n_H = n_a + 2 \, n_m$, because most studies consider cold gas clouds ($T \sim 100 \, \mathrm{K}$) dominated by molecular and atomic hydrogen. In contrast, we consider atomic, molecular, and ionized gas within our gas cell; so, we will use $n_H = n_a + 2\, n_m + n_i$.

We also note that the expression for $\mathrm{d}n_m / \mathrm{d}t$ is only used in equilibrium equations in most of the early works ([Hollenbach1971a](https://doi.org/10.1086/150754), [Hollenbach1971b](https://doi.org/10.1086/150755), [Jura1974](https://doi.org/10.1086/152975), [Jura1975](https://doi.org/10.1086/153545), [Black1987](https://doi.org/10.1086/165740), [Sternberg1988](https://doi.org/10.1086/166664), and [Goldshmidt1995](https://doi.org/10.1086/175168)), while first appearing in an actual differential equation that does not assume equilibrium in [Draine1996](https://doi.org/10.1086/177689).

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
|:-----:|:-------------------------------------------------------------:|
| $1$   | eq. 3 in [Hollenbach1971b](https://doi.org/10.1086/150755)    |
| $5$   | section III in [Jura1974](https://doi.org/10.1086/152975)     |
| $3$   | section V in [Jura1975](https://doi.org/10.1086/153545)       |
| $0.9$ | eq. 6 in [Black1987](https://doi.org/10.1086/165740)          |
| $3$   | section 3 in [Goldshmidt1995](https://doi.org/10.1086/175168) |
| $0.6$ | eq. 18 in [Draine1996](https://doi.org/10.1086/177689)        |
| $3.5$ | section 4.2 in [Wolfire2008](https://doi.org/10.1086/587688)  |

Theoretical and experimental work has shown that $R_d \propto T^{1/2}$ (eq. 6 in [Black1987](https://doi.org/10.1086/165740)) and $R_d \propto n_\mathrm{dust}$ ([Hollenbach1971b](https://doi.org/10.1086/150755)). Assuming the simplest dust model $n_\mathrm{dust} \propto Z$, we have $R_d \propto T^{1/2} \, Z$ ([Pelupessy2006](https://doi.org/10.1086/504366) and [Wolfire2008](https://doi.org/10.1086/587688)).

Following previous prescriptions ([Pelupessy2006](https://doi.org/10.1086/504366), [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55), [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x), [Mollá2017](https://doi.org/10.1093/mnras/stx419), [Millan-Irigoyen2020](https://doi.org/10.1093/mnras/staa635)) we will only scale $R_d$ with the metallicity and use $\sim 100\,\mathrm{K}$ for the temperature of the cold neutral gas.

Our reference value for $R_d$ at $Z = Z_\odot$ is ([Wolfire2008](https://doi.org/10.1086/587688))

$\begin{equation}
    R_\odot := R_d(Z = Z_\odot) = 3.5 \times 10^{-17} \, \mathrm{cm^3 \, s^{-1}} \, ,
\end{equation}$

so, we have

$\begin{equation}
    R_d = Z \, \frac{R_\odot}{Z_\odot} \, .
\end{equation}$

The value of $R_\odot$ within $1 \, \mathrm{dex}$ should not affect significantly the results, but we will use an adjustable global factor (like the clumping factor $C_\rho$ in [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) and [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x)) to account for unresolved sub-cell density structures.

The fact that the final expression for $R_d$ would give $0$ when $Z = 0$ is problematic because $\tau_\mathrm{cond}$ would diverge. The problem comes from the incorrect assumption that the only conversion channel for $\mathrm{HI} \rightarrow \mathrm{H_2}$ is the condensation on the surface of dust grains. To solve this we do the replacement $Z \rightarrow Z + Z_\mathrm{eff}$ with $Z_\mathrm{eff} = 10^{-3} \, Z_\odot$, which eliminates the divergence and takes into account the molecular formation that occurs below $10^{-3} \, Z_\odot$ ([Glover2007](https://doi.org/10.1086/519445)).

So, we finally have

$\begin{align}
    \tau_\mathrm{cond} &= \frac{m_p \, x_a \, Z_\odot}{2 \, R_\odot \, C_\rho \, (Z + Z_\mathrm{eff}) \, \rho_\mathrm{cell} \, x_m \left( \dfrac{f_a}{x_a} + \dfrac{f_m}{x_m} + \dfrac{f_i}{x_i} \right)} \\
	&= \frac{1}{C_\mathrm{cond} \, (Z + Z_\mathrm{eff}) \, \rho_\mathrm{cell} \left( \dfrac{f_a}{x_a} + \dfrac{f_m}{x_m} + \dfrac{f_i}{x_i} \right)}  \\
	&= \frac{1}{C_\mathrm{cond} \, (f_Z + f_d + Z_\mathrm{eff}) \, \rho_\mathrm{cell} \, (f_i + f_a + f_m)} \, ,
\end{align}$

where we use that the total metallicity can be written as $Z_\mathrm{tot} = f_Z + f_d$ and $x_i = x_a = x_i = 1$.

$\begin{equation}
	C_\mathrm{cond} = \frac{2 \, R_\odot \, C_\rho}{m_p \, Z_\odot} \, .
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
    const Rsun    = 3.5e-17u"cm^3 * s^-1"
    const Zsun    = 0.0127
	const Zeff    = 1e-3 * Zsun
	const Cρ      = 100.0
	const C_cond  = (2 * Rsun * Cρ) / (m_u * Zsun)
	const cc_unit = m_u^-1 * t_u^-1 * l_u^3
	const c_cond  = ustrip(Unitful.NoUnits, C_cond / cc_unit)

	τ_cond(ρ_cell, Zt, fg) = 1.0 / (c_cond * ρ_cell * fg * (Zt + Zeff))
end;

# ╔═╡ 6cab6cb7-a432-40b6-9390-ad0083fe486d
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Dust destruction

For the destruction of dust we will summarize several physical processes with a single time scale as in [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635). The physical proceses in question are: sputtering,
collision with cosmic rays, supernova shock waves, and radiative torque of a powerful and anisotropic radiation field. In our model, we will ignore astration (the "consumption" of dust when stars are created).

From [Slavin2015](https://doi.org/10.1088/0004-637X/803/1/7) we have that the time scale of dust destruction is (Table 3):

|     **ISM model**     | **SN mass interval sources** |                 |
|:---------------------:|:----------------------------:|:---------------:|
|                       |           **Local**          |    **Global**   |
| _Carbonaceous grains_ |                              |                 |
|     HIM dominated     |         $1.6 \pm 0.7$        |  $1.2 \pm 0.3$  |
|      WM dominated     |         $3.2 \pm 1.4$        |  $2.6 \pm 0.7$  |
| _Silicate grains_     |                              |                 |
|     HIM dominated     |        $0.92 \pm 0.39$       | $0.72 \pm 0.20$ |
|      WM dominated     |         $2.0 \pm 0.8$        |  $1.5 \pm 0.4$  |

As we can see the time scale scale depends on the dust grain composition and on the environment: hot ionized intercloud medium (HIM) dominated values are for the assumption that the grain destruction occurs in warm clouds embedded in a hot medium, while the warm medium (WM) dominated values use the fiducial model of [Slavin2015](https://doi.org/10.1088/0004-637X/803/1/7).

We will consider this fiducial model and the local measurements for the SN mass interval, so the possible time scales are:

|    **Composition**    |      $\tau_\mathrm{dd}$      |
|:---------------------:|:----------------------------:|
| _Carbonaceous grains_ |         $3.2 \pm 1.4$        |
| _Silicate grains_     |         $2.0 \pm 0.8$        |

Finally, we will not consider different gran compositions, so we will take the weighted average of these values, which is:

$\begin{equation}
	\tau_\mathrm{dd} = 2.3 \, \mathrm{Gyr} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 38e91d0d-1020-44a8-becc-8542fd600104
begin
	times = [3.2 ± 1.4, 2.0 ± 0.8]
	destruction_time = Measurements.value(weightedmean(times)) * u"Gyr"
	
	τ_dd = ustrip(t_u, destruction_time)
	inv_τ_dd = ustrip(t_u^-1, 1.0 / destruction_time)
end;

# ╔═╡ 4b81f302-9637-4161-b745-8ac39b9e31d3
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Dust creation

Following [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635) we will model the dust growth as the accretion of metals from the gas. From [Dwek1998](https://doi.org/10.1086/305829) (eq. 32) and [Hirashita1999](https://doi.org/10.48550/arXiv.astro-ph/9903259) (eq. 2) we can write

$\begin{align}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_d(t)\right|_{\mathrm{acc}} = \frac{f_Z}{f_Z + f_d} \, \frac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m) \, ,
\end{align}$

where $\tau_\mathrm{dg}$ is the time scale of dust growth defined by this expression. Notice how we made the choice $X_\mathrm{cl} = f_a + f_m$, based on the definition of $X_\mathrm{cold}$ in [Hirashita1999](https://doi.org/10.48550/arXiv.astro-ph/9903259) eq.29: ''$X_\mathrm{cold}$ represents the mass fraction of the cold phase to the total mass of ISM''.

The formulation presented here uses $f_Z + f_d$ in the denominator to represent the total mass fraction of metals, which differs from the $f_Z$ used in the denominator of eq. 28 in [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635). The choice here is to align with the definition that the term represents the fraction of all available metal-phase mass that is not currently locked in dust ([Dwek1998](https://doi.org/10.1086/305829), [Hirashita1998](https://doi.org/10.1086/311806), and [Hirashita2019](https://doi.org/10.1093/mnras/sty2838))

$\begin{align}
	\frac{f_Z}{f_Z + f_d} = 1 - \frac{f_d}{f_Z + f_d} \, ,
\end{align}$

is the fraction of metals trap in the dust, so it is indeed f_d / (f_Z + f_d)$ in our notation.

The time scale can be estimated from [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) (eq. 19 and 20)

$\begin{equation}
	\tau_\mathrm{dg} = c_a \, A \, \frac{a}{a^0} \, \frac{Z_\odot^d}{Z} \, \frac{n_H^0}{n_H} \, \sqrt{\frac{T^0}{T}} \, \frac{S^0}{S} \, ,
\end{equation}$
where:

| Parameter   | Description                     |
|:-----------:|:-------------------------------:|
| $c_a$       | Geometric and integral factor   |
| $A$         | Normalization constant          |
| $a$         | Mean grain radius               |
| $Z$         | Metallicity                     |
| $n_H$       | Hydrogen number density         |
| $T$         | Gas temperature                 |
| $S$         | Sticking efficiency             |
| $a^0$       | Fiducial value of $a$           |
| $Z_\odot^d$ | Fiducial solar metallicity      |
| $n_H^0$     | Fiducal hydrogen number density |
| $T^0$       | Fiducial gas temperature        |
| $S^0$       | Fiducial sticking efficiency    |

Of all these parameters only $Z$ and $n_H$ will vary from cell to cell, all others are constants.
"""
  ╠═╡ =#

# ╔═╡ f525b379-b5e2-48bc-af1a-e57436f6176e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### $c_a$ - Geometric and integral factor

Notice how, in contrast with eq. 19 and 20 from [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x), we have the aditional factor $c_a$. To understand where this comes from we have to map the notation in the Hirashita et al. paper to ours.

From [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x) eq. 30 we have that

$\begin{align}
	\tau_\mathrm{grow} \equiv \tau_\mathrm{dg} \, ,
\end{align}$

Then, from eq. 33 we have

$\begin{align}
	\tau_\mathrm{dg} \equiv \tau_\mathrm{grow} \approx \frac{\langle a^3 \rangle_0}{3 \, \langle a^2 \rangle_0 \, a_0} \, \tau \, ,
\end{align}$

where $\langle a^3 \rangle_0$, $\langle a^2 \rangle_0$, and $a_0$ can be taken from Table 2, and $\tau$ from eq. 23 and eq. 24 of [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x) (even though we will use the constants of [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) as explained below).

The $\langle a \rangle$ factors comes from the distict definition of the time scales. $\tau$ in [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x) is a local, microscopic time scale (associated to the growth of a single grain), while $\tau_\mathrm{grow}$ is the time scale for the dust growth in a whole cloud. One and the other will be related by integrals of the chosen density distribution of grains (eq. 26 in [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x)).

In particular, we choose the model C and F (for our porposes they are the same), so we have

| Parameter                                       | Value                      |
|:-----------------------------------------------:|:--------------------------:|
| $\langle a^3 \rangle_0 / \langle a^2 \rangle_0$ | $0.0159 \, \mathrm{\mu m}$ |
| $a_0$                                           | $0.1 \, \mathrm{\mu m}$    |

$\begin{align}
	c_a := \frac{\langle a^3 \rangle_0}{3 \, \langle a^2 \rangle_0 \, a_0} = 0.053 \, ,
\end{align}$

The factor of 3 relates a time scale of radial growth to a time scale of mass growth (see eq. 12 and 13 of [Hirashita2019](https://doi.org/10.1093/mnras/sty2838)), as follows:

If we define $\tau^a$ and $\tau^m$ with

$\begin{align}
	\frac{\mathrm{d} a}{\mathrm{d} t} &= \xi \, \frac{a}{\tau^a} \, , \\
	\frac{\mathrm{d} m}{\mathrm{d} t} &= \xi \, \frac{m}{\tau^m} \, ,
\end{align}$

where $a$ is the radius of a grain and $m$ its mass. We can relate both time scales assuming spherical grains

$\begin{align}
	m = \frac{4}{3} \, \pi \, a^3 \, s \, ,
\end{align}$

where $s$ is the grain density.

If we consider $s = \mathrm{cte.}$ then is easy to show that

$\begin{align}
	\tau^m = \frac{\tau^a}{3} \, .
\end{align}$
"""
  ╠═╡ =#

# ╔═╡ c01cf078-c0cb-4000-bf76-a90cbeb10b89
const ca = 0.053;

# ╔═╡ 46c5bd5a-21b1-4f92-8eb1-a11ec2a0c94a
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### $A$ - Normalization constant

This parameter conceptually depens only on the physical properties of the dust grains, so is a function of the compositions (silicate vs. graphite grains), but in practice (given the particular way $\tau_\mathrm{dg}$ is modelled) its value will depend on the choice of fiducial parameters in the $\tau_\mathrm{dg}$ equation.

From [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x) we have

| Parameter   | Fiducial value             |
|:-----------:|:--------------------------:|
| $a^0$       | $0.1 \, \mathrm{\mu m}$    |
| $Z_\odot^g$ | $0.015$                    |
| $n_H^0$     | $1000 \, \mathrm{cm^{-3}}$ |
| $T^0$       | $50 \, \mathrm{K}$         |
| $S^0$       | $0.3$                      |

So the $A$ is

| Composition | $A$ fiducial value                |
|:-----------:|:---------------------------------:|
| Silicate    | $6.30 \times 10^7 \, \mathrm{yr}$ |
| Graphite    | $5.59 \times 10^7 \, \mathrm{yr}$ |

From [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) and [Hirashita2019](https://doi.org/10.1093/mnras/sty2838) we have

| Parameter   | Fiducial value             |
|:-----------:|:--------------------------:|
| $a^0$       | $0.1 \, \mathrm{\mu m}$    |
| $Z_\odot^d$ | $0.02$                     |
| $n_H^0$     | $1000 \, \mathrm{cm^{-3}}$ |
| $T^0$       | $10 \, \mathrm{K}$         |
| $S^0$       | $0.3$                      |

So the $A$ is

| Composition | $A$ fiducial value                 |
|:-----------:|:----------------------------------:|
| Silicate    | $1.61 \times 10^8 \, \mathrm{yr}$  |
| Graphite    | $0.993 \times 10^8 \, \mathrm{yr}$ |

Comparing with the values in [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x) one would think that the difference is only due to the choice of fiducial parameters in each paper. But doing the quotient between $\tau_\mathrm{dg}^{2011}$ and $\tau_\mathrm{dg}^{2012/19}$ it can be easily checked that it is not the case. The computation of $A$ is different beyond the choices of fiducial parameters in the papers.

We will follow [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) and [Hirashita2019](https://doi.org/10.1093/mnras/sty2838), so we end up with

| Parameter   | Fiducial value                   |
|:-----------:|:--------------------------------:|
| $A$         | $1.3 \times 10^8 \, \mathrm{yr}$ |
| $a^0$       | $0.1 \, \mathrm{\mu m}$          |
| $Z_\odot^d$ | $0.02$                           |
| $n_H^0$     | $1000 \, \mathrm{cm^{-3}}$       |
| $T^0$       | $10 \, \mathrm{K}$               |
| $S^0$       | $0.3$                            |

where we used the mean $A$ between grain compositions.
"""
  ╠═╡ =#

# ╔═╡ 9c5b30d5-08aa-48a8-9ae2-c3b6c432ab89
const A = 1.3e8u"yr";

# ╔═╡ ce2383dc-c34f-4b56-8501-42ce0539c95c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### $a$ - Mean grain radius

This value will be given by the grain-size distribution used. In particular, we will follow [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x) and as before adopt the C and F models

$\begin{equation}
	\langle a \rangle =  0.00167 \, \mathrm{\mu m} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ a3d1e1bf-c513-4d6b-a43b-3dab0106f1a5
begin
	const a = 0.00167u"μm"
	const a0 = 0.1u"μm"
end;

# ╔═╡ 5680e62b-973b-4a60-bb3f-8785ce07e581
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### $n_H$ - Hydrogen number density

The hydrogen number density will depend on which gas phase we consider to be around the formation region of dust (molecular, atomic, molecular + atomic, etc.).

For consistency with our previous choice we will use the neutral fraction,

$\begin{equation}
    n_H = n_a + 2 \, n_m = f_a \, \frac{\rho_\mathrm{cell}}{m_p} + 2 \, f_m \, \frac{\rho_\mathrm{cell}}{2 \, m_p} = f_n \, \frac{\rho_\mathrm{cell}}{m_p} \, ,
\end{equation}$

where $m_p$ is the proton mass and $f_n = f_a + f_m$.
"""
  ╠═╡ =#

# ╔═╡ 47bf94da-1368-41bd-ba46-c9c1f75cf44e
const nH0 = 1000.0u"cm^-3";

# ╔═╡ 9af6ce74-5678-4b48-8758-37b0d5a6f0e4
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### $T$ - Gas temperature

As before, this values will depend on which gas phase we consider. For consistency with $n_H$, we will adopt the molecular gas fiducial temperature,

$\begin{equation}
    T = 10 \, \mathrm{K} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ a9c6a292-086e-4aa0-9856-78d3ab3fbe35
begin
	const T0 = 10.0u"K"
	const T = 10.0u"K"
end;

# ╔═╡ 5890b699-b7de-47a3-bee7-1e7dd7663fbe
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### $S$ - Sticking efficiency

For the sticking efficiency we will use the fit of [Grassi2011](https://doi.org/10.1051/0004-6361/200913779) to the data of [Leitch-Devlin1985](https://doi.org/10.1093/mnras/213.2.295):

$\begin{equation}
    S(T_g, T_d) = 0.019 \, T_g \, (0.0017 \, T_d + 0.4) \, \exp(-0.007 \, T_g) \, .
\end{equation}$

where $T_g$ is the gas temperature and $T_d$ is the dust temperature, both in Kelvin.

Using the simplest option of $T_g = T_d = 10 \, \mathrm{K}$, we get

$\begin{equation}
    S = 0.073 \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ ef79392b-395c-497a-9c0e-dc2cd468f6e1
begin
	const S0 = 0.3
	const S = 0.073
end;

# ╔═╡ 8c9ab125-2acb-4732-a9bf-7838e819e4f7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Now, from

$\begin{equation}
	\tau_\mathrm{dg} = c_a \, A \, \frac{a}{a^0} \, \frac{Z_\odot^d}{Z} \, \frac{n_H^0}{n_H} \, \sqrt{\frac{T^0}{T}} \, \frac{S^0}{S} \, ,
\end{equation}$

we can write

$\begin{equation}
	\tau_\mathrm{dg} = \frac{1}{C_\mathrm{dg} \, Z \, \rho_\mathrm{cell} \, f_n} \, ,
\end{equation}$

where

$\begin{equation}
	C_\mathrm{dg} = \frac{a^0 \, \sqrt{T} \, S}{c_a \, A \, a \, Z_\odot^d \, n_H^0 \, m_p \, \sqrt{T^0} \, S^0} \, .
\end{equation}$

Notice that $Z_\odot^d$ is not the solar metallicity we use everywhere else, but the one used in [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) and [Hirashita2019](https://doi.org/10.1093/mnras/sty2838).

Using $Z = f_Z + f_d$, where we are considering that $Z$ represent the total metallicity in the gas and dust (see eq. 11 in [Hirashita2011](https://doi.org/10.1111/j.1365-2966.2011.19131.x)), the dust creation term in the equations can be rewritten as

$\begin{align}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_d(t)\right|_{\text{accretion}} &= \frac{f_Z}{f_Z + f_d} \, \frac{f_d}{\tau_\mathrm{dg}} \, f_n \\
	&= \frac{f_Z \, f_d \, f_n}{f_Z + f_d} \, C_\mathrm{dg} \, Z \, \rho_\mathrm{cell} \, f_n  \\
	&= C_\mathrm{dg} \, \frac{f_Z \, f_d \, f_n^2}{f_Z + f_d} \, (f_Z + f_d) \, \rho_\mathrm{cell} \\
	&= C_\mathrm{dg} \, f_Z \, f_d \, f_n^2 \, \rho_\mathrm{cell} \, .
\end{align}$
"""
  ╠═╡ =#

# ╔═╡ 2a39d6f8-da49-4976-9aa7-889391e55a5d
begin
	const Zdsun    = 0.02
	const C_dg     = (a0 * sqrt(T) * S) / (ca * A * a * Zdsun * nH0 * Unitful.mp * sqrt(T0) * S0)
	const cdg_unit = m_u^-1 * t_u^-1 * l_u^3
	const c_dg     = ustrip(Unitful.NoUnits, C_dg / cdg_unit)

	τ_dg(ρ_cell, fn, Zt) = 1 / (c_dg * ρ_cell * fn * Zt)
end;

# ╔═╡ ada0a462-e4cf-438d-b2a7-35d7280e955a
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Initial mass functions (IMF)

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
# Implementation of several IMFs
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

# ╔═╡ 93c15a25-51da-4019-b09c-83a95aa5999d
md"""
## Photodissociation and photoionization
"""

# ╔═╡ 43ee281f-1a16-445d-894d-23e0319b1fd0
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
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

where $\dot{N}_\mathrm{diss}$ is the number of photodissociating photons produced per unit time (in the Lyman–Werner band, $912\,\mathrm{Å}$ to $1107\,\mathrm{Å}$), $\dot{N}_\mathrm{ion}$ the number of ionizing photons produced per unit time (between $0$ and $912\,Å$), and $c_\mathrm{diss}$ and $c_\mathrm{ion}$ are the unit conversion factors (proton mass into solar mass).

The unit conversion factors are
"""
  ╠═╡ =#

# ╔═╡ c8efcb3f-cc6f-4c45-a0f3-55cdb73d5195
begin
	const c_diss = ustrip(1.0u"mp" |> u"Msun")
	const c_ion  = ustrip(1.0u"mp" |> u"Msun")
end

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
	Z(t - t') &\approx Z(t) \, ,
\end{align}$

where we are assuming that the $\mathrm{SFR}$ and $Z$ stayed approximately constant since the younger stars where born (time $t - t'$) until today (time $t$).

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

Using the values from $\texttt{PopStar}$ in [Mollá2009](https://doi.org/10.1111/j.1365-2966.2009.15160.x) we compute a table of $Q$ for six IMFs ([Salpeter1955](https://doi.org/10.1086/145971) in two mass ranges: $0.85\,\mathrm{M}_\odot$ to $120\,\mathrm{M}_\odot$ and $0.15\,\mathrm{M}_\odot$ to $100\,\mathrm{M}_\odot$, [Kroupa2001](https://doi.org/10.1046/j.1365-8711.2001.04022.x), [Ferrini1990](https://ui.adsabs.harvard.edu/abs/1990A%26A...231..391F), and [Chabrier2003](https://doi.org/10.1086/374879)), six metallicities (0.0001, 0.0004, 0.004, 0.008, 0.02, and 0.05) and for ages between $0.1\,\mathrm{Myr}$, and $15\,\mathrm{Gyr}$.
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
        r"(?<=_)(.*?)(?=_)",        # mlow [M⊙]
        r"[^_]+$",                  # mup [M⊙]
        r"(?<=z)(.*?)(?=_)",        # metallicity
        r"(?<=t)(.*)",              # log(age [yr])
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
In what follows we will use the IMF of [Chabrier2003](https://doi.org/10.1086/374879)
"""
  ╠═╡ =#

# ╔═╡ 72794eb3-7f09-4445-906a-af9756dceebd
begin
	const IMF      = "Chabrier2003"             # Initial mass function
	const Q_df     = Q_by_imf[IMF]
	const Q_ages   = unique(Q_df[!, :log_age])  # List of stellar ages
    const Q_metals = unique(Q_df[!, :Zmet])     # List of metallicities
end;

# ╔═╡ 221942a6-5e3b-46f2-ac68-ac0a55d17f04
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Young stars contribution
"""
  ╠═╡ =#

# ╔═╡ b96de0fd-a65b-463a-9c91-2504b8427dba
# ╠═╡ skip_as_script = true
#=╠═╡
let
	# Choose a specific metallicity to test the approximation
	aprox_test = @subset(Q_df, :Zmet .== 0.05)

	# Index corresponding to a stellar age of 10Myr
	idx = 39

	max_age = round(
		ustrip(u"Myr", exp10(aprox_test[!, :log_age][idx]) * u"yr");
		sigdigits=4,
	)

	Q_diss_quot = round(
		sum(aprox_test[!, :Q_diss][1:idx]) / sum(aprox_test[!, :Q_diss][idx+1:end]);
		sigdigits=4,
	)
	Q_ion_quot = round(
		sum(aprox_test[!, :Q_diss][1:idx]) / sum(aprox_test[!, :Q_diss][idx+1:end]);
		sigdigits=4,
	)

	Q_diss_frac = round((1.0 / (1.0 + (1.0 / Q_diss_quot))) * 100; sigdigits=4)
	Q_ion_frac = round((1.0 / (1.0 + (1.0 / Q_ion_quot))) * 100; sigdigits=4)

	println("The accumulated Q_diss for stars younger than $(max_age) Myr is $(Q_diss_quot) larger that the accumulated Q_diss produced by all the stars older than that. So they represent $(Q_diss_frac)% of the value of Q_diss.\n")

	println("The accumulated Q_ion for stars younger than $(max_age) Myr is $(Q_ion_quot) larger that the accumulated Q_ion produced by all the stars older than that. So they represent $(Q_ion_frac)% of the value of Q_ion.\n")
end
  ╠═╡ =#

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

# ╔═╡ 9a756404-f35b-40b9-bc88-4b7e3a7df0b6
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Shielding

Not all the photons discussed in the previous section will reach an atom or molecule, due to the shielding effect of dust and molecular gas.

Following [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) and [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x) (which are based on [Draine1996](https://doi.org/10.1086/177689) and [Glover2007](https://doi.org/10.1086/512238)) we have that the fraction of ionizing and dissociation flux that reaches the gas is

$\begin{align}
    e_\mathrm{ion} &= S_d \, , \\
    e_\mathrm{diss} &= S_d \, S_\mathrm{H_2} \, ,
\end{align}$

where $S_d$ is the dust term and $S_\mathrm{H_2}$ the molecular self-shielding term. Both can be modelled as functions of metallicity and density (eq. 8 from [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x))

$\begin{equation}
    S_d = \exp(−\sigma_{d, eff} \, Z^* \, (N_\mathrm{HI} + 2 \, N_\mathrm{H_2})) \, ,
\end{equation}$

where $\sigma_{d, eff} = 4 \times 10^{-21} \, \mathrm{cm^{-2}}$, $Z^*$ is the metallicity in solar units, and $N_i$ is the number column density of the corresponding quantity. Notice that these values are the fits from [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) which are different from the fits in [Draine1996](https://doi.org/10.1086/177689) and [Glover2007](https://doi.org/10.1086/512238).

The column density can be stimated using a characteristic column length ($h$). Following [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x) we will use the smoothing length of the corresponding gas cell (even gas cells have a force softening length in $\texttt{AREPO}$).

Using the relation between number density and fractions, we can write

$\begin{equation}
    N_\mathrm{HI} + 2 \, N_\mathrm{H_2} = (f_a + f_m) \, \frac{\rho_\mathrm{cell} \, h}{m_p} \, .
\end{equation}$

So, the dust term ends up being

$\begin{equation}
    S_d = \exp\left(C_{sd} \, (f_Z + f_d) \, (f_a + f_m) \, \rho_\mathrm{cell} \, h\right) \, ,
\end{equation}$

where we are considering that $Z$ represent the total metallicity in the gas and dust ($Z = f_Z + f_d$) and

$\begin{equation}
    C_{sd} = −\frac{\sigma_{d, eff}}{Z_\odot \, m_p} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 805ef471-1bd2-4ce2-bf50-eeb9c6b10e4d
begin
	const σdeff    = 4.0e-21u"cm^2"
	const C_sd     = -σdeff / (Zsun * u"mp")
	const csd_unit = m_u^-1 * l_u^2
	const c_sd     = ustrip(Unitful.NoUnits, C_sd / csd_unit)

	Sd(h, fa, fm, fZ, fd, ρc) = exp(c_sd * (fZ + fd) * (fa + fm) * ρc * h)
end;

# ╔═╡ 6eef3827-1af6-431d-a163-67bb595c337d
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The molecular term (eq. 9 from [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x)) is

$\begin{equation}
    S_\mathrm{H_2} = \frac{1 − \omega_\mathrm{H_2}}{(1 + x)^2} + \frac{\omega_\mathrm{H_2}}{(1 + x)^{1/2}} \, \exp(−0.00085 \, (1 + x)^{1/2}) \, ,
\end{equation}$

where $\omega_\mathrm{H_2} = 0.2$ and $x = N_\mathrm{H_2} \, / \, 5 \times 10^{14} \, \mathrm{cm^{-2}}$.

Notice that in [Draine1996](https://doi.org/10.1086/177689) and [Glover2007](https://doi.org/10.1086/512238) the first term is

$\begin{equation}
    \frac{1 − \omega_\mathrm{H_2}}{(1 + x / b_5)^2} \, ,
\end{equation}$

where $b_5 = b / 10^5 \, \mathrm{cm \, s^{-1}}$ and $b$ is the Doppler broadening parameter. It is implicit in [Gnedin2009](https://doi.org/10.1088/0004-637X/697/1/55) and [Christensen2012](https://doi.org/10.1111/j.1365-2966.2012.21628.x) that they took $b = 10^5 \, \mathrm{cm \, s^{-1}}$. We will do the same.

We can rewrite $x$ as

$\begin{equation}
    x = \frac{f_m \, \rho_\mathrm{cell} \, h}{2 \, m_p \, x_\mathrm{nf}} = C_{sh2} \, f_m \, \rho_\mathrm{cell} \, h \, ,
\end{equation}$

where $x_\mathrm{nf} = 5 \times 10^{14} \, \mathrm{cm^{-2}}$ and

$\begin{equation}
    C_{sh2} = \frac{1}{2 \, m_p \, x_\mathrm{nf}} \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ e14648ff-69f6-47c2-924c-c5ac69e200ed
begin
	const ωH2      = 0.2
	const xnf      = 5.0e14u"cm^-2"
	const exp_fac  = -8.5e-4
	const C_sh2    = 1.0 / (2u"mp" * xnf)
	const ch2_unit = m_u^-1 * l_u^2
	const c_sh2    = ustrip(Unitful.NoUnits, C_sh2 / ch2_unit)

	xp1(h, ρc, fm)    = 1.0 + c_sh2 * fm * ρc * h
	sq_xpi(h, ρc, fm) = NaNMath.sqrt(xp1(h, ρc, fm))

	Sh2(h, ρc, fm) = ((1 - ωH2) / xp1(h, ρc, fm)^2) + (ωH2 * exp(exp_fac * sq_xpi(h, ρc, fm)) / sq_xpi(h, ρc, fm))
end;

# ╔═╡ aaba63dd-8fc3-405d-ad19-8e8368e70019
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### UVB photoionization

We will consider the ionizing effect of the metagalactic ultraviolet background radiation (UVB) using the Table 3 from [Haardt2012](https://doi.org/10.1088/0004-637X/746/2/125). 

Notice how these values are for the optically thin limit (Section 9.1 of [Haardt2012](https://doi.org/10.1088/0004-637X/746/2/125)), so we have to add the dust shielding factor. And from the definition of $\Gamma_i$ (eq. 18), we have to multiply by the atomic fraction (the cross section make the expresion be _per ion_).

Then, the term for the ODEs is

$\begin{equation}
	\left. \frac{\mathrm{d}}{\mathrm{d}t}f_i(t)\right|_{\mathrm{UVB}} = S_d \, \Gamma_\mathrm{UVB} \, f_a \, .
\end{equation}$
"""
  ╠═╡ =#

# ╔═╡ 03754250-3831-41db-ac11-fe6f5a08f0c1
begin
	const UVB_TABLE_PATH = "./data/photoionization_rate_hm12.txt"
	const UVB_TABLE = readdlm(UVB_TABLE_PATH)
end;

# ╔═╡ 2d27f581-2c53-409e-afde-812e70bba8bf
##################################################################################
# Compute ΓUVB [Myr^(-1)]
#
# a: scale factor [dimensionless]
##################################################################################

function photoionization_UVB(a::Float64)::Float64

	z = (1 / a) - 1

    ΓUVB = linear_interpolation(
		UVB_TABLE[:, 1],
		UVB_TABLE[:, 2],
		extrapolation_bc=Flat(),
	)

	return ΓUVB(z)

end;

# ╔═╡ 6ddf393c-9f37-43b2-b401-9ab0e7a9e88c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### LWB photodissociation

The photodissociation of molecular hydrogen ($H_2$) by Lyman-Werner photons is a two-step process, known as the Solomon process ([Incatasciato2023](https://doi.org/10.1093/mnras/stad1008)).

1.  **Absorption:** An $H_2$ molecule absorbs a UV photon in the Lyman-Werner bands (absorption lines in the energy range of $11 \, \mathrm{eV} < E < 13.6 \, \mathrm{eV}$). This excites the molecule to a higher electronic state.

$H_2 + \gamma_{LW} \rightarrow H_2^*$

2.  **Decay:** The excited molecule ($H_2^*$) quickly decays. About 85% of the time, it decays back to a stable vibrational state of the ground electronic level, emitting a photon. However, about 15% of the time, it decays to the vibrational continuum of the ground state, meaning the molecule dissociates into two hydrogen atoms.

$H_2^* \rightarrow H + H \,\, \mathrm{(dissociation, \,\, \sim 15\% \,\, probability)}$
$H_2^* \rightarrow H_2 + \gamma \,\, \mathrm{(radiative \,\, decay, \,\, \sim 85\% \,\, probability)}$

The overall rate of this dissociation process depends on how many LW photons are available and the probability that a molecule will absorb one.

Then, the photodissosiation rate can be written as ([Abel1997](https://doi.org/10.1016/S1384-1076(97)00010-9))

$\frac{\mathrm{d}}{\mathrm{d}t} \rho_m(t) = \rho_m(t) \int_{LW \text{ bands}} \sigma_{pd}(\nu) \, \frac{F_\nu}{h \, \nu} \, \mathrm{d}\nu$

where 

 *  $F_\nu$ is the energy flux of photons per unit area per unit frequency $[\mathrm{erg \, s^{-1} \, cm^{-2} \, Hz^{-1}}]$,
 *  $h \, \nu$ is the energy of a single photon, $F_\nu / (h \, \nu)$ is therefore the number flux of photons $[\mathrm{photons \, s^{-1} \, cm^{-2} \, Hz^{-1}}]$,
 *  $\sigma_{pd}(\nu)$ is the photodissociation cross-section $[\mathrm{cm^{-2}}]$.

Assuming an homogeneous background field we can write the energy flux as 

$F_\nu = 4\pi \, J_\nu$

where $J_\nu$ is the radiation intensity and is generally parametrized as

$J_\nu(z) = J_{21}(z) \times 10^{-21} \, \mathrm{erg \, s^{-1} \, cm^{-2} \, Hz^{-1} \, sr^{-1}}$

Putting all together we can write the photodissociation rate as

$\Gamma_\mathrm{LWB}(z) = \left(\frac{4\pi}{h} \times 10^{-21}\right) J_{21}(z) \int_{LW \mathrm{bands}} \frac{\sigma_{pd}(\nu)}{\nu} \mathrm{d}\nu$

The integral part contains all the complex quantum mechanics of the $H_2$ molecule's many absorption lines. This integral has been numerically calculated and is treated as a constant value. 

From [Abel1997](https://doi.org/10.1016/S1384-1076(97)00010-9) we have (in the optically thin limit)

$\Gamma_\mathrm{LWB}(z) = (1.38 \times 10^{-12} \, \mathrm{s^{-1}}) \, J_{21}(z)$

Now, If we want to consider the shielding effects of dust and the self-shielding of molecular hydrogen we simply use the term $S_d \, S_\mathrm{H_2}$ previously computed

$k_{pd}(z) = S_d \, S_\mathrm{H_2} \, \Gamma_\mathrm{LWB}$

where $k_{pd}$ is the final photodissociation rate.

For $J_{21}(z)$ we will use the model from eq. 9 in [Incatasciato2023](https://doi.org/10.1093/mnras/stad1008)

$\log_{10} J_{21}(z) = A + B \, (1 + z) + C \, (1 + z)^2$

where $A = 2.119$, $B = -1.117 \times 10^{-1}$, and $C = -2.782 \times 10^{-3}$.
"""
  ╠═╡ =#

# ╔═╡ 5ce798b5-4366-4cc1-8d2d-881a827d6dc9
##################################################################################
# LW background radiation intensity 
# in [10^(-21) erg s^(-1) cm^(-2) Hz^(-1) sr^(-1)]
##################################################################################
begin
	const lwb_A = 2.119
	const lwb_B = −1.117e-1
	const lwb_C = −2.782e-3

	J21(z) = exp10(lwb_A + lwb_B * (1 + z) + lwb_C * (1 + z)^2)
end;

# ╔═╡ 400b6787-f2b5-40ef-be11-cc5d60d19b9f
##################################################################################
# Compute ΓLWB [Myr^(-1)]
#
# a: scale factor [dimensionless]
##################################################################################
begin
	const abel97 = ustrip(t_u^-1, 1.38e-12u"s^-1")
	
	function photodissociation_LWB(a::Float64)::Float64
	
		z = (1 / a) - 1
	
		return abel97 * J21(z)
	
	end
end;

# ╔═╡ 7df63fdf-ef41-44f2-8a55-e5c2c849029c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Mass recycling

There are two mass recycling parameters: $R$ which is defined as the mass fraction of a stellar population that is returned to the ISM under the instantaneous recycling hypothesis (stars under a certain mass live forever and stars above that mass die instantly, first explicitly used by [Schmidt1963](https://doi.org/10.1086/147553)), and $Z_\mathrm{sn}$ which is the fraction of the returned gas that is composed of metals (the rest is assumed to be ionized gas).

Notice that the instantaneous recycling hypothesis can be avoided by considering the lifetimes of the stars (using empirical relations) as it is done in [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635) (sections 2.2.2 and 2.2.3). This would introduce a time-dependence to R, as the integrals would need to be re-evaluated at each time step, significantly increasing computational cost. To avoid this, we will simply assume that $R$ is constant within the ODEs integration time scales.

A stellar yield model gives the amount of each element (as a fraction of the stellar mass) that is returned to the ISM by stars with masses between $m$ and $m + \mathrm{d}m$, so $R$ can be computed as

$\begin{equation}
	R = \dfrac{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} (m - m_\mathrm{rem}(m)) \, \phi(m) \, \mathrm{d}m}{\int_{m_\mathrm{low}}^{m_\mathrm{high}} m \, \phi(m) \, \mathrm{d}m} \, ,
\end{equation}$

where $\phi(m)$ is the ISM, $m_\mathrm{low}$ and $m_\mathrm{high}$ are the extremes in the mass range of the ISM, $m_\mathrm{ir}$ is the mass limit for the instantaneous recycling hypothesis, and $m_\mathrm{rem}(m)$ is the remnant stellar mass given by the yield model ([Pipino2014](https://doi.org/10.1093/mnras/stu579) and [Ascasibar2015](https://doi.org/10.1093/mnras/stv098)).

Notice that the denominator in the expression for $R$ is the total mass of the stellar population modeled by $\phi(m)$, so it is just a normalization factor. This is needed because the IMF is in general defined except for a global constant.

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

# ╔═╡ 137c43b7-a541-4ca0-a6d0-0ad3579b345a
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
In what follows we will use the stellar yield model of [Portinari1998](https://doi.org/10.48550/arXiv.astro-ph/9711337)
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
	const N_EQU = 6  # Number of equations
	const N_PAR = 8  # Number of parameters

	# Index of each phase in the ODE solution  matrix
	const phase_name_to_index = Dict(
		"ionized"   => 1,
		"atomic"    => 2,
		"molecular" => 3,
		"stellar"   => 4,
		"metals"    => 5,
		"dust"      => 6,
	)
end;

# ╔═╡ 9f18a8cb-93a6-4a87-b0c4-db27ab7b80bd
md"""
C_star, C_rec, C_cond, C_dg, C_sd, C_sh2, R, Z_SN, τ_dd
"""

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

# ╔═╡ 9958610a-3384-4410-beae-8f4ce657846d
#####################################################################################
# Personalized method of fma() for the computation of the jacobian
#####################################################################################
function Base.fma(x::Symbolics.Num, y::Symbolics.Num, z::Symbolics.Num)
	return x * y + z
end;

# ╔═╡ bd8743d6-8f21-413d-835a-e543926baa09
#####################################################################################
# System of ODEs, where each equation has units of time^(-1), and
#
# Ionized gas fraction:    fi(t) = Mi(t) / M_cell --> y[1]
# Atomic gas fraction:     fa(t) = Ma(t) / M_cell --> y[2]
# Molecular gas fraction:  fm(t) = Mm(t) / M_cell --> y[3]
# Stellar fraction:        fs(t) = Ms(t) / M_cell --> y[4]
# Metals fraction:         fZ(t) = MZ(t) / M_cell --> y[5]
# Dust fraction:           fd(t) = Md(t) / M_cell --> y[6]
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

	# Dust shielding
	csd_h       = c_sd * h
	dust_shield = csd_h * ρ_cell * fn * Zt
	sd          = exp(dust_shield)
	
	# Molecular self-shielding
	csh2_h  = c_sh2 * h
	xp1     = fma(csh2_h, fm * ρ_cell, 1.0)
	xp1_2   = xp1 * xp1
	sq_xp1  = NaNMath.sqrt(xp1)
	exp_xp1 = exp(exp_fac * sq_xp1)
	sh2     = ((1 - ωH2) / xp1_2) + (ωH2 * exp_xp1 / sq_xp1)

	#################
	# Net ionization
    #################

	# Recombination [Myr^(-1)]
	recomb = c_rec * fi * fi * ρ_cell

	# Stellar ionization [Myr^(-1)]
	s_ion = sd * η_ion * ψs

	# UVB photoionization [Myr^(-1)]
	uvb = sd * ΓUVB * fa
	
	net_ionization = -recomb + uvb + s_ion

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

	# Stellar dissociation [Myr^(-1)]
	e_diss = sd * sh2
	s_diss = e_diss * η_diss * ψs

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
    dydt[4] = fma(ψs, -R, ψs)
	dydt[5] = Zsn * R_ψs - net_dust_growth
	dydt[6] = net_dust_growth

end;

# ╔═╡ 80005099-7154-4306-9172-c9a168336e14
const JACOBIAN_FUNCTION = construct_jacobian(system!);

# ╔═╡ c291700e-3a84-49a7-85d6-592cfb3b1a11
#####################################################################################
# Evaluate the Jacobian
#
# J:          Matrix to save the results, it must have size N_EQU × N_EQU
# ic:         Initial condition, [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
# parameters: Parameters for the ODEs, [ρ, ΓUVB, η_diss, η_ion, R, Zsn, h]
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
# ic:          Initial conditions, [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
# base_params: Parameters for the ODEs, [ρ_cell, Z, a, h]
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
	args::Tuple=(),
    kwargs::NamedTuple=(
		dense=false,
		reltol=1.0e-10,
		verbose=false,
	),
)::Vector{Vector{Float64}}

	ρ_cell = base_params[1]
	Z      = base_params[2]
	a      = base_params[3]
	h      = base_params[4]

	ΓUVB          = photoionization_UVB(a)
	ΓLWB          = photodissociation_LWB(a)
	η_diss, η_ion = photodissociation_efficiency(tspan[2], Z)
	R, Zsn        = recycled_fractions(Z)

	parameters = [ρ_cell, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h]

    sol = solve(ODEProblem(
		ode_function,
		ic,
		tspan,
		parameters,
	), args...; kwargs...)

    return sol(times).u

end;

# ╔═╡ 8fd489b2-620a-40cd-8283-dc88b7f3584f
# ╠═╡ skip_as_script = true
#=╠═╡
md"# ODE coefficients"
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
# ic:         Initial conditions, [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
# base_parms: Parameters for the ODEs,
#             (ρ [cm⁻³], Z [dimensionless], a [dimensionless], h [cm], it [Myr])
#####################################################################################

function stiffness_ratio(
	ic::Vector{Float64},
	base_params::NTuple{5,Float64},
)::Float64

	# Construct the parameters for the ODEs
	ρ  = base_params[1]  # Density [cm⁻³]
	Z  = base_params[2]  # Arepo metallicity [dimensionless]
	a  = base_params[3]  # Scale factor [dimensionless]
	h  = base_params[4]  # Column height [cm]
	it = base_params[5]  # Integration time [Myr]

	ΓUVB          = photoionization_UVB(a)
	ΓLWB          = photodissociation_LWB(a)
	η_diss, η_ion = photodissociation_efficiency(it, Z)
	R, Zsn        = recycled_fractions(Z)

	parameters = [ρ, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h]
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
	[0.7, 0.27, 0.0, 0.0, 0.02, 0.01],
	# (ρ [cm⁻³], Z [dimensionless], a [dimensionless], h [cm], it [Myr])
	(1.0, Zsun, 1.0, 3.0e20, 5.0),
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
# ╠═╡ skip_as_script = true
#=╠═╡
#####################################################################################
# Compute the Lyapunov spectrum
#
# ic:         Initial conditions, [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
# base_parms: Parameters for the ODEs,
#             (ρ [cm⁻³], Z [dimensionless], a [dimensionless], h [cm], it [Myr])
#####################################################################################

function lyapunov_spectrum(
	ic::Vector{Float64},
	base_params::NTuple{5,Float64},
)::Vector{Float64}

	# Construct the parameters for the ODEs
	ρ  = base_params[1]  # Density [cm⁻³]
	Z  = base_params[2]  # Arepo metallicity [dimensionless]
	a  = base_params[3]  # Scale factor [dimensionless]
	h  = base_params[4]  # Column height [cm]
	it = base_params[5]  # Integration time [Myr]

	ΓUVB          = photoionization_UVB(a)
	ΓLWB          = photodissociation_LWB(a)
	η_diss, η_ion = photodissociation_efficiency(it, Z)
	R, Zsn        = recycled_fractions(Z)

	parameters = [ρ, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h]
	ds = CoupledODEs(system!, ic, parameters)

	# Compute the Lyapunov spectrum
	return lyapunovspectrum(ds, 100; Δt=1e-4)

end;
  ╠═╡ =#

# ╔═╡ b1e67fcd-f70b-4dd2-87e5-74bef1fa201e
# ╠═╡ skip_as_script = true
#=╠═╡
lyapunov_spectrum(
	# [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
	[0.7, 0.27, 0.0, 0.0, 0.02, 0.01],
	# (ρ [cm⁻³], Z [dimensionless], a [dimensionless], h [cm], it [Myr])
	(30.0, Zsun, 1.0, 3.0e20, 10.0),
)
  ╠═╡ =#

# ╔═╡ d1f3d076-b9d6-4b92-b25c-4155b5fc464b
# ╠═╡ skip_as_script = true
#=╠═╡
lyapunov_spectrum(
	# [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
	[0.7, 0.27, 0.0, 0.0, 0.02, 0.01],
	# (ρ [cm⁻³], Z [dimensionless], a [dimensionless], h [cm], it [Myr])
	(1000.0, Zsun, 1.0, 3.0e20, 10.0)
)
  ╠═╡ =#

# ╔═╡ f336b083-99c7-4345-8ada-78b166a98abe
# ╠═╡ skip_as_script = true
#=╠═╡
md"The model, for a typical set of initial conditions, is not chaotic (maximum lyapunov exponent < 0)."
  ╠═╡ =#

# ╔═╡ 24920a1b-355c-43e1-aff2-028e6f25584c
md"# Extras"

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

# ╔═╡ 74c82c98-f233-4e72-8f74-17bdfeddf884
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Dust initial condition"
  ╠═╡ =#

# ╔═╡ ab1b7695-3cd6-4e67-9cfd-fbaf2f4d1d15
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Schematic diagram of each phase contribution to the initial condition.
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
Following eq. 18 of [Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x) we have that the inital dust mass density is

$\begin{equation}
	\rho_d(0) = f_d \, \rho_\mathrm{cell} = \frac{m_\mathrm{X_d}}{f_\mathrm{X_d}}  \, (1 - \xi_0) \, \frac{Z}{Z_\odot} \, \left( \mathrm{\frac{X_d}{H}} \right)_\odot \, n_H \, ,
\end{equation}$

where $f_\mathrm{X_d}$ is the mass fraction in the dust of element $\mathrm{X_d}$, $m_\mathrm{X_d}$ is the atomic mass of element $\mathrm{X_d}$, $\left( \mathrm{X_d / H} \right)_\odot$ is the solar abundance of element $\mathrm{X_d}$, $n_H$ is the number density of hydrogen, and $\xi_0$ is the initial value of

$\begin{equation}
	\xi(t) = \frac{n_\mathrm{X_d}^\mathrm{gas}}{n_\mathrm{X_d}^\mathrm{tot}} \, ,
\end{equation}$

where $n_\mathrm{X_d}^\mathrm{gas}$ is the number density in the gas of element $\mathrm{X_d}$ and $n_\mathrm{X_d}^\mathrm{tot}$ is the number density in the gas and dust of element $\mathrm{X_d}$.

Using

|  Species | $\mathrm{X_d}$ | $f_\mathrm{X_d}$ | $m_\mathrm{X_d}$ | $(\mathrm{X_d / H})_\odot$ |
|:--------:|:--------------:|:----------------:|:----------------:|:--------------------------:|
| Silicate |       Si       |      $0.166$     |      $28.1$      |    $3.55 \times 10^{−5}$   |
| Graphite |        C       |        $1$       |       $12$       |    $3.63 \times 10^{−4}$   |

where $m_\mathrm{X_d}$ is in daltons (Da), and with the fiducial value $\xi_0 = 0.3$ ([Hirashita2012](https://doi.org/10.1111/j.1365-2966.2012.20702.x)), we can write

$\begin{equation}
    f_d = C_\mathrm{X_d} \, Z \, f_n \, ,
\end{equation}$

where

$\begin{equation}
    C_\mathrm{X_d} = \frac{m_\mathrm{X_d}}{f_\mathrm{X_d}}  \, (1 - \xi_0) \, \frac{1}{Z_\odot \, m_p} \, \left( \mathrm{\frac{X_d}{H}} \right)_{\!\odot} \, ,
\end{equation}$

where we used

$\begin{align}
    n_H &= n_a + 2 \, n_m = f_n \, \frac{\rho_\mathrm{cell}}{m_p} \, , \\
	Z_\odot &= 0.0127 \, .
\end{align}$

Using the values from the table we have

|  Species | $\mathrm{X_d}$ | $C_\mathrm{X_d}$ |
|:--------:|:--------------:|:----------------:|
| Silicate |       Si       |      $0.329$     |
| Graphite |        C       |      $0.238$     |

So our fiducial value is the average $C_\mathrm{X_d} = 0.283$.
"""
  ╠═╡ =#

# ╔═╡ d7cc8a66-220a-4e9b-b42e-a0e085ed3a0f
begin
	const ξ0 = 0.3

	C_xd(mx, fx, XH) = (mx / fx) * (1.0 - ξ0) * XH / (Zsun * 1.0u"mp")

	silicate = C_xd(28.1u"u", 0.166, 3.55e-5)
	graphite = C_xd(12.0u"u", 1.0, 3.63e-4)

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
As defined in [Krumholz2011](https://iopscience.iop.org/article/10.1088/0004-637X/745/1/69) the efficiency per free-fall time is

$\begin{equation}
	\epsilon_\mathrm{ff} = \mathrm{SFR} \, \frac{t_\mathrm{ff}}{M_g} \, ,
\end{equation}$
where $M_g$ is the mass of gas in consideration, $\mathrm{SFR}$ the star formation rate (in the case of the simulations, the one used in the probability calculations: $\mathrm{SFR}_p$), and

$\begin{equation}
	t_\mathrm{ff} = \sqrt{\frac{3 \, \pi}{32 \, G \, \rho}} \, .
\end{equation}$

So we can write

$\begin{equation}
	\epsilon_\mathrm{ff} = c \, \frac{\mathrm{SFR}_p}{M_\mathrm{cell} \, \sqrt{\rho_\mathrm{cell}}} \, .
\end{equation}$
where $c = \sqrt{3 \, \pi / 32 \, G}$.
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
	\frac{\mathrm{d}}{\mathrm{d}t} f_i &= - \frac{f_i}{\tau_\mathrm{rec}} + S_d \, \eta_\mathrm{ion} \, \psi + S_d \, \Gamma_\mathrm{UVB} \, f_a + R \, (1 - Z_\mathrm{SN}) \, \psi \, , \\
	\frac{\mathrm{d}}{\mathrm{d}t} f_a &= \frac{f_i}{\tau_\mathrm{rec}} - S_d \, \eta_\mathrm{ion} \, \psi - S_d \, \Gamma_\mathrm{UVB} \, f_a - \frac{f_a}{\tau_\mathrm{cond}} + e_\mathrm{diss} \, \eta_\mathrm{diss} \, \psi + e_\mathrm{diss} \, \Gamma_\mathrm{LWB} \, f_m \, , \\
	\frac{\mathrm{d}}{\mathrm{d}t} f_m &= \frac{f_a}{\tau_\mathrm{cond}} - e_\mathrm{diss}  \, \eta_\mathrm{diss} \, \psi - e_\mathrm{diss} \, \Gamma_\mathrm{LWB} \, f_m - \psi \, , \\
	\frac{\mathrm{d}}{\mathrm{d}t} f_s &= \psi - R \, \psi \, , \\
	\frac{\mathrm{d}}{\mathrm{d}t} f_Z &= R \, Z_\mathrm{SN} \, \psi + \frac{f_d}{\tau_\mathrm{dd}} - \frac{f_Z}{f_Z + f_d} \, \frac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m) \, , \\
	\frac{\mathrm{d}}{\mathrm{d}t} f_d &= \frac{f_Z}{f_Z + f_d} \, \frac{f_d}{\tau_\mathrm{dg}} \, (f_a + f_m) - \frac{f_d}{\tau_\mathrm{dd}} \, ,
\end{align}$
and set the derivatives to $0$.
"""
  ╠═╡ =#

# ╔═╡ 245ad10f-4c70-4585-b2e0-dea4cf5d19b8
#####################################################################################
# Copy of the ODE system, for testing the equilibrium condition
#
# fractions:   Gas fractions, [fi, fa, fm, fs, fZ, fd]
# parameters: Parameters, [ρ_cell, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h]
#####################################################################################

function odes(
	fractions::Vector{Float64}, 
	parameters::Vector{Float64},
)::Vector{Float64}

	dydt = Vector{Float64}(undef, N_EQU)

    system!(dydt, fractions, parameters, 1.0)

	return dydt

end;

# ╔═╡ e8c51ae9-ded4-450a-80f9-c484c2b991d1
#####################################################################################
# Compute the derivatives of the ODE system
#
# fractions:   Gas fractions, [fi, fa, fm, fs, fZ, fd]
# base_params: Parameters for the ODEs, 
#              [ρ [cm⁻³], Z [dimensionless], a [dimensionless], h [cm], it [Myr]]
#####################################################################################

function equilibrium(
    fractions::Vector{Float64},
	base_params::Vector{Float64};
)::Vector{Float64}

	ρ_cell = base_params[1]
	Z      = base_params[2]
	a      = base_params[3]
	h      = base_params[4]
	it     = base_params[5]

	ΓUVB          = photoionization_UVB(a)
	ΓLWB          = photodissociation_LWB(a)
	η_diss, η_ion = photodissociation_efficiency(it, Z)
	R, Zsn        = recycled_fractions(Z)

	parameters = [ρ_cell, ΓUVB, ΓLWB, η_diss, η_ion, R, Zsn, h]

    return odes(fractions, parameters)

end;

# ╔═╡ c96ea0ad-47c3-469e-a2f7-0743e681a6d1
#####################################################################################
# Test the equilibrium of the ODE system
#
# fractions:   Gas fractions, [fi, fa, fm, fs, fZ, fd]
# base_params: Parameters for the ODEs, 
#              [ρ [cm⁻³], Z [dimensionless], a [dimensionless], h [cm], it [Myr]]
# limit: Equilibrium criteria.
#####################################################################################

function test_equilibrium(
    fractions::Vector{Float64},
	base_params::Vector{Float64};
	limit::Float64=1e-10,
)::Vector{Bool}

	time_scale = base_params[5]

    return (equilibrium(fractions, base_params) ./ (1 / time_scale)) .< limit

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
Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
Trapz = "592b5752-818d-11e9-1e9a-2b8ca4a44cd1"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"

[compat]
ChaosTools = "~3.3.1"
DataFrames = "~1.7.0"
DataFramesMeta = "~0.15.4"
DifferentialEquations = "~7.16.1"
Interpolations = "~0.16.1"
Measurements = "~2.14.0"
NaNMath = "~1.1.3"
PlutoUI = "~0.7.68"
QuadGK = "~2.11.2"
SpecialFunctions = "~2.5.1"
Symbolics = "~6.44.0"
TikzPictures = "~3.5.0"
Trapz = "~2.0.3"
Unitful = "~1.23.1"
UnitfulAstro = "~1.2.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "6dc642e68a3c77c0cb6e18f94478a445abaa9140"

[[deps.ADTypes]]
git-tree-sha1 = "be7ae030256b8ef14a441726c4c37766b90b93a3"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.15.0"
weakdeps = ["ChainRulesCore", "ConstructionBase", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
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
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

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
git-tree-sha1 = "f7817e2e585aa6d924fd714df1e2a84be7896c60"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.3.0"
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
git-tree-sha1 = "9606d7832795cbef89e06a550475be300364a8aa"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.19.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
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
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "4e25216b8fea1908a0ce0f5d87368587899f75be"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.11.1"
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
git-tree-sha1 = "e35c672b239c5105f597963c33e740eeb46cf0ab"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.9.4"

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"
    CliqueTreesExt = "CliqueTrees"

    [deps.BandedMatrices.weakdeps]
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

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
git-tree-sha1 = "9b302ba0af3e17e8d468ae95af13415016be8ab0"
uuid = "56b672f2-a5fe-4263-ab2d-da677488eb3a"
version = "1.11.0"

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
git-tree-sha1 = "a9014924595b7a2c1dd14aac516e38fa10ada656"
uuid = "70df07ce-3d50-431d-a3e7-ca6ddb60ac1e"
version = "1.3.0"
weakdeps = ["ChainRulesCore", "ForwardDiff"]

    [deps.BracketingNonlinearSolve.extensions]
    BracketingNonlinearSolveChainRulesCoreExt = ["ChainRulesCore", "ForwardDiff"]
    BracketingNonlinearSolveForwardDiffExt = "ForwardDiff"

[[deps.BranchAndPrune]]
deps = ["AbstractTrees"]
git-tree-sha1 = "9a97232c3aab366fc4408ddc2239939b4cad0179"
uuid = "d3bc4f2e-91e6-11e9-365e-cd067da536ce"
version = "0.2.1"

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
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "5a97e67919535d6841172016c9530fd69494e5ec"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.6"

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

[[deps.Chain]]
git-tree-sha1 = "9ae9be75ad8ad9d26395bf625dea9beac6d519f1"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.6.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "06ee8d1aa558d2833aa799f6f0b31b30cada405f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.2"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChaosTools]]
deps = ["Combinatorics", "DSP", "DataStructures", "Distances", "Distributions", "DynamicalSystemsBase", "IntervalRootFinding", "LinearAlgebra", "LombScargle", "Neighborhood", "Optim", "ProgressMeter", "Random", "Reexport", "Roots", "SpecialFunctions", "Statistics", "StatsBase"]
git-tree-sha1 = "67bf00462e4b46235ba83161cc83e1cda7bd5cb8"
uuid = "608a59af-f2a3-5ad4-90b4-758bdf3122a7"
version = "3.3.1"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

    [deps.ColorTypes.weakdeps]
    StyledStrings = "f489334b-da3d-4c2e-b8f0-e476e12c162b"

[[deps.Combinatorics]]
git-tree-sha1 = "8010b6bb3388abe68d95743dcbea77650bb2eddf"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.3"

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
git-tree-sha1 = "3a3dfb30697e96a440e4149c8c51bf32f818c0f3"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.17.0"
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

[[deps.DSP]]
deps = ["Bessels", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "5989debfc3b38f736e69724818210c67ffee4352"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.8.4"
weakdeps = ["OffsetArrays"]

    [deps.DSP.extensions]
    OffsetArraysExt = "OffsetArrays"

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
git-tree-sha1 = "21a4335f249f8b5f311d00d5e62938b50ccace4e"
uuid = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
version = "0.15.4"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqRosenbrock", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack", "SymbolicIndexingInterface"]
git-tree-sha1 = "8b416f6b1f9ef8df4c13dd0fe6c191752722b36f"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.53.1"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "FastPower", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "e9b34e0eb3443492f396c97e7fed08630752a4f2"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.177.2"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDistributionsExt = "Distributions"
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
deps = ["ConcreteStructs", "DataStructures", "DiffEqBase", "DifferentiationInterface", "Functors", "LinearAlgebra", "Markdown", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "80a782f3e65d4900dcf5f2cb71f5e19d9459c04a"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "4.8.0"

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "516d553f5deee7c55b2945b5edf05b6542837887"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.24.1"

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
git-tree-sha1 = "afdc7dfee475828b4f0286d63ffe66b97d7a3fa7"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.16.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "f620da805b82bec64ab4d5f881c7592c82dbc08a"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.3"

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
git-tree-sha1 = "3e6d038b77f22791b8e3472b7c633acea1ecac06"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.120"

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
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "a7e9f13f33652c533d49868a534bfb2050d1365f"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.15"

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
git-tree-sha1 = "98c4bb95af37e5d980129261fdd6dab0392c6607"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.2"

[[deps.DynamicalSystemsBase]]
deps = ["ForwardDiff", "LinearAlgebra", "OrdinaryDiffEqTsit5", "Reexport", "Roots", "SciMLBase", "SparseArrays", "StateSpaceSets", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "fe3a7586545e5e1bb1af4de840cff53395e225fa"
uuid = "6e36e845-645a-534a-86f2-f5d4aa5a06b4"
version = "3.14.0"
weakdeps = ["StochasticDiffEq"]

    [deps.DynamicalSystemsBase.extensions]
    StochasticSystemsBase = "StochasticDiffEq"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.EnzymeCore]]
git-tree-sha1 = "8272a687bca7b5c601c0c24fc0c71bff10aafdfd"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.12"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d55dffd9ae73ff72f1c0482454dcf2ec6c6c4a63"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.5+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "cae251c76f353e32d32d76fae2fea655eab652af"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.27.0"
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

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "797762812ed063b9b94f6cc7742bc8883bb5e69e"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.9.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FastAlmostBandedMatrices]]
deps = ["ArrayInterface", "ArrayLayouts", "BandedMatrices", "ConcreteStructs", "LazyArrays", "LinearAlgebra", "MatrixFactorizations", "PrecompileTools", "Reexport"]
git-tree-sha1 = "9482a2b4face8ade73792c23a54796c79ed1bcbf"
uuid = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
version = "0.1.5"

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
git-tree-sha1 = "fd923962364b645f3719855c88f7074413a6ad92"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "1.0.2"

[[deps.FastPower]]
git-tree-sha1 = "5f7afd4b1a3969dc34d692da2ed856047325b06e"
uuid = "a4df4552-cc26-4903-aec0-212e50a0e84b"
version = "1.1.3"

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
git-tree-sha1 = "f089ab1f834470c525562030c8cfde4025d5e915"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.27.0"

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
git-tree-sha1 = "301b5d5d731a0654825f1f2e906990f7141a106b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.16.0+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"
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

[[deps.Functors]]
deps = ["Compat", "ConstructionBase", "LinearAlgebra", "Random"]
git-tree-sha1 = "60a0339f28a233601cb74468032b5c302d5067de"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.5.2"

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
git-tree-sha1 = "f88e0ba1f6b42121a7c1dfe93a9687d8e164c91b"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.5"

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
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "c5abfa0ae0aaee162a3fbb053c13ecda39be545b"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.0"

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
git-tree-sha1 = "8594fac023c5ce1ef78260f24d1ad18b4327b420"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.4"

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
git-tree-sha1 = "0f14a5456bdc6b9731a5682f439a672750a09e48"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.0.4+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "f2905febca224eade352a573e129ef43aa593354"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.1"
weakdeps = ["ForwardDiff", "Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Random", "RoundingEmulator"]
git-tree-sha1 = "79342df41c3c24664e5bf29395cfdf2f2a599412"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.36"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalRootFinding]]
deps = ["BranchAndPrune", "ForwardDiff", "IntervalArithmetic", "LinearAlgebra", "Reexport", "StaticArrays"]
git-tree-sha1 = "68c9d23b092424df6b66e06cd241d2709f1b430e"
uuid = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
version = "0.6.0"

[[deps.IntervalSets]]
git-tree-sha1 = "5fbb102dcb8b1a858111ae81d56682376130517d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.11"
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
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

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
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eac1206917768cb54957c65a615460d87b455fc1"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.1+0"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "f8da88993c914357031daf0023f18748ff473924"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.16.1"
weakdeps = ["FastBroadcast"]

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "b94257a1a8737099ca40bc7271a8b374033473ed"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.10.1"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

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

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cd10d2cc78d34c0e2a3a36420ab607b611debfbb"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.7"

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
git-tree-sha1 = "866ce84b15e54d758c11946aacd4e5df0e60b7a3"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.6.1"

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
git-tree-sha1 = "a31572773ac1b745e0343fe5e2c8ddda7a37e997"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "321ccef73a96ba828cd51f2ab5b9f917fa73945a"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.0+0"

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "97d502765cc5cf3a722120f50da03c2474efce04"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.4"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "4adee99b7262ad2a1a4bbbc59d993d24e55ea96f"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.4.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "GPUArraysCore", "InteractiveUtils", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "13464637e13bc2a6577c2e456d561c5603ca54c7"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "3.19.1"

    [deps.LinearSolve.extensions]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveEnzymeExt = "EnzymeCore"
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveFastLapackInterfaceExt = "FastLapackInterface"
    LinearSolveForwardDiffExt = "ForwardDiff"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = ["Pardiso", "SparseArrays"]
    LinearSolveRecursiveFactorizationExt = "RecursiveFactorization"
    LinearSolveSparseArraysExt = "SparseArrays"
    LinearSolveSparspakExt = ["SparseArrays", "Sparspak"]

    [deps.LinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    FastLapackInterface = "29a986be-02c6-4525-aec4-84b980013641"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Sparspak = "e56a9233-b9d6-4f03-8d0f-1825330902ac"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

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

[[deps.LombScargle]]
deps = ["FFTW", "LinearAlgebra", "Measurements", "Random", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "d64a0ce7539181136a85fd8fe4f42626387f0f26"
uuid = "fc60dff9-86e7-5f2f-a8a0-edeadbb75bd9"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "5de60bc6cb3899cd318d80d627560fae2e2d99ae"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.0.1+1"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

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
deps = ["Calculus", "LinearAlgebra", "Printf"]
git-tree-sha1 = "030f041d5502dbfa41f26f542aaac32bcbe89a64"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.14.0"

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
version = "2023.12.12"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "fade91fe9bee7b142d332fc6ab3f0deea29f637b"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.9"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "491bdcdc943fcbc4c005900d7463c9f216aabf4c"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.4"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "25a6638571a902ecfb1ae2a18fc1575f86b1d4df"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.10.0"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ca7e18198a166a1f3eb92a3650d53d94ed8ca8a1"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.22"

[[deps.Neighborhood]]
deps = ["Distances", "NearestNeighbors", "Random", "Test"]
git-tree-sha1 = "fdea60ca30d724e76cc3b3d90d7f9d29d3d5cab5"
uuid = "645ca80c-8b79-4109-87ea-e1f58159d116"
version = "0.2.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "NonlinearSolveBase", "NonlinearSolveFirstOrder", "NonlinearSolveQuasiNewton", "NonlinearSolveSpectralMethods", "PrecompileTools", "Preferences", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SparseArrays", "SparseMatrixColorings", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "aeb6fb02e63b4d4f90337ed90ce54ceb4c0efe77"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "4.9.0"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = ["NLsolve", "LineSearches"]
    NonlinearSolvePETScExt = ["PETSc", "MPI"]
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
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"
    Sundials = "c3572dad-4567-51f8-b174-8c6c989267f4"

[[deps.NonlinearSolveBase]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "Compat", "ConcreteStructs", "DifferentiationInterface", "EnzymeCore", "FastClosures", "LinearAlgebra", "Markdown", "MaybeInplace", "Preferences", "Printf", "RecursiveArrayTools", "SciMLBase", "SciMLJacobianOperators", "SciMLOperators", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "404d71dd057759f4d590191a643113485c4a482a"
uuid = "be0214bd-f91f-a760-ac4e-3421ce2b2da0"
version = "1.12.0"
weakdeps = ["BandedMatrices", "DiffEqBase", "ForwardDiff", "LineSearch", "LinearSolve", "SparseArrays", "SparseMatrixColorings"]

    [deps.NonlinearSolveBase.extensions]
    NonlinearSolveBaseBandedMatricesExt = "BandedMatrices"
    NonlinearSolveBaseDiffEqBaseExt = "DiffEqBase"
    NonlinearSolveBaseForwardDiffExt = "ForwardDiff"
    NonlinearSolveBaseLineSearchExt = "LineSearch"
    NonlinearSolveBaseLinearSolveExt = "LinearSolve"
    NonlinearSolveBaseSparseArraysExt = "SparseArrays"
    NonlinearSolveBaseSparseMatrixColoringsExt = "SparseMatrixColorings"

[[deps.NonlinearSolveFirstOrder]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConcreteStructs", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLJacobianOperators", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "9c8cd0a986518ba317af263549b48e34ac8f776d"
uuid = "5959db7a-ea39-4486-b5fe-2dd0bf03d60d"
version = "1.5.0"

[[deps.NonlinearSolveQuasiNewton]]
deps = ["ArrayInterface", "CommonSolve", "ConcreteStructs", "DiffEqBase", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLOperators", "StaticArraysCore"]
git-tree-sha1 = "e3888bdbab6e0bfadbc3164ef4595e40e7b7e954"
uuid = "9a2c21bd-3a47-402d-9113-8faf9a0ee114"
version = "1.6.0"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveQuasiNewton.extensions]
    NonlinearSolveQuasiNewtonForwardDiffExt = "ForwardDiff"

[[deps.NonlinearSolveSpectralMethods]]
deps = ["CommonSolve", "ConcreteStructs", "DiffEqBase", "LineSearch", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "3398222199e4b9ca0b5840907fb509f28f1a2fdc"
uuid = "26075421-4e9a-44e1-8bd1-420ed7ad02b2"
version = "1.2.0"
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
git-tree-sha1 = "567515ca155d0020a45b05175449b499c63e7015"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.29+0"

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
version = "0.8.5+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ad31332567b189f508a3ea8957a2640b1147ab00"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["Compat", "EnumX", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "61942645c38dd2b5b78e2082c9b51ab315315d10"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.13.2"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqAdamsBashforthMoulton", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqExplicitRK", "OrdinaryDiffEqExponentialRK", "OrdinaryDiffEqExtrapolation", "OrdinaryDiffEqFIRK", "OrdinaryDiffEqFeagin", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqHighOrderRK", "OrdinaryDiffEqIMEXMultistep", "OrdinaryDiffEqLinear", "OrdinaryDiffEqLowOrderRK", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqNordsieck", "OrdinaryDiffEqPDIRK", "OrdinaryDiffEqPRK", "OrdinaryDiffEqQPRK", "OrdinaryDiffEqRKN", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqSSPRK", "OrdinaryDiffEqStabilizedIRK", "OrdinaryDiffEqStabilizedRK", "OrdinaryDiffEqSymplecticRK", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "Static", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "1c2b2df870944e0dc01454fd87479847c55fa26c"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.98.0"

[[deps.OrdinaryDiffEqAdamsBashforthMoulton]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqLowOrderRK", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "82f78099ecf4e0fa53545811318520d87e7fe0b8"
uuid = "89bda076-bce5-4f1c-845f-551c83cdda9a"
version = "1.2.0"

[[deps.OrdinaryDiffEqBDF]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqSDIRK", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "9124a686af119063bb4d3a8f87044a8f312fcad9"
uuid = "6ad6398a-0878-4a85-9266-38940aa047c8"
version = "1.6.0"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "FastBroadcast", "FastClosures", "FastPower", "FillArrays", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleUnPack", "Static", "StaticArrayInterface", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "1bd20b621e8dee5f2d170ae31631bf573ab77eec"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "1.26.2"

    [deps.OrdinaryDiffEqCore.extensions]
    OrdinaryDiffEqCoreEnzymeCoreExt = "EnzymeCore"
    OrdinaryDiffEqCoreMooncakeExt = "Mooncake"

    [deps.OrdinaryDiffEqCore.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"

[[deps.OrdinaryDiffEqDefault]]
deps = ["ADTypes", "DiffEqBase", "EnumX", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "PrecompileTools", "Preferences", "Reexport"]
git-tree-sha1 = "7e2f4ec76ebac709401064fd2cf73ad993d1e694"
uuid = "50262376-6c5a-4cf5-baba-aaf4f84d72d7"
version = "1.5.0"

[[deps.OrdinaryDiffEqDifferentiation]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqCore", "SciMLBase", "SciMLOperators", "SparseArrays", "SparseMatrixColorings", "StaticArrayInterface", "StaticArrays"]
git-tree-sha1 = "efecf0c4cc44e16251b0e718f08b0876b2a82b80"
uuid = "4302a76b-040a-498a-8c04-15b101fed76b"
version = "1.10.0"

[[deps.OrdinaryDiffEqExplicitRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "TruncatedStacktraces"]
git-tree-sha1 = "4dbce3f9e6974567082ce5176e21aab0224a69e9"
uuid = "9286f039-9fbf-40e8-bf65-aa933bdc4db0"
version = "1.1.0"

[[deps.OrdinaryDiffEqExponentialRK]]
deps = ["ADTypes", "DiffEqBase", "ExponentialUtilities", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqVerner", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "8d2ab84d7fabdfde995e5f567361f238069497f5"
uuid = "e0540318-69ee-4070-8777-9e2de6de23de"
version = "1.4.0"

[[deps.OrdinaryDiffEqExtrapolation]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastPower", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "80a636aac325c546b04e3bf20f0c80eaa0173dd4"
uuid = "becaefa8-8ca2-5cf9-886d-c06f3d2bd2c4"
version = "1.5.0"

[[deps.OrdinaryDiffEqFIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastGaussQuadrature", "FastPower", "LinearAlgebra", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "0da8ec3491821262a3d2828e6370e76b51a770a3"
uuid = "5960d6e9-dd7a-4743-88e7-cf307b64f125"
version = "1.12.0"

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
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Reexport"]
git-tree-sha1 = "095bab73a3ff185e9ef971fc42ecc93c7824e589"
uuid = "9f002381-b378-40b7-97a6-27a27c83f129"
version = "1.3.0"

[[deps.OrdinaryDiffEqLinear]]
deps = ["DiffEqBase", "ExponentialUtilities", "LinearAlgebra", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "940cef72ec8799d869ff1ba3dcf47cf7758e51cf"
uuid = "521117fe-8c41-49f8-b3b6-30780b3f0fb5"
version = "1.3.0"

[[deps.OrdinaryDiffEqLowOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "d4bb32e09d6b68ce2eb45fb81001eab46f60717a"
uuid = "1344f307-1e59-4825-a18e-ace9aa3fa4c6"
version = "1.2.0"

[[deps.OrdinaryDiffEqLowStorageRK]]
deps = ["Adapt", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "StaticArrays"]
git-tree-sha1 = "52ec7081e65291fa5c19749312df0818db2fa1bc"
uuid = "b0944070-b475-4768-8dec-fb6eb410534d"
version = "1.3.0"

[[deps.OrdinaryDiffEqNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "PreallocationTools", "RecursiveArrayTools", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "StaticArrays"]
git-tree-sha1 = "ffdb0f5207b0e30f8b1edf99b3b9546d9c48ccaf"
uuid = "127b3ac7-2247-4354-8eb6-78cf4e7c58e8"
version = "1.10.0"

[[deps.OrdinaryDiffEqNordsieck]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "ef44754f10e0dfb9bb55ded382afed44cd94ab57"
uuid = "c9986a66-5c92-4813-8696-a7ec84c806c8"
version = "1.1.0"

[[deps.OrdinaryDiffEqPDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "Reexport", "StaticArrays"]
git-tree-sha1 = "ab9897e4bc8e3cf8e15f1cf61dbdd15d6a2341d7"
uuid = "5dd0a6cf-3d4b-4314-aa06-06d4e299bc89"
version = "1.3.1"

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
deps = ["ADTypes", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "1ce0096d920e95773220e818f29bf4b37ea2bb78"
uuid = "43230ef6-c299-4910-a778-202eb28ce4ce"
version = "1.11.0"

[[deps.OrdinaryDiffEqSDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "b3a7e3a2f355d837c823b435630f035aef446b45"
uuid = "2d112036-d095-4a1e-ab9a-08536f3ecdbf"
version = "1.3.0"

[[deps.OrdinaryDiffEqSSPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "StaticArrays"]
git-tree-sha1 = "651756c030df7a1d49ad484288937f8c398e8a08"
uuid = "669c94d9-1f4b-4b64-b377-1aa079aa2388"
version = "1.3.0"

[[deps.OrdinaryDiffEqStabilizedIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "111c23b68ad644b47e38242af920d5805c7bedb1"
uuid = "e3e12d00-db14-5390-b879-ac3dd2ef6296"
version = "1.3.0"

[[deps.OrdinaryDiffEqStabilizedRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "1b0d894c880e25f7d0b022d7257638cf8ce5b311"
uuid = "358294b1-0aab-51c3-aafe-ad5ab194a2ad"
version = "1.1.0"

[[deps.OrdinaryDiffEqSymplecticRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "a13d59a2d6cfb6a3332a7782638ca6e1cb6ca688"
uuid = "fa646aed-7ef9-47eb-84c4-9443fc8cbfa8"
version = "1.3.0"

[[deps.OrdinaryDiffEqTsit5]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "96552f7d4619fabab4038a29ed37dd55e9eb513a"
uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a"
version = "1.1.0"

[[deps.OrdinaryDiffEqVerner]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "08f2d3be30874b6e2e937a06b501fb9811f7d8bd"
uuid = "79d7bb75-1356-48c1-b8c0-6832512096c2"
version = "1.2.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f07c06228a1c670ae4c87d1276b92c7c597fdda0"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.35"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

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
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ec9e63bd098c50e4ad28e7cb95ca7a4860603298"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.68"

[[deps.PoissonRandom]]
deps = ["LogExpFunctions", "Random"]
git-tree-sha1 = "bb178012780b34046c6d1600a315d8dbee89d83d"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.5"

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

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "972089912ba299fba87671b025cd0da74f5f54f7"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.0"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
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
git-tree-sha1 = "2cc315bb7f6e4d59081bad744cdb911d6374fc7f"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.29"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"
    PreallocationToolsSparseConnectivityTracerExt = "SparseConnectivityTracer"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"

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
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "13c5103482a8ed1536a54c08d0e742ae3dca2d42"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.4"

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

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "efc718978d97745c58e69c5115a35c51a080e45e"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.34.1"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsKernelAbstractionsExt = "KernelAbstractions"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStructArraysExt = "StructArrays"
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
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
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
git-tree-sha1 = "668e411c0616a70860249b4c96e5d35296631a1d"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.8"

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
git-tree-sha1 = "86a8a8b783481e1ea6b9c91dd949cb32191f8ab4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.15"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "e6a28a9a2dd9bc3ed46391fa0e6c35839bde4028"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.103.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMLStyleExt = "MLStyle"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    MLStyle = "d8e11817-5142-5d16-987a-aa16d5891078"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLJacobianOperators]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DifferentiationInterface", "FastClosures", "LinearAlgebra", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "7da1216346ad79499d08d7e2a3dbf297dc80c829"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.6"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "3249fe77f322fe539e935ecb388c8290cd38a3fc"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.3.1"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "566c4ed301ccb2a44cbd5a27da5f885e0ed1d5df"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

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
git-tree-sha1 = "7aaa5fe4617271b64fce0466d187f2a72edbd81a"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "2.5.0"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolveDiffEqBaseExt = "DiffEqBase"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveTrackerExt = "Tracker"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffEqBase = "2b5f629d-d688-5b77-993f-72d75c75574e"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

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
git-tree-sha1 = "182990067a09adf950274f97f38f68c76f81d2d0"
uuid = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
version = "0.6.21"

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

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "Random", "SparseArrays"]
git-tree-sha1 = "9de43e0b9b976f1019bf7a879a686c4514520078"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.21"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsCUDAExt = "CUDA"
    SparseMatrixColoringsCliqueTreesExt = "CliqueTrees"
    SparseMatrixColoringsColorsExt = "Colors"

    [deps.SparseMatrixColorings.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "41852b8679f78c8d8961eeadc8f62cef861a52e3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StateSpaceSets]]
deps = ["Distances", "LinearAlgebra", "Neighborhood", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "b39e3bb75852c62514e474131e48607300349cc3"
uuid = "40b095a5-5852-4c12-98c7-d43bf788e795"
version = "2.4.0"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools"]
git-tree-sha1 = "f737d444cb0ad07e61b3c1bef8eb91203c321eff"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.2.0"

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
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"
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
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "b81c5035922cc89c2d9523afc6c54be512411466"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.5"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "8e45cecc66f3b42633b8ce14d431e8e57a3e242e"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.SteadyStateDiffEq]]
deps = ["ConcreteStructs", "DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NonlinearSolveBase", "Reexport", "SciMLBase"]
git-tree-sha1 = "66a028f9a2bb44d0f6de0814a2b9840af548143a"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "2.5.0"

[[deps.StochasticDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FastPower", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SparseArrays", "StaticArrays", "UnPack"]
git-tree-sha1 = "2992af2739fdcf5862b12dcf53a5f6e3e4acd358"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.80.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "f35f6ab602df8413a50c4a25ca14de821e8605fb"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.7"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "725421ae8e530ec29bcbdddbe91ff8053421d023"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.1"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "7c7a7ee705724b3c80d5451ac49779db36c6f758"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.28.0"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "SuiteSparse_jll", "libblastrampoline_jll"]
git-tree-sha1 = "91db7ed92c66f81435fe880947171f1212936b14"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.3+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "PrettyTables", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "658f6d01bfe68d6bf47915bf5d868228138c7d71"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.41"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "fabf4650afe966a2ba646cabd924c3fd43577fc3"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "ExproniconLite", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "TimerOutputs", "Unityper", "WeakValueDicts"]
git-tree-sha1 = "fa63e8f55e99aee528951ba26544403b09645979"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "3.29.0"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "LaTeXStrings", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "OffsetArrays", "PrecompileTools", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "3d7b491b60bf3d24a5c2a74821db4520f0307c72"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.44.0"

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
deps = ["LaTeXStrings", "Poppler_jll", "Requires", "tectonic_jll"]
git-tree-sha1 = "79e2d29b216ef24a0f4f905532b900dcf529aa06"
uuid = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
version = "3.5.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.Trapz]]
git-tree-sha1 = "79eb0ed763084a3e7de81fe1838379ac6a23b6a0"
uuid = "592b5752-818d-11e9-1e9a-2b8ca4a44cd1"
version = "2.0.3"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

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

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d2282232f8a4d71f79e85dc4dd45e5b12a6297fb"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.23.1"
weakdeps = ["ConstructionBase", "ForwardDiff", "InverseFunctions", "Printf"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    PrintfExt = "Printf"

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

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.WeakValueDicts]]
git-tree-sha1 = "98528c2610a5479f091d470967a25becfd83edd0"
uuid = "897b6980-f191-5a31-bcb0-bf3c4585e0c1"
version = "0.1.0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "b8b243e47228b4a3877f1dd6aee0c5d56db7fcf4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.6+1"

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
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "07b6a107d926093898e82b3b1db657ebe33134ec"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.50+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

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
# ╟─800dc762-ce0b-463f-858f-6e8eabbc26b0
# ╟─03943fac-642a-4063-b043-a6b05b417c7d
# ╟─417aad38-6928-4b13-9286-3f11efcffb99
# ╟─b842e98e-34e2-40f2-84b6-c180815c2df3
# ╟─4fc9bda3-1f0f-41c2-8005-5b8613271d4a
# ╟─3c0c1fea-1ae6-4c9d-9db1-764937cd5fbc
# ╟─99d0fdc9-e368-462c-a357-86f07624a52e
# ╟─6d6cbebb-bb3b-437d-9835-0a79f36857f2
# ╟─71b94af7-76ee-4228-987c-2f22e0951552
# ╟─22c37732-1cf7-4d80-a250-9fd5e4a2f88c
# ╟─8da87629-09ff-4677-b5e8-2a3ff5777a8e
# ╟─c5e21675-f120-4555-84be-d99b35592f2d
# ╟─8dcc45e1-192a-4007-9d6e-6e4149b563d7
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
# ╠═cea20293-a562-4b30-94c1-aabad784fbfc
# ╠═5b08171f-5ee3-431c-8fc8-f666a4f6ee5d
# ╠═6c533a69-6d71-4fea-8aad-224c0bf7d53c
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
# ╠═2a39d6f8-da49-4976-9aa7-889391e55a5d
# ╟─ada0a462-e4cf-438d-b2a7-35d7280e955a
# ╠═6a59ed2e-f040-4e10-bc55-91d2c1dcc97e
# ╟─93c15a25-51da-4019-b09c-83a95aa5999d
# ╟─43ee281f-1a16-445d-894d-23e0319b1fd0
# ╠═c8efcb3f-cc6f-4c45-a0f3-55cdb73d5195
# ╟─4d09ed45-423c-4bd6-802d-59389a966d2e
# ╟─35083c71-0a7c-4b97-94ef-4d06ecbd2ee8
# ╟─76a10b7c-1135-49f1-a298-36c59bb94b37
# ╠═72794eb3-7f09-4445-906a-af9756dceebd
# ╟─221942a6-5e3b-46f2-ac68-ac0a55d17f04
# ╟─b96de0fd-a65b-463a-9c91-2504b8427dba
# ╠═1f190b07-164b-4ba9-86f8-edcf6e7ab3a9
# ╟─9a756404-f35b-40b9-bc88-4b7e3a7df0b6
# ╠═805ef471-1bd2-4ce2-bf50-eeb9c6b10e4d
# ╟─6eef3827-1af6-431d-a163-67bb595c337d
# ╠═e14648ff-69f6-47c2-924c-c5ac69e200ed
# ╟─aaba63dd-8fc3-405d-ad19-8e8368e70019
# ╠═03754250-3831-41db-ac11-fe6f5a08f0c1
# ╠═2d27f581-2c53-409e-afde-812e70bba8bf
# ╟─6ddf393c-9f37-43b2-b401-9ab0e7a9e88c
# ╠═5ce798b5-4366-4cc1-8d2d-881a827d6dc9
# ╠═400b6787-f2b5-40ef-be11-cc5d60d19b9f
# ╟─7df63fdf-ef41-44f2-8a55-e5c2c849029c
# ╟─0895d464-a029-410d-9e7d-89cfac2d1615
# ╟─137c43b7-a541-4ca0-a6d0-0ad3579b345a
# ╠═c195ae82-0e44-4181-b863-bdd365059f4b
# ╠═14e326a5-d2bf-4873-9cf7-b57be9416da2
# ╟─0bdf9dbf-479c-46f6-bd86-50576095cba0
# ╟─6fe43e3a-2e8f-4708-a3ec-6f5a8088060e
# ╠═f863d68f-590e-4b96-8433-dc6b5177539f
# ╠═9f18a8cb-93a6-4a87-b0c4-db27ab7b80bd
# ╟─f5a983bf-ef3a-4d5d-928d-da1097b91ee8
# ╠═bd8743d6-8f21-413d-835a-e543926baa09
# ╟─9ab0a10b-8165-401a-a2b6-c9726526a906
# ╠═1b8af600-56eb-4508-bc52-aa4e138b4c7e
# ╠═9958610a-3384-4410-beae-8f4ce657846d
# ╠═80005099-7154-4306-9172-c9a168336e14
# ╠═c291700e-3a84-49a7-85d6-592cfb3b1a11
# ╟─1ec99792-905d-4e1b-a413-ef58143d3c68
# ╠═4e1140be-3e61-47bd-9f9f-6b5dfbff6ae2
# ╟─2afd459c-90e9-4105-9121-27e21bb89eeb
# ╠═4173e91c-c01a-4fa2-ae77-0408bf7d9a1b
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
# ╟─24920a1b-355c-43e1-aff2-028e6f25584c
# ╟─0813f6a3-aadc-491c-ad1e-ec54dcbd0d56
# ╟─de4ba100-5f58-4657-aea0-a3f31eacda65
# ╠═77ca6eab-c5bf-479c-8982-1483933bbb6e
# ╟─d4a1ec85-e6f5-48ed-9724-202c4dae9ae9
# ╠═7f8ebe64-ea86-4eb8-9939-6682180b9dd2
# ╟─aec6e4fc-e496-4add-b982-ab60f9f900a0
# ╠═8940cb2b-2c8c-407e-bc6b-426843cf6125
# ╠═28955b95-df19-403a-bf79-b68e9be8e1dd
# ╟─74c82c98-f233-4e72-8f74-17bdfeddf884
# ╟─ab1b7695-3cd6-4e67-9cfd-fbaf2f4d1d15
# ╟─6475b25b-8711-44cd-bc3b-e3d42681bb93
# ╟─ea4e58e9-d041-4a6e-b0d8-83e3aef7648b
# ╠═d7cc8a66-220a-4e9b-b42e-a0e085ed3a0f
# ╟─477c8f59-97d1-405d-a1c8-64b7e0b9119f
# ╟─08e05c65-06a2-4560-92eb-b014dc7c3d70
# ╟─2235689d-9c83-4907-aa17-c2624fbeb68d
# ╟─df8a9449-851c-4546-97a7-7fd4a270a867
# ╠═245ad10f-4c70-4585-b2e0-dea4cf5d19b8
# ╠═e8c51ae9-ded4-450a-80f9-c484c2b991d1
# ╠═c96ea0ad-47c3-469e-a2f7-0743e681a6d1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
