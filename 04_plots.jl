### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ f79c680a-7c97-4100-97e8-19190707c55b
let
	using CairoMakie, ColorSchemes, Measurements, Libdl, Printf, PlutoLinks, Statistics, StatsBase

	using ChaosTools, DataFrames, CSV, DifferentialEquations, Interpolations, LinearAlgebra, Measurements, NaNMath, PlutoUI, QuadGK, SpecialFunctions, Symbolics, TikzPictures, Unitful, UnitfulAstro
end

# ╔═╡ 37653018-aaf5-42d1-a938-e05a44f18918
# ╠═╡ skip_as_script = true
#=╠═╡
TableOfContents(title="Plots", depth=4)
  ╠═╡ =#

# ╔═╡ 6e5823eb-3614-476b-80e7-66fc27ed0392
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	MARKERS = [
		:circle,
		:rect,
		:dtriangle,
		:utriangle,
		:cross,
		:diamond,
		:ltriangle,
		:rtriangle,
		:pentagon,
		:xcross,
		:hexagon,
	]

	LINE_STYLES=[:solid, :dash, :dot, :dashdot, :dashdotdot]

	CYCLE = Cycle([:color, :linestyle, :marker])

	DEFAULT_THEME = Theme(
	    size=(880, 880),
	    fontsize=35,
	    figure_padding=(1, 15, 5, 15),
	    palette=(color=Makie.wong_colors(), marker=MARKERS, linestyle=LINE_STYLES),
	    CairoMakie=(px_per_unit=2.3622, pt_per_unit=0.283466),
	    Axis=(
	        xlabelpadding=15,
	        xticklabelpad=10,
	        xticksize=7,
	        xgridvisible=false,
	        spinewidth=3,
	        xminorticksvisible=true,
	        xminorticks=IntervalsBetween(5),
	        ylabelpadding=15,
	        yticklabelpad=10,
	        yticksize=7,
	        ygridvisible=false,
	        yminorticksvisible=true,
	        yminorticks=IntervalsBetween(5),
	        aspect=AxisAspect(1),
	    ),
	    Legend=(
	        tellheight=false,
	        tellwidth=false,
	        framevisible=false,
	        colgap=20,
	        halign=:right,
	        valign=:bottom,
	        nbanks=3,
	        titlegap=-5,
	        labelsize=30,
	        linewidth=5,
	        markersize=28,
	        patchsize=(50, 50),
	        linepoints=[Point2f(0.0, 0.5), Point2f(0.9, 0.5)],
	        polypoints=[
	            Point2f(0.15, 0.15),
	            Point2f(0.85, 0.15),
	            Point2f(0.85, 0.85),
	            Point2f(0.15, 0.85),
	        ],
	    ),
	    Lines=(linewidth=5, cycle=CYCLE),
	    VLines=(linewidth=3, cycle=CYCLE),
	    HLines=(linewidth=3, cycle=CYCLE),
	    Heatmap=(colormap=:CMRmap, nan_color=ColorSchemes.CMRmap[1]),
	    Scatter=(markersize=22, cycle=CYCLE),
	    Colorbar=(
	        colormap=:CMRmap,
	        size=25,
	        ticklabelpad=10,
	        minorticksvisible=true,
	        ticksize=7,
	        labelpadding=2,
	    ),
	)

	THEME = merge(DEFAULT_THEME, theme_latexfonts())
end;
  ╠═╡ =#

# ╔═╡ 4b65fd9c-bd69-4e97-a28c-357c6a9c7bda
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Load the model"
  ╠═╡ =#

# ╔═╡ a0ca487a-8a51-48cb-acc0-7b318c793d09
begin
	const TESTING = @ingredients("./03_testing.jl")
	const CODEGEN = TESTING.CODEGEN
	const MODEL   = CODEGEN.MODEL
end;

# ╔═╡ cfacae0b-ee2e-4ba5-b31e-0fc122a3676b
# ╠═╡ skip_as_script = true
#=╠═╡
md"# Plots"
  ╠═╡ =#

# ╔═╡ d5d5c8f9-4eb2-4332-8de5-6e20da7e528a
begin
	# If the plots will be made with the Julia implementation
	const SOVE_WITH_JULIA = true

	if SOVE_WITH_JULIA

		##################################################
		# Solver functions using the Julia implementation 
		##################################################
		
		function solve_ode(
			ic::Vector{Float64},
			base_params::Vector{Float64},
			time_list::Vector{Float64},
		)::Vector{Vector{Float64}}

			return MODEL.integrate_model(
				ic,
				base_params,
				(time_list[1], time_list[end]);
				times=time_list,
			)

		end

		function solve_ode(
			ic::Vector{Float64},
			base_params::Vector{Float64},
			it::Float64,
		)::Vector{Float64}

			return MODEL.integrate_model(ic, base_params, (0.0, it))[1]

		end

	else

		##############################################
		# Solver functions using the C implementation 
		##############################################
		
		c_solver = TESTING.c_solver_generator();

		function solve_ode(
			ic::Vector{Float64},
			base_params::Vector{Float64},
			time_list::Vector{Float64},
		)::Vector{Vector{Float64}}

			return [c_solver(ic, base_params, it) for it in time_list]

		end

		function solve_ode(
			ic::Vector{Float64},
			 base_params::Vector{Float64},
			it::Float64,
		)::Vector{Float64}

			return c_solver(ic, base_params, it)

		end

	end
end;

# ╔═╡ b0c851e1-2cf7-48b7-b045-f5960147fb04
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Physical processes"
  ╠═╡ =#

# ╔═╡ 78385b04-ccee-48ad-8bb3-c6b1af1bfd1e
# ╠═╡ skip_as_script = true
#=╠═╡
diagram = TikzPictures.TikzPicture(
	L"""
		\node[box] (stars) {Stars};
		\node[box, text width=2em, above=2.5cm of stars] (atom) {HI};
		\node[box, text width=2em, right=2cm of atom] (molecule) {\ch{H2}};
		\node[box, text width=2em, left=2cm of atom] (ion) {HII};
		\node[box, text width=2em, below=0.8cm of stars] (metals) {Z};
		\node[box, text width=2em, below=0.8cm of metals] (dust) {Dust};
		\draw[line, ->]
		(ion) edge [bend left, "recombination"] (atom)
		(atom) edge [bend left, "condensation"] (molecule)
		(molecule) edge [bend left,"dissociation"] (atom)
		(atom) edge [bend left,"ionization"] (ion)
		(stars) edge [bend left, "supernova"] (ion)
		(molecule) edge [bend left, "star formation"] (stars)
		(stars) edge ["supernova"] (metals)
		(metals) edge [bend left, "dust growth"] (dust)
		(dust) edge [bend left, "dust destruction"] (metals);
	""",
	width="75em",
	preamble = """
		\\usepackage{chemformula}
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

# ╔═╡ 6162a115-fdb7-43eb-becb-0f25192001e0
# ╠═╡ skip_as_script = true
#=╠═╡
let
    mkpath("./generated_files/plots")
    TikzPictures.save(PDF("./generated_files/plots/physical_processes"), diagram)
end;
  ╠═╡ =#

# ╔═╡ 6775e416-0878-4d82-9b3a-067eeebbe3ad
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Parameters"
  ╠═╡ =#

# ╔═╡ 25ae9eb0-86ce-4e73-8b4b-128c33f15e4e
# ╠═╡ skip_as_script = true
#=╠═╡
let
	resolution = 100

	##################
	# Fiducial values
	##################

	# Neutral fraction
	fn = 0.8
	# Gas fraction (fi + fa + fm)
	fg = 1.0

	# Total cell density [mp * cm^(-3)]
	log_ρc_range = range(-1, 3, resolution)
	ρc_range = exp10.(log_ρc_range)

    # Ionized gas fraction [dimensionless]
	logfi_range = range(-2, 0, resolution)
	fi_range = exp10.(logfi_range)

	# Metallicity [dimensionless]
	solarZt_range = range(0.0, 1.0, resolution)
	Zt_range = solarZt_range .* MODEL.Zsun

	with_theme(THEME) do

		f = Figure(size=(880, 1200), figure_padding=(5, 10, 10, 1))

		#########
		# τ_star
		#########

		# Time scale
		τstar = [MODEL.τ_star(ρc) for ρc in ρc_range]
		ext_τstar = extrema(τstar)
		log_τstar = log10.(τstar)

		# Color variable
		ρ_τs(τs) = 1.0 / (MODEL.c_star * τs)^2
		c_range = [ρ_τs(τs) for τs in τstar, _ in 1:1]
		log_c_range = replace(x -> -1 < x < 3 ? x : NaN, log10.(c_range))

		ax_01 = CairoMakie.Axis(
			f[1, 1];
			limits=(-4, 5, nothing, nothing),
			title=L"\tau_\mathrm{star}",
			backgroundcolor=:gray85,
			aspect=nothing,
		)

		hm = heatmap!(
			ax_01,
			log_τstar,
			[1,],
			log_c_range;
			nan_color=:gray85,
			colorrange=(-1.0, 3.0),
		)

		hidedecorations!(ax_01)

		println("τ_star range: $(round.(ext_τstar; sigdigits=4)) Myr\n")

		###########
		# Colorbar
		###########

        Colorbar(
			f[1:4, 2],
			hm;
			label=L"\log_{10} \, \rho_\mathrm{cell} \,\, [\mathrm{%$(MODEL.l_u)^{-3}}]",
		)

		#########
		# τ_cond
		#########

		# Time scale
		τcond = [MODEL.τ_cond(ρc, Zt, fg - Zt) for ρc in ρc_range, Zt in Zt_range]
		ext_τcond = extrema(τcond)
		log_τcond = range(log10.(ext_τcond)..., resolution)
		liniar_τcond = exp10.(log_τcond)

		# Color variable
		ρ_τc(τc, Zt, fg) = 1.0 / (MODEL.c_cond * τc * fg * (Zt + MODEL.Zeff))
		c_range = [ρ_τc(τc, Zt, fg - Zt) for τc in liniar_τcond, Zt in Zt_range]
		log_c_range = replace(x -> -1 < x < 3 ? x : NaN, log10.(c_range))

		ax_02 = CairoMakie.Axis(
			f[2, 1];
			ylabel=L"Z \, / \, Z_\odot",
			limits=(-4, 5, 0.0, 1.0),
			title=L"\tau_\mathrm{cond}",
			backgroundcolor=:gray85,
			aspect=nothing,
			yticks=[0.0, 0.5, 1.0],
		)

		hidexdecorations!(ax_02)

		heatmap!(
			ax_02,
			log_τcond,
			solarZt_range,
			log_c_range;
			nan_color=:gray85,
			colorrange=(-1.0, 3.0),
		)

		println("τ_cond range: $(round.(ext_τcond; sigdigits=4)) Myr\n")

		########
		# τ_rec
		########

		# Time scale
		τrec = [MODEL.τ_rec(fi, ρc) for ρc in ρc_range, fi in fi_range]
		ext_τrec = extrema(τrec)
		log_τrec = range(log10.(ext_τrec)..., resolution)
		liniar_τrec = exp10.(log_τrec)

		# Color variable
		ρ_τr(fi, τr) = 1.0 / (MODEL.c_rec * fi * τr)
		c_range = [ρ_τr(fi, τr) for τr in liniar_τrec, fi in fi_range]
		log_c_range = replace(x -> -1 < x < 3 ? x : NaN, log10.(c_range))

		ax_03 = CairoMakie.Axis(
			f[3, 1];
			xlabel=L"\log_{10} \, \tau \,\, [\mathrm{%$(MODEL.t_u)}]",
			ylabel=L"f_i",
			limits=(-4, 5, 0.0, 1.0),
			title=L"\tau_\mathrm{rec}",
			backgroundcolor=:gray85,
			aspect=nothing,
			yticks=[0.0, 0.5, 1.0],
		)

		hidexdecorations!(ax_03)

		heatmap!(
			ax_03,
			log_τrec,
			fi_range,
			log_c_range;
			nan_color=:gray85,
			colorrange=(-1.0, 3.0),
		)

		println("τ_rec range: $(round.(ext_τrec; sigdigits=4)) Myr\n")

		#######
		# τ_dg
		#######

		# Ignore the case Z = 0
		p_Zt_range = Zt_range[2:end]
		p_solarZt_range = solarZt_range[2:end]

		# Time scale
		τdg = [MODEL.τ_dg(ρc, fn, Zt) for ρc in ρc_range, Zt in p_Zt_range]
		ext_τdg = extrema(τdg)
		log_τdg = range(log10.(ext_τdg)..., resolution)
		liniar_τdg = exp10.(log_τdg)

		# Color variable
		ρ_τdg(fn, τdg, Zt) = 1 / (MODEL.c_dg * τdg * fn * Zt)
		c_range = [ρ_τdg(fn, τdg, Zt) for τdg in liniar_τdg, Zt in p_Zt_range]
		log_c_range = replace(x -> -1 < x < 3 ? x : NaN, log10.(c_range))

		ax_04 = CairoMakie.Axis(
			f[4, 1];
			xlabel=L"\log_{10} \, \tau \,\, [\mathrm{%$(MODEL.t_u)}]",
			ylabel=L"Z \, / \, Z_\odot",
			limits=(-4, 5, 0.0, 1.0),
			title=L"\tau_\mathrm{dg}",
			backgroundcolor=:gray85,
			aspect=nothing,
			yticks=[0.0, 0.5, 1.0],
		)

		heatmap!(
			ax_04,
			log_τdg,
			p_solarZt_range,
			log_c_range;
			nan_color=:gray85,
			colorrange=(-1.0, 3.0),
		)

		println("τ_dg range: $(round.(ext_τdg; sigdigits=4)) Myr")

		#############################################################################
		# Save the plot
		#############################################################################
		
		rowsize!(f.layout, 1, Relative(1/20))

		mkpath("./generated_files/plots/model/")
		Makie.save("./generated_files/plots/model/timescale_comparison.png", f)

		f

	end
end
  ╠═╡ =#

# ╔═╡ dfb8b71a-b3a3-4ea0-9ea3-213c206117ca
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Photodissociation and mass recycling"
  ╠═╡ =#

# ╔═╡ 4ad4589f-9834-44a5-9c4b-994b7d91e442
# ╠═╡ skip_as_script = true
#=╠═╡
let
	# Stellar ages, as log10(age [Myr])
	log_ages = range(-1, 2, 200)

	# Stellar ages in [Myr]
	ages = exp10.(log_ages)

	# Stellar metalicities [Z⊙]
	metalicities = range(0, 2.0, 200)

	# Photodissociation efficiency
	η_diss(t, Z) = MODEL.compute_η(t, Z * MODEL.Zsun)[1]
	
	# Photoionization efficiency
	η_ion(t, Z) = MODEL.compute_η(t, Z * MODEL.Zsun)[2]
	
	# Mass recycling fraction
	R(Z) = MODEL.compute_recycled_fractions(Z * MODEL.Zsun)[1]
	
	with_theme(THEME) do

		f = Figure(size=(800, 2200))

		##################################
		# Photoionization efficiency plot
		##################################
		
		ax = CairoMakie.Axis(
			f[1, 1],
			xlabel=L"\mathrm{stellar \,\, age \,\, [Myr]}",
			ylabel=L"\eta_\mathrm{ion}",
			xscale=log10,
			aspect=AxisAspect(1),
		)

		for Z in [0.0, 0.5, 1.0]
			label = L"Z \, / \, Z_\odot = %$(Z)"
			lines!(ax, ages, x -> η_ion(x, Z); label)
		end

		axislegend(ax; position=:lt, nbanks=1, labelsize=28)

		####################################
		# Photodissociation efficiency plot
		####################################
		
		ax = CairoMakie.Axis(
			f[2, 1],
			xlabel=L"\mathrm{stellar \,\, age \,\, [Myr]}",
			ylabel=L"\eta_\mathrm{diss}",
			xscale=log10,
			aspect=AxisAspect(1),
		)

		for Z in [0.0, 0.5, 1.0]
			label = L"Z \, / \, Z_\odot = %$(Z)"
			lines!(ax, ages, x -> η_diss(x, Z); label)
		end

		axislegend(ax; position=:lt, nbanks=1, labelsize=28)

		###############################
		# Mass recycling fraction plot
		###############################

		ax = CairoMakie.Axis(
			f[3, 1],
			xlabel=L"Z \, / \, Z_\odot",
			ylabel=L"R",
			aspect=AxisAspect(1),
		)

		lines!(ax, metalicities, R)

		f

	end
end
  ╠═╡ =#

# ╔═╡ b2197f1a-ff9a-429c-8486-dfa27abc880d
# ╠═╡ skip_as_script = true
#=╠═╡
let
	# Stellar ages, as log10(age [Myr])
	log_ages = range(-1, 2, 200)

	# Stellar ages in [Myr]
	ages = exp10.(log_ages)

	# Photodissociation efficiency
	η_diss(t, Z) = MODEL.compute_η(t, Z * MODEL.Zsun)[1]
	
	# Photoionization efficiency
	η_ion(t, Z) = MODEL.compute_η(t, Z * MODEL.Zsun)[2]

	with_theme(THEME) do

		f = Figure(size=(880, 880))

		ax = CairoMakie.Axis(
			f[1,1],
			xlabel=L"\mathrm{stellar \,\, age \,\, [Myr]}",
			ylabel=L"\eta",
			xscale=log10,
			aspect=AxisAspect(1),
		)

		#############################
		# Photoionization efficiency
		#############################

		labels = [
			L"\eta_\mathrm{ion} \,\, (Z = 0.0)",
			L"\eta_\mathrm{ion} \,\, (Z = Z_\odot)",
		]

		for (Z, label) in zip([0.0, 1.0], labels)
			lines!(ax, ages, x -> η_ion(x, Z); label)
		end

		###############################
		# Photodissociation efficiency 
		###############################

		labels = [
			L"\eta_\mathrm{diss} \,\, (Z = 0.0)",
			L"\eta_\mathrm{diss} \,\, (Z = Z_\odot)",
		]

		for (Z, label) in zip([0.0, 1.0], labels)
			lines!(ax, ages, x -> η_diss(x, Z); label)
		end

		axislegend(ax; position=:lt, nbanks=1)

		f

	end
end
  ╠═╡ =#

# ╔═╡ ad34b750-5e81-4cb8-a444-52220c116790
# ╠═╡ skip_as_script = true
#=╠═╡
md"## Integration"
  ╠═╡ =#

# ╔═╡ b31c615a-36b3-46c0-bf8e-f8b5a77296f4
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Fractions vs. time"
  ╠═╡ =#

# ╔═╡ 6f341d1b-8b0f-4bc5-b18c-d1e9dc9281d5
# ╠═╡ skip_as_script = true
#=╠═╡
let
	#################################################################################
	# Metallicity x density grid
	#################################################################################

	frac_labels = [L"f_i", L"f_a", L"f_m", L"f_s", L"f_Z", L"f_d"]

	xaxis_visible = [false, false, false, false, false, false, true, true, true]
	yaxis_visible = [true, false, false, true, false, false, true, false, false]

	# Initial conditions
	# [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
	ic(Z, fi) = [fi - Z, 1.0 - fi, 0.0, 0.0, 0.9 * Z, 0.1 * Z]

	# Total cell density values [mp * cm^(-3)]
	ρ_str  = [L"0.19", L"1", L"10"]
	ρ_list = [0.19, 1.0, 10.0]

	# Metallicity values [dimensionless]
	Z_str  = [L"0.0", L"0.1", L"1.0"]
	Z_list = [0.0, 0.1, 1.0] * MODEL.Zsun

	# Parameters for the ODEs
	# [ρ_cell, Z, a, h]
	print_params = [[ρ, Z] for ρ in ρ_str, Z in Z_str]
	base_params  = [[ρ, Z, 1.0, 3.0e20] for ρ in ρ_list, Z in Z_list]

	# Integration time [Myr]
	time_list  = exp10.(range(-3.0, 1.0, 1000))
	time_range = (time_list[1], time_list[end])

	# Numerical integration
	fractions_05 = [
		solve_ode(ic(base_param[2], 0.5), base_param, time_list)
		for base_param in base_params
	]

	fractions_001 = [
		solve_ode(ic(base_param[2], 0.01), base_param, time_list)
		for base_param in base_params
	]

	iterator = zip(
		print_params,
		fractions_05,
		fractions_001,
		xaxis_visible,
		yaxis_visible,
	)

	scale = 1700 / length(ρ_list)

	with_theme(THEME) do

		f = Figure(
			size=(1700, scale * length(Z_list)),
			figure_padding=(1, 15, 5, 15),
		)

		for (
			idx,
			((ρ_str, Z_str), fraction_05, fraction_001, xaxis_v, yaxis_v),
		) in enumerate(iterator)

			row = ceil(Int, idx / length(ρ_list))
			col = mod1(idx, length(Z_list))

			ax = CairoMakie.Axis(
				f[row, col];
				xlabel=L"\log_{10} \, t \,\, [\mathrm{%$(MODEL.t_u)}]",
		        ylabel=L"f",
				title=L"$Z / Z_\odot =$%$Z_str$\quad \rho_\mathrm{cell} =$%$ρ_str$\, \mathrm{%$(MODEL.l_u)^{-3}}$",
				xminorticksvisible=xaxis_v,
		        xticksvisible=xaxis_v,
		        xlabelvisible=xaxis_v,
		        xticklabelsvisible=xaxis_v,
				yminorticksvisible=yaxis_v,
		        yticksvisible=yaxis_v,
		        ylabelvisible=yaxis_v,
		        yticklabelsvisible=yaxis_v,
				limits=((-3.1, 1.1), (-0.05, 1.02)),
				xticks=[-3, -2, -1, 0, 1],
			)

			for (label, idx, color) in zip(frac_labels, 1:6, Makie.wong_colors())

				lines!(
					ax,
					log10.(time_list),
					getindex.(fraction_001, idx);
					linestyle=:solid,
					label,
					color,
				)

				lines!(
					ax,
					log10.(time_list),
					getindex.(fraction_05, idx);
					linestyle=:dash,
					color,
				)

			end

			if row == 3 && col == 1
				axislegend(ax; position=(0.05, 0.9), nbanks=2, labelsize=32)
			end

		end

		mkpath("./generated_files/plots/model/")
		Makie.save("./generated_files/plots/model/fractions-vs-time.png", f)

		f

	end
end
  ╠═╡ =#

# ╔═╡ 1d870ff9-fade-43c8-af55-6137ce051b8e
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Fractions vs. density"
  ╠═╡ =#

# ╔═╡ c76cf08f-aaef-46df-ae18-9cd018141b70
# ╠═╡ skip_as_script = true
#=╠═╡
let
	#################################################################################
	# ICs x fractions grid
	#################################################################################

	phase_labels = [L"f_i", L"f_a", L"f_m", L"f_s", L"f_Z", L"f_d"]
	phases       = [
		"ionized",
		"atomic",
		"molecular",
		"stellar",
		"metals",
		"dust",
	]

	ylimits = [
		(-0.02, 0.62),
		(-0.04, 1.02),
		(-0.04, 1.02),
		(-0.004, 0.154),
		(-0.0004, 0.0154),
		(-0.0002, 0.0102),
	]

	# Axis visibility
	xaxis_visible = [false, false, false, false, false, true]
	yaxis_visible = [true, false, false]
	title_visible = [true, false, false, false, false, false]

	# Initial conditions
	# [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
	fis = [0.01, 0.25, 0.5]
	ic(Z, fi) = [fi - Z, 1.0 - fi, 0.0, 0.0, 0.9 * Z, 0.1 * Z]

	# Total cell density range [mp * cm^(-3)]
	ρ_exp_list = range(-1, 2, 100)
	ρ_list     = exp10.(ρ_exp_list)

	# Metallicities [dimensionless]
	Z_labels = [
		L"Z / Z_\odot = 0.0",
		L"Z / Z_\odot = 10^{-3}",
		L"Z / Z_\odot = 10^{-2}",
		L"Z / Z_\odot = 10^{-1}",
		L"Z / Z_\odot = 0.5",
		L"Z / Z_\odot = 1.0",
	]
	Z_list = [0.0, 0.001, 0.01, 0.1, 0.5, 1.0] * MODEL.Zsun

	# Integration time [Myr]
	it = 1.0

	# Solver configuration
	args = (TRBDF2(),)
	kwargs = (dense=false, reltol=1.0e-15, verbose=false)

	# Iterators
	col_iterator = enumerate(zip(fis, yaxis_visible))
	row_iterator = enumerate(
		zip(phases, phase_labels, xaxis_visible, title_visible, ylimits),
	)

	scale = 1800 / length(fis)

	with_theme(THEME) do

		f = Figure(size=(1800, scale * length(phases)))

		for (col, (fi, yaxis_v)) in col_iterator

			base_params_01 = [[ρ, Z_list[1], 1.0, 3.0e20] for ρ in ρ_list]
			base_params_02 = [[ρ, Z_list[2], 1.0, 3.0e20] for ρ in ρ_list]
			base_params_03 = [[ρ, Z_list[3], 1.0, 3.0e20] for ρ in ρ_list]
			base_params_04 = [[ρ, Z_list[4], 1.0, 3.0e20] for ρ in ρ_list]
			base_params_05 = [[ρ, Z_list[5], 1.0, 3.0e20] for ρ in ρ_list]
			base_params_06 = [[ρ, Z_list[6], 1.0, 3.0e20] for ρ in ρ_list]

			fractions_01 = [
				solve_ode(ic(Z_list[1], fi), base_param, it)
				for base_param in base_params_01
			]

			fractions_02 = [
				solve_ode(ic(Z_list[2], fi), base_param, it)
				for base_param in base_params_02
			]

			fractions_03 = [
				solve_ode(ic(Z_list[3], fi), base_param, it)
				for base_param in base_params_03
			]

			fractions_04 = [
				solve_ode(ic(Z_list[4], fi), base_param, it)
				for base_param in base_params_04
			]

			fractions_05 = [
				solve_ode(ic(Z_list[5], fi), base_param, it)
				for base_param in base_params_05
			]

			fractions_06 = [
				solve_ode(ic(Z_list[6], fi), base_param, it)
				for base_param in base_params_06
			]

			for (row, (phase, p_label, xaxis_v, title_v, ylimit)) in row_iterator

				ax = CairoMakie.Axis(
					f[row, col];
					xlabel=L"\log_{10} \, \rho_\mathrm{cell} \,\, [\mathrm{%$(MODEL.l_u)^{-3}}]",
				    ylabel=p_label,
					title=L"f_i = %$(string(fi))",
					titlevisible=title_v,
					xminorticksvisible=xaxis_v,
		            xticksvisible=xaxis_v,
		            xlabelvisible=xaxis_v,
		            xticklabelsvisible=xaxis_v,
					yminorticksvisible=yaxis_v,
		            yticksvisible=yaxis_v,
		            ylabelvisible=yaxis_v,
		            yticklabelsvisible=yaxis_v,
					aspect=AxisAspect(1),
					limits=(nothing, ylimit),
				)

				idx = MODEL.phase_name_to_index[phase]

				fraction = getindex.(fractions_01, idx)
				lines!(ax, ρ_exp_list, fraction; label=Z_labels[1])

				fraction = getindex.(fractions_02, idx)
				lines!(ax, ρ_exp_list, fraction; label=Z_labels[2])

				fraction = getindex.(fractions_03, idx)
				lines!(ax, ρ_exp_list, fraction; label=Z_labels[3])

				fraction = getindex.(fractions_04, idx)
				lines!(ax, ρ_exp_list, fraction; label=Z_labels[4])

				fraction = getindex.(fractions_05, idx)
				lines!(ax, ρ_exp_list, fraction; label=Z_labels[5])

				fraction = getindex.(fractions_06, idx)
				lines!(ax, ρ_exp_list, fraction; label=Z_labels[6])

				if row == length(phases) && col == length(fis)
					axislegend(ax; position=:lt, nbanks=1, labelsize=32)
				end

			end

		end

		colgap!(f.layout, 35)

		f

	end
end
  ╠═╡ =#

# ╔═╡ 3bd92a5f-4efd-441e-a7d1-c3597e1f4a39
# ╠═╡ skip_as_script = true
#=╠═╡
let
	#################################################################################
	# ICs x fractions grid
	#################################################################################

	phase_labels = [L"f_i", L"f_a", L"f_m", L"f_s", L"f_Z", L"f_d"]
	phases       = [
		"ionized",
		"atomic",
		"molecular",
		"stellar",
		"metals",
		"dust",
	]

    ylimits = [
		(-0.04, 1.05),
		(-0.04, 1.05),
		(-0.04, 1.05),
		(-0.04, 1.05),
		(-0.001, 0.023),
		(-0.001, 0.041),
	]

	# Axis visibility
	xaxis_visible = [false, false, false, false, false, true]
	yaxis_visible = [true, false, false]
	title_visible = [true, false, false, false, false, false]

	# Initial conditions
	# [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
	f_i = 0.1
	ic(Z, fi) = [fi - Z, 1.0 - fi, 0.0, 0.0, 0.9 * Z, 0.1 * Z]

	# Total cell density range [mp * cm^(-3)]
	ρ_exp_list = range(-1, 2, 100)
	ρ_list     = exp10.(ρ_exp_list)

	# Metallicities [dimensionless]
	Z_labels = [
		L"Z / Z_\odot = 0.0",
		L"Z / Z_\odot = 10^{-3}",
		L"Z / Z_\odot = 10^{-2}",
		L"Z / Z_\odot = 10^{-1}",
		L"Z / Z_\odot = 0.5",
		L"Z / Z_\odot = 1.0",
	]
	Z_list = [0.0, 0.001, 0.01, 0.1, 0.5, 1.0] * MODEL.Zsun

	# Integration time [Myr]
	its = [0.5, 1.0, 10.0]

	# Solver configuration
	args = (TRBDF2(),)
	kwargs = (dense=false, reltol=1.0e-18, verbose=false)

	# Iterators
	col_iterator = enumerate(zip(its, string.(its), yaxis_visible))
	row_iterator = enumerate(
		zip(phases, phase_labels, xaxis_visible, title_visible, ylimits),
	)

	scale = 1700 / length(its)

	with_theme(THEME) do

		f = Figure(size=(1700, scale * length(phases)))

		for (col, (it, fi, yaxis_v)) in col_iterator

			base_params_01 = [[ρ, Z_list[1], 1.0, 3.0e20] for ρ in ρ_list]
			base_params_02 = [[ρ, Z_list[2], 1.0, 3.0e20] for ρ in ρ_list]
			base_params_03 = [[ρ, Z_list[3], 1.0, 3.0e20] for ρ in ρ_list]
			base_params_04 = [[ρ, Z_list[4], 1.0, 3.0e20] for ρ in ρ_list]
			base_params_05 = [[ρ, Z_list[5], 1.0, 3.0e20] for ρ in ρ_list]
			base_params_06 = [[ρ, Z_list[6], 1.0, 3.0e20] for ρ in ρ_list]

			fractions_01 = [
				solve_ode(ic(Z_list[1], f_i), base_param, it)
				for base_param in base_params_01
			]

			fractions_02 = [
				solve_ode(ic(Z_list[2], f_i), base_param, it)
				for base_param in base_params_02
			]

			fractions_03 = [
				solve_ode(ic(Z_list[3], f_i), base_param, it)
				for base_param in base_params_03
			]

			fractions_04 = [
				solve_ode(ic(Z_list[4], f_i), base_param, it)
				for base_param in base_params_04
			]

			fractions_05 = [
				solve_ode(ic(Z_list[5], f_i), base_param, it)
				for base_param in base_params_05
			]

			fractions_06 = [
				solve_ode(ic(Z_list[6], f_i), base_param, it)
				for base_param in base_params_06
			]

			for (row, (phase, p_label, xaxis_v, title_v, ylimit)) in row_iterator

				ax = CairoMakie.Axis(
					f[row, col];
					xlabel=L"\log_{10} \, \rho_\mathrm{cell} \,\, [\mathrm{%$(MODEL.l_u)^{-3}}]",
				    ylabel=p_label,
					title=L"t_f = %$it \, \mathrm{Myr}",
					titlevisible=title_v,
					xminorticksvisible=xaxis_v,
		            xticksvisible=xaxis_v,
		            xlabelvisible=xaxis_v,
		            xticklabelsvisible=xaxis_v,
					yminorticksvisible=yaxis_v,
		            yticksvisible=yaxis_v,
		            ylabelvisible=yaxis_v,
		            yticklabelsvisible=yaxis_v,
					aspect=AxisAspect(1),
					limits=(nothing, ylimit),
				)

				idx = MODEL.phase_name_to_index[phase]

				fraction = getindex.(fractions_01, idx)
				lines!(ax, ρ_exp_list, fraction; label=Z_labels[1])

				fraction = getindex.(fractions_02, idx)
				lines!(ax, ρ_exp_list, fraction; label=Z_labels[2])

				fraction = getindex.(fractions_03, idx)
				lines!(ax, ρ_exp_list, fraction; label=Z_labels[3])

				fraction = getindex.(fractions_04, idx)
				lines!(ax, ρ_exp_list, fraction; label=Z_labels[4])

				fraction = getindex.(fractions_05, idx)
				lines!(ax, ρ_exp_list, fraction; label=Z_labels[5])

				fraction = getindex.(fractions_06, idx)
				lines!(ax, ρ_exp_list, fraction; label=Z_labels[6])

				if row == 1 && col == 1
					axislegend(ax; position=:lt, nbanks=2, labelsize=30)
				end

			end

		end

		colgap!(f.layout, 35)

		mkpath("./generated_files/plots/model/")
		Makie.save("./generated_files/plots/model/fractions-vs-density.png", f)

		f

	end
end
  ╠═╡ =#

# ╔═╡ e99ef380-2646-4948-8673-436baad64d3b
# ╠═╡ skip_as_script = true
#=╠═╡
md"### Long integration time"
  ╠═╡ =#

# ╔═╡ 660aa82c-78db-469d-87cf-9325fa6a30ff
# ╠═╡ skip_as_script = true
#=╠═╡
let
	phase_labels = [L"f_i", L"f_a", L"f_m", L"f_s", L"f_Z", L"f_d"]

	# Metallicities [dimensionless]
	Z = 0.5 * MODEL.Zsun

	# Initial conditions
	# [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
	ic = [0.2 - Z, 0.8, 0.0, 0.0, 0.9 * Z, 0.1 * Z]

	# Total cell density range [mp * cm^(-3)]
	ρ = 1.0

	# Integration time [Myr]
	time_list = exp10.(range(-4, 5, 1000))

	# Integrate the system
	fractions = solve_ode(ic, [ρ, Z, 1.0, 3.0e20], time_list)

	colors = Makie.wong_colors()

	with_theme(THEME) do

		f = Figure(size=(1500, 2100))

		for (idx, (label, color)) in enumerate(zip(phase_labels, colors))

			row = ceil(Int, idx / 2)
			col = mod1(idx, 2)

			ax = CairoMakie.Axis(
				f[row, col];
				xlabel=L"\log_{10} \, t \,\, [\mathrm{%$(MODEL.t_u)}]",
		        ylabel=label,
				aspect=AxisAspect(1),
				xticks=(-4:1:6),
			)

			lines!(ax, log10.(time_list), getindex.(fractions, idx); color)

		end

		f

	end
end
  ╠═╡ =#

# ╔═╡ 05c3a97a-4096-4642-a1ae-c47283666a43
# ╠═╡ skip_as_script = true
#=╠═╡
let
	phase_labels = [L"f_i", L"f_a", L"f_m", L"f_s", L"f_Z", L"f_d"]

	# Metallicities [dimensionless]
	Z = 0.5 * MODEL.Zsun

	# Initial conditions
	# [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
	ic = [0.2 - Z, 0.8, 0.0, 0.0, 0.9 * Z, 0.1 * Z]

	# Total cell density range [mp * cm^(-3)]
	ρ = 1.0

	# Integration time [Myr]
	time_list = exp10.(range(-4, 5, 1000))

	# Integrate the system
	fractions = solve_ode(ic, [ρ, Z, 1.0, 3.0e20], time_list)

	colors = Makie.wong_colors()

	with_theme(THEME) do

		f = Figure(size=(880, 880))

		ax = CairoMakie.Axis(
			f[1,1];
			xlabel=L"\log_{10} \, t \,\, [\mathrm{%$(MODEL.t_u)}]",
		    ylabel=L"\log_{10} f",
			aspect=AxisAspect(1),
			xticks=(-4:1:5),
			yticks=(-8:1:0),
			limits=(nothing, nothing, -8.1, nothing),
		)

		for (idx, (label, color)) in enumerate(zip(phase_labels, colors))
			log_f = log10.(getindex.(fractions, idx))
			lines!(ax, log10.(time_list), log_f; color, label)
		end

		axislegend(ax; position=(0.7, 0.03), nbanks=2, labelsize=32)

		f

	end
end
  ╠═╡ =#

# ╔═╡ a6b4d852-43dd-4fe9-a0d8-f07e1c53c1f1
# ╠═╡ skip_as_script = true
#=╠═╡
let
	phase_labels = [L"f_i", L"f_a", L"f_m", L"f_s", L"f_Z", L"f_d"]

	# Metallicities [dimensionless]
	Z = 0.5 * MODEL.Zsun

	# Initial conditions
	# [fi(0), fa(0), fm(0), fs(0), fZ(0), fd(0)]
	ic = [0.2 - Z, 0.8, 0.0, 0.0, 0.9 * Z, 0.1 * Z]

	# Total cell density range [mp * cm^(-3)]
	ρ = 1.0

	# Integration time [Myr]
	it = 100.0
	time_list = exp10.(range(-4, log10(it), 1000))

	# Integrate the system
	fractions = solve_ode(ic, [ρ, Z, 1.0, 3.0e20], time_list)

	equilibrium = [
		MODEL.equilibrium(fraction, [ρ, Z, 1.0, 3.0e20, it]) for
		fraction in fractions
	]

	# Equilibrium fractions
	eq_fi = getindex.(equilibrium, 1)
	eq_fa = getindex.(equilibrium, 2)
	eq_fm = getindex.(equilibrium, 3)
	eq_fs = getindex.(equilibrium, 4)
	eq_fZ = getindex.(equilibrium, 5)
	eq_fd = getindex.(equilibrium, 6)

	with_theme(THEME) do

		f = Figure(size=(880, 880), figure_padding=(5, 10, 5, 5))

		ax = CairoMakie.Axis(
			f[1,1];
			xlabel=L"\log_{10} \, t \,\, [\mathrm{%$(MODEL.t_u)}]",
		    ylabel=L"|\mathrm{df/dt}|",
			aspect=AxisAspect(1),
			xticks=(-4:1:2),
		)

		lines!(
			ax,
			log10.(time_list),
			abs.(eq_fi);
			color=Makie.wong_colors()[1],
			label=L"\mathrm{HII}",
		)

		lines!(
			ax,
			log10.(time_list),
			abs.(eq_fa);
			color=Makie.wong_colors()[2],
			label=L"\mathrm{HI}",
		)

		lines!(
			ax,
			log10.(time_list),
			abs.(eq_fm);
			color=Makie.wong_colors()[3],
			label=L"\mathrm{H}_2",
		)

		hlines!(ax, [0.0]; color=:gray65)

		axislegend(ax; position=(0.98, 1.0), nbanks=1, labelsize=30)

		f

	end
end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
ChaosTools = "608a59af-f2a3-5ad4-90b4-758bdf3122a7"
ColorSchemes = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
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
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"

[compat]
CSV = "~0.10.15"
CairoMakie = "~0.15.4"
ChaosTools = "~3.3.1"
ColorSchemes = "~3.30.0"
DataFrames = "~1.7.0"
DifferentialEquations = "~7.16.1"
Interpolations = "~0.15.1"
Measurements = "~2.14.0"
NaNMath = "~1.1.3"
PlutoLinks = "~0.1.6"
PlutoUI = "~0.7.68"
QuadGK = "~2.11.2"
SpecialFunctions = "~2.5.1"
StatsBase = "~0.34.5"
Symbolics = "~6.45.0"
TikzPictures = "~3.5.0"
Unitful = "~1.23.1"
UnitfulAstro = "~1.2.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "98d713bca07338544af2309e538096ea9db8c1de"

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

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

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

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

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
deps = ["FillArrays", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "120e392af69350960b1d3b89d41dcc1d66543858"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.11.2"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "SIMD", "TranscodingStreams"]
git-tree-sha1 = "a8f503e8e1a5f583fbef15a8440c8c7e32185df2"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.1.0"

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

[[deps.BaseDirs]]
git-tree-sha1 = "bca794632b8a9bbe159d56bf9e31c422671b35e0"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.3.2"

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

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

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
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "22c63a303e8d6ff1a02a0933d579dbdb6e0f2925"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.15.4"

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

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "5ac098a7c8660e217ffac31dc2af0964a8c3182a"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "2.0.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "e771a63cc8b539eca78c85b0cabd9233d6c8f06f"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "a656525c8b46aa6a1c76891552ed5381bb32ae7b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.30.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

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

[[deps.Compiler]]
git-tree-sha1 = "382d79bfe72a406294faca39ef0c3cef6e6ce1f1"
uuid = "807dbc54-b67e-4c79-8afb-eafe4df6f2e1"
version = "0.1.1"

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

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "e215ba0e9a9e9377f2ed87cf3eb26840c8990585"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.3"

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

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

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

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "5620ff4ee0084a6ab7097a27ba0c19290200b037"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.4"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqRosenbrock", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack", "SymbolicIndexingInterface"]
git-tree-sha1 = "5ed2cc538b174dc2531cec319f294d5feb6b118c"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.54.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "FastPower", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "f069ea960f7a92ef5287a1b9831317dfdc19e1cc"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.179.0"

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
git-tree-sha1 = "a0c11553b6c58bfc0283d8279dd3798150b77ac5"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.12"

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

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

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

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "b3f2ff58735b5f024c392fde763f29b057e4b025"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.8"

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

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

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

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "b66970a70db13f45b7e57fbda1736e1cf72174ea"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

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

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

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

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "1f5a80f4ed9f5a4aada88fc2db456e637676414b"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.10"

    [deps.GeometryBasics.extensions]
    GeometryBasicsGeoInterfaceExt = "GeoInterface"

    [deps.GeometryBasics.weakdeps]
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "b0036b392358c80d2d2124746c2bf3d48d457938"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.4+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

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

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "dc6bed05c15523624909b3953686c5f5ffa10adc"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

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

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

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
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

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
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

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

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eac1206917768cb54957c65a615460d87b455fc1"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.1+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "e09121f4c523d8d8d9226acbed9cb66df515fcf2"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.10.4"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "f8da88993c914357031daf0023f18748ff473924"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.16.1"
weakdeps = ["FastBroadcast"]

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "ba51324b894edaf1df3ab16e2cc6bc3280a2f1a7"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.10"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "b94257a1a8737099ca40bc7271a8b374033473ed"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.10.1"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

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
git-tree-sha1 = "76627adb8c542c6b73f68d4bfd0aa71c9893a079"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.6.2"

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

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

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
git-tree-sha1 = "22e90c33c5297d6162ee54e2383584849379aa53"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "3.23.0"

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

[[deps.LoweredCodeUtils]]
deps = ["CodeTracking", "Compiler", "JuliaInterpreter"]
git-tree-sha1 = "73b98709ad811a6f81d84e105f4f695c229385ba"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.4.3"

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

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "96d73e05b6f3079df0963b66c2844d162263a896"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.4"

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
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "6e64d2321257cc52f47e193407d0659ea1b2b431"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.5"

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

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

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

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "NonlinearSolveBase", "NonlinearSolveFirstOrder", "NonlinearSolveQuasiNewton", "NonlinearSolveSpectralMethods", "PrecompileTools", "Preferences", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SparseArrays", "SparseMatrixColorings", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "d2ec18c1e4eccbb70b64be2435fc3b06fbcdc0a1"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "4.10.0"

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
git-tree-sha1 = "ee395563ae6ffaecbdf86d430440fddc779253a4"
uuid = "be0214bd-f91f-a760-ac4e-3421ce2b2da0"
version = "1.13.0"
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
git-tree-sha1 = "65101a20b135616a13625ae6f84b052ef5780363"
uuid = "5959db7a-ea39-4486-b5fe-2dd0bf03d60d"
version = "1.6.0"

[[deps.NonlinearSolveQuasiNewton]]
deps = ["ArrayInterface", "CommonSolve", "ConcreteStructs", "DiffEqBase", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLOperators", "StaticArraysCore"]
git-tree-sha1 = "3e04c917d4e3cd48b2a5091b6f76c720bb3f7362"
uuid = "9a2c21bd-3a47-402d-9113-8faf9a0ee114"
version = "1.7.0"
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

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "567515ca155d0020a45b05175449b499c63e7015"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.29+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

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

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c392fc5dd032381919e3b22dd32d6443760ce7ea"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.5.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqAdamsBashforthMoulton", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqExplicitRK", "OrdinaryDiffEqExponentialRK", "OrdinaryDiffEqExtrapolation", "OrdinaryDiffEqFIRK", "OrdinaryDiffEqFeagin", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqHighOrderRK", "OrdinaryDiffEqIMEXMultistep", "OrdinaryDiffEqLinear", "OrdinaryDiffEqLowOrderRK", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqNordsieck", "OrdinaryDiffEqPDIRK", "OrdinaryDiffEqPRK", "OrdinaryDiffEqQPRK", "OrdinaryDiffEqRKN", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqSSPRK", "OrdinaryDiffEqStabilizedIRK", "OrdinaryDiffEqStabilizedRK", "OrdinaryDiffEqSymplecticRK", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "Static", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "55c21fdb4626037cdbcb04fec3afa192345a24de"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.101.0"

[[deps.OrdinaryDiffEqAdamsBashforthMoulton]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqLowOrderRK", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "82f78099ecf4e0fa53545811318520d87e7fe0b8"
uuid = "89bda076-bce5-4f1c-845f-551c83cdda9a"
version = "1.2.0"

[[deps.OrdinaryDiffEqBDF]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqSDIRK", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "b0bbc6541ea4a27974bd67e0a10b26211cb95e58"
uuid = "6ad6398a-0878-4a85-9266-38940aa047c8"
version = "1.7.0"

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
git-tree-sha1 = "382a4bf7795eee0298221c37099db770be524bab"
uuid = "4302a76b-040a-498a-8c04-15b101fed76b"
version = "1.10.1"

[[deps.OrdinaryDiffEqExplicitRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "TruncatedStacktraces"]
git-tree-sha1 = "4dbce3f9e6974567082ce5176e21aab0224a69e9"
uuid = "9286f039-9fbf-40e8-bf65-aa933bdc4db0"
version = "1.1.0"

[[deps.OrdinaryDiffEqExponentialRK]]
deps = ["ADTypes", "DiffEqBase", "ExponentialUtilities", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "585f73f10a1b444654d739853a9328d1bb7fce6b"
uuid = "e0540318-69ee-4070-8777-9e2de6de23de"
version = "1.5.0"

[[deps.OrdinaryDiffEqExtrapolation]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastPower", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "80a636aac325c546b04e3bf20f0c80eaa0173dd4"
uuid = "becaefa8-8ca2-5cf9-886d-c06f3d2bd2c4"
version = "1.5.0"

[[deps.OrdinaryDiffEqFIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastGaussQuadrature", "FastPower", "LinearAlgebra", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "d7cbd84ba96a91e765fc20d2c3b3b1702c15eeed"
uuid = "5960d6e9-dd7a-4743-88e7-cf307b64f125"
version = "1.13.0"

[[deps.OrdinaryDiffEqFeagin]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "a7cc74d3433db98e59dc3d58bc28174c6c290adf"
uuid = "101fe9f7-ebb6-4678-b671-3a81e7194747"
version = "1.1.0"

[[deps.OrdinaryDiffEqFunctionMap]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "a4cb67794464352b69331c903cfa91a40e9a79ac"
uuid = "d3585ca7-f5d3-4ba6-8057-292ed1abd90f"
version = "1.2.0"

[[deps.OrdinaryDiffEqHighOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "3466ac9ba5121c700a17e1d5ba42757405d636d4"
uuid = "d28bc4f8-55e1-4f49-af69-84c1a99f0f58"
version = "1.2.0"

[[deps.OrdinaryDiffEqIMEXMultistep]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Reexport"]
git-tree-sha1 = "f2e7decd8b8b92a13e9d48f87780fdfecbc85708"
uuid = "9f002381-b378-40b7-97a6-27a27c83f129"
version = "1.4.0"

[[deps.OrdinaryDiffEqLinear]]
deps = ["DiffEqBase", "ExponentialUtilities", "LinearAlgebra", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "940cef72ec8799d869ff1ba3dcf47cf7758e51cf"
uuid = "521117fe-8c41-49f8-b3b6-30780b3f0fb5"
version = "1.3.0"

[[deps.OrdinaryDiffEqLowOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "2b7a3ef765c5e3e9d8eee6031da143a0b1c0805c"
uuid = "1344f307-1e59-4825-a18e-ace9aa3fa4c6"
version = "1.3.0"

[[deps.OrdinaryDiffEqLowStorageRK]]
deps = ["Adapt", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "StaticArrays"]
git-tree-sha1 = "52ec7081e65291fa5c19749312df0818db2fa1bc"
uuid = "b0944070-b475-4768-8dec-fb6eb410534d"
version = "1.3.0"

[[deps.OrdinaryDiffEqNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "PreallocationTools", "RecursiveArrayTools", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "StaticArrays"]
git-tree-sha1 = "e98aa2f8da8386bc26daeb7c9b161bc351ea6a77"
uuid = "127b3ac7-2247-4354-8eb6-78cf4e7c58e8"
version = "1.11.0"

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
git-tree-sha1 = "6548bff67665b13a60abfb33e95fcf7ae08d186d"
uuid = "af6ede74-add8-4cfd-b1df-9a4dbb109d7a"
version = "1.2.0"

[[deps.OrdinaryDiffEqRosenbrock]]
deps = ["ADTypes", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "a2f83c9b6e977c8dc5f37e0b448ad64f17c3a9c1"
uuid = "43230ef6-c299-4910-a778-202eb28ce4ce"
version = "1.12.0"

[[deps.OrdinaryDiffEqSDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "62327171ae40737b7874d4bdf70e0a213d60086a"
uuid = "2d112036-d095-4a1e-ab9a-08536f3ecdbf"
version = "1.4.0"

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
git-tree-sha1 = "8b54bcaf8634548bd13c0bd9a0522f15e5438b67"
uuid = "358294b1-0aab-51c3-aafe-ad5ab194a2ad"
version = "1.2.0"

[[deps.OrdinaryDiffEqSymplecticRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "0d3e0527149d7ece68850c51de67e99bc4477b1b"
uuid = "fa646aed-7ef9-47eb-84c4-9443fc8cbfa8"
version = "1.4.0"

[[deps.OrdinaryDiffEqTsit5]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "d0b069075f4a5e54b29e412419e5a733a83e6240"
uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a"
version = "1.2.0"

[[deps.OrdinaryDiffEqVerner]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "91f0a004785791c8c538d34a67c19cf3f7776e85"
uuid = "79d7bb75-1356-48c1-b8c0-6832512096c2"
version = "1.3.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f07c06228a1c670ae4c87d1276b92c7c597fdda0"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.35"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9dd97171646850ee607593965ce1f55063d8d3f9"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.0+0"

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
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

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

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "972089912ba299fba87671b025cd0da74f5f54f7"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.0"
weakdeps = ["ChainRulesCore", "FFTW", "Makie", "MutableArithmetics"]

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
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

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "8b3fc30bc0390abdce15f8822c889f669baed73d"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.1"

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
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
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
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "f8726bd5a8b7f5f5d3f6c0ce4793454a599b5243"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.36.0"

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

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

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

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "20ccb7e2501e9da93fe8450d01aeabf16a5f0c82"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.8.1"
weakdeps = ["Distributed"]

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

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

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "fea870727142270bdf7624ad675901a1ee3b4c87"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.1"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "57d3969321911e101c94dde1d1ff54d6308669fd"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.105.0"

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
git-tree-sha1 = "3414071e3458f3065de7fa5aed55283b236b4907"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.8"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "7d3a1519dc4d433a6b20035eaff20bde8be77c66"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.4.0"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "566c4ed301ccb2a44cbd5a27da5f885e0ed1d5df"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

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

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "818554664a2e01fc3784becb2eb3a82326a604b6"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

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
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "09d986e27a606f172c5b6cffbd8b8b2f10bf1c75"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "2.7.0"

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

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

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
git-tree-sha1 = "7bd2b8981cc57adcf5cf1add282aba2713a7058f"
uuid = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
version = "1.0.0"

    [deps.SparseConnectivityTracer.extensions]
    SparseConnectivityTracerLogExpFunctionsExt = "LogExpFunctions"
    SparseConnectivityTracerNNlibExt = "NNlib"
    SparseConnectivityTracerNaNMathExt = "NaNMath"
    SparseConnectivityTracerSpecialFunctionsExt = "SpecialFunctions"

    [deps.SparseConnectivityTracer.weakdeps]
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

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

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
git-tree-sha1 = "cbea8a6bd7bed51b1619658dec70035e07b8502f"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.14"
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
git-tree-sha1 = "c3a55a2a1e180e249a0550d30a58c700487aa7ef"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.81.0"

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

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "8ad2e38cbb812e29348719cc63580ec1dfeb9de4"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.1"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

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
git-tree-sha1 = "59ca6eddaaa9849e7de9fd1153b6faf0b1db7b80"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.42"

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
git-tree-sha1 = "8922f50f13eac24583a445968d3cfbdcc0c621ac"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.45.0"

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

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

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

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "02aca429c9885d1109e58f400c333521c13d48a0"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.4"

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

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

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

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

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

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WeakValueDicts]]
git-tree-sha1 = "98528c2610a5479f091d470967a25becfd83edd0"
uuid = "897b6980-f191-5a31-bcb0-bf3c4585e0c1"
version = "0.1.0"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

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

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4bba74fa59ab0755167ad24f98800fe5d727175b"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.12.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "07b6a107d926093898e82b3b1db657ebe33134ec"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.50+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "ccbb625a89ec6195856a50aa2b668a5c08712c94"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.4.0+0"

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
# ╠═f79c680a-7c97-4100-97e8-19190707c55b
# ╟─37653018-aaf5-42d1-a938-e05a44f18918
# ╟─6e5823eb-3614-476b-80e7-66fc27ed0392
# ╟─4b65fd9c-bd69-4e97-a28c-357c6a9c7bda
# ╠═a0ca487a-8a51-48cb-acc0-7b318c793d09
# ╟─cfacae0b-ee2e-4ba5-b31e-0fc122a3676b
# ╟─d5d5c8f9-4eb2-4332-8de5-6e20da7e528a
# ╟─b0c851e1-2cf7-48b7-b045-f5960147fb04
# ╟─78385b04-ccee-48ad-8bb3-c6b1af1bfd1e
# ╠═6162a115-fdb7-43eb-becb-0f25192001e0
# ╟─6775e416-0878-4d82-9b3a-067eeebbe3ad
# ╟─25ae9eb0-86ce-4e73-8b4b-128c33f15e4e
# ╟─dfb8b71a-b3a3-4ea0-9ea3-213c206117ca
# ╟─4ad4589f-9834-44a5-9c4b-994b7d91e442
# ╟─b2197f1a-ff9a-429c-8486-dfa27abc880d
# ╟─ad34b750-5e81-4cb8-a444-52220c116790
# ╟─b31c615a-36b3-46c0-bf8e-f8b5a77296f4
# ╟─6f341d1b-8b0f-4bc5-b18c-d1e9dc9281d5
# ╟─1d870ff9-fade-43c8-af55-6137ce051b8e
# ╟─c76cf08f-aaef-46df-ae18-9cd018141b70
# ╟─3bd92a5f-4efd-441e-a7d1-c3597e1f4a39
# ╟─e99ef380-2646-4948-8673-436baad64d3b
# ╟─660aa82c-78db-469d-87cf-9325fa6a30ff
# ╟─05c3a97a-4096-4642-a1ae-c47283666a43
# ╟─a6b4d852-43dd-4fe9-a0d8-f07e1c53c1f1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
