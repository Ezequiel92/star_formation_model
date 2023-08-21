<div align="center">
    <h1>✨ Star formation model</h1>
</div>

<p align="center">
    <a href="https://julialang.org"><img src="https://img.shields.io/badge/-Julia-9558B2?style=for-the-badge&logo=julia&logoColor=white"></a>
</p>

<p align="center">
    <a href="https://github.com/ezequiel92/star_formation_model/blob/main/LICENSE"><img src="https://img.shields.io/github/license/ezequiel92/star_formation_model?style=flat&logo=GNU&labelColor=2B2D2F"></a>
</p>

Implementation in [Julia](https://julialang.org) and C of our star formation model.

The goal of the models is to improve the realism of the star formation rate (SFR) in [Arepo](https://arepo-code.org/) simulations.

The notebooks in this repo, can be view as static HTML files [here](https://ezequiel92.github.io/star_formation_model/),

- `001_model.jl`: Description of the model.
- `002_codegen.jl`: Tables and C code generation.
- `003_testing.jl`: Tests and error computation.
- `004_analysis.jl`: Benchmarks.

## ⚠️ Warning

These scripts are written for my personal use and may break at any moment. So, use them at your own risk.
