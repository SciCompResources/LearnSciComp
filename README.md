# LearnSciComp.jl

LearnSciComp is under-construction Julia package, which will have generic tools required for scientific computing by an under-grad or masters student. The purpose of the package is for now purely limited to learning.
The documentation of some example-based blogs based on LearnSciComp tools can be found [here](https://scicompresources.github.io/). 

## What and How do I contribute ?
The guidelines and scope of contributions(pull-requests) that are accepted is also discussed in (https://scicompresources.github.io/). Please feel free to check [Discussions](https://github.com/SciCompResources/LearnSciComp/discussions) for knowing current plans, or feel free to create an [issue](https://github.com/SciCompResources/LearnSciComp/issues)

## How do I work and test LearnSciComp.jl locally on my computer?
1) Clone the repository through terminal with `git clone https://github.com/SciCompResources/LearnSciComp.git`
2) Test cases have been added in `test` folder. You may take a look! 
3) Start Julia in your favourite IDEs(VS Code, Juno to name a few) or inside Julia REPL
4) Change directory by setting path to `LearnSciComp` folder using cd("path").
5) Enter package mode by pressing `]`. Enter `test` to test any newly added tests-cases or newly features added.(And hope they pass!)

## How to start with LearnSciComp.jl features wih julia ?
1) Get Julia installed through suitable Julia binaries at https://julialang.org/downloads/. 
2) Start Julia in your favourite IDEs(VS Code, Juno to name a few) or inside Julia REPL
3) Enter package mode by pressing `]` 
4) Enter `add https://github.com/SciCompResources/LearnSciComp.git` while you are in package mode
5) Enter back to julia mode. Precompile by entering `using LearnSciComp`. This shall enable the features meantioned below

## Features available currently
1) Deirvative tools -  function `fornberg` to calculate weights of finite-difference formulas for arbitrary grid spacing.
```julia
"""
Central difference FD for second derivative, with order of accuracy= 2
    u''(xᵢ) = ( u(xᵢ-₁) - 2.u(xᵢ) + u(xᵢ+₁) ) / (Δx)^2

Consider,            i = 0,  Δx = 1
         therefore,  weights should be {1 , -2, 1}
"""
using LearnSciComp
order = 2;  # order of derivative you wish to approximate     
z = 0;      # location of point at which you wish to approximate the derivative
x = [-1, 0, 1];  # grid points over which the stencil is extended
julia> fornberg(order, z, x)
3-element Vector{Float64}:
  1.0
 -2.0
  1.0
```
You may also find hermite-based finite difference weights for arbitrary grid spacing by providing an optional argument `dfdx = true` as follows
```julia
order = 2;  # order of derivative you wish to approximate     
z = 0;      # location of point at which you wish to approximate the derivative, using weights of `f(x)` and `f'(x)`
x = [-1, 0, 1];  # grid points over which the stencil is extended
julia> d, e = fornberg(order, z, x;dfdx = true)
([2.0, -4.0, 2.0], [0.5, 0.0, -0.5]) 
```
where, `d` consists of weights of `f(x)` and `e` contains weights of `f'(x)`

2) Spectral tools - Discrete fourier transform function `DFT_1`, `DFT_2`, `DIT_FFT_radix2` and `DIT_FFT_radix2_mem` which are basically less-effective but accurate version of Fast Fourier transform `fft` function in [FFTW package](https://github.com/JuliaMath/FFTW.jl). 

**Note:** You may enter help mode in julia by pressing `?` and entering the name of features, say `fornberg` to see the description of the feature
