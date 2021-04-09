# Numerical-methods-HPC-Aspects

This is a repo for sharing ideas/codes/demo's related to Numerical Methods used in CFD, FEA, DL, ML and/or its HPC aspects, with a purpose of academic, professional or lesiure satisfaction. Be organized !!!

## Guidelines for contributors
- For those interested in contributing to the demonstrations, please contribute your ideas [here](https://github.com/yewalenikhil65/Numerical-methods-HPC-aspects/discussions) in Julia/MATLAB like readable pseudocode syntax. (Julia is exactly MATLAB like syntaxwise, so don't be alarmed if you see actual code in it too). Julia and/or MATLAB should be preferable languages that you use for contribution to new ideas/concepts as demonstrations. For HPC related stuff, you can digress later into other languages if needed.
- Request to review your contributions from other contributors
- Get familiar with github first or keep a github [cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet#images) ready with you. It should help you in effectively expressing yourself in `Issues` or `PRs` or `Discussions`

## Trends in Numerical Methods used in CFD/ FEA/ ML and DL these days
Lets start with 
### 1.) Differentiation
A lot depends on how we approximate differential operators or calculate gradients in these numerical methods. To put it succinctly, we have not been taught well how to approximate these differential operators properly. If one keeps up with major ways of approximating differential operators are by following way, 
- **Automatic differentiation** (a numerical approach used in automatic differentiation libraries like Tensorflow/PyTorch/Zygote. An ongoing discussion is [here](https://github.com/yewalenikhil65/Numerical-methods-HPC-aspects/discussions/1)) 
- **Traditional Numerical differentiation** (good old Taylor series, majorly used in CFD softwares) We have been exposed to traditional numerical differentiation a lot times, but a bulk of us don't know how to create automatically stencils and Differential operators matrices and applying boundary conditions on them. (On going discussion of creating general finite difference stencils, [here](https://github.com/yewalenikhil65/Numerical-methods-HPC-aspects/discussions/3))
- **Symbolic differentiation** (used in Maple/Mathematica/Sage like Libraries) 

In coming weeks, a small demonstration/blog of these would be added to the repo, along with useful references


## Brief Intro to High Performance Computing

## Resources: 

### Nice must-watch educational videos/online conferences on YouTube etc. ##

### Nice must-read type books/online resources/papers ##
