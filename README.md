# Numerical-methods-HPC-Aspects

This is a repo for sharing ideas/codes/demo's related to Numerical Methods used in CFD, FEA, DL, ML and/or its HPC aspects, with a purpose of academic, professional
or lesiure satisfaction.
Be organized !!!

## A perspective on trends in Numerical Methods used in CFD, FEA, ML, and DL these days
### Differentiation
A lot depends on how we approximate differential operators or calculate gradients in these numerical methods. To put it succinctly, we have not been taught well how to approximate these differential operators properly. If one keeps up with major ways of approximating differential operators are by following way, 
- **Automatic differentiation** (a numerical approach used in automatic differentiation libraries like Tensorflow/PyTorch/Zygote) 
- **Traditional Numerical differentiation** (good old Taylor series, majorly used in CFD softwares) 
- **Symbolic differentiation** (used in Maple/Mathematica/Sage like Libraries) 

Of these, we have been exposed to traditional numerical differentiation a lot times, but a bulk of us don't know how to create automatically stencils of linear operators and apply boundary conditions on them.

In coming weeks, a small demonstration/blog of these would be added to the repo, along with useful references

For those interested in contributing to the demonstrations, please contribute your ideas [here](https://github.com/yewalenikhil65/Numerical-methods-HPC-aspects/discussions) in Julia/MATLAB like readable pseudocode syntax. (Julia is exactly MATLAB like syntaxwise, so don't be alarmed if you see actual code in it too)


## Brief Intro to High Performance Computing

## Resources: 

### Nice must-watch educational videos/online conferences on YouTube etc. ##

### Nice must-read type books/online resources/papers ##
