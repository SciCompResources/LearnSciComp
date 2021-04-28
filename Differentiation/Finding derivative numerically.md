# Finding derivatives numerically

Lets start computing function derivatives using an age-old method of finite-difference. We might as well start with the definition of derivative we learn in class of calculus. 
$$
	f'(x) = \lim_{h \to 0} \frac{f(x+h) - f(x)}{h}   \tag1
$$
## Finite-difference method
Finite-difference method is based on this very definition of derivative in equation $(1)$. The accuracy in general depends on the choice of $h$, and the choice of grid-points(or also called as choice of stencil)
Some commonly used Finite-difference approximations are as follows(written here only for first order derivatives):
- **Forward difference** 
Uses future/next grid-point data in the stencil
$$
			u'(x_i) = \frac{u(x_{i+1} )- u(x_{i} )	}{h} + O(h)  \tag2
$$
- **Backward difference**
Uses past/previous grid-point data in the stencil
$$
			u'(x_i) = \frac{u(x_{i}) - u(x_{i-1} )	}{h} + O(h)  \tag3
$$
- **Central difference**
Uses future/next and past/previous grid-point data in the stencil
$$
			u'(x_i) = \frac{u(x_{i+1} )- u(x_{i-1} )	}{2h} + O(h^2)  \tag4
$$
- **[Upwind differencing](https://en.wikipedia.org/wiki/Upwind_scheme)** Includes a bias towards choosing either side of the grid-points , depending on the characteristic 'flow' observed.

Above-mentioned schemes tell us about computing errors we face in estimating derivative of a function, namely, `truncation error` and `round-off error`, both of which depend on the choice of $h$
You choose large $h$ and you risk running into `truncation error`.
You choose too small value of $h$, and you risk into subtracting two very closely spaced numbers (also known as Subtractive cancellation). Thus, resulting in a `round-off error`. More about this can be found [here](https://en.wikipedia.org/wiki/Round-off_error#Roundoff_error_caused_by_floating-point_arithmetic)


You may choose , a small enough value of $h$ to strike a balance between `truncation error` and `round-off error`. But this also means that there is an upper-ceiling on the accuracy you can achieve.

Lets revisit what we have stated above in a form of an example below.
Consider a function (example taken from [here](https://hal.archives-ouvertes.fr/hal-01483287/document)) 

$$
	f(x) = \frac{e^{x}}{\sqrt{sin^{3}x + cos^{3}x}}   \tag5
$$
Point of interest for calculation of derivative , $x = 1.5$
Following is the demonstration of forward difference, backward difference and central difference at $x=  1.5$ in Julia.

```julia
f(x) = ℯ^x/√( (sin(x))^3 + (cos(x))^3 )  # a generic julia function

# analytical derivative/ exact derivative..(symobolic differentiation or manual way)
∂f(x) = ℯ^x*(9*sin(x) + sin(3x) + 3*cos(x) + 5*cos(3x))/
            ((√2*sin(π/4- 3x) + 3*√2*sin(x + π/4))^(3/2))

∂f(1.5)   # outputs, 4.05342789389862.. exact derivative at x=1.5
```


Finite-difference methods and steps sizes:
```julia  
# generic functions for forward, backward and central differences
forward_diff(x,h) = (f(x+h) - f(x))/h 
backward_diff(x,h) = (f(x) - f(x-h))/h
central_diff(x,h) = (f(x+h) - f(x-h))/(2h)

# Lets make a choice for h from 0.1 decreasing to 10e-16 as follows
N = 16  # number of choices
h =  fill(10.0,N).^collect(-(1:N))
```
Julia's ability of broadcasting function using `@.` operator makes things ridiculously simple. Lets Pretty our findings in tabular and plots as follows:
```julia
julia> [h  forward_diff.(1.5,h) backward_diff.(1.5,h) central_diff.(1.5,h)]
#=
16×4 Matrix{Float64}:
 0.1      4.58484  3.6303   4.10757
 0.01     4.10128  4.00664  4.05396
 0.001    4.05816  4.0487   4.05343
 0.0001   4.0539   4.05295  4.05343
 1.0e-5   4.05348  4.05338  4.05343
 1.0e-6   4.05343  4.05342  4.05343
 1.0e-7   4.05343  4.05343  4.05343
 1.0e-8   4.05343  4.05343  4.05343
 1.0e-9   4.05343  4.05343  4.05343
 1.0e-10  4.05343  4.05342  4.05343
 1.0e-11  4.05347  4.05338  4.05342
 1.0e-12  4.05365  4.05365  4.05365
 1.0e-13  4.05009  4.05009  4.05009
 1.0e-14  4.08562  3.9968   4.04121
 1.0e-15  4.44089  4.44089  4.44089
 1.0e-16  0.0      0.0      0.0
=#


# Error estimates for finite-difference schemes
err_fd = @. abs(forward_diff(1.5,h) - ∂f(1.5))
err_bd = @. abs(backward_diff(1.5,h) - ∂f(1.5))
err_cd = @. abs(central_diff(1.5,h) - ∂f(1.5))

using Plots
plotly()  # or gr().... use any backend for plot that you like
plot(h, err_fd, xaxis=:log, yaxis=:log, 
            label = "Forward difference",legend=:outertopright,lw=2)
plot!(h, err_bd, xaxis=:log, yaxis=:log, 
            label = "Backward difference",legend=:outertopright,lw=2)
plot!(h, err_cd, xaxis=:log, yaxis=:log, 
            label = "Central difference",legend=:outertopright,lw=2)
plot!(xlabel = " Step-size(h) ", ylabel = " |(Error)| ")
```
![](https://github.com/yewalenikhil65/Numerical-methods-HPC-aspects/blob/master/assets/finite_diff.png)</br>

The figure above verifies our previous statement of choosing smartly the value of $h$ to strike a balance in `truncation error` and `round-off error`
But, we can do better than this, using an equally simple technique for estimating above mentioned first order derivatives.

## Complex-step differentiation
Consider now, a *[complex analytic function](https://en.wikipedia.org/wiki/Analytic_function) $f = u + i v$, dependent on complex variable $z= x+iy$. Since $f$ is an analytic function(i.e.  differentiable function in the complex plane), we can apply the famous Cauchy-Riemann equations as follows:
$$
	\frac{∂u}{∂x} = \frac{∂v}{∂y}	\tag6
$$
$$
		\frac{∂u}{∂y} = -\frac{∂v}{∂x}	  \tag7
$$
Using first-principle definition of derivative in equation $(1)$, we could re-arrange the equation $(6)$ as
$$
	\frac{∂u}{∂x} = \lim_{h \to 0} \frac{v(x + i(y+h)) - v(x+iy)}{h}  \tag8
$$
where h, a step-size is a real number. Considering now the function defined in equation $(5)$, which is a real function (real functions are generally our interest), we can set $y=0$, $u(x) = f(x)$, and $v(x) = 0$.
In this case, we can expand a real-valued $f$ with `a step` on imaginary-axis as using a Taylor series about point $x$
$$
f(x+ih) = f(x) + ih f'(x) + (ih)^2\frac{f''(x)}{2} + ....    \tag9
$$
Taking Imaginary parts on both sides and dividing by $h$, we obtain another form of equation $(8)$ as follows,

$$
	\frac{∂f}{∂x} = \lim_{h \to 0} \frac{Im[f(x+ih)]}{h} + O(h^2)      \tag{10}
$$
Thus, we have avoided the `round-off error` due to subtractive cancellation. Equation $(10)$ is called complex-step differentiation. 
And there is little stopping us from choosing as small $h$ as we want unlike Finite-difference methods.

```julia
# generic functions for complex_step_differentiation
complex_step_diff(x,h) = imag(f(x + im*h))/h  # im is i = √-1 in julia

err_compdiff = @. abs(complex_step_diff(1.5,h) - ∂f(1.5)) # error estimate

# compare finite-diff with complex step derivatives in tabular form
[h complex_step_diff.(1.5,h) forward_diff.(1.5,h) backward_diff.(1.5,h) central_diff.(1.5,h)]
#= Output is
16×5 Matrix{Float64}:
 0.1      4.00033  4.58484  3.6303   4.10757
 0.01     4.05289  4.10128  4.00664  4.05396
 0.001    4.05342  4.05816  4.0487   4.05343
 0.0001   4.05343  4.0539   4.05295  4.05343
 1.0e-5   4.05343  4.05348  4.05338  4.05343
 1.0e-6   4.05343  4.05343  4.05342  4.05343
 1.0e-7   4.05343  4.05343  4.05343  4.05343
 1.0e-8   4.05343  4.05343  4.05343  4.05343
 1.0e-9   4.05343  4.05343  4.05343  4.05343
 1.0e-10  4.05343  4.05343  4.05342  4.05343
 1.0e-11  4.05343  4.05347  4.05338  4.05342
 1.0e-12  4.05343  4.05365  4.05365  4.05365
 1.0e-13  4.05343  4.05009  4.05009  4.05009
 1.0e-14  4.05343  4.08562  3.9968   4.04121
 1.0e-15  4.05343  4.44089  4.44089  4.44089
 1.0e-16  4.05343  0.0      0.0      0.0
=#

```
Lets see visually how complex-step differentiation fares with finite-difference schemes:

```julia
err_compdiff = @. abs(complex_step_diff(1.5,h) - ∂f(1.5))
plot!(h, err_compdiff, xaxis=:log, yaxis=:log,
       label = "complex-step differentiation", lw=2, legend=:outertopright,
       linestyle =:dashdot, ylims=(eps(Float64), 1))

```
![](https://github.com/yewalenikhil65/Numerical-methods-HPC-aspects/blob/master/assets/complex_step.png)</br>

As we can see,unlike finite-difference schemes, smaller choice of $h$ does not decrease the accuracy using complex-step differentiation method.

## Is complex-step differentiation the best method ?
Certainly not ! We would like to highlight some important observations about this method.
- Unlike finite difference methods, the method cannot be iterated to find a second order/higher order derivatives, since our $f$ is a real-valued function. This is self-explanatory from equation $(9)$
- Taylor expansion in equation $(9)$ also hints us towards an important idea of usage of Dual numbers for differentiation, when we consider $i h = \epsilon$ is real numbr and  $\epsilon^2 \to 0$ . We will explore this concept in our next topic.
