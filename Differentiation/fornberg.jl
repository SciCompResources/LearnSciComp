#=
Given:
    z : location at which we want to approximate the derivative (may
    but need not be a grid point),
    x = {x₀ , x₁ , x₂, .... , xₙ} : grid point locations which the stencil is to extend over (distinct,
    otherwise arbitrary),
    n = N - 1  : one less than the number of grid points, and
    m  : highest derivative for which weights are sought,

Output Parameter
    c[0:nd,0:m] weights at grid locations x[0:n] for derivatives of order 0:m, found in c[0:n,0:m]
=#
function fornberg(m::Int, z::S, x::AbstractVector; nd = length(x)) where S <: Real
    nd >= m || error("number of data-points must be larger than order of derivative")
    n = nd - 1;    #
    c1 = one(S);
    c4 = x[1] - z;
    C = zeros(nd, m+1);  # first column corresponds to 0'th derivative..
                        # so m+1 column corresponds to m'th derivative
    C[1,1] = 1;
    for i in 1:n
        c2 = one(S);
        mn = min(i, m);
        c5 = c4;
        c4 = x[i+1] - z;
        for j in 0:i-1
            c3 = x[i+1] - x[j+1];
            c2 = c2*c3;
            if j == i-1
                for k in mn:-1:1
                    C[i+1,k+1] = c1*(k*C[i,k] - c5*C[i,k+1])/c2;
                end
                C[i+1,1] = -c1*c5*C[i,1]/c2;
            end
            for k in mn:-1:1
                C[j+1,k+1] = (c4*C[j+1,k+1] - k*C[j+1,k])/c3;
            end
            C[j+1,1] = c4*C[j+1,1]/c3;
        end
        c1 = c2;
    end
    c = C[:,end]; 	# last column of c gives desired row vector
    c[(n÷2) + 1] -= sum(c)
end
