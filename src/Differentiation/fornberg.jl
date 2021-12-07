""" Given:
    z : location at which we want to approximate the derivative (may
    but need not be a grid point),
    x = {x₀ , x₁ , x₂, .... , xₙ} : grid point locations which the stencil is to extend over (distinct,
    otherwise arbitrary),
    n = N - 1  : one less than the number of grid points, and
    m  : highest derivative for which weights are sought,

Output Parameter
    c[0:nd,0:m] weights at grid locations x[0:n] for derivatives of order 0:m, found in c[0:n,0:m]
"""
@views function fornberg(m::Int, z::S, x::AbstractVector; nd = length(x)) where S <: Real
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
    return c
end

"""
Based on MATLAB routine in paper , An algorithm for calculating Hermite-based finite difference weights by Bengt Fornberg.
Paper found at https://www.colorado.edu/amath/sites/default/files/attached-files/2020_f_hermite-fd_ima_j_num_anal.pdf
"""
@views function fornberg_hermite(m::Int, z::S, x::AbstractVector) where S <: Real
    nd = length(x);
    c = zeros(m+2, nd);
    c[2,1] = 1;
    A = x .- x';
    s = sum(1 ./ (A + I(nd)), dims=1) .- 1;
    b = cumprod([ones(nd) A], dims = 2);
    rmat = repeat(0:m+1, 1,nd-1);
    d0 = diag(b);
    d1 = d0[1:nd-1]./d0[2:nd];
    for i ∈ 2:nd
        mn = min(i, m+1);
        c[2:mn+1,i] .= d1[i-1]*(rmat[1:mn, 1].*c[1:mn,i-1] - (x[i-1]-z)*c[2:mn+1,i-1]);
        c[2:mn+1,1:i-1] .= ((x[i]-z)*c[2:mn+1,1:i-1] - rmat[1:mn,1:i-1].*c[1:mn,1:i-1])./transpose(x[i] .- x[1:i-1]);
    end
    c2 = zeros(m+2, nd)
    cp = cumprod([1; 1:m])
    c = c[setdiff(1:end, 1), :]
    cc = c./cp
    for k=1:m+1
        c2[k+1,:] = sum(cc[1:k,:].*cc[k:-1:1,:], dims = 1);
    end
    e = c2[1:m+1,:] - transpose(x .- z).*c2[2:m+2,:];
    d = c2[2:m+2,:] + 2*s.*e;
    d = d.*cp;
    e = e.*cp;
    return c[end,:], d[end,:] ,e[end,:]
end
