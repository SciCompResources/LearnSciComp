""" Given:
    z : location at which we want to approximate the derivative (may
        but need not be a grid point),
    x = {x₀ , x₁ , x₂, .... , xₙ} : grid point locations which the stencil is to extend over (distinct,
        otherwise arbitrary),
    n = N - 1  : one less than the number of grid points, and
    m  : highest derivative for which weights are sought,
    dfdx :  optional argument to consider weights for first derivative or not
            (dfdx == true , calculates weights using hermite-based algorithm)

Output Parameter
    c : weights to approximate derivative using function values only (dfdx == false)
    d ,e : weights to approximate derivative using function and its first derivative values (dfdx == true)
"""
function fornberg(m::Int, z::S, x::T; dfdx::Bool = false) where {S <: Real, T <: AbstractVector}
    nd = length(x)
    nd > m || error("number of data-points must be larger than order of derivative")
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
    
    if dfdx == false
        return c
    else
        A = x .- x';
        s = sum(1 ./ (A + I(nd)), dims = 1) .- 1;
        cp = factorial.(0:m);
        cc = C./cp'
        c̃ = zeros(nd, m+2);
        for k in 1:m+1
           c̃[:,k+1] = sum(cc[:,1:k].*cc[:,k:-1:1], dims = 2);
        end
        E = c̃[:,1:m+1] - (x .- z).*c̃[:,2:m+2];
        D = c̃[:,2:m+2] + 2*E.*s';
        D = D.*cp';
        E = E.*cp';

        d = D[:,end];   e = E[:,end]
        return d, e
    end
end
