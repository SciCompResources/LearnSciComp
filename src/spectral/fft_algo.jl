using FFTW
using Memoize

"""DFT_1 function, a naive O(𝑁²) method
"""
# DFT_1 is plain definition of DFT... O(𝑁²)
function DFT_1(x::AbstractArray; 𝑁::Int = length(x))
    X::Array{ComplexF64} = zeros(ComplexF64, 𝑁)
    for k ∈ 1:𝑁
        for n ∈ 1:𝑁
            X[k] += x[n]*ℯ^(-im*2π*(k-1)*(n-1)/𝑁)
        end
    end
    return X
end


"""DFT_2 function, a less-naive O(𝑁²/2 + 𝑁) method
"""
# DFT_2 with Little improvement over... O(𝑁²)
function DFT_2(x::AbstractArray; 𝑁::Int = length(x))
    X::Array{ComplexF64} = zeros(ComplexF64, 𝑁)
    for k ∈ 1:(𝑁÷2 + 1)
        for n ∈ 1:𝑁
            X[k] += x[n]*ℯ^(-im*2π*(k-1)*(n-1)/𝑁)
        end
        if k ≠ 1   # symmetric
            X[end - (k-2)] = conj(X[k])
        end
    end
    return X
end


"""DIT_FFT_radix2, a Cooley_Tuckey_FFT radix 2,  O(𝑁 . log(𝑁)) method
"""
# Cooley_Tuckey_FFT radix 2, based on divide and conquer
function DIT_FFT_radix2(x::AbstractArray; 𝑁::Int64 = length(x))
    X₁ = Vector{ComplexF64}()
    X₂ = Vector{ComplexF64}()
    if 𝑁 == 1
        return Array{ComplexF64}(x)
    else
        Xeven = DIT_FFT_radix2(x[1:2:end])  # recursion
        Xodd = DIT_FFT_radix2(x[2:2:end])
        for k ∈ 1:(𝑁÷2)
            push!(X₁, Xeven[k] + ℯ^(-2π*im*(k-1)/𝑁)*Xodd[k])
            push!(X₂, Xeven[k] - ℯ^(-2π*im*(k-1)/𝑁)*Xodd[k])
        end
        return [X₁; X₂]
    end
end


"""memoized a Cooley_Tuckey_FFT radix 2,  O(𝑁 . log(𝑁)) method
"""
# Memoization of Cooley_Tuckey_FFT radix 2, based on divide and conquer
@memoize function DIT_FFT_radix2_mem(x::AbstractArray; 𝑁::Int64 = length(x))
    X₁ = Vector{ComplexF64}()
    X₂ = Vector{ComplexF64}()
    if 𝑁 == 1
        return Array{ComplexF64}(x)
    else
        Xeven = DIT_FFT_radix2(x[1:2:end])  # recursion
        Xodd = DIT_FFT_radix2(x[2:2:end])
        for k ∈ 1:(𝑁÷2)
            push!(X₁, Xeven[k] + ℯ^(-2π*im*(k-1)/𝑁)*Xodd[k])
            push!(X₂, Xeven[k] - ℯ^(-2π*im*(k-1)/𝑁)*Xodd[k])
        end
        return [X₁; X₂]
    end
end
