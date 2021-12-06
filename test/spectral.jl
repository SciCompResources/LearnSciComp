@info "Starting spectral tests"
using LearnSciComp
using FFTW

@info "Testing DFT functions"
x = rand(2^10);
𝑁 = length(x);
@testset "DFT algorithms" begin
    @test DFT_1(x; 𝑁 = 𝑁) ≈ fft(x)
    @test DFT_2(x; 𝑁 = 𝑁) ≈ fft(x)
    @test DIT_FFT_radix2(x; 𝑁 = 𝑁) ≈ fft(x)
    @test DIT_FFT_radix2_mem(x; 𝑁 = 𝑁) ≈ fft(x)
end
