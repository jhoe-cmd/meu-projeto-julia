using FFTW
using Plots

# parâmetros
ħ = 1.0
m = 1.0

Nx = 512
L = 100.0
dx = L / Nx
x = (-Nx/2:Nx/2-1) .* dx

dk = 2π / L
k = (-Nx/2:Nx/2-1) .* dk

dt = 0.05
steps = 100

# pacote inicial (gaussiano)
x0 = -20.0
k0 = 2.0
σ = 2.0

ψ = exp.(-(x .- x0).^2 ./ (2σ^2)) .* exp.(im * k0 .* x)

# potencial zero
V = zeros(Nx)

# operadores
expV = exp.(-im .* V .* dt ./ (2ħ))
expK = exp.(-im .* (ħ^2 .* k.^2) ./ (2m) .* dt ./ ħ)

# evolução temporal
for i in 1:steps
    ψ .= expV .* ψ
    ψ .= ifft(expK .* fft(ψ))
    ψ .= expV .* ψ
end

# resultado
prob = abs.(ψ).^2

# salvar resultado
using DelimitedFiles
writedlm("results/probability.csv", [x prob])

# plot
plot(x, prob, title="Densidade de Probabilidade", xlabel="x", ylabel="|ψ|²")
savefig("results/plot.png")

println("Simulação concluída. Resultados em /results/")
using Pkg
Pkg.add("FFTW")
Pkg.add("Plots")
