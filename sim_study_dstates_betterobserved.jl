include( "algs.jl" )
using BenchmarkTools, JLD2

d = 100

D_loc = 1.0 .+ zeros(d)

P_loc = zeros(Float64, d,d-1)
for i in 1:d
    P_loc[i,:] = repeat([1/(d-1)], d-1)
end

P_loc

Δ = 1.0

n = Int64(1000*D_loc[1]*Δ)

T, X = data_gen( D_loc, P_loc, n)
#@save "T_better_obs.txt" T
#@save "X_better_obs.txt" X

#@load "T_better_obs.txt" T
#@load "X_better_obs.txt" X

T_thin, X_thin = data_thin( Δ, T, X)
Y_thin = convert( Array{Float64}, X_thin)

# Compile
D₀, P₀ = prop_inic(repeat([1.0], d))

metr_hast_chain_D_P_mode( d, X_thin, 1.0, 2, 0.1, 1.0, 1.0, 1.0, 1.0, (D₀,P₀))
fearnhead_gibbs( d, X_thin, T_thin, 2, 1.0, 0.6666, ones(d,d-1), (D₀,P₀))
raotehHMMGibbs( d, 2, T_thin, T_thin[end], X_thin, Array( 1.0 .* I(d)[28,:] ), 1.0, 1.0, ones(d,d-1), (D₀,P₀))

# Run algorithms for k iterations
k = 100

tiemp_metro_stat  = @elapsed chain_metro_stat = metr_hast_chain_D_P_mode( d, X_thin, Δ, k, 0.05, 1.0, 2.0, 1.0, 1.0, (D₀,P₀))

tiemp_fearn_stat  = @elapsed chain_fearn_stat = fearnhead_gibbs( d, X_thin, T_thin, k, 2.0, 1.0, ones(d,d-1), (D₀,P₀))

tiemp_rao_stat  = @elapsed chain_rao_stat = raotehHMMGibbs( d, k, T_thin, T_thin[end], X_thin, Array( 1.0 .* I(d)[28,:] ), 2.0, 1.0, ones(d,d-1), (D₀,P₀))
