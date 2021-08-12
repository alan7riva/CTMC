include( "algs.jl" )
using BenchmarkTools, JLD2

using RCall

R"library(ctmcd)"

d = 5

Δ = 0.5

4/0.5

D_loc = repeat( [4/Δ], d)

P_loc = zeros(Float64, d,d-1)
for i in 1:d
    P_loc[i,:] = repeat([1/(d-1)], d-1)
end

P_loc

n = Int64(10000*D_loc[1]*Δ)

T, X = data_gen( D_loc, P_loc, n)
T_thin, X_thin = data_thin(Δ, T, X)

@save "T_5by5_scale0pt125_ob0pt5.txt" T
@save "X_5by5_scale0pt125_ob0pt5.txt" X

@load "T_5by5_scale0pt125_ob0pt5.txt" T
@load "X_5by5_scale0pt125_ob0pt5.txt" X

T_thin, X_thin = data_thin(Δ, T, X)

D₀, P₀ = prop_inic(repeat([1.0], d))

metr_hast_chain_D_P_mode( 5, X_thin, Δ, 2, 0.1, 1.0, 1.0, 1.0, 1.0, (D_loc,P_loc))

G0 = GeneratorParam( D₀, P₀ )

toy_abs = Int64.( pfeuffer_matrix( T_thin, X_thin, 5) )

@rput G0
@rput toy_abs

tiemp_metro_scale0pt125_ob0pt5  = @elapsed CHAIN_metro_scale0pt125_ob0pt5 = metr_hast_chain_D_P_mode( 5,  X_thin, Δ, 10000, 0.5, 0.5, 1.0, 1.0, 0.25, (D₀,P₀))
@save "CHAIN_metro_scale0pt125_ob0pt5.txt" CHAIN_metro_scale0pt125_ob0pt5
@save "tiemp_metro_scale0pt125_ob0pt5.txt" tiemp_metro_scale0pt125_ob0pt5

tiemp_EM_scale0pt125_ob0pt5 = @elapsed R"gmem <- gm(tm=toy_abs, te=0.5, method='EM', gmguess=G0, niter=10000)"
@save "tiemp_EM_scale0pt125_ob0pt5.txt" tiemp_EM_scale0pt125_ob0pt5

@rget gmem
EM_mat_scale0pt125_ob0pt5 = gmem[:par]

MH_mat_scale0pt125_ob0pt5 = mean( [ GeneratorParam(CHAIN_metro_scale0pt125_ob0pt5[1][i],CHAIN_metro_scale0pt125_ob0pt5[2][i]) for i in 3000:10000] )

@save "EM_mat_scale0pt125_ob0pt5.txt" EM_mat_scale0pt125_ob0pt5
@save "MH_mat_scale0pt125_ob0pt5.txt" MH_mat_scale0pt125_ob0pt5


2/0.5

D_loc = repeat( [2/Δ], d)

P_loc = zeros(Float64, d,d-1)
for i in 1:d
    P_loc[i,:] = repeat([1/(d-1)], d-1)
end

P_loc

n = Int64(10000*D_loc[1]*Δ)

T, X = data_gen( D_loc, P_loc, n)
T_thin, X_thin = data_thin(Δ, T, X)

@save "T_5by5_scale0pt25_ob0pt5.txt" T
@save "X_5by5_scale0pt25_ob0pt5.txt" X

@load "T_5by5_scale0pt25_ob0pt5.txt" T
@load "X_5by5_scale0pt25_ob0pt5.txt" X

T_thin, X_thin = data_thin(Δ, T, X)

D₀, P₀ = prop_inic(repeat([1.0], d))

G0 = GeneratorParam( D₀, P₀ )

toy_abs = Int64.( pfeuffer_matrix( T_thin, X_thin, 5) )

@rput G0
@rput toy_abs

tiemp_metro_scale0pt25_ob0pt5  = @elapsed CHAIN_metro_scale0pt25_ob0pt5 = metr_hast_chain_D_P_mode( 5,  X_thin, Δ, 10000, 0.5, 0.5, 1.0, 1.0, 0.25, (D₀,P₀))
@save "CHAIN_metro_scale0pt25_ob0pt5.txt" CHAIN_metro_scale0pt25_ob0pt5
@save "tiemp_metro_scale0pt25_ob0pt5.txt" tiemp_metro_scale0pt25_ob0pt5

tiemp_EM_scale0pt25_ob0pt5 = @elapsed R"gmem <- gm(tm=toy_abs, te=0.5, method='EM', gmguess=G0, niter=10000)"
@save "tiemp_EM_scale0pt25_ob0pt5.txt" tiemp_EM_scale0pt25_ob0pt5

@rget gmem
EM_mat_scale0pt25_ob0pt5 = gmem[:par]






MH_mat_scale0pt25_ob0pt5 = mean( [ GeneratorParam(CHAIN_metro_scale0pt25_ob0pt5[1][i],CHAIN_metro_scale0pt25_ob0pt5[2][i]) for i in 3000:10000] )

@save "EM_mat_scale0pt25_ob0pt5.txt" EM_mat_scale0pt25_ob0pt5
@save "MH_mat_scale0pt25_ob0pt5.txt" MH_mat_scale0pt25_ob0pt5

1/0.5

D_loc = repeat( [1/Δ], d)

P_loc = zeros(Float64, d,d-1)
for i in 1:d
    P_loc[i,:] = repeat([1/(d-1)], d-1)
end

P_loc

n = Int64(10000*D_loc[1]*Δ)

T, X = data_gen( D_loc, P_loc, n)
T_thin, X_thin = data_thin(Δ, T, X)

@save "T_5by5_scale0pt5_ob0pt5.txt" T
@save "X_5by5_scale0pt5_ob0pt5.txt" X

D₀, P₀ = prop_inic(repeat([1.0], d))

G0 = GeneratorParam( D₀, P₀ )

toy_abs = Int64.( pfeuffer_matrix( T_thin, X_thin, 5) )

@rput G0
@rput toy_abs

tiemp_metro_scale0pt5_ob0pt5  = @elapsed CHAIN_metro_scale0pt5_ob0pt5 = metr_hast_chain_D_P_mode( 5,  X_thin, Δ, 10000, 0.5, 0.5, 1.0, 1.0, 0.25, (D₀,P₀))
@save "CHAIN_metro_scale0pt5_ob0pt5.txt" CHAIN_metro_scale0pt5_ob0pt5
@save "tiemp_metro_scale0pt5_ob0pt5.txt" tiemp_metro_scale0pt5_ob0pt5

tiemp_EM_scale0pt5_ob0pt5 = @elapsed R"gmem <- gm(tm=toy_abs, te=0.5, method='EM', gmguess=G0, niter=10000)"
@save "tiemp_EM_scale0pt5_ob0pt5.txt" tiemp_EM_scale0pt5_ob0pt5

@rget gmem
EM_mat_scale0pt5_ob0pt5 = gmem[:par]






MH_mat_scale0pt5_ob0pt5 = mean( [ GeneratorParam(CHAIN_metro_scale0pt5_ob0pt5[1][i],CHAIN_metro_scale0pt5_ob0pt5[2][i]) for i in 3000:10000] )

@save "EM_mat_scale0pt5_ob0pt5.txt" EM_mat_scale0pt5_ob0pt5
@save "MH_mat_scale0pt5_ob0pt5.txt" MH_mat_scale0pt5_ob0pt5

################################################################################

Δ = 1.0

4/1.0

D_loc = repeat( [4/Δ], d)

P_loc = zeros(Float64, d,d-1)
for i in 1:d
    P_loc[i,:] = repeat([1/(d-1)], d-1)
end

P_loc

n = Int64(5000*D_loc[1]*Δ)

T, X = data_gen( D_loc, P_loc, n)
T_thin, X_thin = data_thin(Δ, T, X)

@save "T_5by5_scale0pt25_ob1.txt" T
@save "X_5by5_scale0pt25_ob1.txt" X

@load "T_5by5_scale0pt25_ob1.txt" T
@load "X_5by5_scale0pt25_ob1.txt" X

T_thin, X_thin = data_thin(Δ, T, X)

D₀, P₀ = prop_inic(repeat([1.0], d))

G0 = GeneratorParam( D₀, P₀ )

toy_abs = Int64.( pfeuffer_matrix( T_thin, X_thin, 5) )

@rput G0
@rput toy_abs

tiemp_metro_scale0pt25_ob1  = @elapsed CHAIN_metro_scale0pt25_ob1 = metr_hast_chain_D_P_mode( 5,  X_thin, Δ, 10000, 0.5, 0.5, 1.0, 1.0, 0.25, (D₀,P₀))
@save "CHAIN_metro_scale0pt25_ob1.txt" CHAIN_metro_scale0pt25_ob1
@save "tiemp_metro_scale0pt25_ob1.txt" tiemp_metro_scale0pt25_ob1

tiemp_EM_scale0pt25_ob1 = @elapsed R"gmem <- gm(tm=toy_abs, te=1.0, method='EM', gmguess=G0, niter=10000)"
@save "tiemp_EM_scale0pt25_ob1.txt" tiemp_EM_scale0pt25_ob1

@rget gmem
EM_mat_scale0pt25_ob1 = gmem[:par]






MH_mat_scale0pt25_ob1 = mean( [ GeneratorParam(CHAIN_metro_scale0pt25_ob1[1][i],CHAIN_metro_scale0pt25_ob1[2][i]) for i in 3000:10000] )

@save "EM_mat_scale0pt25_ob1.txt" EM_mat_scale0pt25_ob1
@save "MH_mat_scale0pt25_ob1.txt" MH_mat_scale0pt25_ob1

2/1.0

D_loc = repeat( [2/Δ], d)

P_loc = zeros(Float64, d,d-1)
for i in 1:d
    P_loc[i,:] = repeat([1/(d-1)], d-1)
end

P_loc

n = Int64(5000*D_loc[1]*Δ)

T, X = data_gen( D_loc, P_loc, n)
T_thin, X_thin = data_thin(Δ, T, X)

@save "T_5by5_scale0pt5_ob1.txt" T
@save "X_5by5_scale0pt5_ob1.txt" X

@load "T_5by5_scale0pt5_ob1.txt" T
@load "X_5by5_scale0pt5_ob1.txt" X

T_thin, X_thin = data_thin(Δ, T, X)

D₀, P₀ = prop_inic(repeat([1.0], d))

G0 = GeneratorParam( D₀, P₀ )

toy_abs = Int64.( pfeuffer_matrix( T_thin, X_thin, 5) )

@rput G0
@rput toy_abs

tiemp_metro_scale0pt5_ob1  = @elapsed CHAIN_metro_scale0pt5_ob1 = metr_hast_chain_D_P_mode( 5,  X_thin, Δ, 10000, 0.5, 0.5, 1.0, 1.0, 0.25, (D₀,P₀))
@save "CHAIN_metro_scale0pt5_ob1.txt" CHAIN_metro_scale0pt5_ob1
@save "tiemp_metro_scale0pt5_ob1.txt" tiemp_metro_scale0pt5_ob1

tiemp_EM_scale0pt5_ob1 = @elapsed R"gmem <- gm(tm=toy_abs, te=1.0, method='EM', gmguess=G0, niter=10000)"
@save "tiemp_EM_scale0pt5_ob1.txt" tiemp_EM_scale0pt5_ob1

@rget gmem
EM_mat_scale0pt5_ob1 = gmem[:par]






MH_mat_scale0pt5_ob1 = mean( [ GeneratorParam(CHAIN_metro_scale0pt5_ob1[1][i],CHAIN_metro_scale0pt5_ob1[2][i]) for i in 3000:10000] )

@save "EM_mat_scale0pt5_ob1.txt" EM_mat_scale0pt5_ob1
@save "MH_mat_scale0pt5_ob1.txt" MH_mat_scale0pt5_ob1

1/1.0

D_loc = repeat( [1/Δ], d)

P_loc = zeros(Float64, d,d-1)
for i in 1:d
    P_loc[i,:] = repeat([1/(d-1)], d-1)
end

P_loc

n = Int64(5000*D_loc[1]*Δ)

T, X = data_gen( D_loc, P_loc, n)
T_thin, X_thin = data_thin(Δ, T, X)

@save "T_5by5_scale1_ob1.txt" T
@save "X_5by5_scale1_ob1.txt" X

@load "T_5by5_scale1_ob1.txt" T
@load "X_5by5_scale1_ob1.txt" X

T_thin, X_thin = data_thin(Δ, T, X)

D₀, P₀ = prop_inic(repeat([1.0], d))

G0 = GeneratorParam( D₀, P₀ )

toy_abs = Int64.( pfeuffer_matrix( T_thin, X_thin, 5) )

@rput G0
@rput toy_abs

tiemp_metro_scale1_ob1  = @elapsed CHAIN_metro_scale1_ob1 = metr_hast_chain_D_P_mode( 5,  X_thin, Δ, 10000, 0.5, 0.5, 1.0, 1.0, 0.25, (D₀,P₀))
@save "CHAIN_metro_scale1_ob1.txt" CHAIN_metro_scale1_ob1
@save "tiemp_metro_scale1_ob1.txt" tiemp_metro_scale1_ob1

tiemp_EM_scale1_ob1  = @elapsed R"gmem <- gm(tm=toy_abs, te=1.0, method='EM', gmguess=G0, niter=10000)"
@save "tiemp_EM_scale1_ob1.txt" tiemp_EM_scale1_ob1

@rget gmem
EM_mat_scale1_ob1 = gmem[:par]






MH_mat_scale1_ob1 = mean( [ GeneratorParam(CHAIN_metro_scale1_ob1[1][i],CHAIN_metro_scale1_ob1[2][i]) for i in 3000:10000] )

@save "EM_mat_scale1_ob1.txt" EM_mat_scale1_ob1
@save "MH_mat_scale1_ob1.txt" MH_mat_scale1_ob1

################################################################################

Δ = 2.0

4/2.0

D_loc = repeat( [4/Δ], d)

P_loc = zeros(Float64, d,d-1)
for i in 1:d
    P_loc[i,:] = repeat([1/(d-1)], d-1)
end

P_loc

n = Int64(5000*D_loc[1]*Δ)

T, X = data_gen( D_loc, P_loc, n)
T_thin, X_thin = data_thin(Δ, T, X)

@save "T_5by5_scale0pt5_ob2.txt" T
@save "X_5by5_scale0pt5_ob2.txt" X

@load "T_5by5_scale0pt5_ob2.txt" T
@load "X_5by5_scale0pt5_ob2.txt" X

T_thin, X_thin = data_thin(Δ, T, X)

D₀, P₀ = prop_inic(repeat([1.0], d))

G0 = GeneratorParam( D₀, P₀ )

toy_abs = Int64.( pfeuffer_matrix( T_thin, X_thin, 5) )

@rput G0
@rput toy_abs

tiemp_metro_scale0pt5_ob2  = @elapsed CHAIN_metro_scale0pt5_ob2 = metr_hast_chain_D_P_mode( 5,  X_thin, Δ, 10000, 0.5, 0.5, 1.0, 1.0, 0.25, (D₀,P₀))
@save "CHAIN_metro_scale0pt5_ob2.txt" CHAIN_metro_scale0pt5_ob2
@save "tiemp_metro_scale0pt5_ob2.txt" tiemp_metro_scale0pt5_ob2

tiemp_EM_scale0pt5_ob2  = @elapsed R"gmem <- gm(tm=toy_abs, te=2.0, method='EM', gmguess=G0, niter=10000)"
@save "tiemp_EM_scale0pt5_ob2.txt" tiemp_EM_scale0pt5_ob2


@rget gmem


EM_mat_scale0pt5_ob2 = gmem[:par]






MH_mat_scale0pt5_ob2 = mean( [ GeneratorParam(CHAIN_metro_scale0pt5_ob2[1][i],CHAIN_metro_scale0pt5_ob2[2][i]) for i in 3000:10000] )

@save "EM_mat_scale0pt5_ob2.txt" EM_mat_scale0pt5_ob2
@save "MH_mat_scale0pt5_ob2.txt" MH_mat_scale0pt5_ob2

2/2.0

D_loc = repeat( [2/Δ], d)

P_loc = zeros(Float64, d,d-1)
for i in 1:d
    P_loc[i,:] = repeat([1/(d-1)], d-1)
end

P_loc

n = Int64(5000*D_loc[1]*Δ)

T, X = data_gen( D_loc, P_loc, n)
T_thin, X_thin = data_thin(Δ, T, X)

@save "T_5by5_scale1_ob2.txt" T
@save "X_5by5_scale1_ob2.txt" X

@load "T_5by5_scale1_ob2.txt" T
@load "X_5by5_scale1_ob2.txt" X

T_thin, X_thin = data_thin(Δ, T, X)

D₀, P₀ = prop_inic(repeat([1.0], d))

G0 = GeneratorParam( D₀, P₀ )

toy_abs = Int64.( pfeuffer_matrix( T_thin, X_thin, 5) )

@rput G0
@rput toy_abs

tiemp_metro_scale1_ob2  = @elapsed CHAIN_metro_scale1_ob2 = metr_hast_chain_D_P_mode( 5,  X_thin, Δ, 10000, 0.5, 0.5, 1.0, 1.0, 0.25, (D₀,P₀))
@save "CHAIN_metro_scale1_ob2.txt" CHAIN_metro_scale1_ob2
@save "tiemp_metro_scale1_ob2.txt" tiemp_metro_scale1_ob2

tiemp_EM_scale1_ob2  = @elapsed R"gmem <- gm(tm=toy_abs, te=2.0, method='EM', gmguess=G0, niter=10000)"
@save "tiemp_EM_scale1_ob2.txt" tiemp_EM_scale1_ob2


@rget gmem


EM_mat_scale1_ob2 = gmem[:par]






MH_mat_scale1_ob2 = mean( [ GeneratorParam(CHAIN_metro_scale1_ob2[1][i],CHAIN_metro_scale1_ob2[2][i]) for i in 3000:10000] )

@save "EM_mat_scale1_ob2.txt" EM_mat_scale1_ob2
@save "MH_mat_scale1_ob2.txt" MH_mat_scale1_ob2

1/2.0

D_loc = repeat( [1/Δ], d)

P_loc = zeros(Float64, d,d-1)
for i in 1:d
    P_loc[i,:] = repeat([1/(d-1)], d-1)
end

P_loc

n = Int64(5000*D_loc[1]*Δ)

T, X = data_gen( D_loc, P_loc, n)
T_thin, X_thin = data_thin(Δ, T, X)

@save "T_5by5_scale2_ob2.txt" T
@save "X_5by5_scale2_ob2.txt" X

@load "T_5by5_scale2_ob2.txt" T
@load "X_5by5_scale2_ob2.txt" X

T_thin, X_thin = data_thin(Δ, T, X)

D₀, P₀ = prop_inic(repeat([1.0], d))

G0 = GeneratorParam( D₀, P₀ )

toy_abs = Int64.( pfeuffer_matrix( T_thin, X_thin, 5) )

@rput G0
@rput toy_abs

tiemp_metro_scale2_ob2  = @elapsed CHAIN_metro_scale2_ob2 = metr_hast_chain_D_P_mode( 5,  X_thin, Δ, 10000, 0.5, 0.5, 1.0, 1.0, 0.25, (D₀,P₀))
@save "CHAIN_metro_scale2_ob2.txt" CHAIN_metro_scale2_ob2
@save "tiemp_metro_scale2_ob2.txt" tiemp_metro_scale2_ob2

tiemp_EM_scale2_ob2  = @elapsed R"gmem <- gm(tm=toy_abs, te=2.0, method='EM', gmguess=G0, niter=10000)"
@save "tiemp_EM_scale2_ob2.txt" tiemp_EM_scale2_ob2


@rget gmem


EM_mat_scale2_ob2 = gmem[:par]






MH_mat_scale2_ob2 = mean( [ GeneratorParam(CHAIN_metro_scale2_ob2[1][i],CHAIN_metro_scale2_ob2[2][i]) for i in 3000:10000] )

@save "EM_mat_scale2_ob2.txt" EM_mat_scale2_ob2
@save "MH_mat_scale2_ob2.txt" MH_mat_scale2_ob2


################################################################################
# Fearnhead
################################################################################

@load "T_5by5_scale0pt125_ob0pt5.txt" T
@load "X_5by5_scale0pt125_ob0pt5.txt" X

Δ = 0.5

T_thin, X_thin = data_thin(Δ, T, X)

@load "CHAIN_metro_scale0pt125_ob0pt5.txt"

D₀ = CHAIN_metro_scale0pt125_ob0pt5[1][1]
P₀ = CHAIN_metro_scale0pt125_ob0pt5[2][1]

fearnhead_gibbs( 5, X_thin, T_thin, 3, 1.0, 1.0, 0.25.*ones(5,4), (D₀,P₀))

tiemp_fearn_scale0pt125_ob0pt5  = @elapsed CHAIN_fearn_scale0pt125_ob0pt5 = fearnhead_gibbs( 5, X_thin, T_thin, 10000, 1.0, 1.0, 0.25.*ones(5,4), (D₀,P₀))
@save "CHAIN_fearn_scale0pt125_ob0pt5.txt" CHAIN_fearn_scale0pt125_ob0pt5
@save "tiemp_fearn_scale0pt125_ob0pt5.txt" tiemp_fearn_scale0pt125_ob0pt5

GF_mat_scale0pt125_ob0pt5 = mean( [ CHAIN_fearn_scale0pt125_ob0pt5[i] for i in 3000:10000] )
@save "GF_mat_scale0pt125_ob0pt5.txt" GF_mat_scale0pt125_ob0pt5

@load "T_5by5_scale0pt25_ob0pt5.txt" T
@load "X_5by5_scale0pt25_ob0pt5.txt" X

T_thin, X_thin = data_thin(Δ, T, X)

@load "CHAIN_metro_scale0pt25_ob0pt5.txt"

D₀ = CHAIN_metro_scale0pt25_ob0pt5[1][1]
P₀ = CHAIN_metro_scale0pt25_ob0pt5[2][1]

tiemp_fearn_scale0pt25_ob0pt5  = @elapsed CHAIN_fearn_scale0pt25_ob0pt5 = fearnhead_gibbs( 5, X_thin, T_thin, 10000, 1.0, 1.0, 0.25.*ones(5,4), (D₀,P₀))
@save "CHAIN_fearn_scale0pt25_ob0pt5.txt" CHAIN_fearn_scale0pt25_ob0pt5
@save "tiemp_fearn_scale0pt25_ob0pt5.txt" tiemp_fearn_scale0pt25_ob0pt5

GF_mat_scale0pt25_ob0pt5 = mean( [ CHAIN_fearn_scale0pt25_ob0pt5[i] for i in 3000:10000] )
@save "GF_mat_scale0pt25_ob0pt5.txt" GF_mat_scale0pt25_ob0pt5

@load "T_5by5_scale0pt5_ob0pt5.txt" T
@load "X_5by5_scale0pt5_ob0pt5.txt" X

T_thin, X_thin = data_thin(Δ, T, X)

@load "CHAIN_metro_scale0pt5_ob0pt5.txt"

D₀ = CHAIN_metro_scale0pt5_ob0pt5[1][1]
P₀ = CHAIN_metro_scale0pt5_ob0pt5[2][1]

tiemp_fearn_scale0pt5_ob0pt5  = @elapsed CHAIN_fearn_scale0pt5_ob0pt5 = fearnhead_gibbs( 5, X_thin, T_thin, 10000, 1.0, 1.0, 0.25.*ones(5,4), (D₀,P₀))
@save "CHAIN_fearn_scale0pt5_ob0pt5.txt" CHAIN_fearn_scale0pt5_ob0pt5
@save "tiemp_fearn_scale0pt5_ob0pt5.txt" tiemp_fearn_scale0pt5_ob0pt5

GF_mat_scale0pt5_ob0pt5 = mean( [ CHAIN_fearn_scale0pt5_ob0pt5[i] for i in 3000:10000] )
@save "GF_mat_scale0pt5_ob0pt5.txt" GF_mat_scale0pt5_ob0pt5

@load "T_5by5_scale0pt25_ob1.txt" T
@load "X_5by5_scale0pt25_ob1.txt" X

Δ = 1.0

T_thin, X_thin = data_thin(Δ, T, X)

@load "CHAIN_metro_scale0pt25_ob1.txt"

D₀ = CHAIN_metro_scale0pt25_ob1[1][1]
P₀ = CHAIN_metro_scale0pt25_ob1[2][1]

tiemp_fearn_scale0pt25_ob1  = @elapsed CHAIN_fearn_scale0pt25_ob1 = fearnhead_gibbs( 5, X_thin, T_thin, 10000, 1.0, 1.0, 0.25.*ones(5,4), (D₀,P₀))
@save "CHAIN_fearn_scale0pt25_ob1.txt" CHAIN_fearn_scale0pt25_ob1
@save "tiemp_fearn_scale0pt25_ob1.txt" tiemp_fearn_scale0pt25_ob1

GF_mat_scale0pt25_ob1 = mean( [ CHAIN_fearn_scale0pt25_ob1[i] for i in 3000:10000] )
@save "GF_mat_scale0pt25_ob1.txt" GF_mat_scale0pt25_ob1

@load "T_5by5_scale0pt5_ob1.txt" T
@load "X_5by5_scale0pt5_ob1.txt" X

T_thin, X_thin = data_thin(Δ, T, X)

@load "CHAIN_metro_scale0pt5_ob1.txt"

D₀ = CHAIN_metro_scale0pt5_ob1[1][1]
P₀ = CHAIN_metro_scale0pt5_ob1[2][1]

tiemp_fearn_scale0pt5_ob1  = @elapsed CHAIN_fearn_scale0pt5_ob1 = fearnhead_gibbs( 5, X_thin, T_thin, 10000, 1.0, 1.0, 0.25.*ones(5,4), (D₀,P₀))
@save "CHAIN_fearn_scale0pt5_ob1.txt" CHAIN_fearn_scale0pt5_ob1
@save "tiemp_fearn_scale0pt5_ob1.txt" tiemp_fearn_scale0pt5_ob1

GF_mat_scale0pt5_ob1 = mean( [ CHAIN_fearn_scale0pt5_ob1[i] for i in 3000:10000] )
@save "GF_mat_scale0pt5_ob1.txt" GF_mat_scale0pt5_ob1

@load "T_5by5_scale1_ob1.txt" T
@load "X_5by5_scale1_ob1.txt" X

T_thin, X_thin = data_thin(Δ, T, X)

@load "CHAIN_metro_scale1_ob1.txt"

D₀ = CHAIN_metro_scale1_ob1[1][1]
P₀ = CHAIN_metro_scale1_ob1[2][1]

tiemp_fearn_scale1_ob1  = @elapsed CHAIN_fearn_scale1_ob1 = fearnhead_gibbs( 5, X_thin, T_thin, 10000, 1.0, 1.0, 0.25.*ones(5,4), (D₀,P₀))
@save "CHAIN_fearn_scale1_ob1.txt" CHAIN_fearn_scale1_ob1
@save "tiemp_fearn_scale1_ob1.txt" tiemp_fearn_scale1_ob1

GF_mat_scale1_ob1 = mean( [ CHAIN_fearn_scale1_ob1[i] for i in 3000:10000] )
@save "GF_mat_scale1_ob1.txt" GF_mat_scale1_ob1

@load "T_5by5_scale0pt5_ob2.txt" T
@load "X_5by5_scale0pt5_ob2.txt" X

Δ = 2.0

T_thin, X_thin = data_thin(Δ, T, X)

@load "CHAIN_metro_scale0pt5_ob2.txt"

D₀ = CHAIN_metro_scale0pt5_ob2[1][1]
P₀ = CHAIN_metro_scale0pt5_ob2[2][1]

tiemp_fearn_scale0pt5_ob2  = @elapsed CHAIN_fearn_scale0pt5_ob2 = fearnhead_gibbs( 5, X_thin, T_thin, 10000, 1.0, 1.0, 0.25.*ones(5,4), (D₀,P₀))
@save "CHAIN_fearn_scale0pt5_ob2.txt" CHAIN_fearn_scale0pt5_ob2
@save "tiemp_fearn_scale0pt5_ob2.txt" tiemp_fearn_scale0pt5_ob2

GF_mat_scale0pt5_ob2 = mean( [ CHAIN_fearn_scale0pt5_ob2[i] for i in 3000:10000] )
@save "GF_mat_scale0pt5_ob2.txt" GF_mat_scale0pt5_ob2

@load "T_5by5_scale1_ob2.txt" T
@load "X_5by5_scale1_ob2.txt" X

T_thin, X_thin = data_thin(Δ, T, X)

@load "CHAIN_metro_scale1_ob2.txt"

D₀ = CHAIN_metro_scale1_ob2[1][1]
P₀ = CHAIN_metro_scale1_ob2[2][1]

tiemp_fearn_scale1_ob2  = @elapsed CHAIN_fearn_scale1_ob2 = fearnhead_gibbs( 5, X_thin, T_thin, 10000, 1.0, 1.0, 0.25.*ones(5,4), (D₀,P₀))
@save "CHAIN_fearn_scale1_ob2.txt" CHAIN_fearn_scale1_ob2
@save "tiemp_fearn_scale1_ob2.txt" tiemp_fearn_scale1_ob2

GF_mat_scale1_ob2 = mean( [ CHAIN_fearn_scale1_ob2[i] for i in 3000:10000] )
@save "GF_mat_scale1_ob2.txt" GF_mat_scale1_ob2

@load "T_5by5_scale2_ob2.txt" T
@load "X_5by5_scale2_ob2.txt" X

T_thin, X_thin = data_thin(Δ, T, X)

@load "CHAIN_metro_scale2_ob2.txt"

D₀ = CHAIN_metro_scale2_ob2[1][1]
P₀ = CHAIN_metro_scale2_ob2[2][1]

tiemp_fearn_scale2_ob2  = @elapsed CHAIN_fearn_scale2_ob2 = fearnhead_gibbs( 5, X_thin, T_thin, 10000, 1.0, 1.0, 0.25.*ones(5,4), (D₀,P₀))
@save "CHAIN_fearn_scale2_ob2.txt" CHAIN_fearn_scale2_ob2
@save "tiemp_fearn_scale2_ob2.txt" tiemp_fearn_scale2_ob2

GF_mat_scale2_ob2 = mean( [ CHAIN_fearn_scale2_ob2[i] for i in 3000:10000] )
@save "GF_mat_scale2_ob2.txt" GF_mat_scale2_ob2
