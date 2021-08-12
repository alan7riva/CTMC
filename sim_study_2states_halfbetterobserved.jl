include( "algs.jl" )
using BenchmarkTools, JLD2

d = 2

D_loc = [2.0,1.0]

P_loc = zeros(Float64, d,d-1)
for i in 1:d
    P_loc[i,:] = repeat([1/(d-1)], d-1)
end

P_loc

Δ = 1.0

n = Int64(1000*D_loc[1]*Δ)

# T, X = data_gen( D_loc, P_loc, n)
# @save "T_better_obs.txt" T
# @save "X_better_obs.txt" X

@load "T_better_obs.txt" T
@load "X_better_obs.txt" X

T_thin, X_thin = data_thin( Δ, T, X)
Y_thin = convert( Array{Float64}, X_thin)

n_ch_1 = Num_visits_from_j_to_k( 1, 2, T_thin, X_thin)
n_ch_2 = Num_visits_from_j_to_k( 2, 1, T_thin, X_thin)
n_eq_1 = Num_visits_from_j_to_k( 1, 1, T_thin, X_thin)
n_eq_2 = Num_visits_from_j_to_k( 2, 2, T_thin, X_thin)

# Compile
D₀, P₀ = prop_inic(repeat([1.0], 2))

metr_hast_chain_D_P_mode( 2, X_thin, 6.0, 2, 0.1, 1.0, 1.0, 1.0, 1.0, (D₀,P₀))
fearnhead_gibbs( 2, X_thin, T_thin, 2, 1.0, 0.6666, ones(2,1), (D₀,P₀))
raotehHMMGibbs( 2, 2, T_thin, T_thin[end], X_thin, [1.0,0.0], 1.0, 1.0, ones(2,1), (D₀,P₀))

################################################################################
########### Before convergence analysis ########################################
################################################################################

k = 5000

D₀vec = [ D_inic(d) for _ in 1:4 ]

# chains_metro, tiemps_metro = stationary_diagnost_metro_modeprop(4,d,3,0.05, 1.0, 2.0, 1.0, 1.0, D₀vec)

#@save "chains_metro.txt" chains_metro
#@save "tiemps_metro.txt" tiemps_metro

@load "chains_metro.txt"
@load "tiemps_metro.txt"

potential_scale_reduction(chains_metro[1])
potential_scale_reduction(chains_metro[2])

# chains_rao, tiemps_rao = stationary_diagnost_raoteh(4,d, 3, 2.0, 1.0, 1.0, D₀vec)

#@save "chains_rao.txt" chains_rao
#@save "tiemps_rao.txt" tiemps_rao

@load "chains_rao.txt"
@load "tiemps_rao.txt"

potential_scale_reduction(chains_rao[1])

potential_scale_reduction(chains_rao[2])

# chains_fearn, tiemps_fearn = stationary_diagnost_fearnhead(4,d,3, 1.0, 1.0, 1.0, D₀vec)

#@save "chains_fearn.txt" chains_fearn
#@save "tiemps_fearn.txt" tiemps_fearn

@load "chains_fearn.txt"
@load "tiemps_fearn.txt"

potential_scale_reduction(chains_fearn[1])
potential_scale_reduction(chains_fearn[2])

length(chains_metro[1][1])
length(chains_rao[1][1])
length(chains_fearn[1][1])

mean(tiemps_metro)
mean(tiemps_rao)
mean(tiemps_fearn)

ESSMetroλ1BeforeStat = [ effective_sample_size(chains_metro[1][i]) for i in 1:4 ]
mean(ESSMetroλ1BeforeStat)
ESSMetroλ2BeforeStat = [ effective_sample_size(chains_metro[2][i]) for i in 1:4 ]
mean(ESSMetroλ2BeforeStat)
ESSRaoλ1BeforeStat = [ effective_sample_size(chains_rao[1][i]) for i in 1:4 ]
mean(ESSRaoλ1BeforeStat)
ESSRaoλ2BeforeStat = [ effective_sample_size(chains_rao[2][i]) for i in 1:4 ]
mean(ESSRaoλ2BeforeStat)
ESSFearnλ1BeforeStat = [ effective_sample_size(chains_fearn[1][i]) for i in 1:4 ]
mean(ESSFearnλ1BeforeStat)
ESSFearnλ2BeforeStat = [ effective_sample_size(chains_fearn[2][i]) for i in 1:4 ]
mean(ESSFearnλ2BeforeStat)

ESSperSecMetroλ1BeforeStat = [ effective_sample_size(chains_metro[1][i])/tiemps_metro[i] for i in 1:4 ]





mean(ESSperSecMetroλ1BeforeStat)

ESSperSecMetroλ2BeforeStat = [ effective_sample_size(chains_metro[2][i])/tiemps_metro[i] for i in 1:4 ]





mean(ESSperSecMetroλ2BeforeStat)


ESSperSecRaoλ1BeforeStat = [ effective_sample_size(chains_rao[1][i])/tiemps_rao[i] for i in 1:4 ]





mean(ESSperSecRaoλ1BeforeStat)

ESSperSecRaoλ2BeforeStat = [ effective_sample_size(chains_rao[2][i])/tiemps_rao[i] for i in 1:4 ]





mean(ESSperSecRaoλ2BeforeStat)

ESSperSecFearnλ1BeforeStat = [ effective_sample_size(chains_fearn[1][i])/tiemps_rao[i] for i in 1:4 ]





mean(ESSperSecFearnλ1BeforeStat)

ESSperSecFearnλ2BeforeStat = [ effective_sample_size(chains_fearn[2][i])/tiemps_rao[i] for i in 1:4 ]





mean(ESSperSecFearnλ2BeforeStat)

################################################################################
########### After convergence analysis #########################################
################################################################################

@load "chains_metro.txt"
@load "chains_fearn.txt"
@load "chains_rao.txt"

rand_chain = rand(Categorical(4))
D₀metro_stat = [ chains_metro[1][rand_chain][end], chains_metro[2][rand_chain][end] ]


tiemp_metro_stat  = @elapsed chain_metro_stat = metr_hast_chain_D_P_mode( 2, X_thin, Δ, 5000, 0.05, 1.0, 1.0, 1.0, 1.0, (D₀metro_stat,P₀))

λ_1_metro_stat_5 = [ chain_metro_stat[1][i][1] for i in 1:length(chain_metro_stat[1]) ]
λ_2_metro_stat_5 = [ chain_metro_stat[1][i][2] for i in 1:length(chain_metro_stat[1]) ]

rand_chain = rand(Categorical(4))
D₀rao_stat = [ chains_rao[1][rand_chain][end], chains_metro[2][rand_chain][end] ]


tiemp_rao_stat  = @elapsed chain_rao_stat = raotehHMMGibbs( 2, k, T_thin, T_thin[end], X_thin, [1.0,0.0], 1.0, 1.0, ones(2,1), (D₀rao_stat,P₀))

λ_1_rao_stat_5 = [ chain_rao_stat[1][i][1] for i in 1:length(chain_rao_stat[1]) ]
λ_2_rao_stat_5 = [ chain_rao_stat[1][i][2] for i in 1:length(chain_rao_stat[1]) ]

rand_chain = rand(Categorical(4))
D₀fearn_stat = [ chains_fearn[1][rand_chain][end], chains_fearn[2][rand_chain][end] ]

tiemp_fearn_stat  = @elapsed chain_fearn_stat = fearnhead_gibbs( 2, X_thin, T_thin, k, 1.0, 1.0, ones(2,1), (D₀fearn_stat,P₀))

λ_1_fearn_stat_5 = [ -chain_fearn_stat[i][1,1] for i in 1:length(chain_fearn_stat) ]
λ_2_fearn_stat_5 = [ -chain_fearn_stat[i][2,2] for i in 1:length(chain_fearn_stat) ]

tiemp_metro_stat
tiemp_fearn_stat
tiemp_rao_stat

effective_sample_size(λ_1_metro_stat_5 )
effective_sample_size(λ_2_metro_stat_5 )
effective_sample_size(λ_1_rao_stat_5 )
effective_sample_size(λ_2_rao_stat_5 )
effective_sample_size(λ_1_fearn_stat_5 )
effective_sample_size(λ_2_fearn_stat_5 )

effective_sample_size(λ_1_metro_stat_5 )/tiemp_metro_stat
effective_sample_size(λ_2_metro_stat_5 )/tiemp_metro_stat
effective_sample_size(λ_1_rao_stat_5 )/tiemp_rao_stat
effective_sample_size(λ_2_rao_stat_5 )/tiemp_rao_stat
effective_sample_size(λ_1_fearn_stat_5 )/tiemp_fearn_stat
effective_sample_size(λ_2_fearn_stat_5 )/tiemp_fearn_stat

################################################################################
########### Initialized in stationarity analysis ###############################
################################################################################

function posterior( λ₁::Float64, λ₂::Float64, Δ::Float64, n_eq_1::Int64, n_eq_2::Int64, n_ch_1::Int64, n_ch_2::Int64, α::Float64, β::Float64)
    return ( λ₂/(λ₁+λ₂)  + λ₁*exp(-(λ₁+λ₂)*Δ)/(λ₁+λ₂) )^n_eq_1 * ( λ₁/(λ₁+λ₂)  + λ₂*exp(-(λ₁+λ₂)*Δ)/(λ₁+λ₂) )^n_eq_2 *
           ( λ₁*( 1.0 - exp(-(λ₁+λ₂)*Δ) )/(λ₁+λ₂) )^n_ch_1 * ( λ₂*( 1.0 - exp(-(λ₁+λ₂)*Δ) )/(λ₁+λ₂) )^n_ch_2 *
           pdf( Gamma(α,1.0/β), λ₁) * pdf( Gamma(α,1.0/β), λ₂)
end


function logposterior( λ₁::Float64, λ₂::Float64, Δ::Float64, n_eq_1::Int64, n_eq_2::Int64, n_ch_1::Int64, n_ch_2::Int64, α::Float64, β::Float64)
    return n_eq_1 * log( λ₂/(λ₁+λ₂)  + λ₁*exp(-(λ₁+λ₂)*Δ)/(λ₁+λ₂) ) + n_eq_2 * log( λ₁/(λ₁+λ₂)  + λ₂*exp(-(λ₁+λ₂)*Δ)/(λ₁+λ₂) ) +
           n_ch_1 * log( λ₁*( 1.0 - exp(-(λ₁+λ₂)*Δ) )/(λ₁+λ₂) ) + n_ch_2 * log( λ₂*( 1.0 - exp(-(λ₁+λ₂)*Δ) )/(λ₁+λ₂) ) +
           logpdf( Gamma(α,1.0/β), λ₁) + logpdf( Gamma(α,1.0/β), λ₂)
end


using HCubature

f(x) = exp( logposterior( x[1], x[2], Δ, n_eq_1, n_eq_2, n_ch_1, n_ch_2, 1.0, 1.0) - logposterior( 2.0, 1.0, Δ, n_eq_1, n_eq_2, n_ch_1, n_ch_2, 1.0, 1.0)  )

norm_c = hcubature(f, [0.0,0.0], [10.0,10.0] )[1]

f_norm(x) = f(x)/norm_c

hcubature(f_norm, [0.0,0.0], [10.0,10.0] )[1]

marg1 = [ hquadrature( i -> f_norm([x,i]), 0.0, 10.0 )[1] for x in collect(1.0:0.01:4.5) ]
marg2 = [ hquadrature( i -> f_norm([i,y]), 0.0, 10.0 )[1] for y in collect(0.0:0.01:2.0) ]

using Plots
#pyplot()

################################################################################
########### Stationarity in plot ##### #########################################
################################################################################

p1_r_stat_5 = histogram( λ_1_rao_stat_5, bins=100, xlim=(1.0,4.5), color=:lightblue, normalize=true, grid=false, label="Gibbs sampler 2 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(1.0:0.01:4.5),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("f_{\\lambda_1 \\, | \\, X }")
title!("\\lambda_1 posterior marginal distribution")
savefig("lambda1histraostat5.png")

p1_m_stat_5 = histogram( λ_1_metro_stat_5, bins=100, xlim=(1.0,4.5), color=:lightblue, normalize=true, grid=false, label="M-H sampler histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(1.0:0.01:4.5),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\\lambda_1 posterior marginal distribution")
savefig("lambda1histmetrostat5.png")

p1_f_stat_5 = histogram( λ_1_fearn_stat_5, bins=100, xlim=(1.0,4.5), color=:lightblue, normalize=true, grid=false, label="Gibbs sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(1.0:0.01:4.5),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\\lambda_1 posterior marginal distribution")
savefig("lambda1histfearnstat5.png")





















################################################################################



f_marg1(x) = hquadrature( i -> f_norm([x,i]), 0.0, 10.0 )[1]

findmax(marg1)

sum( marg1 .<= 2.3 .*pdf.(Normal(collect(0.0:0.01:7.0)[argmax(marg1)],2.0),  collect(0.0:0.01:7.0) ) )

plot(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9)
plot!( collect(0.0:0.01:7.0), 2.3 .*pdf.(Normal(collect(0.0:0.01:7.0)[argmax(marg1)],2.0),  collect(0.0:0.01:7.0) ) )

z₁ = rand(Normal(collect(0.0:0.01:7.0)[argmax(marg1)],2.0))
f_marg1(z₁) / ( 2.3 .*pdf.(Normal(collect(0.0:0.01:7.0)[argmax(marg1)],2.0),z₁) ) > rand(Uniform())

f_cond(y) = f_norm([z₁,y])/f_marg1(z₁)

plot( collect(0.0:0.01:4.0), f_cond.(collect(0.0:0.01:4.0)) )
plot!( collect(0.0:0.01:4.0), 3.25 .* pdf.(Normal(collect(0.0:0.01:4.0)[ findmax( f_cond.(collect(0.0:0.01:4.0)) )[2] ],0.2),collect(0.0:0.01:4.0)) )

collect(0.0:0.01:4.0)[ findmax( f_cond.(collect(0.0:0.01:4.0)) )[2] ]

length(collect(0.0:0.01:4.0))
sum( f_cond.(collect(0.0:0.01:4.0)) .<=  3.25 .* pdf.(Normal( collect(0.0:0.01:4.0)[ findmax( f_cond.(collect(0.0:0.01:4.0)) )[2] ],0.2),collect(0.0:0.01:4.0)) )

z₂ = rand(Normal(1.1,0.2))
f_cond(z₂) / ( 1.8 .*pdf.(Normal(1.1,0.2),z₁) ) > rand(Uniform())

(z₁,z₂)

D_stat = [ z₁, z₂ ]

tiemp_rao_50_stat  = @elapsed  CHAIN_rao_50_stat = raotehHMMGibbs( 2, 50000, T_thin, 1926.0, X_thin, [ [1.0,0.0], [0.0,1.0] ], Categorical, [1.0,0.0], 1.0, 0.6666, ones(2,1),
    (D_stat,P₀))

#@save "CHAIN_rao_50_stat.txt" CHAIN_rao_50_stat
#@save "tiemp_rao_50_stat.txt" tiemp_rao_50_stat
@load "tiemp_rao_50_stat.txt"
(tiemp_rao_50_stat/60)/60


################################################################################
########### Make plots section ##### ###########################################
################################################################################


using Plots
#pyplot()

p1_r_5 = histogram( λ_1_raoteh_5, bins=150, xlim=(1,4), color=:pink, normalize=true, grid=false, label="Gibbs sampler 2 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(1.0:0.01:4.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\\lambda_1 posterior marginal distribution")

p1_m_5 = histogram( λ_1_metro_5[5000:end], color=:pink, normalize=true, grid=false, label="M-H sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(1.0:0.01:4.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\\lambda_1 posterior marginal distribution")


p1_r_50_2 = histogram( λ_1_raoteh_50_2[10000:end], color=:pink, normalize=true, grid=false, label="M-H sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

p1_g = histogram( chains_vec_gibbs_1[4], bins=150, xlim=(0,7), color=:lightblue, normalize=true, grid=false, label="Gibbs sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

p1_m = histogram( chains_vec_metro_1[4], bins=150, xlim=(0,7), color=:lightblue, normalize=true, grid=false, label="M-H sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

p1_r = histogram( λ_1_raoteh_10, bins=150, xlim=(0,7), color=:lightblue, normalize=true, grid=false, label="M-H sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

# Stationary Plots

p1_m_s = histogram!( λ_1_metro_10_stat, bins=150, xlim=(0,7), color=:pink, normalize=true, grid=false, label="M-H sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

p1_r_s = histogram( λ_1_raoteh_10_stat, bins=150, xlim=(0,7), color=:lightblue, normalize=true, grid=false, label="M-H sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

p1_g_s = histogram( λ_1_gibbs_10_stat, bins=150, xlim=(0,7), color=:lightblue, normalize=true, grid=false, label="Gibbs sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

p2=histogram(λ_2_gibbs_10_2, bins=150, xlim=(0,4), color=:lightblue, normalize=true, grid=false, label="Gibbs sampler 1 histogram", dpi=300)
plot!(collect(0.0:0.01:4.0),marg2,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, width=9cm,height=10cm]X } \$")
title!("\$ \\lambda_2 \$ posterior marginal distribution")
plot(p1,p2,size=(800,500),dpi=300)


p1_m_prev = histogram( λ_1_metro_100_prev, bins=150, xlim=(0,7), color=:lightblue, normalize=true, grid=false, label="M-H sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

################################################################################

effective_sample_size(λ_1_metro_10_stat)

effective_sample_size(λ_1_metro_10)


p1_m_50 = histogram!( λ_1_metro_50, bins=150, xlim=(0,7), color=:pink, normalize=true, grid=false, label="M-H sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

p1_r_50_1 = histogram( λ_1_raoteh_50_1[10000:end], bins=150, xlim=(0,7), color=:pink, normalize=true, grid=false, label="M-H sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

p1_r_50_2 = histogram( λ_1_raoteh_50_2[10000:end], bins=150, xlim=(0,7), color=:pink, normalize=true, grid=false, label="M-H sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

p1_r_50_3 = histogram( λ_1_raoteh_50_3[10000:end], bins=150, xlim=(0,7), color=:pink, normalize=true, grid=false, label="M-H sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

p1_r_50_4 = histogram( λ_1_raoteh_50_4[10000:end], bins=150, xlim=(0,7), color=:pink, normalize=true, grid=false, label="M-H sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

p1_r_50_5 = histogram( λ_1_raoteh_50_5[10000:end], bins=150, xlim=(0,7), color=:pink, normalize=true, grid=false, label="M-H sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

p1_r_150 = histogram([λ_1_raoteh_50_1[25000:end]; λ_1_raoteh_50_2[25000:end]; λ_1_raoteh_50_3[25000:end]], bins=150, xlim=(0,7), color=:pink, normalize=true, grid=false, label="M-H sampler 1 histogram", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
plot!(collect(0.0:0.01:7.0),marg1,linewidth=2.5, color=:red, label="True posterior marginal", legendfontsize=12, xtickfontsize=9, ytickfontsize=9, dpi=300)
#ylabel!("\$ f_{\\lambda_1 \\, | \\, X } \$")
title!("\$ \\lambda_1 \$ posterior marginal distribution")

loglikelihood_given_hidden_path( T, T_thin, X_thin, [ [1.0,0.0], [0.0,1.0] ], Categorical)


###############################################################################

using JLD2

@load "CHAIN_rao_5.txt"

s_1000 = [ CHAIN_rao_5[4][i][findlast( CHAIN_rao_5[3][i] .<= 1000.5 )] for i in 1:length(CHAIN_rao_5[3]) ]

using StatsBase

mean(s_1000)


prop_correct_s_estim = [ mean( X[findlast( T .<= j )]  .== [ CHAIN_rao_5[4][i][findlast( CHAIN_rao_5[3][i] .<= j )] for i in 1:length(CHAIN_rao_5[3]) ]) for j in collect( 0.5:1.0:(T_thin[end]-0.5) )]


mid_s = [ X[findlast( T .<= j )] for j in collect( 0.5:1.0:(T_thin[end]-0.5) ) ]

using Plots

scatter( collect( 0.5:1.0:(T_thin[end]-0.5) ), prop_correct_s_estim )
#scatter!( collect( 0.5:1.0:(T_thin[end]-0.5) ), mid_s )

prop_correct_s_estim
