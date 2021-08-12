
#using LinearAlgebra, Distributions, ExpmV, SparseArrays, Distributed, ProgressMeter
using LinearAlgebra, Distributions, SciPy, SparseArrays, Distributed, ProgressMeter
using SpecialFunctions: factorial
using Statistics: mean, var
using StatsBase: mean_and_var
using StatsFuns: logistic
expm = SciPy.sparse.linalg.expm
expm_multiply = SciPy.sparse.linalg.expm_multiply

# Functions for explicit posteriors in Sim Study 1

function Num_visits_from_j_to_k( j::Int64, k::Int64, T::Array{Float64,1}, X::Array{Int64,1})
    return length( findall( i -> ( X[i-1] == j ) && ( X[i] == k ), collect(2:length(X)) ) )
end

function Num_visits_from_j( num_states::Int64, j::Int64, T::Array{Float64,1}, X::Array{Int64,1})
    return mapreduce( k -> Num_visits_from_j_to_k( j, k, T, X), +, collect(1:num_states)[1:end .!= j] )
end

function Time_spent_in_j( j::Int64, T::Array{Float64,1}, X::Array{Int64,1})
    in_j_vec = findall( i -> X[i] == j, collect(1:(length(X)-1) ) )
    return mapreduce( k -> T[k+1] - T[k], +, in_j_vec, init = 0.0)
end

# Functions for D and P initialization

function D_inic(n_states::Int64,hyp_α::Float64=1.0,hyp_β::Float64=1.0)
    D = zeros(Float64,n_states)
    for i in 1:n_states
        D[i] = rand( Gamma(hyp_α,1.0/hyp_β) )
    end
    return D
end

function P_inic(n_states::Int64,hyp_α::Float64=1.0)
    P = zeros(Float64,n_states,n_states-1)
    for i in 1:n_states
        P[i,:] = rand(Dirichlet(n_states-1,hyp_α))
    end
    return P
end

# Proposal for initial value in Markov Chain
function prop_inic( λ::Array{Float64} )
    k = length(λ)
    D_0 = rand.( Exponential.(1.0./λ) )
    P_0 = rand(  Dirichlet( k-1, 1.0), k)'
    P_0 = copy(P_0)
    return D_0, P_0
end


# Function to obtain generator G from D and P

@inline function GeneratorParam( D::Array{Float64,1}, P::Array{Float64,2})
    k = size(D)[1]
    G = zeros( Float64, k, k)
    for i in 1:k
        G[i,i] = -D[i]
        for j in 1:(i-1)
            G[i,j] = P[i,j] * D[i]
        end
        for j in (i+1):k
            G[i,j] = P[i,j-1] * D[i]
        end
    end
    return G
end

################################################################################
###### Metropolis ##############################################################
################################################################################

# Loglikelihood function with usual squaring and scaling exp mat
function parallel_loglik( D::Array{Float64,1}, P::Array{Float64,2},  Δ::Float64, S::Array{Int64,1} )
    G = GeneratorParam( D, P)
    M = exp(Δ*G)
    return @distributed (+) for i in 1:(length(S)-1)
        log( M[ S[i], S[i+1] ] )
    end
end

# Loglikelihood function with Al-Mohy exp mat
# function parallel_loglik_AlMohy( D::Array{Float64,1}, P::Array{Float64,2}, Δ::Float64, S::Array{Int64,1} )
#     G = sparse( GeneratorParam( D, P) )
#     return @distributed (+) for i in 1:(length(S)-1)
#         log(  expmv( Δ, G, Matrix(I,size(G))[:,S[i+1]] )[S[i]] )
#     end
# end

# Loglikelihood function with Al-Mohy of SciPy exp mat
function parallel_loglik_AlMohy( D::Array{Float64,1}, P::Array{Float64,2}, Δ::Float64, S::Array{Int64,1} )
    G = GeneratorParam( D, P)
    return @distributed (+) for i in 1:(length(S)-1)
        log( expm_multiply( Δ*G, Matrix(I,size(G))[:,S[i+1]] )[S[i]]  )
    end
end

# Logarith of prior function with cancelations considered for a step which changes the n-th diagonal elemnt of G
function logprior_diag( n::Int64, D::Array{Float64,1}, α::Float64, β::Float64)
    return logpdf( Gamma(α,1.0/β), D[n] )
end

# Logarith of prior function with cancelations considered for a step which changes the n-th row of G without the diagonal element
function logprior_no_diag( n::Int64, P::Array{Float64,2}, c::Float64, α::Array{Float64,1})
    return logpdf( Dirichlet( c .* α ), P[n,:] )
end

# Proposal for diagonal elements in the n-th row of G given by a
# LogNormal centered in log(D[n]) and with scale var
function prop_diag( n::Int64, D::Array{Float64,1}, var::Float64)
    D_prop = copy(D)
    D_prop[n] = rand( LogNormal( log( D[n] ), var ) )
    return D_prop
end

# Proposal for off-diagonal elements in the n-th row of G given by a
# Dirichlet distribution uniform on the (k-1)-simplex, where G is a k×k
# matrix.
function prop_row_mode( n::Int64, P::Array{Float64,2}, c::Float64)
    P_prop = copy(P)
    P_prop[n,:] = rand( Dirichlet( 1.0 .+ c.*P[n,:] ) )
    return P_prop
end

# Metropolis-Hastings proposal kernel without last row of P as it corresponds to absorbing state
function log_q( m::Int64, P_new::Array{Float64,2}, P_old::Array{Float64,2}, c::Float64)
    return logpdf( Dirichlet( 1.0 .+ c .* P_old[m,:] ), P_new[m,:] )
end

# Metropolis-Hastings algorithm at log scale
@inline function metr_hast_chain_D_P_mode( num_states::Int64, S::Array{Int64,1}, Δ::Float64, l::Int64,
    prop_lognorm_var::Float64, prop_dir_prec::Float64, prior_α::Float64, prior_β::Float64, prior_dir_prec::Float64,
    Inic_set::Tuple{Array{Float64,1},Array{Float64,2}})
    D_chain = [ zeros(Float64, num_states) for i in 1:l ]
    P_chain = [ zeros(Float64, num_states,num_states-1) for i in 1:l]
    D = Inic_set[1]
    P = Inic_set[2]
    log_π_lik = parallel_loglik( D, P, Δ, S)
    prog = Progress( l, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50)
    for i in 1:l
        for m in 1:num_states
            log_π_prior = logprior_diag( m, D, prior_α, prior_β)
            D_prop = prop_diag( m, D, prop_lognorm_var)
            log_π_lik_prop =  parallel_loglik( D_prop, P, Δ, S)
            log_π_prior_prop = logprior_diag( m, D_prop, prior_α, prior_β)
            log_q_given_D = - log( D_prop[m] ) #  Because of exponential transformation Jacobian in LogNormal proposal
            log_q_given_D_prop = - log( D[m] ) #  Because of exponential transformation Jacobian in LogNormal proposal
            u = rand(Uniform())
            if log(u) < min( log_π_lik_prop + log_π_prior_prop - log_π_lik - log_π_prior + log_q_given_D_prop - log_q_given_D, 0.)
                D = D_prop
                log_π_lik = log_π_lik_prop
            end
        end
        for m in 1:num_states
            log_π_prior = logprior_no_diag( m, P, prior_dir_prec, repeat([1.0],num_states-1))
            P_prop = prop_row_mode( m, P, prop_dir_prec )
            log_π_lik_prop = parallel_loglik( D, P_prop, Δ, S)
            log_π_prior_prop = logprior_no_diag( m, P_prop, prior_dir_prec, repeat([1.0],num_states-1))
            log_q_given_P = log_q( m, P_prop, P, prop_dir_prec)
            log_q_given_P_prop = log_q( m, P, P_prop, prop_dir_prec )
            u = rand(Uniform())
            if log(u) < min( log_π_lik_prop + log_π_prior_prop - log_π_lik - log_π_prior + log_q_given_P_prop - log_q_given_P, 0.)
                P = P_prop
                log_π_lik = log_π_lik_prop
            end
        end
        D_chain[i] = D
        P_chain[i] = P
        next!(prog)
    end
    return [ D_chain, P_chain ]
end

# Metropolis-Hastings algorithm at log scale with Al-Mohy expoenntial matrix
@inline function metr_hast_chain_D_P_mode_AlMohy( num_states::Int64, S::Array{Int64,1}, Δ::Float64, l::Int64,
    prop_lognorm_var::Float64, prop_dir_prec::Float64, prior_α::Float64, prior_β::Float64, prior_dir_prec::Float64,
    Inic_set::Tuple{Array{Float64,1},Array{Float64,2}})
    D_chain = [ zeros(Float64, num_states) for i in 1:l ]
    P_chain = [ zeros(Float64, num_states,num_states-1) for i in 1:l]
    D = Inic_set[1]
    P = Inic_set[2]
    log_π_lik = parallel_loglik( D, P, Δ, S)
    prog = Progress( l, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50)
    for i in 1:l
        for m in 1:num_states
            log_π_prior = logprior_diag( m, D, prior_α, prior_β)
            D_prop = prop_diag( m, D, prop_lognorm_var)
            log_π_lik_prop =  parallel_loglik_AlMohy( D_prop, P, Δ, S)
            log_π_prior_prop = logprior_diag( m, D_prop, prior_α, prior_β)
            log_q_given_D = - log( D_prop[m] ) #  Because of exponential transformation Jacobian in LogNormal proposal
            log_q_given_D_prop = - log( D[m] ) #  Because of exponential transformation Jacobian in LogNormal proposal
            u = rand(Uniform())
            if log(u) < min( log_π_lik_prop + log_π_prior_prop - log_π_lik - log_π_prior + log_q_given_D_prop - log_q_given_D, 0.)
                D = D_prop
                log_π_lik = log_π_lik_prop
            end
        end
        for m in 1:num_states
            log_π_prior = logprior_no_diag( m, P, prior_dir_prec, repeat([1.0],num_states-1))
            P_prop = prop_row_mode( m, P, prop_dir_prec )
            log_π_lik_prop = parallel_loglik_AlMohy( D, P_prop, Δ, S)
            log_π_prior_prop = logprior_no_diag( m, P_prop, prior_dir_prec, repeat([1.0],num_states-1))
            log_q_given_P = log_q( m, P_prop, P, prop_dir_prec)
            log_q_given_P_prop = log_q( m, P, P_prop, prop_dir_prec )
            u = rand(Uniform())
            if log(u) < min( log_π_lik_prop + log_π_prior_prop - log_π_lik - log_π_prior + log_q_given_P_prop - log_q_given_P, 0.)
                P = P_prop
                log_π_lik = log_π_lik_prop
            end
        end
        D_chain[i] = D
        P_chain[i] = P
        next!(prog)
    end
    return [ D_chain, P_chain ]
end

################################################################################
###### Fearnhead ###############################################################
################################################################################

# Sampler for number of dominating events
function number_dom_sim( t::Float64, G::Array{Float64,2}, Sₑ::Int64, S₀::Int64)
    U = rand(Uniform())
    r_run = 0
    ρ = maximum( -diag(G) )
    B = exp(t*G)
    M = (1.0/ρ)*G + Matrix{Float64}(I, size(G))
    M_vec = [ Matrix{Float64}(I, size(G)) ]
    M_run = I(size(G)[1])
    ρt_run = 1.0
    r_fact = 1.0
    sum_run = exp(-ρ*t) * ρt_run * (M_run)[S₀,Sₑ] / ( B[S₀,Sₑ] * r_fact )
    while U > sum_run
        M_run *= M
        push!( M_vec, M_run)
        r_run += 1
        ρt_run *= ρ*t
        r_fact *= r_run
        sum_run += exp(-ρ*t) * ρt_run * (M_run)[S₀,Sₑ] / ( B[S₀,Sₑ] * r_fact )
    end
    return r_run, M_vec
end

# Sampler for CTMC bridge between times T₀ and Tₑ with corresponding states S₀
# and Sₑ in the bridge boundaries and G the generator matrix
function sim_bridge( G::Array{Float64,2}, Sₑ::Int64, S₀::Int64, Tₑ::Float64, T₀::Float64=0.0  )
    n, M_vec = number_dom_sim( Tₑ-T₀, G, Sₑ, S₀)
    T = sort( rand( Uniform(T₀,Tₑ), n) )
    d = size(G)[1] # number of hidden states
    M = (1.0/maximum( -diag(G) ))*G + Matrix{Float64}(I, size(G))
    S = zeros( Int64, n+2) # array for hidden Markov porcess states, starting with S₀ and ending with Sₑ
    S[1] = S₀ # initial hiiden state value
    S[n+2] = Sₑ # ending hidden state value
    for i in 2:(n+1)
        S[i] = rand( Categorical( [ M[S[i-1],j] * M_vec[n-i+2][j,Sₑ] /  M_vec[n-i+3][S[i-1],Sₑ] for j in 1:d ] ) )
    end
    T = [ T₀; T; Tₑ ]
    return T, S
end

# Full sample of CTMC given discrete observations X at times T and G the generator matrix
function data_sim_bridge( G::Array{Float64,2}, T::Array{Float64,1}, X::Array{Int64,1})
    m = length(T)
    X_full = Int64[]
    T_full = Float64[]
    for i in 2:m
        T_loc, X_loc = sim_bridge( G, X[i], X[i-1], T[i], T[i-1]  )
        append!( T_full, T_loc[1:(end-1)] )
        append!( X_full, X_loc[1:(end-1)] )
    end
    return [ T_full; T[end] ], [ X_full; X[end] ]
end

# Sample of generator matrix G given a complete observed path given by jumps to states X at times T
function G_given_full_path(  num_states::Int64, prior_α::Float64, prior_β::Float64, prior_dir_α::Array{Float64,2}, T::Array{Float64,1}, X::Array{Int64,1})
    D_new = [ rand( Gamma( prior_α + Num_visits_from_j( num_states, i, T, X), 1.0/( prior_β + Time_spent_in_j( i, T, X) ) ) ) for i in 1:num_states ]
    P_new = zeros(Float64, num_states, num_states -1 )
    for i in 1:num_states
        P_new[i,:] = rand( Dirichlet( prior_dir_α[i,:] + [ Num_visits_from_j_to_k(i,k,T,X) for k in  collect(1:num_states)[1:end .!= i] ]  ) )
    end
    return GeneratorParam( D_new, P_new)
end

# Gibbs sampler using Fearnhead and Sherlock's bridge algorith for simulating missing
# jumps for a discretely observed CTMC
function fearnhead_gibbs( num_states::Int64, X::Array{Int64,1}, T::Array{Float64,1}, l::Int64,
                    prior_α::Float64, prior_β::Float64, prior_dir_α::Array{Float64,2},
                    Inic_set::Tuple{Array{Float64,1},Array{Float64,2}})
    G_chain = [ zeros(Float64, num_states, num_states) for i in 1:l ]
    G_run = GeneratorParam( Inic_set[1], Inic_set[2])
    G_chain[1] = G_run
    prog = Progress( l, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50)
    for i in 1:l
        T_full_run, X_full_run = data_sim_bridge( G_run, T, X)
        G_run = G_given_full_path( num_states, prior_α, prior_β, prior_dir_α, T_full_run, X_full_run)
        G_chain[i] = G_run
        next!(prog)
    end
    return G_chain
end

################################################################################
###### Rao Teh Gibbs ###########################################################
################################################################################

# Sample of generator matrix G given a complete observed path given by jumps to states X at times T
function D_P_given_full_path(  num_states::Int64, prior_α::Float64, prior_β::Float64, prior_dir_α::Array{Float64,2}, T::Array{Float64,1}, X::Array{Int64,1})
    D_new = [ rand( Gamma( prior_α + Num_visits_from_j( num_states, i, T, X), 1.0/( prior_β + Time_spent_in_j( i, T, X) ) ) ) for i in 1:num_states ]
    P_new = zeros(Float64, num_states, num_states -1 )
    for i in 1:num_states
        P_new[i,:] = rand( Dirichlet( prior_dir_α[i,:] + [ Num_visits_from_j_to_k(i,k,T,X) for k in  collect(1:num_states)[1:end .!= i] ]  ) )
    end
    return D_new, P_new
end

function deuniformization(T::Array{Float64,1},S::Array{Int64,1})
    T_deunif = T[1:1]
    S_deunif = S[1:1]
    for i in 2:length(T)
        if S[i] != S_deunif[end]
            append!( T_deunif, T[i] )
            append!( S_deunif, S[i] )
        end
    end
    return T_deunif, S_deunif
end

# Simulation of a piecewise Poisson process which has intensity λ[1] in [t_0,T[1])
# λ[i] in [T[i],T[i+1]) for i in 1:(length(T)-1)
# λ[end] in [T[end],t_e)
# with length(λ)=length(T)+1
function piecewise_poisson_process( λ::Array{Float64,1}, T::Array{Float64,1},
        t_e::Float64)
    U = Float64[]
    T_v = [T;t_e]
    c =  [ 0.0; cumsum( [T_v[j+1] - T_v[j] for j in 1:(length(T)-1)].*λ[1:(end-1)] ) ]
    Γ_run = rand(Exponential())
    i_run = findfirst( Γ_run .<= c )
    if !isnothing(i_run)
        t_run = ( Γ_run - c[i_run] )/λ[i_run] + T[i_run]
        while t_run < t_e
            push!(U, t_run)
            Γ_run += rand(Exponential())
            if  Γ_run > c[i_run]
                if c[end] > Γ_run
                    i_run += findfirst( Γ_run .<= c[(i_run+1):end] ) - 1
                else
                    i_run = length(λ)
                end
            end
            t_run = ( Γ_run - c[i_run] )/λ[i_run] + T[i_run]
        end
    end
    return U
end

# Likelihood given the full hidden CTMC
function loglikelihood_given_hidden_path( T_S::Array{Float64,1}, T_Y::Array{Float64,1}, Y::Array{Int64,1}, n_states::Int64)
    logl = [ zeros(Float64,n_states) for i in 1:length(T_S) ]
    for i in 1:(length(T_S)-1)
        run_loc = ( T_S[i] .<= T_Y .< T_S[i+1] )
        logl[i] = [ mapreduce( x -> j == x ? 0.0 : -Inf, +, Y[run_loc], init=0.0) for j in 1:n_states ]
    end
    run_loc = ( T_S[end] .<= T_Y  )
    logl[end] = [ mapreduce( x -> j == x ? 0.0 : -Inf, +, Y[run_loc], init=0.0) for j in 1:n_states ]
    return logl
end

# Forward Backward algorithm for when only one transition matrix B is used
function ForwBackw_discrete_time( B::Array{Float64,2},
           logl::Array{Array{Float64,1},1}, π₀::Array{Float64,1})
    log_B = log.(B)
    n = length(logl) # Number of observations
    m = length(π₀) # Number of hidden states
    S = zeros( Int64, n)
    # Array of probabilities for forward recurssion
    log_α = zeros(n,m)
    # Initial forward accumulation of probabilities given by prior π₀
    log_α[1,:] = log.(π₀)
    # Loop forward accumulation of probabilities with logsumexp trick
    for t in 2:n
        α_tmp = zeros(m,m)
        b_tmp = zeros(m)
        for j in 1:m
            α_tmp[j,:] = [ log_α[t-1,k] + logl[t-1][k] + log_B[j,k] for k in 1:m ]
            b_tmp[j] = maximum(α_tmp[j,:])
        end
        log_α[t,:] = [ b_tmp[j] + log( sum(  exp.( α_tmp[j,:] .- b_tmp[j] ) ) ) for j in 1:m ]
    end
    # Initial backward simulation of hidden states with Max Gumbel trick
    log_β = logl[n] .+ log_α[n,:]
    S[n] = findmax( rand( Gumbel(), m) .+ log_β )[2]
    # Loop bacward simulation of the hidden states with Max Gumbel trick
    for i in (n-1):-1:1
        log_β = logl[i] .+ log_α[i,:] .+ log_B[S[i+1],:]
        S[i] = findmax( rand( Gumbel(), m) .+ log_β )[2]
    end
    return S
end

# Simulation for hidden trajectory given generator G and previous complete hidden
# trajectory for Rao Teh Gibbs algorithm
function raotehHiddengivenG(  D::Array{Float64,1}, P::Array{Float64,2},
        T_S::Array{Float64,1}, S::Array{Int64,1}, Y::Array{Int64,1},
        T_Y::Array{Float64,1}, Tₑ::Float64, π₀::Array{Float64,1},
        ϵ::Float64=1.0)
        n_states = length(D)
        Ω = (1.0+ϵ) * maximum(D)
        G = GeneratorParam( D, P)
        λ = Ω .- D[S]
        B =  one(G) .+ (G./Ω)
        U = piecewise_poisson_process( λ, T_S, Tₑ)
        completed_T = sort([T_S;U])
        logl = loglikelihood_given_hidden_path( completed_T, T_Y, Y, n_states)
        unif_S = ForwBackw_discrete_time( B, logl, π₀) # uniformized state values so there are jumps to the same state
        new_T_S, new_S = deuniformization( [completed_T; Tₑ], [unif_S; unif_S[end]])
        return new_T_S, new_S
end


function raotehHMMGibbs( n_states::Int64, n_iter::Int64,
        T::Array{Float64}, Tₑ::Float64, Y::Array{Int64}, π₀::Array{Float64,1},
        prior_α::Float64, prior_β::Float64, prior_dir_α::Array{Float64,2},
        Inic_set::Tuple{Array{Float64,1},Array{Float64,2}})
    chain_D = [ zeros(Float64,n_states) for i in 1:n_iter]
    chain_P = [ zeros(Float64,n_states,n_states-1) for i in 1:n_iter]
    chain_T_S = [ Float64[] for i in 1:n_iter ]
    chain_S = [ Int64[] for i in 1:n_iter ]
    chain_D[1] = Inic_set[1]
    chain_P[1] =  Inic_set[2]
    chain_T_S[1], chain_S[1] = deuniformization( data_sim_bridge( GeneratorParam( chain_D[1], chain_P[1]), T, convert( Array{Int64}, Y))... )
    #println("MCMC run has started.")
    prog = Progress( (n_iter-1), dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50)
    for i in 2:n_iter
        chain_D[i], chain_P[i] = D_P_given_full_path( n_states, prior_α, prior_β, prior_dir_α, chain_T_S[i-1], chain_S[i-1])
        chain_T_S[i], chain_S[i] = raotehHiddengivenG( chain_D[i], chain_P[i], chain_T_S[i-1], chain_S[i-1], Y, T, Tₑ, π₀)
        next!(prog)
    end
    return chain_D, chain_P, chain_T_S, chain_S
end

################################################################################
###### Extra functions #########################################################
################################################################################

function data_thin( δ::Float64, T::Array{Float64,1}, X::Array{Int64,1})
    m = Int64( fld(T[end], δ) )
    T_thin = zeros(Float64, m)
    X_thin = zeros(Int64, m )
    for i in 1:m
        index_loc = findlast( T .<= δ*i )
        T_thin[i] = δ*i
        X_thin[i] = X[index_loc]
    end
    return [ T[1]; T_thin], [X[1];X_thin]
end

function data_gen( D::Array{Float64,1}, P::Array{Float64,2}, m::Int64, t_0::Float64=0.0)
    num_states = size(D)[1]
    p_0 = zeros(Float64, num_states-1) .+ 1.0/(num_states-1)
    m_count = 0
    t_run = 0.0
    T = Float64[ t_0 ]
    X = Int64[]
    x_run = Int64( rand( Categorical(p_0) ) )
    push!( X, x_run )
    while m_count < m
        t_run += rand( Exponential( 1.0/D[x_run] ) )
        push!(T,t_run)
        x_run = collect(1:num_states)[1:end .!= x_run ][ rand( Categorical( P[x_run,:] ) ) ]
        push!(X, x_run )
        m_count += 1
    end
    return T, X
end

################################################################################
###### MCMC diagnostic functions ###############################################
################################################################################

function autocorrelation(x::Array{Float64,1}, k::Integer, v::Float64 = var(x))
    x1 = x[1:(end-k)]
    x2 = x[(1+k):end]
    V = sum((x1 .- x2).^2) / length(x1)
    return 1 - V / (2*v)
end

function ess_factor_estimate(x::Array{Float64,1}, v::Float64 = var(x))
    N = length(x)
    τ_inv = 1 + 2 * autocorrelation(x, 1, v)
    K = 2
    while K < N - 2
        Δ = autocorrelation(x, K, v) + autocorrelation(x, K + 1, v)
        if Δ < 0
            break
        else
            τ_inv += 2*Δ
            K += 2
        end
    end
    return 1 / τ_inv
end

function effective_sample_size(x::Array{Float64,1}, v::Float64 = var(x))
    τ = ess_factor_estimate(x, v)
    return τ * length(x)
end

function potential_scale_reduction(chains::Array{Array{Float64,1},1})
    N = length(chains[1])
    mvs = mean_and_var.(chains)
    W = mean(last.(mvs))
    B = var(first.(mvs))
    return √( (N-1)/N + B/W )
end

[ D_inic(2) for _ in 1:4 ]

function stationary_diagnost_metro_modeprop(C::Int64,n_states::Int64,k::Int64,
    prop_lognorm_var::Float64,prop_dir_prec::Float64,prior_α::Float64,prior_β::Float64,prior_dir_prec::Float64,D₀vec::Array{Array{Float64,1},1}=[[]])
    if D₀vec == [[]]
        D₀run = [ D_inic(n_states) for _ in 1:C ]
    else
        D₀run = D₀vec
    end
    D₀run = [ D_inic(n_states) for _ in 1:C ]
    P₀ = hcat([1.0,1.0])
    chains = [ [ [ D₀run[j][i] ] for j in 1:C  ] for i in 1:n_states ]
    tiemp = zeros(C)
    psr = [ Inf for _ in n_states ]
    while ! all( psr .< 1.01 )
        for i in 1:C
            #println(D₀run[i])
            tiemp[i] += @elapsed  chain = metr_hast_chain_D_P_mode( n_states, X_thin, Δ, k, prop_lognorm_var, prop_dir_prec, prior_α, prior_β, prior_dir_prec, (D₀run[i],P₀))
            D₀run[i] = chain[1][end]
            for j in 1:n_states
                append!( chains[j][i], [ chain[1][h][j] for h in 2:length(chain[1]) ] )
            end
        end
        psr = [ potential_scale_reduction(chains[i]) for i in 1:n_states ]
    end
    return chains, tiemp
end

function stationary_diagnost_raoteh(C::Int64,n_states::Int64,k::Int64,
    prior_α::Float64,prior_β::Float64,prior_dir_prec::Float64,D₀vec::Array{Array{Float64,1},1}=[[]])
    if D₀vec == [[]]
        D₀run = [ D_inic(n_states) for _ in 1:C ]
    else
        D₀run = D₀vec
    end
    P₀ = hcat([1.0,1.0])
    chains = [ [ [ D₀run[j][i] ] for j in 1:C  ] for i in 1:n_states ]
    tiemp = zeros(C)
    psr = [ Inf for _ in n_states ]
    while ! all( psr .< 1.01 )
        for i in 1:C
            tiemp[i] += @elapsed chain = raotehHMMGibbs( n_states, k, T_thin, T_thin[end], X_thin, [1.0,0.0], prior_α, prior_β, prior_dir_prec.*ones(2,1), (D₀run[i],P₀))
            D₀run[i] = chain[1][end]
            for j in 1:n_states
                append!( chains[j][i], [ chain[1][h][j] for h in 2:length(chain[1]) ] )
            end
        end
        psr = [ potential_scale_reduction(chains[i]) for i in 1:n_states ]
        D₀run = D₀vec
    end
    return chains, tiemp
end

function stationary_diagnost_fearnhead(C::Int64,n_states::Int64,k::Int64,
    prior_α::Float64,prior_β::Float64,prior_dir_prec::Float64,D₀vec::Array{Array{Float64,1},1}=[[]])
    if D₀vec == [[]]
        D₀run = [ D_inic(n_states) for _ in 1:C ]
    else
        D₀run = D₀vec
    end
    P₀ = hcat([1.0,1.0])
    chains = [ [ [ D₀run[j][i] ] for j in 1:C  ] for i in 1:n_states ]
    tiemp = zeros(C)
    psr = [ Inf for _ in n_states ]
    while ! all( psr .< 1.01 )
        for i in 1:C
            # println( D₀run[i] )
            tiemp[i] += @elapsed chain = fearnhead_gibbs( n_states, X_thin, T_thin, k, prior_α, prior_β, prior_dir_prec.*ones(2,1), (D₀run[i],P₀))
            D₀run[i] = -diag( chain[end] )
            for j in 1:n_states
                append!( chains[j][i], [ -chain[h][j,j] for h in 2:length(chain) ] )
            end
        end
        psr = [ potential_scale_reduction(chains[i]) for i in 1:n_states ]
        D₀run = D₀vec
    end
    return chains, tiemp
end

################################################################################
########### Initialization functions ###########################################
################################################################################

function D_inic(n_states::Int64)
    D = zeros(Float64,n_states)
    for i in 1:n_states
        u = rand(Uniform(-2,2))
        D[i] = exp(u)
    end
    return D
end

function P_inic(n_states::Int64)
    P = zeros(Float64,n_states,n_states-1)
    for i in 1:n_states
        u = rand( Uniform(-2,2), n_states-1)
        z = logistic.( u .+ log.( 1.0./( (n_states-1) .- collect(1:n_states-1) ) ) )
        for j in 1:(n_states-2)
            P[i,j] = ( 1 - sum( P[i,1:(j-1)] ) )*z[j]
        end
        P[i,end] = 1 - sum( P[i,1:(end-1)] )
    end
    return P
end

################################################################################
###### Metropolis with Forward Backward ########################################
################################################################################

# Likelihood given the hidden CTMC at the discrete observation times
function likelihood_given_hidden( Y::Array{Int64,1}, θ::Array{Array{Float64,1},1}, F::UnionAll)
    l = [ zeros( Float64, length(θ)) for i in 1:length(Y) ]
    for i in 1:length(Y)
        l[i] = [ pdf( F( θ[j]... ), Y[i]) for j in 1:length(θ) ]
    end
    return l
end

# Simulation for hidden trajectory given generator G and emitted data Y at time
# intervals Δ
function hiddengivenG(  D::Array{Float64,1}, P::Array{Float64,2},
        Δ::Float64, Y::Array{Float64,1},
        θ::Array{Array{Float64,1},1}, F::UnionAll, π₀::Array{Float64,1})
        G = GeneratorParam( D, P)
        B = exp(Δ*G)
        l = likelihood_given_hidden( Y, θ, F) # Likelihood vector for Forw-Back recurssion
        S = ForwBackw_discrete_time( B, l, π₀) # Forw-Backw recurssion for simulating hiddem values at obsrved times for Y
        return  S
end

# Random walk uniform (0,1) Metropolis-Hastings algorithm at log scale
@inline function metr_hast_chain_D_P_hidden( n_iter::Int64, n_MH_D_P::Int64, n_states::Int64, Y::Array{Float64,1}, Δ::Float64,
                    θ::Array{Array{Float64,1},1}, F::UnionAll, π₀::Array{Float64,1},
                    prop_dir_prec::Float64, prop_lognorm_var::Float64, Gamma_hyp_α::Float64, Gamma_hyp_β::Float64, Dir_hyp_α::Float64,
                    D_inic_α::Float64, D_inic_β::Float64, P_inic_α::Float64)
    chain_D = [ zeros(Float64, n_states) for i in 1:n_iter ]
    chain_P = [ zeros(Float64, n_states, n_states-1) for i in 1:n_iter ]
    chain_S = [ zeros(Int64, length(Y)) for i in 1:n_iter ]
    chain_D[1] = D_inic(n_states,D_inic_α,D_inic_β)
    chain_P[1] = P_inic(n_states,P_inic_α)
    chain_S[1] = rand( Categorical(n_states), length(Y) )
    println("MCMC run has started.")
    prog = Progress( (n_iter-1), dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50)
    for i in 2:n_iter
        chain_D[i], chain_P[i] = metr_hast_chain_D_P_hidden( n_MH_D_P, chain_D[i-1], chain_P[i-1], chain_S[i-1], Δ, prop_lognorm_var, prop_dir_prec, Gamma_hyp_α, Gamma_hyp_β, Dir_hyp_α)
        chain_S[i] = hiddengivenG( chain_D[i], chain_P[i], Δ, Y, θ, F, π₀)
        next!(prog)
    end
    return chain_D, chain_P, chain_S
end
