using LogDensityProblems, Statistics, LogDensityProblemsAD

struct NormalPosterior{T} # contains the summary statistics
    N::Int
    x̄::T
    S::T
end

# calculate summary statistics from a data vector
function NormalPosterior(x::AbstractVector)
    NormalPosterior(length(x), mean(x), var(x; corrected = false))
end

# define a callable that unpacks parameters, and evaluates the log likelihood
function (problem::NormalPosterior)(θ)
    @unpack μ, σ = θ
    @unpack N, x̄, S = problem
    loglikelihood = -N * (log(σ) + (S + abs2(μ - x̄)) / (2 * abs2(σ)))
    logprior = - abs2(σ)/8 - abs2(μ)/50
    loglikelihood + logprior
end

problem = NormalPosterior(randn(100))