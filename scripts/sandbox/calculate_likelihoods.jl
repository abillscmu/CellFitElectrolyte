function chain_to_nt(chain,idx)
    return (ω = chain[:ω].data[idx,1], εₑ⁺=chain[:εₑ⁺].data[idx,1],εₑ⁻ = chain[:εₑ⁻].data[idx,1], frac_sol_am_pos = chain[:frac_sol_am_pos].data[idx,1],frac_sol_am_neg=chain[:frac_sol_am_neg].data[idx,1])
end

function calc_likelihood(chain,sample_size)
    N = length(chain)
    rmse = zeros(N)
    for i in 1:N
        lp = chain[:lp]
        priors = 