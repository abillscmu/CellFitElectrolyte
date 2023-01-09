function plot_sampler(chain; label="")
    # Extract values from chain.
    val = get(chain, [:s², :m, :lp])
    ss = link.(Ref(InverseGamma(2, 3)), val.s²)
    ms = val.m
    lps = val.lp

    # How many surface points to sample.
    granularity = 100

    # Range start/stop points.
    spread = 0.5
    σ_start = minimum(ss) - spread * std(ss); σ_stop = maximum(ss) + spread * std(ss);
    μ_start = minimum(ms) - spread * std(ms); μ_stop = maximum(ms) + spread * std(ms);
    σ_rng = collect(range(σ_start, stop=σ_stop, length=granularity))
    μ_rng = collect(range(μ_start, stop=μ_stop, length=granularity))

    # Make surface plot.
    p = StatsPlots.surface(σ_rng, μ_rng, evaluate,
          camera=(30, 65),
        #   ticks=nothing,
          colorbar=false,
          color=:inferno,
          title=label)

    line_range = 1:length(ms)

    StatsPlots.scatter3d!(ss[line_range], ms[line_range], lps[line_range],
        mc =:viridis, marker_z=collect(line_range), msw=0,
        legend=false, colorbar=false, alpha=0.5,
        xlabel="σ", ylabel="μ", zlabel="Log probability",
        title=label)

    return p
end;