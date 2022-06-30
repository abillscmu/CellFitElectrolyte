struct cache{T}
    A::Array{T,2}
    du_transport::Array{T,1}
    control::Array{T,1}
    controller::Array{T,1}
    mm_cache::Array{T,1}
end

function initialize_cache(T)
    A = zeros(T,7,7)
    du_transport = zeros(T,7)
    control = zeros(T,7)
    controller = T.([-1,0,1,0,-1,1,0]./F)
    mm_cache = zeros(T,7)

    return cache(A,du_transport,control,controller,mm_cache)
end
