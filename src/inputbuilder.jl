project_path(parts...) = normpath(joinpath(@__DIR__, "..", parts...))
function load_airbus_cyclearrays()
    vec = load(project_path("src/cycle_array_vector_new.jld2"))
    return vec
end 

#For now, I will make a simple 5 amp discharge for 10 minutes
function current_profile(current,time)
    #Builds a properly formatted input vector. The input vector should be model-agnostic.
        L=length(current)
        time_vec = time
        classification_vec = 3 .*ones(length(current))
        value_vec = current
        p_load = vcat(L,time_vec,classification_vec,value_vec)
        return p_load
    end

