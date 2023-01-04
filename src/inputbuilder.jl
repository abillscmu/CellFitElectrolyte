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

function truncate_cycle_array(num_segments, cycle_array)
    num_segments_original = Int(cycle_array[1])
    if num_segments_original == -1
        error("can't cut an RPT")
    elseif num_segments_original < num_segments
        error("not enough segments in original cycle array")
    end
    new_times = cycle_array[2:1+num_segments]
    new_types = cycle_array[num_segments_original+2:num_segments_original+num_segments+1]
    new_values = cycle_array[2*num_segments_original+2:2*num_segments_original+num_segments+1]
    new_cycle_array = vcat(num_segments, new_times, new_types, new_values)
    return new_cycle_array
end

