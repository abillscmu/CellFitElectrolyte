function build_nn(num_layers, activation, width, input_dim, output_dim)
    first_layer = FastDense(input_dim, width, activation)
    layer_vec = Any[first_layer,]
    for layer in 2:num_layers-1
        this_layer = FastDense(width, width, activation)
        push!(layer_vec, this_layer)
    end
    last_layer = FastDense(width, output_dim)
    push!(layer_vec, last_layer)
    nn = FastChain(layer_vec...)
    return nn
end


