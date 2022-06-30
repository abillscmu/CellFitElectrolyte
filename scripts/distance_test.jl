N = 10000


x_vec = -1 .+2 .*rand(N)
y_vec = -1 .+2 .*rand(N)

Rs = 1
Rb = 0.5

my_distance = 


inner_circle = sqrt.(x_vec.^2 + y_vec.^2) .< Rb
outer_circle = (sqrt.(x_vec.^2 + y_vec.^2) .> Rb).&(sqrt.(x_vec.^2 + y_vec.^2) .< Rs)


x_inner = x_vec[inner_circle]
y_inner = y_vec[inner_circle]
x_outer = x_vec[outer_circle]
y_outer = y_vec[outer_circle]

N_inner = length(x_inner)
N_outer = length(x_outer)

