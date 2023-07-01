
using PythonPlot, MultiKDE, Turing, ProgressMeter
p(x₁, x₂) = exp(-((100*(x₂ - x₁^2)^2) + (1 - x₁)^2)/20)


x_1r = range(start=-10, stop=10, length=25)
x_2r = range(start=0, stop=60, length=25)

x_g = [[x_1, x_2] for x_1 in x_1r for x_2 in x_2r]


@model function rosenbrock(a, b)
    NormalizingConstant = sqrt(a*b)/π
    x_1 ~ Uniform(-15, 15)
    x_2 ~ Uniform(-10, 165)
    Turing.@addlogprob! log(p(x_1, x_2)/NormalizingConstant)
    return nothing
end

a = 1/20
b = 100/20
NormalizingConstant = sqrt(a*b)/π
model = rosenbrock(a, b)


N = 1000000

chain = sample(model, NUTS(0.65), MCMCSerial(), N, 1; progress=true)


observations = Vector{Vector{Float64}}(undef, N)
@showprogress "getting observations" for i in 1:N
    observations[i] = [chain[:x_1].data[i,1], chain[:x_2].data[i,1]]
end

dims = [ContinuousDim(), ContinuousDim()]
bw = [1, 1]
kde = KDEMulti(dims, bw, observations)



z = zeros(length(x_g))
@showprogress "computing pd on grid" for (i,x) in enumerate(x_g)
    z[i] = MultiKDE.pdf(kde, x)
end
z = reshape(z, 25, 25)

z2 = zeros(length(x_g))
@showprogress "computing true z" for (i,x) in enumerate(x_g)
    z2[i] = p(x[1], x[2])/NormalizingConstant
end
z2 = reshape(z2, 25, 25)

figure(1)
clf()
scatter(chain[:x_1].data[:,1], chain[:x_2].data[:,1])
contour(x_1r, x_2r, z)
xlim([-5,5])
ylim([0, 20])


