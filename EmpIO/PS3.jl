using DataFrames
using FastGaussQuadrature
using JuMP
using Ipopt

data = readtable("E:\\Documents\\GitHub\\NYU\\EmpIO\\data_ps2.txt", header=false, separator=',')
rename!(data, names(data), [:car, :year, :firm, :price, :quantity, :weight, :hp, :ac, :nest3, :nest4]);
data[:share] = data[:quantity] / 1e8;
meanweight = mean(data[:price])
meanhp = mean(data[:price])
meanprice = mean(data[:price])

data[:weight] = data[:weight] / mean(data[:weight])
data[:hp]     = data[:hp] / mean(data[:hp])
data[:price]  = data[:price] / mean(data[:price]);

function competitor_mean_characteristic(good, char::Symbol)
    mean(data[(data[:firm] .!= good[:firm]) .* (data[:year] .== good[:year]), char]) #same year, different firm
end

function firm_mean_characteristic(good, char::Symbol)
    mean(data[(data[:firm] .== good[:firm]) .* (data[:year] .== good[:year]) .*
            (data[:car] .!= good[:car]), char]) #same year, same firm, different product
end

data[:comp_weight] = [competitor_mean_characteristic(good, :weight) for good in eachrow(data)]
data[:comp_hp]     = [competitor_mean_characteristic(good, :hp) for good in eachrow(data)]
data[:comp_ac]     = [competitor_mean_characteristic(good, :ac) for good in eachrow(data)]

data[:firm_weight] = [firm_mean_characteristic(good, :weight) for good in eachrow(data)]
data[:firm_hp]     = [firm_mean_characteristic(good, :hp) for good in eachrow(data)]
data[:firm_ac]     = [firm_mean_characteristic(good, :ac) for good in eachrow(data)];

data[isnan.(data[:firm_weight]),:firm_weight] = 0.0
data[isnan.(data[:firm_hp]),:firm_hp] = 0.0
data[isnan.(data[:firm_ac]),:firm_ac] = 0.0;

year = 1990

w = convert(Array{Float64}, data[data[:year] .== year, [:share, :weight, :hp, :ac, :price]])
z = convert(Array{Float64}, data[data[:year] .== year, [:weight, :hp, :ac, :comp_weight, :comp_hp, :comp_ac,
                                                                      :firm_weight, :firm_hp, :firm_ac]]);
share = w[:, 1]
weight = w[:, 2]
hp = w[:, 3]
ac = w[:, 4]
price = w[:, 5];

m = Model(solver=IpoptSolver(print_level=1))

@variable(m, β[1:3], start = 0.0)
@variable(m, μy, start = -0.1)
@variable(m, σy >= 0, start = 0.1)

ξ = @variable(m, [i=1:131], start=share[i])
#@variable(m, int[1:131], start=0)

@variable(m, denom[1:30])
@variable(m, integrand[1:131, 1:30])

(e0, q0) = gausshermite(30)
#@NLparameter(m, e[k=1:30] == e0[k])

@objective(m, Min, (z' * ξ)' * (z' * ξ))

@constraint(m, integral_constraint[i=1:131], share[i] == dot(integrand[i, :], q0))
#dot((exp(ξ[i]) + weight[i] * β[1] + hp[i] * β[2] + ac[i] * β[3] - exp(-μy -(√2) * σy * e) * price[i]) / (1 + sum((exp(ξ[j]) + weight[j] * β[1] + hp[j] * β[2] + ac[j] * β[3] - exp(-μy -(√2) * σy * e) * price[j]) for j in 1:131)), q)


for k in 1:30
    @NLconstraint(m, denom[k] == (1 + sum((exp(ξ[j]) + weight[j] * β[1] + hp[j] * β[2] + ac[j] * β[3] - exp(-μy - sqrt(2) * σy * e0[k]) * price[j]) for j in 1:131)))
end

for i in 1:131, k in 1:30
    @NLconstraint(m, integrand[i,k] == (exp(ξ[i]) + weight[i] * β[1] + hp[i] * β[2] + ac[i] * β[3] - exp(-μy -sqrt(2) * σy * e0[k]) * price[i]) / denom[k])
end
