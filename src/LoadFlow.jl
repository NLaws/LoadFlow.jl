module LoadFlow

using NLPModels, NLPModelsJuMP, JuMP, CommonOPF

export
  add_variables_bfm,
  define_power_with_admittance_bfm,
  set_loads_bfm,
  add_variables_bim,
  define_power_with_admittance_bim,
#=
No timesteps, simplest load flow model possible
See if NLPModelsJuMP works with complex variables: nope :/
    # ERROR: Unexpected object im of type Complex{Bool} in nonlinear expression.

=#

include("bfm.jl")
include("bim.jl")

netdict = Dict(
    :network => Dict(:substation_bus => "1"),
    :conductors => [
        Dict(
            :busses => ("1", "2"),
            :r1 => 0.1,
            :x1 => 0.1,
            :length => 100
        ),
        Dict(
            :busses => ("2", "3"),
            :r1 => 0.1,
            :x1 => 0.1,
            :length => 100
        ),
        Dict(
            :busses => ("3", "1"),
            :r1 => 0.1,
            :x1 => 0.1,
            :length => 100
        ),
    ],
    :loads => [
        Dict(
            :bus => "3",
            :kws1 => [5],
            :kvars1 => [0.5]
        )
    ]
)

net = Network(netdict)

#=
using Ipopt, LoadFlow, JuMP, NLPModelsJuMP, CommonOPF

m = JuMP.Model(Ipopt.Optimizer)
set_attribute(m, "max_iter", 10_000)

add_variables_bim(m, net)
define_power_with_admittance_bim(m, net)
set_knowns_bim(m, net)

julia> all_variables(m)
12-element Vector{VariableRef}:
 pj[1]
 pj[2]
 pj[3]
 qj[1]
 qj[2]
 qj[3]
 v_mag[1]
 v_mag[2]
 v_mag[3]
 v_ang[1]
 v_ang[2]
 v_ang[3]

nlp = MathOptNLPModel(m)

x0 = copy(nlp.meta.x0)
for i=7:9
    x0[i] = 1.0
end
x0[3] = -5
x0[6] = -0.5
 
Jv = jprod(nlp, x0, x0)

Jx = jac(nlp, x0)


add_variables_bfm(m, net)

define_power_with_admittance_bfm(m, net)

set_loads_bfm(m, net)

optimize!(m)

julia> all_variables(m)
20-element Vector{VariableRef}:
 pij[("1", "2")]
 pij[("1", "3")]
 pij[("2", "3")]
 pij[("2", "1")]
 pij[("3", "1")]
 pij[("3", "2")]
 qij[("1", "2")]
 qij[("1", "3")]
 qij[("2", "3")]
 qij[("2", "1")]
 qij[("3", "1")]
 qij[("3", "2")]
 v_real[1]
 v_real[2]
 v_real[3]
 v_imag[1]
 v_imag[2]
 v_imag[3]
 p0
 q0

r = Dict()
for varref in all_variables(m)
    r[replace(string(varref), "\"" => "")] = value(varref)
end
r



nlp = MathOptNLPModel(m)


x0 = copy(nlp.meta.x0)
for i=13:15
    x0[i] = 1.0
end
x0[17] = 5
x0[18] = 0.5



f(x) = obj(nlp, x)
g(x) = grad(nlp, x)
H(x) = hess(nlp, x)

d = -H(x0) \ g(x0)

# ERROR: LinearAlgebra.SingularException(0)


for i = 1:5
    global x
    x = x - H(x) \ g(x)
    println("x = $x")
end



julia> typeof(net.graph.graph) <: SimpleGraph
true

julia> typeof(net.graph.graph) <: SimpleDiGraph
false

julia> SimpleDiGraph <: SimpleGraph
false


=#
end # module LoadFlow
