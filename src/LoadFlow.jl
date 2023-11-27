module LoadFlow

using NLPModels, NLPModelsJuMP, JuMP, CommonOPF

export
  add_variables,
  define_power,
  set_loads
#=
No timesteps, simplest load flow model possible
See if NLPModelsJuMP works with complex variables: nope :/
    # ERROR: Unexpected object im of type Complex{Bool} in nonlinear expression.

=#


function add_variables(m::JuMP.AbstractModel, net::Network{SinglePhase})
    edges_both_ways = vcat(collect(edges(net)), reverse.(edges(net)))
    all_busses = Set(busses(net))
    b0 = net.substation_bus
    busses_except_sub = setdiff(all_busses, [b0])
    @variables(m, begin
        pij[edges_both_ways]  # should undirected net return duplicate edges? i.e. (i,j) and (j,i) ?
        qij[edges_both_ways]
        # 1.2 >= v_real[busses_except_sub] >= 0.8, (start = 1.0, container = Dict{String, Union{VariableRef, Real}})
        # 1.2 >= v_imag[busses_except_sub] >= 0.8, (start = 1.0, container = Dict{String, Union{VariableRef, Real}})
        v_real[busses(net)] >= 0.0, (start = 1.0)
        v_imag[busses(net)], (start = 0.0)
        p0
        q0
    end)
    # we define the container type for voltages so that we can store substation voltage in the model
    # like a variable; but doing this makes the Ipopt problem infeasible
    # m[:v_real][b0] = net.v0
    # m[:v_imag][b0] = 0.0

    n_edge_vars = 4*length(edges(net))
    n_bus_vars = 2*length(busses_except_sub)

    @info "Variables count:\n"*
        "\t$(n_edge_vars) edge variables\n"*
        "\t$(n_bus_vars) bus variables\n"*
        "\t2 slack bus variables\n"*
        "\t-------------------------\n"*
        "\t$(n_edge_vars + n_bus_vars + 2) total variables"
    nothing
end


function define_power_with_impedance(m::JuMP.AbstractModel, net::Network{SinglePhase})
    pij = m[:pij]
    qij = m[:qij]
    v_real = m[:v_real]
    v_imag = m[:v_imag]
    edges_both_ways = vcat(collect(edges(net)), reverse.(edges(net)))
    # note cannot put functions in NLconstraint unless it is truly part of the math
    @NLconstraint(m, [(i,j) in edges_both_ways],
        pij[(i,j)] * net[(i,j)][:Conductor][:r1] * net[(i,j)][:Conductor][:length] / net.Zbase
        + qij[(i,j)] * net[(i,j)][:Conductor][:x1] * net[(i,j)][:Conductor][:length] / net.Zbase
        ==
        v_real[i]^2 - v_real[i] * v_real[j] - v_imag[i] * v_imag[j] + v_imag[i]^2
    )

    @NLconstraint(m, [(i,j) in edges_both_ways],
        qij[(i,j)] * net[(i,j)][:Conductor][:r1] * net[(i,j)][:Conductor][:length] / net.Zbase
        + pij[(i,j)] * net[(i,j)][:Conductor][:x1] * net[(i,j)][:Conductor][:length] / net.Zbase
        == 
        v_real[i] * v_imag[j] - v_imag[i] * v_real[j]
    )
    nothing
end


function define_power_with_admittance(m::JuMP.AbstractModel, net::Network{SinglePhase})
    pij = m[:pij]
    qij = m[:qij]
    v_real = m[:v_real]
    v_imag = m[:v_imag]
    edges_both_ways = vcat(collect(edges(net)), reverse.(edges(net)))
    # note cannot put functions in NLconstraint unless it is truly part of the math
    @NLconstraint(m, [(i,j) in edges_both_ways],
        pij[(i,j)] ==
        (v_real[i]^2 - v_real[i] * v_real[j] - v_imag[i] * v_imag[j] + v_imag[i]^2) *
        net[(i,j)][:Conductor][:r1] / (net[(i,j)][:Conductor][:r1]^2 + net[(i,j)][:Conductor][:x1]^2) *
        net.Zbase / net[(i,j)][:Conductor][:length] 
        +
        (v_real[i] * v_imag[j] - v_imag[i] * v_real[j]) * 
        -net[(i,j)][:Conductor][:x1] / (net[(i,j)][:Conductor][:r1]^2 + net[(i,j)][:Conductor][:x1]^2) *
        net.Zbase / net[(i,j)][:Conductor][:length] 
    )

    @NLconstraint(m, [(i,j) in edges_both_ways],
        qij[(i,j)] ==
        (v_real[i] * v_imag[j] - v_imag[i] * v_real[j]) *
        net[(i,j)][:Conductor][:r1] / (net[(i,j)][:Conductor][:r1]^2 + net[(i,j)][:Conductor][:x1]^2) *
        net.Zbase / net[(i,j)][:Conductor][:length] 
        -
        (v_real[i]^2 - v_real[i] * v_real[j] - v_imag[i] * v_imag[j] + v_imag[i]^2) *
        -net[(i,j)][:Conductor][:x1] / (net[(i,j)][:Conductor][:r1]^2 + net[(i,j)][:Conductor][:x1]^2) *
        net.Zbase / net[(i,j)][:Conductor][:length] 
    )
    nothing
end


function set_loads(m::JuMP.AbstractModel, net::Network{SinglePhase})
    pij = m[:pij]
    qij = m[:qij]
    # assuming that substation_bus is not in load busses for now
    @constraint(m, [i in real_load_busses(net)],
        -net[i][:Load][:kws1][1] * 1e3 / net.Sbase ==
        sum(
            pij[(i,j)] for j in j_to_k(i, net)
        )
    )
    @constraint(m, [i in reactive_load_busses(net)],
        -net[i][:Load][:kvars1][1] * 1e3 / net.Sbase ==
        sum(
            qij[(i,j)] for j in j_to_k(i, net)
        )
    )
    all_busses = Set(busses(net))
    b0 = net.substation_bus
    non_load_busses_real = setdiff(all_busses, [b0, collect(real_load_busses(net))...])
    non_load_busses_reactive = setdiff(all_busses, [b0, collect(reactive_load_busses(net))...])
    # TODO put the sets in CommonOPF
    @constraint(m, [i in non_load_busses_real],
        0 ==
        sum(
            pij[(i,j)] for j in j_to_k(i, net)
        )
    )
    @constraint(m, [i in non_load_busses_reactive],
        0 ==
        sum(
            qij[(i,j)] for j in j_to_k(i, net)
        )
    )
    @constraint(m, m[:p0] == sum( pij[(b0,j)] for j in j_to_k(b0, net) ) )
    @constraint(m, m[:q0] == sum( qij[(b0,j)] for j in j_to_k(b0, net) ) )
    # NOTE for undirected graph i_to_j == j_to_k == outneighbors == inneighbors (TODO test this)
    nothing
end


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

m = JuMP.Model(Ipopt.Optimizer)
set_attribute(m, "max_iter", 10_000)

add_variables(m, net)

define_power_with_admittance(m, net)

set_loads(m, net)

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
