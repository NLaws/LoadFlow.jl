


function add_variables_bfm(m::JuMP.AbstractModel, net::Network{SinglePhase})
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


function define_power_with_impedance_bfm(m::JuMP.AbstractModel, net::Network{SinglePhase})
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


function define_power_with_admittance_bfm(m::JuMP.AbstractModel, net::Network{SinglePhase})
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


function set_loads_bfm(m::JuMP.AbstractModel, net::Network{SinglePhase})
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
