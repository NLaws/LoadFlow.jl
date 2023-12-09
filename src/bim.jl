


function add_variables_bim(m::JuMP.AbstractModel, net::Network{SinglePhase})
    @variables(m, begin
        pj[busses(net)]  # we specify all PQ busses values later
        qj[busses(net)]
        v_mag[busses(net)] >= 0.0, (start = 1.0)
        v_ang[busses(net)], (start = 0.0)
    end)
end


function define_power_with_admittance_bim(m::JuMP.AbstractModel, net::Network{SinglePhase})
    pj = m[:pj]
    qj = m[:qj]
    v_mag = m[:v_mag]
    v_ang = m[:v_ang]

    # note cannot put functions in NLconstraint unless it is truly part of the math
    @NLconstraint(m, [j in busses(net)],
        pj[j] == v_mag[j] * sum( v_mag[i] * (
                net[(i,j)][:Conductor][:r1] / (net[(i,j)][:Conductor][:r1]^2 + net[(i,j)][:Conductor][:x1]^2) *
                net.Zbase / net[(i,j)][:Conductor][:length] # conductance per-unit
            *   cos(v_ang[j] - v_ang[i])
            +   net[(i,j)][:Conductor][:x1] / (net[(i,j)][:Conductor][:r1]^2 + net[(i,j)][:Conductor][:x1]^2) *
                net.Zbase / net[(i,j)][:Conductor][:length] # susceptance per-unit
            *   sin(v_ang[j] - v_ang[i])
            )
        for i in union(i_to_j(j, net), j_to_k(j, net))
        )
    )

    @NLconstraint(m, [j in busses(net)],
        qj[j] == v_mag[j] * sum( v_mag[i] * (
                net[(i,j)][:Conductor][:r1] / (net[(i,j)][:Conductor][:r1]^2 + net[(i,j)][:Conductor][:x1]^2) *
                net.Zbase / net[(i,j)][:Conductor][:length] # conductance per-unit
            *   sin(v_ang[j] - v_ang[i])
            -   net[(i,j)][:Conductor][:x1] / (net[(i,j)][:Conductor][:r1]^2 + net[(i,j)][:Conductor][:x1]^2) *
                net.Zbase / net[(i,j)][:Conductor][:length] # susceptance per-unit
            *   cos(v_ang[j] - v_ang[i])
            )
        for i in union(i_to_j(j, net), j_to_k(j, net))
        )
    )
    nothing
end


function set_knowns_bim(m::JuMP.AbstractModel, net::Network{SinglePhase})
    p = m[:pj]
    q = m[:qj]
    # assuming that substation_bus is not in load busses for now
    @constraint(m, [j in real_load_busses(net)],
        -net[j][:Load][:kws1][1] * 1e3 / net.Sbase == p[j]
    )
    @constraint(m, [j in reactive_load_busses(net)],
        -net[j][:Load][:kvars1][1] * 1e3 / net.Sbase == q[j]
    )
    all_busses = Set(busses(net))
    b0 = net.substation_bus
    non_load_busses_real = setdiff(all_busses, [b0, collect(real_load_busses(net))...])
    non_load_busses_reactive = setdiff(all_busses, [b0, collect(reactive_load_busses(net))...])
    # TODO put the sets in CommonOPF
    @constraint(m, [j in non_load_busses_real],
        0 == p[j]
    )
    @constraint(m, [j in non_load_busses_reactive],
        0 == q[j]
    )

    # substation / slack bus
    @constraint(m, m[:v_mag][b0] == 1.0)
    @constraint(m, m[:v_ang][b0] == 0)
    
    nothing
end
