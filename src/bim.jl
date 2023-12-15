


function add_variables_bim(m::JuMP.AbstractModel, net::Network{SinglePhase})
    # assuming all PQ busses for now
    m[:v_mag] = Dict{String, Any}()
    m[:v_ang] = Dict{String, Any}()

    b0 = net.substation_bus
    m[:v_mag][b0] = net.v0
    m[:v_ang][b0] = 0.0
    vreg_busses = voltage_regulator_busses(net)
    # vreg_busses have v_mag and p specified (or p gets set to zero)

    for b in setdiff(busses(net), [b0])
        if !(b in vreg_busses)
            m[:v_mag][b] = @variable(
                m, lower_bound=0.5, upper_bound=1.5, 
                start = 1.0, 
                base_name="v_mag[" * string(b) * "]"
            )
        else
            m[:v_mag][b] = net[b][:VoltageRegulator][:vreg_pu]
        end
        m[:v_ang][b] = @variable(
            m, lower_bound=0.0, upper_bound=2π, 
            start = 0.0, 
            base_name="v_ang[" * string(b) * "]"
        )
    end
end


"""
    set_loads_bim(m::JuMP.AbstractModel, net::Network{SinglePhase})

Fill in `m[:p]` and `m[:q]` with one of 
    - a provided value,
    - zero,
    - or a variable.
"""
function set_loads_bim(m::JuMP.AbstractModel, net::Network{SinglePhase})
    p = m[:p] = Dict{String, Number}()
    q = m[:q] = Dict{String, Union{Number, JuMP.AbstractVariableRef}}()
    # assuming that substation_bus is not in load busses for now
    for j in real_load_busses(net)
        p[j] = -net[j][:Load][:kws1][1] * 1e3 / net.Sbase
    end
    for j in reactive_load_busses(net)
        q[j] = -net[j][:Load][:kvars1][1] * 1e3 / net.Sbase
    end
    all_busses = Set(busses(net))
    b0 = net.substation_bus
    non_load_busses_real = setdiff(all_busses, [b0, collect(real_load_busses(net))...])
    non_load_busses_reactive = setdiff(all_busses, [b0, collect(reactive_load_busses(net))...])
    # TODO put the sets in CommonOPF
    for j in non_load_busses_real
        p[j] = 0
    end
    vreg_busses = voltage_regulator_busses(net)
    for j in non_load_busses_reactive
        if j in vreg_busses
            q[j] = @variable(m, base_name="q[" * string(j) * "]")
        else
            q[j] = 0
        end
    end
    
    nothing
end


function define_power_with_admittance_bim(m::JuMP.AbstractModel, net::Network{SinglePhase})
    p = m[:p]
    q = m[:q]
    v_mag = m[:v_mag]
    v_ang = m[:v_ang]
    busses_no_sub = setdiff(busses(net), [net.substation_bus])

    G, B = build_Y(net)

    # NOTE cannot put functions in NLconstraint unless it is truly part of the math
    # NOTE NLPModelsJuMP does not work with complex variables
    @NLexpression(m, deltaP[j in busses_no_sub],
        -p[j] + (
            v_mag[j] 
            * sum( v_mag[i] * (
                    G[i,j] # conductance
                *   cos(v_ang[j] - v_ang[i])
                +   B[i,j] # susceptance
                *   sin(v_ang[j] - v_ang[i])
            ) for i in union(i_to_j(j, net), j_to_k(j, net), [j])
            )
        )
    )

    @NLexpression(m, deltaQ[j in busses_no_sub],
        -q[j] + (
            v_mag[j] 
            * sum( v_mag[i] * (
                    G[i,j] # conductance
                *   sin(v_ang[j] - v_ang[i])
                -   B[i,j] # susceptance
                *   cos(v_ang[j] - v_ang[i])
            ) for i in union(i_to_j(j, net), j_to_k(j, net), [j])
            )
        )
    )

    @NLconstraint(m, [j in busses_no_sub], deltaP[j] == 0)
    @NLconstraint(m, [j in busses_no_sub], deltaQ[j] == 0)
    nothing
end
