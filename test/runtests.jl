using LoadFlow
using CommonOPF
using Ipopt, JuMP
using NLPModelsJuMP, NLPModels


@testset "" begin


# we need to solve:
# -J Δxⁱ = f(xⁱ)
# for Δxⁱ to get
# x⁽ⁱ⁺⁾ = xⁱ + Δxⁱ
# until f(xⁱ) is within some tolerance (or Δxⁱ). f(xⁱ) contains the deltaP and deltaQ so it should go
# to zero along with Δxⁱ 
    
netdict = Dict(
    :network => Dict(:substation_bus => "1", :Sbase => 1e3),
    :conductors => [
        Dict(
            :busses => ("1", "2"),
            :r1 => 0.0,
            :x1 => 0.1,
            :length => 1
        ),
        Dict(
            :busses => ("2", "3"),
            :r1 => 0.0,
            :x1 => 0.1,
            :length => 1
        ),
        Dict(
            :busses => ("3", "1"),
            :r1 => 0.0,
            :x1 => 0.1,
            :length => 1
        ),
    ],
    :loads => [
        Dict(
            :bus => "3",
            :kws1 => [.0028653],
            :kvars1 => [.0012244]
        ),
        Dict(
            :bus => "2",
            :kws1 => [-0.0006661],
            :kvars1 => [-.0016395]
        )
    ]
)

net = Network(netdict)


m = JuMP.Model(Ipopt.Optimizer)
set_attribute(m, "max_iter", 10_000)

add_variables_bim(m, net)
set_loads_bim(m, net)
define_power_with_admittance_bim(m, net)


nlp = MathOptNLPModel(m)

busses_no_sub = setdiff(busses(net), [net.substation_bus])


x = Vector{Vector{Real}}()
push!(x, copy(nlp.meta.x0))

for i = 1:5
    J = jac(nlp, x[i])
    fxi = Vector{Real}()
    inp = Dict(
        var_ref => var_val for (var_ref, var_val) in zip(all_variables(m), x[i])
    )
    for b in busses_no_sub
        push!(fxi, value(k -> inp[k], m[:deltaP][b]))
    end
    for b in busses_no_sub
        push!(fxi, value(k -> inp[k], m[:deltaQ][b]))
    end
    Δxⁱ = J \ -fxi
    push!(x, x[i] + Δxⁱ)
end


## another way:

# for i = 1:5
#     J = jac(nlp, x[i])
#     fxi = Vector{Real}()
#     for b in busses_no_sub
#         push!(fxi, value(k -> get(inp, k, 0.0), m[:deltaP][b]))
#     end
#     for b in busses_no_sub
#         push!(fxi, value(k -> get(inp, k, 0.0), m[:deltaQ][b]))
#     end


#     w = lu!(J)

#     Δxⁱ = ldiv!(w, fxi)
#     push!(x, x[i] - Δxⁱ)
# end


end