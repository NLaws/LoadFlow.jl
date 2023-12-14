using LoadFlow
using Ipopt, JuMP
using NLPModelsJuMP, NLPModels
using Test
# TODO will not need all these dependencies in test env eventually

@testset "McCalley ISU example" begin


# we need to solve:
# -J Δxⁱ = f(xⁱ)
# for Δxⁱ to get
# x⁽ⁱ⁺⁾ = xⁱ + Δxⁱ
# until f(xⁱ) is within some tolerance (or Δxⁱ). f(xⁱ) contains the deltaP and deltaQ so it should go
# to zero along with Δxⁱ 
    
netdict = Dict(
    :network => Dict(:substation_bus => "1", :Sbase => 1),
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
        ) # this should be a "generator" with P and V specified
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
x[1][1] = 1.05  # v_mag["2"]

change_tol = 0.0001
change = 1
i = 1
while change > change_tol
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
    change = maximum( Δxⁱ )
    push!(x, x[i] + Δxⁱ)
    i += 1
    x[end][1] = 1.05
end
# @test x[end][1] ≈ 1.05 rtol = change_tol * 10 # v_mag["2"]
@test x[end][2] * 180/pi ≈ -3 rtol = change_tol * 10    # v_ang["2"]
@test x[end][3] ≈ 0.9499 rtol = change_tol * 10   # v_mag["3"]
@test x[end][4] * 180/pi ≈ -10.01 rtol = change_tol * 10   # v_ang["3"]



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