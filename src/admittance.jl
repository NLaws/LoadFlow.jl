

"""

line admittance
"""
function yij(i::AbstractString, j::AbstractString, net::CommonOPF.Network)::Tuple{Real, Real}
    r = rij(i, j, net)
    x = xij(i, j, net)
    denom = r^2 + x^2
    return r/denom, -x/denom
end


"""

entry of admittance matrix at i,j
"""
function Yij(i::AbstractString, j::AbstractString, net::CommonOPF.Network)::Tuple{Real, Real}
    if i != j
        return -1 .* yij(i, j, net)
    end
    # TODO shunt impedance
    # sum up the y_jk where k is connected to j
    g, b = 0, 0
    for k in union(i_to_j(j, net), j_to_k(j, net))
        (g, b) = (g, b) .+ yij(i, k, net)
    end
    return g, b
end


"""

NLPModelsJuMP does not allow or functions in NLexpression, etc
so we have to build the admittance matrix rather than look up values as needed.

"""
function build_Y(net::CommonOPF.Network)
    # I = Int[]; J = Int[]; V = Real[];
    # ibmap = net.graph.graph_data[:int_bus_map]
    # TODO sparse array?
    bs = collect(busses(net))
    N = length(bs)
    G = AxisArray(spzeros(N, N), bs, bs)
    B = AxisArray(spzeros(N, N), bs, bs)
    for j in bs
        (G[j, j], B[j, j]) = Yij(j, j, net)
        for i in i_to_j(j, net)
            (G[i, j], B[i, j]) = Yij(i, j, net)
        end
    end
    return G, B
end
