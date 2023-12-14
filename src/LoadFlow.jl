module LoadFlow

using NLPModels, NLPModelsJuMP, JuMP, CommonOPF
using AxisArrays, SparseArrays

export
  Network,
  busses,
  # end CommonOPF stuff
  add_variables_bfm,
  define_power_with_admittance_bfm,
  set_loads_bfm,
  add_variables_bim,
  set_loads_bim,
  define_power_with_admittance_bim

include("admittance.jl")
include("bfm.jl")
include("bim.jl")

end # module LoadFlow
