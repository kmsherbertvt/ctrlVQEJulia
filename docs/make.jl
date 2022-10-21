include("../src/utils.jl")
include("../src/device.jl")
include("../src/pulse.jl")
include("../src/evolve.jl")
using Documenter, .Devices, .Evolutions, .Pulses, .Utils

makedocs(sitename="ctrl-VQE Julia m1")
