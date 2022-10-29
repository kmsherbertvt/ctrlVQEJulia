#= Run sequence of tests to ensure time evolution of a basic square pulse is functional.

NOTE: This is not (yet) meant to be a fully robust unit-test framework.
=#

include("../src/utils.jl")
include("../src/device.jl")
include("../src/pulse.jl")
include("../src/evolve.jl")


module UnitTesting
using Test

import ..Utils
import ..Pulses
import ..Devices
import ..Evolutions

@testset "Utils" begin
    @test Utils.a_matrix()  == [0 1; 0 0]
    @test Utils.a_matrix(3) == [0 1 0; 0 0 √2; 0 0 0]
    a = Utils.a_matrix()
    @test Utils.on(a, 1, 2) == [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
    @test Utils.on(a, 2, 2) == [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0]
end

@testset "Basic Square Pulse" begin
    pt = Pulses.BasicSquarePulse(10, 20, [1, 2, 3], [2.5, 6.5])
    @test length(pt) == 10
    @test Pulses.frequency(pt, 0)   == 20
    @test Pulses.frequency(pt, 2.5) == 20
    @test Pulses.frequency(pt, 5.0) == 20
    @test Pulses.frequency(pt, 7.5) == 20
    @test Pulses.amplitude(pt, 0)   == 0
    @test Pulses.amplitude(pt, 2.5) == 2
    @test Pulses.amplitude(pt, 5.0) == 2
    @test Pulses.amplitude(pt, 7.5) == 3
end

@testset "Transmon Device" begin
    device = Devices.Transmon([1,2], [.1, .2], Dict(Devices.QubitCouple(1,2)=>.5))
    @test length(device) == 2
    @test size(Devices.static_hamiltonian(device, 3)) == (9, 9)
    HD = Devices.static_hamiltonian(device, 3)
    @test Devices.static_hamiltonian(device) == [0 0 0 0; 0 2 0.5 0; 0 0.5 1 0; 0 0 0 3]
    # Technically we're never testing δ terms, but I'm not inclined to enter a 9x9 matrix...
end

# TODO: Port evolvetest.jl into a concise(?) test script here.

end
