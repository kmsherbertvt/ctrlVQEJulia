#= I'd like to know if multiplying scalar by vector should have . or no.

Also a good excuse to experiment with time testing in Julia. ^_^
=#


N = 9999999
T = 12

function multiply(scalar, vector)
    return scalar * vector
end

function dotmultiply(scalar, vector)
    return scalar .* vector
end

scalar = (-im*10/1999)
vector = [i for i ∈ 1:N]

# COMPILATION RUNS
@time multiply(scalar, vector)
@time dotmultiply(scalar, vector)


println("----------------------------------------------------")
for trial ∈ 1:T
    @time multiply(scalar, vector)
end
println("----------------------------------------------------")
for trial ∈ 1:T
    @time multiply(scalar, vector)
end


# CONCLUSION: Naturally it doesn't matter. ^_^
