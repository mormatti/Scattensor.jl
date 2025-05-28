function local_expvals(ψ::MPS, A₀::MPO, d::dType) where {dType <: Integer}
    if length(A₀) < length(ψ)
        L₀ = length(A₀)
        L = length(ψ)
        Lc = L - L₀
        vals = []
        for j in 1:Lc
            Aext = insert_local(j, A₀, Lc - j, d)
            substitute_siteinds!(Aext, ψ)
            push!(vals, real(inner(ψ', Aext, ψ)))
        end
        return vals
    elseif length(A₀) == length(ψ)
        return real(inner(mps', A₀, mps))
    else
        error("The length of the local operator cannot be greater than the one of the .")
    end
end

function local_expvals(ψ::Vector{MPS}, A₀::MPO, d::dType) where {dType <: Integer}
    # We assert that all the MPS inside ψ has the same length
    @assert all(mps -> length(mps) == length(first(ψ)), ψ) "All the MPS must have the same length"
    println("")

    timelapsed = 0
    matrix = zeros(length(ψ), length(ψ[1]) - length(A₀))
    N = length(ψ)
    for n in 1:N
        print("Step $n / $N. Estimated remaining time: $(timelapsed * (N - n)).")
        timelapsed = @elapsed begin
        matrix[n,:] = local_expvals(ψ[n], A₀, d)
        end
        print("\r\u001b[2K")
    end
    return matrix
end

export local_expvals