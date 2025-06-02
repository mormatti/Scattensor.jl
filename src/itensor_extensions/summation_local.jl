function summation_local(mpo::MPO, d::Integer, L::Integer; func::Function = identity; kwargs...)
    L0 = length(mpo)
    Lc = L - L0
    if func == identity
        f(j) = 1
        func = f
    end
    finalsum = func(0) * insert_local(0, mpo, Lc, d)
    for j in 1:Lc
        tosum = func(j) * insert_local(j, mpo, Lc - j, d)
        substitute_siteinds!(tosum, finalsum)
        finalsum += tosum
    end
    finalsum = truncate(finalsum, kwargs...)
end

export summation_local