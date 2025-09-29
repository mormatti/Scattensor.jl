# TODO write documentation

function adjoint_mpo(A::MPO)::MPO
    Ad = dag(A)
    return swapprime(Ad, 0 => 1; tags = "Site")
end

export adjoint_mpo