function blind_product_matricial(mpo::MPO, mps::MPS)
    substitute_siteinds!(mpo, mps)
    return apply(mpo, mps)
end

function blind_product_matrcial(mpo1::MPO, mpo2::MPS)
    substitute_siteinds!(mpo1, mpo2)
    return apply(mpo1, mpo2)
end

export blind_product_matricial