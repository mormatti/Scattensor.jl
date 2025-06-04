function blind_product_inner(mps1::MPS, mps2::MPS)
    replace_siteinds!(mps2, siteinds(mps1))
    return inner(mps1', mps2)
end

export blind_product_inner