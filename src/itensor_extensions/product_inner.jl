function product_inner(mps1::MPS, mps2::MPS)
    replace_siteinds!(mps2, siteinds(mps1))
    return inner(mps1', mps2)
end

export product_inner