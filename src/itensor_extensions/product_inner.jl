function product_inner(mps1::MPS, mps2::MPS)
    return inner(mps1', replace_siteinds(mps2, siteinds(mps1)))
end

export product_inner