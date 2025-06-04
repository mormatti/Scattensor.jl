function blind_product_outer(mps1::MPS, mps2::MPS)
    return outer(mps1', replace_siteinds(mps2, siteinds(mps1)))
end

export blind_product_outer