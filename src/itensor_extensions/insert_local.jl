function insert_local(L1::L1Type, localmpo::MPO, L2::L2Type, d::dType) where {L1Type <: Integer, L2Type <: Integer, dType <: Integer}
    idm(L) = operator_identity(MPO, d, L)
    if L1 == 0 && L2 == 0
        return c
    elseif L1 == 0
        return product_direct(localmpo, idm(L2))
    elseif L2 == 0
        return product_direct(idm(L1), localmpo)
    else
        return product_direct(idm(L1), product_direct(localmpo, idm(L2)))
    end
end 

export insert_local