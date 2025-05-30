function logdim(hs::UniformChain)
    return hs.length * log(hs.localdim)
end