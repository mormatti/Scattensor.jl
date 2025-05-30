function log10dim(hs::UniformChain)
    return hs.length * log10(hs.localdim)
end