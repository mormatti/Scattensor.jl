
function wannier(args...)
    # Def. The (unique) groundstate
    groundstate = states[0 // 1][1]

    # Def. The single-particle dispersion relation, a dictionary momentum => state
    singleparticle = Dict{Rational, Vector{Complex}}()
    
    for (kl, st) in states
        if kl != 0 // 1
            singleparticle[kl] = st[1]
        end
    end
end

# Step 3: l'utente conferma che la Wannier va bene, e da questa calcola l'operatore di
# creazione della particella, lo applica al vuoto DMRG e calcola le funzioni di
# correlazione da cui si possono estrarre le simulazioni di scattering.