# Electric Hamiltonian Opsum
function OsHE()::OpSum
    H = OpSum()
    for j in 1:L
        H += -1,"Sz",j
    end
    for j in 1:(L-1)
        H += -2,"Sz",j,"Sz",j+1
    end
    return H
end

# Magnetic Hamiltonian Opsum
function OsHB()::OpSum
    H = OpSum()
    for j in 1:L
        H += -√8,"Sx",j
    end
    return H
end

# Electric impurity Opsum
function OsHEimp(p::Integer)::OpSum
    H = OpSum()
    H += -1,"Sz",p
    H += -1,"Sz",p+1
    H += 5/4,"Sz",p,"Sz",p+1
    return H
end

# Magnetic impurity Opsum
function OsHBimp(p::Integer)::OpSum
    H = OpSum()
    c1 = √6 + 2 - √8
    c2 = √6 + 2
    H += -c1,"Sx",p
    H += -c1,"Sx",p+1
    H += -c2,"Sz",p,"Sx",p+1
    H += -c2,"Sx",p,"Sz",p+1
    return H
end

# Local electric Hamiltonian Opsum
function OsHE(j::Integer)::OpSum
    H = OpSum()
    H += -1,"Sz",j
    H += -1,"Sz",j-1,"Sz",j
    H += -1,"Sz",j,"Sz",j+1
    return H
end

# Local magnetic Hamiltonian Opsum
function OsHB(j::Integer)::OpSum
    H = OpSum()
    H += -√8,"Sx",j
    return H
end

# Naive Wannier creation operator OpSum
function OsNW(j₀::Integer)::OpSum
    H = OpSum()
    H += "Sz",j₀
    H += "Sx",j₀
    return H
end

# Naive wavepacket creation operator OpSum
function OsNWp(j₀::Integer, k::Float64, σ::Float64, L::Integer, ϵ::Float64)::OpSum
    H = OpSum()
    f(j) = exp(im * k * j) * exp(-(j - j₀)^2 / (4 * σ^2))
    for j in 1:L
        if abs(f(j)) > ϵ
            H += f(j),"Sz",j
            H += f(j),"Sx",j
        end
    end
    return H
end

function testTDVP()
    # System parameters
    L = 100 # Number of sites
    λ = 2 # Coupling constant
    jⁱᵐᵖ = Integer(round(L/2)) # Impurity position

    # Wavepacket parameters
    jᵂᴾ = Integer(round(L/4)) # Wavevector initial average position
    kᵂᴾ = π/2 # Wavevector initial average momentum
    σᵂᴾ = 2.0 # Wavevector initial standard deviation
    ϵᵂᴾ = 10^(-10) # Wavevector cutoff for the creation operator

    # DMRG parameters
    Nᴰᴹᴿᴳ = 10 # Number of sweeps
    χᴰᴹᴿᴳ = 100 # Maximum bond dimension
    ϵᴰᴹᴿᴳ = 1e-10 # Cutoff

    # Simulation parameters
    χᵀᴰⱽᴾ = 100 # Maximum bond dimension for the time evolution
    ϵᵀᴰⱽᴾ = 1e-10 # Cutoff for the time evolution
    dt = 0.5 # Time step
    Δt = 130.0 # Total time of simulation

    # We create the sites and the Hamiltonian
    os = λ * OsHE() + (1 - λ) * OsHB() + (λ * OsHEimp(jⁱᵐᵖ) + (1 - λ) * OsHBimp(jⁱᵐᵖ))
    sites = siteinds("S=1/2", L)
    H = MPO(os, sites)

    # We find the groundstate
    ℰ₀,ψ₀ = dmrg(H, randomMPS(sites); 
                nsweeps = Nᴰᴹᴿᴳ, 
                maxdim = χᴰᴹᴿᴳ, 
                cutoff = ϵᴰᴹᴿᴳ)

    # We generate the creation operator and we apply it to the groundstate
    Loc = MPO(OsNWp(jᵂᴾ, kᵂᴾ, σᵂᴾ, L, ϵᵂᴾ), sites)
    ψ = noprime(Loc * ψ₀)
    normalize!(ψ)

    # We time-evolve the system
    function tdvp_time_evolution(
        H::MPO,
        ψ::MPS,
        ψ₀::MPS,
        dt::Float64,
        Δt::Float64
        )

        N = Integer(Δt / dt) # Number of steps

        # We create the list of times 
        times = [n * dt for n in 1:N]
        positions = [x for x in 2:L-1]

        # We create a 2D array E of dimensions (N, L-2)
        En = zeros(N, L-2)
        EnLog = zeros(N, L-2)
        Linkdims = zeros(N, L-2)
        Maxlinkdim = zeros(N)

        for n in 1:N
            ψ = tdvp(H, ψ, -im * dt, maxlinkdim = 100, cutoff = 1e-10)
            normalize!(ψ)

            timeElapsed = @elapsed begin
            for j in 2:(L-1)
                locEn = λ * OsHE(j) + (1 - λ) * OsHB(j)
                A = MPO(locEn, sites)
                ComputedLocEn = real(inner(ψ', A, ψ)) - real(inner(ψ₀', A, ψ₀))
                En[n,j-1] = ComputedLocEn
                EnLog[n,j-1] = log10(abs(ComputedLocEn))
                Linkdims[n,j-1] = linkdim(ψ,j)
                Maxlinkdim[n] = maxlinkdim(ψ)
            end
            end

            # Print information
            println("Step ", n, " of ", N, " completed.")
            println("Maximum link dimension of MPS: ", maxlinkdim(ψ))
            println("Relative Energy: ", real(inner(ψ', H, ψ) - inner(ψ₀', H, ψ₀)))
            println("Estimated time remained: ", timeElapsed * (N - n), " seconds.")
        end

        # heatmap(positions, times, En)
        # savefig("simulation_output/En.png")

        # heatmap(positions, times, EnLog)
        # savefig("simulation_output/EnLog.png")

        # heatmap(positions, times, Linkdims)
        # savefig("simulation_output/Linkdims.png")

        titleString = "L = $L, λ = $λ, jⁱᵐᵖ = $jⁱᵐᵖ, jᵂᴾ = $jᵂᴾ, kᵂᴾ = $kᵂᴾ,
                        σᵂᴾ = $σᵂᴾ, ϵᵂᴾ = $ϵᵂᴾ, Nᴰᴹᴿᴳ = $Nᴰᴹᴿᴳ, χᴰᴹᴿᴳ ≤ $χᴰᴹᴿᴳ,
                        ϵᴰᴹᴿᴳ = $ϵᴰᴹᴿᴳ, χᵀᴰⱽᴾ ≤ $χᵀᴰⱽᴾ, ϵᵀᴰⱽᴾ = $ϵᵀᴰⱽᴾ, dt = $dt, Δt = $Δt"

        l = @layout [grid(2,2); b{0.2h}]

        title = plot(title = titleString, grid = false, showaxis = false)
        p1 = heatmap(positions, times, En, title = "Energy density")
        p2 = heatmap(positions, times, EnLog, title = "Log10 of Energy dens.")
        p3 = heatmap(positions, times, Linkdims, title = "Link dims")
        p4 = plot(times, Maxlinkdim, title = "Max link dims")
        plot(p1, p2, p3, p4, title, layout = l, dpi = 300, size=(1000, 1000))
        savefig("simulation_output/All.png")
    end

    tdvp_time_evolution(H, ψ, ψ₀, dt, Δt)
end
