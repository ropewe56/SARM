σ  = 5.670374419e-8
SC = 1367.0

function DT(a, ΔF, T)
    # from sun = out travelling
    F0 = SC/4*(1-a)
    local Te
    if T  <= 0.0
        # effective radiation temperature
        Te  = (F0/σ)^(1/4)
    else
        Te = T
    end
    # upward from earth surface
    F0  = σ*(Te)^4
    # necessary ΔT to compensate ΔF
    ΔT = ((F0+ΔF)/σ)^(1/4) - Te

    ΔT2 = ((F0 + ΔF)/(σ*Te^4) - 1) * Te/4.0

    SC/4*(1-a),F0,T,ΔT,ΔT2
end

# (F0+DF)/ (σ * T^4) T/4 = (1 + 4 D/T)
# T^4 (1 + 4 D/T)

F0 = SC/4*(1-a)
Fs  = σ*(289)^4
Fs-F0
a  = 0.3
ΔF = 5.0
T  = -289.0 
DT(a, ΔF, T)