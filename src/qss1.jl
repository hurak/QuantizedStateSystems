"""
    t,x = qss1(f,x₀,tspan,Δq)

Solve the initial value problem (IVP) for a given first-order ODE using the method of hysteretically quantized state system (QSS) or order 1.

Consider a single first-order explicit ODEs (aka state equation) `ẋ = f(x)`, with the initial value `x₀` specified at an initial time `t₀`. Solution is to be found on a time span `tspan=(t₀,t₁)`. The only parameter is the quantum `Δq`. The hysteresis band has the same width as the quantum.

The outputs are two arrays of arrays: the 1D array `t` of 1D arrays (vectors) of times and the correspoding 1D array `x` of the 1D arrays corresponding to the evolution of individual quantized state variables.

# References

1. E. Kofman and S. Junco, "Quantized-state systems: a DEVS Approach for continuous system simulation", Trans. Soc. Comput. Simul. Int., vol. 18, no. 3, pp. 123–132, Sep. 2001.
2. F. E. Cellier and E. Kofman, Continuous System Simulation, Springer, 2010.
3. B. P. Zeigler, A. Muzy, and E. Kofman, Theory of Modeling and Simulation: Discrete Event and Iterative System Computational Foundations, 3rd ed. Academic Press, 2018.

# Example

```julia
julia> f(x) = -1.0*x
f (generic function with 1 method)

julia> x₀ = 10.0
10.0

julia> tspan = (0.0,1.0)
(0.0, 1.0)

julia> Δq = 0.1
0.1

julia> T,X = qss1(f,x₀,tspan,Δq)
([0.0, 0.01, 0.020101010101010102, 0.030305091733663164, 0.040614370084178626, 0.05103103675084529, 0.06155735254031897, 0.07219565041265939, 0.08294833858470241, 0.0938179038020937  …  0.7924294020882965, 0.8146516243105186, 0.8373788970377912, 0.8606347109912795, 0.8844442348008033, 0.9088344787032422, 0.9338344787032421, 0.9594755043442676, 0.9857912938179517, 1.0128183208449786], [10.0, 9.9, 9.8, 9.700000000000001, 9.600000000000001, 9.500000000000002, 9.400000000000002, 9.300000000000002, 9.200000000000003, 9.100000000000003  …  4.5000000000000195, 4.40000000000002, 4.30000000000002, 4.200000000000021, 4.100000000000021, 4.000000000000021, 3.9000000000000212, 3.800000000000021, 3.700000000000021, 3.600000000000021])
```
"""
function qss1(f,x₀,tspan,Δq)
    (t₀,t₁) = tspan             # Initial and final times.
    tarray = [t₀,]              # Initializing an array of times. The length not known apriori.
    xarray = [x₀,]              # Initializing an array of state vectors.
    x = x₀                      # Initializing the state vector.
    t = t₀                      # Initializing the current time.
    q = Δq*floor(x/Δq)          # Quantized initial state.
    while t <= t₁
        k = f(q)                # Derivative corresponding to the quantized state.
        Δt = abs(Δq/k)          # Time increment.
        t += Δt                 # Next time.
        q += sign(k)*Δq         # Next quantized state.
        push!(tarray,t)         # Append the time to the array of times.
        push!(xarray,q)         # Append the quantized state to the array of (samples of) states.
    end
    return (tarray,xarray)
end
