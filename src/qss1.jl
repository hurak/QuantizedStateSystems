"""
    txarray,xarray,tqarray,qarray = qss1(f,x₀,tspan,Δq,tuarray,uarray)

Solve the initial value problem (IVP) for a given first-order ODE with inputs using the method of hysteretically quantized state system (QSS) or order 1.

For a single first-order explicit ODEs (aka state equation) `ẋ = f(x,u)`, with the initial value `x₀` specified at an initial time `t₀`, and possibly a piecewise (control) input `u` given by a pair of vectors `tuarray` and `uarray`of times and values, find the solution on the time span `tspan=(t₀,t₁)`. The only parameter for the method is the quantum `Δq`, while the hysteresis band has the same width as the quantum.

# Arguments
- `f`: function defining the right hand side of the differential equation.
- `x₀`: initial condition.
- `tspan`: time span.
- `Δq`: quantum.
- `tuarray`: array of times at which the the input changes to a new constant value.
- `uarray`: array of values of the piecewise constant input.

# Outputs
- `txarray`: a vector of times when the value of `x` is computed, between these values it evolves linearly.
- `xarray`: a vector of values of `x` corresponding to the times.
- `tqarray`: a vector of times when the value of piecewise constant `q` changes.
- `qarray`: a vector of values of the piecewise constant quantized `q`.

Note that for a first-order system with no external events the `x` and `q` arrays are identical up to the initial values (`x₀` before vs. after quantization).
"""
function qss1(f,x₀,tspan,Δq,tuarray,uarray)
    (t₀,t₁) = tspan             # Initial and final times for the simulation.
    if tuarray[1]==t₀           # Setting the initial value of the input:
        u = uarray[1]           # - either to the first value of the (piecewise) input,
        popfirst!(uarray)
        popfirst!(tuarray)
    else
        u = 0.0                 # - or to zero.
    end
    txarray = [t₀,]             # Initializing an array of times. The length not known apriori.
    xarray = [x₀,]              # Initializing an array of state variables.
    t = t₀                      # Initializing the current time.
    x = x₀                      # Initializing the current state.
    q = Δq*floor(x₀/Δq)         # Setting the quantized initial state.
    tqarray = [t,]              # Initializing an array of times for the quantized state variable.
    qarray = [q,]               # Initializing an array of (samples of) the quant. state variable.
    while t < t₁
        k = f(q,u)              # Derivative corresponding to the quantized state.
        Δt = abs(Δq/k)          # Time increment based only on f and not u.
        if isempty(tuarray) && t+Δt > t₁ # Predicted crossing of the quant. level beyond time span.
            Δt = t₁-t           # The time from now till the end.
            t = t₁              # The last time x is evaluated is at the very end of the time span.
            x += k*Δt           # The state x updated at the very end of the timespan.
            push!(txarray,t)    # Store the time when x experiences a change in derivative and
            push!(xarray,x)     # the actual value of x when that happens. No change to q and its t.
        elseif (~isempty(tuarray) && t+Δt < tuarray[1]) || (isempty(tuarray) && t+Δt <= t₁)   # If the next crossing of quant. level is before the change of u or before the end.
            t += Δt             # Set the next time.
            q += sign(k)*Δq     # Set the next quantized state.
            x = q               # State variable x is following q at the moment q changes.
            push!(txarray,t)    # Append the time to the array of times.
            push!(xarray,q)     # Append the quantized state to the array of (samples of) states.
            push!(tqarray,t)
            push!(qarray,q)     # The same for the quantized states (up to the first value).
        elseif ~isempty(tuarray) && t+Δt >= tuarray[1] # If the change in u arrives before x crosses the quant. level.
            Δt = tuarray[1]-t   # The time from now till the moment u changes.
            t = tuarray[1]      # The moment when u changes.
            u = uarray[1]       # The new value of u.
            popfirst!(tuarray)  # Discard the time so that it is not used in the next iteration.
            popfirst!(uarray)   # Discard also the value of control for the same reason.
            x += k*Δt           # The state x gets updated but q does not.
            push!(txarray,t)    # Store the time when x experiences a change in derivative and
            push!(xarray,x)     # the actual value of x when that happens. No change to q and its t.
        end
    end
    return (txarray,xarray,tqarray,qarray)
end
