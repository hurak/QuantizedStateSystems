using QuantizedStateSystems
using LinearAlgebra
using Test

@testset "QSS1 for a scalar systems" begin
    @test begin
        a = -1.0
        f(x) = a*x
        x₀ = 10.0
        tspan = (0.0,1.0)
        Δq = 0.01
        tsol,xsol = qss1(f,x₀,tspan,Δq)
        xafun(t) = x₀*exp(a*t)
        xasol = xafun.(tsol)
        relerrsol = (xsol-xasol)./xasol
        norm(relerrsol)<0.01
    end
end
