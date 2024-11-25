# function BL_eqn(  x₁, x₂, x₃, h₁, h₂, h₃, vars )
using DifferentialEquations;
using Symbolics;
using Latexify;
using SymbolicUtils;

@variables  t x₁ x₂ x₃ h₁(x₁,x₃) h₃(x₁,x₃) U₁(x₁,x₂,x₃) U₂(x₂) U₃(x₁,x₂,x₃) ρ(x₁,x₂,x₃) P(x₁,x₃) T(x₁,x₂,x₃);
@variables  κ(T), μ(x₂), Cₚ Rgas λ Re γ Mach Pᵣ; 
    


h = [h₁ 1 h₃];
x = [x₁ x₂ x₃];
U = [U₁ U₂ U₃];

Dx₁ = Differential(x₁);
Dx₂ = Differential(x₂);
Dx₃ = Differential(x₃);
Dt = Differential(t);
DT = Differential(T);


vars_base = [ρ, U₁, U₂, U₃, T, P, κ, μ];

Eqns_base = NS_eqn_for_BL( x₁, x₂, x₃, h₁, 1, h₃, vars_base );

Ceqn = (Symbolics.value(( ( Differential(x[1])(h[2]*h[3]*(ρ*U[1])) + Differential(x[2])(h[1]*h[3]*(ρ*U[2])) + Differential(x[3])(h[1]*h[2]*(ρ*U[3])) ) ) ));
Meqn =  ρ*U₁/h₁*Differential(x₁)(U₁) + ρ*U₂*Differential(x₂)(U₁) + + ρ*U₃/h₃*Differential(x₃)(U₁) + ρ*U₁*U₃*( Differential(x₃)(h₁)/(h₁*h₃) )
Eneqn = Symbolics.value(expand_derivatives(ET1)) - Symbolics.value(expand_derivatives(ET2)) - Symbolics.value(expand_derivatives(ET3)) - Symbolics.value(expand_derivatives(ET4));


tc1,tc2 = Symbolics.arguments(Symbolics.value(simplify(Eqns_base[2])));
tcc1 =  Symbolics.arguments(Symbolics.value(simplify(tc1)));

Meqn2 = Symbolics.diff2term((((tcc1./tc2))));

K = Symbolics.coeff((sum(Meqn2)), Differential(x₁)(h₃) )

latexify((:Ceqn_tr~(( Symbolics.value( ( expand_derivatives( K ) ))))),render=true) |> print
