using DifferentialEquations;
using Symbolics;
using Latexify;
using SymbolicUtils;

@variables t x₁ x₂ x₃ h₁(x₁,x₂) h₂(x₁,x₂) h₃(x₁,x₂) U₁(x₁,x₂,t) U₂(x₁,x₂,t) U₃(x₁,x₂,t) ρ(x₁,x₂,t); 
@variables T(x₁,x₂,t) P(x₁,t) Cₚ Rgas κ μ λ Re γ Mach Pᵣ; 

h₂ = 1;
h₃ = rₘ;
U₃= 0;

h = [h₁ h₂ h₃];
x = [x₁ x₂ x₃];
U = [U₁ U₂ U₃];

Dx₁ = Differential(x₁);
Dx₂ = Differential(x₂);
Dx₃ = Differential(x₃);

Dt = Differential(t);
Dx₁x₁ = Dx₁^2; Dx₂x₂ = Dx₂^2; Dx₃x₃ = Dx₃^2;
Dx₁x₂ = Dx₁*Dx₂; Dx₂x₃ = Dx₂*Dx₃; Dx₁x₃ = Dx₁*Dx₃;

### Functions to calculate gradient, divergence and curl with curvilinear parameters

function divergence_curv(f,x,h)
    q = 1/(prod(h))*( Differential(x[1])(h[2]*h[3]*f[1]) + Differential(x[2])(h[1]*h[3]*f[2]) + Differential(x[3])(h[1]*h[2]*f[3]) );
    return q;
end;

Ceqn =  Dt(ρ) + expand_derivatives( divergence_curv(ρ.*U,x,h) ) ;

# @variables  h₃;

κ = μ;

h = [h₁ h₂ h₃];

###### SUBSTITUTE DOESN'T WORK PROPERLY ###########

#=      
Momentum equation = ρ*( Dt(V) + V.ΔV ) - 1/Re*(ρ*f - Δ.Pᵢⱼ)
=#

# V_ΔV_inT2_x₁ = U₁/h₁*Dx₁(U₁) + U₂/h₂*Dx₂(U₁) + U₃/h₃*Dx₃(U₁) + U₁*U₂/(h₁*h₂)*Dx₂(h₁) + U₁*U₃/(h₁*h₃)*Dx₃(h₁) - U₂^2/(h₁*h₂)*Dx₁(h₂) - U₃^2/(h₁*h₃)*Dx₁(h₃);
# V_ΔV_inT2_x₂ = substitute(V_ΔV_inT2_x₁,(Dict(U₁=>U₂,h₁=>h₂)));
# V_ΔV_inT2_x₃ = substitute(V_ΔV_inT2_x₁,(Dict(U₁=>U₃,h₁=>h₃)));

T1 = ρ*[    simplify( Symbolics.diff2term(expand_derivatives( U₁/h₁*Dx₁(U₁) + U₂/h₂*Dx₂(U₁) + U₃/h₃*Dx₃(U₁) + U₁*U₂/(h₁*h₂)*Dx₂(h₁) + U₁*U₃/(h₁*h₃)*Dx₃(h₁) - U₂^2/(h₁*h₂)*Dx₁(h₂) - U₃^2/(h₁*h₃)*Dx₁(h₃) + Dt(U₁) )) );
            simplify( Symbolics.diff2term(expand_derivatives( U₁/h₁*Dx₁(U₂) + U₂/h₂*Dx₂(U₂) + U₃/h₃*Dx₃(U₂) + U₁*U₂/(h₁*h₂)*Dx₁(h₂) + U₂*U₃/(h₂*h₃)*Dx₃(h₂) - U₁^2/(h₁*h₂)*Dx₂(h₁) - U₃^2/(h₂*h₃)*Dx₂(h₃) + Dt(U₂) )) );
            simplify( Symbolics.diff2term(expand_derivatives( U₁/h₁*Dx₁(U₃) + U₂/h₂*Dx₂(U₃) + U₃/h₃*Dx₃(U₃) + U₂*U₃/(h₂*h₃)*Dx₂(h₃) + U₁*U₃/(h₁*h₃)*Dx₁(h₃) - U₂^2/(h₂*h₃)*Dx₃(h₂) - U₁^2/(h₁*h₃)*Dx₃(h₁) + Dt(U₃) )) )];


ex₁x₁ = simplify( Symbolics.diff2term(expand_derivatives( 1/h₁*Dx₁(U[1]) + U[2]/(h₁*h₂)*Dx₂(h₁) + U[3]/(h₁*h₃)*Dx₃(h₁) )));
ex₂x₂ = simplify( Symbolics.diff2term(expand_derivatives( 1/h₂*Dx₂(U[2]) + U[1]/(h₁*h₂)*Dx₁(h₂) + U[3]/(h₂*h₃)*Dx₃(h₂) )));
ex₃x₃ = simplify( Symbolics.diff2term(expand_derivatives( 1/h₃*Dx₃(U[3]) + U[1]/(h₁*h₃)*Dx₁(h₃) + U[2]/(h₂*h₃)*Dx₂(h₃) )));
# ex₂x₂ = substitute(ex₁x₁,(Dict(h₁=>h₂,h₂=>h₁,U₁=>U₂,U₂=>U₁,x₁=>x₂,x₂=>x₁)));
# ex₃x₃ = substitute(ex₁x₁,(Dict(h₁=>h₃,h₃=>h₁,U₁=>U₃,U₃=>U₁,x₁=>x₃,x₃=>x₁)));
ex₁x₂ = simplify( Symbolics.diff2term(expand_derivatives( h[2]/h[1]*Dx₁(U[2]/h[2]) + h[1]/h[2]*Dx₂(U[1]/h[1]) )));
ex₂x₃ = simplify( Symbolics.diff2term(expand_derivatives( h[3]/h[2]*Dx₂(U[3]/h[3]) + h[2]/h[3]*Dx₃(U[2]/h[2]) )));
ex₁x₃ = simplify( Symbolics.diff2term(expand_derivatives( h[3]/h[1]*Dx₁(U[3]/h[3]) + h[1]/h[3]*Dx₃(U[1]/h[1]) )));
# ex₁x₃ = substitute(ex₂x₃,(Dict(h₁=>h₂,h₂=>h₁,U₁=>U₂,U₂=>U₁,x₁=>x₂,x₂=>x₁)));
# ex₁x₂ = substitute(ex₂x₃,(Dict(h₂=>h₁,h₃=>h₂,U₂=>U₁,U₃=>U₂,x₂=>x₁,x₃=>x₂)));

Πx₁ = simplify( Symbolics.diff2term(expand_derivatives( -P + 2/3*μ/Re*( 2*ex₁x₁ -ex₂x₂ -ex₃x₃ ) )));
Πx₂ = simplify( Symbolics.diff2term(expand_derivatives( -P + 2/3*μ/Re*( 2*ex₂x₂ -ex₁x₁ -ex₃x₃ ) )));
Πx₃ = simplify( Symbolics.diff2term(expand_derivatives( -P + 2/3*μ/Re*( 2*ex₃x₃ -ex₁x₁ -ex₂x₂ ) )));
Πx₁x₂ = simplify( Symbolics.diff2term(expand_derivatives( μ/Re*ex₁x₂ )));
Πx₂x₃ = simplify( Symbolics.diff2term(expand_derivatives( μ/Re*ex₂x₃ )));
Πx₁x₃ = simplify( Symbolics.diff2term(expand_derivatives( μ/Re*ex₁x₃ )));

# Δ_Πᵢⱼ_1 = ( 1/(prod(h)).*( Dx₁(h₂*h₃*Πx₁) + Dx₂(h₁*h₃*Πx₁x₂) + Dx₃(h₁*h₂*Πx₁x₃) ) ) + ( Πx₁x₂./(h₁*h₂)*Dx₂(h₁) + Πx₁x₃./(h₁*h₃)*Dx₃(h₁) - Πx₂./(h₁*h₂)*Dx₁(h₂) - Πx₃./(h₁*h₃)*Dx₁(h₃) );
# # Δ_Πᵢⱼ_2 = substitute(Δ_Πᵢⱼ_1,(Dict(x₁=>x₂,x₂=>x₁,h₁=>h₂,h₂=>h₁)));
# Δ_Πᵢⱼ_2 = ( 1/(prod(h)).*( Dx₁(h₂*h₃*Πx₁x₂) + Dx₂(h₁*h₃*Πx₂) + Dx₃(h₁*h₂*Πx₂x₃) ) ) + ( Πx₁x₂./(h₁*h₂)*Dx₁(h₂) + Πx₂x₃./(h₂*h₃)*Dx₃(h₂) - Πx₁./(h₁*h₂)*Dx₂(h₁) - Πx₃./(h₂*h₃)*Dx₂(h₃) );
# Δ_Πᵢⱼ_3 = ( 1/(prod(h)).*( Dx₁(h₂*h₃*Πx₁x₃) + Dx₂(h₁*h₃*Πx₂) + Dx₃(h₁*h₂*Πx₃) ) ) + ( Πx₁x₃./(h₁*h₃)*Dx₁(h₃) + Πx₂x₃./(h₂*h₃)*Dx₂(h₃) - Πx₁./(h₁*h₃)*Dx₃(h₁) - Πx₂./(h₂*h₃)*Dx₃(h₂) );

T2 = [  Symbolics.diff2term(expand_derivatives( ( 1/(prod(h)).*( Dx₁(h₂*h₃*Πx₁) + Dx₂(h₁*h₃*Πx₁x₂) + Dx₃(h₁*h₂*Πx₁x₃) ) ) )) + Symbolics.diff2term(expand_derivatives( ( Πx₁x₂./(h₁*h₂)*Dx₂(h₁) + Πx₁x₃./(h₁*h₃)*Dx₃(h₁) - Πx₂./(h₁*h₂)*Dx₁(h₂) - Πx₃./(h₁*h₃)*Dx₁(h₃) ) ));
        ( 1/(prod(h)).*( Dx₁(h₂*h₃*Πx₁x₂) + Dx₂(h₁*h₃*Πx₂) + Dx₃(h₁*h₂*Πx₂x₃) ) ) + ( Πx₁x₂./(h₁*h₂)*Dx₁(h₂) + Πx₂x₃./(h₂*h₃)*Dx₃(h₂) - Πx₁./(h₁*h₂)*Dx₂(h₁) - Πx₃./(h₂*h₃)*Dx₂(h₃) );
        ( 1/(prod(h)).*( Dx₁(h₂*h₃*Πx₁x₃) + Dx₂(h₁*h₃*Πx₂) + Dx₃(h₁*h₂*Πx₃) ) ) + ( Πx₁x₃./(h₁*h₃)*Dx₁(h₃) + Πx₂x₃./(h₂*h₃)*Dx₂(h₃) - Πx₁./(h₁*h₃)*Dx₃(h₁) - Πx₂./(h₂*h₃)*Dx₃(h₂) )];

Meqn = (Symbolics.diff2term(expand_derivatives(T1))) - (Symbolics.diff2term(expand_derivatives(T2))) ;


#= Energy equation
En.eqn = ρ*( Dt(T) + V.ΔT ) -( γ-1 )*Mach^2*( Dt(P) + V.ΔP ) - 1/(Pᵣ*Re)*Δ²T  - (γ-1)*Mach^2/Re*( ϕ )

=#


ET1 = ρ*( Dt(T) + U[1]/h[1]*Dx₁(T) + U[2]/h[2]*Dx₂(T) + U[3]/h[3]*Dx₃(T)  );

ET2 = (γ-1)*Mach^2*( Dt(P) + U[1]/h[1]*Dx₁(P) + U[2]/h[2]*Dx₂(P) + U[3]/h[3]*Dx₃(P)  );

ET3 = 1/(Pᵣ*Re)*1/(prod(h))*( Dx₁(κ*h[2]*h[3]/h[1]*Dx₁(T)) + Dx₂(κ*h[1]*h[3]/h[2]*Dx₂(T)) + Dx₃(κ*h[1]*h[2]/h[3]*Dx₃(T)) );

ET4 = (γ-1)*Mach^2/Re*μ*( 2*( ex₁x₁.^2 + ex₂x₂.^2 + ex₃x₃.^2 ) + ex₁x₂.^2 + ex₂x₃.^2 + ex₁x₃.^2 -2/3*( ex₁x₁ + ex₂x₂ + ex₃x₃ ).^2 );

En_eqn = ET1 - ET2 - ET3 - ET4;

print(latexify((:En_eqn~Symbolics.diff2term(expand_derivatives( Symbolics.value( ( En_eqn ) )))),render = true,cdot=true,adjustment=:c)) 