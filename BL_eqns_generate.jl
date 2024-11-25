using DifferentialEquations;
using Symbolics;
using Latexify;
using SymbolicUtils;
include("NS_eqn_for_BL.jl");

@variables t ξ η ζ h₁(ξ,η,ζ) h₂(ξ,η,ζ) h₃(ξ,η,ζ) U₁(ξ,η,ζ) U₂(ξ,η,ζ) U₃(ξ,η,ζ) ρ(ξ,η,ζ) ρref Lref Uref Tref κref μref;
@variables T(ξ,η,ζ) P(ξ,η,ζ) Cₚ Rgas κ(T) μ(T) λ Re γ Mach Pᵣ; 

@variables h₁ξ h₂ξ h₃ξ h₁η h₂η h₃η h₁ζ h₂ζ h₃ζ;
@variables h₁ξξ h₂ξξ h₃ξξ h₁ηη h₂ηη h₃ηη h₁ζζ h₂ζζ h₃ζζ;
@variables h₁ξη h₂ξη h₃ξη h₁ηζ h₂ηζ h₃ηζ h₁ξζ h₂ξζ h₃ξζ;

@variables U₁ξ U₂ξ U₃ξ U₁η U₂η U₃η U₁ζ U₂ζ U₃ζ;
@variables U₁ξξ U₂ξξ U₃ξξ U₁ηη U₂ηη U₃ηη U₁ζζ U₂ζζ U₃ζζ;
@variables U₁ξη U₂ξη U₃ξη U₁ηζ U₂ηζ U₃ηζ U₁ξζ U₂ξζ U₃ξζ;

@variables Tξ Tη Tζ;
@variables Tξξ Tηη Tζζ;
@variables Tξη Tηζ Tξζ;

@variables ρξ ρη ρζ;
@variables ρξξ ρηη ρζζ;
@variables ρξη ρηζ ρξζ;

@variables μ_T κ_T μ_TT κ_TT;

# h₁ = 1; h₂ = 1; h₃ = 1;
# h = [h₁ h₂ h₃];
# x = [ξ η ζ];
# U = [U₁ U₂ U₃];

Dξ = Differential(ξ);
Dη = Differential(η);
Dζ = Differential(ζ);
Dt = Differential(t);
DT = Differential(T);
Dξξ = Dξ^2; Dηη = Dη^2; Dζζ = Dζ^2;
Dξη = Dξ*Dη; Dηζ = Dη*Dζ; Dξζ = Dξ*Dζ;
Dηξ = Dη*Dξ; Dζη = Dζ*Dη; Dζξ = Dζ*Dξ;
DTT = DT^2;



@variables δu(ξ,η,ζ,t) δv(ξ,η,ζ,t) δw(ξ,η,ζ,t) δT(ξ,η,ζ,t) δρ(ξ,η,ζ,t) δP(ξ,η,ζ,t) ϵ i ω;
A = Symbolics.variables(:A,1:5,1:5);
B = Symbolics.variables(:B,1:5,1:5);
C = Symbolics.variables(:C,1:5,1:5);
D = Symbolics.variables(:D,1:5,1:5);
E = Symbolics.variables(:E,1:5,1:5);
F = Symbolics.variables(:F,1:5,1:5);
G = Symbolics.variables(:G,1:5,1:5);
H = Symbolics.variables(:H,1:5,1:5);
I = Symbolics.variables(:I,1:5,1:5);
J = Symbolics.variables(:J,1:5,1:5);
K = Symbolics.variables(:K,1:5,1:5);

vars_base = [ρ, U₁, U₂, U₃, T, P, κ, μ];

Eqns_base = NS_eqn_for_BL( ξ, η, ζ, h₁, h₂, h₃, vars_base );

@variables r θ z;
@variables ρcyl(r,θ,z), U₁cyl(r,θ,z), U₂cyl(r,θ,z), U₃cyl(r,θ,z), T(r,θ,z), Pcyl(r,θ,z), κcyl(r,θ,z), μcyl(r,θ,z);

vars_base_cyl = [ρcyl, U₁cyl, U₂cyl, U₃cyl, T, Pcyl, κcyl, μcyl];
Eqns_base_cyl = NS_eqn_for_BL( r, θ, z, 1, r, 1, vars_base_cyl );


##### Non dimensionalization

llist = [ ρ, U₁, U₂, U₃, T, P, κ, μ, 
          Dξ(h₁), Dξ(h₂), Dξ(h₃), Dη(h₁), Dη(h₂), Dη(h₃), Dζ(h₁), Dζ(h₂), Dζ(h₃),
          Dξξ(h₁), Dξξ(h₂), Dξξ(h₃), Dηη(h₁), Dηη(h₂), Dηη(h₃), Dζζ(h₁), Dζζ(h₂), Dζζ(h₃),
          Dξη(h₁), Dξη(h₂), Dξη(h₃), Dηζ(h₁), Dηζ(h₂), Dηζ(h₃), Dξζ(h₁), Dξζ(h₂), Dξζ(h₃),
          Dηξ(h₁), Dηξ(h₂), Dηξ(h₃), Dζη(h₁), Dζη(h₂), Dζη(h₃), Dζξ(h₁), Dζξ(h₂), Dζξ(h₃),
          Dξ(U₁), Dξ(U₂), Dξ(U₃), Dη(U₁), Dη(U₂), Dη(U₃), Dζ(U₁), Dζ(U₂), Dζ(U₃),
          Dξξ(U₁), Dξξ(U₂), Dξξ(U₃), Dηη(U₁), Dηη(U₂), Dηη(U₃), Dζζ(U₁), Dζζ(U₂), Dζζ(U₃),
          Dξη(U₁), Dξη(U₂), Dξη(U₃), Dηζ(U₁), Dηζ(U₂), Dηζ(U₃), Dξζ(U₁), Dξζ(U₂), Dξζ(U₃),
          Dηξ(U₁), Dηξ(U₂), Dηξ(U₃), Dζη(U₁), Dζη(U₂), Dζη(U₃), Dζξ(U₁), Dζξ(U₂), Dζξ(U₃),
          Dξ(T), Dη(T), Dζ(T), Dξ(ρ), Dη(ρ), Dζ(ρ),
          Dξξ(T), Dηη(T), Dζζ(T), Dξξ(ρ), Dηη(ρ), Dζζ(ρ),
          Dξη(T), Dηζ(T), Dξζ(T), Dξη(ρ), Dηζ(ρ), Dξζ(ρ),
          Dηξ(T), Dζη(T), Dζξ(T), Dηξ(ρ), Dζη(ρ), Dζξ(ρ),
          DT(μ), DT(κ), DTT(μ), DTT(κ) ];


rlist2 = [ ρ*ρref, U₁*Uref, U₂*Uref/sqrt(Re), U₃*Uref, T*Tref, P*ρref*Uref^2, κ*κref, μ*μref,
          h₁ξ/Lref, h₂ξ/Lref, h₃ξ/Lref, h₁η/Lref*sqrt( Re ), h₂η/Lref*sqrt( Re ), h₃η/Lref*sqrt( Re ), h₁ζ/Lref, h₂ζ/Lref, h₃ζ/Lref,
          h₁ξξ/Lref^2, h₂ξξ/Lref^2, h₃ξξ/Lref^2, h₁ηη/Lref^2*( Re ), h₂ηη/Lref^2*( Re ), h₃ηη/Lref^2*( Re ), h₁ζζ/Lref^2, h₂ζζ/Lref^2, h₃ζζ/Lref^2,
          h₁ξη/Lref^2*sqrt( Re ), h₂ξη/Lref^2*sqrt( Re ), h₃ξη/Lref^2*sqrt( Re ), h₁ηζ/Lref^2*sqrt( Re ), h₂ηζ/Lref^2*sqrt( Re ), h₃ηζ/Lref^2*sqrt( Re ), h₁ξζ/Lref^2, h₂ξζ/Lref^2, h₃ξζ/Lref^2,
          h₁ξη/Lref^2*sqrt( Re ), h₂ξη/Lref^2*sqrt( Re ), h₃ξη/Lref^2*sqrt( Re ), h₁ηζ/Lref^2*sqrt( Re ), h₂ηζ/Lref^2*sqrt( Re ), h₃ηζ/Lref^2*sqrt( Re ), h₁ξζ/Lref^2, h₂ξζ/Lref^2, h₃ξζ/Lref^2,
          U₁ξ*Uref/Lref, U₂ξ*Uref/sqrt(Re)/Lref, U₃ξ*Uref/Lref, U₁η*Uref/Lref*sqrt( Re ), U₂η*Uref/sqrt(Re)/Lref*sqrt( Re ), U₃η*Uref/Lref*sqrt( Re ), U₁ζ*Uref/Lref, U₂ζ*Uref/sqrt(Re)/Lref, U₃ζ*Uref/Lref,
          U₁ξξ*Uref/Lref^2, U₂ξξ*Uref/sqrt(Re)/Lref^2, U₃ξξ*Uref/Lref^2, U₁ηη*Uref/Lref^2*( Re ), U₂ηη*Uref/sqrt(Re)/sqrt(Re)/Lref^2*( Re ), U₃ηη*Uref/Lref^2*( Re ), U₁ζζ*Uref/Lref^2, U₂ζζ*Uref/sqrt(Re)/Lref^2, U₃ζζ*Uref/Lref^2,
          U₁ξη*Uref/Lref^2*sqrt( Re ), U₂ξη*Uref/sqrt(Re)/Lref^2*sqrt( Re ), U₃ξη*Uref/Lref^2*sqrt( Re ), U₁ηζ*Uref/Lref^2*sqrt( Re ), U₂ηζ*Uref/sqrt(Re)/Lref^2*sqrt( Re ), U₃ηζ*Uref/Lref^2*sqrt( Re ), U₁ξζ*Uref/Lref^2, U₂ξζ*Uref/sqrt(Re)/Lref^2, U₃ξζ*Uref/Lref^2,
          U₁ξη*Uref/Lref^2*sqrt( Re ), U₂ξη*Uref/sqrt(Re)/Lref^2*sqrt( Re ), U₃ξη*Uref/Lref^2*sqrt( Re ), U₁ηζ*Uref/Lref^2*sqrt( Re ), U₂ηζ*Uref/sqrt(Re)/Lref^2*sqrt( Re ), U₃ηζ*Uref/Lref^2*sqrt( Re ), U₁ξζ*Uref/Lref^2, U₂ξζ*Uref/sqrt(Re)/Lref^2, U₃ξζ*Uref/Lref^2,
          Tξ*Tref/Lref, Tη*Tref/Lref*sqrt( Re ), Tζ*Tref/Lref, ρξ*ρref/Lref, ρη*ρref/Lref*sqrt( Re ), ρζ*ρref/Lref,
          Tξξ*Tref/Lref^2, Tηη*Tref/Lref^2*( Re ), Tζζ*Tref/Lref^2, ρξξ*ρref/Lref^2, ρηη*ρref/Lref^2*( Re ), ρζζ*ρref/Lref^2,
          Tξη*Tref/Lref^2*sqrt( Re ), Tηζ*Tref/Lref^2*sqrt( Re ), Tξζ*Tref/Lref^2, ρξη*ρref/Lref^2*sqrt( Re ), ρηζ*ρref/Lref^2*sqrt( Re ), ρξζ*ρref/Lref^2,
          Tξη*Tref/Lref^2*sqrt( Re ), Tηζ*Tref/Lref^2*sqrt( Re ), Tξζ*Tref/Lref^2, ρξη*ρref/Lref^2*sqrt( Re ), ρηζ*ρref/Lref^2*sqrt( Re ), ρξζ*ρref/Lref^2,
          μ_T*μref/Tref, κ_T*κref/Tref, μ_TT*μref/Tref^2, κ_TT*κref/Tref^2 ];



# vars_base_non = [ρ*ρref, U₁*Uref, U₂*Uref/sqrt(Re), U₃*Uref, T*Tref, P*ρref*Uref^2, κ*κref, μ*μref];


Eqns_base = simplify(expand_derivatives(Eqns_base));


# Eqns_base = substitute(( Symbolics.value(simplify(expand_derivatives(Eqns_base)))),Dict( llist.=>rlist2 ));
# Eqns_base = substitute(( Symbolics.value(simplify(expand_derivatives(Eqns_base)))),Dict( rlist.=>rlist2 ));
# Eqns_base = substitute(( Symbolics.value(simplify(expand_derivatives(Eqns_base)))),Dict( vars_base.=>vars_base_non ));
Eqns_base = substitute(( Symbolics.value(simplify(expand_derivatives(Eqns_base)))),Dict( llist.=>rlist2 ));


tc1 = Symbolics.arguments(Symbolics.value((Eqns_base_cyl[2])));

tcc1,tcc2 = Symbolics.arguments(Symbolics.value(simplify(tc1)));
tccc1 = Symbolics.arguments(Symbolics.value(simplify(tcc1)));


ceqn = ( (((tccc1./(tcc2*tc2  ))))  );
latexify(:Ceqn_tr~(( ( ( Symbolics.value( Eqns_base_cyl[2] ) )))),render = true) |> print;


#### mom Eqn_n


tc1 = Symbolics.arguments(Symbolics.value(expand_derivatives(Eqns_base[2])));

tcc1,tcc2 = Symbolics.arguments(Symbolics.value(simplify(tc1[1])));
tccc1 = Symbolics.arguments(Symbolics.value(simplify(tcc1)));


meqn = ( ((sum(tccc1./(tcc2*(1)  ))))  );
latexify(:Ceqn_tr~(( ( ( Symbolics.value( ((meqn) ) ) )))),render = true) |> print;


# simplify(Differential(ζ)(ρ*ρref))
