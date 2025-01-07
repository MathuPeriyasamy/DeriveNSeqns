# function NS_eqn_v2(  r, θ, z, h1, h2, h3, vars )
using DifferentialEquations;
using Symbolics;
using Latexify;
using SymbolicUtils;
# @variables t r θ z h1(r,θ,z) h2(r,θ,z) h3(r,θ,z) U1(r,θ,z) U2(r,θ,z) U3(r,θ,z) ρ(r,θ,z) P(r,θ,z) T(r,θ,z) P(r,θ,z) Cp Rgas κ(T) μ(T) λ Re γ Mach Pr ϵ;
# @variables δU1(r,θ,z)
# h = [h1 h2 h3];
# x = [r θ z];
# U = [U1 U2 U3];

## Cylindrical CS
@variables t r θ z h1(r,θ,z) h2(r,θ,z) h3(r,θ,z) U1(r,θ,z) U2(r,θ,z) U3(r,θ,z) Ur(r,θ,z) Uth(r,θ,z) Uz(r,θ,z) ρ(r,θ,z) T(r,θ,z) P(r,θ,z) Cp Rgas κ(T) μ λ Re γ Mach Pr ϵ;
@variables δU1(r,θ,z)
h = [1 r 1];
x = [r θ z];
U = [U1 U2 U3];
h1 = 1; h2 = r; h3 = 1;
# U1 = Ur; U2 = Uth; U3 = Uz;
## The derivative variables
@variables h1r h2r h3r h1θ h2θ h3θ h1z h2z h3z;
@variables h1rr h2rr h3rr h1θθ h2θθ h3θθ h1zz h2zz h3zz;
@variables h1rθ h2rθ h3rθ h1θz h2θz h3θz h1rz h2rz h3rz;
# @variables h1θr h2θr h3θr h1zθ h2zθ h3zθ h1zr h2zr h3zr;

@variables U1r U2r U3r U1θ U2θ U3θ U1z U2z U3z;
@variables U1rr U2rr U3rr U1θθ U2θθ U3θθ U1zz U2zz U3zz;
@variables U1rθ U2rθ U3rθ U1θz U2θz U3θz U1rz U2rz U3rz;
# @variables U1θr U2θr U3θr U1zθ U2zθ U3zθ U1zr U2zr U3zr;

@variables Tr Tθ Tz;
@variables Trr Tθθ Tzz;
@variables Trθ Tθz Trz;
# @variables Tθr Tzθ Tzr;

@variables ρr ρθ ρz;
@variables ρrr ρθθ ρzz;
@variables ρrθ ρθz ρrz;
# @variables ρθr ρzθ ρzr;

@variables μ_T κ_T μ_TT κ_TT;

# Dr = Differential(r);
# Dθ = Differential(θ);
# Dz = Differential(z);

Dr = Differential(r);
Dθ = Differential(θ);
Dz = Differential(z);
Dt = Differential(t);
DT = Differential(T);
Drr = Dr^2; Dθθ = Dθ^2; Dzz = Dz^2;
Drθ = Dr*Dθ; Dθz = Dθ*Dz; Drz = Dr*Dz;
Dθr = Dθ*Dr; Dzθ = Dz*Dθ; Dzr = Dz*Dr;
DTT = DT^2;

# P = ρ*T/(γ*Mach^2);

T1 = ρ*[    ( (( U1/h1*Dr(U1) + U2/h2*Dθ(U1) + U3/h3*Dz(U1) + U1*U2/(h1*h2)*Dθ(h1) + U1*U3/(h1*h3)*Dz(h1) - U2^2/(h1*h2)*Dr(h2) - U3^2/(h1*h3)*Dr(h3) + Dt(U1) )) );
            ( (( U1/h1*Dr(U2) + U2/h2*Dθ(U2) + U3/h3*Dz(U2) + U1*U2/(h1*h2)*Dr(h2) + U2*U3/(h2*h3)*Dz(h2) - U1^2/(h1*h2)*Dθ(h1) - U3^2/(h2*h3)*Dθ(h3) + Dt(U2) )) );
            ( (( U1/h1*Dr(U3) + U2/h2*Dθ(U3) + U3/h3*Dz(U3) + U2*U3/(h2*h3)*Dθ(h3) + U1*U3/(h1*h3)*Dr(h3) - U2^2/(h2*h3)*Dz(h2) - U1^2/(h1*h3)*Dz(h1) + Dt(U3) )) )];


err = ( expand_derivatives(( 1/h1*Dr(U[1]) + U[2]/(h1*h2)*Dθ(h1) + U[3]/(h1*h3)*Dz(h1) )));
eθθ = ( expand_derivatives(( 1/h2*Dθ(U[2]) + U[1]/(h1*h2)*Dr(h2) + U[3]/(h2*h3)*Dz(h2) )));
ezz = ( expand_derivatives(( 1/h3*Dz(U[3]) + U[1]/(h1*h3)*Dr(h3) + U[2]/(h2*h3)*Dθ(h3) )));

erθ = ( expand_derivatives(( h[2]/h[1]*Dr(U[2]/h[2]) + h[1]/h[2]*Dθ(U[1]/h[1]) )));
eθz = ( expand_derivatives(( h[3]/h[2]*Dθ(U[3]/h[3]) + h[2]/h[3]*Dz(U[2]/h[2]) )));
erz = ( expand_derivatives(( h[3]/h[1]*Dr(U[3]/h[3]) + h[1]/h[3]*Dz(U[1]/h[1]) )));

Πrr = ( expand_derivatives(( 2//3*μ/Re*( 2*err -eθθ -ezz ) )));
Πθθ = ( expand_derivatives(( 2//3*μ/Re*( 2*eθθ -err -ezz ) )));
Πzz = ( expand_derivatives(( 2//3*μ/Re*( 2*ezz -err -eθθ ) )));

Πrθ = ( expand_derivatives(( μ/Re*erθ )));
Πθz = ( expand_derivatives(( μ/Re*eθz )));
Πrz = ( expand_derivatives(( μ/Re*erz )));

T2 = [  (( ( 1/(prod(h)).*( Dr(h2*h3*Πrr) + Dθ(h1*h3*Πrθ) + Dz(h1*h2*Πrz) ) ) )) + (( ( Πrθ./(h1*h2)*Dθ(h1) + Πrz./(h1*h3)*Dz(h1) - Πθθ./(h1*h2)*Dr(h2) - Πzz./(h1*h3)*Dr(h3) ) )) - 1/h1*Dr(P) ;
        (( ( 1/(prod(h)).*( Dr(h2*h3*Πrθ) + Dθ(h1*h3*Πθθ) + Dz(h1*h2*Πθz) ) ) )) + (( ( Πrθ./(h1*h2)*Dr(h2) + Πθz./(h2*h3)*Dz(h2) - Πrr./(h1*h2)*Dθ(h1) - Πzz./(h2*h3)*Dθ(h3) ) )) - 1/h1*Dθ(P) ;
        (( ( 1/(prod(h)).*( Dr(h2*h3*Πrz) + Dθ(h1*h3*Πθz) + Dz(h1*h2*Πzz) ) ) )) + (( ( Πrz./(h1*h3)*Dr(h3) + Πθz./(h2*h3)*Dθ(h3) - Πrr./(h1*h3)*Dz(h1) - Πθθ./(h2*h3)*Dz(h2) ) )) - 1/h1*Dz(P) ];

Cp = 1; # non dimensional

ET1 = (( ρ*Cp*( Dt(T) + U[1]/h[1]*Dr(T) + U[2]/h[2]*Dθ(T) + U[3]/h[3]*Dz(T)  )  ));

ET2 = (( (γ-1)*Mach^2*( Dt(P) + U[1]/h[1]*Dr(P) + U[2]/h[2]*Dθ(P) + U[3]/h[3]*Dz(P)  )  )); # non dimensional

ET3 =  (( 1/(Pr*Re)*( 1/(prod(h))*( Dr(κ*h[2]*h[3]/h[1]*Dr(T)) + Dθ(κ*h[1]*h[3]/h[2]*Dθ(T)) + Dz(κ*h[1]*h[2]/h[3]*Dz(T)) ) )  )); # non dimensional

ET4 = (( (γ-1)*Mach^2/Re*μ*( 2*( err.^2 + eθθ.^2 + ezz.^2 ) + erθ.^2 + eθz.^2 + erz.^2 -2//3*( err + eθθ + ezz ).^2 ) )); # non dimensional

Ceqn = ((( Dt(ρ) + ( 1/(prod(h))*( Differential(x[1])(h[2]*h[3]*(ρ*U[1])) + Differential(x[2])(h[1]*h[3]*(ρ*U[2])) + Differential(x[3])(h[1]*h[2]*(ρ*U[3])) )) ) ));
Meqn = ( (((T1))) - (((T2))) ); #Non dimensional form
Eneqn = ((ET1)) - ((ET2)) - ((ET3)) - ((ET4));

llist = [ Dr(h1), Dr(h2), Dr(h3), Dθ(h1), Dθ(h2), Dθ(h3), Dz(h1), Dz(h2), Dz(h3),
          Drr(h1), Drr(h2), Drr(h3), Dθθ(h1), Dθθ(h2), Dθθ(h3), Dzz(h1), Dzz(h2), Dzz(h3),
          Drθ(h1), Drθ(h2), Drθ(h3), Dθz(h1), Dθz(h2), Dθz(h3), Drz(h1), Drz(h2), Drz(h3),
          Dθr(h1), Dθr(h2), Dθr(h3), Dzθ(h1), Dzθ(h2), Dzθ(h3), Dzr(h1), Dzr(h2), Dzr(h3),
          Dr(U1), Dr(U2), Dr(U3), Dθ(U1), Dθ(U2), Dθ(U3), Dz(U1), Dz(U2), Dz(U3),
          Drr(U1), Drr(U2), Drr(U3), Dθθ(U1), Dθθ(U2), Dθθ(U3), Dzz(U1), Dzz(U2), Dzz(U3),
          Drθ(U1), Drθ(U2), Drθ(U3), Dθz(U1), Dθz(U2), Dθz(U3), Drz(U1), Drz(U2), Drz(U3),
          Dθr(U1), Dθr(U2), Dθr(U3), Dzθ(U1), Dzθ(U2), Dzθ(U3), Dzr(U1), Dzr(U2), Dzr(U3),
          Dr(T), Dθ(T), Dz(T), Dr(ρ), Dθ(ρ), Dz(ρ),
          Drr(T), Dθθ(T), Dzz(T), Drr(ρ), Dθθ(ρ), Dzz(ρ),
          Drθ(T), Dθz(T), Drz(T), Drθ(ρ), Dθz(ρ), Drz(ρ),
          Dθr(T), Dzθ(T), Dzr(T), Dθr(ρ), Dzθ(ρ), Dzr(ρ),
          DT(μ), DT(κ), DTT(μ), DTT(κ) ];
 
rlist = [ h1r, h2r, h3r, h1θ, h2θ, h3θ, h1z, h2z, h3z,
          h1rr, h2rr, h3rr, h1θθ, h2θθ, h3θθ, h1zz, h2zz, h3zz,
          h1rθ, h2rθ, h3rθ, h1θz, h2θz, h3θz, h1rz, h2rz, h3rz,
          h1rθ, h2rθ, h3rθ, h1θz, h2θz, h3θz, h1rz, h2rz, h3rz,
          U1r, U2r, U3r, U1θ, U2θ, U3θ, U1z, U2z, U3z,
          U1rr, U2rr, U3rr, U1θθ, U2θθ, U3θθ, U1zz, U2zz, U3zz,
          U1rθ, U2rθ, U3rθ, U1θz, U2θz, U3θz, U1rz, U2rz, U3rz,
          U1rθ, U2rθ, U3rθ, U1θz, U2θz, U3θz, U1rz, U2rz, U3rz,
          Tr, Tθ, Tz, ρr, ρθ, ρz,
          Trr, Tθθ, Tzz, ρrr, ρθθ, ρzz,
          Trθ, Tθz, Trz, ρrθ, ρθz, ρrz,
          Trθ, Tθz, Trz, ρrθ, ρθz, ρrz,
          μ_T, κ_T, μ_TT, κ_TT ];

Eqns_base = [ expand_derivatives(Ceqn);expand_derivatives(Meqn[1]);expand_derivatives(Meqn[2]);expand_derivatives(Meqn[3]);expand_derivatives(Eneqn) ];
Eqns_base = substitute(( (Eqns_base)),Dict( llist.=>rlist ));
open("Test.txt","w") do file    
    println(file,"eqns = ",replace(string(simplify(Symbolics.coeff( Eqns_base[4],μ )) ),"μ(T)"=>"μ","(r, θ, z)"=>"") );    
end         
eqn = replace(string(simplify(Symbolics.coeff( Eqns_base[5],μ )) ),"μ(T)"=>"μ","(r, θ, z)"=>"");
open("Test_latexify.txt","w") do file    
    println(file,latexify(eqn) );    
end   

