# function NS_eqn_v2(  x1, x2, x3, h1, h2, h3, vars )
using DifferentialEquations;
using Symbolics;
using Latexify;
using SymbolicUtils;
@variables t x1 x2 x3 h1(x1,x2,x3) h2(x1,x2,x3) h3(x1,x2,x3) U1(x1,x2,x3) U2(x1,x2,x3) U3(x1,x2,x3) ρ(x1,x2,x3) P(x1,x2,x3) T(x1,x2,x3) P(x1,x2,x3) Cp Rgas κ(T) μ(T) λ Re γ Mach Pr ϵ;

h = [h1 h2 h3];
x = [x1 x2 x3];
U = [U1 U2 U3];

## The derivative variables
@variables h1x1 h2x1 h3x1 h1x2 h2x2 h3x2 h1x3 h2x3 h3x3;
@variables h1x1x1 h2x1x1 h3x1x1 h1x2x2 h2x2x2 h3x2x2 h1x3x3 h2x3x3 h3x3x3;
@variables h1x1x2 h2x1x2 h3x1x2 h1x2x3 h2x2x3 h3x2x3 h1x1x3 h2x1x3 h3x1x3;
# @variables h1x2x1 h2x2x1 h3x2x1 h1x3x2 h2x3x2 h3x3x2 h1x3x1 h2x3x1 h3x3x1;

@variables U1x1 U2x1 U3x1 U1x2 U2x2 U3x2 U1x3 U2x3 U3x3;
@variables U1x1x1 U2x1x1 U3x1x1 U1x2x2 U2x2x2 U3x2x2 U1x3x3 U2x3x3 U3x3x3;
@variables U1x1x2 U2x1x2 U3x1x2 U1x2x3 U2x2x3 U3x2x3 U1x1x3 U2x1x3 U3x1x3;
# @variables U1x2x1 U2x2x1 U3x2x1 U1x3x2 U2x3x2 U3x3x2 U1x3x1 U2x3x1 U3x3x1;

@variables Tx1 Tx2 Tx3;
@variables Tx1x1 Tx2x2 Tx3x3;
@variables Tx1x2 Tx2x3 Tx1x3;
# @variables Tx2x1 Tx3x2 Tx3x1;

@variables ρx1 ρx2 ρx3;
@variables ρx1x1 ρx2x2 ρx3x3;
@variables ρx1x2 ρx2x3 ρx1x3;
# @variables ρx2x1 ρx3x2 ρx3x1;

@variables μ_T κ_T μ_TT κ_TT;

Dx1 = Differential(x1);
Dx2 = Differential(x2);
Dx3 = Differential(x3);
Dt = Differential(t);
DT = Differential(T);
Dx1x1 = Dx1^2; Dx2x2 = Dx2^2; Dx3x3 = Dx3^2;
Dx1x2 = Dx1*Dx2; Dx2x3 = Dx2*Dx3; Dx1x3 = Dx1*Dx3;
Dx2x1 = Dx2*Dx1; Dx3x2 = Dx3*Dx2; Dx3x1 = Dx3*Dx1;
DTT = DT^2;

# P = ρ*T/(γ*Mach^2);

T1 = ρ*[    ( (( U1/h1*Dx1(U1) + U2/h2*Dx2(U1) + U3/h3*Dx3(U1) + U1*U2/(h1*h2)*Dx2(h1) + U1*U3/(h1*h3)*Dx3(h1) - U2^2/(h1*h2)*Dx1(h2) - U3^2/(h1*h3)*Dx1(h3) + Dt(U1) )) );
            ( (( U1/h1*Dx1(U2) + U2/h2*Dx2(U2) + U3/h3*Dx3(U2) + U1*U2/(h1*h2)*Dx1(h2) + U2*U3/(h2*h3)*Dx3(h2) - U1^2/(h1*h2)*Dx2(h1) - U3^2/(h2*h3)*Dx2(h3) + Dt(U2) )) );
            ( (( U1/h1*Dx1(U3) + U2/h2*Dx2(U3) + U3/h3*Dx3(U3) + U2*U3/(h2*h3)*Dx2(h3) + U1*U3/(h1*h3)*Dx1(h3) - U2^2/(h2*h3)*Dx3(h2) - U1^2/(h1*h3)*Dx3(h1) + Dt(U3) )) )];


ex1x1 = ( (( 1/h1*Dx1(U[1]) + U[2]/(h1*h2)*Dx2(h1) + U[3]/(h1*h3)*Dx3(h1) )));
ex2x2 = ( (( 1/h2*Dx2(U[2]) + U[1]/(h1*h2)*Dx1(h2) + U[3]/(h2*h3)*Dx3(h2) )));
ex3x3 = ( (( 1/h3*Dx3(U[3]) + U[1]/(h1*h3)*Dx1(h3) + U[2]/(h2*h3)*Dx2(h3) )));

ex1x2 = ( (( h[2]/h[1]*Dx1(U[2]/h[2]) + h[1]/h[2]*Dx2(U[1]/h[1]) )));
ex2x3 = ( (( h[3]/h[2]*Dx2(U[3]/h[3]) + h[2]/h[3]*Dx3(U[2]/h[2]) )));
ex1x3 = ( (( h[3]/h[1]*Dx1(U[3]/h[3]) + h[1]/h[3]*Dx3(U[1]/h[1]) )));

Πx1x1 = ( (( 2//3*μ/Re*( 2*ex1x1 -ex2x2 -ex3x3 ) )));
Πx2x2 = ( (( 2//3*μ/Re*( 2*ex2x2 -ex1x1 -ex3x3 ) )));
Πx3x3 = ( (( 2//3*μ/Re*( 2*ex3x3 -ex1x1 -ex2x2 ) )));

Πx1x2 = ( (( μ/Re*ex1x2 )));
Πx2x3 = ( (( μ/Re*ex2x3 )));
Πx1x3 = ( (( μ/Re*ex1x3 )));

T2 = [  (( ( 1/(prod(h)).*( Dx1(h2*h3*Πx1x1) + Dx2(h1*h3*Πx1x2) + Dx3(h1*h2*Πx1x3) ) ) )) + (( ( Πx1x2./(h1*h2)*Dx2(h1) + Πx1x3./(h1*h3)*Dx3(h1) - Πx2x2./(h1*h2)*Dx1(h2) - Πx3x3./(h1*h3)*Dx1(h3) ) )) - 1/h1*Dx1(P) ;
        (( ( 1/(prod(h)).*( Dx1(h2*h3*Πx1x2) + Dx2(h1*h3*Πx2x2) + Dx3(h1*h2*Πx2x3) ) ) )) + (( ( Πx1x2./(h1*h2)*Dx1(h2) + Πx2x3./(h2*h3)*Dx3(h2) - Πx1x1./(h1*h2)*Dx2(h1) - Πx3x3./(h2*h3)*Dx2(h3) ) )) - 1/h1*Dx2(P) ;
        (( ( 1/(prod(h)).*( Dx1(h2*h3*Πx1x3) + Dx2(h1*h3*Πx2x3) + Dx3(h1*h2*Πx3x3) ) ) )) + (( ( Πx1x3./(h1*h3)*Dx1(h3) + Πx2x3./(h2*h3)*Dx2(h3) - Πx1x1./(h1*h3)*Dx3(h1) - Πx2x2./(h2*h3)*Dx3(h2) ) )) - 1/h1*Dx3(P) ];

Cp = 1; # non dimensional

ET1 = (( ρ*Cp*( Dt(T) + U[1]/h[1]*Dx1(T) + U[2]/h[2]*Dx2(T) + U[3]/h[3]*Dx3(T)  )  ));

ET2 = (( (γ-1)*Mach^2*( Dt(P) + U[1]/h[1]*Dx1(P) + U[2]/h[2]*Dx2(P) + U[3]/h[3]*Dx3(P)  )  )); # non dimensional

ET3 =  (( 1/(Pr*Re)*( 1/(prod(h))*( Dx1(κ*h[2]*h[3]/h[1]*Dx1(T)) + Dx2(κ*h[1]*h[3]/h[2]*Dx2(T)) + Dx3(κ*h[1]*h[2]/h[3]*Dx3(T)) ) )  )); # non dimensional

ET4 = (( (γ-1)*Mach^2/Re*μ*( 2*( ex1x1.^2 + ex2x2.^2 + ex3x3.^2 ) + ex1x2.^2 + ex2x3.^2 + ex1x3.^2 -2//3*( ex1x1 + ex2x2 + ex3x3 ).^2 ) )); # non dimensional

Ceqn = (expand_derivatives(( Dt(ρ) + ( 1/(prod(h))*( Differential(x[1])(h[2]*h[3]*(ρ*U[1])) + Differential(x[2])(h[1]*h[3]*(ρ*U[2])) + Differential(x[3])(h[1]*h[2]*(ρ*U[3])) )) ) ));
Meqn = expand_derivatives( (((T1))) - (((T2))) ); #Non dimensional form
Eneqn = ((ET1)) - ((ET2)) - ((ET3)) - ((ET4));

open("Test.txt","w") do file    
    println(file,"Meqn[1] = ",replace(string(Meqn[1]), "(x1, x2, x3)"=>"") );    
end         

