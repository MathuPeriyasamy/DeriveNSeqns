
function NS_eqn_for_BL(  x₁, x₂, x₃, h₁, h₂, h₃, vars )

# @variables t x₁ x₂ x₃ h₁(x₁,x₂,x₃) h₂(x₁,x₂,x₃) h₃(x₁,x₂,x₃) U₁(x₁,x₂,x₃) U₂(x₁,x₂,x₃) U₃(x₁,x₂,x₃) ρ(x₁,x₂,x₃) ;
@variables  Cₚ Rgas λ Re γ Mach Pᵣ; 
    
ρ = vars[1];
U₁ = vars[2];
U₂ = vars[3];
U₃ = vars[4];
T = vars[5];
P = vars[6];
κ = vars[7];
μ = vars[8];

h = [h₁ h₂ h₃];
x = [x₁ x₂ x₃];
U = [U₁ U₂ U₃];

Dx₁ = Differential(x₁);
Dx₂ = Differential(x₂);
Dx₃ = Differential(x₃);
Dt = Differential(t);
DT = Differential(T);

# P = ρ*T/(γ*Mach^2);

T1 = ρ*[    simplify( (( U₁/h₁*Dx₁(U₁) + U₂/h₂*Dx₂(U₁) + U₃/h₃*Dx₃(U₁) + U₁*U₂/(h₁*h₂)*Dx₂(h₁) + U₁*U₃/(h₁*h₃)*Dx₃(h₁) - U₂^2/(h₁*h₂)*Dx₁(h₂) - U₃^2/(h₁*h₃)*Dx₁(h₃)  )) );
            simplify( (( 0 )) );
            simplify( (( U₁/h₁*Dx₁(U₃) + U₂/h₂*Dx₂(U₃) + U₃/h₃*Dx₃(U₃) + U₂*U₃/(h₂*h₃)*Dx₂(h₃) + U₁*U₃/(h₁*h₃)*Dx₁(h₃) - U₂^2/(h₂*h₃)*Dx₃(h₂) - U₁^2/(h₁*h₃)*Dx₃(h₁)  )) )];


ex₁x₁ = simplify( ((  0  )));
ex₂x₂ = simplify( ((  0   )));
ex₃x₃ = simplify( ((  0 )));

ex₁x₂ = simplify( (( h[1]/h[2]*Dx₂(U[1]/h[1]) )));
ex₂x₃ = simplify( (( h[3]/h[2]*Dx₂(U[3]/h[3])  )));
ex₁x₃ = simplify( (( 0 )));



Πx₁x₁ = simplify( (( 2//3*μ*( 2*ex₁x₁ -ex₂x₂ -ex₃x₃ ) ))) ;
Πx₂x₂ = simplify( (( 2//3*μ*( 2*ex₂x₂ -ex₁x₁ -ex₃x₃ ) ))) ;
Πx₃x₃ = simplify( (( 2//3*μ*( 2*ex₃x₃ -ex₁x₁ -ex₂x₂ ) ))) ;

Πx₁x₂ = simplify( (( μ*ex₁x₂ )));
Πx₂x₃ = simplify( (( μ*ex₂x₃ )));
Πx₁x₃ = simplify( (( μ*ex₁x₃ )));

<<<<<<< HEAD
T2 = [  - 1/h₁*Dx₁(P) + (( ( 1/(prod(h)).*(  Dx₂(h₁*h₃*Πx₁x₂)  ) ) )) + (( ( Πx₁x₂./(h₁*h₂)*Dx₂(h₁)  ) ));
        - 1/h₂*Dx₂(P) ;
        - 1/h₃*Dx₃(P) + (( ( 1/(prod(h)).*(  Dx₂(h₁*h₃*Πx₂x₃)  ) ) )) + (( (  Πx₂x₃./(h₂*h₃)*Dx₂(h₃)  ) ))];
=======
T2 = [  - 1/h₁*Dx₁(P) + Symbolics.diff2term(( ( 1/(prod(h)).*(  Dx₂(h₁*h₃*Πx₁x₂)  ) ) )) + Symbolics.diff2term(( ( Πx₁x₂./(h₁*h₂)*Dx₂(h₁)  ) ));
        - 1/h₂*Dx₂(P) ;
        - 1/h₃*Dx₃(P) + Symbolics.diff2term(( ( 1/(prod(h)).*(  Dx₂(h₁*h₃*Πx₂x₃)  ) ) )) + Symbolics.diff2term(( (  Πx₂x₃./(h₂*h₃)*Dx₂(h₃)  ) ))];
>>>>>>> 6af94779961f0ef42c01941e8163efe01e139261


# κ = Cₚ*μ/Pᵣ;
# κ = μ; # non dimensional
Cₚ = 1; # non dimensional

ET1 = simplify(( ρ*Cₚ*(  U[1]/h[1]*Dx₁(T) + U[2]/h[2]*Dx₂(T) + U[3]/h[3]*Dx₃(T)  )  ));

# ET2 = ( Dt(P) + U[1]/h[1]*Dx₁(P) + U[2]/h[2]*Dx₂(P) + U[3]/h[3]*Dx₃(P)  );
ET2 = simplify(( (  U[1]/h[1]*Dx₁(P) + U[2]/h[2]*Dx₂(P) + U[3]/h[3]*Dx₃(P)  )  )); # non dimensional

# ET3 = 1/(prod(h))*( Dx₁(κ*h[2]*h[3]/h[1]*Dx₁(T)) + Dx₂(κ*h[1]*h[3]/h[2]*Dx₂(T)) + Dx₃(κ*h[1]*h[2]/h[3]*Dx₃(T)) );
ET3 =  simplify(( 1/(Pᵣ)*( 1/(prod(h))*(   Dx₂(κ*h[1]*h[3]/h[2]*Dx₂(T))  ) )  )); # non dimensional

# ET4 = μ*( 2*( ex₁x₁.^2 + ex₂x₂.^2 + ex₃x₃.^2 ) + ex₁x₂.^2 + ex₂x₃.^2 + ex₁x₃.^2 -2/3*( ex₁x₁ + ex₂x₂ + ex₃x₃ ).^2 );

ET4 = simplify(( μ*( 2*( ex₁x₁.^2 + ex₂x₂.^2 + ex₃x₃.^2 ) + ex₁x₂.^2 + ex₂x₃.^2 + ex₁x₃.^2 -2//3*( ex₁x₁ + ex₂x₂ + ex₃x₃ ).^2 ) )); # non dimensional

Ceqn = (Symbolics.value((  + ( 1/(prod(h))*( Differential(x[1])(h[2]*h[3]*(ρ*U[1])) + Differential(x[2])(h[1]*h[3]*(ρ*U[2])) + Differential(x[3])(h[1]*h[2]*(ρ*U[3])) )) ) ));
Meqn = ( (((T1))) - (((T2))) ); #Non dimensional form
Eneqn = Symbolics.value((ET1)) - Symbolics.value((ET2)) - Symbolics.value((ET3)) - Symbolics.value((ET4));

# Ceqn = Symbolics.substitute((Symbolics.value((Ceqn))),Dict([x₁=>x₁ᵢ,x₂=>x₂ᵢ,x₃=>x₃ᵢ,h₁=>h₁ᵢ,h₂=>h₂ᵢ,h₃=>h₃ᵢ,U₁=>U₁ᵢ,U₂=>U₂ᵢ,U₃=>U₃ᵢ,T=>Tᵢ,P=>Pᵢ,κ=>κᵢ,μ=>μᵢ]));
# Meqn = Symbolics.substitute((Symbolics.value((Meqn))),Dict([x₁=>x₁ᵢ,x₂=>x₂ᵢ,x₃=>x₃ᵢ,h₁=>h₁ᵢ,h₂=>h₂ᵢ,h₃=>h₃ᵢ,U₁=>U₁ᵢ,U₂=>U₂ᵢ,U₃=>U₃ᵢ,T=>Tᵢ,P=>Pᵢ,κ=>κᵢ,μ=>μᵢ]));
# Eneqn = Symbolics.substitute((Symbolics.value((Eneqn))),Dict([x₁=>x₁ᵢ,x₂=>x₂ᵢ,x₃=>x₃ᵢ,h₁=>h₁ᵢ,h₂=>h₂ᵢ,h₃=>h₃ᵢ,U₁=>U₁ᵢ,U₂=>U₂ᵢ,U₃=>U₃ᵢ,T=>Tᵢ,P=>Pᵢ,κ=>κᵢ,μ=>μᵢ]));

return [Ceqn;Meqn[1];Meqn[2];Meqn[3];Eneqn]


end