using DifferentialEquations;
using Symbolics;
using Latexify;
using SymbolicUtils;
include("NS_eqn.jl");

@variables t x₁ x₂ x₃ h₁(x₁,x₂,x₃) h₂(x₁,x₂,x₃) h₃(x₁,x₂,x₃) U₁(x₁,x₂,x₃) U₂(x₁,x₂,x₃) U₃(x₁,x₂,x₃) ρ(x₁,x₂,x₃) ;
@variables T(x₁,x₂,x₃) P(x₁,x₂,x₃) Cₚ Rgas κ(T) μ(T) λ Re γ Mach Pᵣ; 

@variables h₁x₁ h₂x₁ h₃x₁ h₁x₂ h₂x₂ h₃x₂ h₁x₃ h₂x₃ h₃x₃;
@variables h₁x₁x₁ h₂x₁x₁ h₃x₁x₁ h₁x₂x₂ h₂x₂x₂ h₃x₂x₂ h₁x₃x₃ h₂x₃x₃ h₃x₃x₃;
@variables h₁x₁x₂ h₂x₁x₂ h₃x₁x₂ h₁x₂x₃ h₂x₂x₃ h₃x₂x₃ h₁x₁x₃ h₂x₁x₃ h₃x₁x₃;

@variables U₁x₁ U₂x₁ U₃x₁ U₁x₂ U₂x₂ U₃x₂ U₁x₃ U₂x₃ U₃x₃;
@variables U₁x₁x₁ U₂x₁x₁ U₃x₁x₁ U₁x₂x₂ U₂x₂x₂ U₃x₂x₂ U₁x₃x₃ U₂x₃x₃ U₃x₃x₃;
@variables U₁x₁x₂ U₂x₁x₂ U₃x₁x₂ U₁x₂x₃ U₂x₂x₃ U₃x₂x₃ U₁x₁x₃ U₂x₁x₃ U₃x₁x₃;

@variables Tx₁ Tx₂ Tx₃;
@variables Tx₁x₁ Tx₂x₂ Tx₃x₃;
@variables Tx₁x₂ Tx₂x₃ Tx₁x₃;

@variables ρx₁ ρx₂ ρx₃;
@variables ρx₁x₁ ρx₂x₂ ρx₃x₃;
@variables ρx₁x₂ ρx₂x₃ ρx₁x₃;

@variables μ_T κ_T μ_TT κ_TT;

# h₁ = 1; h₂ = 1; h₃ = 1;
# h = [h₁ h₂ h₃];
# x = [x₁ x₂ x₃];
# U = [U₁ U₂ U₃];

Dx₁ = Differential(x₁);
Dx₂ = Differential(x₂);
Dx₃ = Differential(x₃);
Dt = Differential(t);
DT = Differential(T);
Dx₁x₁ = Dx₁^2; Dx₂x₂ = Dx₂^2; Dx₃x₃ = Dx₃^2;
Dx₁x₂ = Dx₁*Dx₂; Dx₂x₃ = Dx₂*Dx₃; Dx₁x₃ = Dx₁*Dx₃;
Dx₂x₁ = Dx₂*Dx₁; Dx₃x₂ = Dx₃*Dx₂; Dx₃x₁ = Dx₃*Dx₁;
DTT = DT^2;

llist = [ Dx₁(h₁), Dx₁(h₂), Dx₁(h₃), Dx₂(h₁), Dx₂(h₂), Dx₂(h₃), Dx₃(h₁), Dx₃(h₂), Dx₃(h₃),
          Dx₁x₁(h₁), Dx₁x₁(h₂), Dx₁x₁(h₃), Dx₂x₂(h₁), Dx₂x₂(h₂), Dx₂x₂(h₃), Dx₃x₃(h₁), Dx₃x₃(h₂), Dx₃x₃(h₃),
          Dx₁x₂(h₁), Dx₁x₂(h₂), Dx₁x₂(h₃), Dx₂x₃(h₁), Dx₂x₃(h₂), Dx₂x₃(h₃), Dx₁x₃(h₁), Dx₁x₃(h₂), Dx₁x₃(h₃),
          Dx₂x₁(h₁), Dx₂x₁(h₂), Dx₂x₁(h₃), Dx₃x₂(h₁), Dx₃x₂(h₂), Dx₃x₂(h₃), Dx₃x₁(h₁), Dx₃x₁(h₂), Dx₃x₁(h₃),
          Dx₁(U₁), Dx₁(U₂), Dx₁(U₃), Dx₂(U₁), Dx₂(U₂), Dx₂(U₃), Dx₃(U₁), Dx₃(U₂), Dx₃(U₃),
          Dx₁x₁(U₁), Dx₁x₁(U₂), Dx₁x₁(U₃), Dx₂x₂(U₁), Dx₂x₂(U₂), Dx₂x₂(U₃), Dx₃x₃(U₁), Dx₃x₃(U₂), Dx₃x₃(U₃),
          Dx₁x₂(U₁), Dx₁x₂(U₂), Dx₁x₂(U₃), Dx₂x₃(U₁), Dx₂x₃(U₂), Dx₂x₃(U₃), Dx₁x₃(U₁), Dx₁x₃(U₂), Dx₁x₃(U₃),
          Dx₂x₁(U₁), Dx₂x₁(U₂), Dx₂x₁(U₃), Dx₃x₂(U₁), Dx₃x₂(U₂), Dx₃x₂(U₃), Dx₃x₁(U₁), Dx₃x₁(U₂), Dx₃x₁(U₃),
          Dx₁(T), Dx₂(T), Dx₃(T), Dx₁(ρ), Dx₂(ρ), Dx₃(ρ),
          Dx₁x₁(T), Dx₂x₂(T), Dx₃x₃(T), Dx₁x₁(ρ), Dx₂x₂(ρ), Dx₃x₃(ρ),
          Dx₁x₂(T), Dx₂x₃(T), Dx₁x₃(T), Dx₁x₂(ρ), Dx₂x₃(ρ), Dx₁x₃(ρ),
          Dx₂x₁(T), Dx₃x₂(T), Dx₃x₁(T), Dx₂x₁(ρ), Dx₃x₂(ρ), Dx₃x₁(ρ),
          DT(μ), DT(κ), DTT(μ), DTT(κ) ];

rlist = [ h₁x₁, h₂x₁, h₃x₁, h₁x₂, h₂x₂, h₃x₂, h₁x₃, h₂x₃, h₃x₃,
          h₁x₁x₁, h₂x₁x₁, h₃x₁x₁, h₁x₂x₂, h₂x₂x₂, h₃x₂x₂, h₁x₃x₃, h₂x₃x₃, h₃x₃x₃,
          h₁x₁x₂, h₂x₁x₂, h₃x₁x₂, h₁x₂x₃, h₂x₂x₃, h₃x₂x₃, h₁x₁x₃, h₂x₁x₃, h₃x₁x₃,
          h₁x₁x₂, h₂x₁x₂, h₃x₁x₂, h₁x₂x₃, h₂x₂x₃, h₃x₂x₃, h₁x₁x₃, h₂x₁x₃, h₃x₁x₃,
          U₁x₁, U₂x₁, U₃x₁, U₁x₂, U₂x₂, U₃x₂, U₁x₃, U₂x₃, U₃x₃,
          U₁x₁x₁, U₂x₁x₁, U₃x₁x₁, U₁x₂x₂, U₂x₂x₂, U₃x₂x₂, U₁x₃x₃, U₂x₃x₃, U₃x₃x₃,
          U₁x₁x₂, U₂x₁x₂, U₃x₁x₂, U₁x₂x₃, U₂x₂x₃, U₃x₂x₃, U₁x₁x₃, U₂x₁x₃, U₃x₁x₃,
          U₁x₁x₂, U₂x₁x₂, U₃x₁x₂, U₁x₂x₃, U₂x₂x₃, U₃x₂x₃, U₁x₁x₃, U₂x₁x₃, U₃x₁x₃,
          Tx₁, Tx₂, Tx₃, ρx₁, ρx₂, ρx₃,
          Tx₁x₁, Tx₂x₂, Tx₃x₃, ρx₁x₁, ρx₂x₂, ρx₃x₃,
          Tx₁x₂, Tx₂x₃, Tx₁x₃, ρx₁x₂, ρx₂x₃, ρx₁x₃,
          Tx₁x₂, Tx₂x₃, Tx₁x₃, ρx₁x₂, ρx₂x₃, ρx₁x₃,
          μ_T, κ_T, μ_TT, κ_TT ];


@variables δu(x₁,x₂,x₃,t) δv(x₁,x₂,x₃,t) δw(x₁,x₂,x₃,t) δT(x₁,x₂,x₃,t) δρ(x₁,x₂,x₃,t) δP(x₁,x₂,x₃,t) ϵ i ω;
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

Eqns_base = NS_eqn( x₁, x₂, x₃, h₁, h₂, h₃, vars_base );
dist_pars = [δρ, δu, δv, δw, δT];
vars_dist = [ρ + ϵ*δρ, U₁ + ϵ*δu, U₂ + ϵ*δv, U₃ + ϵ*δw, T + ϵ*δT, P + ϵ*δP, κ + DT(κ)*(ϵ*δT) , μ + DT(μ)*(ϵ*δT)];

Eqns_dist = NS_eqn( x₁, x₂, x₃, h₁, h₂, h₃, vars_dist );

Eqns_base[1] = expand_derivatives(Eqns_base[1]);
Eqns_base[2] = expand_derivatives(Eqns_base[2]);
Eqns_base[3] = expand_derivatives(Eqns_base[3]);
Eqns_base[4] = expand_derivatives(Eqns_base[4]);
Eqns_base[5] = expand_derivatives(Eqns_base[5]);

Eqns_dist[1] = expand_derivatives(Eqns_dist[1]);
Eqns_dist[2] = expand_derivatives(Eqns_dist[2]);
Eqns_dist[3] = expand_derivatives(Eqns_dist[3]);
Eqns_dist[4] = expand_derivatives(Eqns_dist[4]);
Eqns_dist[5] = expand_derivatives(Eqns_dist[5]);

Eqns_base = substitute(( Symbolics.value(Eqns_base)),Dict( llist.=>rlist ));
Eqns_dist = substitute(( Symbolics.value(Eqns_dist)),Dict( llist.=>rlist ));

Eqns_dist = simplify(expand_derivatives(Eqns_dist))-simplify(expand_derivatives(Eqns_base));
Eqns_dist = simplify(expand_derivatives(Eqns_dist))./ϵ;
Eqns_dist = Symbolics.substitute(simplify.(expand_derivatives.(Eqns_dist)),Dict([ϵ=>0]));

writeoutput = 1;
# n,d = (Symbolics.arguments(simplify(Symbolics.value(Eqns_base[5]))))
# terms = Symbolics.arguments(simplify(n));
# open("debug.txt","a") do file
#     for ii = 1:length(terms)
#         println(file,simplify(terms[ii]./d))
#     end 
# end
for ii = 1:5
    local tc1,tc2
    if length(Symbolics.arguments(Symbolics.value(simplify(Eqns_dist[ii])))) > 2 
        tc1 = Symbolics.arguments(Symbolics.value(simplify(Eqns_dist[ii])));
        tc1 = sum(tc1);
        tc2 = 1;
    else
        tc1,tc2 =  Symbolics.arguments(Symbolics.value(simplify(Eqns_dist[ii])));
    end;
        # println( "Writing the Eqn_no  :  ",ii,"  onto the output file " );
    for jj = 1:5
        A[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),(dist_pars[jj]))./tc2);
        B[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁(dist_pars[jj]))./tc2);
        C[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₁(dist_pars[jj]))./tc2);
        D[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂(dist_pars[jj]))./tc2);
        E[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂x₂(dist_pars[jj]))./tc2);
        F[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃(dist_pars[jj]))./tc2);
        G[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃x₃(dist_pars[jj]))./tc2);
        H[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₂(dist_pars[jj]))./tc2 + Symbolics.coeff((expand_derivatives(tc1)),Dx₂x₁(dist_pars[jj]))./tc2);
        I[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂x₃(dist_pars[jj]))./tc2 + Symbolics.coeff((expand_derivatives(tc1)),Dx₃x₂(dist_pars[jj]))./tc2);
        J[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃x₁(dist_pars[jj]))./tc2 + Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₃(dist_pars[jj]))./tc2);
        K[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(dist_pars[jj]))./tc2);

        # if writeoutput == 1
        #     open("Stability_Equations.txt","a") do file
        #         println(file,"A",ii,jj," = " , A[ii,jj],"  ;  ",'\n' );
        #         println(file,"B",ii,jj," = " , B[ii,jj],"  ;  ",'\n' );
        #         println(file,"C",ii,jj," = " , C[ii,jj],"  ;  ",'\n' ); 
        #         println(file,"D",ii,jj," = " , D[ii,jj],"  ;  ",'\n' );
        #         println(file,"E",ii,jj," = " , E[ii,jj],"  ;  ",'\n' );
        #         println(file,"F",ii,jj," = " , F[ii,jj],"  ;  ",'\n' );
        #         println(file,"G",ii,jj," = " , G[ii,jj],"  ;  ",'\n' );
        #         println(file,"H",ii,jj," = " , H[ii,jj],"  ;  ",'\n' );
        #         println(file,"I",ii,jj," = " , I[ii,jj],"  ;  ",'\n' );
        #         println(file,"J",ii,jj," = " , J[ii,jj],"  ;  ",'\n' );
        #         println(file,"K",ii,jj," = " , K[ii,jj],"  ;  ",'\n' );
        #     end
        # end

    end

end

if writeoutput == 1
    open("Stability_Equations_curv.txt","a") do file
        ## A
        for ii = 1:5
            for jj = 1:5                
                println(file,"A",ii,jj," = " , A[ii,jj],"  ;  ",'\n' );                
            end
        end
        ## B
        for ii = 1:5
            for jj = 1:5                
                println(file,"B",ii,jj," = " , B[ii,jj],"  ;  ",'\n' );                
            end
        end
        ## C
        for ii = 1:5
            for jj = 1:5                
                println(file,"C",ii,jj," = " , C[ii,jj],"  ;  ",'\n' );                
            end
        end
        ## D
        for ii = 1:5
            for jj = 1:5                
                println(file,"D",ii,jj," = " , D[ii,jj],"  ;  ",'\n' );                
            end
        end
        ## E
        for ii = 1:5
            for jj = 1:5                
                println(file,"E",ii,jj," = " , E[ii,jj],"  ;  ",'\n' );                
            end
        end
        ## F
        for ii = 1:5
            for jj = 1:5                
                println(file,"F",ii,jj," = " , F[ii,jj],"  ;  ",'\n' );                
            end
        end
        ## G
        for ii = 1:5
            for jj = 1:5                
                println(file,"G",ii,jj," = " , G[ii,jj],"  ;  ",'\n' );                
            end
        end
        ## H
        for ii = 1:5
            for jj = 1:5                
                println(file,"H",ii,jj," = " , H[ii,jj],"  ;  ",'\n' );                
            end
        end
        ## I
        for ii = 1:5
            for jj = 1:5                
                println(file,"I",ii,jj," = " , I[ii,jj],"  ;  ",'\n' );                
            end
        end
        ## J
        for ii = 1:5
            for jj = 1:5                
                println(file,"J",ii,jj," = " , J[ii,jj],"  ;  ",'\n' );                
            end
        end
        ## K
        for ii = 1:5
            for jj = 1:5                
                println(file,"K",ii,jj," = " , K[ii,jj],"  ;  ",'\n' );
            end
        end
    end
end



####### Continuity Equation ########
# Ceqn_dist = simplify(expand_derivatives(Eqns_dist[1]))
# Ceqn_dist = simplify(expand_derivatives(Ceqn_dist))-simplify(expand_derivatives(Ceqn_base));
# Ceqn_dist = simplify(expand_derivatives(Ceqn_dist))./(ϵ);
# Ceqn_dist = Symbolics.substitute(simplify(expand_derivatives(Ceqn_dist)),Dict([ϵ=>0]));

# if length(Symbolics.arguments(Symbolics.value(simplify(Ceqn_dist)))) > 2 
#     tc1 = Symbolics.arguments(Symbolics.value(simplify(Ceqn_dist)));
#     tc1 = sum(tc1);
#     tc2 = 1;
# else
#     tc1,tc2 =  Symbolics.arguments(Symbolics.value(simplify(Ceqn_dist)));
# end;

#= Formulate the equations on the following format 

A + B*Dx₁ + C*Dx₁x₁ + D̂*Dx₂ + E*Dx₂x₂ + F*Dx₃ + G*Dx₃ + H*Dx₁x₂ + I*Dx₂x₃ + J*Dx₁x₃ + K*(-1i*ω*t) = 0;

# =#
# A11 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δρ))./tc2)
# A12 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δu))./tc2)
# A13 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δv))./tc2)
# A14 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δw))./tc2)
# A15 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δT))./tc2)
# B11 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁(δρ))./tc2)
# B12 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁(δu))./tc2)
# B13 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁(δv))./tc2)
# B14 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁(δw))./tc2)
# B15 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁(δT))./tc2)
# D11 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂(δρ))./tc2)
# D12 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂(δu))./tc2)
# D13 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂(δv))./tc2)
# D14 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂(δw))./tc2)
# D15 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂(δT))./tc2)
# F11 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃(δρ))./tc2)
# F12 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃(δu))./tc2)
# F13 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃(δv))./tc2)
# F14 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃(δw))./tc2)
# F15 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃(δT))./tc2)
# K11 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δρ))./tc2)
# K12 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δu))./tc2)
# K13 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δv))./tc2)
# K14 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δw))./tc2)
# K15 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δT))./tc2)

# ####### momentum Equation ########

# mom_eqn_dist = simplify(expand_derivatives(Meqn_dist))-((Meqn_base));
# mom_eqn_dist = simplify(expand_derivatives(mom_eqn_dist))./(ϵ);
# mom_eqn_dist_mod = mom_eqn_dist;
# for ii = 1:3
#     mom_eqn_dist_mod[ii] = Symbolics.substitute(simplify(mom_eqn_dist[ii]),Dict([ϵ=>0]));
# end
# # x₁_mom_eqn_dist = Symbolics.substitute(simplify(mom_eqn_dist[1]),Dict([ϵ=>0]));
# # x₂_mom_eqn_dist = Symbolics.substitute(simplify(mom_eqn_dist[2]),Dict([ϵ=>0]));
# # x₃_mom_eqn_dist = Symbolics.substitute(simplify(mom_eqn_dist[3]),Dict([ϵ=>0]));

# if length(Symbolics.arguments(Symbolics.value(simplify(mom_eqn_dist_mod[1])))) > 2 
#     tc1 = Symbolics.arguments(Symbolics.value(simplify(mom_eqn_dist_mod[1])));
#     tc1 = sum(tc1);
#     tc2 = 1;
# else
#     tc1,tc2 =  Symbolics.arguments(Symbolics.value(simplify(x₁_mom_eqn_dist)));
# end;

# #= Formulate the equations on the following format 

# A + B*Dx₁ + C*Dx₁x₁ + D̂*Dx₂ + E*Dx₂x₂ + F*Dx₃ + G*Dx₃ + H*Dx₁x₂ + I*Dx₂x₃ + J*Dx₁x₃ + K*(-1i*ω*t) = 0;

# =#

# A21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δρ))./tc2);
# A22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δu))./tc2);
# A23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δv))./tc2);
# A24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δw))./tc2);
# A25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δT))./tc2);

# B21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁(δρ))./tc2);
# B22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁(δu))./tc2);
# B23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁(δv))./tc2);
# B24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁(δw))./tc2);
# B25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁(δT))./tc2);

# C21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₁(δρ))./tc2);
# C22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₁(δu))./tc2);
# C23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₁(δv))./tc2);
# C24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₁(δw))./tc2);
# C25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₁(δT))./tc2);

# D21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂(δρ))./tc2);
# D22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂(δu))./tc2);
# D23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂(δv))./tc2);
# D24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂(δw))./tc2);
# D25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂(δT))./tc2);

# E21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂x₂(δρ))./tc2);
# E22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂x₂(δu))./tc2);
# E23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂x₂(δv))./tc2);
# E24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂x₂(δw))./tc2);
# E25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂x₂(δT))./tc2);

# F21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃(δρ))./tc2);
# F22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃(δu))./tc2);
# F23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃(δv))./tc2);
# F24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃(δw))./tc2);
# F25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃(δT))./tc2);

# G21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃x₃(δρ))./tc2);
# G22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃x₃(δu))./tc2);
# G23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃x₃(δv))./tc2);
# G24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃x₃(δw))./tc2);
# G25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₃x₃(δT))./tc2);

# H21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₂(δρ))./tc2);
# H22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₂(δu))./tc2);
# H23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₂(δv))./tc2);
# H24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₂(δw))./tc2);
# H25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₂(δT))./tc2);

# I21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂x₃(δρ))./tc2);
# I22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂x₃(δu))./tc2);
# I23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂x₃(δv))./tc2);
# I24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂x₃(δw))./tc2);
# I25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₂x₃(δT))./tc2);

# J21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₃(δρ))./tc2);
# J22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₃(δu))./tc2);
# J23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₃(δv))./tc2);
# J24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₃(δw))./tc2);
# J25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx₁x₃(δT))./tc2);

# K21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δρ))./tc2);
# K22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δu))./tc2);
# K23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δv))./tc2);
# K24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δw))./tc2);
# K25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δT))./tc2);

# (latexify(simplify(:Ceqn_tr~(( Symbolics.value( ( ( expand_derivatives(T2[1]) ) )))),expand=true),render = true,cdot=true,adjustment=:c,snakecase=true));



