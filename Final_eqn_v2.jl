using DifferentialEquations;
using Symbolics;
using Latexify;
using SymbolicUtils;

include("NS_eqn_v2.jl");

writeoutput = 1;
pressureform = 1;
densityform = 1;
#### Verification ==> uncomment depending upon the Coordinate system
## Cartesian system 
# @variables x y z;
# x1 = x; x2 = y; x3 = z;

##  Cylindrical system
@variables r θ z;
x1 = r; x2 = θ; x3 = z;

## Sperical system
# @variables r θ ϕ
# x1 = r; x2 = θ; x3 = ϕ;

@variables t x1 x2 x3 h1(x1,x2,x3) h2(x1,x2,x3) h3(x1,x2,x3) U1(x1,x2,x3) U2(x1,x2,x3) U3(x1,x2,x3) ρ(x1,x2,x3) ;
@variables T(x1,x2,x3) P(x1,x2,x3) Cp Rgas κ(T) μ(T) λ Re γ Mach Pr Ec; 

@variables h1x1 h2x1 h3x1 h1x2 h2x2 h3x2 h1x3 h2x3 h3x3;
@variables h1x1x1 h2x1x1 h3x1x1 h1x2x2 h2x2x2 h3x2x2 h1x3x3 h2x3x3 h3x3x3;
@variables h1x1x2 h2x1x2 h3x1x2 h1x2x3 h2x2x3 h3x2x3 h1x1x3 h2x1x3 h3x1x3;

@variables U1x1 U2x1 U3x1 U1x2 U2x2 U3x2 U1x3 U2x3 U3x3;
@variables U1x1x1 U2x1x1 U3x1x1 U1x2x2 U2x2x2 U3x2x2 U1x3x3 U2x3x3 U3x3x3;
@variables U1x1x2 U2x1x2 U3x1x2 U1x2x3 U2x2x3 U3x2x3 U1x1x3 U2x1x3 U3x1x3;

@variables Tx1 Tx2 Tx3;
@variables Tx1x1 Tx2x2 Tx3x3;
@variables Tx1x2 Tx2x3 Tx1x3;

@variables ρx1 ρx2 ρx3 Px1 Px2 Px3;
@variables ρx1x1 ρx2x2 ρx3x3 Px1x1 Px2x2 Px3x3; 
@variables ρx1x2 ρx2x3 ρx1x3 Px1x2 Px2x3 Px1x3;

@variables μ_T κ_T μ_TT κ_TT;

@variables δu(x1,x2,x3,t) δv(x1,x2,x3,t) δw(x1,x2,x3,t) δT(x1,x2,x3,t) δρ(x1,x2,x3,t) δP(x1,x2,x3,t) ϵ i ω;
###### uncomment corresponding CS
## Cartesian system
# h1 = 1; h2 = 1; h3 = 1;

##  Cylindrical system
h1 = 1; h2 = r; h3 = 1;

## Sperical system

# h1 = 1; h2 = r; h3 = r*sin(θ);

hi = [h1 h2 h3];
xi = [x1 x2 x3];

Dx1 = Differential(x1);
Dx2 = Differential(x2);
Dx3 = Differential(x3);
Dt = Differential(t);
DT = Differential(T);
Dx1x1 = Dx1^2; Dx2x2 = Dx2^2; Dx3x3 = Dx3^2;
Dx1x2 = Dx1*Dx2; Dx2x3 = Dx2*Dx3; Dx1x3 = Dx1*Dx3;
Dx2x1 = Dx2*Dx1; Dx3x2 = Dx3*Dx2; Dx3x1 = Dx3*Dx1;
DTT = DT^2;

if densityform == 1;

    P = ρ*T/(γ*Mach^2);
    δP = (δρ*T+ρ*δT)./(γ*Mach^2);
    dist_pars = [δρ, δu, δv, δw, δT];

    llist = [ Dx1(h1), Dx1(h2), Dx1(h3), Dx2(h1), Dx2(h2), Dx2(h3), Dx3(h1), Dx3(h2), Dx3(h3),
            Dx1x1(h1), Dx1x1(h2), Dx1x1(h3), Dx2x2(h1), Dx2x2(h2), Dx2x2(h3), Dx3x3(h1), Dx3x3(h2), Dx3x3(h3),
            Dx1x2(h1), Dx1x2(h2), Dx1x2(h3), Dx2x3(h1), Dx2x3(h2), Dx2x3(h3), Dx1x3(h1), Dx1x3(h2), Dx1x3(h3),
            Dx2x1(h1), Dx2x1(h2), Dx2x1(h3), Dx3x2(h1), Dx3x2(h2), Dx3x2(h3), Dx3x1(h1), Dx3x1(h2), Dx3x1(h3),
            Dx1(U1), Dx1(U2), Dx1(U3), Dx2(U1), Dx2(U2), Dx2(U3), Dx3(U1), Dx3(U2), Dx3(U3),
            Dx1x1(U1), Dx1x1(U2), Dx1x1(U3), Dx2x2(U1), Dx2x2(U2), Dx2x2(U3), Dx3x3(U1), Dx3x3(U2), Dx3x3(U3),
            Dx1x2(U1), Dx1x2(U2), Dx1x2(U3), Dx2x3(U1), Dx2x3(U2), Dx2x3(U3), Dx1x3(U1), Dx1x3(U2), Dx1x3(U3),
            Dx2x1(U1), Dx2x1(U2), Dx2x1(U3), Dx3x2(U1), Dx3x2(U2), Dx3x2(U3), Dx3x1(U1), Dx3x1(U2), Dx3x1(U3),
            Dx1(T), Dx2(T), Dx3(T), 
            Dx1(ρ), Dx2(ρ), Dx3(ρ), Dx1(P), Dx2(P), Dx3(P),
            Dx1x1(T), Dx2x2(T), Dx3x3(T), 
            Dx1x1(ρ), Dx2x2(ρ), Dx3x3(ρ), Dx1x1(P), Dx2x2(P), Dx3x3(P),
            Dx1x2(T), Dx2x3(T), Dx1x3(T), 
            Dx1x2(ρ), Dx2x3(ρ), Dx1x3(ρ), Dx1x2(P), Dx2x3(P), Dx1x3(P),
            Dx2x1(T), Dx3x2(T), Dx3x1(T), 
            Dx2x1(ρ), Dx3x2(ρ), Dx3x1(ρ), Dx2x1(P), Dx3x2(P), Dx3x1(P),
            DT(μ), DT(κ), DTT(μ), DTT(κ) ];

    rlist = [ h1x1, h2x1, h3x1, h1x2, h2x2, h3x2, h1x3, h2x3, h3x3,
            h1x1x1, h2x1x1, h3x1x1, h1x2x2, h2x2x2, h3x2x2, h1x3x3, h2x3x3, h3x3x3,
            h1x1x2, h2x1x2, h3x1x2, h1x2x3, h2x2x3, h3x2x3, h1x1x3, h2x1x3, h3x1x3,
            h1x1x2, h2x1x2, h3x1x2, h1x2x3, h2x2x3, h3x2x3, h1x1x3, h2x1x3, h3x1x3,
            U1x1, U2x1, U3x1, U1x2, U2x2, U3x2, U1x3, U2x3, U3x3,
            U1x1x1, U2x1x1, U3x1x1, U1x2x2, U2x2x2, U3x2x2, U1x3x3, U2x3x3, U3x3x3,
            U1x1x2, U2x1x2, U3x1x2, U1x2x3, U2x2x3, U3x2x3, U1x1x3, U2x1x3, U3x1x3,
            U1x1x2, U2x1x2, U3x1x2, U1x2x3, U2x2x3, U3x2x3, U1x1x3, U2x1x3, U3x1x3,
            Tx1, Tx2, Tx3, 
            ρx1, ρx2, ρx3, Px1, Px2, Px3,
            Tx1x1, Tx2x2, Tx3x3, 
            ρx1x1, ρx2x2, ρx3x3, Px1x1, Px2x2, Px3x3,
            Tx1x2, Tx2x3, Tx1x3, 
            ρx1x2, ρx2x3, ρx1x3, Px1x2, Px2x3, Px1x3,
            Tx1x2, Tx2x3, Tx1x3, 
            ρx1x2, ρx2x3, ρx1x3, Px1x2, Px2x3, Px1x3,
            μ_T, κ_T, μ_TT, κ_TT ];

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

    vars_base = [ρ, U1, U2, U3, T, P, κ, μ];

    Eqns_base = NS_eqn_v2( xi, hi, vars_base );

    vars_dist = [ρ + ϵ*δρ, U1 + ϵ*δu, U2 + ϵ*δv, U3 + ϵ*δw, T + ϵ*δT, P + ϵ*δP, κ + DT(κ)*(ϵ*δT) , μ + DT(μ)*(ϵ*δT)];

    Eqns_dist = NS_eqn_v2( xi, hi, vars_dist );

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
            B[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1(dist_pars[jj]))./tc2);
            C[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x1(dist_pars[jj]))./tc2);
            D[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2(dist_pars[jj]))./tc2);
            E[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2x2(dist_pars[jj]))./tc2);
            F[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3(dist_pars[jj]))./tc2);
            G[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3x3(dist_pars[jj]))./tc2);
            H[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x2(dist_pars[jj]))./tc2 + Symbolics.coeff((expand_derivatives(tc1)),Dx2x1(dist_pars[jj]))./tc2);
            I[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2x3(dist_pars[jj]))./tc2 + Symbolics.coeff((expand_derivatives(tc1)),Dx3x2(dist_pars[jj]))./tc2);
            J[ii,jj] = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3x1(dist_pars[jj]))./tc2 + Symbolics.coeff((expand_derivatives(tc1)),Dx1x3(dist_pars[jj]))./tc2);
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
        open("Stability_Equations_dens_form.txt","a") do file
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

A + B*Dx1 + C*Dx1x1 + D̂*Dx2 + E*Dx2x2 + F*Dx3 + G*Dx3 + H*Dx1x2 + I*Dx2x3 + J*Dx1x3 + K*(-1i*ω*t) = 0;

# =#
# A11 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δρ))./tc2)
# A12 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δu))./tc2)
# A13 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δv))./tc2)
# A14 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δw))./tc2)
# A15 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δT))./tc2)
# B11 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1(δρ))./tc2)
# B12 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1(δu))./tc2)
# B13 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1(δv))./tc2)
# B14 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1(δw))./tc2)
# B15 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1(δT))./tc2)
# D11 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2(δρ))./tc2)
# D12 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2(δu))./tc2)
# D13 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2(δv))./tc2)
# D14 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2(δw))./tc2)
# D15 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2(δT))./tc2)
# F11 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3(δρ))./tc2)
# F12 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3(δu))./tc2)
# F13 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3(δv))./tc2)
# F14 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3(δw))./tc2)
# F15 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3(δT))./tc2)
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
# # x1_mom_eqn_dist = Symbolics.substitute(simplify(mom_eqn_dist[1]),Dict([ϵ=>0]));
# # x2_mom_eqn_dist = Symbolics.substitute(simplify(mom_eqn_dist[2]),Dict([ϵ=>0]));
# # x3_mom_eqn_dist = Symbolics.substitute(simplify(mom_eqn_dist[3]),Dict([ϵ=>0]));

# if length(Symbolics.arguments(Symbolics.value(simplify(mom_eqn_dist_mod[1])))) > 2 
#     tc1 = Symbolics.arguments(Symbolics.value(simplify(mom_eqn_dist_mod[1])));
#     tc1 = sum(tc1);
#     tc2 = 1;
# else
#     tc1,tc2 =  Symbolics.arguments(Symbolics.value(simplify(x1_mom_eqn_dist)));
# end;

# #= Formulate the equations on the following format 

# A + B*Dx1 + C*Dx1x1 + D̂*Dx2 + E*Dx2x2 + F*Dx3 + G*Dx3 + H*Dx1x2 + I*Dx2x3 + J*Dx1x3 + K*(-1i*ω*t) = 0;

# =#

# A21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δρ))./tc2);
# A22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δu))./tc2);
# A23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δv))./tc2);
# A24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δw))./tc2);
# A25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),(δT))./tc2);

# B21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1(δρ))./tc2);
# B22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1(δu))./tc2);
# B23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1(δv))./tc2);
# B24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1(δw))./tc2);
# B25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1(δT))./tc2);

# C21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x1(δρ))./tc2);
# C22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x1(δu))./tc2);
# C23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x1(δv))./tc2);
# C24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x1(δw))./tc2);
# C25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x1(δT))./tc2);

# D21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2(δρ))./tc2);
# D22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2(δu))./tc2);
# D23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2(δv))./tc2);
# D24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2(δw))./tc2);
# D25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2(δT))./tc2);

# E21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2x2(δρ))./tc2);
# E22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2x2(δu))./tc2);
# E23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2x2(δv))./tc2);
# E24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2x2(δw))./tc2);
# E25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2x2(δT))./tc2);

# F21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3(δρ))./tc2);
# F22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3(δu))./tc2);
# F23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3(δv))./tc2);
# F24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3(δw))./tc2);
# F25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3(δT))./tc2);

# G21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3x3(δρ))./tc2);
# G22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3x3(δu))./tc2);
# G23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3x3(δv))./tc2);
# G24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3x3(δw))./tc2);
# G25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx3x3(δT))./tc2);

# H21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x2(δρ))./tc2);
# H22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x2(δu))./tc2);
# H23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x2(δv))./tc2);
# H24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x2(δw))./tc2);
# H25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x2(δT))./tc2);

# I21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2x3(δρ))./tc2);
# I22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2x3(δu))./tc2);
# I23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2x3(δv))./tc2);
# I24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2x3(δw))./tc2);
# I25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx2x3(δT))./tc2);

# J21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x3(δρ))./tc2);
# J22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x3(δu))./tc2);
# J23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x3(δv))./tc2);
# J24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x3(δw))./tc2);
# J25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dx1x3(δT))./tc2);

# K21 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δρ))./tc2);
# K22 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δu))./tc2);
# K23 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δv))./tc2);
# K24 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δw))./tc2);
# K25 = simplify(Symbolics.coeff((expand_derivatives(tc1)),Dt(δT))./tc2);

# (latexify(simplify(:Ceqn_tr~(( Symbolics.value( ( ( expand_derivatives(T2[1]) ) )))),expand=true),render = true,cdot=true,adjustment=:c,snakecase=true));



