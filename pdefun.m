%% Code Equations
% Formulate a system of 12 PDE equations. 8 partial differential equations (CO2, DOC, O2, NH4, NO2, NO3, N2O, N2) and 3 odinary differential equations (b1,b2,b3).
% PDE: 2 independent varibales x and t, dependent variable u. ODE: t and u.
% x [m], t [day], u [mmol/L]
% indexes of vector u. 1: CO2, 2: DOC (mmolC/L), 3: O2, 4: b_AER, 5: NH4, 6:
% NO2, 7: NO3, 8: N2O, 9: b_AOB, 10: N2, 11: b_DEN, 12: b_NOB.
% INPUTS: 
% x,t,u,dudx: internal handles for pdepe
% msInfo: measurement information
% par: parameters 
% Use interp1 to interpolate data_thetaR(x) and data_Deff(x) to the grid point x for which the solver requires information.
%%
function [c,f,s] = pdefun(x,t,u,dudx,msInfo,par) % Equation to solve
u = real(u);
u(u<0.000)=0; % avoid negative concentrations
if u(10)<33.000
    u(10)=33;
end

kC_CO2_r = par(1); % mmol C/L
kO2_CO2_r = par(2); % mmol O2/L in air
kNH4_NO2_n = par(3); % mmol N/L
kO2_NO2_n = par(4); % mmol N/L
kNO2_NO3_n = par(5); % mmol N/L
kO2_NO3_n = par(6); % mmol N/L
kNH4_N2O_n = par(7); % mmol N/L
kO2_N2O_n = par(8); % mmol N/L
kNO2_N2O_nd = par(9); % mmol N/L

kNH4_N2O_nd = par(10); % mmol N/L
kO2_N2O_nd = par(11); % mmol N/L
kNO3_NO2_dn = par(12); % mmol N/L
kC_NO2_dn = par(13); % mmol C/L
kNO2_N2O_dn = par(14); % mmol N/L
kC_N2O_dn = par(15); % mmol C/L
kN2O_N2_dn = par(16); % mmol N2O/L
kC_N2_dn = par(17); % mmol C/L
kI_NO2_dn = par(18); % mmol O2/L
kI_N2O_dn = par(19); % mmol O2/L
kI_N2O_nd = par(20); % mmol O2/L
kI_N2_dn = par(21); % mmol O2/L
u_CO2_r = par(22);  % mmol/g biomass/day
u_NO2_n = par(23); % mmol/g biomass/day
u_NO3_n = par(24); % mmol/g biomass/day
u_N2O_n = par(25); % mmol/g biomass/day
u_N2O_nd = par(26); % mmol/g biomass/day
u_NO2_dn = par(27); % mmol/g biomass/day
u_N2O_dn = par(28); % mmol N2O/g biomass/day
u_N2_dn = par(29); % mmol N2O/g biomass/day
y_AER = par(30); % g C/g C, yield coefficient
y_AOB = par(31); % g N/g N
y_NOB = par(32); % g N/g N
y_DEN = par(33); % g C/g C
a_AER = par(34); % 1/day, decay rate
a_AOB = par(35); % 1/day
a_NOB = par(36); % 1/day
a_DEN = par(37); % 1/day
SOC_alpha = par(38); %1/day, Mass transfer coefficient of SOC desorption
POC_alpha = par(39); %1/day, Mass transfer coefficient of POC desorption

f_Cbio = 0.53; % C content in microbial biomass
f_Nbio = 0.066; % N content in microbial biomass

icMesh = msInfo.icMesh;
thetaMat = msInfo.thetaMat;
N_NH4 = msInfo.N_NH4;
KF_NH4 = msInfo.KF_NH4;
rhob = msInfo.rhob; % bulk density [kg/L soil]
rhob_gL = 1000*rhob; % bulk density [g/L soil]
DeffMat = msInfo.DeffMat;
theta_tot = msInfo.theta_tot;
restSOC_t0 = msInfo.restSOC_t0;
POC_t0 = msInfo.POC_t0;
biomass_base = msInfo.biomass_base; %[g bio/gC]
b_AER_0 = biomass_base(1);
b_AOB_0 = biomass_base(2);
b_NOB_0 = biomass_base(3);
b_DEN_0 = biomass_base(4);



theta = interp1(icMesh,thetaMat',x,'pchip')';% a function of x, returning interpolated values of theta at x for 11 components
FreuC_NH4 = KF_NH4 * N_NH4 * (18*u(5)).^(N_NH4-1); %KF[mg/kg][L/mg]^N, NH4 molar conc: u(5)[mmol/L], NH4 molar mass: 18[mg/mmol], NH4 mass conc:(18*u(5))[mg/L]
R_NH4 = 1 + rhob./theta(5).*FreuC_NH4; %rhob = 1.4[kg/L soil], FreuC_NH4 [L/kg]
c = theta.*[1 1 1 1 R_NH4 1 1 1 1 1 1 1]';
Deff = interp1(icMesh,DeffMat',x,'pchip')';% [m2/d], a function of x, returning interpolated values of Deff at x for 11 components
f = Deff.*dudx;

ix_DOC = theta_tot.^(-3) .* theta(2).^3; %DOC--2, water-filled porosity
ix_O2 = theta_tot.^(-4/3) .* theta(3).^(4/3); %O2--3, AIR-filled porosity
ix_NH4 = theta_tot.^(-3) .* theta(5).^3; %NH4--5, water-filled porosity
ix_NO2 = theta_tot.^(-3) .* theta(6).^3; %NO2--6, water-filled porosity
ix_NO3 = theta_tot.^(-3) .* theta(7).^3; %NO3--7, water-filled porosity
ix_N2O = theta_tot.^(-4/3) .* theta(8).^(4/3); %N2O--8, AIR-filled porosity
%%%%%% ix_N2O = theta_tot.^(-3) .* theta(2).^3; %N2O--8, water-filled porosity

% !!! ALERT!!!NEED TO Multiply rhob on RHS (Done. 1400 g/L), Need to add Vmax on RHS below (Done)
% 4: b_AER, 9: b_AOB, 11: b_DEN, 12: b_NOB. % source and sink terms S, [mmol/L soil/day]
S_CO2pr_r = rhob_gL.*(u(4)+b_AER_0).*u_CO2_r.*ix_DOC.*u(2)./(kC_CO2_r+ix_DOC.*u(2)).*ix_O2.*u(3)./(kO2_CO2_r+ix_O2.*u(3)); % CO2 production rate, [mmol CO2 produced/L soil/day]
S_NO2pr_n = rhob_gL.*(u(9)+b_AOB_0).*u_NO2_n.*ix_NH4.*u(5)./(kNH4_NO2_n+ix_NH4.*u(5)).*ix_O2.*u(3)./(kO2_NO2_n+ix_O2.*u(3)); % NO2 production rate from N, [mmol NO2 produced/L soil/day]
S_NO3pr_n = rhob_gL.*(u(12)+b_NOB_0).*u_NO3_n.*ix_NO2.*u(6)./(kNO2_NO3_n+ix_NO2.*u(6)).*ix_O2.*u(3)./(kO2_NO3_n+ix_O2.*u(3)); % NO3 production rate from N, [mmol NO3 produced/L soil/day]
S_N2Opr_nn = rhob_gL.*(u(9)+b_AOB_0).*u_N2O_n.*ix_NH4.*u(5)./(kNH4_N2O_n+ix_NH4.*u(5)).*ix_O2.*u(3)./(kO2_N2O_n+ix_O2.*u(3)); % N2O production rate from NN, [mmol N2O produced/L soil/day]
S_N2Opr_nd = rhob_gL.*(u(9)+b_AOB_0).*u_N2O_nd.*ix_NO2.*u(6)./(kNO2_N2O_nd+ix_NO2.*u(6)).*ix_NH4.*u(5)./(kNH4_N2O_nd+ix_NH4.*u(5)).*ix_O2.*u(3)./(kO2_N2O_nd+ix_O2.*u(3)).*kI_N2O_nd./(kI_N2O_nd+ix_O2.*u(3)); % N2O production rate from ND, [mmol N2O produced/L soil/day]
S_NO2pr_dn  = rhob_gL.*(u(11)+b_DEN_0).*u_NO2_dn.*ix_NO3.*u(7)./(kNO3_NO2_dn+ix_NO3.*u(7)).*ix_DOC.*u(2)./(kC_NO2_dn+ix_DOC.*u(2)).*kI_NO2_dn./(kI_NO2_dn+ix_O2.*u(3)); % NO2 production rate from DN, [mmol NO2 produced/L soil/day]
S_N2Opr_dn  = rhob_gL.*(u(11)+b_DEN_0).*u_N2O_dn.*ix_NO2.*u(6)./(kNO2_N2O_dn+ix_NO2.*u(6)).*ix_DOC.*u(2)./(kC_N2O_dn+ix_DOC.*u(2)).*kI_N2O_dn./(kI_N2O_dn+ix_O2.*u(3)); % N2O production rate from DN, [mmol N2O produced/L soil/day]
S_N2pr_dn  = rhob_gL.*(u(11)+b_DEN_0).*u_N2_dn.*ix_N2O.*u(8)./(kN2O_N2_dn+ix_N2O.*u(8)).*ix_DOC.*u(2)./(kC_N2_dn+ix_DOC.*u(2)).*kI_N2_dn./(kI_N2_dn+ix_O2.*u(3)); % N2 production rate from DN, [mmol N2 produced/L soil/day]

S_CO2pr_r(S_CO2pr_r<0.000)=0;
S_NO2pr_n(S_NO2pr_n<0.000)=0;
S_NO3pr_n(S_NO3pr_n<0.000)=0;
S_N2Opr_nn(S_N2Opr_nn<0.000)=0;
S_N2Opr_nd(S_N2Opr_nd<0.000)=0;
S_NO2pr_dn(S_NO2pr_dn<0.000)=0;
S_N2Opr_dn(S_N2Opr_dn<0.000)=0;
S_N2pr_dn(S_N2pr_dn<0.000)=0;

S_CO2pr_dn = 0.5*S_NO2pr_dn + S_N2Opr_dn + 0.5*S_N2pr_dn; %[mmol CO2 produced/L soil/day]

% to calculate SOC decomposed to DOC along the soil profile
S_restSOC_cons = -SOC_alpha.*restSOC_t0.*exp(-SOC_alpha.*t); % gC/g dw/day 
S_DOC_pr_soc = -rhob_gL.*S_restSOC_cons.*1000./12; % mmolC/L soil/day

% to calculate POC decomposed to DOC in the 2 mm center of the soil core
if (x>=0.049) && (x<=0.051)
    S_POC_cons = -POC_alpha.*POC_t0.*exp(-POC_alpha.*t);% gC/g dw/day 
    S_DOC_pr_poc = -rhob_gL.*S_POC_cons.*1000./12; % mmolC/L soil/day
else
    S_DOC_pr_poc = 0;
end

p_a_AER = 1;
p_a_AOB = 1;
p_a_DEN = 1;
p_a_NOB = 1;

s = [S_CO2pr_r + S_CO2pr_dn;... %CO2--1 [mmol CO2 produced/L soil/day]
    -S_CO2pr_r/(1-y_AER) - S_CO2pr_dn/(1-y_DEN) + S_DOC_pr_soc + S_DOC_pr_poc;... %DOC--2
    -S_CO2pr_r - 1.5*S_NO2pr_n - 0.5*S_NO3pr_n - 2.75*S_N2Opr_nn - 0.5*S_N2Opr_nd;... %O2--3
    y_AER/(1-y_AER)/f_Cbio *0.001*12*S_CO2pr_r/rhob_gL - a_AER*u(4)*p_a_AER;... %b_AER, heterotroph--4, y_r[g bio/gC] * 0.001[mol/mmol] * 12 [gC/mol] * S_CO2pr[mmol/L soil/day] / rho_b[g dw/L soil] --> [g bio/g dw/day]
    -(S_NO2pr_n + 2.5*S_N2Opr_nn + S_N2Opr_nd)/(1-y_AOB);... %NH4--5, substrate of Bnit
    S_NO2pr_n + S_NO2pr_dn + 0.5*S_N2Opr_nn - S_NO3pr_n/(1-y_NOB) - S_N2Opr_nd - 2*S_N2Opr_dn;... %NO2--6
    S_NO3pr_n - S_NO2pr_dn;... %NO3--7
    S_N2Opr_nn + S_N2Opr_nd + S_N2Opr_dn - S_N2pr_dn;... %N2O--8
    y_AOB/(1-y_AOB)/f_Nbio *0.001*14*(S_NO2pr_n + 2.5*S_N2Opr_nn + S_N2Opr_nd)/rhob_gL - a_AOB.*u(9)*p_a_AOB;...%b_AOB, nit--9, y_n[g bio/gN] * 0.001[mol/mmol] * 14 [gN/mol] * S_NH4cons[mmol/L soil/day] / rho_b[g dw/L soil] --> [g bio/g dw/day]
    S_N2pr_dn;...%N2--10
    y_DEN/(1-y_DEN)/f_Cbio *0.001*12*S_CO2pr_dn/rhob_gL - a_DEN*u(11)*p_a_DEN;...%b_DEN, denit--11, y_dn[g bio/gN] * 0.001[mol/mmol] * 14 [gN/mol] * S_NO3cons[mmol/L soil/day] / rho_b[g dw/L soil] --> [g bio/g dw/day]
    y_NOB/(1-y_NOB)/f_Nbio *0.001*14*S_NO3pr_n/rhob_gL - a_NOB*u(12)*p_a_NOB...%b_NOB, nit--12, y_n[g bio/gN] * 0.001[mol/mmol] * 14 [gN/mol] * S_NH4cons[mmol/L soil/day] / rho_b[g dw/L soil] --> [g bio/g dw/day]
    ];
end
