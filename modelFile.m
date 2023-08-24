%%
% SCRIPT FILE FOR RUNNING THE REACTIVE TRANSPORT MODEL (a.k.a. the respiration-nitrification-denitrification model)
% AUTHOR: JIE ZHANG jiezh@agro.au.dk
close all;clear;clc
%% Define the model domain and load experiment information
x = linspace(0,0.1,101); % x[m]
finex = 0.04: 2e-3/20: 0.06;  %%%%% 
x = unique(round(sort([x,finex]),5));
t = 0:6/24:28;  % t[day] output time interval
expInfo = load('Info_30hPa.mat'); %measInfo

parTab = readtable('par.xlsx', 'Sheet', 'Sheet1');% Load parameters
par = parTab.Fitted';
%% Model comuputation: solve the pde system
tic;
opts = odeset('RelTol',1e-3,'AbsTol',1e-6); % default tolerance
m = 0;
sol = pdepe(m,@(x,t,u,dudx)pdefun(x,t,u,dudx,expInfo,par),@(x)pdeic(x,expInfo),@pdebc,x,t,opts);
sol = real(sol);
toc;
%% Postprocessing
%------ Extract concentrations for 12 components-----------------------
u1 = sol(:,:,1); %CO2--1, time steps * depth
u2 = sol(:,:,2);
u3 = sol(:,:,3); %O2--3
u4 = sol(:,:,4);
u5 = sol(:,:,5);
u6 = sol(:,:,6);
u7 = sol(:,:,7);
u8 = sol(:,:,8); %N2O--8
u9 = sol(:,:,9);
u10 = sol(:,:,10); %N2--10
u11 = sol(:,:,11);
u12 = sol(:,:,12);
%------ Calculate gas fluxes--------------------------------------------
x_up = x(2); %0.01; %[m]
x_down = x(end-1);%0.09; %[m]
dx = x_up - 0; % [m]
x_index = find(ismember(x, [x_up,x_down]));
D = interp1(expInfo.icMesh,expInfo.DeffMat',[x_up, x_down],'pchip')'; %gas diff. at the depth of 1 cm to the border, heter. measurements 
% D = interp1(scenexpInfo.icMesh,scenexpInfo.DeffMat',[x_up, x_down],'pchip')'; %gas diff. at the depth of 1 cm to the border, homo. scenrios 
D_CO2 = D(1,:); % [m2/d] 1. CO2
D_N2O = D(8,:); % [m2/d] 8. N2O
D_N2 = D(10,:); % [m2/d] 10. N2
J_CO2_mol = (D_CO2(1)*(u1(:,x_index(1)) - u1(:,1))+D_CO2(2)*(u1(:,x_index(2)) - u1(:,end)))/dx; % [mol CO2/m2/day]
J_N2O_mol = (D_N2O(1)*(u8(:,x_index(1)) - u8(:,1))+D_N2O(2)*(u8(:,x_index(2)) - u8(:,end)))/dx; % [mol N2O/m2/day]
u10(u10<33)=33;
J_N2_mol = (D_N2(1)*(u10(:,x_index(1)) - u10(:,1))+D_N2(2)*(u10(:,x_index(2)) - u10(:,end)))/dx; % [mol N2/m2/day]
J_CO2 = J_CO2_mol*1e6*12/24 ; % [ug C/m2/h], converted from [mol CO2/m2/day]
J_N2O = J_N2O_mol*1e6*28/24 ; % [ug N/m2/h]
J_N2 = J_N2_mol*1e6*28/24 ; % [ug N/m2/h]