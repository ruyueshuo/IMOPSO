%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MATLAB Code for                                                  %
%                                                                   %
%  Multi-Objective Particle Swarm Optimization (MOPSO)              %
%  Version 1.0 - Feb. 2011                                          %
%                                                                   %
%  According to:                                                    %
%  Carlos A. Coello Coello et al.,                                  %
%  "Handling Multiple Objectives with Particle Swarm Optimization," %
%  IEEE Transactions on Evolutionary Computation, Vol. 8, No. 3,    %
%  pp. 256-279, June 2004.                                          %
%                                                                   %
%  Developed Using MATLAB R2009b (Version 7.9)                      %
%                                                                   %
%  Programmed By: S. Mostapha Kalami Heris                          %
%                                                                   %
%         e-Mail: sm.kalami@gmail.com                               %
%                 kalami@ee.kntu.ac.ir                              %
%                                                                   %
%       Homepage: http://www.kalami.ir                              %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

function z=MyCost3(x)

%     n=numel(x);
% 
%     z=[0 0];
%     
%     z(1)=1-exp(-sum((x-1/sqrt(n)).^2));
%     
%     z(2)=1-exp(-sum((x+1/sqrt(n)).^2));
    
%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus; 
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%Input data of IEEE 30 bus test system
%[baseMVA, bus, gen, branch, areas, gencost]=loadcase('case_ieee30'); 
%[mpc]=loadcase('case118'); 
global mpc_data
global mpc_change
mpc = mpc_data;

nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ng   = size(mpc.gen, 1);    %% number of dispatchable injections

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(mpc.bus, mpc.gen);

%% generator info
on = mpc.gen(:, GEN_STATUS) > 0;      %% which generators are on?
gbus = mpc.gen(on , GEN_BUS);                %% what buses are they at?

pv = intersect(pv, gbus);
% npv = size(pv , 1);
[pvnum, ipv, ion] = intersect( pv , mpc_data.gen(: , GEN_BUS));
%%
% for i = 1:ng
%     mpc.gen(i,PG) = x(i);
%     mpc.gen(i,QG) = x(ng+i);
% end 
npv =length(ion);
for i = 1:npv
    mpc.gen(ion(i),PG) = x(i);
    mpc.gen(ion(i),QG) = x(npv+i);
end 
%savecase('C:\Users\Feng Da\Documents\MATLAB\AGCAVC\case_test', mpc);
%save mpc
mpc_change = mpc;

%% Power Flow Calculation
%[result, success] = runpf('case_test');
%[baseMVA, bus, gen, branch, success, et] = runpf(...)
[result, success] = runpf(mpc_change);

%% objective
z = [];
%% Objective function one-system losses calculation
if success == 1
	total_system_real_losses = sum(real(get_losses(result)));
else 
	total_system_real_losses = 1000000;
end
z(1) = total_system_real_losses;

%% Objective function two-reactive power equilibrium degree calculation
reactive_power_output_level = zeros(ng,1);
if success == 1
    for i = 1:ng
	reactive_power_output_level(i) = (result.gen(i,QG) - result.gen(i,QMIN))/...
        (result.gen(i,QMAX) - result.gen(i,QMIN));
    end
    reactive_power_equilibrium_degree = var(reactive_power_output_level);
else
    reactive_power_equilibrium_degree = 1000000;
end
z(2) = reactive_power_equilibrium_degree;

%% Objective function three-quality index
active_power_output_level = zeros(nb,1);
if success == 1
    for i = 1:nb
        active_power_output_level(i) = result.bus(i,VM)-...
            (result.bus(i,VMAX) + result.bus(i,VMIN))/2;
    end
    active_power_equilibrium_degree = sum(abs((active_power_output_level)));
else
    active_power_equilibrium_degree = 1000000;
end
z(3) = active_power_equilibrium_degree;

end