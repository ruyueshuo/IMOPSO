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

clc;
clear;
close all;
tic

global mpc_data 
[mpc_data]=loadcase('case39'); 
% mpc_data.bus(4,3) = mpc_data.bus(4,3) + 20;
 
%% Problem Definition

TestProblem=4;   % Set to 1, 2,  3

switch TestProblem
    case 1
        CostFunction=@(x) MyCost1(x);
        nVar=5;
        VarMin=-4;
        VarMax=4;
        
    case 2
        CostFunction=@(x) MyCost2(x);
        nVar=3;
        VarMin=-5;
        VarMax=5;
        
    case 3
        CostFunction=@(x) MyCost3(x);
        nVar=2;
        VarMin=0;
        VarMax=1;
    case 4
        CostFunction=@(x) MyCost4(x);

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
       
        %% get bus index lists of each type of bus
        [ref, pv, pq] = bustypes(mpc_data.bus, mpc_data.gen);
        
        %% generator info
        on = mpc_data.gen(:, GEN_STATUS) > 0;      %% which generators are on?
        gbus = mpc_data.gen(on , GEN_BUS);                %% what buses are they at?
        pv = intersect(pv, gbus);
%         [pvnum, ipv, ion] = intersect(pv, on);
        [pvnum, ipv, ion] = intersect( pv , mpc_data.gen(: , GEN_BUS));
        
        %% set up initial variables and bounds
%         Pg   = mpc_data.gen(:, PG);
        Pg   = mpc_data.gen(ion, PG);
%         Qg   = mpc_data.gen(ion, QG);
        Pmin = mpc_data.gen(ion, PMIN);
        Pmax = mpc_data.gen(ion, PMAX);
%         Qmin = mpc_data.gen(ion, QMIN);
%         Qmax = mpc_data.gen(ion, QMAX);
%         Pg   = mpc_data.gen(:, PG);
%         Qg   = mpc_data.gen(:, QG);
%         Pmin = mpc_data.gen(:, PMIN);
%         Pmax = mpc_data.gen(:, PMAX);
%         Qmin = mpc_data.gen(:, QMIN);
%         Qmax = mpc_data.gen(:, QMAX);
        Vg   = mpc_data.gen(ion, VG);
        for i = 1:length(ion)
            Vmin(i) = mpc_data.bus(ion(i), VMIN);
            Vmax(i) = mpc_data.bus(ion(i), VMAX);
        end
        
%         Vmin = mpc_data.bus(gbus(1):gbus(end), VMIN);
%         Vmax = mpc_data.bus(gbus(1):gbus(end), VMAX);
        %%
        % Obtain the number of decision variables
%         x = [Pg ; Qg];
        x = [Pg ; Vg];
        nVar = length(x);
        % Obtain the minimum possible value for each decision variable
%        xmin = [Pmin ; Qmin];
        xmin = [Pmin ; Vmin'];
        VarMin = xmin';
        % Obtain the maximum possible value for each decision variable
%        xmax = [Pmax ; Qmax]; 
        xmax = [Pmax ; Vmax']; 
        VarMax = xmax';
end

VarSize=[1 nVar];

VelMax=(VarMax-VarMin)/10;

%% MOPSO Settings

nPop=100;   % Population Size

nRep=500;   % Repository Size

MaxIt=200;  % Maximum Number of Iterations

phi1=2.05;
phi2=2.05;
phi=phi1+phi2;
chi=2/(phi-2+sqrt(phi^2-4*phi));

w=chi;              % Inertia Weight   %%% need to change %%%
wdamp=1;            % Inertia Weight Damping Ratio
c1=chi*phi1;        % Personal Learning Coefficient
c2=chi*phi2;        % Global Learning Coefficient

alpha=0.1;  % Grid Inflation Parameter

nGrid=10;   % Number of Grids per each Dimension

beta=4;     % Leader Selection Pressure Parameter

gamma=2;    % Extra (to be deleted) Repository Member Selection Pressure

%% Initialization

particle=CreateEmptyParticle(nPop);

for i=1:nPop
    particle(i).Velocity=0;
    particle(i).Position=unifrnd(VarMin,VarMax,VarSize);

    particle(i).Cost=CostFunction(particle(i).Position);

    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
end

particle=DetermineDomination(particle);

rep=GetNonDominatedParticles(particle);

rep_costs=GetCosts(rep);
G=CreateHypercubes(rep_costs,nGrid,alpha);

for i=1:numel(rep)
    [rep(i).GridIndex,rep(i).GridSubIndex]=GetGridIndex(rep(i),G);
end
    
%% MOPSO Main Loop

for it=1:MaxIt
    for i=1:nPop
        rep_h=SelectLeader(rep,beta);

        particle(i).Velocity=w*particle(i).Velocity ...
                             +c1*rand*(particle(i).Best.Position - particle(i).Position) ...
                             +c2*rand*(rep_h.Position -  particle(i).Position);

        particle(i).Velocity=min(max(particle(i).Velocity,-VelMax),+VelMax);

        particle(i).Position=particle(i).Position + particle(i).Velocity;

        flag=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(flag)=-particle(i).Velocity(flag);
        
        particle(i).Position=min(max(particle(i).Position,VarMin),VarMax);

        particle(i).Cost=CostFunction(particle(i).Position);

        if Dominates(particle(i),particle(i).Best)
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
        elseif ~Dominates(particle(i).Best,particle(i))
            if rand<0.5
                particle(i).Best.Position=particle(i).Position;
                particle(i).Best.Cost=particle(i).Cost;
            end
        end

    end
    
    particle=DetermineDomination(particle);
    nd_particle=GetNonDominatedParticles(particle);
    
    rep=[rep
         nd_particle];
    
    rep=DetermineDomination(rep);
    rep=GetNonDominatedParticles(rep);
    
    for i=1:numel(rep)
        [rep(i).GridIndex,rep(i).GridSubIndex]=GetGridIndex(rep(i),G);
    end
    
    if numel(rep)>nRep
        EXTRA=numel(rep)-nRep;
        rep=DeleteFromRep(rep,EXTRA,gamma);
        
        rep_costs=GetCosts(rep);
        G=CreateHypercubes(rep_costs,nGrid,alpha);
        
    end
   
    disp(['Iteration ' num2str(it) ': Number of Repository Particles = ' num2str(numel(rep))]);
    
    w=w*wdamp;
end

%% Results

disp(['Calculation Used ',num2str(toc),' s']);

costs=GetCosts(particle);
rep_costs=GetCosts(rep);
M = size(rep_costs);

%% fuzzy clustering 
for i = 1:M(1)
    for j = 1:M(2)
%         u(i,j) = (rep_costs(i,j)-min(rep_costs(i,:)))/(max(rep_costs(i,:))-min(rep_costs(i,:)));
         u(i,j) = (max(rep_costs(i,:))-rep_costs(i,j))/(max(rep_costs(i,:))-min(rep_costs(i,:)));
    end
end

for j = 1:M(2)
    w(j) = sum(u(:,j))/(sum(sum(u)));  %% need to change namda
end

idx = find(w==(max(w)));
rep_costs(:,idx)


%% TOPSIS
%  f = rep_costs./sqrt(sum(rep_costs.^2,2));
%  for i = 1:M(1)
%      fmax(i) = max(f(i,:));
%      fmin(i) = min(f(i,:));
%  end
%  fmax = fmax';
%  fmin = fmin';
%  namda = ones(M(1),1)./M(1);  %% need to change
%  dmax = sqrt(sum(namda.*(f-fmax)).^2);
%  dmin = sqrt(sum(namda.*(f-fmin)).^2);
%  d = dmax./(dmax+dmin);
%  idx = find(d==(max(d)));
%  rep_costs(:,idx)

%% Visualize
% The following is used to visualize the result if objective space
% dimension is visualizable.

if M(1) == 2
%     plot(costs(1,:),costs(2,:),'b.');
%     hold on
    plot(rep_costs(1,:),rep_costs(2,:),'r*');
    hold on
%     plot(rep_costs(1,idx),rep_costs(2,idx),'ks');
%     legend('Main Population','Repository');
    x1=xlabel('Power Loss');        
    x2=ylabel('Reactive Power Equilibrium Degree');     
elseif M(1) ==3
    plot3(rep_costs(1,:),rep_costs(2,:),rep_costs(3,:),'r*');
    hold on
    plot3(rep_costs(1,idx),rep_costs(2,idx),rep_costs(3,idx),'bs');
    x1=xlabel('Economy Index');        
    x2=ylabel('Safety Index');      
    x3=zlabel('Quality Index'); 
end
hold on

% figure;      

save 'C:\Users\Feng Da\Documents\MATLAB\AGCAVC-MOPSO\results\rep_costs.mat' rep_costs 
save 'C:\Users\Feng Da\Documents\MATLAB\AGCAVC-MOPSO\results\rep_costs.txt' rep_costs -ASCII
rep_cell = struct2cell(rep);
rep_position = rep_cell(1,:);
save 'C:\Users\Feng Da\Documents\MATLAB\AGCAVC-MOPSO\results\rep_position.mat' rep_position 
% x1=xlabel('NSGA-ii');        
% x2=ylabel('MOPSO'); 
%%
cost_best_position = cell2mat(rep_position(idx));

for i = 1:nVar/2
    mpc.gen(i,PG) = cost_best_position(i);
    mpc.gen(i,VG) = cost_best_position(nVar/2+i);
end 
%savecase('C:\Users\Feng Da\Documents\MATLAB\AGCAVC\case_test', mpc);
mpc_best  = mpc;
%%
%[result, success] = runpf('case_test');
%[baseMVA, bus, gen, branch, success, et] = runpf(...)
[result_best, success_best] = runpf(mpc_best);