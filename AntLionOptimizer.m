%___________________________________________________________________%
%      Ant Lion Optimizer (ALO) source codes demo version 1.0       %
%                                                                   %
%                Developed in MATLAB R2011b(7.13)                   %
%                                                                   %
%           Author and programmer: Seyedali Mirjalili               %
%                                                                   %
%                   Edithor by Riki Khomarudin                      %
%              Muhammadiyah University Of Yogyakarta                %
%                                                                   %
%                   e-Mail:ali.mirjalili@gmail.com                  %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%   Homepage: http://www.alimirjalili.com                           %
%                                                                   %
%   Main paper:                                                     %
%                                                                   %
%   S. Mirjalili, The Ant Lion Optimizer                            %
%   Advances in Engineering Software , in press,2015                %
%   DOI: http://dx.doi.org/10.1016/j.advengsoft.2015.01.010         %
%                                                                   %
%___________________________________________________________________%

% You can simply define your cost in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of generations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers

% To run ALO: [Best_score,Best_pos,cg_curve]=ALO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)

%function [Elite_antlion_fitness,Elite_antlion_position,Convergence_curve]=ALO(N,Max_iter,lb,ub,dim,fobj)
clear;
%clc;
tic;
global scenario alpha_t CENS ICdg MCdg minVoltage maxVoltage

scenario_=[1 2 3];
scenario=scenario_(3);

minVoltage_=[0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99];
minVoltage=minVoltage_(6);

maxVoltage_=[1.00 1.05];
maxVoltage=maxVoltage_(2);

alpha_t_=[0.7 1 1.05]; %Faktor keekonomian yang mengubah laju aliran kas masa depan untuk nilai sekarang
alpha_t=alpha_t_(3);

CENS_=[0 8 10 20 30 40 50 60 78 90 100];
CENS=CENS_(4);

 %C_Purc=[1700 1800 1900 2000 2100 2200 2300 2400 2500];
 %CP=C_Purc(4);

MCdg=[26.04 33.00 36.48];
%MCdg=MCdg_(3);

%Ct_=[0.1 0.2 0.3]; %Electricity market price
%Ct=Ct_(1);

ICdg_=[100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500];
ICdg=ICdg_(6);

%Run base load flow to read data, create distribution structure, and finevoltage
BaseLoadFlow;
br;
no;
LF;
brlenght=1;
PF=0.85;        % Power factor, option: 0.85; 0.95; 1.
%PF2=0.9;
%PF3=1;
R;
X;
P;
Q;
MVAb;
Vb;
adjcb;
adjb;
oldloss=PL;
%priceENS=CENS;
minVoltage;
maxVoltage;
plot(Voltage,'--rs','LineWidth',2,...
                'MarkerEdgeColor','g',...
                'MarkerFaceColor','g',...
                'MarkerSize',4)
set(gca,'XTick',1:1:no);
title('Voltage profile');
xlabel('Bus');
ylabel('Voltage magnitude (p.u.)');
hold all

% Dimension of the search variables
%ObjectiveFunction=@func;

lb=[1 40];
ub=[33  16800];
dimension=2;
obj_no=2;
scenario=2;               % Pilih scenario 1,2,atau 3

%lb&&ub(2)=populasi_DG;
dgmax=ub(2);
dgmin=lb(2);
populasi_DG=(dgmax-dgmin)*rand()+dgmin;

Elite_antlions_position=populasi_DG;
%iterasi parameter
Max_iter=99;  % Jumlah penrulangan
N=100;    % jumlah populasi

% Initialize the positions of antlions and ants
antlion_position=initialization(N,dimension,ub,lb);
ant_position=initialization(N,dimension,ub,lb);

% Initialize variables to save the position of elite, sorted antlions, 
% convergence curve, antlions fitness, and ants fitness
Sorted_antlions=zeros(N,dimension);

Elite_antlion_position=zeros(1,dimension);
Elite_antlion_fitness=inf;

Convergence_curve=zeros(1,Max_iter);

antlions_fitness=zeros(1,N);
ants_fitness=zeros(1,N);

% Calculate the fitness of initial antlions and sort them
for i=1:size(antlion_position,1)
antlions_fitness(1,i)=BaseLoadFlow2(antlion_position(i,:),br,no,brlenght,PF,R,X,P,Q,MVAb,Vb,adjcb,adjb,oldloss); 
end

[sorted_antlion_fitness,sorted_indexes]=sort(antlions_fitness);
    
for newindex=1:N
     Sorted_antlions(newindex,:)=antlion_position(sorted_indexes(newindex),:);
end
    
Elite_antlion_position=Sorted_antlions(1,:);
Elite_antlion_fitness=sorted_antlion_fitness(1);


% Main loop start from the second iteration since the first iteration 
% was dedicated to calculating the fitness of antlions
Current_iter=2; 
while Current_iter<Max_iter+1
    
Convergence_curve(1,1)=Elite_antlion_fitness;

    % This for loop simulate random walks
    for i=1:size(ant_position,1)
        % Select ant lions based on their fitness (the better anlion the higher chance of catching ant)
        Rolette_index=RouletteWheelSelection(1./sorted_antlion_fitness);
        if Rolette_index==-1  
            Rolette_index=1;
        end
      
        % RA is the random walk around the selected antlion by rolette wheel
        RA=Random_walk_around_antlion(dimension,Max_iter,lb,ub, Sorted_antlions(Rolette_index,:),Current_iter);
        
        % RA is the random walk around the elite (best antlion so far)
        [RE]=Random_walk_around_antlion(dimension,Max_iter,lb,ub, Elite_antlion_position(1,:),Current_iter);
        
        ant_position(i,:)= (RA(Current_iter,:)+RE(Current_iter,:))/2; % Equation (2.13) in the paper          
    end
    
    for i=1:size(ant_position,1)  
        
        % Boundar checking (bring back the antlions of ants inside search
        % space if they go beyoud the boundaries
        Flag4ub=ant_position(i,:)>ub;
        Flag4lb=ant_position(i,:)<lb;

ant_position(i,:)=(ant_position(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
ants_fitness(1,i)=BaseLoadFlow2(ant_position(i,:),br,no,brlenght,PF,R,X,P,Q,MVAb,Vb,adjcb,adjb,oldloss);        
       
    end
    
    % Update antlion positions and fitnesses based of the ants (if an ant 
    % becomes fitter than an antlion we assume it was cought by the antlion  
    % and the antlion update goes to its position to build the trap)
    double_population=[Sorted_antlions;ant_position];
    double_fitness=[sorted_antlion_fitness ants_fitness];
        
    [double_fitness_sorted, I]=sort(double_fitness);
    double_sorted_population=double_population(I,:);
        
    antlions_fitness=double_fitness_sorted(1:N);
    Sorted_antlions=double_sorted_population(1:N,:);
        
    % Update the position of elite if any antlinons becomes fitter than it
    if antlions_fitness(1)<Elite_antlion_fitness 
        Elite_antlion_position=Sorted_antlions(1,:);
        Elite_antlion_fitness=antlions_fitness(1);
    end
      
    % Keep the elite in the population
    Sorted_antlions(1,:)=Elite_antlion_position;
    antlions_fitness(1)=Elite_antlion_fitness;
  
    % Update the convergence curve
    Convergence_curve(Current_iter)=Elite_antlion_fitness;

    % Display the iteration and best optimum obtained so far
    if mod(Current_iter,50)==0
        display(['At iteration ', num2str(Current_iter), ' the elite fitness is ', num2str(Elite_antlion_fitness)]);
        if round(Current_iter/100)==Current_iter/100,
        dispElite_antlion_fitness=[round(Elite_antlion_fitness(1)) Elite_antlion_fitness(2)];
         disp(' ');
         disp(['Iteration: ',num2str(Current_iter)]);
         disp(['Best solution: ',num2str(dispElite_antlion_fitness)]);
         disp(['Elite_antlion_fitness: ',num2str(Elite_antlion_position)]);
        end
    end
    %Convergence_curve(Current_iter+1,1)=Elite_antlion_fitness;
    Current_iter=Current_iter+1; 
end

% Output/display
dispElite_antlion_position=[round(Elite_antlion_position(1)) (Elite_antlion_position(2))];
disp(' ');
disp(['Total number of evaluations: ',num2str(Max_iter*N)]);
disp(['Best solution= ',num2str(dispElite_antlion_position)]);
disp(['Elite_antlion_fitness: ',num2str(Elite_antlion_fitness)]);

% Run Load Flow with Optimal DG size and Placement

OV=Elite_antlion_position;
OptimalLoadFlow;

PV=load('Photovoltaic2.m');

nobus=round(Elite_antlion_position(1,1));
total_unit_PV_based_DG1=Elite_antlion_position(1,2)/PV(nobus,2);

total_PV_based_DG=ceil(total_unit_PV_based_DG1);
vkec=13;


plot(Voltage,'--ks','LineWidth',2,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','b',...
                'MarkerSize',3);
                 legend('before DG','after DG',1);
                 grid on;
                 
                                  
figure
plot(Convergence_curve)
set(gca,'XTick',0:10:Max_iter)
title('Convergence curve');
xlabel('Number of iteration');
ylabel('Rugi Daya (kW)');

toc






