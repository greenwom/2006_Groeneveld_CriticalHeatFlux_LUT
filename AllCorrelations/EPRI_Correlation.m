function [CHF_EPRI,CHF_EPRI_Avg,L_EPRI,x_EPRI_local] = ...
    EPRI_Correlation(G,A_heated,A_test,L_heated,x_in,...
              h_fg,Pr,K_g,cwall,nu,toggle_vis,toggle_plotsave,plotfilepath)
% Prediction of the critical heat flux using the EPRI correlation as
% presented in 1982 EPRI Parametric Study of CHF Data Vol. 1-3
%
%       "Outputs"
% CHF_EPRI => Critical heat flux predicted by EPRI              [kW/m^2]
% CHF_EPRI_Avg => Average heat flux based on total heater power [kW/m^2]
% L_EPRI => Location of CHF event                               [m]
% x_EPRI_local => quality at location of predicted CHF          [-]
%
%       "Inputs"
% G => mass flux per subchannel                     [kg/m^2s]
% A_heated => total heated area per heater element  [m^2]
% A_test => flow area per subchannel                [m^2]
% L_heated => heated length per heater element      [m]
% x_in => inlet quality                             [-]
% h_fg => latent heat of vaporization               [J/kg]
% Pr => reduced pressure                            [-]
% K_g => grid spacer pressure loss coefficient [-]. Set default K_g = 1.
% cwall toggels cold wall effect correction factor -> 1/0 = on/off
% nu toggels nonuniform heat flux effect correction factor -> 1/0 = on/off
% toggle_vis => toggle visibility of the convergence plot. 'on'/'off'

if nu == 0 || nu == 1
else
    fprintf('Check input for nu. Must be 1 or 0.')
    return;
end
if cwall == 0 || cwall == 1
else
    fprintf('Check input for cwall. Must be 1 or 0.')
    return;
end

%EPRI Regressed Constants
P1 = 0.5328;
P2 = 0.1212;
P3 = 1.6151;
P4 = 1.4066;
P5 = -0.3040;
P6 = 0.4843;
P7 = -0.3285;
P8 = -2.0749;

%Convert from SI to English Units
G = G*7.3733812*10^(-4);                %[kg/m^2-s] to [Mlbm/ft^2-hr]
A_heated = A_heated*(100/(2.54*12))^2;  %[m^2] to [ft^2]
A_test = A_test*(100/(2.54*12))^2;      %[m^2] to [ft^2]
L_heated = L_heated*100/(2.54*12);      %[m] to [ft]
h_fg = h_fg*4.2992261*10^(-4);          %[J/kg] to [MBTU/Mlbm]
    
%EPRI Defined Variables
A = P1*Pr^(P2)*G^(P5+P7*Pr);
C = P3*Pr^(P4)*G^(P6+P8*Pr);

%Grid spacer effect correction factor
F_g = 1.3 - 0.3*K_g;

%Cold wall effect correction factors
if cwall == 1
    F_A = G^(0.1);
    F_C = 1.183*G^(0.1);
elseif cwall == 0
    F_A = 1;
    F_C = 1;
end

%Constants in Stern Document 'Technical Specification SLTS-76'
theta_0 = 0.8187458177;
theta_1 = 0.6812541823;
theta_2 = 2.436354311;

%Define position vector
x_n = 200;
xx = linspace(0,L_heated,x_n); %[ft]

%Initial Guess
P_CHF_Avg(1,:) = [1,2].*3.41214163*10^(-3);  %Heater power [kW] to [MBTU/hr]

Xratio=xx./L_heated;                 

%Solution convergence criteria
Error = 1;
tol = 0.0001;
itmax = 100;
iter = 0;
kk = 0;

fig_title = 'EPRI Correlation Convergence Plot';
figure_EPRIconv = figure('Name',fig_title,'NumberTitle','off','Visible',toggle_vis);
hold on;
while (Error > tol && itmax > iter)
    
    q_flux_local = zeros(x_n,2);
    q_CHF = zeros(x_n,2);
    iter = iter + 1;
    kk = kk + 1;
    
    for j = 1:2
    for i=1:x_n

        %Average heat flux per heater [MBTU/ft^2-hr]
        q_Flux_Avg = P_CHF_Avg(kk,j)/A_heated;
        
        %Definate integral (relative) of the chopped cosine power profile from
        %Xratio = 0..Xratio
        Power_Profile = theta_0 + theta_1*cos(2*theta_2*(Xratio(i) - 0.5));
        relative = theta_0*Xratio(i)+...
            0.5*theta_1/theta_2*sin(2*theta_2*(Xratio(i) - 0.5))+...
            0.5*theta_1/theta_2*sin(theta_2);
        
        %Integrated heat flux up to location xx [MBTU/ft^2-hr]
        q_flux_integrated = q_Flux_Avg*relative;
        
        %Local heat flux at location xx [MBTU/ft^2-hr]
        q_flux_local(i,j) = q_Flux_Avg*Power_Profile;
        
        %Nonuniform effect correction factor
        if nu == 1
            q_average = q_flux_integrated/xx(i);
            Y = q_average/q_flux_local(i,j);
            F_nu = 1 + (1-Y)/(1+G);
        elseif nu == 0
            F_nu = 1;
        end
        
        %Increase in quality due to heating up to location xx
        dx = q_flux_integrated*A_heated/(G*A_test*h_fg);
        
        %Local quality at location xx
        x_local(i) = x_in + dx;
        
        %EPRI Correlation for Critical Heat Flux
        q_CHF(i,j) = (A*F_A-x_in)/(C*F_g*F_C*F_nu + (x_local(i) - x_in)/q_flux_local(i,j));
    end 
    end
    
    error_1 = (q_CHF(:,1) - q_flux_local(:,1))./q_CHF(:,1);
    error_2 = (q_CHF(:,2) - q_flux_local(:,2))./q_CHF(:,2);
    [min_error_1,index_1] = min(error_1);
    [min_error_2,index_2] = min(error_2);
    
    %Calculate error and next step (Secant Method)
    gp = (min_error_2-min_error_1) / (P_CHF_Avg(kk,2)-P_CHF_Avg(kk,1));
    dPower = -min_error_2 / gp;
    P_CHF_Avg_New = P_CHF_Avg(kk,2) + dPower;
    
    % Update the error
    Error = abs(P_CHF_Avg_New-P_CHF_Avg(kk,2));
    P_CHF_Avg(kk+1,1) = P_CHF_Avg(kk,2);
    P_CHF_Avg(kk+1,2) = P_CHF_Avg_New;
        
    %[ft] to [m]
    xx_si = xx.*2.54*12/100;
    L_heated_si = L_heated*2.54*12/100;
    
    %[BTU/ft^2-hr] to [kW/m^2]
    q_EPRI = q_CHF(:,2)./(3.169983306*10^(-4));
    q_profile = q_flux_local(:,2)./(3.169983306*10^(-4));
    CHF_EPRI = q_EPRI(index_2);
    L_EPRI = xx_si(index_2);
        
    %Define plot limits
    minx = 0;
    maxx = L_heated_si;
    miny = 0;
    maxy = max(max(q_EPRI,q_profile));

    %Plot actual vs predicted heat flux as a function of position
    plot(xx_si,q_profile,'-k',xx_si,q_EPRI,'-r',[L_EPRI,L_EPRI],[miny,CHF_EPRI],'-k',[0 L_EPRI],[CHF_EPRI,CHF_EPRI],'-k');
    legend('Heater Flux','EPRI Flux');
    axis([minx maxx miny maxy]);
    xlabel('Heater Element Location [m]');
    ylabel('Heat Flux [kW/m^2]');
end
hold off;

if toggle_plotsave == 1
    FigSave = strcat(plotfilepath,'\',fig_title);
    print('-dpng','-r400',FigSave);
end

x_EPRI_local = x_local(index_2);
%Average heat flux [MBTU/ft^2-hr] to [kW/m^2]
CHF_EPRI_Avg = (P_CHF_Avg(kk,2)/(3.41214163*10^(-3)))/(A_heated/(100/(2.54*12))^2);

if iter == itmax
    fprintf('The maximum # of iterations was reached for EPRI\n');
end
end

