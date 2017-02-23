function [CHF_CISEGE,CHF_CISEGE_Avg,L_CISEGE,x_CISEGE_local] = ...
    CISEGE_Correlation(A_heated,A_test,G,h_fg,Inlet_Press,L_heated,x_in,...
    toggle_vis,toggle_plotsave,plotfilepath)
% Prediction of the critical heat flux using the CISE-GE Critical Quatlity
% Correlation as presented in the TRACE Theory Manual V5.0 by the USNRC.
% The correlation was taken from J. A. Borkowski and N. L. Wade,
% 'TRAC-BF1/MOD1 Models and Correlations', NUREG/CR-4391, Section 4.2.12,
% 1992. The mentioned version is a modified version of the CISE-GE
% correlation. The original CISE-GE correlation was provided by GE as a
% public correlation by GE for convienience for outside design reviewers
% (GE report 1975). The actual GEXL correlation is proprietary.
%
%       "Outputs"
% CHF_CISEGE => Critical heat flux predicted b the LUT             [kW/m^2]
% CHF_CISEGE_Avg => Average heat flux based on total heater power  [kW/m^2]
% L_CISEGE => Location of CHF event                                [m]
% x_CISEGE_local => quality at location of predicted CHF           [-]
%
%       "Inputs"
% A_heated => total heated area per heater element  [m^2]
% A_test => flow area per subchannel                [m^2]
% G => mass flux per subchannel                     [kg/m^2s]
% h_fg => latent heat of vaporization               [J/kg]
% Inlet_Press => inlet pressure                     [Pa]
% L_heated => heated length per heater element      [m]
% x_in => inlet quality                             [-]

h_fg = h_fg/1000;   % [J/kg] to [kJ/kg]

% Define Correlation parameters
%7.37*10^(-4) is a correction factor for SI to English units
%Constants in B have been converted from [in] to [m]
A = 1.055 - 0.013*((Inlet_Press-4.137*10^6)/(2.758*10^6))^2 ...
    - 1.233*(7.37*10^(-4)*G) + 0.907*(7.37*10^(-4)*G)^2  ...
    - 0.285*(7.37*10^(-4)*G)^3;
    
B = 0.457 + 2.003*(7.37*10^(-4)*G) - 0.901*(7.37*10^(-4)*G)^2;
R_f = 1;    %Radial peaking factor

%Constants in Stern Document 'Technical Specification SLTS-76'
theta_0 = 0.8187458177;
theta_1 = 0.6812541823;
theta_2 = 2.436354311;

%Define position vector
x_n = 200;
xx = linspace(0,L_heated,x_n); %[m]

%Initialize arrays
x_local = zeros(x_n,1);

%Initial Guess
P_CHF_Avg(1,:) = [30,31];  %Heater power [kW]

Xratio=xx./L_heated;                 

%Solution convergence criteria
Error = 1;
tol = 0.001;
itmax = 100;
iter = 0;
kk = 0;

fig_title = 'CISE-GE Correlation Convergence Plot';
figure_CISEGEconv = figure('Name',fig_title,'NumberTitle','off','Visible',toggle_vis);
hold on;
while (Error > tol && itmax > iter)
    
    x_local = zeros(x_n,2);
    x_CHF = zeros(x_n,2);
    
    iter = iter + 1;
    kk = kk + 1;
    
    for j = 1:2
    for i=1:x_n

        %Average heat flux per heater [kW/m^2]
        q_Flux_Avg = P_CHF_Avg(kk,j)/A_heated;
        
        %Definate integral (relative) of the chopped cosine power profile
        %from Xratio = 0..Xratio
        Power_Profile = theta_0 + theta_1*cos(2*theta_2*(Xratio(i) - 0.5));
        relative = theta_0*Xratio(i)+...
            0.5*theta_1/theta_2*sin(2*theta_2*(Xratio(i) - 0.5))+...
            0.5*theta_1/theta_2*sin(theta_2);
        
        %Integrated heat flux up to location xx [kW/m^2]
        q_flux_integrated = q_Flux_Avg*relative;
                       
        %Increase in quality due to heating up to location xx
        dx = q_flux_integrated*A_heated/(G*A_test*h_fg);
        
        %Local quality at location xx
        x_local(i,j) = x_in + dx;
        
        %CISE-GE correlation
        if x_local(i,j) > 0
            L_B = BoilingLength(A_test,G,h_fg,L_heated,P_CHF_Avg(kk,j),x_in,Xratio(i));
            x_CHF(i,j) = A*L_B/(B + L_B)*(1.24/R_f);
        else
            L_B = 0;
            x_CHF(i,j) = .25 + j/4;
        end
            
    end 
    end
    
    error_1 = (x_CHF(:,1) - x_local(:,1))./x_CHF(:,1);
    error_2 = (x_CHF(:,2) - x_local(:,2))./x_CHF(:,2);
    [min_error_1,index_1] = min(error_1);
    [min_error_2,index_2] = min(error_2);
    
    %Calculate error and next step (Secant Method)
    gp = (min_error_2-min_error_1) / (P_CHF_Avg(kk,2)-P_CHF_Avg(kk,1));
    dPower = -min_error_2 / gp;
    P_CHF_Avg_New = P_CHF_Avg(kk,2) + dPower/4;
    
    % Update the error
    Error = abs(P_CHF_Avg_New-P_CHF_Avg(kk,2));
    P_CHF_Avg(kk+1,1) = P_CHF_Avg(kk,2);
    P_CHF_Avg(kk+1,2) = P_CHF_Avg_New;
            
    i=1:x_n;
    Power_Profile(i) = theta_0 + theta_1*cos(2*theta_2*(Xratio(i) - 0.5));
    
    %Local heat flux at location xx [kW/m^2]
    q_flux_local = q_Flux_Avg.*Power_Profile;
    
    %[kW/m^2],[m] and store information at CHF location
    x_CISEGE = x_CHF(:,2);
    x_profile = x_local(:,2);
    CHF_CISEGE = q_flux_local(index_2);
    L_CISEGE = xx(index_2);
        
    %Define plot limits
    minx = 0;
    maxx = L_heated;
    miny = 0;
    x_CISEGE(x_CISEGE==0.25) = 0;
    x_CISEGE(x_CISEGE==0.75) = 0;
    if x_CISEGE(:,:) <= 0
        maxy = 0.25;
    else
    maxy = ceil(max(max(x_CISEGE,x_profile))*10)/10;
    end

    %Plot actual vs predicted heat flux as a function of position
    plot(xx,x_profile,'-k',xx,x_CISEGE,'-r');%,[L_CISEGE,L_CISEGE],[miny,CHF_CISEGE],'-k',[0 L_CISEGE],[CHF_CISEGE,CHF_CISEGE],'-k');
    legend('Profile Quality','CISE-GE Quality');
    axis([minx maxx miny maxy]);
    xlabel('Heater Element Location [m]');
    ylabel('Quality [-]');
end
hold off;

if toggle_plotsave == 1
    FigSave = strcat(plotfilepath,'\',fig_title);
    print('-dpng','-r400',FigSave);
end

x_CISEGE_local = x_local(index_2,2);

%Average heat flux [kW/m^2]
CHF_CISEGE_Avg = P_CHF_Avg(kk,2)/A_heated;
if iter == itmax
    fprintf('The maximum # of iterations was reached for CISE-GE\n');
end
end

function L_Boiling = BoilingLength(A_test,G,h_fg,L_htd,P_rod_avg,x_in,XLcurrent)
%       "Outputs"
% L_Boiling => boiling length                   [m]
%
%       "Inputs"
% A_heated => element heated area               [m^2]
% A_test => subchannel flow area                [m^2]
% G => mass flux                                [kg/m^2-s]
% h_fg => inlet subcooling                      [kJ/kg]
% L_htd => total heated length                  [m]
% P_rod_avg => total heater power               [kW]
% x_in => inlet quality                         [-]
% XLcurrent => current location xx(i)/L_htd     [-]

h_sub = -x_in*h_fg; %inlet subcooling [kJ/kg]
%Power required to bring subcooled inlet to a x_quality = 0
P_zero = h_sub.*G.*A_test;

%Constants in Stern Document 'Technical Specification SLTS-76'
theta_0 = 0.8187458177;
theta_1 = 0.6812541823;
theta_2 = 2.436354311;

L_zero = [0.25 0.3]; %[m]
Error = 1;
tol = 0.001;
itmax = 100;
iter = 0;
Xratio = L_zero./L_htd;
while (Error > tol && iter < itmax)
    iter = iter + 1;
    
    %Guess location 1
    %Integral of Profile from 0 to L_CHF and normalized to 1
    rel_zero_1 = theta_0*Xratio(1)+...
        0.5*theta_1/theta_2*sin(2*theta_2*(Xratio(1) - 0.5))+...
        0.5*theta_1/theta_2*sin(theta_2);

    %Integrated Power up to location of x_quality = 0
    P_int_zero_1 = P_rod_avg*rel_zero_1;
    Error_1 = P_zero - P_int_zero_1;
    
    %Guess location 2
    rel_zero_2 = theta_0*Xratio(2)+...
        0.5*theta_1/theta_2*sin(2*theta_2*(Xratio(2) - 0.5))+...
        0.5*theta_1/theta_2*sin(theta_2);
    
    %Integrated Power up to location of x_quality = 0
    P_int_zero_2 = P_rod_avg*rel_zero_2;
    Error_2 = P_zero - P_int_zero_2;
    
    %Calculate error and next step
    gp = (Error_2-Error_1) / (L_zero(2)-L_zero(1));
    dx = -Error_2 / gp;
    L_zero(3) = L_zero(2) + dx;
    
    % Update the error
    Error = abs(L_zero(3)-L_zero(2));
    L_zero(1) = L_zero(2);
    L_zero(2) = L_zero(3);
    Xratio = L_zero./L_htd;
end

XLratio = [L_zero(3)/L_htd XLcurrent];  %Boiling Limits
L_Boiling = L_htd*(XLratio(2)-XLratio(1)); %Boiling length

end