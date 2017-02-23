function [CHF_LUT,CHF_LUT_Avg,L_LUT,x_LUT_local,K_Factors] = ...
    Groeneveld_LUT(G,Inlet_Press,x_in,A_heated,A_test,L_heated,h_fg,...
    toggle_K,D_hyd,D_htr,Pitch,K_g,rho_fsat,rho_gsat,toggle_vis,...
    toggle_plotsave,plotfilepath)
% Prediction of the critical heat flux using the 2006 Groeneveld Look-Up
% Table as presented in 2006 CHF Look-Up Table, Nuclear Engineering and
% Design 237, pp. 190-1922.
%
% This file requires an additional file 2006LUTdata.txt
%
%       "Outputs"
% CHF_LUT => Critical heat flux predicted b the LUT             [kW/m^2]
% CHF_LUT_Avg => Average heat flux based on total heater power  [kW/m^2]
% L_LUT => Location of CHF event                                [m]
% x_LUT_local => quality at location of predicted CHF           [-]
% K => the final K factors used to obtain CHF_LUT               [-]
%
%       "Inputs"
% G => mass flux per subchannel                     [kg/m^2s]
% Inlet_Press => inlet pressure                     [Pa]
% x_in => inlet quality                             [-]
% A_heated => total heated area per heater element  [m^2]
% A_test => flow area per subchannel                [m^2]
% L_heated => heated length per heater element      [m]
% h_fg => latent heat of vaporization               [J/kg]
%
% LUT Specific Inputs
% toggle_K => 8x1 vector of 1(on) or 0(off) indicating which correction
%            factors are to be used to calculate the 'corrected' heat flux.
% D_hyd => hydraulic diameter of subchannel         [m]
% D_htr => heater element diameter                  [m]
% Pitch => center distance between heating elements [m]
% K_g => pressure loss coefficient of spacer        [-]
% rho_fsat => fluid saturation density              [kg/m^3]
% rho_gsat => vapor saturation density              [kg/m^3]
% toggle_vis => toggle visibility of convergence plot ['on'/'off']

%       "Inputs for Correction Factors - K(1-8)"
% Correction factors definitions can be found in IAEA TecDoc 1203.
%
% If a K factor is not to be used a value of 1 can simply be defined for
% the applicable variables associated with that K factor.
%
% !!NOTE!! - If K(3) is used K(4) = 1. This is because K(3) already
% includes the heated length effect. This condition is handled in the code.
%
% Definition of correction factors K:
%   K(1) => Rod Diameter Factor
%   K(2) => Bundle Geometry Factor
%   K(3) => Mid-Plane Spacer Factor
%   K(4) => Heated-Length Factor
%   K(5) => Axial Flux Distribution Factor
%   K(6) => Radial/Circumferential Flux Distribution Factor
%   K(7) => Flow-Orientation Factor
%   K(8) => Vertical Low-Flow Factor

if toggle_K(3) == 1;
    toggle_K(4) = 0;
end

% **********************

% Load 2006 LUT for interpolation based on iteration

% Pressure range [MPa] from 2006 LUT
P_lut = [0.10,0.30,0.50,1.0,2.0,3.0,5.0,7.0,...
    10.0,12.0,14.0,16.0,18.0,20.0,21.0];

% Mass Flux range [kg/m^2-s] from 2006 LUT;
G_lut = [0,50,100,300,500,750,1000,1500,2000,2500,3000,3500,4000,4500,...
    5000,5500,6000,6500,7000,7500,8000];

% Quality range from 2006 LUT
x_lut = [-0.50,-0.40,-0.30,-0.20,-0.15,-0.10,-0.05,0.00,0.05,0.10,0.15,...
    0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80,0.90,1.00];

% Import the 2006 LUT from the text file
LUT_array = importdata('2006LUTdata.txt','\t');

% Convert the imported array into a (MxNxQ) where:
%   M is number of mass flux divisions
%   N is number of quality divisions
%   Q is number of pressure divisions

i_stop = 0;

lenG = length(G_lut);
lenx = length(x_lut);
lenP = length(P_lut);
LUT_Data = zeros(lenG,lenx,lenP);
for i=1:lenP
    i_start = i_stop + 1;
    i_stop = i_stop + lenG;
    index = i_start:i_stop;
    LUT_Data(:,:,i) = LUT_array(index,:);
end
% The 2006LUT is now ready to be called as necessary using the MATLAB
% built-in function => CHF = interpn(G_lut,x_lut,P_lut,LUT_Data,G,x,P)

% **********************

Inlet_Press = Inlet_Press/10^6; % Convert [Pa] to [MPa]
h_fg = h_fg/1000;               % Convert [J/kg] to [kJ/kg]

%Constants in Stern Document 'Technical Specification SLTS-76'
theta_0 = 0.8187458177;
theta_1 = 0.6812541823;
theta_2 = 2.436354311;

%Define position vector
x_n = 200;
xx = linspace(0,L_heated,x_n); %[m]

%Initialize arrays
K = ones(8,x_n);
x_local = zeros(x_n,1);

%Initial Guess
P_CHF_Avg(1,:) = [30,31];  %Heater power [kW]

Xratio=xx./L_heated;                 

%Solution convergence criteria
Error = 1;
tol = 0.0001;
itmax = 100;
iter = 0;
kk = 0;

fig_title = '2006 LUT Correlation Convergence Plot';
figure_LUTconv = figure('Name',fig_title,'NumberTitle','off','Visible',toggle_vis);
hold on;
while (Error > tol && itmax > iter)
    
    q_flux_local = zeros(x_n,2);
    q_CHF = zeros(x_n,2);       q_CHF_base = zeros(x_n,2);
    
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
        
        %Local heat flux at location xx [kW/m^2]
        q_flux_local(i,j) = q_Flux_Avg*Power_Profile;
                
        %Increase in quality due to heating up to location xx
        dx = q_flux_integrated*A_heated/(G*A_test*h_fg);
        
        %Local quality at location xx
        x_local(i) = x_in + dx;
               
        q_CHF_base(i,j) = interpn(G_lut,x_lut,P_lut,LUT_Data,...
                                  G,x_local(i),Inlet_Press);
   
        if toggle_K(1) == 1;
            K(1,i) = K_1_Factor(D_hyd);
        end
        if toggle_K(2) == 1;
            K(2,i) = K_2_Factor(D_htr,Pitch,x_local(i));
        end
        if toggle_K(3) == 1;
            Grid_Location = [1.0,1.5,2.0]; %Location of grid spacers [m]

            %Define distance from upstream grid spacer L_sp [m]
            if xx(i) >= Grid_Location(1) && xx(i) < Grid_Location(2)
                L_sp = xx(i) - Grid_Location(1);
            elseif xx(i) >= Grid_Location(2) && xx(i) < Grid_Location(3)
                L_sp = xx(i) - Grid_Location(2);
            elseif xx(i) >= Grid_Location(3)
                L_sp = xx(i) - Grid_Location(3);
            else
                L_sp = 100;
            end

            K(3,i) = K_3_Factor(D_hyd,G,K_g,L_sp);
        end
        if toggle_K(4) == 1;
            K(4,i) = K_4_Factor(D_hyd,L_heated,rho_fsat,rho_gsat,x_local(i));
        end
        if toggle_K(5) == 1;
            if x_local(i) > 0
            q_local = q_flux_local(i,j);
            q_BLA = BoilingLengthAvg(A_heated,A_test,G,h_fg,L_heated,...
                                     P_CHF_Avg(kk,j),x_in,Xratio(i));
            K(5,i) = K_5_Factor(q_BLA,q_local,x_local(i));
            else
            K(5,i) = 1;
            end
        end
%         if toggle_K(6) == 1;
%             %Not used in this research.
%             K(6,i) = K_6_Factor(q_rc_avg,q_rc_max,x_local);
%         end
%         if toggle_K(7) == 1;
%             %Not used in this research.
%             K(7,i) = K_7_Factor(f_L,G,grav,rho_fsat,rho_gsat,x_local);
%         end
%         if toggle_K(8) == 1;
%             %Not used in this research.
%             K(8,i) = K_8_Factor(G,rho_fsat,rho_gsat,x_local);
%         end

        q_CHF(i,j) = q_CHF_base(i,j)*prod(K(:,i));
        
    end 
    end
    
    error_1 = (q_CHF(:,1) - q_flux_local(:,1))./q_CHF(:,1);
    error_2 = (q_CHF(:,2) - q_flux_local(:,2))./q_CHF(:,2);
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
            
    %[kW/m^2]
    q_LUT = q_CHF(:,2);
    q_profile = q_flux_local(:,2);
    CHF_LUT = q_LUT(index_2);
    L_LUT = xx(index_2);
        
    %Define plot limits
    minx = 0;
    maxx = L_heated;
    miny = 0;
    maxy = max(max(q_LUT,q_profile));

    %Plot actual vs predicted heat flux as a function of position
    plot(xx,q_profile,'-k',xx,q_LUT,'-r',[L_LUT,L_LUT],[miny,CHF_LUT],'-k',[0 L_LUT],[CHF_LUT,CHF_LUT],'-k');
    legend('Heater Flux','LUT Flux');
    axis([minx maxx miny maxy]);
    xlabel('Heater Element Location [m]');
    ylabel('Heat Flux [kW/m^2]');
end
hold off;

if toggle_plotsave == 1
    FigSave = strcat(plotfilepath,'\',fig_title);
    print('-dpng','-r400',FigSave);
end
    
x_LUT_local = x_local(index_2);
K_Factors = K(:,index_2);
%Average heat flux [kW/m^2]
CHF_LUT_Avg = P_CHF_Avg(kk,2)/A_heated;

if iter == itmax
    fprintf('The maximum # of iterations was reached for 2006 LUT\n');
end
end

%****************************************%

function K_1 = K_1_Factor(D_hyd)
% K(1) - Diameter Effect Factor
% D_hyd => hydraulic diameter of subchannel         [m]

if D_hyd < 0.002
    K_1 = 1;
elseif D_hyd > 0.025
    K_1 = 0.57;
else
    K_1 = sqrt(0.008/D_hyd);
end
end

%****************************************%

function K_2 = K_2_Factor(D_htr,Pitch,x_quality)
% K(2) - Bundle Geometry Factor
% D_htr => heater element diameter                  [m]
% Pitch => center distance between heating elements [m]
% x_quality => same as under Primary Input          [-]

if x_quality < 0
    x_quality = 0;
end

delta = Pitch - D_htr;
arg(1) = 1;
arg(2) = (0.5 + 2*delta/D_htr)*exp(-x_quality^(1/3)/2);
K_2 = min(arg);
end

%****************************************%

function K_3 = K_3_Factor(D_hyd,G,K_g,L_sp)
% K(3) - Mid-Plane Spacer Factor
% D_hyd => same as under K(1)                       [m]
% G => same as under Primary Inputs                 [kg/m^2-s]
% K_g => pressure loss coefficient of spacer        [-]
% L_sp => distance between mid-plane of spacers     [m]

A = 1.5*sqrt(K_g)*(G/1000)^(0.2);
B = 0.10;
K_3 = 1+A*exp(-B*L_sp/D_hyd);
end

%****************************************%

function K_4 = K_4_Factor(D_hyd,L_htd,rho_fsat,rho_gsat,x_quality)
% K(4) - Heated Length Factor
% D_hyd => same as under K(1)                       [m]
% L_htd => heated length                            [m]
% rho_fsat => fluid saturation density              [kg/m^3]
% rho_gsat => vapor saturation density              [kg/m^3]
% x_quality => same as under Primary Input          [-]

if x_quality < 0
    alpha_h = 0;
else
    alpha_h = x_quality*rho_fsat/(x_quality*rho_fsat + (1-x_quality)*rho_gsat);
end

if L_htd/D_hyd >= 5
    K_4 = exp(D_hyd/L_htd*exp(2*alpha_h));
else
    K_4 = 1;
end
end

%****************************************%

function K_5 = K_5_Factor(q_BLA,q_local,x_quality)
% K(5) - Axial Flux Distribution Factor
% q_BLA => boling length average heat flux          [kW/m^2]
% q_local => heat flux at CHF location              [kW/m^2]
% x_quality => same as under Primary Inputs         [-]

if x_quality > 0
    K_5 = q_local/q_BLA;
else
    K_5 = 1;
end
end

%****************************************%

function K_6 = K_6_Factor(q_rc_avg,q_rc_max,x_quality)
% K(6) - Radial/Circumferential (R/C) Flux Distribution Factor
% q_rc_avg => average R/C flux at a height z        [kW/m^2]
% q_rc_max => maximum R/C flux at a height z        [kW/m^2]
% x_quality => same as under Primary Inputs         [-]

if x_quality > 0
    K_6 = q_rc_avg/q_rc_max;
else
    K_6 = 1;
end
end

%****************************************%

function K_7 = K_7_Factor(f_L,G,grav,rho_fsat,rho_gsat,x_quality)
% K(7) - Flow-Orientation Factor
% f_L => friction factor of the channel             [-]
% G => same as under Primary Inputs                 [kg/m^2-s]
% grav => gravity coefficient (e.g. 9.81 m/s^2)     [m/s^2]
% rho_fsat => same as under K(4)                    [kg/m^3]
% rho_gsat => same as under K(4)                    [kg/m^3]
% x_quality => same as under Primary Inputs         [-]

% This alpha should be changed to the correclation of Premoli et al. (1970)
if x_quality < 0
    alpha = 0;
else
    alpha = x_quality*rho_fsat/(x_quality*rho_fsat + (1-x_quality)*rho_gsat);
end

T_1 = ((1-x_quality)/(1-alpha))^2*f_L*G^2/...
      (grav*D_hyd*rho_fsat*(rho_fsat-rho_gsat)*sqrt(alpha));
K_7 = 1 - exp(-sqrt(T_1/3));
end

%****************************************%

function K_8 = K_8_Factor(G,rho_fsat,rho_gsat,x_quality)
% K(8) - Vertical Low-Flow Factor - Aminus sign refers to downward flow
% G => same as under Primary Inputs                 [kg/m^2-s]
% rho_fsat => same as under K(4)                    [kg/m^3]
% rho_gsat => same as under K(4)                    [kg/m^3]
% x_quality => same as under Primary Inputs         [-]

if x_quality < 0
    alpha = 0;
else
    alpha = x_quality*rho_fsat/(x_quality*rho_fsat + (1-x_quality)*rho_gsat);
end

if alpha_h < 0.8
    C_1 = 1.0;
else
    C_1 = (0.8 + 0.2*rho_fsat/rho_gsat)/...
          (alpha_h + (1-alpha_h)*rho_fsat/rho_gsat);
end

if G > -400 && G < 0
    K_8 = (1-alpha_h)*C_1;
else
    K_8 = 1;    
end
end

%****************************************%

function q_BLA = BoilingLengthAvg(A_heated,A_test,G,h_fg,L_htd,P_rod_avg,x_in,XLcurrent)
%       "Outputs"
% q_BLA => boiling length average heat flux     [kW/m^2]
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

XLratio = [L_zero(3)/L_htd XLcurrent];  %Integration limits for q_BLA
L_Boiling = L_htd*(XLratio(2)-XLratio(1)); %Boiling length
relative = theta_0*(XLratio(2)-XLratio(1)) + theta_1/(2*theta_2)*...
    (sin(2*theta_2*(XLratio(2)-0.5)) - sin(2*theta_2*(XLratio(1)-0.5)));

% Boiling length average heat flux
q_BLA = (P_rod_avg/A_heated)*relative/(XLratio(2)-XLratio(1));

end