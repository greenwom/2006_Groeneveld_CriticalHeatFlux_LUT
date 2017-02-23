function [CHF_W3,CHF_W3_Avg,L_W3,x_W3_local,F_C] = ...
    W3_Correlation(A_heated,A_test,G,h_fg,Inlet_Press,L_heated,x_in,...
    h_fsat,D_ee,D_h,cw_toggle,toggle_nonuni,toggle_vis,toggle_plotsave,plotfilepath)
% Prediction of the critical heat flux using the W-3 (Tong F-Factor)
% Correlation as presented in 'Prediction of Departure from Nucleate
% Boiling for an Axially Non-Uniform Heat Flux Distribution', Journal of
% Nuclear Energy 1967, pp. 241-248.
%
%
%       "Outputs"
% CHF_W3 => Critical heat flux predicted b the LUT             [kW/m^2]
% CHF_W3_Avg => Average heat flux based on total heater power  [kW/m^2]
% L_W3 => Location of CHF event                                [m]
% x_W3_local => quality at location of predicted CHF           [-]
% F_C => non-uniform correction factor at CHF location         [-]
%       "Inputs"
% A_heated => total heated area per heater element  [m^2]
% A_test => flow area per subchannel                [m^2]
% D_ee => subchannel wetted perimeter                [m]
% D_h => equivalent dimaeter based on heated perimeter [m]
% G => mass flux per subchannel                     [kg/m^2s]
% h_fg => latent heat of vaporization               [J/kg]
% h_fsat => saturate liquid enthalpy                [J/kg]
% Inlet_Press => inlet pressure                     [Pa]
% L_heated => heated length per heater element      [m]
% x_in => inlet quality                             [-]

A_heated = A_heated*(100/2.54)^2;   % [m^2] to [in^2]
A_test = A_test*(100/2.54)^2;       % [m^2] to [in^2]
D_ee = D_ee*100/2.54;               % [m] to [in]
D_h = D_h*100/2.54;                 % [m] to [in]
G = G*737.338;                      % [kg/m^2-s] to [lbm/ft^2-hr]
h_fg = h_fg/2326;                   % [J/kg] to [BTU/lbm]
h_fsat = h_fsat/2326;               % [J/kg] to [BTU/lbm]
h_in = x_in*h_fg + h_fsat;          % Inlet enthalphy [BTU/lbm]
Inlet_Press = Inlet_Press/6894.757; % [Pa] to [psia]
L_heated = L_heated*100/2.54;       % [m] to [in]

% Define Correlation parameters


if cw_toggle == 1
    D_e = D_h;
    Ra = 1-(D_ee/D_h);
else
    D_e = D_ee;
    Ra = 0;
end

%Constants in Stern Document 'Technical Specification SLTS-76'
theta_0 = 0.8187458177;
theta_1 = 0.6812541823;
theta_2 = 2.436354311;

%Define position vector
x_n = 200;
xx = linspace(0,L_heated,x_n); %[in]

%Initialize arrays
x_local = zeros(x_n,1);

%Initial Guess
P_CHF_Avg(1,:) = [30,31];             %Heater power [kW]
P_CHF_Avg = P_CHF_Avg.*3412.1416;   % [kW] to [BTU/hr]

Xratio=xx./L_heated;                 

%Solution convergence criteria
Error = 1;
tol = 0.0001;
itmax = 100;
iter = 0;
kk = 0;

fig_title = 'W-3 (Tong F-Factor) Correlation Convergence Plot';
figure_W3conv = figure('Name',fig_title,'NumberTitle','off','Visible',toggle_vis);
hold on;
while (Error > tol && itmax > iter)
    
    q_flux_local = zeros(x_n,2);
    q_CHF = zeros(x_n,2);
    
    iter = iter + 1;
    kk = kk + 1;
    
    for j = 1:2
    for i=1:x_n

        %Average heat flux per heater [BTU/in^2-hr]
        q_Flux_Avg = P_CHF_Avg(kk,j)/A_heated;
        
        %Definate integral (relative) of the chopped cosine power profile
        %from Xratio = 0..Xratio
        Power_Profile = theta_0 + theta_1*cos(2*theta_2*(Xratio(i) - 0.5));
        relative = theta_0*Xratio(i)+...
            0.5*theta_1/theta_2*sin(2*theta_2*(Xratio(i) - 0.5))+...
            0.5*theta_1/theta_2*sin(theta_2);
        
        %Integrated heat flux up to location xx [BTU/in^2-hr]
        q_flux_integrated = q_Flux_Avg*relative;
               
        %Local heat flux at location xx [BTU/in^2-hr]
        q_flux_local(i,j) = q_Flux_Avg*Power_Profile;
        
        %Increase in quality due to heating up to location xx
        dx = q_flux_integrated*A_heated/(G/12^2*A_test*h_fg);
        
        %Local quality at location xx
        x_local(i) = x_in + dx;
        
        %Pressure Factor (F_P)
        F_P(i) = ((2.022 - 0.0004302*Inlet_Press)+...
                  (0.1722 - 0.0000984*Inlet_Press)*...
                   exp((18.177 - 0.004129*Inlet_Press)*x_local(i)))*...
                 (1.157 - 0.869*x_local(i));
        
        %Flux Factor (F_G)
        F_G(i) = (0.1484 - 1.596*x_local(i) + 0.1729*x_local(i)*...
                  abs(x_local(i)))*G/10^6 + 1.037;
        
        %Diameter Factor (F_De)
        F_De(i) = 0.2664 + 0.8357*exp(-3.151*D_e);
        
        %Inlet Subcooling Factor (F_Hin)
        F_Hin(i) = 0.8258 + 0.000794*(h_fsat - h_in);

        %Cold wall Facot (cwf)
        cwf(i) = 1 - Ra*(13.76-1.372*exp(x_local(i)) ...
                  - 4.732*(G/10^6)^(-0.0535) ...
                  - 0.0619*(Inlet_Press/1000)^(0.14) ...
                  - 8.509*D_h^(0.107));
        
        %W-3 Correlation Uniform Heat Flux [BTU/in^2-hr]
        q_CHFu(i,j) = F_P(i)*F_G(i)*F_De(i)*F_Hin(i)*10^6/12^2;
        
        if toggle_nonuni == 1
            if xx(i) > 0
                L_DNB_non = xx(i);
                x_DNB = x_local(i);
                F_c(i) = NonUniFac(x_DNB,G,q_Flux_Avg,theta_0,theta_1,...
                    theta_2,L_DNB_non,L_heated,q_flux_local(i,j));
            else
                F_c(i) = 1;
            end
        else
            F_c(i) = 1;
        end
        
        q_CHF(i,j) = q_CHFu(i,j)*cwf(i)/F_c(i);
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
    
    % [BTU/in^2-hr] and [in] to [kW/m^2] and [m]
    q_W3 = q_CHF(:,2).*0.4542611;
    q_profile = q_flux_local(:,2).*0.4542611;
    CHF_W3 = q_W3(index_2);
    L_W3 = xx(index_2)*2.54/100;
        
    %Define plot limits
    minx = 0;
    maxx = L_heated*2.54/100;
    miny = 0;
    maxy = 2000;
    
    %Plot actual vs predicted heat flux as a function of position
    plot(xx*2.54/100,q_profile,'-k',xx*2.54/100,q_W3,'-r',[L_W3,L_W3],[miny,CHF_W3],'-k',[0 L_W3],[CHF_W3,CHF_W3],'-k');
    legend('Heater Flux','W3 Flux');
    axis([minx maxx miny maxy]);
    xlabel('Heater Element Location [m]');
    ylabel('Heat Flux [kW/m^2]');

end
hold off;

if toggle_plotsave == 1
    FigSave = strcat(plotfilepath,'\',fig_title);
    print('-dpng','-r400',FigSave);
end

%Quality at CHF location
x_W3_local = x_local(index_2);

%F_c factor at CHF location
F_C = F_c(index_2);

%Average heat flux [BTU/in^2-hr] to [kW/m^2]
CHF_W3_Avg = P_CHF_Avg(kk,2)/A_heated*0.4542611;

if iter == itmax
    fprintf('The maximum # of iterations was reached for W-3\n');
end
end

function [F_c] = NonUniFac(x_DNB,G,q_Flux_Avg,theta_0,theta_1,theta_2,...
                           L_DNB_non,L_heated,qfl)
%Correction Factor C
Cc = 0.44;
C = Cc*(1-x_DNB)^(7.9)/(G/10^6)^(1.72);

%Calculate Integral
qfloc = @(xloc) q_Flux_Avg.*(theta_0 + theta_1.*cos(2.*theta_2.*(xloc/L_heated - 0.5)));
integrand = @(xloc) qfloc(xloc).*exp(-C.*(L_DNB_non-xloc));
Int = integral(integrand,0,L_DNB_non);

%Non-Uniform Correction Factor
F_c = C/(qfl*(1-exp(-C*L_DNB_non)))*Int;
end