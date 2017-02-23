function [CHF_LM,CHF_LM_Avg,L_LM,x_LM_local] = LeeMudawwar()%...
%     LeeMudawwar(G,Inlet_Press,A_heated,A_test,D_hydw,...
%                h_in,L_heated,toggle_vis)
% Prediction of the critical heat flux using the Lee and Mudawwar model as
% presented in 'Mechanistic Critical Heat Flux Model for Subcooled Flow
% Boiling based on Local Bulk Flow Conditions' from the International
% Journal Multiphase Flow Vol. 14 No. 6 pp. 711-728, 1988.
%
%       "Outputs"
% CHF_LM => Critical heat flux predicted by EPRI              [kW/m^2]
% CHF_LM_Avg => Average heat flux based on total heater power [kW/m^2]
% L_LM => Location of CHF event                               [m]
% x_LM_local => quality at location of predicted CHF          [-]
%
%       "Inputs"
% G => test section mass flux                       [kg/m^2s]
% Inlet_Press => inlet Pressure                     [Pa]
% A_heated => total heated area per heater element  [m^2]
% A_test => flow area per subchannel                [m^2]
% D_hyd => hydraulic diameter                       [m]
% h_in => inlet subcooling                          [J/kg]
% L_heated => heated length per heater element      [m]
% toggle_vis => toggle visibility of the convergence plot. 'on'/'off'

%Initialize procedure with geomerty and inlet conditions.

Inlet_Press = 12e6;
G = 2000;
h_in = 1.427e6;
D_hydw = 0.008;  
A_test = D_hydw^2*pi/4;
A_heated = D_hydw*pi;

%Inlet Parameters
% G = 585.6838;               %total axial mass velocity [kg/m^2-s]
% G = 1511;
% Inlet_Press = 12.555*10^6;  %Inlet Pressure [Pa]
% Inlet_Press = 11.295e6;
% A_heated = 0.0601;          %heated area of heater lement [m^2]
% A_test = 0.000087311;       %total cross sectional flow area [m^2]
% D_hydw = 0.0117;            %wetted hydraulic diameter [m]
% h_in = 1.1067*10^6;         %inlet [J/kg]
% h_in = 1.067e6;
L_heated = 2.0076;          %heated lenth [m]
toggle_vis = 'on';

%Pressure Dependent Thermodynamic Properties
Inlet_Press = Inlet_Press/(10^(6));             %[Pa] to [MPa]
T_sat = IAPWS_IF97('Tsat_p',Inlet_Press);       %saturation temperature [K]

mu_fsat = IAPWS_IF97('mu_pT',Inlet_Press,T_sat-0.1);%viscosity of saturated liquid [kg/m-s]
mu_gsat = IAPWS_IF97('mu_pT',Inlet_Press,T_sat);    %viscosity of saturated vapor [kg/m-s]

rho_fsat = 1/IAPWS_IF97('vL_p',Inlet_Press);    %saturated liquid density [kg/m^3]
rho_gsat = 1/IAPWS_IF97('vV_p',Inlet_Press);    %vapor density (assumed at saturation) [kg/m^3]

h_fsat = IAPWS_IF97('hL_p',Inlet_Press).*1000;  %saturated liquid enthalpy [J/kg]
h_gsat = IAPWS_IF97('hV_p',Inlet_Press).*1000;  %saturated vapor enthalpy [J/kg]
h_fg = h_gsat-h_fsat;                           %latent heat of vaporization [J/kg]

P_crit = 22.064;                        %Critical Pressure [MPa]
T_crit = IAPWS_IF97('Tsat_p',P_crit);   %Critical Temperature [K]

grav = 9.81;    %Gravitational coefficient [m/s^2]

%Discritize the axial space and intialize guess values and power profile

n_xx = 200;
xx = linspace(0,L_heated,n_xx); %axial position [m]
dx = L_heated/(n_xx-1);         %axial spacing [m]

%Define the heat flux profile to test
theta_0 = 0.8187458177;
theta_1 = 0.6812541823;
theta_2 = 2.436354311;

PowProf = @(z) theta_0+theta_1*cos(2*theta_2*(z/L_heated -0.5)); %[-]

%Initialize iteration and set convergence criteria
Pow_htr = [90,91].*1000;          %initial average heater power guess [W]
itmax = 30;
iter = 0;
Error = 1;
tol = 0.01;

fig_title = 'Lee and Mudawwar Model Convergence Plot';
figure_LMconv = figure('Name',fig_title,'NumberTitle','off','Visible',toggle_vis);
hold on;

while (Error>tol && itmax>iter)
    q_pp_local = zeros(n_xx,2); q_CHF = zeros(n_xx,2);
    iter = iter + 1;
    
    for ii = 1:2
        %Calculate average and local properties based on 1-D energy equation
        
        %Average Heater Flux profile [W/m^2]
        q_pp_avg = Pow_htr(ii)/A_heated;
        
        q_pp_net = zeros(n_xx,1); qneti = zeros(n_xx-1,1);
        
        for ix = 2:n_xx
            %Net Heat Flux Input [W/m^2]
            q_pp_net(ix) = q_pp_avg*integral(PowProf,0,xx(ix))./2;
            
            %Heat Flux Associated with Each Interval [W/m^2]
            qneti(ix) =q_pp_avg*integral(PowProf,xx(ix-1),xx(ix));
        end
        
        %Local Heat Flux [W/m^2]
        for ix = 1:n_xx
            q_pp_local(ix,ii) = q_pp_avg*PowProf(xx(ix));
        end
        
        %Average enthalpy as a function of position [J/kg]
        h_avg = (h_in + q_pp_net.*A_heated./(G.*A_test));
        x_local = (h_avg - h_fsat)./h_fg; %Local quality [-]
        
        h_Lavg = h_avg;
        h_Lavg(h_avg>h_fsat)=h_fsat; %pure liquid enthalpy [J/kg]
        
        %Energy/Temp Dependent nodal/Local Thermodynamic Properties
        T_L = IAPWS_IF97('T1_ph',Inlet_Press,h_Lavg./1000); %liquid temperature [K]
        mu_L = IAPWS_IF97('mu_pT',Inlet_Press,T_L-0.1);     %viscosity of liquid [kg/m-s]
        sigma = IAPWS_IF97('sigma_T',T_L);                  %surface tension [N/m]
        C_pL = IAPWS_IF97('cp1_pT',Inlet_Press,T_L).*1000;  %specific heat of liquid [J/kg-K]
        k_L = IAPWS_IF97('k_pT',Inlet_Press,T_L-0.1);       %liquid thermal conductivity [W/m-K]
        
        %Calculate single phase convection coefficient for each node
        Pr = C_pL.*mu_L./k_L;             %Prandtl number for liquid [-]
        Re = G.*D_hydw./mu_L;             %Reynolds number [-]
        
        %Dittus-Boelter single phase heat transfer coefficient [W/m^2-K]
        %Equation
        Nu = 0.023.*Re.^(0.8).*Pr.^(0.4); %Nusselt Number [-]
        h_1phi = Nu.*k_L./D_hydw;
        
        %Calculate the bubble diameter
        %Equation 8b
        D_b = 0.00015.*sqrt(sigma./(grav.*(rho_fsat - rho_gsat))).*...
            (rho_fsat.*C_pL.*T_sat./(rho_gsat.*h_fg)).^(1.25);

         A1 = 0.35;  A2 = 240;   A3 = -0.8;
         
        %Equation 19
        h_sc = 230.*q_pp_local(:,ii).*sqrt(q_pp_local(:,ii)./(G.*h_fg)) ...
            ./(A1.*(T_sat-T_L).*(230.*sqrt(q_pp_local(:,ii)./(G.*h_fg))-1) ...
            +q_pp_local(:,ii)./h_1phi);
        
        %Guess initial microlayer thickness and iterate
%         T_m = T_sat - A1.*(T_sat-T_L);
%         
%         q_B = q_pp_local(:,ii) - h_sc.*(T_sat - T_m);
%         q_B(q_B<0) = 0;
%         delta_m_g = pi.*(0.0584).^2./2.*(rho_gsat./rho_fsat).^(0.2)....
%             .*(1+rho_gsat./rho_fsat).*sigma./rho_gsat...
%             .*(rho_gsat.*h_fg./q_B).^2;

        C_D = 2./3.*D_b./sqrt(sigma./(grav.*(rho_fsat-rho_gsat)));
        
        delta_m_g = 1*10^-6*ones(n_xx,1);
        delta_m = zeros(n_xx,1); L_m = D_b; 
        U_bL = zeros(n_xx,1); U_b = zeros(n_xx,1); G_m = zeros(n_xx,1);
        
        for kk = 1:n_xx
            
            itmax_m = 20;
            iter_m = 0;
            Error_m = 1;
            tol_m = 0.001;
            
            while(iter_m < itmax_m && Error_m > tol_m)
                iter_m = iter_m + 1;
                
                %Equation 26
                delta_m(kk) = 0.421*A2*Re(kk)^(A3-0.1)*G*rho_gsat*h_fg^2*D_b(kk)...
                    /(q_pp_local(kk,ii)-A1*h_sc(kk)*(T_sat-T_L(kk)))^2 ...
                    *(1+delta_m_g(kk)/(delta_m_g(kk)+D_b(kk))) ...
                    *sqrt(L_m(kk)*grav*(rho_fsat-rho_gsat)/(rho_fsat*C_D(kk)));
                
                %Liquid velocity at y = delta_m + D_b/2 [m/s]
                %Equation 11
                U_bL(kk) = 0.758*Re(kk)^(-0.1)*G/rho_fsat*(log(0.152*Re(kk)^(-0.1) ...
                    *G*(delta_m(kk)+D_b(kk)/2)/mu_L(kk))-0.61);
                
                %Vapor blanket velocity [m/s]
                %Equation 9
                U_b(kk) = sqrt(2*L_m(kk)*grav*(rho_fsat-rho_gsat)/(rho_fsat*C_D(kk)))+U_bL(kk);
                
                %Length of the vapor blanket [m]
                %Equation 15
                L_m(kk) = 2*pi*sigma(kk)*(rho_fsat + ...
                    rho_gsat)/(rho_fsat*rho_gsat*U_b(kk)^2);
            
                %Equation 7b
                C_D(kk) = 48*mu_L(kk)/(rho_fsat*D_b(kk)*(U_b(kk)-U_bL(kk)));
                               
                Error_m = abs(delta_m(kk) - delta_m_g(kk))/delta_m_g(kk);
                delta_m_g(kk) = delta_m(kk);
            end
            
            if iter_m == itmax_m
                fprintf('Micro layer thickness max iterations reached.\n');
            end
                        
            %Liquid mass flux flowing into the microlayer kg/m^2-s]
            %Equation 14
            G_m(kk) = rho_fsat*U_b(kk);
            
            %Predicted CHF [W/m^2]
            %Equation 3
            q_CHF(kk,ii) = G_m(kk)*delta_m(kk)*(h_fg + A1*C_pL(kk)*(T_sat-T_L(kk)))/L_m(kk);
                        
        end
    end
    
    
    % Test for convergence and update the solution.
    error_1 = (q_CHF(:,1) - q_pp_local(:,1))./q_CHF(:,1);
    error_2 = (q_CHF(:,2) - q_pp_local(:,2))./q_CHF(:,2);
    [min_error_1,index_1] = min(abs(error_1));
    [min_error_2,index_2] = min(abs(error_2));
    
    %Calculate error and next step (Secant Method)
    gp = (min_error_2-min_error_1) / (Pow_htr(2)-Pow_htr(1));
    dPower = -min_error_2 / gp;
    Pow_htr_new = Pow_htr(2) + dPower/2;
    
    % Update the error
    Error = abs(Pow_htr_new-Pow_htr(2));
    Pow_htr(1) = Pow_htr(2);
    Pow_htr(2) = Pow_htr_new;
    
    plot(xx,q_pp_local(:,2)./1000,'-k',xx,q_CHF(:,2)./1000,'-b')
    legend('Heater Flux','LLP Flux');
    axis([0 L_heated 0 6000]);
    xlabel('Heater Element Location [m]');
    ylabel('Heat Flux [kW/m^2]');
    
end
hold off;

% Define the output variables based on converged solutions

if iter == itmax
    fprintf('The maximum # of iterations was reached for LeeMudawwar.\n');
end

CHF_LM = q_CHF(index_2,2)/1000;         %[kW/m^2]
CHF_LM_Avg = Pow_htr(2)/1000/A_heated;  %[kW/m^2]
L_LM = xx(index_2);                     %[m]
x_LM_local = x_local(index_2);          %[-]

end


