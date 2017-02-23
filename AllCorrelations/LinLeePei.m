function [CHF_LLP,CHF_LLP_Avg,L_LLP,x_LLP_local] = LinLeePei()%...
%     LuitjensWu(G,Inlet_Press,A_heated,A_test,D_hydw,...
%                h_in,L_heated,toggle_vis)
% Prediction of the critical heat flux using the Lin, Lee, and Pei model as
% presented in 'An Improved Theoretical Critical HEat Flux Model for
% Low-Quality Flow' from Nuclear Technology Vol. 88, Dec. 1989.
%
%       "Outputs"
% CHF_LLP => Critical heat flux predicted by EPRI              [kW/m^2]
% CHF_LLP_Avg => Average heat flux based on total heater power [kW/m^2]
% L_LLP => Location of CHF event                               [m]
% x_LLP_local => quality at location of predicted CHF          [-]
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

%% Step 1: 
%Initialize procedure with geomerty and inlet conditions.

%Inlet Parameters
G = 585.6838;               %total axial mass velocity [kg/m^2-s]
% G = 1511;
Inlet_Press = 12.555*10^6;  %Inlet Pressure [Pa]
% Inlet_Press = 11.295e6;
A_heated = 0.0601;          %heated area of heater lement [m^2]
A_test = 0.000087311;       %total cross sectional flow area [m^2]
D_hydw = 0.0117;            %wetted hydraulic diameter [m]
h_in = 1.1067*10^6;         %inlet [J/kg]
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

%% Step 2:
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
Pow_htr = [45,42].*1000;          %initial average heater power guess [W]
itmax = 20;
iter = 0;
Error = 1;
tol = 0.01;

fig_title = 'Lin, Lee, and Pei Model Convergence Plot';
figure_LLPconv = figure('Name',fig_title,'NumberTitle','off','Visible',toggle_vis);
hold on;

%%
while (Error>tol && itmax>iter)
    q_pp_local = zeros(n_xx,2); q_CHF = zeros(n_xx,2);
    iter = iter + 1;
    
    for ii = 1:2
 %% Step 3:
%Calcualte average and local properties based on 1-D energy equation
      
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
        X_e = x_local; %equilibrium local quality
        
        h_Lavg = h_avg;
        h_Lavg(h_avg>h_fsat)=h_fsat; %pure liquid enthalpy [J/kg]
    
        %Energy/Temp Dependent nodal/Local Thermodynamic Properties
        T_L = IAPWS_IF97('T1_ph',Inlet_Press,h_Lavg./1000); %liquid temperature [K]
        mu_L = IAPWS_IF97('mu_pT',Inlet_Press,T_L-0.1);     %viscosity of liquid [kg/m-s]
        sigma = IAPWS_IF97('sigma_T',T_L);                  %surface tension [N/m]
        C_pL = IAPWS_IF97('cp1_pT',Inlet_Press,T_L).*1000;  %specific heat of liquid [J/kg-K]
        k_L = IAPWS_IF97('k_pT',Inlet_Press,T_L-0.1);       %liquid thermal conductivity [W/m-K]
        rho_mix = 1./IAPWS_IF97('v_ph',Inlet_Press,h_avg./1000); %vapor/liquid mixture density [kg/m^3]

%% Step 4:
%Calculate single phase convection coefficient for each node
        Pr = C_pL.*mu_L./k_L;             %Prandtl number for liquid [-]
        Re = G.*D_hydw./mu_L;             %Reynolds number [-]
        
        %Dittus-Boelter single phase heat transfer coefficient [W/m^2-K]
        %Equation 32
        Nu = 0.023.*Re.^(0.8).*Pr.^(0.4); %Nusselt Number [-]
        h_1phi = Nu.*k_L./D_hydw;         

%Calculate the bulk equilibrium quality at the location of onset of
%net vapor generation (ONVG - start of bulk saturated conditions - x=0)
      
        d_index = find(x_local >= 0);
        if isempty(d_index)
            q_pp_d = 0;
            Delta_T_d = 0;
            X_d = 0;
        else
            %Index of ONVG
            d = d_index(1);
            
            %Single-phase wall shear stress [Pa]
            %Equation 20
            tau_w0 = 0.5*(0.046*Re(d)^(-0.2))*G^2/rho_fsat;
        
            %Friction velocity for single-phase flow [m/s]
            %Equation 19
            U_t0 = sqrt(tau_w0/rho_fsat);
        
            %Distance from the wall [m]
            %Equation 35
            Y_b = 0.015*sqrt(sigma(d)*D_hydw/tau_w0);
        
            %Dimensionless distance from the wall [-]
            %Equation 34
            Y_bplus = Y_b*rho_fsat*U_t0/mu_fsat;
        
            %Dimensionless bulk temperature at the location of ONVG [-]
            %Equation 33a-c
            if Y_bplus < 5
                T_bplus = Pr(d)*Y_bplus;
            elseif Y_bplus >= 5 && Y_bplus < 30
                T_bplus = 5*(Pr(d) + ...
                                log(1 + Pr(d)*(Y_bplus/5 - 1)));
            else
                T_bplus = 5*(Pr(d) + log(1 + 5*Pr(d)) + ...
                                0.5*log(Y_bplus/30));
            end
       
            %Local heat flux at location of ONVG [W/m^2]
            q_pp_d = q_pp_local(d);
            
            %Subcooling at location of ONVG [K]
            %Equation 31
            Delta_T_d = q_pp_d*(1/h_1phi(d) - ...
                                T_bplus/(C_pL(d)*rho_fsat.*U_t0));
                            
            %Bulk equilibrium quality at ONVG [-]
            %Equation 30
            X_d = -C_pL(d).*Delta_T_d/h_fg;
        end
     
%% Step 5:
%Calculate the bubble diameter
        D_b = 0.00015.*sqrt(sigma./(grav.*(rho_fsat - rho_gsat))).*...
                       (rho_fsat.*C_pL.*T_sat./(rho_gsat.*h_fg)).^(1.25);
                   
%% Step 6:           
        %Bulk true quality [-]
        %Equation 29
        X_t = zeros(n_xx,1);
        for kk = 1:n_xx
            if X_e(kk) < X_d
                X_t(kk) = X_e(kk);
            else
                X_t(kk) = X_e(kk) - X_d*exp(X_e(kk)/X_d - 1);
            end
        end
        %Effective quality to simplify two-phase property coding
        X_teff = X_t;
        X_teff(X_t<0) = 0; 
        
        %Void fraction based on homogeneous two-phase flow model [-]
        %Equation 36
        alpha = X_teff./(X_teff + (1 - X_teff).*rho_gsat./rho_fsat);
                
        %Homogeneous fluid density [kg/m^3]
        %Equation 27
        rho_2phi = rho_fsat.*(1 - alpha) + rho_gsat.*alpha;

        %Mean two-phase viscosity [kg/m-s]
        %Equation 28
        mu_2phi = rho_2phi.*(X_teff.*mu_gsat./rho_gsat + ...
                             (1-X_teff).*mu_fsat./rho_fsat);
                         
        %Two-phase Reynolds number [-]
        %Between Equation 26 and 27
        Re_2phi = G.*D_hydw./mu_2phi;
        
        %Wall shear stress [Pa]
        %Equation 24
        tau_w = 0.5.*(0.046.*Re_2phi.^(-0.2)).*G.^2./rho_fsat;
        
        %True flow velocity [m/s]
        %Equation 23
        U_t = sqrt(tau_w./rho_2phi);
        
        %Liquid velocity perpendicular to heater surface (y-dir) [m/s]
        %Equation 25
        n_yy = 100;
        y = linspace(1e-5,0.00155,n_yy);
        U_L = zeros(n_xx,n_yy);
        for kk = 1:n_xx
            U_L(kk,:) = U_t(kk).*(5.*...
                    log(rho_2phi(kk).*y.*U_t(kk)./mu_2phi(kk)) - 3.05);
        end
                
%% Step 7:
%Guess initial microlayer thickness and iterate
        delta_m = 0.000001*ones(n_xx,1); %[m] 
        t_index = find(X_t>0);
        delta_m(1:t_index(1)-1) = 0;
        alpha_index = find(alpha>=0.7); %Model requires alpha < 0.7
        alpha_index(isempty(alpha_index))=n_xx;
        q_CHF = 1e7*ones(n_xx,2);
        node_error = zeros(n_xx,1);
        
        itmax_m = 20;
        iter_m = 0;
        Error_m = 1;
        tol_m = 0.001;
        
        while(iter_m < itmax_m && Error_m > tol_m)
            iter_m = iter_m + 1;
            
            for kk = t_index(1):alpha_index(1)
%% Step 8:
            a_1 = 7000; a_2 = -0.4; a_3 = -0.8; a_4 = 0.5;
            %Equation 47
            Y_r = (delta_m(kk) + D_b(kk))/delta_m(kk);
            %Equation 46
            C = a_1*Y_r^a_2*Re(kk)^(a_3 - a_4*alpha(kk)/(1-X_teff(kk)));

            %Effective velocity gradient
            %Equation 44
            dUL_dy = 2.5*U_t(kk)*(1/delta_m(kk) + ...
                            1/(D_b(kk) + delta_m(kk)));

            %Liquid velocity at y = delta_m + D_b/2 [m/s]
            %Equation 26
            U_bL = U_t(kk)*(5*log(rho_2phi(kk)*(delta_m(kk) + ...
                            D_b(kk)/2).*U_t(kk)/mu_2phi(kk)) - 3.05);
U_bL2(kk) = U_bL;
            %Vapor blanket velocity [m/s]
            %Equations 13, 14, 15 and 16
            S_3 = pi*sigma(kk)*grav*D_b(kk)*(rho_fsat^2 - ...
                            rho_gsat^2)/(12*rho_fsat*rho_gsat*mu_fsat);
            S_2 = (S_3/2 + U_bL^3/27 - ...
                    sqrt(S_3^2/4 + S_3*U_bL^3/27))^(1/3);
            S_1 = (S_3/2 + U_bL^3/27 + ...
                    sqrt(S_3^2/4 + S_3*U_bL^3/27))^(1/3);
            U_b = S_1 + S_2 + U_bL/S_3;
        
%% Step 9:
            %Length of the vapor blanket [m]
            %Equation 6
            L_m = 2*pi*sigma(kk)*(rho_fsat + ...
                        rho_gsat)/(rho_fsat*rho_gsat*U_b^2);

            %Liquid mass flux flowing into the microlayer kg/m^2-s]
            %Equation 7
            G_m = rho_fsat*U_b;
        
%% Step 10:  
            %Predicted CHF [W/m^2]
            %Equation 4
            q_CHF(kk,ii) = G_m*delta_m(kk)*h_fg/L_m;

            %Vapor momentum transfer rate into the vapor blanket
            %Equation 41
            F_I = q_CHF(kk,ii)^2*D_b(kk)*L_m/(rho_gsat*h_fg^2);
        
%% Step 11:
            %Equation 52
            Delta_U = F_I/(C*rho_fsat*dUL_dy*pi/4*D_b(kk)^2*L_m);
        
%% Step 12:
            %Equation 53
            U_bL = U_b - Delta_U;
U_bL3(kk) = U_bL;            
%% Step 13:
%Update micro layer thickness
            %Equation 54
            delta_m_new = mu_2phi(kk)/(rho_2phi(kk)*U_t(kk))*...
                          exp((U_bL/U_t(kk) + 3.05)/5) - 0.5*D_b(kk);
                      
%% Step 14:
%Convergence is based on micro layer thickness error
            node_error(kk) = abs(delta_m(kk) - delta_m_new)/delta_m(kk);
            delta_m(kk) = delta_m_new;
            end
            Error_m = max(node_error);
        end
       
        if iter_m == itmax_m
            fprintf('Micro layer thickness max iterations reached.\n');
        end
    end
    
%% Test for convergence and update the solution.
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
    axis([0 L_heated 0 2000]);
    xlabel('Heater Element Location [m]');
    ylabel('Heat Flux [kW/m^2]');   
       
end
hold off;

%% Define the output variables based on converged solutions

if iter == itmax
    fprintf('The maximum # of iterations was reached for LinLeePei.\n');
end

CHF_LLP = q_CHF(index_2,2)/1000;         %[kW/m^2]
CHF_LLP_Avg = Pow_htr(2)/1000/A_heated;  %[kW/m^2]
L_LLP = xx(index_2);                     %[m]
x_LLP_local = x_local(index_2);          %[-]

end
