function [CHF_LuWu,CHF_LuWu_Avg,L_LuWu,x_LuWu_local] = ...
    LuitjensWu(G,Inlet_Press,A_heated,A_test,D_hydw,...
               h_in,L_heated,toggle_vis,Factor_toggle,Identifier)
% Prediction of the critical heat flux using the Luitjens and Wu (OSU) as
% presented in 'Mechanistic CHF Modeling for Natural Circulation
% Applications in SMRs' from the 2015 US-Japan Thermal-Hydraulics Seminar,
% Purdue University, West Lafayette, Indiana, the paper 'Development of a
% Mechanistic Critical Heat Flux Correlation' from the NURETH-16, Chicago,
% Il, 2015. The 'Steps' are in reference to the NURETH-16 conference.
%
%       "Outputs"
% CHF_LuWu => Critical heat flux predicted by EPRI              [kW/m^2]
% CHF_LuWu_Avg => Average heat flux based on total heater power [kW/m^2]
% L_LuWu => Location of CHF event                               [m]
% x_LuWu_local => quality at location of predicted CHF          [-]
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

%% Step 1-2: 
%Specify specific channel geometry/Specify fluid boundary
%conditions (system pressure, inlet mass flux, inlet temperature/enthalpy)

%Inlet Parameters
% G = 585.6838;               %total axial mass velocity [kg/m^2-s]
% Inlet_Press = 12.555*10^6;  %Inlet Pressure [Pa]
% A_heated = 0.0601;          %heated area of heater lement [m^2]
% A_test = 0.000087311;       %total cross sectional flow area [m^2]
% D_hydw = 0.0117;            %wetted hydraulic diameter [m]
% h_in = 1.1067*10^6;         %inlet [J/kg]
% L_heated = 2.0076;          %heated lenth [m]
% toggle_vis = 'on';

%Pressure Dependent Thermodynamic Properties
Inlet_Press = Inlet_Press/(10^(6));             %[Pa] to [MPa]
T_sat = IAPWS_IF97('Tsat_p',Inlet_Press);       %saturation temperature [K]

rho_fsat = 1/IAPWS_IF97('vL_p',Inlet_Press);    %saturated liquid density [kg/m^3]
rho_gsat = 1/IAPWS_IF97('vV_p',Inlet_Press);    %vapor density (assumed at saturation) [kg/m^3]

h_fsat = IAPWS_IF97('hL_p',Inlet_Press).*1000;  %saturated liquid enthalpy [J/kg]
h_gsat = IAPWS_IF97('hV_p',Inlet_Press).*1000;  %saturated vapor enthalpy [J/kg]
h_fg = h_gsat-h_fsat;                           %latent heat of vaporization [J/kg]

P_crit = 22.064;                        %Critical Pressure [MPa]
T_crit = IAPWS_IF97('Tsat_p',P_crit);   %Critical Temperature [K]

%% Step 3: 
%Define resolution of axial nodalization

n_xx = 200;
xx = linspace(0,L_heated,n_xx); %axial position [m]
dx = L_heated/(n_xx-1);         %axial spacing [m]

%% Step 4: 
%Define the heat flux profile to test

theta_0 = 0.8187458177;
theta_1 = 0.6812541823;
theta_2 = 2.436354311;

PowProf = @(z) theta_0+theta_1*cos(2*theta_2*(z/L_heated -0.5)); %[-]

%Experimentally regressed fitting parameters
if Factor_toggle == 0
    a_1 = 26.55;  a_2 = -0.5631;  a_3 = -1.472e-3; %EPRI
elseif Factor_toggle == 1
    a_1 = 2.151; a_2 = -0.09458; a_3 = -5.594e-3; %UW
end

%Initialize iteration and set convergence criteria
Pow_htr = [1,2].*1000;          %initial average heater power guess [W]
itmax = 30;
iter = 0;
Error = 1;
tol = 0.005;

fig_title = sprintf('Luitjens and Wu Model Convergence Plot %s',Identifier);
figure_LuWuconv = figure('Name',fig_title,'NumberTitle','off','Visible',toggle_vis);
hold on;

%%
while (Error>tol && itmax>iter)
    q_pp_local = zeros(n_xx,2); q_CHF = zeros(n_xx,2);
    iter = iter + 1;
    
    for ii = 1:2
        
        lambda = zeros(n_xx,1); S = zeros(n_xx,1); Na = zeros(n_xx,1);
       
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
        
%% Step 5:
%Step 5a: Solve the one-dimensional energy equation

        %Average enthalpy as a function of position [J/kg]
        h_avg = (h_in + q_pp_net.*A_heated./(G.*A_test)); 
        x_local = (h_avg - h_fsat)./h_fg; %Local quality [-]
        
        h_Lavg = h_avg;
        h_Lavg(h_avg>h_fsat)=h_fsat; %pure liquid enthalpy [J/kg]
        
%Step 5b: Determine the each axial node fluid properties

        %Energy/Temp Dependent Thermodynamic Properties
        T_L = IAPWS_IF97('T1_ph',Inlet_Press,h_Lavg./1000); %local liquid temperature [K]
        mu_L = IAPWS_IF97('mu_pT',Inlet_Press,T_L-0.1);     %viscosity of liquid [kg/m-s]
        mu_fsat = IAPWS_IF97('mu_pT',Inlet_Press,T_sat-0.1);%viscosity of saturated liquid [kg/m-s]
        mu_gsat = IAPWS_IF97('mu_pT',Inlet_Press,T_sat);    %viscosity of saturated vapor [kg/m-s]
        sigma = IAPWS_IF97('sigma_T',T_L);                  %surface tension [N/m]
        C_pL = IAPWS_IF97('cp1_pT',Inlet_Press,T_L).*1000;  %specific heat of liquid [J/kg-K]
        k_L = IAPWS_IF97('k_pT',Inlet_Press,T_L-0.1);       %liquid thermal conductivity [W/m-K]
        rho_mix = 1./IAPWS_IF97('v_ph',Inlet_Press,h_avg./1000); %vapor/liquid mixture density [kg/m^3]
       
        Pr = C_pL.*mu_L./k_L;             %Prandtl number for liquid [-]
        Re = G.*D_hydw./mu_L;             %Reynolds number [-]
        Nu = 0.023.*Re.^(0.8).*Pr.^(0.4); %Nusselt Number [-]
        h_conv = Nu.*k_L./D_hydw;         %Single phase heat transfer coefficient [W/m^2-K]
        U_l = G./rho_mix;                 %Bulk Velocity [m/s]
        
        for kk = 1:n_xx
            if Re(kk) < 2320
                lambda(kk) = 64/Re(kk); %Laminar
            elseif Re(kk) > 1*10^5
                lambda(kk) = 0.0032 + 0.221*Re(kk)^(-0.237);
            else
                lambda(kk) = 0.3164/Re(kk)^(0.25);
            end
        end
                
%% Step 6, 10, & 11: 
%Determine the effective active cavity size from the calibarion curve 
%and the associated bubble dimensions
                
        %Minimum/Average active nucleation cavity size [m]
        r_c = (a_1*G^a_2 + a_3*10*Inlet_Press)/10^6;
        
        %Bubble Contact Angle [rad]
        theta_contact = pi/4;
        
        %Bubble Radius [m]
        r_b = r_c/sin(theta_contact);
                
        %Bubble Height [m]
        H_b = r_b*(1+cos(theta_contact));
        
        %Bubble Interfacial Surface Area [m^2]
        A_i = 2*pi*H_b*r_b;

%% Step 7:
%Determine the wall temperature using the Clausius-Claperyeon Eqn.

        %Internal Bubble Pressure [MPa]
        P_bubble = 2.*sigma./r_b./10^(6) + Inlet_Press;
        if P_bubble > P_crit
            P_bubble = P_crit;
        end
        
        %Wall Temperature [K]
%         T_w2 = IAPWS_IF97('Tsat_p',P_bubble);
        
        %Wall Temperature using Clausius-Claperyon [K]
        T_w = T_sat.*(1+2*sigma./(r_c.*rho_gsat.*h_fg));
        
%% Step 8: 
%Step 8a: Determine the boiling suppression factor (S)
        
        %Chen Boiling Model for Suppression Factor (S)
        %Martinelli Parameter for Turbulent-Turbulen Flow
        
        for kk = 1:n_xx
            
            if x_local(kk) <= 0
                Ff = 1.0;
            elseif x_local(kk) > 1
                Ff = 2.35*0.213^(0.736);
            else
                Xtt = ((1 - x_local(kk))./x_local(kk)).^(0.9)...
                    .*(rho_gsat./rho_fsat).^(0.5)...
                    .*(mu_fsat./mu_gsat).^(0.1);
                if 1/Xtt <= 0.1
                    Ff = 1.0;
                else
                    Ff = 2.35*(1/Xtt + 0.213)^(0.736);
                end
            end
            
            if x_local(kk) >= 0 && x_local(kk) <= 1
                x_Temp = x_local(kk);
            elseif x_local(kk) < 0
                x_Temp = 0.0;
            elseif x_local(kk) > 1.0
                x_Temp = 1.0;
            end
            
            Re_TP = G.*(1-x_Temp).*D_hydw./mu_L(kk).*Ff.^(1.25);
            S(kk) = 1./(1 + 1.4*10^(-5).*Re_TP);
        end

%Step 8b: Calculate the suppression corrected wall temperature

        %Effective Wall Superheat [K]        
        Delta_T_w = (T_w - T_sat).*S;
        
        %Effective Wall Temperature [K]
        T_w = Delta_T_w + T_sat;
        T_w(T_w>T_crit) = T_crit; %Limit maximum wall temperature
        
%% Step 9:
%Solve the near wall thermal fluid properties using the wall temperature
        
        P_satw = IAPWS_IF97('psat_T',T_w);          %Pressure at Wall [MPa]
        mu_Lw = IAPWS_IF97('mu_pT',P_bubble,T_w-0.1);%viscosity of liquid [kg/m-s]
        sigma_w = IAPWS_IF97('sigma_T',T_w);        %surface tension [N/m]
        rho_fsatw = 1./IAPWS_IF97('vL_p',P_satw);    %saturated liquid density [kg/m^3]
        rho_gsatw = 1./IAPWS_IF97('vV_p',P_satw);    %vapor density (assumed at saturation) [kg/m^3]
        h_fsatw = IAPWS_IF97('hL_p',P_satw).*1000;  %saturated liquid enthalpy [J/kg]
        h_gsatw = IAPWS_IF97('hV_p',P_satw).*1000;  %saturated vapor enthalpy [J/kg]
        h_fgw = h_gsatw-h_fsatw;                    %latent heat of vaporization [J/kg]
        
%% Step 12:
%Calculate the nucleation site density based on the wall superheat and the
%contact angle.
        
        %Nucleation Site Density [sites/m^2]
        for kk = 1:n_xx
            if Delta_T_w(kk) < 15
                Na(kk) = 3400.*(1 - cos(theta_contact)).*Delta_T_w(kk).^2;
            else
                Na(kk) = 0.34.*(1 - cos(theta_contact)).*Delta_T_w(kk).^(5.3);
            end
        end
        
%% Step 13:
%Calculate the surface boiling fraction
        
        %Average departure radius [m]
        ravg = (0.5)*2.54/100;
        
        %Fractional Boiling Area [-]. Limit to <= 1
        f_boil = Na.*pi.*ravg.^2;
        f_boil(f_boil>1) = 1;
  
%% Step 14:
%Step 14a: Individual terms of the force balance equation are solved
        
        %Shear Stress [kg/m-s^2]
        C_f = lambda./4; %Skin Friction Coefficient
        tau_w = C_f./2.*rho_fsatw.*U_l.^2;
        
        %Bubble Center Velocity [m/s]
        u_Lw = tau_w./mu_Lw.*H_b./2;
        
        %Bubble Reynolds Number [-]
        Re_b = rho_fsatw.*u_Lw.*r_b./mu_Lw;
        
        %Shear Lift Coefficient
        Gs = 1.0;
        C_L = 3.877.*Gs.^(0.5).*(Re_b.^(-2) + 0.014.*Gs.^2).^(0.25);
               
        %Surface Tension Force Component
        F_sigma = 2.*pi.*sigma_w.*r_c*sin(theta_contact); % O(-8)
        
        %Bubble Lift Force Component
        F_lift = 0.5*C_L.*rho_fsatw.*u_Lw.^2.*pi.*r_b^2;
        
        %Bubble Vapor Inertical Force Component
        F_i = A_i./(rho_gsatw.*h_fgw.^2);        
        
%Step 14b: Determine the maximum evaporative heat flux
        
        %Evaporative Interfacial Energy Transfer [W/m^2]
        q_evap = sqrt((F_sigma-F_lift)./F_i);
        for kk = 1:n_xx
            if imag(q_evap(kk))~=0 || isnan(q_evap(kk))
                q_evap(kk) = 0; % means no bubble can be generated
            end
        end
        
%% Step 15:
%Calcualte the single phase convective heat flux
        
        %Single Phase Convection Heat Flux [W/m^2]
        q_conv = h_conv.*(Delta_T_w);
        
%% Step 16 & 17:
%Calcuale the CHF for each node based on the sum of the boiling
%fraction corrected evaporative and convective heat fluxes
        
        %Luitjens and Wu Critical Heat Flux [W/m^2]
        q_CHF(:,ii) = q_conv.*(1-f_boil) + q_evap.*f_boil;
        
    end

%% Test for convergence and update the solution.

    error_1 = (q_CHF(:,1) - q_pp_local(:,1))./q_CHF(:,1);
    error_2 = (q_CHF(:,2) - q_pp_local(:,2))./q_CHF(:,2);
    [min_error_1,index_1] = min(abs(error_1));
    [min_error_2,index_2] = min(abs(error_2));

    %Calculate error and next step (Secant Method)
    gp = (min_error_2-min_error_1) / (Pow_htr(2)-Pow_htr(1));
    dPower = -min_error_2 / gp;
    
    if dPower > 8e4
        dPower = 20e3;
    end
    if dPower < 10e3
        Pow_htr_new = Pow_htr(2) + dPower/2;
    else
        Pow_htr_new = Pow_htr(2) + dPower/6;
    end

    % Update the error
    Error = abs(Pow_htr_new-Pow_htr(2))./Pow_htr(2);
    Pow_htr(1) = Pow_htr(2);
    Pow_htr(2) = Pow_htr_new;

    plot(xx,q_pp_local(:,2)./1000,'-k',xx,q_CHF(:,2)./1000,'-b')
    legend('Heater Flux','LuWu Flux');
    axis([0 L_heated 0 2000]);
    xlabel('Heater Element Location [m]');
    ylabel('Heat Flux [kW/m^2]');   
       
end
hold off;

%% Define the output variables based on converged solutions

if iter == itmax
    fprintf('The maximum # of iterations was reached for LuitjensWu\n');
end

CHF_LuWu = q_CHF(index_2,2)/1000;         %[kW/m^2]
CHF_LuWu_Avg = Pow_htr(2)/1000/A_heated;  %[kW/m^2]
L_LuWu = xx(index_2);                     %[m]
x_LuWu_local = x_local(index_2);          %[-]

end
