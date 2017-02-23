function [CHF_WP,CHF_WP_Avg,L_WP,x_WP_local] = WeismanPei(G,Inlet_Press,...
    A_test,p_h,D_hyd,L_heated,h_in,A_heated,toggle_vis,Identifier)
% Prediction of the critical heat flux using the Weisman and Pei model as
% presented in 1983 'Prediction of Critical Heat Flux in Flow Boiling at
% Low Qualities' by Weisman and Pei and 1985 'A Theoretically based 
% Critical Heat Flux Prediction for Rod Bundles at PWR Condtions' by 
% Wesiman and Ying. The model also employs the use of the 1967 paper
% 'Forced Convection Subcooled Boiling - Prediction of Vapor Volumetric
% Fraction' by Levy.
%
%       "Outputs"
% CHF_WP => Critical heat flux predicted by EPRI              [kW/m^2]
% CHF_WP_Avg => Average heat flux based on total heater power [kW/m^2]
% L_WP => Location of CHF event                               [m]
% x_WP_local => quality at location of predicted CHF          [-]
%
%       "Inputs"
% G => mass flux per subchannel                     [kg/m^2s]
% Inlet_Press => inlet Pressure                     [Pa]
% A_heated => total heated area per heater element  [m^2]
% A_test => flow area per subchannel                [m^2]
% L_heated => heated length per heater element      [m]
% p_h => heated perimeter                           [m]
% D_hyd => hydraulic diameter                       [m]
% h_in => inlet subcooling                          [J/kg]
% toggle_vis => toggle visibility of the convergence plot. 'on'/'off'

theta_0 = 0.8187458177;
theta_1 = 0.6812541823;
theta_2 = 2.436354311;

PowProf = @(z) theta_0+theta_1*cos(2*theta_2*(z/L_heated -0.5)); %[-]

%Pressure Dependent Thermodynamic Properties
Inlet_Press = Inlet_Press/(10^(6));             %[Pa] to [MPa]
T_sat = IAPWS_IF97('Tsat_p',Inlet_Press);       %saturation temperature [K]

rho_fsat = 1/IAPWS_IF97('vL_p',Inlet_Press);    %saturated liquid density [kg/m^3]
rho_gsat = 1/IAPWS_IF97('vV_p',Inlet_Press);    %vapor density (assumed at saturation) [kg/m^3]
nu_fg = 1/rho_gsat - 1/rho_fsat;                %change in specific volume on vaporization [m^3/kg]

h_fsat = IAPWS_IF97('hL_p',Inlet_Press).*1000;  %saturated liquid enthalpy [J/kg]
h_gsat = IAPWS_IF97('hV_p',Inlet_Press).*1000;  %saturated vapor enthalpy [J/kg]
h_fg = h_gsat-h_fsat;                           %latent heat of vaporization [J/kg]

%*****************************

Pow_htr = [1,2].*1000;          %initial average heater power guess [W]

n_xx = 200;
xx = linspace(0,L_heated,n_xx); %axial position [m]
dx = L_heated/(n_xx-1);         %axial spacing [m]

itmax = 40;
iter = 0;
Error = 1;
tol = 0.001;

fig_title = sprintf('Weisman and Pei Model Convergence Plot %0s',Identifier);
figure_WPconv = figure('Name',fig_title,'NumberTitle','off','Visible',toggle_vis);
hold on;
while (Error>tol && itmax>iter)
    q_pp_local = zeros(n_xx,2); q_CHF = zeros(n_xx,2);
    iter = iter + 1;
    
    for ii = 1:2
       
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
        
        %Energy/Temp Dependent Thermodynamic Properties
        T_L = IAPWS_IF97('T1_ph',Inlet_Press,h_Lavg./1000); %local liquid temperature [K]
        mu_L = IAPWS_IF97('mu_pT',Inlet_Press,T_L-.1);      %viscosity of liquid [kg/m-s]
        sigma = IAPWS_IF97('sigma_T',T_L);                  %surface tension [N/m]
        C_pL = IAPWS_IF97('cp1_pT',Inlet_Press,T_L).*1000;  %specific heat of liquid [J/kg-K]
        k_L = IAPWS_IF97('k_pT',Inlet_Press,T_L-.1);        %liquid thermal conductivity [W/m-K]
        rho_L = 1./IAPWS_IF97('v1_pT',Inlet_Press,T_L);     %liquid density at bulk temperature [kg/m^3]
        
        Pr = C_pL.*mu_L./k_L;             %Prandtl number for liquid [-]
        Re = G.*D_hyd./mu_L;              %Reynolds number [-]
        Nu = 0.023.*Re.^(0.8).*Pr.^(0.4); %Nusselt Number [-]
        H_10 = Nu.*k_L./D_hyd;            %single phase heat transfer coefficient [W/m^2-K]
        
        %Calculate the Darcy Friction Factor [-]
        epsDh = 0.015*10^(-3)/D_hyd;      %Relative roughness [-]
        ffn = @(f,Re) 1/sqrt(f) + 2*log10(epsDh/3.7 + 2.51/(Re*sqrt(f)));
        f = zeros(n_xx,1);
        for ix = 1:n_xx
            f(ix) = fzero(@(f)ffn(f,Re(ix)),0.1);
        end
      
        %Non-dimensional distance from the heated wall to the tip of the
        %bubble from Levy (1967 - Forced Convective...)
        y_b_plus = 0.01.*sqrt(sigma.*D_hyd.*rho_fsat)./mu_L; % [-]
        
        %h_Ld - enthalpy of liquid at point of bubble departure [J/kg]
        h_Ld = zeros(n_xx,1);
        for ix = 1:n_xx
            if y_b_plus(ix) <= 5.0
                h_Ld(ix) = h_fsat - C_pL(ix)*q_pp_local(ix)/H_10(ix) + q_pp_local(ix)*Pr(ix)*y_b_plus(ix)/(G*sqrt(f(ix)/8));
            elseif y_b_plus(ix) <= 30.0 && y_b_plus(ix) > 5.0
                h_Ld(ix) = h_fsat - C_pL(ix)*q_pp_local(ix)/H_10(ix) + 5.0*q_pp_local(ix)/(G*sqrt(f(ix)/8))*...
                    (Pr(ix)+log(1+Pr(ix)*(y_b_plus(ix)/5-1)));
            elseif y_b_plus(ix) > 30
                h_Ld(ix) = h_fsat - C_pL(ix)*q_pp_local(ix)/H_10(ix) + 5.0*q_pp_local(ix)/(G*sqrt(f(ix)/8))*...
                    (Pr(ix)+log(1+5*Pr(ix))+0.5*log(y_b_plus(ix)/30));
            else
                fprintf('y_b_plus does not fit any conditions\n')
                return;
            end
        end

        %*****************************
        
        % Iterative section for x_avg and h_L
         x_avg  = zeros(n_xx,1); rho_avg   = zeros(n_xx,1);
         h_L    = zeros(n_xx,1); alpha_avg = zeros(n_xx,1);
         q_pp_b = zeros(n_xx,1); q_pp_cond = zeros(n_xx,1);
         
        for ix = 1:n_xx
                        
            if h_avg(ix) < h_fsat
                
                x_avg(ix) = 0;
                rho_avg(ix) = rho_L(ix);
                h_L(ix) = h_avg(ix);
                alpha_avg(ix) = 0;
                q_pp_b(ix) = 0;
                q_pp_cond(ix) = 0;
                
            else
                dstart = find(h_avg>=h_Ld); %Index of detachment point
                z_d = xx(dstart(1));        %Bubble detachment point [m]
                z_CHF = xx(ix);             %axial distance at which CHF occurs [m]
                
                x_avgg = x_local(ix); %guess - average quality [-]
                
                itxmax = 25;
                iterx = 0;
                xError = 1;
                tolx = 0.001;
                
                while (xError>tolx && itxmax>iterx)
                    
                    iterx = iterx + 1;
                    
                    fn1 = 0;
                    fn2 = 0;
                    
                    %Evaluate q_pp_b/q_pp_cond integrals from z_d to z_CHF
                    for id = dstart(1):ix
                        
                        %enthalpy of liquid [J/kg]
                        h_Lg = (h_avg(id)-h_gsat*x_avgg)/(1-x_avgg);             
                        
                        %average density [kg/m^3]
                        rho_avgg = rho_gsat/(rho_gsat/rho_L(id)*(1-x_avgg)+x_avgg); 
                        
                        %average void coefficient [-]
                        alpha_avgg = x_avgg*rho_avgg/rho_gsat;
                        
                        %?? h_fsat here may be h_fg but proly not [-]
                        epsilon = rho_L(id)*(h_fsat-h_Lg)/(rho_gsat*h_fg);  
                        
                        %q_pp_b - portion of total heat flux effective in
                        %         generating vapor [W/m^2]
                        if h_Lg <= h_Ld(id)
                            q_pp_bg = 0;
                        else
                            q_pp_bg = qneti(id)*(h_Lg-h_Ld(id))/(h_fsat-h_Ld(id));
                        end
                        
                        %q_pp_cond - heat flux due to bubble condensation [W/m^2]
                        if h_Lg <= h_Ld(id)
                            q_pp_condg = 0;
                        else
                            q_pp_condg = H_10(id)*(h_fg/nu_fg)*(A_test/p_h)*alpha_avgg*(T_sat-T_L(id));
                        end
                        
                        fn1 = fn1 + q_pp_bg/(1+epsilon)*dx;
                        fn2 = fn2 + q_pp_condg*dx;
                        
                    end
                    
                    %global quality [-]
                    x_avgg_new = p_h/(G*A_test*h_fg)*(fn1-fn2);
                    
                    xError = x_avgg_new - x_avgg;
                    
                    x_avgg = x_avgg_new;
                    
                end
                if iterx == itxmax
                    fprintf('Max iterations in quality was reached.\n');
                end

                x_avg(ix) = x_avgg;
                rho_avg(ix) = rho_avgg;
                h_L(ix) = h_Lg;
                alpha_avg(ix) = alpha_avgg;
                q_pp_b(ix) = q_pp_bg;
                q_pp_cond(ix) = q_pp_condg;
            end
        end
        
        %*****************************
        
        for ix = 1:n_xx
            
            %average viscosity [kg/m-s]
            mu_avg = mu_L(ix)*exp(2.5*alpha_avg(ix)/(1-39/64*alpha_avg(ix)));
            
            %Reynolds number based on average viscosity [-]
            Re_avg = G*D_hyd/mu_avg;
            
            %a - two-phase multiplier in turbulent intensity equation [-]
            if G <= 2694.44 %G in [kg/m^2-s]
                a = 0.135;
            else
                a = 0.135*(G((161+2/3)*10^(3)))^(-0.3);
            end
            
            %average bubble diameter [m]
            D_p = 0.015*sqrt(8*sigma(ix)*D_hyd*rho_avg(ix)/(f(ix)*G^2));
            
            %turbulent intensity [-] - 0.79 = 0.462(k)^0.6 where k is exp. regressed
            i_b = 0.79*(Re_avg)^(-0.1)*(D_p/D_hyd)^(0.6)*(1+a*(rho_L(ix)-rho_gsat)/rho_gsat);
            
            %standard deviation of nu_p [-]
            sigma_nu_p = i_b*G/rho_avg(ix);
            
            %radial velocity created by vapor generation [m/s]
            nu_1L = q_pp_b(ix)/(rho_gsat*h_fg);
            
            nu_sig = nu_1L/sigma_nu_p;
            psi(ix) = 1/sqrt(2*pi)*exp(-0.5*nu_sig^2)-0.5*nu_sig*erfc(0.5*nu_sig);
            
            r_0 = D_hyd/2;      %outer radius of tube [m]
            t = 5.5*D_p;        %bubble layer thickness [m]
            
            rho_1 = rho_avg(ix)*r_0^2/(r_0 - t)^2;      
            rho_2 = rho_avg(ix)*r_0^2/(2*t*(r_0-0.5*t));
            
            alpha_1 = (rho_L(ix)-rho_1)/(rho_L(ix)-rho_gsat);   %average void fraction in core [-]
            alpha_2 = 0.82;                                     %bubbly layer void fraction [-]
            
            x_1 = alpha_1*rho_gsat/rho_L(ix);   %average quality core region [-]
            x_2 = alpha_2*rho_gsat/rho_2;   %average quality in bubbly layer [-]
                        
            %Weisman & Pei model heat flux [W/m^2]
            q_CHF(ix,ii) = (x_2-x_1)*psi(ix)*i_b*h_fg*G*(h_fsat-h_Ld(ix))/(h_L(ix)-h_Ld(ix));

        end
    end
    
    ceiling = 2e3;
    for ix = 1:n_xx
        if q_CHF(ix,1) <= 0
            q_CHF(ix,1) = ceiling;
        end
        if q_CHF(ix,2) <= 0
            q_CHF(ix,2) = ceiling;
        end
    end
    
    if q_CHF(end,1) == ceiling
        Pow_htr(1) = Pow_htr(1) + 2000;
        Pow_htr(2) = Pow_htr(2) + 2000;
        Error = 1;
    else
        error_1 = (q_CHF(:,1) - q_pp_local(:,1))./q_CHF(:,1);
        error_2 = (q_CHF(:,2) - q_pp_local(:,2))./q_CHF(:,2);
        [min_error_1,index_1] = min(abs(error_1));
        [min_error_2,index_2] = min(abs(error_2));
        
        %Calculate error and next step (Secant Method)
        gp = (min_error_2-min_error_1) / (Pow_htr(2)-Pow_htr(1));
        dPower = -min_error_2 / gp;
        Pow_htr_new = Pow_htr(2) + dPower/4;
        
        % Update the error
        Error = abs(Pow_htr_new-Pow_htr(2))/Pow_htr(2);
        Pow_htr(1) = Pow_htr(2);
        Pow_htr(2) = Pow_htr_new;
        
    end
    
    plot(xx,q_pp_local(:,2)./1000,'-k',xx,q_CHF(:,2)./1000,'-b')
    legend('Heater Flux','WP Flux');
    axis([0 L_heated 0 2000]);
    xlabel('Heater Element Location [m]');
    ylabel('Heat Flux [kW/m^2]');      
    
end
if iter == itmax
    fprintf('Max iterations in WeismanPei.\n');
end
hold off;
CHF_WP = q_CHF(index_2,2)/1000;         %[kW/m^2]
CHF_WP_Avg = Pow_htr(2)/1000/A_heated;  %[kW/m^2]
L_WP = xx(index_2);                     %[m]
x_WP_local = (h_in + CHF_WP_Avg*1000*integral(PowProf,0,L_WP)/2*A_heated/(G*A_test) - h_fsat)/h_fg;              %[-]

end