function [CHF_Biasi,CHF_Biasi_Avg,L_Biasi,x_Biasi_local] = ...
    Biasi_Correlation(A_test,D_htr,D_hyd,G,h_fg,Inlet_Press,L_heated,...
    x_in,toggle_vis,toggle_plotsave,plotfilepath)
% Prediction of the critical heat flux using the Biasi Correlation as
% presented in 'A new correlation for round duct and uniform heating
% comparison with world data' European Atomic Energy Community 1967
%
%       "Outputs"
% CHF_Biasi => Critical heat flux predicted b the LUT             [kW/m^2]
% CHF_Biasi_Avg => Average heat flux based on total heater power  [kW/m^2]
% L_Biasi => Location of CHF event                                [m]
% x_Biasi_local => quality at location of predicted CHF           [-]
%
%       "Inputs"
% A_test => flow area per subchannel                [m^2]
% D_htr => heater element diameter                  [m]
% D_hyd => hydraulic diameter of subchannel         [m]
% G => mass flux per subchannel                     [kg/m^2s]
% h_fg => latent heat of vaporization               [J/kg]
% Inlet_Press => inlet pressure                     [Pa]
% L_heated => heated length per heater element      [m]
% x_in => inlet quality                             [-]

%Convert from SI to correlation specified units
A_test = A_test*100^2;              % [m^2] to [cm^2]
D_htr = D_htr*100;                  % [m] to [cm]
D_hyd = D_hyd*100;                  % [m] to [cm]  
G = G*1000/100^2;                   % [kg/m^2-s] to [g/cm^2-s]
h_fg = h_fg/1000;                   % [J/kg] to [J/g]
Inlet_Press = Inlet_Press/101325;   % [Pa] to [atm or ATA]
L_heated = L_heated*100;            % [m] to [cm]

% Define Correlation parameters
hp = -1.159 + 0.149*Inlet_Press*exp(-0.019*Inlet_Press) +...
    8.99*Inlet_Press/(10 + Inlet_Press^2);

yp = 0.7249 + 0.099*Inlet_Press*exp(-0.032*Inlet_Press);

if D_hyd >= 1
    alpha = 0.4;
else
    alpha = 0.6;
end

%Heated area of heater element [cm^2]
A_heated = pi*D_htr*L_heated;

%Constants in Stern Document 'Technical Specification SLTS-76'
theta_0 = 0.8187458177;
theta_1 = 0.6812541823;
theta_2 = 2.436354311;

%Define position vector
x_n = 200;
xx = linspace(0,L_heated,x_n); %[cm]

%Initialize arrays
x_local = zeros(x_n,1);

%Initial Guess
P_CHF_Avg(1,:) = [1,2].*1000;  %Heater power [W]

Xratio=xx./L_heated;                 

%Solution convergence criteria
Error = 1;
tol = 0.0001;
itmax = 100;
iter = 0;
kk = 0;

fig_title = 'Biasi Correlation Convergence Plot';
figure_Biasiconv = figure('Name',fig_title,'NumberTitle','off','Visible',toggle_vis);
hold on;
while (Error > tol && itmax > iter)
    
    q_flux_local = zeros(x_n,2);
    q_CHF = zeros(x_n,2);
    
    iter = iter + 1;
    kk = kk + 1;
    
    for j = 1:2
    for i=1:x_n

        %Average heat flux per heater [W/cm^2]
        q_Flux_Avg = P_CHF_Avg(kk,j)/A_heated;
        
        %Definate integral (relative) of the chopped cosine power profile
        %from Xratio = 0..Xratio
        Power_Profile = theta_0 + theta_1*cos(2*theta_2*(Xratio(i) - 0.5));
        relative = theta_0*Xratio(i)+...
            0.5*theta_1/theta_2*sin(2*theta_2*(Xratio(i) - 0.5))+...
            0.5*theta_1/theta_2*sin(theta_2);
        
        %Integrated heat flux up to location xx [W/cm^2]
        q_flux_integrated = q_Flux_Avg*relative;
        
        %Local heat flux at location xx [W/cm^2]
        q_flux_local(i,j) = q_Flux_Avg*Power_Profile;
                
        %Increase in quality due to heating up to location xx
        dx = q_flux_integrated*A_heated/(G*A_test*h_fg);
        
        %Local quality at location xx
        x_local(i) = x_in + dx;

        % Low quality region
        q_CHF_low_quality = 1883/(D_hyd^(alpha)*G^(1/6))*(yp/G^(1/6)-x_local(i));
        
        % High quality region        
        q_CHF_high_quality = 3780*hp/(D_hyd^(alpha)*G^(0.6))*(1-x_local(i));

        %Biasi Prediction is the maximum of the two quality regions
        q_CHF(i,j) = max(q_CHF_low_quality,q_CHF_high_quality);

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
            
    %[W/cm^2] to [kW/m^2],[cm] to [m] and store information at CHF location
    q_Biasi = q_CHF(:,2)./1000.*100.^2;
    q_profile = q_flux_local(:,2)./1000.*100.^2;
    CHF_Biasi = q_Biasi(index_2);
    L_Biasi = xx(index_2)/100;
        
    %Define plot limits
    minx = 0;
    maxx = L_heated/100;
    miny = 0;
    maxy = max(max(q_Biasi,q_profile));

    %Plot actual vs predicted heat flux as a function of position
    plot(xx./100,q_profile,'-k',xx./100,q_Biasi,'-r',[L_Biasi,L_Biasi],[miny,CHF_Biasi],'-k',[0 L_Biasi],[CHF_Biasi,CHF_Biasi],'-k');
    legend('Heater Flux','Biasi Flux');
    axis([minx maxx miny maxy]);
    xlabel('Heater Element Location [m]');
    ylabel('Heat Flux [kW/m^2]');
end
hold off;

if toggle_plotsave == 1
    FigSave = strcat(plotfilepath,'\',fig_title);
    print('-dpng','-r400',FigSave);
end

x_Biasi_local = x_local(index_2);

%Average heat flux [W/cm^2] to [kW/m^2]
CHF_Biasi_Avg = P_CHF_Avg(kk,2)/A_heated./1000*100.^2;
if iter == itmax
    fprintf('The maximum # of iterations was reached for Biasi\n');
end
end