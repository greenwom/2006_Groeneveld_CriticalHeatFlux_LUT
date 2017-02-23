function Predictive_Correlations()
% This file imports data from the Analyzed Data folder and uses results to
% evaluate the predicted CHF for various prediction methods detailed below.
% -CHF Prediction Methods
%       a) The 2006 Groeneveld CHF Look-Up Tables
%       b) 1967 Biasi
%       c) EPRI Correlation
%       d) CISE-GE Critical Quality
%       e) The W-3 or Tong 'F-Factor' Method
%       f) Biasi Critical Quality Correlation
%       g) Weisman & Pei Model
%       h) Luitjens & Wu Model
%
%Procedure:
%   1. Open raw Excel data and resave file (can save as xls/xlsx).
%
%   2. Make sure all .M, .txt, and .xlsm files are located in the in the current 
%       directory. If a given correlation is not being used the file is not
%       required, just ensure the the correlation toggle has been set to 0.
%       a) Groeneveld_LUT.m, 
%               -2006LUTdata.txt
%               -The 2006 CHF Groeneveld Look Up Table.xlsm
%       b) Biasi_Correlation.m
%       c) EPRI Correlation.m
%       d) CISEGE_Correlation.m
%       e) W3_Correlation.m
%       f) BiasiX_Correlation.m
%       g) WeismanPei.m
%       h) LuitjensWu.m
%       i) IAPWS_IF97 folder
%
%   3. Make sure the Groeneveld LUT Excel file is CLOSED and that the 
%       2006LUTdata.txt file is located in the current directory. The excel
%       file must be closed only if the 'ThreatLevel' variable (color 
%       coding of interpolated values) is to be returned.
%       -The program will NOT function if the Excel LUT is OPEN.
%       -Add/Remove correction factors (K)
%       -K factors 6-8 are not used currently.
%
%   4. 'Current Folder' can be any but all necessary files should be 
%       located in that same folder.
%
%   5. Toggle On/Off desired correlations
%
%   6. Update range of tests of desired condtions under variable 'Range'
%
%   7. Run Program.
%
%   8. Open COPYME_Correlations to view results.
%
%   9. Manually copy data from COPYME_Correlations to desired destination.
%
% TROUBLE SHOOTING/ERROR CHECK SECTION
%   1. ??? Error using ==> xlswrite => Reopen a separate copy of the
%   LUT with the exact same filename and then close it. Will not work if
%   their are any other Excel files currently open.%   

%=========================================================================%

IAPWS_Folder = strcat(cd,'\IAPWS_IF97');
addpath(IAPWS_Folder);

%Specify subchannel type for which to base calculations. 'channel' or 'rod'
toggle_chantype = 'channel';  
% toggle_chantype = 'rod';  

%Alter this section to turn on/off (1/0) desired correlations, plots, etc. 
toggle_vis = 'off'; %Convergence Plots
toggle_plotsave = 0;%Save plots to file

K_g = 1;            %Grid spacer pressure loss coefficient (1 is default)

Groeneveld_Toggle   = 0; %Groeneveld 2006 Look-Up Table
Biasi_Toggle        = 0; %1967 Biasi Correlation
EPRI_Toggle         = 0; %EPRI Correlation
CISEGE_Toggle       = 0; %CISE-GE Correlation
W3_Toggle           = 1; %W-3 or Tong F-Factor Correlation
BIASIX_Toggle       = 0; %Biasi Critical Quality Correlation
WP_Toggle           = 0; %Weisman and Pei Model
LuWu_Toggle         = 0; %Luitjens and Wu Model

%Groeneveld Correction Factors K(6-8) are not used
toggle_K = [0,0,0,0,0,0,0,0];
% Definition of correction factors K:
%   K(1) => Rod Diameter Factor
%   K(2) => Bundle Geometry Factor
%   K(3) => Mid-Plane Spacer Factor
%   K(4) => Heated-Length Factor
%   K(5) => Axial Flux Distribution Factor
%   K(6) => Radial/Circumferential Flux Distribution Factor
%   K(7) => Flow-Orientation Factor
%   K(8) => Vertical Low-Flow Factor

%This is color of the interpolated CHF cell. White is good, green is ok,
%red is don't use, blue is limited quality region.
toggle_ThreatLevel = 'on'; %Toggle ThreatLevel determination

%EPRI Parameters
cwall = 1;  %Toggle cold wall effects
nu = 1;     %Toggle nonuniform heat flux

%W-3 Nonuniform Toggle
toggle_nonuni = 1;
cw_toggle = 0;

%LuitjensWu EPRI or UW Factor Toggle. 0 = EPRI, 1 = UW
Factor_toggle = 1;
%=========================================================================%

%Perform a check of the toggle values for accidental user input error
Toggle_Values = [Groeneveld_Toggle,Biasi_Toggle,EPRI_Toggle,...
                 CISEGE_Toggle,W3_Toggle,cwall,nu,WP_Toggle,LuWu_Toggle];
Above_1_Check = any(Toggle_Values(:)>1);
Below_0_Check = any(Toggle_Values(:)<0);
if Above_1_Check == 1
    fprintf('Toggle value(s) are greater than 1.\n')
    fprintf('Please change toggle value(s) to either 0 or 1.\n');
    return;
elseif Below_0_Check == 1
    fprintf('Toggle value(s) are less than 0.\n');
    fprintf('Please change toggle value(s) to either 0 or 1.\n');
    return;
end

%=========================================================================%

%Import Required Parameters from Analyzed Data File
File='C:\Users\msg839\Desktop\HPCHF Data\HPCHF Data Analysis Summary.xlsx';
Sheet = 'Full Data Summary';
Range = 'A2:AK5';
% Range = 'A2:AK97'; %This is the full range of data
[Num,Txt] = xlsread(File,Sheet,Range);

Identifier_s = Txt(:,1);  %Test identifying name
Inlet_Press_s = Num(:,7); %Inlet Pressure [Pa]
Pr_s = Num(:,8);          %Reduced Pressure [-]
x_in_s = Num(:,9);        %Inlet thermodynamic quality [-]
G_flux_s = Num(:,11);     %Mass flux [kg/m^2-s]
h_fg_s = Num(:,18);       %Latent heat of Vaporization [J/kg]
h_fsat_s = Num(:,19);     %Saturated enthalpy [J/kg]
Delta_L_s = Num(:,20);    %Thermal Expansion Coefficient [m/m]
rho_fsat_s = Num(:,21);   %Saturated liquid density [kg/m^3]
rho_gsat_s = Num(:,22);   %Saturated vapor density [kg/m^3]

%Write selected data to Excel file
File_Save = strcat(pwd,'\COPYME_Correlations.xls');
fid = fopen(File_Save,'a');

%=========================================================================%

%%
%Gather data for a run, evaluate correlations, and write to excel file.
for test_num = 1:length(Identifier_s)
Identifier = char(Identifier_s(test_num));
Inlet_Press = Inlet_Press_s(test_num);
Pr = Pr_s(test_num);
x_in = x_in_s(test_num);
G_flux = G_flux_s(test_num);
h_fg = h_fg_s(test_num);
h_fsat = h_fsat_s(test_num);
Delta_L = Delta_L_s(test_num);
rho_fsat = rho_fsat_s(test_num);
rho_gsat = rho_gsat_s(test_num);
h_in = x_in*h_fg + h_fsat; %Inlet enthalpy [J/kg]

%Create Plot Directory if Plots are to be Saved
if toggle_plotsave == 1
    plotname = char(strcat('\Convergence Plots ',Identifier));
    plotfilepath = strcat(cd,plotname);
    if exist(plotfilepath,'dir') == 0
        mkdir(cd,plotname);
    end
else
    plotfilepath = '';
end
%=========================================================================%

%Fixed parameters
D_rod = 0.374;      %Rod Diameter [in]
Pitch = 0.496;      %Rod pitch [in]
Sq_width = 1.053;   %Width of square test section [in]
N_rods = 4;         %Number of rods
L_heated = 2;       %Heated length of rods [m]
pi = 4*atan(1.d0);  %The value of pi

%Convert to SI units
D_rod = D_rod*2.54/100;         %[m]
Pitch = Pitch*2.54/100;         %[m]
Sq_width = Sq_width*2.54/100;   %[m]

%Account for thermal expansion
% D_rod = D_rod*(1+Delta_L);          %Rod Diameter [m]
% Sq_width = Sq_width*(1+Delta_L);    %Width of square test section [m]
% L_heated = L_heated*(1+Delta_L);    %Heated length [m]

%Define key values based on defined variables.
p_h = pi*D_rod;                     %Heated perimeter [m]
A_heated = p_h*L_heated;            %Heated area per rod [m^2]
A_rod = D_rod^2/4*pi;               %Cross sectional area of one rod [m^2]
A_test = Sq_width^2 - N_rods*A_rod; %Total flow area of test section [m^2]

switch toggle_chantype
    case 'rod'
        A_subch = A_test/N_rods;    %Flow area of a single subchannel
        D_e = 2*Sq_width+pi*D_rod;  %Subchannel wetted perimeter
        D_sub = 4*A_subch/D_e;      %Subchannel hydraulic diameter
        D_h = 4*A_subch/(pi*D_rod); %equivalent diameter based on heated perimeter
    case 'channel'
        A_subch = Pitch^2 - A_rod;  %Flow area of a single subchannel
        D_e = pi*D_rod;             %Subchannel wetted perimeter
        D_sub = 4*A_subch/D_e;      %Subchannel hydraulic diameter
        D_h = 4*A_subch/(pi*D_rod); %equivalent diameter based on heated perimeter
    otherwise
        fprintf('Unrecognized subchannel type. Analysis aborted.\n');
        return;
end

%%
%=========================================================================%
%----------------------8. CHF CORRELATION SECTION-------------------------%
%=========================================================================%
%*Note - if more than one TC trips, these correlations, as applicable, will
%only calculate CHF based on the location, etc. of the first TC. This is
%appropriate as most CHF will occur near the same location and there is
%much uncertainty in the correlations already present. If CHF at additional
%TCs is required, run the files separately with the necessary
%information.
%
% The outputs are recorded on in the excel file COPYME on a new line
%
%    *Outputs*
% CHF_**** => Critical heat flux predicted by ****              [kW/m^2]
% CHF_****_Avg => Average heat flux based on total heater power [kW/m^2]
% L_**** => Location of CHF event                               [m]
% x_****_local => quality at location of predicted CHF          [-]
% other => addtional method specific variables may be returned 

%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/%

%Write the name of the test
fprintf(fid,'\n\r%s\t\n',Identifier);

%a) The 2006 Groeneveld CHF Look-Up Tables

if Groeneveld_Toggle == 1
% Prediction of the critical heat flux using the 2006 Groeneveld Look-Up
% Table as presented in 2006 CHF Look-Up Table, Nuclear Engineering and
% Design 237, pp. 190-1922.
       
    [CHF_LUT,CHF_LUT_Avg,L_LUT,x_LUT_local,K_Factors] = Groeneveld_LUT(...
        G_flux,Inlet_Press,x_in,A_heated,A_subch,L_heated,h_fg,...
        toggle_K,D_sub,D_rod,Pitch,K_g,rho_fsat,rho_gsat,toggle_vis,...
        toggle_plotsave,plotfilepath);
    
    ThreatLevel = 'Turned Off'; %Default label
    if strcmp(toggle_ThreatLevel,'on')
        
    % *** Obtain color coding of interpolated values ***
    %This is color of the CHF cell. White is good, green is ok,
    %red is don't use, blue is limited quality region.
    
    %String complete filename together
    filename = strcat(pwd,'\The 2006 CHF Groeneveld Look Up Table.xlsm');
    
    % Variables to be passed to spreadsheet
    Pressure = Inlet_Press/10^6; %Convert [Pa] to [MPa]
    Variables = [Pressure;G_flux;x_LUT_local];
    xlswrite('The 2006 CHF Groeneveld Look Up Table.xlsm',...
        Variables,'2006 LUT','P4');
    
    ExcelApplication = actxserver('Excel.Application'); %Starts ActXserver
    ExcelApplication.Visible = 0;                       %Toggle visibility
    ExcelApplication.Workbooks.Open(filename);          %Opens file
    AllSheets = ExcelApplication.ActiveWorkbook.Sheets; %Store Sheet Names
    
    LUT_Sheet = get(AllSheets,'Item','2006 LUT');
    LUT_Sheet.Activate;
    ExcelApplication.Run('CHF_Lin_Terp');               %Runs Macro

    Color = ExcelApplication.get('Range','P8');
    ThreatLevel = num2str(Color.value);
    ExcelApplication.DisplayAlerts = 0;     %Save file without prompt
    ExcelApplication.Quit;                  %Quits the spreadsheet
    ExcelApplication.release;               %Releases the spreadsheet
    end
    
    % Record results of Prediction on a new line    
    fprintf(fid,'\n');
    
    LUT_Labels = {'CHF_LUT','CHF_LUT_Avg','L_LUT','x_LUT_local',...
                  'K_1','K_2','K_3','K_4','K_5','K_6','K_7','K_8',...
                  'ThreatLevel'};
              
    for i =1:size(LUT_Labels,2)
        fprintf(fid,'%s\t',LUT_Labels{1,i});
    end
    
    fprintf(fid,'\n');
    
    LUT_Data = [CHF_LUT,CHF_LUT_Avg,L_LUT,x_LUT_local,K_Factors'];
    for i = 1:size(LUT_Data,2)
        fprintf(fid,'%f\t',LUT_Data(1,i));
    end
    fprintf(fid,'%s\t',ThreatLevel);   
end

%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/%

%b) 1967 Biasi Correlation

if Biasi_Toggle == 1
%'A new correlation for round duct and uniform heating
%comparison with world data' (1967)
%Valid from P = 2.7-140 bar and 10-6000 kg/m^2-s
    
    [CHF_Biasi,CHF_Biasi_Avg,L_Biasi,x_Biasi_local] = ...
        Biasi_Correlation(A_subch,D_rod,D_sub,G_flux,h_fg,Inlet_Press,...
        L_heated,x_in,toggle_vis,toggle_plotsave,plotfilepath);
    
    % Record results of Prediction on a new line
    fprintf(fid,'\n');
    
    Biasi_Labels = {'CHF_Biasi','CHF_Biasi_Avg','L_Biasi','x_Biasi_local'};
    
    for i =1:size(Biasi_Labels,2)
        fprintf(fid,'%s\t',Biasi_Labels{1,i});
    end
    
    fprintf(fid,'\n');
    
    Biasi_Data = [CHF_Biasi,CHF_Biasi_Avg,L_Biasi,x_Biasi_local];
    for i = 1:size(Biasi_Data,2)
        fprintf(fid,'%f\t',Biasi_Data(1,i));
    end
end

%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/%
%c) EPRI Correlation
if EPRI_Toggle == 1
% Prediction of the critical heat flux using the EPRI correlation as
% presented in 1982 EPRI Parametric Study of CHF Data Vol. 1-3
    
    [CHF_EPRI,CHF_EPRI_Avg,L_EPRI,x_EPRI_local] = ...
        EPRI_Correlation(G_flux,A_heated,A_subch,L_heated,x_in,...
             h_fg,Pr,K_g,cwall,nu,toggle_vis,toggle_plotsave,plotfilepath);
    
    % Record results of Prediction on a new line
    fprintf(fid,'\n');
    
    EPRI_Labels = {'CHF_EPRI','CHF_EPRI_Avg','L_EPRI','x_EPRI_local'};
    for i =1:size(EPRI_Labels,2)
        fprintf(fid,'%s\t',EPRI_Labels{1,i});
    end
    
    fprintf(fid,'\n');
    
    EPRI_Data = [CHF_EPRI,CHF_EPRI_Avg,L_EPRI,x_EPRI_local];
    for i = 1:size(EPRI_Data,2)
        fprintf(fid,'%f\t',EPRI_Data(1,i));
    end
end

%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/%

%d) CISE-GE Critical Quality Correlation
if CISEGE_Toggle == 1
% Prediction of the critical heat flux using the CISE-GE
% Correlation as presented in the TRACE vs. 5.0 USNRC Theory Manual.

    [CHF_CISEGE,CHF_CISEGE_Avg,L_CISEGE,x_CISEGE_local] = ...
        CISEGE_Correlation(A_heated,A_subch,G_flux,h_fg,Inlet_Press,...
                    L_heated,x_in,toggle_vis,toggle_plotsave,plotfilepath);
                       
    % Record results of Prediction on a new line
    fprintf(fid,'\n');
    
    CISEGE_Labels = {'CHF_CISEGE','CHF_CISEGE_Avg',...
                     'L_CISEGE','x_CISEGE_local'};
    for i =1:size(CISEGE_Labels,2)
        fprintf(fid,'%s\t',CISEGE_Labels{1,i});
    end
    
    fprintf(fid,'\n');
    
    CISEGE_Data = [CHF_CISEGE,CHF_CISEGE_Avg,...
                   L_CISEGE,x_CISEGE_local];
    for i = 1:size(CISEGE_Data,2)
        fprintf(fid,'%f\t',CISEGE_Data(1,i));
    end
end

%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/%

%e) The W-3 Correlation or Tong 'F-factor' Method
if W3_Toggle == 1
% Prediction of the critical heat flux using the W-3 (Tong F-Factor)
% Correlation as presented in 'Prediction of Departure from Nucleate
% Boiling for an Axially Non-Uniform Heat Flux Distribution', Journal of
% Nuclear Energy 1967, pp. 241-248.

    [CHF_W3,CHF_W3_Avg,L_W3,x_W3_local,F_c] = ...
    W3_Correlation(A_heated,A_subch,G_flux,h_fg,Inlet_Press,L_heated,...
                  x_in,h_fsat,D_e,D_h,cw_toggle,toggle_nonuni,toggle_vis,...
                  toggle_plotsave,plotfilepath);

    % Record results of Prediction on a new line
    fprintf(fid,'\n');
    
    W3_Labels = {'CHF_W3','CHF_W3_Avg',...
                     'L_W3','x_W3_local','F_c'};
    for i =1:size(W3_Labels,2)
        fprintf(fid,'%s\t',W3_Labels{1,i});
    end
    
    fprintf(fid,'\n');
    
    W3_Data = [CHF_W3,CHF_W3_Avg,...
                   L_W3,x_W3_local,F_c];
    for i = 1:size(W3_Data,2)
        fprintf(fid,'%f\t',W3_Data(1,i));
    end

end

%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/%

%f) BiasiX Critical Quality Correlation
if BIASIX_Toggle == 1
% Prediction of the critical heat flux using the Biasi Critical Quality (X)
% Correlation as presented in the TRACE vs. 5.0 USNRC Theory Manual.

    [CHF_BIASIX,CHF_BIASIX_Avg,L_BIASIX,x_BIASIX_local] = ...
        BiasiX_Correlation(A_heated,A_subch,G_flux,h_fg,Inlet_Press,...
                    L_heated,x_in,toggle_vis,toggle_plotsave,plotfilepath);
                       
    % Record results of Prediction on a new line
    fprintf(fid,'\n');
    
    BIASIX_Labels = {'CHF_BIASIX','CHF_BIASIX_Avg',...
                     'L_BIASIX','x_BIASIX_local'};
    for i =1:size(BIASIX_Labels,2)
        fprintf(fid,'%s\t',BIASIX_Labels{1,i});
    end
    
    fprintf(fid,'\n');
    
    BIASIX_Data = [CHF_BIASIX,CHF_BIASIX_Avg,...
                   L_BIASIX,x_BIASIX_local];
    for i = 1:size(BIASIX_Data,2)
        fprintf(fid,'%f\t',BIASIX_Data(1,i));
    end
end

%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/%

%g) The Weisman and Pei Model
if WP_Toggle == 1
% Prediction of the critical heat flux using the Weisman and Pei model as
% presented in 1983 'Prediction of Critical Heat Flux in Flow Boiling at
% Low Qualities' by Wesiman and Pei and 1985 'A Theoretically based 
% Critical Heat Flux Prediction for Rod Bundles at PWR Condtions' by 
% Wesiman and Ying. The model also employs the use of the 1967 paper
% 'Forced Convection Subcooled Boiling - Prediction of Vapor Volumetric
% Fraction' by Levy.

    [CHF_WP,CHF_WP_Avg,L_WP,x_WP_local] = WeismanPei(G_flux,Inlet_Press,...
        A_subch,p_h,D_sub,L_heated,h_in,A_heated,toggle_vis,Identifier);

    % Record results of Prediction on a new line
    fprintf(fid,'\n');
    
    WP_Labels = {'CHF_WP','CHF_WP_Avg',...
                     'L_WP','x_WP_local'};
    for i =1:size(WP_Labels,2)
        fprintf(fid,'%s\t',WP_Labels{1,i});
    end
    
    fprintf(fid,'\n');
    
    WP_Data = [CHF_WP,CHF_WP_Avg,...
                   L_WP,x_WP_local];
    for i = 1:size(WP_Data,2)
        fprintf(fid,'%f\t',WP_Data(1,i));
    end
end

%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/%

%h) The Luitjens and Wu Model
if LuWu_Toggle == 1
% Prediction of the critical heat flux using the Luitjens and Wu (OSU) as
% presented in 'Mechanistic CHF Modeling for Natural Circulation
% Applications in SMRs' from the 2015 US-Japan Thermal-Hydraulics Seminar,
% Purdue University, West Lafayette, Indiana.

    [CHF_LuWu,CHF_LuWu_Avg,L_LuWu,x_LuWu_local] = ...
        LuitjensWu(G_flux,Inlet_Press,A_heated,A_subch,D_sub,...
        h_in,L_heated,toggle_vis,Factor_toggle,Identifier);

    % Record results of Prediction on a new line
    fprintf(fid,'\n');
    
    LuWu_Labels = {'CHF_LuWu','CHF_LuWu_Avg',...
                     'L_LuWu','x_LuWu_local'};
    for i =1:size(LuWu_Labels,2)
        fprintf(fid,'%s\t',LuWu_Labels{1,i});
    end
    
    fprintf(fid,'\n');
    
    LuWu_Data = [CHF_LuWu,CHF_LuWu_Avg,...
                   L_LuWu,x_LuWu_local];
    for i = 1:size(LuWu_Data,2)
        fprintf(fid,'%f\t',LuWu_Data(1,i));
    end
end
fprintf(sprintf('%s prediction complete\n',Identifier))
end

fclose(fid);
end