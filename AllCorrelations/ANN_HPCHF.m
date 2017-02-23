%%
% Solve an Input-Output Fitting problem with a Neural Network

%To Use Import G, Inlet Pressure, Inleat Subheating, CHF local, and
% CHF average from the Full Data Summary Sheet using the Import Data tool.
% Use the default titels assigned to the variables.

File='C:\Users\msg839\Desktop\HPCHF Data\Neural Network Look Up Table.xlsx';

%%
% === Get Data to Train ANN === %

Sheet = 'Exp Data';
% Range = 'B2:J71';
Range = 'B2:J97'; %This is the full range parameters
[Exp_Num,~] = xlsread(File,Sheet,Range);

Exp_Sub = Exp_Num(:,1);     %Inlet Pressure [kJ/kg]
Exp_Press = Exp_Num(:,2);   %Inlet Pressure [MPa]
Exp_G = Exp_Num(:,3);       %Inlet Mass Flux [kg/m2-s]
Exp_Qlocal = Exp_Num(:,8);  %Local CHF [kW/m2]
Exp_Qavg = Exp_Num(:,9);    %Average CHF [kW/m2]

inputs = [Exp_Sub Exp_Press Exp_G]';

%Switch between these for local vs average CHF
% targets = Exp_Qlocal';
targets = Exp_Qavg'; 

% +++ End +++ %

%%
% === Get Parameters for ANN Prediction === %

Sheet = 'LUT Matrix';
% Range = 'B2:D1174';
Range = 'B2:D6257'; %This is the full range parameters
[Num,~] = xlsread(File,Sheet,Range);

LUT_Sub = Num(:,1); %Inlet Pressure [kJ/kg]
LUT_Press = Num(:,2); %Inlet Pressure [MPa]
LUT_G = Num(:,3); %Inlet Mass Flux [kg/m2-s]

LUT_Inputs = [LUT_Sub LUT_Press LUT_G]';

% +++ End +++ %

%%
% Create a Fitting Network
hiddenLayerSize = 5;

num_iters = 1;
for j=1:num_iters
    net = fitnet(hiddenLayerSize);
    
    % Setup Division of Data for Training, Validation, Testing
    net.divideParam.trainRatio = 70/100;
    net.divideParam.valRatio = 15/100;
    net.divideParam.testRatio = 15/100;
    
    net.trainParam.showWindow = false;
    
    % Train the Network
    [net,tr] = train(net,inputs,targets);
    
    % Test the Network
    outputs = net(inputs);
    q_pred(:,j) = outputs';
    
    %Determine error from non-trained data (maximum error)
    nontr_ind = sort([tr.valInd, tr.testInd]);
        
    errors = (outputs(nontr_ind) - targets(nontr_ind))./targets(nontr_ind);
    RMS(j) = sqrt(sum(errors.^2)/length(errors));
    x_bar(j) = mean(errors);
    stdv(j) = std(errors);
 
    %%
    % === Error plots based on input parameters === %
%     figure;
%     subplot(2,2,1);
%     plot(inputs(1,nontr_ind),errors.*100,'ko');
%     xlabel('Inlet Subcooling [kJ/kg]');
%     ylabel('% Error');
%     axis([200 550 -50 50]);
%     
%     subplot(2,2,2);
%     plot(inputs(2,nontr_ind),errors.*100,'ko');
%     xlabel('Inlet Pressure [MPa]');
%     ylabel('% Error');
%     axis([6 18 -50 50]);
%     
%     subplot(2,2,3);
%     plot(inputs(3,nontr_ind),errors.*100,'ko');
%     xlabel('Inlet Mass Flux [kg/m^2]');
%     ylabel('% Error');
%     axis([300 1700 -50 50]);
    % +++ End +++ %
    
    %%
%     performance = perform(net,targets,outputs)

    % View the Network
%     view(net);

    %%
    % === Training Performance Plots === %
%     figure, plotperform(tr)
%     figure, plottrainstate(tr)
%     figure, plotregression(targets,outputs)
%     figure, ploterrhist(errors)
    % +++ End +++ %
    
    %%
    % === Extract Weights and Biases === %
%     weights=getwb(net);
%     [bb,iww,lww]=separatewb(net,weights);
    
%     b(:,j) = cell2mat(bb);
%     iw = cell2mat(iww);
%     lw = cell2mat(lww);
%     
%     iw_sub(:,j) = iw(:,1);
%     iw_press(:,j) = iw(:,2);
%     iw_flux(:,j) = iw(:,3);
%    
%     lw_output(:,j) = lw';
    % +++ End +++ %
    
    %%
    % === Use Trained Network For Prediction === %
    LUT_Outputs(:,j) = sim(net,LUT_Inputs)';
    % +++ End +++ %
    
    clearvars -except j q_pred RMS b iw_sub iw_press iw_flux lw_output...
                    inputs targets LUT_Inputs hiddenLayerSize x_bar stdv...
                    Inlet_SubheatJkg Inlet_PressurePa G_fluxkgm2s ...
                    q_CHFkWm2 q_CHF_Flux_AvgkWm2 LUT_Outputs
end
%%    
i=1;
q_pred_avg(:,i) = mean(q_pred,2);
RMS_avg(i) = mean(RMS)
Avg(i) = mean(x_bar)
StandDev(i) = mean(stdv)

LUT_Outputs_avg(:,i) = mean(LUT_Outputs,2);

% plot(RMS)
% b_avg = mean(b,2);
% iw_sub_avg = [mean(iw_sub,2);0];
% iw_press_avg = [mean(iw_press,2);0];
% iw_flux_avg = [mean(iw_flux,2);0];
% lw_output_avg = [mean(lw_output,2);0];
% wb_avg = [b_avg,iw_sub_avg,iw_press_avg,iw_flux_avg,lw_output_avg];
% dlmwrite('wb_avg.txt',wb_avg);

% dlmwrite('q_ANN_avg.txt',q_pred_avg);
% dlmwrite('q_ANN.txt',q_pred);
% dlmwrite('ErrorStats.txt',[Avg RMS_avg StandDev]);

dlmwrite('LUT_CHF_avg.txt',LUT_Outputs_avg);
% dlmwrite('LUT_CHF.txt',LUT_Outputs);

clearvars -except Inlet_SubheatJkg Inlet_PressurePa G_fluxkgm2s ....
                  q_CHFkWm2 q_CHF_Flux_AvgkWm2