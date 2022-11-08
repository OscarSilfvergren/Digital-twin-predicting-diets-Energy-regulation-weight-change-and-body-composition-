load('Taylor1996_paramsCalibrated')

[sNames, ODEs, initialvalues] = IQMstates(model);
[pNames, param] = IQMparameters(model);

[row column] = size(Taylor1996_paramsCalibrated);
Taylor1996_paramsCalibrated = sortrows(Taylor1996_paramsCalibrated,column);

%% Cost to full estimation data

for j = 1:row
    
    TaylorTemp              = Taylor1996_paramsCalibrated(j,:);
    PersonSpecificParams(1)  = TaylorTemp(ismember(pNames,"carb_flow"));
    PersonSpecificParams(2)  = TaylorTemp(ismember(pNames,"protein_flow"));
    PersonSpecificParams(3)  = TaylorTemp(ismember(pNames,"lipids_flow"));
    PersonSpecificParams(4)  = TaylorTemp(column-4);
    PersonSpecificParams(5)  = TaylorTemp(column-3);
    PersonSpecificParams(6)  = TaylorTemp(column-2);
    PersonSpecificParams(7)  = TaylorTemp(column-1);
    costValidationTaylor = Taylor1996_costfunction(Taylor1996_data,9,13,13,14,log(PersonSpecificParams),TaylorTemp(1:AmountParameterOpti),model,initialvalues, pNames,sNames);
    
    Taylor1996_paramsCalibrated(j,column) = costValidationTaylor;
end

Taylor1996_paramsCalibrated = sortrows(Taylor1996_paramsCalibrated,column);

disp(' ')
if Taylor1996_paramsCalibrated(1,column) < chi2inv(0.95,517)
    disp('----- Fit to Estimation data is below statistical treshold -----')
end
fprintf('Best prediction of Rothman data: %.2f, Statistical Limit: %.2f (dgf = %i)', Taylor1996_paramsCalibrated(1,column), chi2inv(0.95,11), 11) % Finddatapoints
disp(' ')


%% Simulate Taylor1996


for j = 1:row
    
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5.6,initialvalues,sNames,Taylor1996_paramsCalibrated(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
    % InsulinResistancePersonalScaling
    optimizedParamTemp(ismember(pNames,"Vm_id"))                     = optimizedParamTemp(ismember(pNames,'Vm_id'))                   * Taylor1996_paramsCalibrated(j,end-1);
    optimizedParamTemp(ismember(pNames,"Vf_id"))                     = optimizedParamTemp(ismember(pNames,'Vf_id'))                   * Taylor1996_paramsCalibrated(j,end-1);
    optimizedParamTemp(ismember(pNames,"InsulinResponseFATP1_K"))   = optimizedParamTemp(ismember(pNames,'InsulinResponseFATP1_K')) * Taylor1996_paramsCalibrated(j,end-4);
    optimizedParamTemp(ismember(pNames,"InsulinResponseLiver"))     = optimizedParamTemp(ismember(pNames,'InsulinResponseLiver'))   * Taylor1996_paramsCalibrated(j,end-4);
    % InsulinClearancePersonalScaling
    optimizedParamTemp(ismember(pNames,"InsulinDegradation_PlasmaK"))    = optimizedParamTemp(ismember(pNames,'InsulinDegradation_PlasmaK'))    * Taylor1996_paramsCalibrated(j,end-3);
    optimizedParamTemp(ismember(pNames,"HE_positiveK"))                  = optimizedParamTemp(ismember(pNames,'HE_positiveK'))                  * Taylor1996_paramsCalibrated(j,end-3);
    % InsulinClearancePersonalScaling
    optimizedParamTemp(ismember(pNames,"InsulinStabilisationMagnitude"))    = optimizedParamTemp(ismember(pNames,'InsulinStabilisationMagnitude'))    * Taylor1996_paramsCalibrated(j,end-2);
    optimizedParamTemp(ismember(pNames,"InsulinProductionK"))               = optimizedParamTemp(ismember(pNames,'InsulinProductionK'))               * Taylor1996_paramsCalibrated(j,end-2);
    
    optimizedParamTempMeal               = Mealdeclaration(138.64,29.25,16.9377777777778,15,optimizedParamTemp,pNames);         %Fett är fel
    
    SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = 0;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = 0;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate      = model(0:480:480,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimTaylorToMeal      = model(1:1:360,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimTaylor_meal       = model(360:1:375,SimTaylorToMeal.statevalues(end,:),optimizedParamTempMeal);
    SimTaylorAfterMeal   = model(375:1:1440,SimTaylor_meal.statevalues(end,:),optimizedParamTemp);
    
    simTaylor1996.time           = [SimTaylorToMeal.time,           SimTaylor_meal.time(2:end),           SimTaylorAfterMeal.time(2:end)]-360;
    simTaylor1996.variablevalues = [SimTaylorToMeal.variablevalues; SimTaylor_meal.variablevalues(2:end,:); SimTaylorAfterMeal.variablevalues(2:end,:)];
    simTaylor1996.reactionvalues = [SimTaylorToMeal.reactionvalues; SimTaylor_meal.reactionvalues(2:end,:); SimTaylorAfterMeal.reactionvalues(2:end,:)];
    simTaylor1996.variables      = [SimTaylor_meal.variables];
    simTaylor1996.reactions      = [SimTaylor_meal.reactions];
    
    if j == 1
        simBest = simTaylor1996;
        
        maxGlycogenConcentrationliver1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'GlycogenConcentrationliver'));
        maxGlycogenConcentrationliver2 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver2 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'GlycogenConcentrationliver'));
        
        % GlucoseConcentrationPlasma
        maxGlucoseConcentrationPlasma1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'GlucoseConcentrationPlasma'));
        minGlucoseConcentrationPlasma1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'GlucoseConcentrationPlasma'));
        maxGlucoseConcentrationPlasma2 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'GlucoseConcentrationPlasma'));
        minGlucoseConcentrationPlasma2 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'GlucoseConcentrationPlasma'));
        
        % InsulinConcentrationPlasma
        maxInsulinConcentrationPlasma1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'InsulinConcentrationPlasma'));
        minInsulinConcentrationPlasma1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'InsulinConcentrationPlasma'));
        maxInsulinConcentrationPlasma2 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'InsulinConcentrationPlasma'));
        minInsulinConcentrationPlasma2 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'InsulinConcentrationPlasma'));
        
                % FFAConcentrationPlasma
        maxFFAConcentrationPlasma1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'FFAConcentrationPlasma'));
        minFFAConcentrationPlasma1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'FFAConcentrationPlasma'));
        maxFFAConcentrationPlasma2 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'FFAConcentrationPlasma'));
        minFFAConcentrationPlasma2 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'FFAConcentrationPlasma'));
        
        
    end
    
    %GlycogenConcentrationliver
    maxGlycogenConcentrationliver1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'GlycogenConcentrationliver'));
    minGlycogenConcentrationliver1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'GlycogenConcentrationliver'));
    maxGlycogenConcentrationliver2 = max(maxGlycogenConcentrationliver2,maxGlycogenConcentrationliver1);
    minGlycogenConcentrationliver2 = min(minGlycogenConcentrationliver2,minGlycogenConcentrationliver1);
    
    % GlucoseConcentrationPlasma
    maxGlucoseConcentrationPlasma1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'GlucoseConcentrationPlasma'));
    minGlucoseConcentrationPlasma1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'GlucoseConcentrationPlasma'));
    maxGlucoseConcentrationPlasma2 = max(maxGlucoseConcentrationPlasma2,maxGlucoseConcentrationPlasma1);
    minGlucoseConcentrationPlasma2 = min(minGlucoseConcentrationPlasma2,minGlucoseConcentrationPlasma1);
    
    % InsulinConcentrationPlasma
    maxInsulinConcentrationPlasma1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'InsulinConcentrationPlasma'));
    minInsulinConcentrationPlasma1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'InsulinConcentrationPlasma'));
    maxInsulinConcentrationPlasma2 = max(maxInsulinConcentrationPlasma2,maxInsulinConcentrationPlasma1);
    minInsulinConcentrationPlasma2 = min(minInsulinConcentrationPlasma2,minInsulinConcentrationPlasma1);
    
        % FFAConcentrationPlasma
    maxFFAConcentrationPlasma1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'FFAConcentrationPlasma'));
    minFFAConcentrationPlasma1 = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'FFAConcentrationPlasma'));
    maxFFAConcentrationPlasma2 = max(maxFFAConcentrationPlasma2,maxFFAConcentrationPlasma1);
    minFFAConcentrationPlasma2 = min(minFFAConcentrationPlasma2,minFFAConcentrationPlasma1);
    
    
end

timeFill = [simTaylor1996.time, fliplr(simTaylor1996.time)];

%% core variables Gly, I, G & EGP

figure('Name', "Taylor1996 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Glycogen
hold on
set(gca,'xtick',[0,5,10,15],'ytick',[0,100,200,300,400],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Hepatic Glycogen' ;'(mmol/L)'});
xlabel('Time (h)')
fill(timeFill/60,[maxGlycogenConcentrationliver2', fliplr(minGlycogenConcentrationliver2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBest.time/60, simBest.variablevalues(:,ismember(simBest.variables,'GlycogenConcentrationliver')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
errorbar(Taylor1996_data.time_glycogen_healthy(1:2)/60,Taylor1996_data.glycogen_healthy(1:2),Taylor1996_data.glycogenSEM_healthy(1:2),' k x','MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
errorbar(Taylor1996_data.time_glycogen_healthy(1:11)/60,Taylor1996_data.glycogen_healthy(1:11),Taylor1996_data.glycogenSEM_healthy(1:11),' k .','MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
xlim([-1 18])
ylim([0 400])
legend('Model uncertainty','Best prediction','Calibration Data','Validation data','Mixed Meal','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize+3)
hold off
grid on


%%
figure('Name', "Taylor1996 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Glucose
hold on
set(gca,'xtick',[0,5,10,15],'ytick',[0,3,6,9,12],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Plasma Glucose' ; '(mM)'});
xlabel('Time (h)')
fill(timeFill/60,[maxGlucoseConcentrationPlasma2', fliplr(minGlucoseConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot((simBest.time)/60, simBest.variablevalues(:,ismember(simBest.variables,'GlucoseConcentrationPlasma')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
errorbar(Taylor1996_data.time_glucose_healthy(1:2)/60,Taylor1996_data.glucose_healthy(1:2),Taylor1996_data.glucoseSEM_healthy(1:2),' k x','MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
errorbar(Taylor1996_data.time_glucose_healthy(1:end)/60,Taylor1996_data.glucose_healthy(1:end),Taylor1996_data.glucoseSEM_healthy(1:end),' k .','MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
xlim([-1 18])
ylim([0 12])
hold off
grid on

%%
figure('Name', "Taylor1996 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Insulin
hold on
set(gca,'xtick',[0,5,10,15],'ytick',[0,200,400,600,800,1000],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Plasma Insulin' ;'(pmol/L)'});
xlabel('Time (h)')
fill(timeFill/60,[maxInsulinConcentrationPlasma2', fliplr(minInsulinConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot((simBest.time)/60, simBest.variablevalues(:,ismember(simBest.variables,'InsulinConcentrationPlasma')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
errorbar(Taylor1996_data.time_Insulin_healthy(1:2)/60,Taylor1996_data.Insulin_healthy(1:2),Taylor1996_data.InsulinSEM_healthy(1:2),' k x','MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
errorbar(Taylor1996_data.time_Insulin_healthy(1:end)/60,Taylor1996_data.Insulin_healthy(1:end),Taylor1996_data.InsulinSEM_healthy(1:end),' k .','MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
xlim([-1 18])
ylim([0 1000])
hold off
grid on

%%
figure('Name', "Taylor1996 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% FFA
hold on
set(gca,'xtick',[0,5,10,15],'ytick',[0,0.3,0.6,0.9,1.2],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'FFA'; '(mM)'});
xlabel('Time (h)')
fill(timeFill/60,[maxFFAConcentrationPlasma2', fliplr(minFFAConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot((simBest.time)/60, simBest.variablevalues(:,ismember(simBest.variables,'FFAConcentrationPlasma')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
errorbar(Taylor1996_data.time_FFA_healthy(1:2)/60,Taylor1996_data.FFA_healthy(1:2),Taylor1996_data.FFASEM_healthy(1:2),' k x','MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
errorbar(Taylor1996_data.time_FFA_healthy(1:end)/60,Taylor1996_data.FFA_healthy(1:end),Taylor1996_data.FFASEM_healthy(1:end),' k .','MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
ylim([0 1.2])
xlim([-1 18])
hold off
grid on
