load('Firth1986_paramsCalibrated')

[sNames, ODEs, initialvalues] = IQMstates(model);
[pNames, param] = IQMparameters(model);

[row column] = size(Firth1986_paramsCalibrated);
Firth1986_paramsCalibrated = sortrows(Firth1986_paramsCalibrated,column);

%% Cost to full estimation data

for j = 1:row
    
    FirthTemp              = Firth1986_paramsCalibrated(j,:);
    PersonSpecificParams(1)  = FirthTemp(ismember(pNames,"carb_flow"));
    PersonSpecificParams(2)  = FirthTemp(ismember(pNames,"protein_flow"));
    PersonSpecificParams(3)  = FirthTemp(ismember(pNames,"lipids_flow"));
    PersonSpecificParams(4)  = FirthTemp(column-4);
    PersonSpecificParams(5)  = FirthTemp(column-3);
    PersonSpecificParams(6)  = FirthTemp(column-2);
    PersonSpecificParams(7)  = FirthTemp(column-1);
    costValidationFirth = Firth1986_costfunction(Firth1986_data,14,15,15,log(PersonSpecificParams),FirthTemp(1:AmountParameterOpti),model,initialvalues, pNames,sNames);
    Firth1986_paramsCalibrated(j,column) = costValidationFirth;
end

Firth1986_paramsCalibrated = sortrows(Firth1986_paramsCalibrated,column);

disp(' ')
if Firth1986_paramsCalibrated(1,column) < chi2inv(0.95,44)
    disp('----- Fit to Estimation data is below statistical treshold -----')
end
fprintf('Best prediction of Rothman data: %.2f, Statistical Limit: %.2f (dgf = %i)', Firth1986_paramsCalibrated(1,column), chi2inv(0.95,44), 44) % Finddatapoints
disp(' ')


%% Simulate Firth1986


for j = 1:row
    
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5.6,initialvalues,sNames,Firth1986_paramsCalibrated(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
    % InsulinResistancePersonalScaling
    optimizedParamTemp(ismember(pNames,"Vm_id"))                     = optimizedParamTemp(ismember(pNames,'Vm_id'))                   * Firth1986_paramsCalibrated(j,column-1);
    optimizedParamTemp(ismember(pNames,"Vf_id"))                     = optimizedParamTemp(ismember(pNames,'Vf_id'))                   * Firth1986_paramsCalibrated(j,column-1);
    optimizedParamTemp(ismember(pNames,"InsulinResponseFATP1_K"))   = optimizedParamTemp(ismember(pNames,'InsulinResponseFATP1_K')) * Firth1986_paramsCalibrated(j,column-4);
    optimizedParamTemp(ismember(pNames,"InsulinResponseLiver"))     = optimizedParamTemp(ismember(pNames,'InsulinResponseLiver'))   * Firth1986_paramsCalibrated(j,column-4);
    % InsulinClearancePersonalScaling
    optimizedParamTemp(ismember(pNames,"InsulinDegradation_PlasmaK"))    = optimizedParamTemp(ismember(pNames,'InsulinDegradation_PlasmaK'))    * Firth1986_paramsCalibrated(j,column-3);
    optimizedParamTemp(ismember(pNames,"HE_positiveK"))                  = optimizedParamTemp(ismember(pNames,'HE_positiveK'))                  * Firth1986_paramsCalibrated(j,column-3);
    % InsulinClearancePersonalScaling
    optimizedParamTemp(ismember(pNames,"InsulinStabilisationMagnitude"))    = optimizedParamTemp(ismember(pNames,'InsulinStabilisationMagnitude'))    * Firth1986_paramsCalibrated(j,column-2);
    optimizedParamTemp(ismember(pNames,"InsulinProductionK"))               = optimizedParamTemp(ismember(pNames,'InsulinProductionK'))               * Firth1986_paramsCalibrated(j,column-2);
    
    optimizedParamTempMeal               = Mealdeclaration(50,0,0,15,optimizedParamTemp,pNames);
    
    SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = 0;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = 0;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate      = model(0:480:480,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimFirthToMeal      = model(1:1:360,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimFirth_meal       = model(360:1:375,SimFirthToMeal.statevalues(end,:),optimizedParamTempMeal);
    SimFirthAfterMeal   = model(375:1:1440,SimFirth_meal.statevalues(end,:),optimizedParamTemp);
    
    simFirth1986.time           = [SimFirthToMeal.time,           SimFirth_meal.time(2:end),           SimFirthAfterMeal.time(2:end)]-360;
    simFirth1986.variablevalues = [SimFirthToMeal.variablevalues; SimFirth_meal.variablevalues(2:end,:); SimFirthAfterMeal.variablevalues(2:end,:)];
    simFirth1986.reactionvalues = [SimFirthToMeal.reactionvalues; SimFirth_meal.reactionvalues(2:end,:); SimFirthAfterMeal.reactionvalues(2:end,:)];
    simFirth1986.variables      = [SimFirth_meal.variables];
    simFirth1986.reactions      = [SimFirth_meal.reactions];
    
    if j == 1
        simBest = simFirth1986;
        
        
        % GlucoseConcentrationPlasma
        maxGlucoseConcentrationPlasma1 = simFirth1986.variablevalues(:,ismember(simFirth1986.variables,'GlucoseConcentrationPlasma'));
        minGlucoseConcentrationPlasma1 = simFirth1986.variablevalues(:,ismember(simFirth1986.variables,'GlucoseConcentrationPlasma'));
        maxGlucoseConcentrationPlasma2 = simFirth1986.variablevalues(:,ismember(simFirth1986.variables,'GlucoseConcentrationPlasma'));
        minGlucoseConcentrationPlasma2 = simFirth1986.variablevalues(:,ismember(simFirth1986.variables,'GlucoseConcentrationPlasma'));
        
        % InsulinConcentrationPlasma
        maxInsulinConcentrationPlasma1 = simFirth1986.variablevalues(:,ismember(simFirth1986.variables,'InsulinConcentrationPlasma'));
        minInsulinConcentrationPlasma1 = simFirth1986.variablevalues(:,ismember(simFirth1986.variables,'InsulinConcentrationPlasma'));
        maxInsulinConcentrationPlasma2 = simFirth1986.variablevalues(:,ismember(simFirth1986.variables,'InsulinConcentrationPlasma'));
        minInsulinConcentrationPlasma2 = simFirth1986.variablevalues(:,ismember(simFirth1986.variables,'InsulinConcentrationPlasma'));
        
        % EGP
        maxEGP1 = simFirth1986.reactionvalues(:,ismember(simFirth1986.reactions,'EGP'));
        minEGP1 = simFirth1986.reactionvalues(:,ismember(simFirth1986.reactions,'EGP'));
        maxEGP2 = simFirth1986.reactionvalues(:,ismember(simFirth1986.reactions,'EGP'));
        minEGP2 = simFirth1986.reactionvalues(:,ismember(simFirth1986.reactions,'EGP'));
        
        
    end
    
    % EGP
    maxEGP1 = simFirth1986.reactionvalues(:,ismember(simFirth1986.reactions,'EGP'));
    minEGP1 = simFirth1986.reactionvalues(:,ismember(simFirth1986.reactions,'EGP'));
    maxEGP2 = max(maxEGP2,maxEGP1);
    minEGP2 = min(minEGP2,minEGP1);
    
    % GlucoseConcentrationPlasma
    maxGlucoseConcentrationPlasma1 = simFirth1986.variablevalues(:,ismember(simFirth1986.variables,'GlucoseConcentrationPlasma'));
    minGlucoseConcentrationPlasma1 = simFirth1986.variablevalues(:,ismember(simFirth1986.variables,'GlucoseConcentrationPlasma'));
    maxGlucoseConcentrationPlasma2 = max(maxGlucoseConcentrationPlasma2,maxGlucoseConcentrationPlasma1);
    minGlucoseConcentrationPlasma2 = min(minGlucoseConcentrationPlasma2,minGlucoseConcentrationPlasma1);
    
    % InsulinConcentrationPlasma
    maxInsulinConcentrationPlasma1 = simFirth1986.variablevalues(:,ismember(simFirth1986.variables,'InsulinConcentrationPlasma'));
    minInsulinConcentrationPlasma1 = simFirth1986.variablevalues(:,ismember(simFirth1986.variables,'InsulinConcentrationPlasma'));
    maxInsulinConcentrationPlasma2 = max(maxInsulinConcentrationPlasma2,maxInsulinConcentrationPlasma1);
    minInsulinConcentrationPlasma2 = min(minInsulinConcentrationPlasma2,minInsulinConcentrationPlasma1);
    
    
    
end

timeFill = [simFirth1986.time, fliplr(simFirth1986.time)];

%% core variables Gly, I, G & EGP

figure('Name', "Firth1986 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% EGP
hold on
set(gca,'xtick',[0,5,10,15],'ytick',[0,0.1,0.2,0.3],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'EGP'; '(g/min)'});
xlabel('Time (h)')
fill(timeFill/60,[maxEGP2', fliplr(minEGP2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot((simBest.time)/60, simBest.reactionvalues(:,ismember(simBest.reactions,'EGP')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
errorbar(Firth1986_data.time_EGP_healthy(1:2)/60,Firth1986_data.EGP_healthy(1:2),Firth1986_data.EGPSEM_healthy(1:2),' k x','MarkerSize',60,'LineWidth',10,'CapSize',40');
errorbar(Firth1986_data.time_EGP_healthy(3:14)/60,Firth1986_data.EGP_healthy(3:14),Firth1986_data.EGPSEM_healthy(3:14),' k .','MarkerSize',1,'LineWidth',10,'CapSize',40');
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
xlim([-1 10])
ylim([0 0.3])
legend('Model uncertainty','Best prediction','Calibration Data','Validation data','OGTT','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize+3)
hold off
grid on

%%
figure('Name', "Firth1986 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Glucose
hold on
set(gca,'xtick',[0,5,10,15],'ytick',[0,3,6,9,12],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Plasma Glucose'; '(mM)'});
xlabel('Time (h)')
fill(timeFill/60,[maxGlucoseConcentrationPlasma2', fliplr(minGlucoseConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot((simBest.time)/60, simBest.variablevalues(:,ismember(simBest.variables,'GlucoseConcentrationPlasma')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
errorbar(Firth1986_data.time_glucose_healthy(1:2)/60,Firth1986_data.glucose_healthy(1:2),Firth1986_data.glucoseSEM_healthy(1:2),' k x','MarkerSize',60,'LineWidth',10,'CapSize',40');
errorbar(Firth1986_data.time_glucose_healthy(3:end)/60,Firth1986_data.glucose_healthy(3:end),Firth1986_data.glucoseSEM_healthy(3:end),' k .','MarkerSize',1,'LineWidth',10,'CapSize',40');
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
xlim([-1 10])
ylim([0 12])
hold off
grid on

%%
figure('Name', "Firth1986 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Insulin
hold on
set(gca,'xtick',[0,5,10,15],'ytick',[0,200,400,600],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Plasma Insulin'; '(pM)'});
xlabel('Time (h)')
fill(timeFill/60,[maxInsulinConcentrationPlasma2', fliplr(minInsulinConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot((simBest.time)/60, simBest.variablevalues(:,ismember(simBest.variables,'InsulinConcentrationPlasma')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
errorbar(Firth1986_data.time_insulin_healthy(1:2)/60,Firth1986_data.insulin_healthy(1:2),Firth1986_data.insulinSEM_healthy(1:2),' k x','MarkerSize',60,'LineWidth',10,'CapSize',40');
errorbar(Firth1986_data.time_insulin_healthy(3:end)/60,Firth1986_data.insulin_healthy(3:end),Firth1986_data.insulinSEM_healthy(3:end),' k .','MarkerSize',1,'LineWidth',10,'CapSize',40');
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
xlim([-1 10])
ylim([0 600])
hold off
grid on
