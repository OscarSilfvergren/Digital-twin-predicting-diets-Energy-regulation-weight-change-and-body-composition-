load('Lerche2009_paramsCalibrated')

[sNames, ODEs, initialvalues] = IQMstates(model);
[pNames, param] = IQMparameters(model);

[row column] = size(Lerche2009_paramsCalibrated);
Lerche2009_paramsCalibrated = sortrows(Lerche2009_paramsCalibrated,column);

%% Cost to full estimation data

for j = 1:row
    
    LercheTemp              = Lerche2009_paramsCalibrated(j,:);
    PersonSpecificParams(1)  = LercheTemp(ismember(pNames,"carb_flow"));
    PersonSpecificParams(2)  = LercheTemp(ismember(pNames,"protein_flow"));
    PersonSpecificParams(3)  = LercheTemp(ismember(pNames,"lipids_flow"));
    PersonSpecificParams(4)  = LercheTemp(column-3);
    PersonSpecificParams(5)  = LercheTemp(column-2);
    PersonSpecificParams(6)  = LercheTemp(column-1);
    costValidationLerche = Lerche2009_costfunction(Lerche2009_data,56,21,log(PersonSpecificParams),LercheTemp(1:AmountParameterOpti),model,initialvalues, pNames,sNames);
    
    Lerche2009_paramsCalibrated(j,column) = costValidationLerche;
end

Lerche2009_paramsCalibrated = sortrows(Lerche2009_paramsCalibrated,column);

disp(' ')
if Lerche2009_paramsCalibrated(1,column) < chi2inv(0.95,517)
    disp('----- Fit to Estimation data is below statistical treshold -----')
end
fprintf('Best prediction of Rothman data: %.2f, Statistical Limit: %.2f (dgf = %i)', Lerche2009_paramsCalibrated(1,column), chi2inv(0.95,11), 11) % Finddatapoints
disp(' ')


%% Simulate Lerche2009


for j = 1:row
    
    
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,75,10,initialvalues,sNames,Lerche2009_paramsCalibrated(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
    % InsulinResistancePersonalScaling
    optimizedParamTemp(ismember(pNames,"Vf_id"))                     = optimizedParamTemp(ismember(pNames,'Vf_id'))                 * Lerche2009_paramsCalibrated(j,end-3);
    optimizedParamTemp(ismember(pNames,"Vm_id"))                     = optimizedParamTemp(ismember(pNames,'Vm_id'))                 * Lerche2009_paramsCalibrated(j,end-3);
    optimizedParamTemp(ismember(pNames,"InsulinResponseFATP1_K"))   = optimizedParamTemp(ismember(pNames,'InsulinResponseFATP1_K')) * Lerche2009_paramsCalibrated(j,end-3);
    optimizedParamTemp(ismember(pNames,"InsulinResponseLiver"))     = optimizedParamTemp(ismember(pNames,'InsulinResponseLiver'))   * Lerche2009_paramsCalibrated(j,end-3);
    % InsulinClearancePersonalScaling
    optimizedParamTemp(ismember(pNames,"InsulinDegradation_PlasmaK"))    = optimizedParamTemp(ismember(pNames,'InsulinDegradation_PlasmaK'))    * Lerche2009_paramsCalibrated(j,end-2);
    optimizedParamTemp(ismember(pNames,"HE_positiveK"))                  = optimizedParamTemp(ismember(pNames,'HE_positiveK'))                  * Lerche2009_paramsCalibrated(j,end-2);
    % InsulinClearancePersonalScaling
    optimizedParamTemp(ismember(pNames,"InsulinStabilisationMagnitude"))    = optimizedParamTemp(ismember(pNames,'InsulinStabilisationMagnitude'))    * Lerche2009_paramsCalibrated(j,end-1);
    optimizedParamTemp(ismember(pNames,"InsulinProductionK"))               = optimizedParamTemp(ismember(pNames,'InsulinProductionK'))               * Lerche2009_paramsCalibrated(j,end-1);
    
    optimizedParamTempMeal               = Mealdeclaration(70,0,0,10,optimizedParamTemp,pNames);
    
    SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = 0;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = 0;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate      = model(0:480:480,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimLercheToMeal      = model(1:1:3000,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimLerche_meal       = model(3000:1:3010,SimLercheToMeal.statevalues(end,:),optimizedParamTempMeal);
    SimLercheAfterMeal   = model(3010:1:3500,SimLerche_meal.statevalues(end,:),optimizedParamTemp);
    
    simLerche2009.time           = [SimLercheToMeal.time,           SimLerche_meal.time(2:end),           SimLercheAfterMeal.time(2:end)]-3000;
    simLerche2009.variablevalues = [SimLercheToMeal.variablevalues; SimLerche_meal.variablevalues(2:end,:); SimLercheAfterMeal.variablevalues(2:end,:)];
    simLerche2009.reactionvalues = [SimLercheToMeal.reactionvalues; SimLerche_meal.reactionvalues(2:end,:); SimLercheAfterMeal.reactionvalues(2:end,:)];
    simLerche2009.variables      = [SimLerche_meal.variables];
    simLerche2009.reactions      = [SimLerche_meal.reactions];
    
    if j == 1
        simBest = simLerche2009;
        
        % GlucoseConcentrationPlasma
        maxGlucoseConcentrationPlasma1 = simLerche2009.variablevalues(:,ismember(simLerche2009.variables,'GlucoseConcentrationPlasma'));
        minGlucoseConcentrationPlasma1 = simLerche2009.variablevalues(:,ismember(simLerche2009.variables,'GlucoseConcentrationPlasma'));
        maxGlucoseConcentrationPlasma2 = simLerche2009.variablevalues(:,ismember(simLerche2009.variables,'GlucoseConcentrationPlasma'));
        minGlucoseConcentrationPlasma2 = simLerche2009.variablevalues(:,ismember(simLerche2009.variables,'GlucoseConcentrationPlasma'));
        
        % InsulinConcentrationPlasma
        maxInsulinConcentrationPlasma1 = simLerche2009.variablevalues(:,ismember(simLerche2009.variables,'InsulinConcentrationPlasma'));
        minInsulinConcentrationPlasma1 = simLerche2009.variablevalues(:,ismember(simLerche2009.variables,'InsulinConcentrationPlasma'));
        maxInsulinConcentrationPlasma2 = simLerche2009.variablevalues(:,ismember(simLerche2009.variables,'InsulinConcentrationPlasma'));
        minInsulinConcentrationPlasma2 = simLerche2009.variablevalues(:,ismember(simLerche2009.variables,'InsulinConcentrationPlasma'));
        
        
        
    end
    
    % GlucoseConcentrationPlasma
    maxGlucoseConcentrationPlasma1 = simLerche2009.variablevalues(:,ismember(simLerche2009.variables,'GlucoseConcentrationPlasma'));
    minGlucoseConcentrationPlasma1 = simLerche2009.variablevalues(:,ismember(simLerche2009.variables,'GlucoseConcentrationPlasma'));
    maxGlucoseConcentrationPlasma2 = max(maxGlucoseConcentrationPlasma2,maxGlucoseConcentrationPlasma1);
    minGlucoseConcentrationPlasma2 = min(minGlucoseConcentrationPlasma2,minGlucoseConcentrationPlasma1);
    
    % InsulinConcentrationPlasma
    maxInsulinConcentrationPlasma1 = simLerche2009.variablevalues(:,ismember(simLerche2009.variables,'InsulinConcentrationPlasma'));
    minInsulinConcentrationPlasma1 = simLerche2009.variablevalues(:,ismember(simLerche2009.variables,'InsulinConcentrationPlasma'));
    maxInsulinConcentrationPlasma2 = max(maxInsulinConcentrationPlasma2,maxInsulinConcentrationPlasma1);
    minInsulinConcentrationPlasma2 = min(minInsulinConcentrationPlasma2,minInsulinConcentrationPlasma1);
    
    
end

timeFill = [simLerche2009.time, fliplr(simLerche2009.time)];

%%
figure('Name', "Lerche2009 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Glucose
hold on
set(gca,'xtick',[-2,-1,0],'ytick',[0,3,4,5,6,9,12],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Plasma Glucose (mmol/L)'});
xlabel('Time (days)')
fill(timeFill/1440,[maxGlucoseConcentrationPlasma2', fliplr(minGlucoseConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
errorbar(Lerche2009_data.time_glucose_healthy(1:5)/1440,Lerche2009_data.glucose_healthy(1:5),Lerche2009_data.glucoseSEM_healthy(1:5),' k x','MarkerSize',60,'LineWidth',10,'CapSize',40');
errorbar(Lerche2009_data.time_glucose_healthy(6:end)/1440,Lerche2009_data.glucose_healthy(6:end),Lerche2009_data.glucoseSEM_healthy(6:end),' k .','MarkerSize',1,'LineWidth',10,'CapSize',40');
plot((simBest.time)/1440, simBest.variablevalues(:,ismember(simBest.variables,'GlucoseConcentrationPlasma')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
legend('Model uncertainty','Best prediction','Calibration Data','Validation data','OGTT','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize+3)
line([0.0 0.02], [0 0],'Color','k','LineWidth',LineWidthMeal);
xlim([-2 0.3])
ylim([0 12])
hold off
grid on

%%
figure('Name', "Lerche2009 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Insulin
hold on
set(gca,'xtick',[-2,-1,0],'ytick',[0,200,400,600],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Plasma Insulin (pmol/L)'});
xlabel('Time (days)')
fill(timeFill/1440,[maxInsulinConcentrationPlasma2', fliplr(minInsulinConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
errorbar(Lerche2009_data.time_insulin_healthy(1:2)/1440,Lerche2009_data.insulin_healthy(1:2),Lerche2009_data.insulinSEM_healthy(1:2),' k x','MarkerSize',60,'LineWidth',10,'CapSize',40');
errorbar(Lerche2009_data.time_insulin_healthy(3:end)/1440,Lerche2009_data.insulin_healthy(3:end),Lerche2009_data.insulinSEMFIXED_healthy(3:end),' k .','MarkerSize',1,'LineWidth',10,'CapSize',40');
plot((simBest.time)/1440, simBest.variablevalues(:,ismember(simBest.variables,'InsulinConcentrationPlasma')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
line([0.0 0.02], [0 0],'Color','k','LineWidth',LineWidthMeal);
xlim([-2 0.3])
ylim([0 600])
hold off
grid on
