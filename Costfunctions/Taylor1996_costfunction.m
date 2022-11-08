function Cost_Model = Taylor1996_costfunction(Taylor1996_data,GlyDatapoints,GDatapoints,IDatapoints,FFADatapoints,param,optimizedParamTemp,model, initialvalues, pNames,sNames)

param = exp(param);

optimizedParamTemp(ismember(pNames,"carb_flow"))    = param(1);
optimizedParamTemp(ismember(pNames,"protein_flow")) = param(2);
optimizedParamTemp(ismember(pNames,"lipids_flow"))  = param(3);
% InsulinResistancePersonalScaling
optimizedParamTemp(ismember(pNames,"Vm_id"))                     = optimizedParamTemp(ismember(pNames,'Vm_id'))                   * param(7);
optimizedParamTemp(ismember(pNames,"Vf_id"))                     = optimizedParamTemp(ismember(pNames,'Vf_id'))                   * param(7);
optimizedParamTemp(ismember(pNames,"InsulinResponseFATP1_K"))   = optimizedParamTemp(ismember(pNames,'InsulinResponseFATP1_K')) * param(4);
optimizedParamTemp(ismember(pNames,"InsulinResponseLiver"))     = optimizedParamTemp(ismember(pNames,'InsulinResponseLiver'))   * param(4);
% InsulinClearancePersonalScaling
optimizedParamTemp(ismember(pNames,"InsulinDegradation_PlasmaK"))    = optimizedParamTemp(ismember(pNames,'InsulinDegradation_PlasmaK'))    * param(5);
optimizedParamTemp(ismember(pNames,"HE_positiveK"))                  = optimizedParamTemp(ismember(pNames,'HE_positiveK'))                  * param(5);
% InsulinClearancePersonalScaling
optimizedParamTemp(ismember(pNames,"InsulinStabilisationMagnitude"))    = optimizedParamTemp(ismember(pNames,'InsulinStabilisationMagnitude'))    * param(6);
optimizedParamTemp(ismember(pNames,"InsulinProductionK"))               = optimizedParamTemp(ismember(pNames,'InsulinProductionK'))               * param(6);

%% Declare Taylor

[optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5.6,initialvalues,sNames,optimizedParamTemp,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames

%% Simulate Krssak

optimizedParamTempMeal               = Mealdeclaration(138.64,29.25,16.9377777777778,15,optimizedParamTemp,pNames);         %Fett är fel


try
    
    SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = 0;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = 0;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate      = model(0:480:480,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimTaylorToMeal      = model(1:1:360,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimTaylor_meal       = model(360:1:375,SimTaylorToMeal.statevalues(end,:),optimizedParamTempMeal);
    SimTaylorAfterMeal   = model(375:1:1440,SimTaylor_meal.statevalues(end,:),optimizedParamTemp);
    
catch error
    Cost_Model = 1e60;
    return
end

simTaylor1996.time           = 1:1:1440;
simTaylor1996.variablevalues = [SimTaylorToMeal.variablevalues; SimTaylor_meal.variablevalues(2:end,:); SimTaylorAfterMeal.variablevalues(2:end,:)];
simTaylor1996.reactionvalues = [SimTaylorToMeal.reactionvalues; SimTaylor_meal.reactionvalues(2:end,:); SimTaylorAfterMeal.reactionvalues(2:end,:)];
simTaylor1996.variables      = [SimTaylor_meal.variables];
simTaylor1996.reactions      = [SimTaylor_meal.reactions];

%% Calculate cost

% Gly
simGly = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'GlycogenConcentrationliver'));
simGly = simGly(round(Taylor1996_data.time_glycogen_healthy(1:GlyDatapoints)+360,0));
CostGly = ((Taylor1996_data.glycogen_healthy(1:GlyDatapoints)- simGly).^2)./(Taylor1996_data.glycogenSEM_healthy(1:GlyDatapoints).^2);
Cost_Model =  sum(CostGly,'omitnan');

%G
simG   = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'GlucoseConcentrationPlasma'));
simG   = simG(Taylor1996_data.time_glucose_healthy(1:GDatapoints)+360);
CostG   = ((Taylor1996_data.glucose_healthy(1:GDatapoints) - simG).^2)  ./(Taylor1996_data.glucoseSEM_healthy(1:GDatapoints).^2);
Cost_Model =  Cost_Model + sum(CostG,'omitnan');

%I
simI   = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'InsulinConcentrationPlasma'));
simI   = simI(Taylor1996_data.time_Insulin_healthy(1:IDatapoints)+360);
CostI   = ((Taylor1996_data.Insulin_healthy(1:IDatapoints) - simI).^2)  ./((Taylor1996_data.InsulinSEM_healthy(1:IDatapoints)).^2);
Cost_Model =  Cost_Model + sum(CostI,'omitnan');

%FFA
simFFA   = simTaylor1996.variablevalues(:,ismember(simTaylor1996.variables,'FFAConcentrationPlasma'));
simFFA   = simFFA(Taylor1996_data.time_FFA_healthy(1:FFADatapoints)+360);
CostFFA   = ((Taylor1996_data.FFA_healthy(1:FFADatapoints) - simFFA).^2)  ./((Taylor1996_data.FFASEM_healthy(1:FFADatapoints)).^2);
Cost_Model =  Cost_Model + sum(CostFFA,'omitnan');

%% Print Cost
% if Cost_Model < chi2inv(0.95,11) %&& nargin==4
% fprintf(fid, [sprintf('%.52f, ',params) sprintf('%.52f\n',Cost_Model)]);
% end

end
