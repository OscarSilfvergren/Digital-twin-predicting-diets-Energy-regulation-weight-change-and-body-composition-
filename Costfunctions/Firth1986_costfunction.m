function Cost_Model = Firth1986_costfunction(Firth1986_data,EGPDatapoints,GDatapoints,IDatapoints,param,optimizedParamTemp,model, initialvalues, pNames,sNames)

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

%% Declare Firth

[optimizedParamTemp,initialvalues]   = Persondeclaration(175, 0, 1,75,25,initialvalues,sNames,optimizedParamTemp,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames

%% Simulate Krssak

optimizedParamTempMeal               = Mealdeclaration(50,0,0,15,optimizedParamTemp,pNames);         %Fett är fel


try
    
    SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = 0;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = 0;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate      = model(0:480:480,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimFirthToMeal      = model(1:1:360,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimFirth_meal       = model(360:1:375,SimFirthToMeal.statevalues(end,:),optimizedParamTempMeal);
    SimFirthAfterMeal   = model(375:1:1440,SimFirth_meal.statevalues(end,:),optimizedParamTemp);
    
catch error
    Cost_Model = 1e60;
    return
end

simFirth1986.time           = 1:1:1440;
simFirth1986.variablevalues = [SimFirthToMeal.variablevalues; SimFirth_meal.variablevalues(2:end,:); SimFirthAfterMeal.variablevalues(2:end,:)];
simFirth1986.reactionvalues = [SimFirthToMeal.reactionvalues; SimFirth_meal.reactionvalues(2:end,:); SimFirthAfterMeal.reactionvalues(2:end,:)];
simFirth1986.variables      = [SimFirth_meal.variables];
simFirth1986.reactions      = [SimFirth_meal.reactions];

%% Calculate cost

% EGP
simEGP = simFirth1986.reactionvalues(:,ismember(simFirth1986.reactions,'EGP'));
simEGP = simEGP(round(Firth1986_data.time_EGP_healthy(1:EGPDatapoints)+360,0));
CostGly = ((Firth1986_data.EGP_healthy(1:EGPDatapoints)- simEGP).^2)./(Firth1986_data.EGPSEM_healthy(1:EGPDatapoints).^2);
Cost_Model =  sum(CostGly,'omitnan');

%G
simG   = simFirth1986.variablevalues(:,ismember(simFirth1986.variables,'GlucoseConcentrationPlasma'));
simG   = simG(round(Firth1986_data.time_glucose_healthy(1:GDatapoints)+360,0));
CostG   = ((Firth1986_data.glucose_healthy(1:GDatapoints) - simG).^2)  ./(Firth1986_data.glucoseSEM_healthy(1:GDatapoints).^2);
Cost_Model =  Cost_Model + sum(CostG,'omitnan');

%I
simI   = simFirth1986.variablevalues(:,ismember(simFirth1986.variables,'InsulinConcentrationPlasma'));
simI   = simI(round(Firth1986_data.time_insulin_healthy(1:IDatapoints)+360,0));
CostI   = ((Firth1986_data.insulin_healthy(1:IDatapoints) - simI).^2)  ./((Firth1986_data.insulinSEM_healthy(1:IDatapoints)).^2);
Cost_Model =  Cost_Model + sum(CostI,'omitnan');

%% Print Cost
% if Cost_Model < chi2inv(0.95,11) %&& nargin==4
% fprintf(fid, [sprintf('%.52f, ',params) sprintf('%.52f\n',Cost_Model)]);
% end

end
