function Cost_Model = Bray2012_costfunction(Bray2012_data,DatapointsToConcider,param,optimizedParamTemp,model, initialvalues, pNames,sNames)

param = exp(param);

% Meal
optimizedParamTemp(ismember(pNames,"carb_flow"))    = param(1);
optimizedParamTemp(ismember(pNames,"protein_flow")) = param(2);
optimizedParamTemp(ismember(pNames,"lipids_flow"))  = param(3);

% BMR
optimizedParamTemp(ismember(pNames,"BMRDynamicBW"))    = optimizedParamTemp(ismember(pNames,'BMRDynamicBW'))  * param(4);
optimizedParamTemp(ismember(pNames,"BMRDynamicGly"))   = optimizedParamTemp(ismember(pNames,'BMRDynamicGly')) * param(4);
optimizedParamTemp(ismember(pNames,"BMRbasal"))        = optimizedParamTemp(ismember(pNames,'BMRbasal'))      * param(4);

%% Declare Krssak

[optimizedParamTemp,initialvalues]    = Persondeclaration(166, 0, 1,73.5,31.5,initialvalues,sNames,optimizedParamTemp,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames

%% Simulate Krssak

  try
      
SimSteadystate     = model(0:282240:282240,initialvalues,optimizedParamTemp);
BMR                = SimSteadystate.reactionvalues(:,ismember(SimSteadystate.reactions,'BMR'));
IntakeFoodCalories = BMR(end)-600;
optimizedParamTemp(ismember(pNames,"carb_flow"))    = ((IntakeFoodCalories*0.5)/4)/1440;
optimizedParamTemp(ismember(pNames,"protein_flow")) = ((IntakeFoodCalories*0.25)/4)/1440;
optimizedParamTemp(ismember(pNames,"lipids_flow"))  = ((IntakeFoodCalories*0.25)/9)/1440;
SimSteadystate     = model(0:2880:2880,SimSteadystate.statevalues(end,:),optimizedParamTemp);
SimBray2012        = model(0:1440:241920,SimSteadystate.statevalues(end,:),optimizedParamTemp); % 24 weeks

catch error
    Cost_Model = 1e60;
    return
  end

  %% Extract measured variables from Simulation

simBW = SimBray2012.variablevalues(:,ismember(SimBray2012.variables,'BW'));
simBW = simBW(Bray2012_data.time);
CostBW = ((Bray2012_data.weight(1:DatapointsToConcider)- simBW(1:DatapointsToConcider)).^2)./((Bray2012_data.weight(1:DatapointsToConcider)*0.15).^2);
Cost_Model =  sum(CostBW,'omitnan');

end
