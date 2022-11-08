function Cost_Model = Rothman1991Personspecific_costfunction(Rothman1991_glycogen,Rothman1991_time,DatapointsConcidered,param,optimizedParamTemp,model, initialvalues, pNames,sNames)

param = exp(param);

optimizedParamTemp(ismember(pNames,"carb_flow"))    = param(1);
optimizedParamTemp(ismember(pNames,"protein_flow")) = param(2);
optimizedParamTemp(ismember(pNames,"lipids_flow"))  = param(3);
optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK"))  = optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK")) * param(4);

%% Declare Rothman

[optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,optimizedParamTemp,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames

%% Simulate Krssak

try
    
    SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = 0;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = 0;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate      = model(0:480:480,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimRothman1991      = model(1:1:4320,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    
catch error
    Cost_Model = 1e60;
    return
end

%% Calculate cost

simGly = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
simGly = simGly(round(Rothman1991_time,0));
SEM = (Rothman1991_glycogen(1:DatapointsConcidered)*0.1)+10;
CostGly = ((Rothman1991_glycogen(1:DatapointsConcidered)- simGly(1:DatapointsConcidered)).^2)./(SEM.^2);
Cost_Model =  sum(CostGly,'omitnan');

%% Print Cost
% if Cost_Model < chi2inv(0.95,11) %&& nargin==4
% fprintf(fid, [sprintf('%.52f, ',params) sprintf('%.52f\n',Cost_Model)]);
% end

end
