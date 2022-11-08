function [minCostPS] = Bray2012Optimization(seed,Bray2012_data, model,ParametersInformation,AmountParameterOpti,options2)
load('Bray2012_data');

rng('shuffle') % set the stochastic seed
s = rng;
rng(s.Seed+seed);

%% Extract optimized parameters

EstimationData_params = readtable('FullParameterEstimation_parameters.csv');
EstimationData_params = table2array(EstimationData_params);
[row column] = size(EstimationData_params);
EstimationData_params = sortrows(EstimationData_params,column);
%% Pre-Optimization

[pNames, param]               = IQMparameters(model);
[sNames, ODEs, initialvalues] = IQMstates(model);

%% Clean data

Bray2012_data.time = Bray2012_data.time/1440;
Bray2012_data.time(1) = 1;
EstimationData_params(ismember(pNames,"EnergyExpandature"))       = 0;
%% Bounds

param = 1:4;

% Meals
lb(1)         = log(ParametersInformation.LowerBound(ismember(pNames,'carb_flow')));    %Param 1
ub(1)         = log(ParametersInformation.UpperBound(ismember(pNames,'carb_flow')));
lb(end+1)     = log(ParametersInformation.LowerBound(ismember(pNames,'protein_flow')));  %Param 2
ub(end+1)     = log(ParametersInformation.UpperBound(ismember(pNames,'protein_flow')));
lb(end+1)     = log(ParametersInformation.LowerBound(ismember(pNames,'lipids_flow')));   %Param 3
ub(end+1)     = log(ParametersInformation.UpperBound(ismember(pNames,'lipids_flow')));

%BMR: Personal Scaling -  Param 4
lb(end+1)     = log(0.5);
ub(end+1)     = log(1.2);

%% Optimization

for i = 1:row
    
    optimizedParamTemp = EstimationData_params(i,1:AmountParameterOpti);
    
    func =@(param)Bray2012_costfunction(Bray2012_data,2,param,optimizedParamTemp,model, initialvalues, pNames,sNames);
    
    [Bray2012_params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Bray2012_params = exp(Bray2012_params);
    
    % Meal
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = Bray2012_params(1);
    optimizedParamTemp(ismember(pNames,"protein_flow")) = Bray2012_params(2);
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = Bray2012_params(3);
    % BMR
    optimizedParamTemp(ismember(pNames,"BMRDynamicBW"))    = optimizedParamTemp(ismember(pNames,'BMRDynamicBW'))  * Bray2012_params(4);
    optimizedParamTemp(ismember(pNames,"BMRDynamicGly"))   = optimizedParamTemp(ismember(pNames,'BMRDynamicGly')) * Bray2012_params(4);
    optimizedParamTemp(ismember(pNames,"BMRbasal"))        = optimizedParamTemp(ismember(pNames,'BMRbasal'))      * Bray2012_params(4);
    [Bray2012_params,initialvalues]    = Persondeclaration(166, 0, 1,73.5,31.5,initialvalues,sNames,optimizedParamTemp,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    Bray2012_params(end+1) = Bray2012_params(4);
    Bray2012_params(end+1) = minCostPS;
    Bray2012_paramsCalibrated(i,:) = Bray2012_params;
    i = i
    
end

save(['Bray2012_paramsCalibrated'],'Bray2012_paramsCalibrated');

end

