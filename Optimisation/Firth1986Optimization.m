function [minCostPS] = Firth1986Optimization(seed,Firth1986_data, model,ParametersInformation,AmountParameterOpti,options2)

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

%% Bounds

param = 1:1:7;

% Meals
lb(1)         = log(ParametersInformation.LowerBound(ismember(pNames,'carb_flow')));    %Param 1
ub(1)         = log(ParametersInformation.UpperBound(ismember(pNames,'carb_flow')));
lb(end+1)     = log(ParametersInformation.LowerBound(ismember(pNames,'protein_flow')));  %Param 2
ub(end+1)     = log(ParametersInformation.UpperBound(ismember(pNames,'protein_flow')));
lb(end+1)     = log(ParametersInformation.LowerBound(ismember(pNames,'lipids_flow')));   %Param 3
ub(end+1)     = log(ParametersInformation.UpperBound(ismember(pNames,'lipids_flow')));

% Insulin resistance: Personal Scaling -  Param 4
lb(end+1)     = log(0.7);
ub(end+1)     = log(1.5);

% Insulin Clearance: Personal Scaling-  Param 5
lb(end+1)     = log(0.7);
ub(end+1)     = log(1.5);

% Insulin Repesponse: Personal Scaling-  Param 6
lb(end+1)     = log(0.7);
ub(end+1)     = log(1.5);

% Insulin response glucose utliziation: Personal Scaling -  Param 7
lb(end+1)     = log(1);
ub(end+1)     = log(2);


%% Optimization
for i = 1:row

    optimizedParamTemp = EstimationData_params(i,1:AmountParameterOpti);
    
    func =@(param)Firth1986_costfunction(Firth1986_data,2,2,2,param,optimizedParamTemp,model,initialvalues, pNames,sNames);
    
    [Firth1986_params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Firth1986_params = exp(Firth1986_params);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = Firth1986_params(1);
    optimizedParamTemp(ismember(pNames,"protein_flow")) = Firth1986_params(2);
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = Firth1986_params(3);
    [optimizedParamTemp,initialvalues]   = Persondeclaration(175, 0, 1,75,25,initialvalues,sNames,optimizedParamTemp,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    optimizedParamTemp(end+1)                           = Firth1986_params(4);
    optimizedParamTemp(end+1)                           = Firth1986_params(5);
    optimizedParamTemp(end+1)                           = Firth1986_params(6);
    optimizedParamTemp(end+1)                           = Firth1986_params(7);
    optimizedParamTemp(end+1) = minCostPS;
    Firth1986_paramsCalibrated(i,:) = optimizedParamTemp;
    
    i = i
end

save(['Firth1986_paramsCalibrated'],'Firth1986_paramsCalibrated');

end

