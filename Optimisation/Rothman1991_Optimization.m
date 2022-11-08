function [minCostPS] = Rothman1991_Optimization(seed,Rothman1991_data, model,ParametersInformation,AmountParameterOpti,options2)

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

    param = 1:1:3;
    
    % Meals
    lb(1)         = log(ParametersInformation.LowerBound(ismember(pNames,'carb_flow')));    %Param 1
    ub(1)         = log(ParametersInformation.UpperBound(ismember(pNames,'carb_flow')));
    lb(end+1)     = log(ParametersInformation.LowerBound(ismember(pNames,'protein_flow')));  %Param 2
    ub(end+1)     = log(ParametersInformation.UpperBound(ismember(pNames,'protein_flow')));
    lb(end+1)     = log(ParametersInformation.LowerBound(ismember(pNames,'lipids_flow')));   %Param 3
    ub(end+1)     = log(ParametersInformation.UpperBound(ismember(pNames,'lipids_flow')));
    
    %     % BMR: Personal Scaling -  Param 4
%     lb(end+1)     = log(0.6);
%     ub(end+1)     = log(1.4);
%     
%     % Insulin resistance: Personal Scaling -  Param 5
%     lb(end+1)     = log(0.6);
%     ub(end+1)     = log(1.4);
%     
%     % Insulin Clearance: Personal Scaling-  Param 6
%     lb(end+1)     = log(0.6);
%     ub(end+1)     = log(1.4);
%     
%     % Insulin Repesponse: Personal Scaling-  Param 7
%     lb(end+1)     = log(0.6);
%     ub(end+1)     = log(1.4);

    
%% Optimization
for i = 1:row
    
    optimizedParamTemp = EstimationData_params(i,1:AmountParameterOpti);
    
    func =@(param)Rothman1991_costfunction(Rothman1991_data,1,param,optimizedParamTemp,model,initialvalues, pNames,sNames);
    
    [Rothman1991_params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_params = exp(Rothman1991_params);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = Rothman1991_params(1);
    optimizedParamTemp(ismember(pNames,"protein_flow")) = Rothman1991_params(2);
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = Rothman1991_params(3);
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,optimizedParamTemp,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    optimizedParamTemp(end+1) = minCostPS;
    Rothman1991_paramsCalibrated(i,:) = optimizedParamTemp;
    
    i = i
end

save(['Rothman1991_paramsCalibrated'],'Rothman1991_paramsCalibrated');

end

