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

param = 1:1:4;

% Meals
lb(1)         = log(ParametersInformation.LowerBound(ismember(pNames,'carb_flow')));    %Param 1
ub(1)         = log(ParametersInformation.UpperBound(ismember(pNames,'carb_flow')));
lb(end+1)     = log(ParametersInformation.LowerBound(ismember(pNames,'protein_flow')));  %Param 2
ub(end+1)     = log(ParametersInformation.UpperBound(ismember(pNames,'protein_flow')));
lb(end+1)     = log(ParametersInformation.LowerBound(ismember(pNames,'lipids_flow')));   %Param 3
ub(end+1)     = log(ParametersInformation.UpperBound(ismember(pNames,'lipids_flow')));
lb(end+1)     = log(0.8);                                                                %  Param 4
ub(end+1)     = log(1.2);


%% Optimization
for i = 1:row
    
    % P1
    optimizedParamTemp = EstimationData_params(i,1:AmountParameterOpti);
    func =@(param)Rothman1991Personspecific_costfunction(Rothman1991_data.glycogenp1(1:9),Rothman1991_data.timep1(1:9),1,param,optimizedParamTemp,model,initialvalues, pNames,sNames);
    [Rothman1991_P1params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_P1params = exp(Rothman1991_P1params);
    optimizedParamTemp(ismember(pNames,"carb_flow"))              = Rothman1991_P1params(1);
    optimizedParamTemp(ismember(pNames,"protein_flow"))           = Rothman1991_P1params(2);
    optimizedParamTemp(ismember(pNames,"lipids_flow"))            = Rothman1991_P1params(3);
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,optimizedParamTemp,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    optimizedParamTemp(end+1)                                     = Rothman1991_P1params(4);
    optimizedParamTemp(end+1)                                     = minCostPS;
    Rothman1991_P1paramsCalibrated(i,:) = optimizedParamTemp;
    
    % P2
    optimizedParamTemp = EstimationData_params(i,1:AmountParameterOpti);
    func =@(param)Rothman1991Personspecific_costfunction(Rothman1991_data.glycogenp2(1:11),Rothman1991_data.timep2(1:11),1,param,optimizedParamTemp,model,initialvalues, pNames,sNames);
    [Rothman1991_P1params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_P1params = exp(Rothman1991_P1params);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = Rothman1991_P1params(1);
    optimizedParamTemp(ismember(pNames,"protein_flow")) = Rothman1991_P1params(2);
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = Rothman1991_P1params(3);
    optimizedParamTemp(end+1)                           = Rothman1991_P1params(4);
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,optimizedParamTemp,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    optimizedParamTemp(end+1)                           = Rothman1991_P1params(4);
    optimizedParamTemp(end+1) = minCostPS;
    Rothman1991_P2paramsCalibrated(i,:) = optimizedParamTemp;
    
    % P3
    optimizedParamTemp = EstimationData_params(i,1:AmountParameterOpti);
    func =@(param)Rothman1991Personspecific_costfunction(Rothman1991_data.glycogenp3(1:10),Rothman1991_data.timep3(1:10),1,param,optimizedParamTemp,model,initialvalues, pNames,sNames);
    [Rothman1991_P1params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_P1params = exp(Rothman1991_P1params);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = Rothman1991_P1params(1);
    optimizedParamTemp(ismember(pNames,"protein_flow")) = Rothman1991_P1params(2);
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = Rothman1991_P1params(3);
    optimizedParamTemp(end+1)                           = Rothman1991_P1params(4);
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,optimizedParamTemp,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    optimizedParamTemp(end+1)                           = Rothman1991_P1params(4);
    optimizedParamTemp(end+1) = minCostPS;
    Rothman1991_P3paramsCalibrated(i,:) = optimizedParamTemp;
    
    % P4
    optimizedParamTemp = EstimationData_params(i,1:AmountParameterOpti);
    func =@(param)Rothman1991Personspecific_costfunction(Rothman1991_data.glycogenp4(1:10),Rothman1991_data.timep4(1:10),1,param,optimizedParamTemp,model,initialvalues, pNames,sNames);
    [Rothman1991_P1params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_P1params = exp(Rothman1991_P1params);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = Rothman1991_P1params(1);
    optimizedParamTemp(ismember(pNames,"protein_flow")) = Rothman1991_P1params(2);
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = Rothman1991_P1params(3);
    optimizedParamTemp(end+1)                           = Rothman1991_P1params(4);
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,optimizedParamTemp,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    optimizedParamTemp(end+1)                           = Rothman1991_P1params(4);
    optimizedParamTemp(end+1) = minCostPS;
    Rothman1991_P4paramsCalibrated(i,:) = optimizedParamTemp;
    
    % P5
    optimizedParamTemp = EstimationData_params(i,1:AmountParameterOpti);
    func =@(param)Rothman1991Personspecific_costfunction(Rothman1991_data.glycogenp5(1:8),Rothman1991_data.timep5(1:8),1,param,optimizedParamTemp,model,initialvalues, pNames,sNames);
    [Rothman1991_P1params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_P1params = exp(Rothman1991_P1params);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = Rothman1991_P1params(1);
    optimizedParamTemp(ismember(pNames,"protein_flow")) = Rothman1991_P1params(2);
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = Rothman1991_P1params(3);
    optimizedParamTemp(end+1)                           = Rothman1991_P1params(4);
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,optimizedParamTemp,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    optimizedParamTemp(end+1)                           = Rothman1991_P1params(4);
    optimizedParamTemp(end+1) = minCostPS;
    Rothman1991_P5paramsCalibrated(i,:) = optimizedParamTemp;
    
    % P6
    optimizedParamTemp = EstimationData_params(i,1:AmountParameterOpti);
    func =@(param)Rothman1991Personspecific_costfunction(Rothman1991_data.glycogenp6(1:11),Rothman1991_data.timep6(1:11),1,param,optimizedParamTemp,model,initialvalues, pNames,sNames);
    [Rothman1991_P1params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_P1params = exp(Rothman1991_P1params);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = Rothman1991_P1params(1);
    optimizedParamTemp(ismember(pNames,"protein_flow")) = Rothman1991_P1params(2);
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = Rothman1991_P1params(3);
    optimizedParamTemp(end+1)                           = Rothman1991_P1params(4);
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,optimizedParamTemp,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    optimizedParamTemp(end+1)                           = Rothman1991_P1params(4);
    optimizedParamTemp(end+1) = minCostPS;
    Rothman1991_P6paramsCalibrated(i,:) = optimizedParamTemp;
    
    % P7
    optimizedParamTemp = EstimationData_params(i,1:AmountParameterOpti);
    func =@(param)Rothman1991Personspecific_costfunction(Rothman1991_data.glycogenp7(1:9),Rothman1991_data.timep7(1:9),1,param,optimizedParamTemp,model,initialvalues, pNames,sNames);
    [Rothman1991_P1params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
    Rothman1991_P1params = exp(Rothman1991_P1params);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = Rothman1991_P1params(1);
    optimizedParamTemp(ismember(pNames,"protein_flow")) = Rothman1991_P1params(2);
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = Rothman1991_P1params(3);
    optimizedParamTemp(end+1)                           = Rothman1991_P1params(4);
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,optimizedParamTemp,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    optimizedParamTemp(end+1)                           = Rothman1991_P1params(4);
    optimizedParamTemp(end+1) = minCostPS;
    Rothman1991_P7paramsCalibrated(i,:) = optimizedParamTemp;
    
    i = i
end

save(['Rothman1991_P1paramsCalibrated'],'Rothman1991_P1paramsCalibrated');
save(['Rothman1991_P2paramsCalibrated'],'Rothman1991_P2paramsCalibrated');
save(['Rothman1991_P3paramsCalibrated'],'Rothman1991_P3paramsCalibrated');
save(['Rothman1991_P4paramsCalibrated'],'Rothman1991_P4paramsCalibrated');
save(['Rothman1991_P5paramsCalibrated'],'Rothman1991_P5paramsCalibrated');
save(['Rothman1991_P6paramsCalibrated'],'Rothman1991_P6paramsCalibrated');
save(['Rothman1991_P7paramsCalibrated'],'Rothman1991_P7paramsCalibrated');

end

