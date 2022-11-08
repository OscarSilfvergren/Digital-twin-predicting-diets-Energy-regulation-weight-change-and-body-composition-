function [minCostPS] = FullParameterEstimationOptimization(seed,Silfvergren2022Weight_data,Silfvergren2022GlucoseBloodSensor_data,Silfvergren2022Ketones_data,Silfvergren2022_meals,Krssak2004Healthy_data,DallaMan2007_data,model,ParametersInformation,AmountParameterOpti,options2,fid);

rng('shuffle') % set the stochastic seed
s = rng;
rng(s.Seed+seed);

%% Clean data for lower dataresolution

Silfvergren2022_meals.MealStart           = round(Silfvergren2022_meals.MealStart/120,0)*120;
Silfvergren2022_meals.MealDuration(1:end) = [120];
Silfvergren2022_meals.MealEnd             = Silfvergren2022_meals.MealStart + 120;
Silfvergren2022Weight_data.Time                       = round(Silfvergren2022Weight_data.Time/120,0);
Silfvergren2022Ketones_data.Time_minutesKetones(1:66) = round(Silfvergren2022Ketones_data.Time_minutesKetones(1:66)/120,0);
Silfvergren2022GlucoseBloodSensor_data.glucoseTime    = round(Silfvergren2022GlucoseBloodSensor_data.glucoseTime/120,0);

%% Pre-Optimization
[pNames, params]               = IQMparameters(model);
[sNames, ODEs, initialvalues] = IQMstates(model);

params = 1:1:AmountParameterOpti;
params = params';
lb     = log(ParametersInformation.LowerBound(1:AmountParameterOpti));                                       
ub     = log(ParametersInformation.UpperBound(1:AmountParameterOpti));

%% Add parameters for multiple starts for 3 studies

AmountExtraParametersOptimized = 7;
AmountExtraStudies             = 2;


for i = 2:AmountExtraStudies+1
params(end+1:end+AmountExtraParametersOptimized)  = [i];

% Meals
lb(end+1)     = lb(ismember(ParametersInformation.ParameterName,'carb_flow'));    %Param 1                                                                
ub(end+1)     = ub(ismember(ParametersInformation.ParameterName,'carb_flow')); 
lb(end+1)     = lb(ismember(ParametersInformation.ParameterName,'protein_flow'));          %Param 2                                                           
ub(end+1)     = ub(ismember(ParametersInformation.ParameterName,'protein_flow')); 
lb(end+1)     = lb(ismember(ParametersInformation.ParameterName,'lipids_flow'));                    %Param 3                                                 
ub(end+1)     = ub(ismember(ParametersInformation.ParameterName,'lipids_flow')); 

% BMR: Personal Scaling -  Param 4
lb(end+1)     = log(0.7);                                                                    
ub(end+1)     = log(1.5); 

% Insulin resistance: Personal Scaling -  Param 5
lb(end+1)     = log(0.7);                                                                    
ub(end+1)     = log(1.5); 

% Insulin Clearance: Personal Scaling-  Param 6
lb(end+1)     = log(0.7);                                                                    
ub(end+1)     = log(1.5); 

% Insulin Repesponse: Personal Scaling-  Param 7
lb(end+1)     = log(0.7);                                                                    
ub(end+1)     = log(1.5); 

end

func =@(params)EstimationData_costfunction(Silfvergren2022Weight_data,Silfvergren2022GlucoseBloodSensor_data,Silfvergren2022Ketones_data,Silfvergren2022_meals,Krssak2004Healthy_data,DallaMan2007_data,AmountParameterOpti,params,model,initialvalues, pNames,sNames,fid);

%% Optimization

[EstimationData_params, minCostPS] = particleswarm(func, length(lb), lb, ub, options2);
EstimationData_params = exp(EstimationData_params);
save(['EstimationData_params' datestr(now, 'yymmdd-HHMMSS')],'EstimationData_params');

end

