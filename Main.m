LoadAll

%% Load all settings, parameters and model

modelName = 'ModelMetabolismIndividualmeals';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

ConfigureParameterInformation

%% -------------------------------- Estimation -------------------------------- %%

%% Full parameter estimation

%    fid = fopen('FullParameterEstimation_parameters3.csv','w+');
%    [minCostPS] = FullParameterEstimationOptimization(i,Silfvergren2022Weight_data,Silfvergren2022GlucoseBloodSensor_data,Silfvergren2022Ketones_data,Silfvergren2022_meals,Krssak2004Healthy_data,DallaMan2007_data, model,ParametersInformation,AmountParameterOpti,options,fid);
%    fclose(fid);

% PlotFitToEstimationData % Plot script takes about 35 seconds.

%% -------------------------------- Validation --------------------------------- %%

%%  Silfvergren2022

% PlotPredictionSilfvergren2022 % Plot script takes about 65 seconds.

%% Rothman1991

% [minCostPS] = Rothman1991_Optimization(i,Rothman1991_data, model,ParametersInformation,AmountParameterOpti,optionsParalell);
% [minCostPS] = Rothman1991PS_Optimization(i,Rothman1991_data, model,ParametersInformation,AmountParameterOpti,optionsParalell);

%  PlotRothman1991Validation % Plot script takes about 1 second.
%  PlotRothman1991PSValidation % Plot script takes about 6 seconds.

%% Bray 2012

% [minCostPS] = Bray2012Optimization(i,Bray2012_data, model,ParametersInformation,AmountParameterOpti,optionsParalell);
% PlotBray2012Validation % Plot script takes about 1 second.
 
%% Taylor1996_data

% [minCostPS] = Taylor1996Optimization(i,Taylor1996_data, model,ParametersInformation,AmountParameterOpti,optionsParalell);
% PlotTaylor1996_dataValidation % Plot script takes about 1 second.

%% Firth1986_data

%  [minCostPS] = Firth1986Optimization(i,Firth1986_data, model,ParametersInformation,AmountParameterOpti,optionsParalell);
%  PlotFirth1986Validation % Plot script takes about 1 second.
%% -------------------------------- Prediction plots --------------------------------- %%

% WeightChangeDiets    % Plot script takes about 1 second.
% BMRDifferentCalories % Plot script takes about 1 second.
% BMRDifferentPersons  % Plot script takes about 1 second.

%% Kill all paralell workers and clusters

poolobj = gcp('nocreate');
myCluster = parcluster('local');
delete(myCluster.Jobs)
delete(poolobj)