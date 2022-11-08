%% Pre-calc

[sNames, ODEs, initialvalues] = IQMstates(model);
[pNames, param] = IQMparameters(model);

[param,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,param,pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
paramMeal               = Mealdeclaration(87,23,24,15,param,pNames);                           % CarbohydratesAmount, ProteinAmount, LipidsAmount, MealLength,


SimSteadystate       = model(0:282240:282240,initialvalues,param);    
param(ismember(pNames,"carb_flow"))    = 0;
param(ismember(pNames,"protein_flow")) = 0;
param(ismember(pNames,"lipids_flow"))  = 0;
SimSteadystate       = model(0:360:360,SimSteadystate.statevalues(end,:),param);  
SimRandom_UpToMeal   = model(1:1:450,SimSteadystate.statevalues(end,:),param);
SimRandom_Meal       = model(450:1:465,SimRandom_UpToMeal.statevalues(end,:),paramMeal);
SimRandom_AfterMeal  = model(465:1:1440,SimRandom_Meal.statevalues(end,:),param);

% MergeSimulations
SimRandom.time           = [SimRandom_UpToMeal.time,           SimRandom_Meal.time(2:end),           SimRandom_AfterMeal.time(2:end)];
SimRandom.statevalues    = [SimRandom_UpToMeal.statevalues; SimRandom_Meal.statevalues(2:end,:)   ; SimRandom_AfterMeal.statevalues(2:end,:)];
SimRandom.variablevalues = [SimRandom_UpToMeal.variablevalues; SimRandom_Meal.variablevalues(2:end,:); SimRandom_AfterMeal.variablevalues(2:end,:)];
SimRandom.reactionvalues = [SimRandom_UpToMeal.reactionvalues; SimRandom_Meal.reactionvalues(2:end,:); SimRandom_AfterMeal.reactionvalues(2:end,:)];
SimRandom.states         = [SimRandom_UpToMeal.states];
SimRandom.variables      = [SimRandom_UpToMeal.variables];
SimRandom.reactions      = [SimRandom_UpToMeal.reactions];
SimRandom.pNames         = [pNames];
time = SimRandom.time';
%% Stomach and intestines
figure('Name', "Stomach and intestines", 'units', 'normalized', 'outerposition', [0 0 1 1])
figureData.Position = [10 10 1300 800]; 

% Stomach_Carbohydrates
subplot(2,3,1)
hold on
title('Stomach Carbohydrates','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60,SimRandom.statevalues(:,ismember(SimRandom.states,'Stomach_Carbohydrates')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')

hold off

% Stomach_Protein
subplot(2,3,2)
hold on
title('Stomach Protein','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Stomach_Protein')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% Stomach_Protein
subplot(2,3,3)
hold on
title('Stomach Lipid','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Stomach_Lipid')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% Intestines_Carbohydrates
subplot(2,3,4)
hold on
title('Intestines Carbohydrates','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Intestines_Carbohydrates')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% Intestines_Protein
subplot(2,3,5)
hold on
title('Intestines Protein','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Intestines_Protein')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% Intestines_Lipid
subplot(2,3,6)
hold on
title('Intestines Lipid','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Intestines_Lipid')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

%% Plasma
figure('Name', "Plasma", 'units', 'normalized', 'outerposition', [0 0 1 1])
figureData.Position = [10 10 1300 800]; 

% Glucose_Plasma
subplot(2,3,1)
hold on
title('Glucose Plasma','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Glucose_Plasma')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% AminoAcids_Plasma
subplot(2,3,2)
hold on
title('AminoAcids Plasma','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'AminoAcids_Plasma')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% Lipids_Plasma
subplot(2,3,3)
hold on
title('Lipids Plasma','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Lipids_Plasma')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% Insulin_Plasma
subplot(2,3,4)
hold on
title('Insulin Plasma','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Insulin_Plasma')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% Ketones_Plasma
subplot(2,3,5)
hold on
title('Ketones Plasma','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Ketones_Plasma')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% EGP 
subplot(2,3,6)  
hold on
a = gca;
set(a,'xtick',[0,2,4,6,8,10,12],'ytick',[0,0.1,0.2,0.3],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'EGP (g/min)'});
xlabel('Time [h]')
line([0.0 0.2], [0 0],'Color','k','LineWidth',5);
errorbar(round(DallaMan2007_data.EGP_time,0)/60,DallaMan2007_data.EGP,DallaMan2007_data.EGPSEM, " .k",'MarkerSize',15,'LineWidth',4);
plot((time-450)/60, SimRandom.reactionvalues(:,ismember(SimRandom.reactions,'EGP')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlim([-1 8])
% ylim([0 0.3])
hold off

%% Liver
figure('Name', "Liver", 'units', 'normalized', 'outerposition', [0 0 1 1])
figureData.Position = [10 10 1300 800]; 

% Glucose_Liver
subplot(3,4,1)
hold on
title('Glucose Liver','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Glucose_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% Insulin_Liver
subplot(3,4,2)
hold on
title('Insulin Liver','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Insulin_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% Glycogen_Liver
subplot(3,4,3)
hold on
title('Glycogen Liver','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Glycogen_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% PEP_Liver
subplot(3,4,4)
hold on
title('PEP Liver','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'PEP_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% Pyruvate_Liver
subplot(3,4,5)
hold on
title('Pyruvate Liver','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Pyruvate_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% AcetylCoA_Liver
subplot(3,4,6)
hold on
title('AcetylCoA Liver','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'AcetylCoA_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% Citrate_Liver
subplot(3,4,7)
hold on
title('Citrate Liver','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Citrate_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% TCA_Liver
subplot(3,4,8)
hold on
title('TCA Liver','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'TCA_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% OAA_Liver
subplot(3,4,9)
hold on
title('OAA Liver','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'OAA_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% FatPool_Liver
subplot(3,4,10)
hold on
title('FatPool Liver','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'FatPool_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% TGA_Liver
subplot(3,4,11)
hold on
title('TGA Liver','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'TGA_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% Ketones_Liver
subplot(3,4,12)
hold on
title('Ketones Liver','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Ketones_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

%% Pancrease
figure('Name', "Pancrease", 'units', 'normalized', 'outerposition', [0 0 1 1])
figureData.Position = [10 10 1300 800]; 

% Insulin_Pancrease
hold on
title('Insulin Pancrease','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Insulin_Pancrease')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

%% Adipocytes
figure('Name', "Adipocytes", 'units', 'normalized', 'outerposition', [0 0 1 1])
figureData.Position = [10 10 1300 800]; 

% Glucose_Adipocytes
subplot(2,2,1)
hold on
title('Glucose Adipocytes','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Glucose_Adipocytes')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% MalonylCoA_Adipocytes
subplot(2,2,2)
hold on
title('MalonylCoA Adipocytes','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'MalonylCoA_Adipocytes')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% FatPool_Adipocytes
subplot(2,2,3)
hold on
title('FatPool Adipocytes','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'FatPool_Adipocytes')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

% TGA_Adipocyte
subplot(2,2,4)
hold on
title('TGA Adipocyte','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'TGA_Adipocyte')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[g]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

%% LeanBodyWeight
figure('Name', "Leanweight", 'units', 'normalized', 'outerposition', [0 0 1 1])
figureData.Position = [10 10 1300 800]; 

% LeanBodyWeight
hold on
title('LeanBodyWeight','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'LeanBodyWeight')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
 
hold off

%% Variables
figure('Name', "Cortesol", 'units', 'normalized', 'outerposition', [0 0 1 1])
figureData.Position = [10 10 1300 800]; 

% GlycogenConcentrationliver
subplot(3,2,1)
hold on
title('GlycogenConcentrationliver','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.variablevalues(:,ismember(SimRandom.variables,'GlycogenConcentrationliver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% InsulinConcentrationPlasma
subplot(3,2,2)
hold on
title('InsulinConcentrationPlasma','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.variablevalues(:,ismember(SimRandom.variables,'InsulinConcentrationPlasma')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% GlucoseConcentrationPlasma
subplot(3,2,3)
hold on
title('GlucoseConcentrationPlasma','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.variablevalues(:,ismember(SimRandom.variables,'GlucoseConcentrationPlasma')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% BW
subplot(3,2,4)
hold on
title('BW','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.variablevalues(:,ismember(SimRandom.variables,'BW')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% FFA
subplot(3,2,5)
hold on
title('FFAConcentrationPlasma','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.variablevalues(:,ismember(SimRandom.variables,'FFAConcentrationPlasma')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

%% Delays
figure('Name', "Cortesol", 'units', 'normalized', 'outerposition', [0 0 1 1])
figureData.Position = [10 10 1300 800]; 

subplot(2,2,1)
hold on
title('InsulinResponseFATP1','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'InsulinResponseFATP1_Delay')),'LineWidth',LineWidthValue)
plot((time-450)/60 ,SimRandom.statevalues(:,ismember(SimRandom.states,'Insulin_Plasma')),'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
legend( 'Delay','Insulin_Plasma')
hold off

subplot(2,2,2)
hold on
title('InsulinResponseScaled','FontSize', 15,'FontSmoothing','on')
plot((time-450)/60 ,SimRandom.reactionvalues(:,ismember(SimRandom.reactions,'InsulinResponseLiver')),'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
legend( 'InsulinResponseLiver')
hold off
