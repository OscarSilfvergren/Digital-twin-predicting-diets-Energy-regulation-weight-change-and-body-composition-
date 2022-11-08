%% Load parameters
% load('EstimationData_params2')
load('Silfvergren2022GlucoseBloodSensor_data');
load('Silfvergren2022Food_data');
load('Silfvergren2022Weight_data');
load('Silfvergren2022Ketones_data');
load('Silfvergren2022_meals');

[sNames, ODEs, initialvalues] = IQMstates(model);
[pNames, param] = IQMparameters(model);

% Clean data
Silfvergren2022_meals.MealStart           = round(Silfvergren2022_meals.MealStart/120,0)*120;
Silfvergren2022_meals.MealDuration(1:end) = [120];
Silfvergren2022_meals.MealEnd             = Silfvergren2022_meals.MealStart + 120;
Silfvergren2022Weight_data.Time                       = round(Silfvergren2022Weight_data.Time/120,0);
Silfvergren2022Ketones_data.Time_minutesKetones(1:66) = round(Silfvergren2022Ketones_data.Time_minutesKetones(1:66)/120,0);
Silfvergren2022GlucoseBloodSensor_data.glucoseTime    = round(Silfvergren2022GlucoseBloodSensor_data.glucoseTime/120,0);

% Extract parameters
EstimationData_params = readtable('FullParameterEstimation_parameters.csv');
EstimationData_params = table2array(EstimationData_params);
[row column] = size(EstimationData_params);
EstimationData_params = sortrows(EstimationData_params,column);

%% Cost to full estimation data
costEstimationData = EstimationData_costfunction(Silfvergren2022Weight_data,Silfvergren2022GlucoseBloodSensor_data,Silfvergren2022Ketones_data,Silfvergren2022_meals,Krssak2004Healthy_data,DallaMan2007_data,AmountParameterOpti,log(EstimationData_params(1,1:column-1)),model,initialvalues, pNames,sNames);
disp(' ')
if costEstimationData < chi2inv(0.95,800)
    disp('----- Fit to Estimation data is below statistical treshold -----')
end
fprintf('Best prediction of Estimation data: %.2f, Statistical Limit: %.2f (dgf = %i)', costEstimationData, chi2inv(0.95,800), 800) % Finddatapoints
disp(' ')

%% Simulate Silfvergren

[row column] = size(EstimationData_params);

for j = 1:row
    
    [Silfvergren2022_params,start]   = Persondeclaration(180, 1, 0,75,10,initialvalues,sNames,EstimationData_params(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
    % Simulate Silfvergren2022
    SimSteadystate       = model(0:282240:282240,start,Silfvergren2022_params);
    Silfvergren2022_params(ismember(pNames,"carb_flow"))    = 0;
    Silfvergren2022_params(ismember(pNames,"protein_flow")) = 0;
    Silfvergren2022_params(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate       = model(0:480:480,SimSteadystate.statevalues(end,:),Silfvergren2022_params);
    Silfvergren2022_paramsMeal            = Mealdeclaration(Silfvergren2022_meals.CarbohydratesAmount(1),Silfvergren2022_meals.ProteinAmount(1),Silfvergren2022_meals.LipidsAmount(1),Silfvergren2022_meals.MealDuration(1),Silfvergren2022_params,pNames);       % CarbohydratesAmount, ProteinAmount, LipidsAmount, MealLength,
    SimUpToMeal        = model(0:120:Silfvergren2022_meals.MealStart(1),SimSteadystate.statevalues(end,:),Silfvergren2022_params);
    SimMeal            = model(Silfvergren2022_meals.MealStart(1):120:Silfvergren2022_meals.MealEnd(1),SimUpToMeal.statevalues(end,:),Silfvergren2022_paramsMeal);
    Sim.time           = [SimUpToMeal.time          , SimMeal.time(2:end)          ];
    Sim.variablevalues = [SimUpToMeal.variablevalues; SimMeal.variablevalues(2:end,:)];
    Sim.reactionvalues = [SimUpToMeal.reactionvalues; SimMeal.reactionvalues(2:end,:)];
    Sim.statevalues    = [SimUpToMeal.statevalues   ; SimMeal.statevalues(2:end,:)   ];
    Sim.states         = [SimUpToMeal.states];
    Sim.variables      = [SimUpToMeal.variables];
    Sim.reactions      = [SimUpToMeal.reactions];
    [AmountOfMeals Start] = size(Silfvergren2022_meals.MealStart);
    for i=Start+1:AmountOfMeals-1
        Silfvergren2022_paramsMeal   = Mealdeclaration(Silfvergren2022_meals.CarbohydratesAmount(i),Silfvergren2022_meals.ProteinAmount(i),Silfvergren2022_meals.LipidsAmount(i),Silfvergren2022_meals.MealDuration(i),Silfvergren2022_params,pNames);       % CarbohydratesAmount, ProteinAmount, LipidsAmount, MealLength,
        
        SimUpToMeal = model(Silfvergren2022_meals.MealEnd(i-1):120:Silfvergren2022_meals.MealStart(i),Sim.statevalues(end,:),Silfvergren2022_params);
        SimMeal     = model(Silfvergren2022_meals.MealStart(i):120:Silfvergren2022_meals.MealEnd(i),SimUpToMeal.statevalues(end,:),Silfvergren2022_paramsMeal);
        
        Sim.time           = [Sim.time          , SimUpToMeal.time(2:end)            , SimMeal.time(2:end)];
        Sim.variablevalues = [Sim.variablevalues; SimUpToMeal.variablevalues(2:end,:); SimMeal.variablevalues(2:end,:)];
        Sim.reactionvalues = [Sim.reactionvalues; SimUpToMeal.reactionvalues(2:end,:); SimMeal.reactionvalues(2:end,:)];
        Sim.statevalues    = [Sim.statevalues   ; SimUpToMeal.statevalues(2:end,:)   ; SimMeal.statevalues(2:end,:)];
        Sim.states         = [Sim.states];
        Sim.variables      = [Sim.variables];
        Sim.reactions      = [Sim.reactions];
    end
    Silfvergren2022_paramsMeal = Mealdeclaration(Silfvergren2022_meals.CarbohydratesAmount(end),Silfvergren2022_meals.ProteinAmount(end),Silfvergren2022_meals.LipidsAmount(end),Silfvergren2022_meals.MealDuration(end),Silfvergren2022_params,pNames);       % CarbohydratesAmount, ProteinAmount, LipidsAmount, MealLength,
    SimUpToMeal          = model(Silfvergren2022_meals.MealEnd(end-1):120:Silfvergren2022_meals.MealStart(end),Sim.statevalues(end,:),Silfvergren2022_params);
    SimMeal              = model(Silfvergren2022_meals.MealStart(end):120:Silfvergren2022_meals.MealEnd(end)+480,SimUpToMeal.statevalues(end,:),Silfvergren2022_paramsMeal);
    simSilfvergren2022.time             = [Sim.time          , SimUpToMeal.time(2:end)            , SimMeal.time(2:end)];
    simSilfvergren2022.variablevalues   = [Sim.variablevalues; SimUpToMeal.variablevalues(2:end,:); SimMeal.variablevalues(2:end,:)];
    simSilfvergren2022.reactionvalues   = [Sim.reactionvalues; SimUpToMeal.reactionvalues(2:end,:); SimMeal.reactionvalues(2:end,:)];
    simSilfvergren2022.statevalues      = [Sim.statevalues   ; SimUpToMeal.statevalues(2:end,:)   ; SimMeal.statevalues(2:end,:)];
    simSilfvergren2022.states           = [Sim.states];
    simSilfvergren2022.variables        = [Sim.variables];
    simSilfvergren2022.reactions        = [Sim.reactions];
    
    if j == 1
        simBest = simSilfvergren2022;
        
        maxBW1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'BW'));
        minBW1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'BW'));
        maxBW2 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'BW'));
        minBW2 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'BW'));
        
        maxLeanBodyWeightKg1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'LeanBodyWeightKg'));
        minLeanBodyWeightKg1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'LeanBodyWeightKg'));
        maxLeanBodyWeightKg2 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'LeanBodyWeightKg'));
        minLeanBodyWeightKg2 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'LeanBodyWeightKg'));
        
        maxTGA_AdipocyteKg1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'TGA_AdipocyteKg'));
        minTGA_AdipocyteKg1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'TGA_AdipocyteKg'));
        maxTGA_AdipocyteKg2 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'TGA_AdipocyteKg'));
        minTGA_AdipocyteKg2 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'TGA_AdipocyteKg'));
        
        maxGlucoseConcentrationPlasma1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'GlucoseConcentrationPlasma'));
        minGlucoseConcentrationPlasma1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'GlucoseConcentrationPlasma'));
        maxGlucoseConcentrationPlasma2 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'GlucoseConcentrationPlasma'));
        minGlucoseConcentrationPlasma2 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'GlucoseConcentrationPlasma'));
        
        maxKetoneConcentrationPlasma1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'KetoneConcentrationPlasma'));
        minKetoneConcentrationPlasma1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'KetoneConcentrationPlasma'));
        maxKetoneConcentrationPlasma2 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'KetoneConcentrationPlasma'));
        minKetoneConcentrationPlasma2 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'KetoneConcentrationPlasma'));
    end
    
    
    %Glucose
    maxBW1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'BW'));
    minBW1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'BW'));
    maxBW2 = max(maxBW2,maxBW1);
    minBW2 = min(minBW2,minBW1);
    
    % Insulin
    maxLeanBodyWeightKg1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'LeanBodyWeightKg'));
    minLeanBodyWeightKg1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'LeanBodyWeightKg'));
    maxLeanBodyWeightKg2 = max(maxLeanBodyWeightKg2,maxLeanBodyWeightKg1);
    minLeanBodyWeightKg2 = min(minLeanBodyWeightKg2,minLeanBodyWeightKg1);
    
    
    % TGA_AdipocyteKg
    maxTGA_AdipocyteKg1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'TGA_AdipocyteKg'));
    minTGA_AdipocyteKg1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'TGA_AdipocyteKg'));
    maxTGA_AdipocyteKg2 = max(maxTGA_AdipocyteKg2,maxTGA_AdipocyteKg1);
    minTGA_AdipocyteKg2 = min(minTGA_AdipocyteKg2,minTGA_AdipocyteKg1);
    
    % GlucoseConcentrationPlasma
    maxGlucoseConcentrationPlasma1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'GlucoseConcentrationPlasma'));
    minGlucoseConcentrationPlasma1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'GlucoseConcentrationPlasma'));
    maxGlucoseConcentrationPlasma2 = max(maxGlucoseConcentrationPlasma2,maxGlucoseConcentrationPlasma1);
    minGlucoseConcentrationPlasma2 = min(minGlucoseConcentrationPlasma2,minGlucoseConcentrationPlasma1);
    
    
    % KetoneConcentrationPlasma
    maxKetoneConcentrationPlasma1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'KetoneConcentrationPlasma'));
    minKetoneConcentrationPlasma1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'KetoneConcentrationPlasma'));
    maxKetoneConcentrationPlasma2 = max(maxKetoneConcentrationPlasma2,maxKetoneConcentrationPlasma1);
    minKetoneConcentrationPlasma2 = min(minKetoneConcentrationPlasma2,minKetoneConcentrationPlasma1);
end

timeFill = [simSilfvergren2022.time, fliplr(simSilfvergren2022.time)];
%%

figure('Name', "Silfvergren et al. 2022 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Weight
hold on
grid on
set(gca,'ytick',[80,84,88,92],'xtick',[0,2,4,6,8],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Weight'; '(kg)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[maxBW2', fliplr(minBW2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBest.time/10080 ,simBest.variablevalues(:,ismember(simBest.variables,'BW')),'Color', BlueColor,'LineWidth',LineWidthValue-2)
plot(Silfvergren2022Weight_data.Time(1:113)*120/10080,Silfvergren2022Weight_data.Weight(1:113),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
line([Silfvergren2022_meals.MealStart Silfvergren2022_meals.MealStart+Silfvergren2022_meals.MealDuration+100]/10080, [80 80],'Color','k','LineWidth',LineWidthMeal);
legend('Model uncertainty','Best fit to all estimation data','Data','Meals','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize+3)
ylim([80 92])
xlim([0 8])
hold off

%%
figure('Name', "Silfvergren et al. 2022 - WeightMuscle", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Body composition
subplot(1,2,1)
hold on
grid on
set(gca,'ytick',[70,75,80],'xtick',[0,2,4,6,8],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Lean mass';'(kg)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[maxLeanBodyWeightKg2', fliplr(minLeanBodyWeightKg2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBest.time/10080,simBest.variablevalues(:,ismember(simBest.variables,'LeanBodyWeightKg')),'Color', BlueColor,'LineWidth',LineWidthValue-2)
plot(Silfvergren2022Weight_data.Time(1:113)*120/10080,Silfvergren2022Weight_data.Leanmass(1:113),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize,'LineWidth',LineWidthValue)
line([Silfvergren2022_meals.MealStart Silfvergren2022_meals.MealStart+Silfvergren2022_meals.MealDuration+100]/10080, [70 70],'Color','k','LineWidth',LineWidthMeal);
ylim([70 80])
xlim([0 8])
hold off

%Body composition
subplot(1,2,2)
hold on
grid on
set(gca,'ytick',[10,12,14,16],'xtick',[0,2,4,6,8],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Fat mass';'(kg)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[maxTGA_AdipocyteKg2', fliplr(minTGA_AdipocyteKg2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBest.time/10080,simBest.variablevalues(:,ismember(simBest.variables,'TGA_AdipocyteKg')),'Color', BlueColor,'LineWidth',LineWidthValue-2)
plot(Silfvergren2022Weight_data.Time(1:113)*120/10080,Silfvergren2022Weight_data.Fatmass(1:113),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize,'LineWidth',LineWidthValue)
line([Silfvergren2022_meals.MealStart Silfvergren2022_meals.MealStart+Silfvergren2022_meals.MealDuration+100]/10080, [10 10],'Color','k','LineWidth',LineWidthMeal);
ylim([10 16])
xlim([0 8])
hold off
grid on
%%
figure('Name', "Silfvergren et al. 2022", 'units', 'normalized', 'outerposition', [0 0 1 1])

timeFillBasalGlucose            = [simSilfvergren2022.time(Silfvergren2022GlucoseBloodSensor_data.glucoseTime), fliplr(simSilfvergren2022.time(Silfvergren2022GlucoseBloodSensor_data.glucoseTime))];
maxGlucose2  = maxGlucoseConcentrationPlasma2(Silfvergren2022GlucoseBloodSensor_data.glucoseTime);
minGlucose2  = minGlucoseConcentrationPlasma2(Silfvergren2022GlucoseBloodSensor_data.glucoseTime);
yFill                           = [maxGlucose2', fliplr(minGlucose2')];
xBest                           = simBest.time(Silfvergren2022GlucoseBloodSensor_data.glucoseTime);
yBest                           = simBest.variablevalues(:,ismember(simBest.variables,'GlucoseConcentrationPlasma'));

%Glucose
subplot(1,2,1)
hold on
grid on
set(gca,'ytick',[0,2,4,6,8,10],'xtick',[0,2,4,6,8],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Basal plasma glucose';'(mM)'});
xlabel({'Time (weeks)'});
fill(timeFillBasalGlucose/10080,yFill,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(xBest/10080,yBest(Silfvergren2022GlucoseBloodSensor_data.glucoseTime),'Color', BlueColor,'LineWidth',LineWidthValue-2)
plot(Silfvergren2022GlucoseBloodSensor_data.glucoseTime(1:56)*120/10080,Silfvergren2022GlucoseBloodSensor_data.bloodvalue(1:56),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize,'LineWidth',LineWidthValue)
line([Silfvergren2022_meals.MealStart Silfvergren2022_meals.MealStart+Silfvergren2022_meals.MealDuration+100]/10080, [0 0],'Color','k','LineWidth',LineWidthMeal);
ylim([0 10])
xlim([0 8])
hold off
grid on

%Ketones
subplot(1,2,2)
hold on
grid on
set(gca,'ytick',[0,0.5,1,1.5,2],'xtick',[0,2,4,6,8],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Plasma ketones'; '(mM)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[maxKetoneConcentrationPlasma2', fliplr(minKetoneConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBest.time/10080,simBest.variablevalues(:,ismember(simBest.variables,'KetoneConcentrationPlasma')),'Color', BlueColor,'LineWidth',LineWidthValue-2)
plot(Silfvergren2022Ketones_data.Time_minutesKetones(1:66)*120/10080,Silfvergren2022Ketones_data.Ketones(1:66),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize,'LineWidth',LineWidthValue)
line([Silfvergren2022_meals.MealStart Silfvergren2022_meals.MealStart+Silfvergren2022_meals.MealDuration+100]/10080, [0 0],'Color','k','LineWidth',LineWidthMeal);
ylim([0 2])
xlim([0 8])
hold off
grid on

%%
figure('Name', "Liver", 'units', 'normalized', 'outerposition', [0 0 1 1])
figureData.Position = [10 10 1300 800]; 

% Glucose_Liver
subplot(3,4,1)
hold on
grid on
title('Glucose Liver','FontSize', 15,'FontSmoothing','on')
plot(simBest.time/10080,simBest.statevalues(:,ismember(simBest.states,'Glucose_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% Insulin_Liver
subplot(3,4,2)
hold on
grid on
title('Insulin Liver','FontSize', 15,'FontSmoothing','on')
plot(simBest.time/10080,simBest.statevalues(:,ismember(simBest.states,'Insulin_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% Glycogen_Liver
subplot(3,4,3)
hold on
title('Glycogen Liver','FontSize', 15,'FontSmoothing','on')
plot(simBest.time/10080,simBest.statevalues(:,ismember(simBest.states,'Glycogen_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% PEP_Liver
subplot(3,4,4)
hold on
title('PEP Liver','FontSize', 15,'FontSmoothing','on')
plot(simBest.time/10080,simBest.statevalues(:,ismember(simBest.states,'PEP_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% Pyruvate_Liver
subplot(3,4,5)
hold on
title('Pyruvate Liver','FontSize', 15,'FontSmoothing','on')
plot(simBest.time/10080,simBest.statevalues(:,ismember(simBest.states,'Pyruvate_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% AcetylCoA_Liver
subplot(3,4,6)
hold on
title('AcetylCoA Liver','FontSize', 15,'FontSmoothing','on')
plot(simBest.time/10080,simBest.statevalues(:,ismember(simBest.states,'AcetylCoA_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% Citrate_Liver
subplot(3,4,7)
hold on
title('Citrate Liver','FontSize', 15,'FontSmoothing','on')
plot(simBest.time/10080,simBest.statevalues(:,ismember(simBest.states,'Citrate_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% TCA_Liver
subplot(3,4,8)
hold on
title('TCA Liver','FontSize', 15,'FontSmoothing','on')
plot(simBest.time/10080,simBest.statevalues(:,ismember(simBest.states,'TCA_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% OAA_Liver
subplot(3,4,9)
hold on
title('OAA Liver','FontSize', 15,'FontSmoothing','on')
plot(simBest.time/10080,simBest.statevalues(:,ismember(simBest.states,'OAA_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% FatPool_Liver
subplot(3,4,10)
hold on
title('FatPool Liver','FontSize', 15,'FontSmoothing','on')
plot(simBest.time/10080,simBest.statevalues(:,ismember(simBest.states,'FatPool_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% TGA_Liver
subplot(3,4,11)
hold on
title('TGA Liver','FontSize', 15,'FontSmoothing','on')
plot(simBest.time/10080,simBest.statevalues(:,ismember(simBest.states,'TGA_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off

% Ketones_Liver
subplot(3,4,12)
hold on
title('Ketones Liver','FontSize', 15,'FontSmoothing','on')
plot(simBest.time/10080,simBest.statevalues(:,ismember(simBest.states,'Ketones_Liver')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlabel("Time[h]",'FontSize', 15,'FontSmoothing','on');
y_label = ylabel("[mg/kg/min]");
set(y_label,'FontSize', 15,'FontSmoothing','on')
hold off
%%

%% Krssak2004 - Study 2

[row column] = size(EstimationData_params);

for j = 1:row
    
    % Study 2
Study2_SteadystateMeal                   = [EstimationData_params(j,AmountParameterOpti+1) , EstimationData_params(j,AmountParameterOpti+2),  EstimationData_params(j,AmountParameterOpti+3)];
Study2_BMRPersonalScaling                = EstimationData_params(j,AmountParameterOpti+4);
Study2_InsulinResistancePersonalScaling  = EstimationData_params(j,AmountParameterOpti+5);
Study2_InsulinClearancePersonalScaling   = EstimationData_params(j,AmountParameterOpti+6);
Study2_InsulinResponsePersonalScaling    = EstimationData_params(j,AmountParameterOpti+7);

    [Krssak2004_params,start]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,EstimationData_params(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
% Declare personal scaling
%Meals
Krssak2004_params(ismember(pNames,"carb_flow"))    = Study2_SteadystateMeal(1);
Krssak2004_params(ismember(pNames,"protein_flow")) = Study2_SteadystateMeal(2);
Krssak2004_params(ismember(pNames,"lipids_flow"))  = Study2_SteadystateMeal(3);
% BMR
Krssak2004_params(ismember(pNames,"BMRDynamicBW"))    = Krssak2004_params(ismember(pNames,'BMRDynamicBW'))  * Study2_BMRPersonalScaling;
Krssak2004_params(ismember(pNames,"BMRDynamicGly"))   = Krssak2004_params(ismember(pNames,'BMRDynamicGly')) * Study2_BMRPersonalScaling;
Krssak2004_params(ismember(pNames,"BMRbasal"))        = Krssak2004_params(ismember(pNames,'BMRbasal'))      * Study2_BMRPersonalScaling;
% InsulinResistancePersonalScaling
Krssak2004_params(ismember(pNames,"Vf_id"))                     = Krssak2004_params(ismember(pNames,'Vf_id'))                   * Study2_InsulinResistancePersonalScaling;
Krssak2004_params(ismember(pNames,"Vm_id"))                     = Krssak2004_params(ismember(pNames,'Vm_id'))                   * Study2_InsulinResistancePersonalScaling;
Krssak2004_params(ismember(pNames,"InsulinResponseFATP1_K"))   = Krssak2004_params(ismember(pNames,'InsulinResponseFATP1_K')) * Study2_InsulinResistancePersonalScaling;
Krssak2004_params(ismember(pNames,"InsulinResponseLiverK"))     = Krssak2004_params(ismember(pNames,'InsulinResponseLiverK'))   * Study2_InsulinResistancePersonalScaling;
% InsulinClearancePersonalScaling
Krssak2004_params(ismember(pNames,"InsulinDegradation_PlasmaK"))    = Krssak2004_params(ismember(pNames,'InsulinDegradation_PlasmaK'))    * Study2_InsulinClearancePersonalScaling;
Krssak2004_params(ismember(pNames,"HE_positiveK"))                  = Krssak2004_params(ismember(pNames,'HE_positiveK'))                  * Study2_InsulinClearancePersonalScaling;
% InsulinClearancePersonalScaling
Krssak2004_params(ismember(pNames,"InsulinStabilisationMagnitude"))    = Krssak2004_params(ismember(pNames,'InsulinStabilisationMagnitude'))    * Study2_InsulinResponsePersonalScaling;
Krssak2004_params(ismember(pNames,"InsulinProductionK"))               = Krssak2004_params(ismember(pNames,'InsulinProductionK'))               * Study2_InsulinResponsePersonalScaling;

Krssak2004_paramsMeal               = Mealdeclaration(87,23,24,15,Krssak2004_params,pNames);                           % CarbohydratesAmount, ProteinAmount, LipidsAmount, MealLength,

    % Simultate
    SimSteadystate       = model(0:282240:282240,start,Krssak2004_params);
    Krssak2004_params(ismember(pNames,"carb_flow"))    = 0;
    Krssak2004_params(ismember(pNames,"protein_flow")) = 0;
    Krssak2004_params(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate       = model(0:360:360,SimSteadystate.statevalues(end,:),Krssak2004_params);
    SimKrssakToMeal      = model(1:1:360,SimSteadystate.statevalues(end,:),Krssak2004_params);
    SimKrssak_meal       = model(360:1:375,SimKrssakToMeal.statevalues(end,:),Krssak2004_paramsMeal);
    SimKrssakAfterMeal   = model(375:1:1260,SimKrssak_meal.statevalues(end,:),Krssak2004_params);
    SimKrssakAfterday    = model(0:2725:2725,SimKrssakAfterMeal.statevalues(end,:),Krssak2004_params);
    simKrssak2004Healthy.time           = [SimKrssakToMeal.time,           SimKrssak_meal.time(2:end),           SimKrssakAfterMeal.time(2:end)]-360;
    simKrssak2004Healthy.statevalues    = [SimKrssakToMeal.statevalues; SimKrssak_meal.statevalues(2:end,:)   ; SimKrssakAfterMeal.statevalues(2:end,:)];
    simKrssak2004Healthy.variablevalues = [SimKrssakToMeal.variablevalues; SimKrssak_meal.variablevalues(2:end,:); SimKrssakAfterMeal.variablevalues(2:end,:)];
    simKrssak2004Healthy.reactionvalues = [SimKrssakToMeal.reactionvalues; SimKrssak_meal.reactionvalues(2:end,:); SimKrssakAfterMeal.reactionvalues(2:end,:)];
    simKrssak2004Healthy.states         = [SimKrssak_meal.states];
    simKrssak2004Healthy.variables      = [SimKrssak_meal.variables];
    simKrssak2004Healthy.reactions      = [SimKrssak_meal.reactions];
    simKrssak2004Healthy.pNames         = [pNames];
    
    if j == 1
        simBest = simKrssak2004Healthy;
        
        %GlycogenConcentrationliver
        maxGlycogenConcentrationliver1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'GlycogenConcentrationliver'));
        maxGlycogenConcentrationliver2 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver2 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'GlycogenConcentrationliver'));
        
        % EGP
        maxEGP1 = simKrssak2004Healthy.reactionvalues(:,ismember(simKrssak2004Healthy.reactions,'EGP'));
        minEGP1 = simKrssak2004Healthy.reactionvalues(:,ismember(simKrssak2004Healthy.reactions,'EGP'));
        maxEGP2 = simKrssak2004Healthy.reactionvalues(:,ismember(simKrssak2004Healthy.reactions,'EGP'));
        minEGP2 = simKrssak2004Healthy.reactionvalues(:,ismember(simKrssak2004Healthy.reactions,'EGP'));
        
        % GlucoseConcentrationPlasma
        maxGlucoseConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'GlucoseConcentrationPlasma'));
        minGlucoseConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'GlucoseConcentrationPlasma'));
        maxGlucoseConcentrationPlasma2 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'GlucoseConcentrationPlasma'));
        minGlucoseConcentrationPlasma2 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'GlucoseConcentrationPlasma'));
        
        % InsulinConcentrationPlasma
        maxInsulinConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'InsulinConcentrationPlasma'));
        minInsulinConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'InsulinConcentrationPlasma'));
        maxInsulinConcentrationPlasma2 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'InsulinConcentrationPlasma'));
        minInsulinConcentrationPlasma2 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'InsulinConcentrationPlasma'));
        
        % KetoneConcentrationPlasma
        maxKetoneConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'KetoneConcentrationPlasma'));
        minKetoneConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'KetoneConcentrationPlasma'));
        maxKetoneConcentrationPlasma2 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'KetoneConcentrationPlasma'));
        minKetoneConcentrationPlasma2 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'KetoneConcentrationPlasma'));
        
        % FFAConcentrationPlasma
        maxFFAConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'FFAConcentrationPlasma'));
        minFFAConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'FFAConcentrationPlasma'));
        maxFFAConcentrationPlasma2 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'FFAConcentrationPlasma'));
        minFFAConcentrationPlasma2 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'FFAConcentrationPlasma'));
    end
    
    
    %GlycogenConcentrationliver
    maxGlycogenConcentrationliver1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'GlycogenConcentrationliver'));
    minGlycogenConcentrationliver1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'GlycogenConcentrationliver'));
    maxGlycogenConcentrationliver2 = max(maxGlycogenConcentrationliver2,maxGlycogenConcentrationliver1);
    minGlycogenConcentrationliver2 = min(minGlycogenConcentrationliver2,minGlycogenConcentrationliver1);
    
    % EGP
    maxEGP1 = simKrssak2004Healthy.reactionvalues(:,ismember(simKrssak2004Healthy.reactions,'EGP'));
    minEGP1 = simKrssak2004Healthy.reactionvalues(:,ismember(simKrssak2004Healthy.reactions,'EGP'));
    maxEGP2 = max(maxEGP2,maxEGP1);
    minEGP2 = min(minEGP2,minEGP1);
    
    % GlucoseConcentrationPlasma
    maxGlucoseConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'GlucoseConcentrationPlasma'));
    minGlucoseConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'GlucoseConcentrationPlasma'));
    maxGlucoseConcentrationPlasma2 = max(maxGlucoseConcentrationPlasma2,maxGlucoseConcentrationPlasma1);
    minGlucoseConcentrationPlasma2 = min(minGlucoseConcentrationPlasma2,minGlucoseConcentrationPlasma1);
    
    % InsulinConcentrationPlasma
    maxInsulinConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'InsulinConcentrationPlasma'));
    minInsulinConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'InsulinConcentrationPlasma'));
    maxInsulinConcentrationPlasma2 = max(maxInsulinConcentrationPlasma2,maxInsulinConcentrationPlasma1);
    minInsulinConcentrationPlasma2 = min(minInsulinConcentrationPlasma2,minInsulinConcentrationPlasma1);
    
    % KetoneConcentrationPlasma
    maxKetoneConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'KetoneConcentrationPlasma'));
    minKetoneConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'KetoneConcentrationPlasma'));
    maxKetoneConcentrationPlasma2 = max(maxKetoneConcentrationPlasma2,maxKetoneConcentrationPlasma1);
    minKetoneConcentrationPlasma2 = min(minKetoneConcentrationPlasma2,minKetoneConcentrationPlasma1);
    
    % FFAConcentrationPlasma
    maxFFAConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'FFAConcentrationPlasma'));
    minFFAConcentrationPlasma1 = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'FFAConcentrationPlasma'));
    maxFFAConcentrationPlasma2 = max(maxFFAConcentrationPlasma2,maxFFAConcentrationPlasma1);
    minFFAConcentrationPlasma2 = min(minFFAConcentrationPlasma2,minFFAConcentrationPlasma1);
end

timeFill = [simKrssak2004Healthy.time, fliplr(simKrssak2004Healthy.time)];

%%

figure('Name', "Krssak et al. 2004 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Glycogen
hold on
grid on
set(gca,'xtick',[0,5,10,15],'ytick',[150,200,250,300,350,400],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Hepatic Glycogen'; '(mmol/L)'});
xlabel('Time [h]')
fill(timeFill/60,[maxGlycogenConcentrationliver2', fliplr(minGlycogenConcentrationliver2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot((simBest.time)/60, simBest.variablevalues(:,ismember(simBest.variables,'GlycogenConcentrationliver')),'Color', BlueColor,'LineWidth',LineWidthValue-2)
errorbar(Krssak2004Healthy_data.timeGlycogen(1:22)/60,Krssak2004Healthy_data.glycogen(1:22),Krssak2004Healthy_data.glycogenSEM(1:22), ". ",'Color', BlackColor,'MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
line([0.0 0.2], [180 180],'Color','k','LineWidth',LineWidthMeal);
legend('Model uncertainty','Best fit to all estimation data','Data','Mixed Meal','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize+3)
ylim([180 400])
xlim([-5 16])
hold off

%%
figure('Name', "Krssak et al. 2004 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

%EGP
hold on
grid on
set(gca,'xtick',[0,5,10,15],'ytick',[0,0.05,0.1,0.15],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'EGP'; '(g/min)'});
xlabel('Time [h]')
fill(timeFill/60,[maxEGP2', fliplr(minEGP2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
%ylim([150 400])
errorbar(Krssak2004Healthy_data.timeEGP(1:20)/60,Krssak2004Healthy_data.EGP(1:20),Krssak2004Healthy_data.EGPSEM(1:20), ".",'Color', BlackColor,'MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
plot((simBest.time)/60, simBest.reactionvalues(:,ismember(simBest.reactions,'EGP')),'Color', BlueColor,'LineWidth',LineWidthValue-2)
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
xlim([-5 16])
ylim([0 0.15])
hold off
%%
figure('Name', "Krssak et al. 2004 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Glucose
hold on
grid on
set(gca,'xtick',[0,5,10,15],'ytick',[2,4,6,8,10,12],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Plasma glucose'; '(mmol/L)'});
xlabel('Time [h]')
fill(timeFill/60,[maxGlucoseConcentrationPlasma2', fliplr(minGlucoseConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
errorbar(Krssak2004Healthy_data.timeGlucose(1:25)/60,Krssak2004Healthy_data.glucose(1:25),Krssak2004Healthy_data.glucoseSEM(1:25), " . ",'Color', BlackColor,'MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
plot((simBest.time)/60, simBest.variablevalues(:,ismember(simBest.variables,'GlucoseConcentrationPlasma')),'Color', BlueColor,'LineWidth',LineWidthValue-2)
line([0.0 0.2], [2 2],'Color','k','LineWidth',LineWidthMeal);
xlim([-5 16])
ylim([2 12])
hold off

%%
figure('Name', "Krssak et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Insulin
hold on
grid on
set(gca,'xtick',[0,5,10,15],'ytick',[0,200,400,600,800],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Insulin Blood' ;'(pM)'});
xlabel('Time [h]')
fill(timeFill/60,[maxInsulinConcentrationPlasma2', fliplr(minInsulinConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
errorbar(Krssak2004Healthy_data.timeInsulin(1:10)/60,Krssak2004Healthy_data.insulin(1:10),Krssak2004Healthy_data.insulinSEM(1:10), " . ",'Color', BlackColor,'MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
plot((simBest.time)/60, simBest.variablevalues(:,ismember(simBest.variables,'InsulinConcentrationPlasma')),'Color', BlueColor,'LineWidth',LineWidthValue-2)
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
ylim([0 800])
xlim([-5 16])
hold off
grid on

%%
figure('Name', "Krssak et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% FFA
hold on
grid on
set(gca,'xtick',[0,5,10,15],'ytick',[0,0.2,0.4,0.6,0.8,1],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'FFA'; '(mM)'});
xlabel('Time [h]')
fill(timeFill/60,[maxFFAConcentrationPlasma2', fliplr(minFFAConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
errorbar(Krssak2004Healthy_data.timeFFA(1:10)/60,Krssak2004Healthy_data.FFA(1:10),Krssak2004Healthy_data.FFASEM(1:10), ". ",'Color', BlackColor,'MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
plot((simBest.time)/60, simBest.variablevalues(:,ismember(simBest.variables,'FFAConcentrationPlasma')),'Color', BlueColor,'LineWidth',LineWidthValue-2)
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
ylim([0 1])
xlim([-5 16])
hold off

%%

%% DallaMan2007 - study 3

[row column] = size(EstimationData_params);

for j = 1:row
    
    % Study 3
Study3_SteadystateMeal                   = [EstimationData_params(j,AmountParameterOpti+8) , EstimationData_params(j,AmountParameterOpti+9),  EstimationData_params(j,AmountParameterOpti+10)];
Study3_BMRPersonalScaling                = EstimationData_params(j,AmountParameterOpti+11);
Study3_InsulinResistancePersonalScaling  = EstimationData_params(j,AmountParameterOpti+12);
Study3_InsulinClearancePersonalScaling   = EstimationData_params(j,AmountParameterOpti+13);
Study3_InsulinResponsePersonalScaling    = EstimationData_params(j,AmountParameterOpti+14);

    
    [DallaMan2007_params,start]   = Persondeclaration(170, 1, 0,70,8,initialvalues,sNames,EstimationData_params(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
    % Meals
    DallaMan2007_params(ismember(pNames,"carb_flow"))    = Study3_SteadystateMeal(1);
    DallaMan2007_params(ismember(pNames,"protein_flow")) = Study3_SteadystateMeal(2);
    DallaMan2007_params(ismember(pNames,"lipids_flow"))  = Study3_SteadystateMeal(3);
    % BMR
    DallaMan2007_params(ismember(pNames,"BMRDynamicBW"))    = DallaMan2007_params(ismember(pNames,'BMRDynamicBW'))  * Study3_BMRPersonalScaling;
    DallaMan2007_params(ismember(pNames,"BMRDynamicGly"))   = DallaMan2007_params(ismember(pNames,'BMRDynamicGly')) * Study3_BMRPersonalScaling;
    DallaMan2007_params(ismember(pNames,"BMRbasal"))        = DallaMan2007_params(ismember(pNames,'BMRbasal'))      * Study3_BMRPersonalScaling;
    % InsulinResistancePersonalScaling
    DallaMan2007_params(ismember(pNames,"Vf_id"))                          = DallaMan2007_params(ismember(pNames,'Vf_id'))                   * Study3_InsulinResistancePersonalScaling;
    DallaMan2007_params(ismember(pNames,"Vm_id"))                          = DallaMan2007_params(ismember(pNames,'Vm_id'))                   * Study3_InsulinResistancePersonalScaling;
    DallaMan2007_params(ismember(pNames,"InsulinResponseFATP1_K"))        = DallaMan2007_params(ismember(pNames,'InsulinResponseFATP1_K')) * Study3_InsulinResistancePersonalScaling;
    DallaMan2007_params(ismember(pNames,"InsulinResponseLiverK"))          = DallaMan2007_params(ismember(pNames,'InsulinResponseLiverK'))   * Study3_InsulinResistancePersonalScaling;
    % InsulinClearancePersonalScaling
    DallaMan2007_params(ismember(pNames,"InsulinDegradation_PlasmaK"))    = DallaMan2007_params(ismember(pNames,'InsulinDegradation_PlasmaK'))    * Study3_InsulinClearancePersonalScaling;
    DallaMan2007_params(ismember(pNames,"HE_positiveK"))                  = DallaMan2007_params(ismember(pNames,'HE_positiveK'))                  * Study3_InsulinClearancePersonalScaling;
    % InsulinClearancePersonalScaling
    DallaMan2007_params(ismember(pNames,"InsulinStabilisationMagnitude")) = DallaMan2007_params(ismember(pNames,'InsulinStabilisationMagnitude'))    * Study3_InsulinResponsePersonalScaling;
    DallaMan2007_params(ismember(pNames,"InsulinProductionK"))            = DallaMan2007_params(ismember(pNames,'InsulinProductionK'))               * Study3_InsulinResponsePersonalScaling;
    DallaMan2007_paramsMeal               = Mealdeclaration(78,0,0,1,DallaMan2007_params,pNames);                           % CarbohydratesAmount, ProteinAmount, LipidsAmount, MealLength,
    
    % Simultate DallaMan2007
    SimSteadystate       = model(0:282240:282240,start,DallaMan2007_params);
    DallaMan2007_params(ismember(pNames,"carb_flow"))    = 0;
    DallaMan2007_params(ismember(pNames,"protein_flow")) = 0;
    DallaMan2007_params(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate       = model(0:360:360,SimSteadystate.statevalues(end,:),DallaMan2007_params);
    SimDallaMan2007ToMeal      = model(1:1:360,SimSteadystate.statevalues(end,:),DallaMan2007_params);
    SimDallaMan2007_meal       = model(360:1:361,SimDallaMan2007ToMeal.statevalues(end,:),DallaMan2007_paramsMeal);
    SimDallaMan2007AfterMeal   = model(361:1:1440,SimDallaMan2007_meal.statevalues(end,:),DallaMan2007_params);
    SimDallaMan2007.time           = [SimDallaMan2007ToMeal.time,           SimDallaMan2007_meal.time(2:end),           SimDallaMan2007AfterMeal.time(2:end)]-360;
    SimDallaMan2007.statevalues    = [SimDallaMan2007ToMeal.statevalues; SimDallaMan2007_meal.statevalues(2:end,:)   ; SimDallaMan2007AfterMeal.statevalues(2:end,:)];
    SimDallaMan2007.variablevalues = [SimDallaMan2007ToMeal.variablevalues; SimDallaMan2007_meal.variablevalues(2:end,:); SimDallaMan2007AfterMeal.variablevalues(2:end,:)];
    SimDallaMan2007.reactionvalues = [SimDallaMan2007ToMeal.reactionvalues; SimDallaMan2007_meal.reactionvalues(2:end,:); SimDallaMan2007AfterMeal.reactionvalues(2:end,:)];
    SimDallaMan2007.states         = [SimDallaMan2007_meal.states];
    SimDallaMan2007.variables      = [SimDallaMan2007_meal.variables];
    SimDallaMan2007.reactions      = [SimDallaMan2007_meal.reactions];
    SimDallaMan2007.pNames         = [pNames];
    
    if j == 1
        simBest = SimDallaMan2007;
        
        % GlucoseConcentrationPlasma
        maxGlucoseConcentrationPlasma1 = SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'GlucoseConcentrationPlasma'));
        minGlucoseConcentrationPlasma1 = SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'GlucoseConcentrationPlasma'));
        maxGlucoseConcentrationPlasma2 = SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'GlucoseConcentrationPlasma'));
        minGlucoseConcentrationPlasma2 = SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'GlucoseConcentrationPlasma'));
        
        % InsulinConcentrationPlasma
        maxInsulinConcentrationPlasma1 = SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'InsulinConcentrationPlasma'));
        minInsulinConcentrationPlasma1 = SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'InsulinConcentrationPlasma'));
        maxInsulinConcentrationPlasma2 = SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'InsulinConcentrationPlasma'));
        minInsulinConcentrationPlasma2 = SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'InsulinConcentrationPlasma'));
        
        % EGP
        maxEGP1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'EGP'));
        minEGP1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'EGP'));
        maxEGP2 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'EGP'));
        minEGP2 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'EGP'));
        
        % Ra
        maxRa1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'Carbohydrates_Intestines_out'));
        minRa1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'Carbohydrates_Intestines_out'));
        maxRa2 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'Carbohydrates_Intestines_out'));
        minRa2 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'Carbohydrates_Intestines_out'));
        
                % U
        maxU1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'U'));
        minU1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'U'));
        maxU2 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'U'));
        minU2 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'U'));
        
        % S
        maxS1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'S'));
        minS1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'S'));
        maxS2 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'S'));
        minS2 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'S'));
    end
    
    % GlucoseConcentrationPlasma
    maxGlucoseConcentrationPlasma1 = SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'GlucoseConcentrationPlasma'));
    minGlucoseConcentrationPlasma1 = SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'GlucoseConcentrationPlasma'));
    maxGlucoseConcentrationPlasma2 = max(maxGlucoseConcentrationPlasma2,maxGlucoseConcentrationPlasma1);
    minGlucoseConcentrationPlasma2 = min(minGlucoseConcentrationPlasma2,minGlucoseConcentrationPlasma1);
    
    % InsulinConcentrationPlasma
    maxInsulinConcentrationPlasma1 = SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'InsulinConcentrationPlasma'));
    minInsulinConcentrationPlasma1 = SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'InsulinConcentrationPlasma'));
    maxInsulinConcentrationPlasma2 = max(maxInsulinConcentrationPlasma2,maxInsulinConcentrationPlasma1);
    minInsulinConcentrationPlasma2 = min(minInsulinConcentrationPlasma2,minInsulinConcentrationPlasma1);
    
    % EGP
    maxEGP1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'EGP'));
    minEGP1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'EGP'));
    maxEGP2 = max(maxEGP2,maxEGP1);
    minEGP2 = min(minEGP2,minEGP1);
    
        % Ra
    maxRa1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'Carbohydrates_Intestines_out'));
    minRa1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'Carbohydrates_Intestines_out'));
    maxRa2 = max(maxRa2,maxRa1);
    minRa2 = min(minRa2,minRa1);
    
        % U
    maxU1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'U'));
    minU1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'U'));
    maxU2 = max(maxU2,maxU1);
    minU2 = min(minU2,minU1);
    
        % S
    maxS1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'S'));
    minS1 = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'S'));
    maxS2 = max(maxS2,maxS1);
    minS2 = min(minS2,minS1);

end

timeFill = [SimDallaMan2007.time, fliplr(SimDallaMan2007.time)];


% Clean data
Glucose_time = round(DallaMan2007_data.Glucose_time,0);
Insulin_time = round(DallaMan2007_data.Insulin_time,0);
EGP_time = round(DallaMan2007_data.EGP_time,0);
Ra_time = round(DallaMan2007_data.Ra_time,0);
U_time = round(DallaMan2007_data.U_time,0);
S_time = round(DallaMan2007_data.S_time,0);

%%
figure('Name', "Krssak et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Glucose
hold on
grid on
a = gca;
set(a,'xtick',[0,2,4,6,8,10,12],'ytick',[4,6,8,10,12],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Plasma Glucose'; '(mM)'});
xlabel('Time [h]')
fill(timeFill/60,[maxGlucoseConcentrationPlasma2', fliplr(minGlucoseConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
errorbar(Glucose_time/60,DallaMan2007_data.Glucose,DallaMan2007_data.GlucoseSEM, " .k ",'MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
plot((SimDallaMan2007.time)/60, SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'GlucoseConcentrationPlasma')),'Color', BlueColor,'LineWidth',LineWidthValue)
line([0.0 0.2], [4 4],'Color','k','LineWidth',LineWidthMeal);
legend('Model uncertainty','Best fit to all estimation data','Data','OGTT','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize+3)
xlim([-1 8])
ylim([4 12])
hold off

%%
figure('Name', "Krssak et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Insulin
hold on
grid on
a = gca;
set(a,'xtick',[0,2,4,6,8,10,12],'ytick',[0,200,400,600],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Plasma Insulin' ;'(pM)'});
xlabel('Time [h]')
fill(timeFill/60,[maxInsulinConcentrationPlasma2', fliplr(minInsulinConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
errorbar(Insulin_time/60,DallaMan2007_data.Insulin,DallaMan2007_data.InsulinSEM, " .k",'MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
plot((SimDallaMan2007.time)/60, SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'InsulinConcentrationPlasma')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlim([-1 8])
ylim([0 600])
hold off
%%
figure('Name', "Krssak et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% EGP
hold on
grid on
a = gca;
set(a,'xtick',[0,2,4,6,8,10,12],'ytick',[0,0.1,0.2,0.3],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'EGP' ; '(g/min)'});
xlabel('Time [h]')
fill(timeFill/60,[maxEGP2', fliplr(minEGP2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
errorbar(EGP_time/60,DallaMan2007_data.EGP,DallaMan2007_data.EGPSEM, " .k",'MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
plot((SimDallaMan2007.time)/60, SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'EGP')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlim([-1 8])
ylim([0 0.3])
hold off

%%
figure('Name', "Krssak et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Ra
hold on
grid on
a = gca;
set(a,'xtick',[0,2,4,6,8,10,12],'ytick',[0,0.4,0.8,1.2],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Ra' ; '(g/min)'});
xlabel('Time [h]')
fill(timeFill/60,[maxRa2', fliplr(minRa2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
errorbar(Ra_time/60,DallaMan2007_data.Ra,DallaMan2007_data.RaSEM, " .k",'MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
plot((SimDallaMan2007.time)/60, SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'Carbohydrates_Intestines_out')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlim([-1 8])
ylim([0 1.2])
hold off

%%
figure('Name', "Krssak et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% U
hold on
grid on
a = gca;
set(a,'xtick',[0,2,4,6,8,10,12],'ytick',[0,0.3,0.6,0.9],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Glucose U' ; '(g/min)'});
xlabel('Time [h]')
fill(timeFill/60,[maxU2', fliplr(minU2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
errorbar(U_time/60,DallaMan2007_data.U,DallaMan2007_data.USEM, " .k",'MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
plot((SimDallaMan2007.time)/60, SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'U')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlim([-1 8])
ylim([0 0.9])
hold off

%%
figure('Name', "Krssak et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% S
hold on
grid on
a = gca;
set(a,'xtick',[0,2,4,6,8,10,12],'ytick',[0,400,800,1200],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Insulin secration' ; '(pM/min)'});
xlabel('Time [h]')
fill(timeFill/60,[maxS2', fliplr(minS2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
line([0.0 0.2], [0 0],'Color','k','LineWidth',LineWidthMeal);
errorbar(S_time/60,DallaMan2007_data.S,DallaMan2007_data.SSEM, "  .k",'MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
plot((SimDallaMan2007.time)/60, SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'S')),'Color', BlueColor,'LineWidth',LineWidthValue)
xlim([-1 8])
ylim([0 1200])
hold off