%% Load parameters
load('EstimationData_params')
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

for j = 1:row
SilfvergrenTemp  = EstimationData_params(j,1:AmountParameterOpti);

costSilfvergrenFull = Silfvergren2022_costfunction(Silfvergren2022Weight_data,Silfvergren2022GlucoseBloodSensor_data,Silfvergren2022Ketones_data,Silfvergren2022_meals,AmountParameterOpti,log(SilfvergrenTemp),model,initialvalues, pNames,sNames);

EstimationData_params(j,column) = costSilfvergrenFull;
end

EstimationData_params = sortrows(EstimationData_params,column);

disp(' ')
if EstimationData_params(1,column) < chi2inv(0.95,517)
    disp('----- Fit to Estimation data is below statistical treshold -----')
end
fprintf('Best prediction of Estimation data: %.2f, Statistical Limit: %.2f (dgf = %i)', costSilfvergrenFull, chi2inv(0.95,517), 517) % Finddatapoints
disp(' ')


%% Simulate Silfvergren

[row column] = size(EstimationData_params);

for j = 1:row
    
    [Silfvergren2022_params,initialvalues]   = Persondeclaration(180, 1, 0,75,10,initialvalues,sNames,EstimationData_params(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
    % Simulate Silfvergren2022
    SimSteadystate       = model(0:282240:282240,initialvalues,Silfvergren2022_params);
    Silfvergren2022_params(ismember(pNames,"carb_flow"))    = 0;
    Silfvergren2022_params(ismember(pNames,"protein_flow")) = 0;
    Silfvergren2022_params(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate       = model(0:360:360,SimSteadystate.statevalues(end,:),Silfvergren2022_params);
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
    
    
    % Glycogen
    maxTGA_AdipocyteKg1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'TGA_AdipocyteKg'));
    minTGA_AdipocyteKg1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'TGA_AdipocyteKg'));
    maxTGA_AdipocyteKg2 = max(maxTGA_AdipocyteKg2,maxTGA_AdipocyteKg1);
    minTGA_AdipocyteKg2 = min(minTGA_AdipocyteKg2,minTGA_AdipocyteKg1);
    
    % EGP
    maxGlucoseConcentrationPlasma1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'GlucoseConcentrationPlasma'));
    minGlucoseConcentrationPlasma1 = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'GlucoseConcentrationPlasma'));
    maxGlucoseConcentrationPlasma2 = max(maxGlucoseConcentrationPlasma2,maxGlucoseConcentrationPlasma1);
    minGlucoseConcentrationPlasma2 = min(minGlucoseConcentrationPlasma2,minGlucoseConcentrationPlasma1);
    
    
    % EGP
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
set(gca,'ytick',[80,85,90,95],'xtick',[0,2,4,6,8,10,12,14,16],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Weight'; '(kg)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[maxBW2', fliplr(minBW2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBest.time/10080 ,simBest.variablevalues(:,ismember(simBest.variables,'BW')),'--','Color', BlueColor,'LineWidth',LineWidthValue-1)
plot(Silfvergren2022Weight_data.Time(1:56)*120/10080,Silfvergren2022Weight_data.Weight(1:56),'Color', BlackColor,'LineStyle', 'n', 'Marker','x','MarkerSize',MarkerSize-15,'LineWidth',LineWidthValue-2)
plot(Silfvergren2022Weight_data.Time(57:113)*120/10080,Silfvergren2022Weight_data.Weight(57:113),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize,'LineWidth',LineWidthValue)
line([Silfvergren2022_meals.MealStart Silfvergren2022_meals.MealStart+Silfvergren2022_meals.MealDuration+100]/10080, [78 78],'Color','k','LineWidth',LineWidthMeal);
legend('Model uncertainty','Best prediction','Estimation data','Validation data','Meals','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize+3)
ylim([78 95])
xlim([0 16])
hold off
grid on

%%
figure('Name', "Silfvergren et al. 2022 - WeightMuscle", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Body composition
hold on
set(gca,'ytick',[70,75,80],'xtick',[0,2,4,6,8,10,12,14,16],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Lean mass';'(kg)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[maxLeanBodyWeightKg2', fliplr(minLeanBodyWeightKg2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBest.time/10080,simBest.variablevalues(:,ismember(simBest.variables,'LeanBodyWeightKg')),'--','Color', BlueColor,'LineWidth',LineWidthValue-1)
plot(Silfvergren2022Weight_data.Time(1:56)*120/10080,Silfvergren2022Weight_data.Leanmass(1:56),'Color', BlackColor,'LineStyle', 'n', 'Marker','x','MarkerSize',MarkerSize-15,'LineWidth',LineWidthValue-2)
plot(Silfvergren2022Weight_data.Time(57:end)*120/10080,Silfvergren2022Weight_data.Leanmass(57:end),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize,'LineWidth',LineWidthValue)
line([Silfvergren2022_meals.MealStart Silfvergren2022_meals.MealStart+Silfvergren2022_meals.MealDuration+100]/10080, [70 70],'Color','k','LineWidth',LineWidthMeal);
ylim([70 80])
xlim([0 16])
hold off
grid on

%%
figure('Name', "Silfvergren et al. 2022 - WeightFat", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Body composition
hold on
set(gca,'ytick',[10,12,14,16],'xtick',[0,2,4,6,8,10,12,14,16],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Fat mass';'(kg)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[maxTGA_AdipocyteKg2', fliplr(minTGA_AdipocyteKg2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBest.time/10080,simBest.variablevalues(:,ismember(simBest.variables,'TGA_AdipocyteKg')),'--','Color', BlueColor,'LineWidth',LineWidthValue-1)
plot(Silfvergren2022Weight_data.Time(1:56)*120/10080,Silfvergren2022Weight_data.Fatmass(1:56),'Color', BlackColor,'LineStyle', 'n', 'Marker','x','MarkerSize',MarkerSize-15,'LineWidth',LineWidthValue-2)
plot(Silfvergren2022Weight_data.Time(57:end)*120/10080,Silfvergren2022Weight_data.Fatmass(57:end),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize,'LineWidth',LineWidthValue)
line([Silfvergren2022_meals.MealStart Silfvergren2022_meals.MealStart+Silfvergren2022_meals.MealDuration+100]/10080, [6 6],'Color','k','LineWidth',LineWidthMeal);
ylim([6 16])
xlim([0 16])
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
%subplot(1,1,1)
hold on
set(gca,'ytick',[0,2,4,6,8,10],'xtick',[0,2,4,6,8,10,12,14,16],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Basal plasma glucose';'(mM)'});
xlabel({'Time (weeks)'});
fill(timeFillBasalGlucose/10080,yFill,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(xBest/10080,yBest(Silfvergren2022GlucoseBloodSensor_data.glucoseTime),'--','Color', BlueColor,'LineWidth',LineWidthValue-1)
plot(Silfvergren2022GlucoseBloodSensor_data.glucoseTime(1:56)*120/10080,Silfvergren2022GlucoseBloodSensor_data.bloodvalue(1:56),'Color', BlackColor,'LineStyle', 'n', 'Marker','x','MarkerSize',MarkerSize-15,'LineWidth',LineWidthValue-2)
plot(Silfvergren2022GlucoseBloodSensor_data.glucoseTime(56:end)*120/10080,Silfvergren2022GlucoseBloodSensor_data.bloodvalue(56:end),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize,'LineWidth',LineWidthValue)
line([Silfvergren2022_meals.MealStart Silfvergren2022_meals.MealStart+Silfvergren2022_meals.MealDuration+100]/10080, [0 0],'Color','k','LineWidth',LineWidthMeal);
ylim([0 10])
xlim([0 16])
hold off
grid on
%%
 figure('Name', "Silfvergren et al. 2022", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Ketones
%subplot(1,2,2)
hold on
set(gca,'ytick',[0,0.5,1,1.5,2],'xtick',[0,2,4,6,8,10,12,14,16],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Plasma ketones'; '(mM)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[maxKetoneConcentrationPlasma2', fliplr(minKetoneConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBest.time/10080,simBest.variablevalues(:,ismember(simBest.variables,'KetoneConcentrationPlasma')),'--','Color', BlueColor,'LineWidth',LineWidthValue-1)
plot(Silfvergren2022Ketones_data.Time_minutesKetones(1:33)*120/10080,Silfvergren2022Ketones_data.Ketones(1:33),'Color', BlackColor,'LineStyle', 'n', 'Marker','x','MarkerSize',MarkerSize-15,'LineWidth',LineWidthValue-2)
plot(Silfvergren2022Ketones_data.Time_minutesKetones(33:end)*120/10080,Silfvergren2022Ketones_data.Ketones(33:end),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize+10,'LineWidth',LineWidthValue)
line([Silfvergren2022_meals.MealStart Silfvergren2022_meals.MealStart+Silfvergren2022_meals.MealDuration+100]/10080, [0 0],'Color','k','LineWidth',LineWidthMeal);
ylim([0 2])
xlim([0 16])
hold off
grid on