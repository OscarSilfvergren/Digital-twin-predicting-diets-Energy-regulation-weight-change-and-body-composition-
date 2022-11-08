%% Load parameters

[sNames, ODEs, initialvalues] = IQMstates(model);
[pNames, param] = IQMparameters(model);

% Extract parameters
EstimationData_params = readtable('FullParameterEstimation_parameters.csv');
EstimationData_params = table2array(EstimationData_params);
[row column] = size(EstimationData_params);
EstimationData_params = sortrows(EstimationData_params,column);

%% Declare diets

Calories     = 2200; %kcal/day

%% Simulate Silfvergren

[row column] = size(EstimationData_params);

for j = 1:row
    
    [optimizedParamTemp,initialvaluesLow]    = Persondeclaration(180, 0, 1,55,10,initialvalues,sNames,EstimationData_params(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    [optimizedParamTemp,initialvaluesHigh]    = Persondeclaration(180, 0, 1,90,10,initialvalues,sNames,EstimationData_params(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
StandardDiet = optimizedParamTemp;

StandardDiet(ismember(pNames,"carb_flow"))    = Calories*0.5/(1440*4);
StandardDiet(ismember(pNames,"protein_flow")) = Calories*0.25/(1440*4);
StandardDiet(ismember(pNames,"lipids_flow"))  = Calories*0.25/(1440*9);


    SimLowSteadystate         = model(0:94080:94080,initialvaluesLow,StandardDiet);
    SimHighSteadystate        = model(0:94080:94080,initialvaluesHigh,StandardDiet);
    SimLowDiet                = model(0:1440:524160*5,SimLowSteadystate.statevalues(end,:),StandardDiet); % 24 weeks
    SimRegular                = model(0:1440:524160*5,SimHighSteadystate.statevalues(end,:),StandardDiet); % 24 weeks
    
    if j == 1
        simBestKeto = SimLowDiet;
        simBestRegular = SimRegular;
        
        maxKetoBW1 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'BW'));
        minKetoBW1 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'BW'));
        maxKetoBW2 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'BW'));
        minKetoBW2 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'BW'));
        
        maxRegularBW1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'BW'));
        minRegularBW1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'BW'));
        maxRegularBW2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'BW'));
        minRegularBW2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'BW'));
        
        maxKetoLeanBodyWeightKg1 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'LeanBodyWeightKg'));
        minKetoLeanBodyWeightKg1 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'LeanBodyWeightKg'));
        maxKetoLeanBodyWeightKg2 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'LeanBodyWeightKg'));
        minKetoLeanBodyWeightKg2 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'LeanBodyWeightKg'));
        
        
        maxRegularLeanBodyWeightKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'LeanBodyWeightKg'));
        minRegularLeanBodyWeightKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'LeanBodyWeightKg'));
        maxRegularLeanBodyWeightKg2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'LeanBodyWeightKg'));
        minRegularLeanBodyWeightKg2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'LeanBodyWeightKg'));
        
        maxKetoTGA_AdipocyteKg1 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'TGA_AdipocyteKg'));
        minKetoTGA_AdipocyteKg1 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'TGA_AdipocyteKg'));
        maxKetoTGA_AdipocyteKg2 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'TGA_AdipocyteKg'));
        minKetoTGA_AdipocyteKg2 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'TGA_AdipocyteKg'));
        
        maxRegularTGA_AdipocyteKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'TGA_AdipocyteKg'));
        minRegularTGA_AdipocyteKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'TGA_AdipocyteKg'));
        maxRegularTGA_AdipocyteKg2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'TGA_AdipocyteKg'));
        minRegularTGA_AdipocyteKg2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'TGA_AdipocyteKg'));
        
        maxKetoBMR1 = SimLowDiet.reactionvalues(:,ismember(SimLowDiet.reactions,'BMR'));
        minKetoBMR1 = SimLowDiet.reactionvalues(:,ismember(SimLowDiet.reactions,'BMR'));
        maxKetoBMR2 = SimLowDiet.reactionvalues(:,ismember(SimLowDiet.reactions,'BMR'));
        minKetoBMR2 = SimLowDiet.reactionvalues(:,ismember(SimLowDiet.reactions,'BMR'));
        
        maxRegularBMR1 = SimRegular.reactionvalues(:,ismember(SimRegular.reactions,'BMR'));
        minRegularBMR1 = SimRegular.reactionvalues(:,ismember(SimRegular.reactions,'BMR'));
        maxRegularBMR2 = SimRegular.reactionvalues(:,ismember(SimRegular.reactions,'BMR'));
        minRegularBMR2 = SimRegular.reactionvalues(:,ismember(SimRegular.reactions,'BMR'));
        
    end
    
    %Glucose
    maxKetoBW1 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'BW'));
    minKetoBW1 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'BW'));
    maxKetoBW2 = max(maxKetoBW2,maxKetoBW1);
    minKetoBW2 = min(minKetoBW2,minKetoBW1);
    
    %Glucose
    maxRegularBW1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'BW'));
    minRegularBW1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'BW'));
    maxRegularBW2 = max(maxRegularBW2,maxRegularBW1);
    minRegularBW2 = min(minRegularBW2,minRegularBW1);
    
    % LeanBodyWeightKg
    maxKetoLeanBodyWeightKg1 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'LeanBodyWeightKg'));
    minKetoLeanBodyWeightKg1 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'LeanBodyWeightKg'));
    maxKetoLeanBodyWeightKg2 = max(maxKetoLeanBodyWeightKg2,maxKetoLeanBodyWeightKg1);
    minKetoLeanBodyWeightKg2 = min(minKetoLeanBodyWeightKg2,minKetoLeanBodyWeightKg1);
    
    % LeanBodyWeightKg
    maxRegularLeanBodyWeightKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'LeanBodyWeightKg'));
    minRegularLeanBodyWeightKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'LeanBodyWeightKg'));
    maxRegularLeanBodyWeightKg2 = max(maxRegularLeanBodyWeightKg2,maxRegularLeanBodyWeightKg1);
    minRegularLeanBodyWeightKg2 = min(minRegularLeanBodyWeightKg2,minRegularLeanBodyWeightKg1);
    
    % TGA_AdipocyteKg
    maxKetoTGA_AdipocyteKg1 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'TGA_AdipocyteKg'));
    minKetoTGA_AdipocyteKg1 = SimLowDiet.variablevalues(:,ismember(SimLowDiet.variables,'TGA_AdipocyteKg'));
    maxKetoTGA_AdipocyteKg2 = max(maxKetoTGA_AdipocyteKg2,maxKetoTGA_AdipocyteKg1);
    minKetoTGA_AdipocyteKg2 = min(minKetoTGA_AdipocyteKg2,minKetoTGA_AdipocyteKg1);
    
    % TGA_AdipocyteKg
    maxRegularTGA_AdipocyteKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'TGA_AdipocyteKg'));
    minRegularTGA_AdipocyteKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'TGA_AdipocyteKg'));
    maxRegularTGA_AdipocyteKg2 = max(maxRegularTGA_AdipocyteKg2,maxRegularTGA_AdipocyteKg1);
    minRegularTGA_AdipocyteKg2 = min(minRegularTGA_AdipocyteKg2,minRegularTGA_AdipocyteKg1);
    
    % BMR
    maxKetoBMR1 = SimLowDiet.reactionvalues(:,ismember(SimLowDiet.reactions,'BMR'));
    minKetoBMR1 = SimLowDiet.reactionvalues(:,ismember(SimLowDiet.reactions,'BMR'));
    maxKetoBMR2 = max(maxKetoBMR2,maxKetoBMR1);
    minKetoBMR2 = min(minKetoBMR2,minKetoBMR1);
    
    % BMR
    maxRegularBMR1 = SimRegular.reactionvalues(:,ismember(SimRegular.reactions,'BMR'));
    minRegularBMR1 = SimRegular.reactionvalues(:,ismember(SimRegular.reactions,'BMR'));
    maxRegularBMR2 = max(maxRegularBMR2,maxRegularBMR1);
    minRegularBMR2 = min(minRegularBMR2,minRegularBMR1);
    
end

timeFill = [SimLowDiet.time, fliplr(SimLowDiet.time)];

%%
figure('Name', "Silfvergren et al. 2022 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Weight
hold on
set(gca,'ytick',[0,40,80,120,150],'xtick',[0,1,2,3,4,5],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Total body weight'; '(kg)'});
xlabel({'Time (years)'});
fill(timeFill/525948.766,[maxKetoBW2', fliplr(minKetoBW2')],'magenta','FaceAlpha',0.2,'EdgeAlpha',0);
fill(timeFill/525948.766,[maxRegularBW2', fliplr(minRegularBW2')],'k','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestKeto.time/525948.766 ,simBestKeto.variablevalues(:,ismember(simBestKeto.variables,'BW')),'m--','LineWidth',LineWidthValue-1)
plot(simBestRegular.time/525948.766,simBestRegular.variablevalues(:,ismember(simBestRegular.variables,'BW')),'k --','LineWidth',LineWidthValue-1)
%legend('Healthy BMI (20): Model uncertainty','Obese BMI (30): Model uncertainty','Healthy BMI (20): Best fit to estimation data','Obese BMI (30): Best fit to estimation data','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize+3,'NumColumns',2)
ylim([0 150])
xlim([0 5])
hold off
grid on

%%
figure('Name', "Silfvergren et al. 2022 - WeightMuscle", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Body composition
hold on
set(gca,'ytick',[0,1000,2000,3000],'xtick',[0,1,2,3,4,5],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'BMR';'(kcal/day)'});
xlabel({'Time (years)'});
fill(timeFill/525948.766,[minKetoBMR2', fliplr(maxKetoBMR2')],'magenta','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestKeto.time/525948.766,simBestKeto.reactionvalues(:,ismember(simBestKeto.reactions,'BMR')),'magenta --','Color', BlueColor,'LineWidth',LineWidthValue-1)
fill(timeFill/525948.766,[minRegularBMR2', fliplr(maxRegularBMR2')],'k','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestRegular.time/525948.766,simBestRegular.reactionvalues(:,ismember(simBestRegular.reactions,'BMR')),'k--','LineWidth',LineWidthValue-1)
ylim([0 3000])
xlim([0 5])
hold off
grid on

%%
WeightObese = simBestRegular.variablevalues(:,ismember(simBestRegular.variables,'BW'));
WeightHealthy   = simBestKeto.variablevalues(:,ismember(simBestKeto.variables,'BW'));

WeightDifHealthy = WeightHealthy(end) - WeightHealthy(1);
WeightDifObese   = WeightObese(end)   - WeightObese(1);

figure('Name', "Silfvergren et al. 2022 - WeightMuscle", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Weight Chang
X = categorical({'Healthy BMI (20)', 'Obese BMI (30)'});
X = reordercats(X,{'Healthy BMI (20)', 'Obese BMI (30)'});

subplot(1,2,1)
hold on
grid on
set(gca,'ytick',[-20,-15,-10,-5,0],'xtick',[0,1,2],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Weight loss'; 'from same diet (kg)'});
bar(X(1),WeightDifHealthy(1),'magenta');
bar(X(2),WeightDifObese(1),'k','FaceAlpha',0.6);
ylim([-20 0])
hold off

