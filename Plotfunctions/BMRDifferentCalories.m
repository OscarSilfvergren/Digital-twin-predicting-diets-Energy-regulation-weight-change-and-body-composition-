%% Load parameters

[sNames, ODEs, initialvalues] = IQMstates(model);
[pNames, param] = IQMparameters(model);

% Extract parameters
EstimationData_params = readtable('FullParameterEstimation_parameters.csv');
EstimationData_params = table2array(EstimationData_params);
[row column] = size(EstimationData_params);
EstimationData_params = sortrows(EstimationData_params,column);

%% Declare diets

CaloriesHigh     = 2500; %kcal/day
CaloriesLow      = 1500; %kcal/day

%% Simulate Silfvergren

[row column] = size(EstimationData_params);

for j = 1:row
    
    [optimizedParamTemp,initialvalues]    = Persondeclaration(180, 0, 1,80,10,initialvalues,sNames,EstimationData_params(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
    LowDietParams = optimizedParamTemp;
    LowDietParams(ismember(pNames,"carb_flow"))    = CaloriesLow/(1440*4)*0.5;
    LowDietParams(ismember(pNames,"protein_flow")) = CaloriesLow/(1440*4)*0.25;
    LowDietParams(ismember(pNames,"lipids_flow"))  = CaloriesLow/(1440*9)*0.75;
    
    StandardDietParams  = optimizedParamTemp;
    StandardDietParams(ismember(pNames,"carb_flow"))    = CaloriesHigh/(1440*4)*0.5;
    StandardDietParams(ismember(pNames,"protein_flow")) = CaloriesHigh/(1440*4)*0.25;
    StandardDietParams(ismember(pNames,"lipids_flow"))  = CaloriesHigh/(1440*9)*0.25;
    
    SimSteadystate = model(0:20160:20160,initialvalues,optimizedParamTemp);
    SimLowDiet     = model(0:60:524160,SimSteadystate.statevalues(end,:),LowDietParams); % 24 weeks
    SimRegular     = model(0:60:524160,SimSteadystate.statevalues(end,:),StandardDietParams); % 24 weeks
    
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

OrangeColor  = [0.8500, 0.3250, 0.0980];
BlackColor   = [0 0 0]; 
GreenColor   = [0.4660, 0.6740, 0.1880];
YellowColor  = [0.9290, 0.6940, 0.1250];

%%
figure('Name', "Silfvergren et al. 2022 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Weight
hold on
set(gca,'ytick',[60,80,100,120],'xtick',[0,4,8,12,16,20,24],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Total body weight'; '(kg)'});
xlabel({'Time (months)'});
fill(timeFill/43829.0639,[maxKetoBW2', fliplr(minKetoBW2')],'g','FaceAlpha',0.3,'EdgeAlpha',0.1);
fill(timeFill/43829.0639,[maxRegularBW2', fliplr(minRegularBW2')],'c','FaceAlpha',0.4,'EdgeAlpha',0.1);
plot(simBestKeto.time/43829.0639 ,simBestKeto.variablevalues(:,ismember(simBestKeto.variables,'BW')),'g--','LineWidth',LineWidthValue-1)
plot(simBestRegular.time/43829.0639 ,simBestRegular.variablevalues(:,ismember(simBestRegular.variables,'BW')),'c--','LineWidth',LineWidthValue-1)
%legend('1500 kcal diet: Model uncertainty','2500 kcal diet: Model uncertainty','1500 diet: Best fit to estimation data','2500 kcal diet: Best fit to estimation data','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize+3,'NumColumns',2)
ylim([60 120])
xlim([0 12])
hold off
grid on

%%
figure('Name', "Silfvergren et al. 2022 - WeightMuscle", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Body composition
hold on
set(gca,'ytick',[0,1000,2000,3000],'xtick',[0,2,4,6,8,10,12],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'BMR';'(kcal/day)'});
xlabel({'Time (months)'});
fill(timeFill/43829.0639,[minKetoBMR2', fliplr(maxKetoBMR2')],'g','FaceAlpha',0.3,'EdgeAlpha',0.1);
plot(simBestKeto.time/43829.0639,simBestKeto.reactionvalues(:,ismember(simBestKeto.reactions,'BMR')),'g--','LineWidth',LineWidthValue-1)
fill(timeFill/43829.0639,[minRegularBMR2', fliplr(maxRegularBMR2')],'c','FaceAlpha',0.4,'EdgeAlpha',0.1);
plot(simBestRegular.time/43829.0639,simBestRegular.reactionvalues(:,ismember(simBestRegular.reactions,'BMR')),'c--','LineWidth',LineWidthValue-1)
ylim([0 3000])
xlim([0 12])
hold off
grid on

%%
BMRHealthy     = simBestRegular.reactionvalues(:,ismember(simBestRegular.reactions,'BMR'));
BMRObese   = simBestKeto.reactionvalues(:,ismember(simBestKeto.reactions,'BMR'));

BMRDifHealthy = BMRHealthy(1) - BMRHealthy(end);
BMRDifObese   = BMRObese(1) - BMRObese(end);

figure('Name', "Silfvergren et al. 2022 - WeightMuscle", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Weight Chang
X = categorical({'1500 kcal/day', '2500 kcal/day'});
X = reordercats(X,{'1500 kcal/day', '2500 kcal/day'});

subplot(1,2,1)
hold on
grid on
set(gca,'ytick',[-20,-10,0,10,20,30],'xtick',[0,1,2],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'BMR change'; 'from different diets (kcal/day)'});
bar(X(1),BMRDifHealthy(1),'g');
bar(X(2),BMRDifObese(1),'c');
ylim([-10 30])
hold off
