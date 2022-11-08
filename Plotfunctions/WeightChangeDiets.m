%% Load parameters

[sNames, ODEs, initialvalues] = IQMstates(model);
[pNames, param] = IQMparameters(model);

% Extract parameters
EstimationData_params = readtable('FullParameterEstimation_parameters.csv');
EstimationData_params = table2array(EstimationData_params);
[row column] = size(EstimationData_params);
EstimationData_params = sortrows(EstimationData_params,column);

%% Declare diets

Calories     = 2500; %kcal/day

Carb_gram    = Calories/(1440*4); %g/min
Protein_gram = Calories/(1440*4); %g/min
Lipids_gram  = Calories/(1440*9); %g/min


%% Simulate Silfvergren

[row column] = size(EstimationData_params);

for j = 1:row
    
    [optimizedParamTemp,initialvalues]    = Persondeclaration(180, 0, 1,80,10,initialvalues,sNames,EstimationData_params(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
    KetogenicDietParams = optimizedParamTemp;
    KetogenicDietParams(ismember(pNames,"carb_flow"))    = Carb_gram*0.1;
    KetogenicDietParams(ismember(pNames,"protein_flow")) = Protein_gram*0.1;
    KetogenicDietParams(ismember(pNames,"lipids_flow"))  = Lipids_gram*0.8;
    
    StandardDietParams  = optimizedParamTemp;
    StandardDietParams(ismember(pNames,"carb_flow"))    = Carb_gram*0.5;
    StandardDietParams(ismember(pNames,"protein_flow")) = Protein_gram*0.25;
    StandardDietParams(ismember(pNames,"lipids_flow"))  = Lipids_gram*0.25;
    
    SimSteadystate         = model(0:282240:282240,initialvalues,optimizedParamTemp);
    SimKeto     = model(0:60:241920,SimSteadystate.statevalues(end,:),KetogenicDietParams); % 24 weeks
    SimRegular  = model(0:60:241920,SimSteadystate.statevalues(end,:),StandardDietParams); % 24 weeks
    
    if j == 1
        simBestKeto = SimKeto;
        simBestRegular = SimRegular;
        
        maxKetoBW1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'BW'));
        minKetoBW1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'BW'));
        maxKetoBW2 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'BW'));
        minKetoBW2 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'BW'));
        
        maxKetoLeanBodyWeightKg1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'LeanBodyWeightKg'));
        minKetoLeanBodyWeightKg1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'LeanBodyWeightKg'));
        maxKetoLeanBodyWeightKg2 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'LeanBodyWeightKg'));
        minKetoLeanBodyWeightKg2 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'LeanBodyWeightKg'));
        
        maxKetoTGA_AdipocyteKg1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'TGA_AdipocyteKg'));
        minKetoTGA_AdipocyteKg1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'TGA_AdipocyteKg'));
        maxKetoTGA_AdipocyteKg2 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'TGA_AdipocyteKg'));
        minKetoTGA_AdipocyteKg2 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'TGA_AdipocyteKg'));
        
        maxKetoGlucoseConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'GlucoseConcentrationPlasma'));
        minKetoGlucoseConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'GlucoseConcentrationPlasma'));
        maxKetoGlucoseConcentrationPlasma2 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'GlucoseConcentrationPlasma'));
        minKetoGlucoseConcentrationPlasma2 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'GlucoseConcentrationPlasma'));
        
        maxKetoKetoneConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'KetoneConcentrationPlasma'));
        minKetoKetoneConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'KetoneConcentrationPlasma'));
        maxKetoKetoneConcentrationPlasma2 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'KetoneConcentrationPlasma'));
        minKetoKetoneConcentrationPlasma2 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'KetoneConcentrationPlasma'));
        
        maxRegularBW1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'BW'));
        minRegularBW1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'BW'));
        maxRegularBW2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'BW'));
        minRegularBW2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'BW'));
        
        maxRegularLeanBodyWeightKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'LeanBodyWeightKg'));
        minRegularLeanBodyWeightKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'LeanBodyWeightKg'));
        maxRegularLeanBodyWeightKg2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'LeanBodyWeightKg'));
        minRegularLeanBodyWeightKg2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'LeanBodyWeightKg'));
        
        maxRegularTGA_AdipocyteKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'TGA_AdipocyteKg'));
        minRegularTGA_AdipocyteKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'TGA_AdipocyteKg'));
        maxRegularTGA_AdipocyteKg2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'TGA_AdipocyteKg'));
        minRegularTGA_AdipocyteKg2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'TGA_AdipocyteKg'));
        
        maxRegularGlucoseConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'GlucoseConcentrationPlasma'));
        minRegularGlucoseConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'GlucoseConcentrationPlasma'));
        maxRegularGlucoseConcentrationPlasma2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'GlucoseConcentrationPlasma'));
        minRegularGlucoseConcentrationPlasma2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'GlucoseConcentrationPlasma'));
        
        maxRegularKetoneConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'KetoneConcentrationPlasma'));
        minRegularKetoneConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'KetoneConcentrationPlasma'));
        maxRegularKetoneConcentrationPlasma2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'KetoneConcentrationPlasma'));
        minRegularKetoneConcentrationPlasma2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'KetoneConcentrationPlasma'));
        
        maxKetoFFAConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'FFAConcentrationPlasma'));
        minKetoFFAConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'FFAConcentrationPlasma'));
        maxKetoFFAConcentrationPlasma2 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'FFAConcentrationPlasma'));
        minKetoFFAConcentrationPlasma2 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'FFAConcentrationPlasma'));
        
        maxRegularFFAConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'FFAConcentrationPlasma'));
        minRegularFFAConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'FFAConcentrationPlasma'));
        maxRegularFFAConcentrationPlasma2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'FFAConcentrationPlasma'));
        minRegularFFAConcentrationPlasma2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'FFAConcentrationPlasma'));
        
        maxKetoInsulinConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'InsulinConcentrationPlasma'));
        minKetoInsulinConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'InsulinConcentrationPlasma'));
        maxKetoInsulinConcentrationPlasma2 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'InsulinConcentrationPlasma'));
        minKetoInsulinConcentrationPlasma2 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'InsulinConcentrationPlasma'));
        
        maxRegularInsulinConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'InsulinConcentrationPlasma'));
        minRegularInsulinConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'InsulinConcentrationPlasma'));
        maxRegularInsulinConcentrationPlasma2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'InsulinConcentrationPlasma'));
        minRegularInsulinConcentrationPlasma2 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'InsulinConcentrationPlasma'));
        
        maxKetoOAA1 = SimKeto.statevalues(:,ismember(SimKeto.states,'OAA_Liver'));
        minKetoOAA1 = SimKeto.statevalues(:,ismember(SimKeto.states,'OAA_Liver'));
        maxKetoOAA2 = SimKeto.statevalues(:,ismember(SimKeto.states,'OAA_Liver'));
        minKetoOAA2 = SimKeto.statevalues(:,ismember(SimKeto.states,'OAA_Liver'));
        
        maxRegularOAA1 = SimRegular.statevalues(:,ismember(SimRegular.states,'OAA_Liver'));
        minRegularOAA1 = SimRegular.statevalues(:,ismember(SimRegular.states,'OAA_Liver'));
        maxRegularOAA2 = SimRegular.statevalues(:,ismember(SimRegular.states,'OAA_Liver'));
        minRegularOAA2 = SimRegular.statevalues(:,ismember(SimRegular.states,'OAA_Liver'));
        
    end
    
    
    %Glucose
    maxKetoBW1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'BW'));
    minKetoBW1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'BW'));
    maxKetoBW2 = max(maxKetoBW2,maxKetoBW1);
    minKetoBW2 = min(minKetoBW2,minKetoBW1);
    
    % Insulin
    maxKetoLeanBodyWeightKg1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'LeanBodyWeightKg'));
    minKetoLeanBodyWeightKg1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'LeanBodyWeightKg'));
    maxKetoLeanBodyWeightKg2 = max(maxKetoLeanBodyWeightKg2,maxKetoLeanBodyWeightKg1);
    minKetoLeanBodyWeightKg2 = min(minKetoLeanBodyWeightKg2,minKetoLeanBodyWeightKg1);
    
    
    % Glycogen
    maxKetoTGA_AdipocyteKg1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'TGA_AdipocyteKg'));
    minKetoTGA_AdipocyteKg1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'TGA_AdipocyteKg'));
    maxKetoTGA_AdipocyteKg2 = max(maxKetoTGA_AdipocyteKg2,maxKetoTGA_AdipocyteKg1);
    minKetoTGA_AdipocyteKg2 = min(minKetoTGA_AdipocyteKg2,minKetoTGA_AdipocyteKg1);
    
    % EGP
    maxKetoGlucoseConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'GlucoseConcentrationPlasma'));
    minKetoGlucoseConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'GlucoseConcentrationPlasma'));
    maxKetoGlucoseConcentrationPlasma2 = max(maxKetoGlucoseConcentrationPlasma2,maxKetoGlucoseConcentrationPlasma1);
    minKetoGlucoseConcentrationPlasma2 = min(minKetoGlucoseConcentrationPlasma2,minKetoGlucoseConcentrationPlasma1);
    
    
    % EGP
    maxKetoKetoneConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'KetoneConcentrationPlasma'));
    minKetoKetoneConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'KetoneConcentrationPlasma'));
    maxKetoKetoneConcentrationPlasma2 = max(maxKetoKetoneConcentrationPlasma2,maxKetoKetoneConcentrationPlasma1);
    minKetoKetoneConcentrationPlasma2 = min(minKetoKetoneConcentrationPlasma2,minKetoKetoneConcentrationPlasma1);
    
    %Glucose
    maxRegularBW1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'BW'));
    minRegularBW1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'BW'));
    maxRegularBW2 = max(maxRegularBW2,maxRegularBW1);
    minRegularBW2 = min(minRegularBW2,minRegularBW1);
    
    % Insulin
    maxRegularLeanBodyWeightKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'LeanBodyWeightKg'));
    minRegularLeanBodyWeightKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'LeanBodyWeightKg'));
    maxRegularLeanBodyWeightKg2 = max(maxRegularLeanBodyWeightKg2,maxRegularLeanBodyWeightKg1);
    minRegularLeanBodyWeightKg2 = min(minRegularLeanBodyWeightKg2,minRegularLeanBodyWeightKg1);
    
    
    % Glycogen
    maxRegularTGA_AdipocyteKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'TGA_AdipocyteKg'));
    minRegularTGA_AdipocyteKg1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'TGA_AdipocyteKg'));
    maxRegularTGA_AdipocyteKg2 = max(maxRegularTGA_AdipocyteKg2,maxRegularTGA_AdipocyteKg1);
    minRegularTGA_AdipocyteKg2 = min(minRegularTGA_AdipocyteKg2,minRegularTGA_AdipocyteKg1);
    
    % EGP
    maxRegularGlucoseConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'GlucoseConcentrationPlasma'));
    minRegularGlucoseConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'GlucoseConcentrationPlasma'));
    maxRegularGlucoseConcentrationPlasma2 = max(maxRegularGlucoseConcentrationPlasma2,maxRegularGlucoseConcentrationPlasma1);
    minRegularGlucoseConcentrationPlasma2 = min(minRegularGlucoseConcentrationPlasma2,minRegularGlucoseConcentrationPlasma1);
    
    
    % EGP
    maxRegularKetoneConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'KetoneConcentrationPlasma'));
    minRegularKetoneConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'KetoneConcentrationPlasma'));
    maxRegularKetoneConcentrationPlasma2 = max(maxRegularKetoneConcentrationPlasma2,maxRegularKetoneConcentrationPlasma1);
    minRegularKetoneConcentrationPlasma2 = min(minRegularKetoneConcentrationPlasma2,minRegularKetoneConcentrationPlasma1);
    
    % FFAConcentrationPlasma
    maxKetoFFAConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'FFAConcentrationPlasma'));
    minKetoFFAConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'FFAConcentrationPlasma'));
    maxKetoFFAConcentrationPlasma2 = max(maxKetoFFAConcentrationPlasma2,maxKetoFFAConcentrationPlasma1);
    minKetoFFAConcentrationPlasma2 = min(minKetoFFAConcentrationPlasma2,minKetoFFAConcentrationPlasma1);
    
    
    % FFAConcentrationPlasma
    maxRegularFFAConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'FFAConcentrationPlasma'));
    minRegularFFAConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'FFAConcentrationPlasma'));
    maxRegularFFAConcentrationPlasma2 = max(maxRegularFFAConcentrationPlasma2,maxRegularFFAConcentrationPlasma1);
    minRegularFFAConcentrationPlasma2 = min(minRegularFFAConcentrationPlasma2,minRegularFFAConcentrationPlasma1);
    
    % FFAConcentrationPlasma
    maxKetoInsulinConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'InsulinConcentrationPlasma'));
    minKetoInsulinConcentrationPlasma1 = SimKeto.variablevalues(:,ismember(SimKeto.variables,'InsulinConcentrationPlasma'));
    maxKetoInsulinConcentrationPlasma2 = max(maxKetoInsulinConcentrationPlasma2,maxKetoInsulinConcentrationPlasma1);
    minKetoInsulinConcentrationPlasma2 = min(minKetoInsulinConcentrationPlasma2,minKetoInsulinConcentrationPlasma1);
    
    
    % InsulinConcentrationPlasma
    maxRegularInsulinConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'InsulinConcentrationPlasma'));
    minRegularInsulinConcentrationPlasma1 = SimRegular.variablevalues(:,ismember(SimRegular.variables,'InsulinConcentrationPlasma'));
    maxRegularInsulinConcentrationPlasma2 = max(maxRegularInsulinConcentrationPlasma2,maxRegularInsulinConcentrationPlasma1);
    minRegularInsulinConcentrationPlasma2 = min(minRegularInsulinConcentrationPlasma2,minRegularInsulinConcentrationPlasma1);
    
    
        % FFAConcentrationPlasma
    maxKetoOAA1 = SimKeto.statevalues(:,ismember(SimKeto.states,'OAA_Liver'));
    minKetoOAA1 = SimKeto.statevalues(:,ismember(SimKeto.states,'OAA_Liver'));
    maxKetoOAA2 = max(maxKetoOAA2,maxKetoOAA1);
    minKetoOAA2 = min(minKetoOAA2,minKetoOAA1);
    
    % OAA
    maxRegularOAA1 = SimRegular.statevalues(:,ismember(SimRegular.states,'OAA_Liver'));
    minRegularOAA1 = SimRegular.statevalues(:,ismember(SimRegular.states,'OAA_Liver'));
    maxRegularOAA2 = max(maxRegularOAA2,maxRegularOAA1);
    minRegularOAA2 = min(minRegularOAA2,minRegularOAA1);
    
    
end

timeFill = [SimKeto.time, fliplr(SimKeto.time)];

%%
figure('Name', "Silfvergren et al. 2022 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Weight
hold on
set(gca,'ytick',[70,80,90,100,110,120],'xtick',[0,4,8,12,16,20,24],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Total body weight'; '(kg)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[maxKetoBW2', fliplr(minKetoBW2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
fill(timeFill/10080,[maxRegularBW2', fliplr(minRegularBW2')],'r','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestKeto.time/10080 ,simBestKeto.variablevalues(:,ismember(simBestKeto.variables,'BW')),'--','Color', BlueColor,'LineWidth',LineWidthValue-1)
plot(simBestRegular.time/10080 ,simBestRegular.variablevalues(:,ismember(simBestRegular.variables,'BW')),'r--','LineWidth',LineWidthValue-1)
legend('Ketogenic diet: Model uncertainty','Standardized diet: Model uncertainty','Ketogenic diet: Best fit to estimation data','Standardized diet: Best fit to estimation data','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize+3,'NumColumns',2)
ylim([70 120])
xlim([0 24])
hold off
grid on

%%
figure('Name', "Silfvergren et al. 2022 - WeightMuscle", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Body composition
hold on
set(gca,'ytick',[70,80,90,100],'xtick',[0,4,8,12,16,20,24],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Lean mass';'(kg)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[maxKetoLeanBodyWeightKg2', fliplr(minKetoLeanBodyWeightKg2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestKeto.time/10080,simBestKeto.variablevalues(:,ismember(simBestKeto.variables,'LeanBodyWeightKg')),'--','Color', BlueColor,'LineWidth',LineWidthValue-1)
fill(timeFill/10080,[maxRegularLeanBodyWeightKg2', fliplr(minRegularLeanBodyWeightKg2')],'r','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestRegular.time/10080,simBestRegular.variablevalues(:,ismember(simBestRegular.variables,'LeanBodyWeightKg')),'r--','LineWidth',LineWidthValue-1)
ylim([70 100])
xlim([0 24])
hold off
grid on

%%
figure('Name', "Silfvergren et al. 2022", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Ketones
%subplot(1,2,2)
hold on
set(gca,'ytick',[0,0.2,0.4,0.6,0.8,1],'xtick',[0,4,8,12,16,20,24],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Plasma ketones'; '(mM)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[maxKetoKetoneConcentrationPlasma2', fliplr(minKetoKetoneConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestKeto.time/10080,simBestKeto.variablevalues(:,ismember(simBestKeto.variables,'KetoneConcentrationPlasma')),'--','Color', BlueColor,'LineWidth',LineWidthValue-1)
fill(timeFill/10080,[maxRegularKetoneConcentrationPlasma2', fliplr(minRegularKetoneConcentrationPlasma2')],'r','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestRegular.time/10080,simBestRegular.variablevalues(:,ismember(simBestRegular.variables,'KetoneConcentrationPlasma')),'r--','LineWidth',LineWidthValue-1)
ylim([0 1])
xlim([0 24])
hold off
grid on

%%
figure('Name', "Silfvergren et al. 2022", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Ketones
%subplot(1,2,2)
hold on
set(gca,'ytick',[0,2,4,6,8,10],'xtick',[0,4,8,12,16,20,24],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Plasma glucose'; '(mM)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[maxKetoGlucoseConcentrationPlasma2', fliplr(minKetoGlucoseConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestKeto.time/10080,simBestKeto.variablevalues(:,ismember(simBestKeto.variables,'GlucoseConcentrationPlasma')),'--','Color', BlueColor,'LineWidth',LineWidthValue-1)
fill(timeFill/10080,[maxRegularGlucoseConcentrationPlasma2', fliplr(minRegularGlucoseConcentrationPlasma2')],'r','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestRegular.time/10080,simBestRegular.variablevalues(:,ismember(simBestRegular.variables,'GlucoseConcentrationPlasma')),'r--','LineWidth',LineWidthValue-1)
ylim([0 10])
xlim([0 24])
hold off
grid on

%%
figure('Name', "Silfvergren et al. 2022", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Ketones
%subplot(1,2,2)
hold on
set(gca,'ytick',[0,1,2,3],'xtick',[0,4,8,12,16,20,24],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'FFA'; '(mM)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[maxKetoFFAConcentrationPlasma2', fliplr(minKetoFFAConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestKeto.time/10080,simBestKeto.variablevalues(:,ismember(simBestKeto.variables,'FFAConcentrationPlasma')),'--','Color', BlueColor,'LineWidth',LineWidthValue-1)
fill(timeFill/10080,[minRegularFFAConcentrationPlasma2', fliplr(maxRegularFFAConcentrationPlasma2')],'r','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestRegular.time/10080,simBestRegular.variablevalues(:,ismember(simBestRegular.variables,'FFAConcentrationPlasma')),'r--','LineWidth',LineWidthValue-1)
ylim([0 3])
xlim([0 24])
hold off
grid on

%%
figure('Name', "Silfvergren et al. 2022", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Ketones
%subplot(1,2,2)
hold on
set(gca,'ytick',[0,100,200,300,400,500],'xtick',[0,4,8,12,16,20,24],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Plasma insulin'; '(pM)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[minKetoInsulinConcentrationPlasma2', fliplr(maxKetoInsulinConcentrationPlasma2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestKeto.time/10080,simBestKeto.variablevalues(:,ismember(simBestKeto.variables,'InsulinConcentrationPlasma')),'--','Color', BlueColor,'LineWidth',LineWidthValue-1)
fill(timeFill/10080,[minRegularInsulinConcentrationPlasma2', fliplr(minRegularInsulinConcentrationPlasma1')],'r','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestRegular.time/10080,simBestRegular.variablevalues(:,ismember(simBestRegular.variables,'InsulinConcentrationPlasma')),'r--','LineWidth',LineWidthValue-1)
ylim([0 300])
xlim([0 24])
hold off
grid on

%%
figure('Name', "Silfvergren et al. 2022", 'units', 'normalized', 'outerposition', [0 0 1 1])

%Ketones
%subplot(1,2,2)
hold on
set(gca,'ytick',[0,0.05,0.1,0.15],'xtick',[0,4,8,12,16,20,24],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Hepatic OOA'; '(g)'});
xlabel({'Time (weeks)'});
fill(timeFill/10080,[maxKetoOAA2', fliplr(minKetoOAA2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestKeto.time/10080,simBestKeto.statevalues(:,ismember(simBestKeto.states,'OAA_Liver')),'--','Color', BlueColor,'LineWidth',LineWidthValue-1)
fill(timeFill/10080,[maxRegularOAA2', fliplr(minRegularOAA2')],'r','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBestRegular.time/10080,simBestRegular.statevalues(:,ismember(simBestRegular.states,'OAA_Liver')),'r--','LineWidth',LineWidthValue-1)
%ylim([0 300])
xlim([0 24])
hold off
grid on
