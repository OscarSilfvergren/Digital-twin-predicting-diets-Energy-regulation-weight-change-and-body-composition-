load('Bray2012_paramsCalibrated')
load('Bray2012_data');

[sNames, ODEs, initialvalues] = IQMstates(model);
[pNames, param] = IQMparameters(model);

%% Clean data

Bray2012_data.time = Bray2012_data.time/1440;
Bray2012_data.time(1) = 1;
%% Cost to full estimation data

[row column] = size(Bray2012_paramsCalibrated);

for j = 1:row
BrayParamTemp                 = Bray2012_paramsCalibrated(j,:);
PersonSpecificParams(1)       = BrayParamTemp(ismember(pNames,"carb_flow"));
PersonSpecificParams(2)       = BrayParamTemp(ismember(pNames,"protein_flow"));
PersonSpecificParams(3)       = BrayParamTemp(ismember(pNames,"lipids_flow"));
PersonSpecificParams(4)       = BrayParamTemp(end-1);

costValidationBray = Bray2012_costfunction(Bray2012_data,13,log(PersonSpecificParams),BrayParamTemp(1:AmountParameterOpti),model,initialvalues, pNames,sNames);
Bray2012_paramsCalibrated(j,column) = costValidationBray;
end

Bray2012_paramsCalibrated = sortrows(Bray2012_paramsCalibrated,column);

disp(' ')
if Bray2012_paramsCalibrated(1,column) < chi2inv(0.95,13)
    disp('----- Prediction of Bray is below statistical treshold -----')
end
fprintf('Best prediction of Bray data: %.2f, Statistical Limit: %.2f (dgf = %i)', Bray2012_paramsCalibrated(1,column), chi2inv(0.95,13), 13) % Finddatapoints
disp(' ')

%% Simulate Rothman1991


%%
%%
[row column] = size(Bray2012_paramsCalibrated);

for j = 1:row
    
    [optimizedParamTemp,initialvalues]    = Persondeclaration(166, 0, 1,73.5,31.5,initialvalues,sNames,Bray2012_paramsCalibrated(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
    SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    BMR                = SimSteadystate.reactionvalues(:,ismember(SimSteadystate.reactions,'BMR'));
    IntakeFoodCalories = BMR(end)-600;
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = ((IntakeFoodCalories*0.5)/4)/1440;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = ((IntakeFoodCalories*0.25)/4)/1440;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = ((IntakeFoodCalories*0.25)/9)/1440;
    SimSteadystate     = model(0:2880:2880,SimSteadystate.statevalues(end,:),optimizedParamTemp);   
    SimBray2012          = model(0:1440:241920,SimSteadystate.statevalues(end,:),optimizedParamTemp); % 24 weeks
    
    if j == 1
        simBest = SimBray2012;
        
        maxBW1 = SimBray2012.variablevalues(:,ismember(SimBray2012.variables,'BW'));
        minBW1 = SimBray2012.variablevalues(:,ismember(SimBray2012.variables,'BW'));
        maxBW2 = SimBray2012.variablevalues(:,ismember(SimBray2012.variables,'BW'));
        minBW2 = SimBray2012.variablevalues(:,ismember(SimBray2012.variables,'BW'));
        
    end
    
    %Glucose
    maxBW1 = SimBray2012.variablevalues(:,ismember(SimBray2012.variables,'BW'));
    minBW1 = SimBray2012.variablevalues(:,ismember(SimBray2012.variables,'BW'));
    maxBW2 = max(maxBW2,maxBW1);
    minBW2 = min(minBW2,minBW1);
    
    
end

timeFill = [SimBray2012.time, fliplr(SimBray2012.time)];

%% Weight
figure('Name', "Bray et al. 2012", 'units', 'normalized', 'outerposition', [0 0 1 1])

hold on
grid on
set(gca,'xtick',[0,4,8,12,16,20,24],'ytick',[50,60,70,80,90,100,110],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Body Weight' ; '(kg)'});
xlabel('Time (weeks)')
y = [maxBW2', fliplr(minBW2')];
fill(timeFill/10080,y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBest.time/10080,simBest.variablevalues(:,ismember(simBest.variables,'BW')),'--','Color', BlueColor,'LineWidth',LineWidthValue-1)
plot(Bray2012_data.time(1:2)/7,Bray2012_data.weight(1:2),'Color', BlackColor,'LineStyle', 'n', 'Marker','x','MarkerSize',MarkerSize-15,'LineWidth',LineWidthValue)
plot(Bray2012_data.time(1:end)/7,Bray2012_data.weight(1:end),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize+10,'LineWidth',LineWidthValue)
legend('Model uncertainty','Best prediction','Calibration data','Validation data','Meals','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize+3)
ylim([80 110])
xlim([0 24])
hold off
