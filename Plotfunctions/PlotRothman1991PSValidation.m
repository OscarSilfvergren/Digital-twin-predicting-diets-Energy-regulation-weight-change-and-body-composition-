
[sNames, ODEs, initialvalues] = IQMstates(model);
[pNames, param] = IQMparameters(model);

%% P1
load('Rothman1991_P1paramsCalibrated')
[row column] = size(Rothman1991_P1paramsCalibrated);

for j = 1:row
    
RothmanTemp              = Rothman1991_P1paramsCalibrated(j,1:AmountParameterOpti);
PersonSpecificParams(1)  = RothmanTemp(ismember(pNames,"carb_flow"));
PersonSpecificParams(2)  = RothmanTemp(ismember(pNames,"protein_flow"));
PersonSpecificParams(3)  = RothmanTemp(ismember(pNames,"lipids_flow"));
PersonSpecificParams(4)  = RothmanTemp(end-1);
costValidationRothman    = Rothman1991Personspecific_costfunction(Rothman1991_data.glycogenp1(1:9),Rothman1991_data.timep1(1:9),9,log(PersonSpecificParams),RothmanTemp,model,initialvalues, pNames,sNames);

Rothman1991_P1paramsCalibrated(j,column) = costValidationRothman;
end

Rothman1991_P1paramsCalibrated = sortrows(Rothman1991_P1paramsCalibrated,column);
[row column] = size(Rothman1991_P1paramsCalibrated);

for j = 1:row
    
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,Rothman1991_P1paramsCalibrated(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
    optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK"))  = optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK")) * Rothman1991_P1paramsCalibrated(j,end-1);
    SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = 0;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = 0;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate      = model(0:480:480,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimRothman1991      = model(1:1:4320,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    
    if j == 1
        simBest = SimRothman1991;
        
        maxGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        maxGlycogenConcentrationliver2 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver2 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        
    end
    
    %GlycogenConcentrationliver
    maxGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
    minGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
    maxGlycogenConcentrationliver2 = max(maxGlycogenConcentrationliver2,maxGlycogenConcentrationliver1);
    minGlycogenConcentrationliver2 = min(minGlycogenConcentrationliver2,minGlycogenConcentrationliver1);
    
end

timeFill = [SimRothman1991.time, fliplr(SimRothman1991.time)];

figure('Name', "Rothman1991 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Glycogen
hold on
set(gca,'xtick',[0,1,2,3],'ytick',[0,100,200,300,400,500],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Hepatic Glycogen' ;'(mmol/L)'});
xlabel('Time (days)')
title('Participant 1')
fill(timeFill/1440,[maxGlycogenConcentrationliver2', fliplr(minGlycogenConcentrationliver2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(simBest.time/1440, simBest.variablevalues(:,ismember(simBest.variables,'GlycogenConcentrationliver')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep1(1:1)/1440,Rothman1991_data.glycogenp1(1:1),'Color', BlackColor,'LineStyle', 'n', 'Marker','x','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep1(1:end)/1440,Rothman1991_data.glycogenp1(1:end),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
xlim([0 3])
ylim([0 530])
legend('Model uncertainty','Best prediction','Calibration data','Person-specific data','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize+3)
hold off
grid on

%% P2
load('Rothman1991_P2paramsCalibrated')
[row column] = size(Rothman1991_P2paramsCalibrated);

for j = 1:row
    
RothmanTemp              = Rothman1991_P2paramsCalibrated(j,1:AmountParameterOpti);
PersonSpecificParams(1)  = RothmanTemp(ismember(pNames,"carb_flow"));
PersonSpecificParams(2)  = RothmanTemp(ismember(pNames,"protein_flow"));
PersonSpecificParams(3)  = RothmanTemp(ismember(pNames,"lipids_flow"));
costValidationRothman    = Rothman1991Personspecific_costfunction(Rothman1991_data.glycogenp2(1:11),Rothman1991_data.timep2(1:11),11,log(PersonSpecificParams),RothmanTemp,model,initialvalues, pNames,sNames);

Rothman1991_P2paramsCalibrated(j,column) = costValidationRothman;
end

Rothman1991_P2paramsCalibrated = sortrows(Rothman1991_P2paramsCalibrated,column);
[row column] = size(Rothman1991_P2paramsCalibrated);

for j = 1:row
    
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,Rothman1991_P2paramsCalibrated(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
    optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK"))  = optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK")) * Rothman1991_P2paramsCalibrated(j,end-1);
    SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = 0;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = 0;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate      = model(0:480:480,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimRothman1991      = model(1:1:4320,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    
    if j == 1
        simBest = SimRothman1991;
        
        maxGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        maxGlycogenConcentrationliver2 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver2 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        
    end
    
    %GlycogenConcentrationliver
    maxGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
    minGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
    maxGlycogenConcentrationliver2 = max(maxGlycogenConcentrationliver2,maxGlycogenConcentrationliver1);
    minGlycogenConcentrationliver2 = min(minGlycogenConcentrationliver2,minGlycogenConcentrationliver1);
    
end

timeFill = [SimRothman1991.time, fliplr(SimRothman1991.time)];

figure('Name', "Rothman1991 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Glycogen
hold on
set(gca,'xtick',[0,1,2,3],'ytick',[0,100,200,300,400,500],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Hepatic Glycogen' ;'(mmol/L)'});
xlabel('Time (days)')
title('Participant 2')
fill(timeFill/1440,[maxGlycogenConcentrationliver2', fliplr(minGlycogenConcentrationliver2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(Rothman1991_data.timep2(1:1)/1440,Rothman1991_data.glycogenp2(1:1),'Color', BlackColor,'LineStyle', 'n', 'Marker','x','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep2(1:end)/1440,Rothman1991_data.glycogenp2(1:end),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
plot(simBest.time/1440, simBest.variablevalues(:,ismember(simBest.variables,'GlycogenConcentrationliver')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
xlim([0 3])
ylim([0 530])
hold off
grid on

%% p3
load('Rothman1991_P3paramsCalibrated')
[row column] = size(Rothman1991_P3paramsCalibrated);

for j = 1:row
    
RothmanTemp              = Rothman1991_P3paramsCalibrated(j,1:AmountParameterOpti);
PersonSpecificParams(1)  = RothmanTemp(ismember(pNames,"carb_flow"));
PersonSpecificParams(2)  = RothmanTemp(ismember(pNames,"protein_flow"));
PersonSpecificParams(3)  = RothmanTemp(ismember(pNames,"lipids_flow"));
costValidationRothman    = Rothman1991Personspecific_costfunction(Rothman1991_data.glycogenp3(1:10),Rothman1991_data.timep3(1:10),10,log(PersonSpecificParams),RothmanTemp,model,initialvalues, pNames,sNames);

Rothman1991_P3paramsCalibrated(j,column) = costValidationRothman;
end

Rothman1991_P3paramsCalibrated = sortrows(Rothman1991_P3paramsCalibrated,column);
[row column] = size(Rothman1991_P3paramsCalibrated);

for j = 1:row
    
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,Rothman1991_P3paramsCalibrated(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
        optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK"))  = optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK")) * Rothman1991_P3paramsCalibrated(j,end-1);
    SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = 0;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = 0;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate      = model(0:480:480,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimRothman1991      = model(1:1:4320,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    
    if j == 1
        simBest = SimRothman1991;
        
        maxGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        maxGlycogenConcentrationliver2 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver2 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        
    end
    
    %GlycogenConcentrationliver
    maxGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
    minGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
    maxGlycogenConcentrationliver2 = max(maxGlycogenConcentrationliver2,maxGlycogenConcentrationliver1);
    minGlycogenConcentrationliver2 = min(minGlycogenConcentrationliver2,minGlycogenConcentrationliver1);
    
end

timeFill = [SimRothman1991.time, fliplr(SimRothman1991.time)];

figure('Name', "Rothman1991 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Glycogen
hold on
set(gca,'xtick',[0,1,2,3],'ytick',[0,100,200,300,400,500],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Hepatic Glycogen' ;'(mmol/L)'});
xlabel('Time (days)')
title('Participant 3')
fill(timeFill/1440,[maxGlycogenConcentrationliver2', fliplr(minGlycogenConcentrationliver2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(Rothman1991_data.timep3(1:1)/1440,Rothman1991_data.glycogenp3(1:1),'Color', BlackColor,'LineStyle', 'n', 'Marker','x','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep3(1:end)/1440,Rothman1991_data.glycogenp3(1:end),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
plot(simBest.time/1440, simBest.variablevalues(:,ismember(simBest.variables,'GlycogenConcentrationliver')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
xlim([0 3])
ylim([0 530])
hold off
grid on

%% P4
load('Rothman1991_P4paramsCalibrated')
[row column] = size(Rothman1991_P4paramsCalibrated);

for j = 1:row
    
RothmanTemp              = Rothman1991_P4paramsCalibrated(j,1:AmountParameterOpti);
PersonSpecificParams(1)  = RothmanTemp(ismember(pNames,"carb_flow"));
PersonSpecificParams(2)  = RothmanTemp(ismember(pNames,"protein_flow"));
PersonSpecificParams(3)  = RothmanTemp(ismember(pNames,"lipids_flow"));
costValidationRothman    = Rothman1991Personspecific_costfunction(Rothman1991_data.glycogenp4(1:10),Rothman1991_data.timep4(1:10),10,log(PersonSpecificParams),RothmanTemp,model,initialvalues, pNames,sNames);

Rothman1991_P4paramsCalibrated(j,column) = costValidationRothman;
end

Rothman1991_P4paramsCalibrated = sortrows(Rothman1991_P4paramsCalibrated,column);
[row column] = size(Rothman1991_P4paramsCalibrated);

for j = 1:row
    
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,Rothman1991_P4paramsCalibrated(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
    optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK"))  = optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK")) * Rothman1991_P4paramsCalibrated(j,end-1);
    SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = 0;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = 0;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate      = model(0:480:480,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimRothman1991      = model(1:1:4320,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    
    if j == 1
        simBest = SimRothman1991;
        
        maxGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        maxGlycogenConcentrationliver2 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver2 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        
    end
    
    %GlycogenConcentrationliver
    maxGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
    minGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
    maxGlycogenConcentrationliver2 = max(maxGlycogenConcentrationliver2,maxGlycogenConcentrationliver1);
    minGlycogenConcentrationliver2 = min(minGlycogenConcentrationliver2,minGlycogenConcentrationliver1);
    
end

timeFill = [SimRothman1991.time, fliplr(SimRothman1991.time)];

figure('Name', "Rothman1991 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Glycogen
hold on
set(gca,'xtick',[0,1,2,3],'ytick',[0,100,200,300,400,500],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Hepatic Glycogen' ;'(mmol/L)'});
xlabel('Time (days)')
title('Participant 4')
fill(timeFill/1440,[maxGlycogenConcentrationliver2', fliplr(minGlycogenConcentrationliver2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(Rothman1991_data.timep4(1:1)/1440,Rothman1991_data.glycogenp4(1:1),'Color', BlackColor,'LineStyle', 'n', 'Marker','x','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep4(1:end)/1440,Rothman1991_data.glycogenp4(1:end),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
plot(simBest.time/1440, simBest.variablevalues(:,ismember(simBest.variables,'GlycogenConcentrationliver')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
xlim([0 3])
ylim([0 530])
hold off
grid on

%% P5
load('Rothman1991_P5paramsCalibrated')
[row column] = size(Rothman1991_P5paramsCalibrated);

for j = 1:row
    
RothmanTemp              = Rothman1991_P5paramsCalibrated(j,1:AmountParameterOpti);
PersonSpecificParams(1)  = RothmanTemp(ismember(pNames,"carb_flow"));
PersonSpecificParams(2)  = RothmanTemp(ismember(pNames,"protein_flow"));
PersonSpecificParams(3)  = RothmanTemp(ismember(pNames,"lipids_flow"));
costValidationRothman    = Rothman1991Personspecific_costfunction(Rothman1991_data.glycogenp5(1:8),Rothman1991_data.timep5(1:8),8,log(PersonSpecificParams),RothmanTemp,model,initialvalues, pNames,sNames);

Rothman1991_P5paramsCalibrated(j,column) = costValidationRothman;
end

Rothman1991_P5paramsCalibrated = sortrows(Rothman1991_P5paramsCalibrated,column);
[row column] = size(Rothman1991_P5paramsCalibrated);

for j = 1:row
    
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,Rothman1991_P5paramsCalibrated(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
        optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK"))  = optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK")) * Rothman1991_P5paramsCalibrated(j,end-1);
    SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = 0;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = 0;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate      = model(0:480:480,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimRothman1991      = model(1:1:4320,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    
    if j == 1
        simBest = SimRothman1991;
        
        maxGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        maxGlycogenConcentrationliver2 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver2 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        
    end
    
    %GlycogenConcentrationliver
    maxGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
    minGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
    maxGlycogenConcentrationliver2 = max(maxGlycogenConcentrationliver2,maxGlycogenConcentrationliver1);
    minGlycogenConcentrationliver2 = min(minGlycogenConcentrationliver2,minGlycogenConcentrationliver1);
    
end

timeFill = [SimRothman1991.time, fliplr(SimRothman1991.time)];

figure('Name', "Rothman1991 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Glycogen
hold on
set(gca,'xtick',[0,1,2,3],'ytick',[0,100,200,300,400,500],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Hepatic Glycogen' ;'(mmol/L)'});
xlabel('Time (days)')
title('Participant 5')
fill(timeFill/1440,[maxGlycogenConcentrationliver2', fliplr(minGlycogenConcentrationliver2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(Rothman1991_data.timep5(1:1)/1440,Rothman1991_data.glycogenp5(1:1),'Color', BlackColor,'LineStyle', 'n', 'Marker','x','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep5(1:end)/1440,Rothman1991_data.glycogenp5(1:end),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
plot(simBest.time/1440, simBest.variablevalues(:,ismember(simBest.variables,'GlycogenConcentrationliver')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
xlim([0 3])
ylim([0 530])
hold off
grid on

%% P6
load('Rothman1991_P6paramsCalibrated')
[row column] = size(Rothman1991_P6paramsCalibrated);

for j = 1:row
    
RothmanTemp              = Rothman1991_P6paramsCalibrated(j,1:AmountParameterOpti);
PersonSpecificParams(1)  = RothmanTemp(ismember(pNames,"carb_flow"));
PersonSpecificParams(2)  = RothmanTemp(ismember(pNames,"protein_flow"));
PersonSpecificParams(3)  = RothmanTemp(ismember(pNames,"lipids_flow"));
costValidationRothman    = Rothman1991Personspecific_costfunction(Rothman1991_data.glycogenp6(1:11),Rothman1991_data.timep6(1:11),11,log(PersonSpecificParams),RothmanTemp,model,initialvalues, pNames,sNames);

Rothman1991_P6paramsCalibrated(j,column) = costValidationRothman;
end

Rothman1991_P6paramsCalibrated = sortrows(Rothman1991_P6paramsCalibrated,column);
[row column] = size(Rothman1991_P6paramsCalibrated);

for j = 1:row
    
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,Rothman1991_P6paramsCalibrated(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
            optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK"))  = optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK")) * Rothman1991_P6paramsCalibrated(j,end-1);
    SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = 0;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = 0;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate      = model(0:480:480,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimRothman1991      = model(1:1:4320,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    
    if j == 1
        simBest = SimRothman1991;
        
        maxGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        maxGlycogenConcentrationliver2 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver2 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        
    end
    
    %GlycogenConcentrationliver
    maxGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
    minGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
    maxGlycogenConcentrationliver2 = max(maxGlycogenConcentrationliver2,maxGlycogenConcentrationliver1);
    minGlycogenConcentrationliver2 = min(minGlycogenConcentrationliver2,minGlycogenConcentrationliver1);
    
end

timeFill = [SimRothman1991.time, fliplr(SimRothman1991.time)];

figure('Name', "Rothman1991 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Glycogen
hold on
set(gca,'xtick',[0,1,2,3],'ytick',[0,100,200,300,400,500],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Hepatic Glycogen' ;'(mmol/L)'});
xlabel('Time (days)')
title('Participant 6')
fill(timeFill/1440,[maxGlycogenConcentrationliver2', fliplr(minGlycogenConcentrationliver2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(Rothman1991_data.timep6(1:1)/1440,Rothman1991_data.glycogenp6(1:1),'Color', BlackColor,'LineStyle', 'n', 'Marker','x','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep6(1:end)/1440,Rothman1991_data.glycogenp6(1:end),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
plot(simBest.time/1440, simBest.variablevalues(:,ismember(simBest.variables,'GlycogenConcentrationliver')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
xlim([0 3])
ylim([0 530])
hold off
grid on

%% P7
load('Rothman1991_P7paramsCalibrated')
[row column] = size(Rothman1991_P7paramsCalibrated);

for j = 1:row
    
RothmanTemp              = Rothman1991_P7paramsCalibrated(j,1:AmountParameterOpti);
PersonSpecificParams(1)  = RothmanTemp(ismember(pNames,"carb_flow"));
PersonSpecificParams(2)  = RothmanTemp(ismember(pNames,"protein_flow"));
PersonSpecificParams(3)  = RothmanTemp(ismember(pNames,"lipids_flow"));
costValidationRothman    = Rothman1991Personspecific_costfunction(Rothman1991_data.glycogenp7(1:9),Rothman1991_data.timep7(1:9),9,log(PersonSpecificParams),RothmanTemp,model,initialvalues, pNames,sNames);

Rothman1991_P7paramsCalibrated(j,column) = costValidationRothman;
end

Rothman1991_P7paramsCalibrated = sortrows(Rothman1991_P7paramsCalibrated,column);
[row column] = size(Rothman1991_P7paramsCalibrated);

for j = 1:row
    
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,Rothman1991_P7paramsCalibrated(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
                optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK"))  = optimizedParamTemp(ismember(pNames,"InsulinResponseLiverK")) * Rothman1991_P7paramsCalibrated(j,end-1);
    SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = 0;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = 0;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate      = model(0:480:480,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    SimRothman1991      = model(1:1:4320,SimSteadystate.statevalues(end,:),optimizedParamTemp);
    
    if j == 1
        simBest = SimRothman1991;
        
        maxGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        maxGlycogenConcentrationliver2 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        minGlycogenConcentrationliver2 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
        
    end
    
    %GlycogenConcentrationliver
    maxGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
    minGlycogenConcentrationliver1 = SimRothman1991.variablevalues(:,ismember(SimRothman1991.variables,'GlycogenConcentrationliver'));
    maxGlycogenConcentrationliver2 = max(maxGlycogenConcentrationliver2,maxGlycogenConcentrationliver1);
    minGlycogenConcentrationliver2 = min(minGlycogenConcentrationliver2,minGlycogenConcentrationliver1);
    
end

timeFill = [SimRothman1991.time, fliplr(SimRothman1991.time)];

figure('Name', "Rothman1991 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])
%%
% Glycogen
hold on
set(gca,'xtick',[0,1,2,3],'ytick',[0,100,200,300,400,500],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Hepatic Glycogen' ;'(mmol/L)'});
xlabel('Time (days)')
title('Participant 7')
fill(timeFill/1440,[maxGlycogenConcentrationliver2', fliplr(minGlycogenConcentrationliver2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(Rothman1991_data.timep7(1:1)/1440,Rothman1991_data.glycogenp7(1:1),'Color', BlackColor,'LineStyle', 'n', 'Marker','x','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep7(1:end)/1440,Rothman1991_data.glycogenp7(1:end),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
plot(simBest.time/1440, simBest.variablevalues(:,ismember(simBest.variables,'GlycogenConcentrationliver')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
xlim([0 3])
ylim([0 530])
hold off
grid on

