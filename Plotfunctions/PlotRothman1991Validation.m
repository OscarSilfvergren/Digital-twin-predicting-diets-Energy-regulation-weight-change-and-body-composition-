load('Rothman1991_paramsCalibrated')

[sNames, ODEs, initialvalues] = IQMstates(model);
[pNames, param] = IQMparameters(model);

[row column] = size(Rothman1991_paramsCalibrated);
Rothman1991_paramsCalibrated = sortrows(Rothman1991_paramsCalibrated,column);

%% Cost to full estimation data

for j = 1:row
    
RothmanTemp              = Rothman1991_paramsCalibrated(j,1:AmountParameterOpti);
PersonSpecificParams(1)  = RothmanTemp(ismember(pNames,"carb_flow"));
PersonSpecificParams(2)  = RothmanTemp(ismember(pNames,"protein_flow"));
PersonSpecificParams(3)  = RothmanTemp(ismember(pNames,"lipids_flow"));

costValidationRothman = Rothman1991_costfunction(Rothman1991_data,11,log(PersonSpecificParams),RothmanTemp,model,initialvalues, pNames,sNames);

Rothman1991_paramsCalibrated(j,column) = costValidationRothman;
end

Rothman1991_paramsCalibrated = sortrows(Rothman1991_paramsCalibrated,column);

disp(' ')
if Rothman1991_paramsCalibrated(1,column) < chi2inv(0.95,11)
    disp('----- Fit to Estimation data is below statistical treshold -----')
end
fprintf('Best prediction of Rothman data: %.2f, Statistical Limit: %.2f (dgf = %i)', Rothman1991_paramsCalibrated(1,column), chi2inv(0.95,11), 11) % Finddatapoints
disp(' ')


%% Simulate Rothman1991

%%
[row column] = size(Rothman1991_paramsCalibrated);

for j = 1:row
    
    
    [optimizedParamTemp,initialvalues]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,Rothman1991_paramsCalibrated(j,1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames
    
SimSteadystate       = model(0:282240:282240,initialvalues,optimizedParamTemp);
    optimizedParamTemp(ismember(pNames,"carb_flow"))    = 0;
    optimizedParamTemp(ismember(pNames,"protein_flow")) = 0;
    optimizedParamTemp(ismember(pNames,"lipids_flow"))  = 0;
SimSteadystate      = model(0:360:360,SimSteadystate.statevalues(end,:),optimizedParamTemp);
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

%% core variables Gly, I, G & EGP 

figure('Name', "Rothman1991 et al. 2044 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

% Glycogen
hold on
grid on
set(gca,'xtick',[0,1,2,3],'ytick',[0,100,200,300,400,500],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Hepatic glycogen' ; '(mM)'});
xlabel('Time (days)')
fill(timeFill/1440,[maxGlycogenConcentrationliver2', fliplr(minGlycogenConcentrationliver2')],'b','FaceAlpha',0.2,'EdgeAlpha',0);
GlycogenSEM = Rothman1991_data.glycogenpSEM;
plot(simBest.time/1440, simBest.variablevalues(:,ismember(simBest.variables,'GlycogenConcentrationliver')),'--','Color', BlueColor,'LineWidth',LineWidthValue)
errorbar(Rothman1991_data.timepALL(1:1)/1440,Rothman1991_data.glycogenpALL(1:1),GlycogenSEM(1:1),' k x','MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
errorbar(Rothman1991_data.timepALL(1:11)/1440,Rothman1991_data.glycogenpALL(1:11),GlycogenSEM(1:11),' k .','MarkerSize',MarkerSizeErrorbars,'CapSize',CapSizeValue,'LineWidth',8);
xlim([0 3])
ylim([0 500])
legend('Model uncertainty','Best prediction','Calibration Data','Data','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize+3)
hold off

