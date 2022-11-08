function Cost_Model = Silfvergren2022_costfunction(Silfvergren2022Weight_data,Silfvergren2022GlucoseBloodSensor_data,Silfvergren2022Ketones_data,Silfvergren2022_meals,AmountParameterOpti,params,model,initialvalues, pNames,sNames)

params=exp(params);

%% Simulate Silfvergren2022 - Study 1

[Silfvergren2022_params,start]   = Persondeclaration(180, 1, 0,75,10,initialvalues,sNames,params(1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames

try
    SimSteadystate       = model(0:282240:282240,start,Silfvergren2022_params);
    Silfvergren2022_params(ismember(pNames,"carb_flow"))    = 0;
    Silfvergren2022_params(ismember(pNames,"protein_flow")) = 0;
    Silfvergren2022_params(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate       = model(0:480:480,SimSteadystate.statevalues(end,:),Silfvergren2022_params);
    
    Silfvergren2022_paramsMeal  = Mealdeclaration(Silfvergren2022_meals.CarbohydratesAmount(1),Silfvergren2022_meals.ProteinAmount(1),Silfvergren2022_meals.LipidsAmount(1),Silfvergren2022_meals.MealDuration(1),Silfvergren2022_params,pNames);       % CarbohydratesAmount, ProteinAmount, LipidsAmount, MealLength,
    SimUpToMeal                 = model(0:120:Silfvergren2022_meals.MealStart(1),SimSteadystate.statevalues(end,:),Silfvergren2022_params);
    SimMeal                     = model(Silfvergren2022_meals.MealStart(1):120:Silfvergren2022_meals.MealEnd(1),SimUpToMeal.statevalues(end,:),Silfvergren2022_paramsMeal);
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
catch error
    Cost_Model = 1e60;
    return
end

% Silfvergren2022 - Cost
Weight_Simulated    = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'BW'));
Weight_Simulated    = Weight_Simulated(Silfvergren2022Weight_data.Time);
CostWeight          = ((Silfvergren2022Weight_data.Weight(1:113)    - Weight_Simulated(1:113)).^2) ./((Silfvergren2022Weight_data.Weight(1:113)*0.05).^2);
CostWeight          = sum(CostWeight,'omitnan') + ((Silfvergren2022Weight_data.Weight(113)    - Weight_Simulated(113)).^2) ./((Silfvergren2022Weight_data.Weight(113)*0.01).^2);
TotalFat_Simulated   = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'TGA_AdipocyteKg'));
TotalFat_Simulated   = TotalFat_Simulated(Silfvergren2022Weight_data.Time);
CostTotalFat         = ((Silfvergren2022Weight_data.Fatmass(1:113)    - TotalFat_Simulated(1:113)).^2) ./((Silfvergren2022Weight_data.Fatmass(1:113)*0.15).^2);
CostTotalFat         = sum(CostTotalFat,'omitnan') + ((Silfvergren2022Weight_data.Fatmass(113)    - TotalFat_Simulated(113)).^2) ./((Silfvergren2022Weight_data.Fatmass(113)*0.1).^2);
TotalLeanmass_Simulated    = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'LeanBodyWeightKg'));
TotalLeanmass_Simulated    = TotalLeanmass_Simulated(Silfvergren2022Weight_data.Time);
CostTotalLeanmass          = ((Silfvergren2022Weight_data.Leanmass(1:113)    - TotalLeanmass_Simulated(1:113)).^2) ./((Silfvergren2022Weight_data.Leanmass(1:113)*0.05).^2);
CostTotalLeanmass          = sum(CostTotalLeanmass,'omitnan') + ((Silfvergren2022Weight_data.Leanmass(113)    - TotalLeanmass_Simulated(113)).^2) ./((Silfvergren2022Weight_data.Leanmass(113)*0.01).^2);
Ketones_Simulated   = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'KetoneConcentrationPlasma'));
Ketones_Simulated   = Ketones_Simulated(Silfvergren2022Ketones_data.Time_minutesKetones);
KetonesSEM          = Silfvergren2022Ketones_data.Ketones*0.2 + 0.01;
CostKetones         = ((Silfvergren2022Ketones_data.Ketones(1:66)    - Ketones_Simulated(1:66)).^2) ./(KetonesSEM(1:66).^2);
CostKetones         = sum(CostKetones,'omitnan');
Glucose_Simulated   = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'GlucoseConcentrationPlasma'));
Glucose_Simulated   = Glucose_Simulated(Silfvergren2022GlucoseBloodSensor_data.glucoseTime);
CostGlucose         = ((Silfvergren2022GlucoseBloodSensor_data.bloodvalue(1:112)    - Glucose_Simulated(1:112)).^2) ./((Silfvergren2022GlucoseBloodSensor_data.bloodvalue(1:112)*0.15).^2);
CostGlucose         = sum(CostGlucose,'omitnan');

Cost_Model = CostKetones + CostWeight + CostTotalFat + CostTotalLeanmass + CostGlucose; %517 datapoints

