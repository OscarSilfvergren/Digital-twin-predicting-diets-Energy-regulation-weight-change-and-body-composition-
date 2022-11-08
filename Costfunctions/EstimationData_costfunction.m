function Cost_Model = EstimationData_costfunction(Silfvergren2022Weight_data,Silfvergren2022GlucoseBloodSensor_data,Silfvergren2022Ketones_data,Silfvergren2022_meals,Krssak2004Healthy_data,DallaMan2007_data,AmountParameterOpti,params,model,initialvalues, pNames,sNames,fid)

params=exp(params);

%% Declare individual params

% Study 2
Study2_SteadystateMeal                   = [params(AmountParameterOpti+1) , params(AmountParameterOpti+2),  params(AmountParameterOpti+3)];
Study2_BMRPersonalScaling                = params(AmountParameterOpti+4);
Study2_InsulinResistancePersonalScaling  = params(AmountParameterOpti+5);
Study2_InsulinClearancePersonalScaling   = params(AmountParameterOpti+6);
Study2_InsulinResponsePersonalScaling    = params(AmountParameterOpti+7);

% Study 3
Study3_SteadystateMeal                   = [params(AmountParameterOpti+8) , params(AmountParameterOpti+9),  params(AmountParameterOpti+10)];
Study3_BMRPersonalScaling                = params(AmountParameterOpti+11);
Study3_InsulinResistancePersonalScaling  = params(AmountParameterOpti+12);
Study3_InsulinClearancePersonalScaling   = params(AmountParameterOpti+13);
Study3_InsulinResponsePersonalScaling    = params(AmountParameterOpti+14);

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
CostWeight          = ((Silfvergren2022Weight_data.Weight(1:56)    - Weight_Simulated(1:56)).^2) ./((Silfvergren2022Weight_data.Weight(1:56)*0.05).^2);
CostWeight          = sum(CostWeight,'omitnan');
TotalFat_Simulated   = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'TGA_AdipocyteKg'));
TotalFat_Simulated   = TotalFat_Simulated(Silfvergren2022Weight_data.Time);
CostTotalFat         = ((Silfvergren2022Weight_data.Fatmass(1:56)    - TotalFat_Simulated(1:56)).^2) ./((Silfvergren2022Weight_data.Fatmass(1:56)*0.15).^2);
CostTotalFat         = sum(CostTotalFat,'omitnan');
TotalLeanmass_Simulated    = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'LeanBodyWeightKg'));
TotalLeanmass_Simulated    = TotalLeanmass_Simulated(Silfvergren2022Weight_data.Time);
CostTotalLeanmass          = ((Silfvergren2022Weight_data.Leanmass(1:56)    - TotalLeanmass_Simulated(1:56)).^2) ./((Silfvergren2022Weight_data.Leanmass(1:56)*0.05).^2);
CostTotalLeanmass          = sum(CostTotalLeanmass,'omitnan');
Ketones_Simulated   = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'KetoneConcentrationPlasma'));
Ketones_Simulated   = Ketones_Simulated(Silfvergren2022Ketones_data.Time_minutesKetones);
KetonesSEM          = Silfvergren2022Ketones_data.Ketones*0.2 + 0.01;
CostKetones         = ((Silfvergren2022Ketones_data.Ketones(1:33)    - Ketones_Simulated(1:33)).^2) ./(KetonesSEM(1:33).^2);
CostKetones         = sum(CostKetones,'omitnan');
Glucose_Simulated   = simSilfvergren2022.variablevalues(:,ismember(simSilfvergren2022.variables,'GlucoseConcentrationPlasma'));
Glucose_Simulated   = Glucose_Simulated(Silfvergren2022GlucoseBloodSensor_data.glucoseTime);
CostGlucose         = ((Silfvergren2022GlucoseBloodSensor_data.bloodvalue(1:56)    - Glucose_Simulated(1:56)).^2) ./((Silfvergren2022GlucoseBloodSensor_data.bloodvalue(1:56)*0.15).^2);
CostGlucose         = sum(CostGlucose,'omitnan');

Cost_Model = CostKetones + CostWeight + CostTotalFat + CostTotalLeanmass + CostGlucose; %258 datapoints

%% Krssak2004 - Study 2

[Krssak2004_params,start]   = Persondeclaration(180, 1, 0,70,5,initialvalues,sNames,params(1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames

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
Krssak2004_params(ismember(pNames,"Vm_id"))                     = Krssak2004_params(ismember(pNames,'Vm_id'))                   * Study2_InsulinResistancePersonalScaling;
Krssak2004_params(ismember(pNames,"Vf_id"))                     = Krssak2004_params(ismember(pNames,'Vf_id'))                   * Study2_InsulinResistancePersonalScaling;
Krssak2004_params(ismember(pNames,"InsulinResponseFATP1_K"))   = Krssak2004_params(ismember(pNames,'InsulinResponseFATP1_K')) * Study2_InsulinResistancePersonalScaling;
Krssak2004_params(ismember(pNames,"InsulinResponseLiverK"))     = Krssak2004_params(ismember(pNames,'InsulinResponseLiverK'))   * Study2_InsulinResistancePersonalScaling;
% InsulinClearancePersonalScaling
Krssak2004_params(ismember(pNames,"InsulinDegradation_PlasmaK"))    = Krssak2004_params(ismember(pNames,'InsulinDegradation_PlasmaK'))    * Study2_InsulinClearancePersonalScaling;
Krssak2004_params(ismember(pNames,"HE_positiveK"))                  = Krssak2004_params(ismember(pNames,'HE_positiveK'))                  * Study2_InsulinClearancePersonalScaling;
% InsulinClearancePersonalScaling
Krssak2004_params(ismember(pNames,"InsulinStabilisationMagnitude"))    = Krssak2004_params(ismember(pNames,'InsulinStabilisationMagnitude'))    * Study2_InsulinResponsePersonalScaling;
Krssak2004_params(ismember(pNames,"InsulinProductionK"))               = Krssak2004_params(ismember(pNames,'InsulinProductionK'))               * Study2_InsulinResponsePersonalScaling;

Krssak2004_paramsMeal               = Mealdeclaration(87,23,24,15,Krssak2004_params,pNames);                           % CarbohydratesAmount, ProteinAmount, LipidsAmount, MealLength,

try
    SimSteadystate       = model(0:282240:282240,start,Krssak2004_params);
    Krssak2004_params(ismember(pNames,"carb_flow"))    = 0;
    Krssak2004_params(ismember(pNames,"protein_flow")) = 0;
    Krssak2004_params(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate       = model(0:360:360,SimSteadystate.statevalues(end,:),Krssak2004_params);
    SimKrssakToMeal      = model(1:1:360,SimSteadystate.statevalues(end,:),Krssak2004_params);
    SimKrssak_meal       = model(360:1:375,SimKrssakToMeal.statevalues(end,:),Krssak2004_paramsMeal);
    SimKrssakAfterMeal   = model(375:1:1260,SimKrssak_meal.statevalues(end,:),Krssak2004_params);

catch error
    Cost_Model = 1e60;
    return
end
simKrssak2004Healthy.time           = 1:1:1260;
simKrssak2004Healthy.variablevalues = [SimKrssakToMeal.variablevalues; SimKrssak_meal.variablevalues(2:end,:); SimKrssakAfterMeal.variablevalues(2:end,:)];
simKrssak2004Healthy.reactionvalues = [SimKrssakToMeal.reactionvalues; SimKrssak_meal.reactionvalues(2:end,:); SimKrssakAfterMeal.reactionvalues(2:end,:)];
simKrssak2004Healthy.variables      = [SimKrssak_meal.variables];
simKrssak2004Healthy.reactions      = [SimKrssak_meal.reactions];


% Krssak2004 - Cost
simGly = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'GlycogenConcentrationliver'));
simG   = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'GlucoseConcentrationPlasma'));
simI   = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'InsulinConcentrationPlasma'));
simEGP = simKrssak2004Healthy.reactionvalues(:,ismember(simKrssak2004Healthy.reactions,'EGP'));
simFFA = simKrssak2004Healthy.variablevalues(:,ismember(simKrssak2004Healthy.variables,'FFAConcentrationPlasma'));
simGly = simGly(Krssak2004Healthy_data.timeGlycogen(1:22)+360);
simG   = simG(Krssak2004Healthy_data.timeGlucose(1:25)+360);
simI   = simI(Krssak2004Healthy_data.timeInsulin(1:10)+360);
simEGP = simEGP(Krssak2004Healthy_data.timeEGP(1:20)+360);
simFFA = simFFA(Krssak2004Healthy_data.timeFFA(1:10)+360);
CostGly = ((Krssak2004Healthy_data.glycogen(1:22)- simGly).^2)./(Krssak2004Healthy_data.glycogenSEM(1:22).^2);
CostG   = ((Krssak2004Healthy_data.glucose(1:25) - simG).^2)  ./(Krssak2004Healthy_data.glucoseSEM(1:25).^2);
CostI   = ((Krssak2004Healthy_data.insulin(1:10) - simI).^2)  ./((Krssak2004Healthy_data.insulinSEM(1:10)).^2);
CostEGP = ((Krssak2004Healthy_data.EGP(1:20)     - simEGP).^2)./(Krssak2004Healthy_data.EGPSEM(1:20).^2);
CostFFA = ((Krssak2004Healthy_data.FFA(1:10)     - simFFA).^2)./(Krssak2004Healthy_data.FFASEM(1:10).^2);

Cost_Model =  Cost_Model + sum(CostG,'omitnan') + sum(CostI,'omitnan') + sum(CostGly,'omitnan') + sum(CostEGP,'omitnan') + sum(CostFFA,'omitnan'); %87 datapoints

%% DallaMan2007 - Study 3

[DallaMan2007_params,start]   = Persondeclaration(170, 1, 0,70,8,initialvalues,sNames,params(1:AmountParameterOpti),pNames); % Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames

% Declare personal scaling
% Meals
DallaMan2007_params(ismember(pNames,"carb_flow"))    = Study3_SteadystateMeal(1);
DallaMan2007_params(ismember(pNames,"protein_flow")) = Study3_SteadystateMeal(2);
DallaMan2007_params(ismember(pNames,"lipids_flow"))  = Study3_SteadystateMeal(3);
% BMR
DallaMan2007_params(ismember(pNames,"BMRDynamicBW"))    = DallaMan2007_params(ismember(pNames,'BMRDynamicBW'))  * Study3_BMRPersonalScaling;
DallaMan2007_params(ismember(pNames,"BMRDynamicGly"))   = DallaMan2007_params(ismember(pNames,'BMRDynamicGly')) * Study3_BMRPersonalScaling;
DallaMan2007_params(ismember(pNames,"BMRbasal"))        = DallaMan2007_params(ismember(pNames,'BMRbasal'))      * Study3_BMRPersonalScaling;
% InsulinResistancePersonalScaling
DallaMan2007_params(ismember(pNames,"Vm_id"))                          = DallaMan2007_params(ismember(pNames,'Vm_id'))                   * Study3_InsulinResistancePersonalScaling;
DallaMan2007_params(ismember(pNames,"Vf_id"))                          = DallaMan2007_params(ismember(pNames,'Vf_id'))                   * Study3_InsulinResistancePersonalScaling;
DallaMan2007_params(ismember(pNames,"InsulinResponseFATP1_K"))        = DallaMan2007_params(ismember(pNames,'InsulinResponseFATP1_K')) * Study3_InsulinResistancePersonalScaling;
DallaMan2007_params(ismember(pNames,"InsulinResponseLiverK"))          = DallaMan2007_params(ismember(pNames,'InsulinResponseLiverK'))   * Study3_InsulinResistancePersonalScaling;
% InsulinClearancePersonalScaling
DallaMan2007_params(ismember(pNames,"InsulinDegradation_PlasmaK"))    = DallaMan2007_params(ismember(pNames,'InsulinDegradation_PlasmaK'))    * Study3_InsulinClearancePersonalScaling;
DallaMan2007_params(ismember(pNames,"HE_positiveK"))                  = DallaMan2007_params(ismember(pNames,'HE_positiveK'))                  * Study3_InsulinClearancePersonalScaling;
% Study3_InsulinResponsePersonalScaling
DallaMan2007_params(ismember(pNames,"InsulinStabilisationMagnitude")) = DallaMan2007_params(ismember(pNames,'InsulinStabilisationMagnitude'))    * Study3_InsulinResponsePersonalScaling;
DallaMan2007_params(ismember(pNames,"InsulinProductionK"))            = DallaMan2007_params(ismember(pNames,'InsulinProductionK'))               * Study3_InsulinResponsePersonalScaling;

DallaMan2007_paramsMeal               = Mealdeclaration(78,0,0,5,DallaMan2007_params,pNames);                           % CarbohydratesAmount, ProteinAmount, LipidsAmount, MealLength,

try
    
    SimSteadystate       = model(0:282240:282240,start,DallaMan2007_params);
    DallaMan2007_params(ismember(pNames,"carb_flow"))    = 0;
    DallaMan2007_params(ismember(pNames,"protein_flow")) = 0;
    DallaMan2007_params(ismember(pNames,"lipids_flow"))  = 0;
    SimSteadystate       = model(0:360:360,SimSteadystate.statevalues(end,:),DallaMan2007_params);
    SimDallaMan2007ToMeal      = model(1:1:360,SimSteadystate.statevalues(end,:),DallaMan2007_params);
    SimDallaMan2007_meal       = model(360:1:365,SimDallaMan2007ToMeal.statevalues(end,:),DallaMan2007_paramsMeal);
    SimDallaMan2007AfterMeal   = model(365:1:840,SimDallaMan2007_meal.statevalues(end,:),DallaMan2007_params);
catch error
    Cost_Model = 1e60;
    return
end
SimDallaMan2007.time           = 1:1:840;
SimDallaMan2007.statevalues    = [SimDallaMan2007ToMeal.statevalues; SimDallaMan2007_meal.statevalues(2:end,:)   ; SimDallaMan2007AfterMeal.statevalues(2:end,:)];
SimDallaMan2007.variablevalues = [SimDallaMan2007ToMeal.variablevalues; SimDallaMan2007_meal.variablevalues(2:end,:); SimDallaMan2007AfterMeal.variablevalues(2:end,:)];
SimDallaMan2007.reactionvalues = [SimDallaMan2007ToMeal.reactionvalues; SimDallaMan2007_meal.reactionvalues(2:end,:); SimDallaMan2007AfterMeal.reactionvalues(2:end,:)];
SimDallaMan2007.states         = [SimDallaMan2007_meal.states];
SimDallaMan2007.variables      = [SimDallaMan2007_meal.variables];
SimDallaMan2007.reactions      = [SimDallaMan2007_meal.reactions];
SimDallaMan2007.pNames         = [pNames];

% DallaMan2007 - Cost
simG = SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'GlucoseConcentrationPlasma'));
simI   = SimDallaMan2007.variablevalues(:,ismember(SimDallaMan2007.variables,'InsulinConcentrationPlasma'));
simEGP = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'EGP'));
simRa = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'Carbohydrates_Intestines_out'));
simU = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'U'));
simS = SimDallaMan2007.reactionvalues(:,ismember(SimDallaMan2007.reactions,'S'));
Glucose_time = round(DallaMan2007_data.Glucose_time(1:21),0)+360;
Insulin_time = round(DallaMan2007_data.Insulin_time(1:21),0)+360;
EGP_time = round(DallaMan2007_data.EGP_time(1:22),0)+360;
Ra_time = round(DallaMan2007_data.Ra_time(1:20),0)+360;
U_time = round(DallaMan2007_data.U_time(1:21),0)+360;
S_time = round(DallaMan2007_data.S_time(1:21),0)+360;
simG   = simG(Glucose_time);
simI   = simI(Insulin_time);
simEGP = simEGP(EGP_time);
simRa  = simRa(Ra_time);
simU   = simU(U_time);
simS   = simS(S_time);
CostG = (DallaMan2007_data.Glucose(1:21) - simG).^2./(DallaMan2007_data.GlucoseSEM(1:21).^2);
CostI = (DallaMan2007_data.Insulin(1:21) - simI).^2./(DallaMan2007_data.InsulinSEM(1:21).^2);
CostEGP = (DallaMan2007_data.EGP(1:22) - simEGP).^2./(DallaMan2007_data.EGPSEM(1:22).^2);
CostRa = (DallaMan2007_data.Ra(1:20) - simRa).^2./(DallaMan2007_data.RaSEM(1:20).^2);
CostU = (DallaMan2007_data.U(1:21) - simU).^2./(DallaMan2007_data.USEM(1:21).^2);
CostS = (DallaMan2007_data.S(1:21) - simS).^2./(DallaMan2007_data.SSEM(1:21).^2);

Cost_Model =  Cost_Model + sum(CostRa,'omitnan') + sum(CostG,'omitnan') + sum(CostI,'omitnan') + sum(CostEGP,'omitnan') + sum(CostU,'omitnan')+ sum(CostS,'omitnan'); %126 datapoints

%% Save cost
% if Cost_Model < chi2inv(0.95,471)
%     fprintf(fid, [sprintf('%.30f, ',params) sprintf('%.30f\n',Cost_Model)]);
% end
end

