
%% Silfvergren et al. 2021
figure('Name', "Silfvergren et al. 2021", 'units', 'normalized', 'outerposition', [0 0 1 1])
 

% p1
subplot(2,1,1)
hold on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[2,4,6],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Glucose plasma (mM)'});
xlabel('Time [days]')
plot((Silfvergren2021_data.time-7709)/1440,Silfvergren2021_data.p1_glucoseCalibrated,'.--','MarkerSize',20)
line([0 0.02], [2 2],'Color','k','LineWidth',4);
line([2 2.02], [2 2],'Color','k','LineWidth',4);
xlim([0 2.2])
ylim([2 6])
legend('P1 Data','Meal');
hold off

% p2
subplot(2,1,2)
hold on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[2,4,6],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Glucose plasma (mM)'});
xlabel('Time [days]')
plot((Silfvergren2021_data.time(2:end-2)-7709)/1440,Silfvergren2021_data.p2_glucoseCalibrated(2:end-2),'.--','MarkerSize',20)
line([0 0.02], [2 2],'Color','k','LineWidth',4);
line([2 2.02], [2 2],'Color','k','LineWidth',4);
xlim([0 2.2])
ylim([2 6])
legend('P2 Data','Meal');
hold off
%%
figure('Name', "Silfvergren et al. 2021", 'units', 'normalized', 'outerposition', [0 0 1 1])
 

%unfed state
subplot(2,2,1)  
hold on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[2,4,6],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Glucose plasma (mM)'});
xlabel('Time [days]')
plot(Silfvergren2021_data.time_fasted_p1/60,Silfvergren2021_data.Value_fasted_p1,'.--','MarkerSize',20)
line([0 0.1], [2 2],'Color','k','LineWidth',4);
legend('P1 Data','Meal');
xlim([0 2])
ylim([2, 6])
hold off

%unfed state
subplot(2,2,2)  
hold on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[2,4,6],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Glucose plasma (mM)'});
xlabel('Time [days]')
plot(Silfvergren2021_data.time_fasted_p2/60,Silfvergren2021_data.Value_fasted_p2,'.--','MarkerSize',20)
line([0 0.1], [2 2],'Color','k','LineWidth',4);
legend('P1 Data','Meal');
legend('P2 Data','Meal');
xlim([0 2])
ylim([2, 6])
hold off

% fed state 
 subplot(2,2,3)  
hold on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[2,4,6],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Glucose plasma (mM)'});
xlabel('Time [days]')
plot(Silfvergren2021_data.time_fed_p1/60-128,Silfvergren2021_data.value_fed_p1,'.--','MarkerSize',20)
line([0 0.1], [2 2],'Color','k','LineWidth',4);
legend('P1 Data','Meal');
xlim([0 2])
ylim([2, 6])
hold off

% fed state 
 subplot(2,2,4)  
hold on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[2,4,6],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Glucose plasma (mM)'});
xlabel('Time [days]')
plot(Silfvergren2021_data.time_fed_p2/60-128,Silfvergren2021_data.value_fed_p2,'.--','MarkerSize',20)
line([0 0.1], [2 2],'Color','k','LineWidth',4);
legend('P2 Data','Meal');
xlim([0 2])
ylim([2, 6])
hold off


%% Qualitative constraints
figure('Name', "Qualitative constraints", 'units', 'normalized', 'outerposition', [0 0 1 1])
 

% EGP vs kidneys  
subplot(1,3,1)
hold on
a = gca;
set(a,'ytick',[0,20,40,60,80,100],'FontSize', 15,'fontname','Arial','FontSmoothing','on')
X = categorical({'EGP kidneys','EGP Liver'});
X = reordercats(X,{'EGP kidneys','EGP Liver'});
bar(X,[20,80]);
errorbar(X,[20,80],[10,10], " k . ",'MarkerSize',10,'LineWidth',2);
title('Mather et al. 2001','FontSize', 15,'FontSmoothing','on');
ylabel('Glucose production (%)','FontSize', 25,'FontSmoothing','on')
hold off
    
 % Fasting; U    
subplot(1,3,2)
hold on
a = gca;
set(a,'ytick',[0,20,40,60,80,100],'FontSize', 15,'fontname','Arial','FontSmoothing','on')
X = categorical({'fat','muscle','brain','liver'});
X = reordercats(X,{'fat','muscle','brain','liver'});
bar(X,[17.5,22.5, 40, 20]);
errorbar(X,[17.5,22.5, 40, 20],[5,5,5,5], " k . ",'MarkerSize',10,'LineWidth',2);
title('Gerich et al. 2010','FontSize', 15,'FontSmoothing','on');
ylabel('Organ Glucose uptake fasting (%)','FontSize', 25,'FontSmoothing','on')
ylim([0 60])
hold off

% Insulin degrad   
subplot(1,3,3)
hold on
a = gca;
set(a,'ytick',[0,20,40,60,80,100],'FontSize', 15,'fontname','Arial','FontSmoothing','on')
X = categorical({'Degradation Liver','Degradation rest of body'});
X = reordercats(X,{'Degradation Liver','Degradation rest of body'});
bar(X,[65,25]);
errorbar(X,[65,25],[10,10], " k . ",'MarkerSize',10,'LineWidth',2);
title('Najjar et al. 2019','FontSize', 15,'FontSmoothing','on');
ylabel('Insulin clearance (%)','FontSize', 25,'FontSmoothing','on')
hold off

%% Mottalib 2018
figure('Name', "Mottalib et al. 2018 ", 'units', 'normalized', 'outerposition', [0 0 1 1])
 

% HbA1c
subplot(2,3,1)
hold on
a = gca;
set(a,'xtick',[0,10,20],'ytick',[7,8,9],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'HbA1c (%)'});
xlabel('Time [weeks]')
errorbar(Mottalib2018_data.time/10080,Mottalib2018_data.GroupA_HbA1c,Mottalib2018_data.GroupA_HbA1cSEM, "b  -- ",'MarkerSize',15,'LineWidth',4);
errorbar(Mottalib2018_data.time/10080,Mottalib2018_data.GroupB_HbA1c,Mottalib2018_data.GroupB_HbA1cSEM, "r  -- ",'MarkerSize',15,'LineWidth',4);
errorbar(Mottalib2018_data.time/10080,Mottalib2018_data.GroupC_HbA1c,Mottalib2018_data.GroupC_HbA1cSEM, "k  -- ",'MarkerSize',15,'LineWidth',4);
ylim([7 9])
legend('Group A','Group B', 'Group C')
hold off

% BW
subplot(2,3,2)
hold on
a = gca;
set(a,'xtick',[0,10,20],'ytick',[90,100,110],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Body weight (kg)'});
xlabel('Time [weeks]')
errorbar(Mottalib2018_data.time/10080,Mottalib2018_data.GroupA_BW,Mottalib2018_data.GroupA_BWSEM, "b -- ",'MarkerSize',15,'LineWidth',4);
errorbar(Mottalib2018_data.time/10080,Mottalib2018_data.GroupB_BW,Mottalib2018_data.GroupB_BWSEM, "r -- ",'MarkerSize',15,'LineWidth',4);
errorbar(Mottalib2018_data.time/10080,Mottalib2018_data.GroupC_BW,Mottalib2018_data.GroupC_BWSEM, "k -- ",'MarkerSize',15,'LineWidth',4);
ylim([90 110])
hold off

% Fasting glucose
subplot(2,3,3)
hold on
a = gca;
set(a,'xtick',[0,10,20],'ytick',[8,9,10,11],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Fasting glucose (mM)'});
xlabel('Time [weeks]')
plot(Mottalib2018_data.timeStartEnd(1:2)/10080,Mottalib2018_data.GroupA_FastingGlucose(1:2), "b --o ",'MarkerSize',15,'LineWidth',4)
plot(Mottalib2018_data.timeStartEnd(1:2)/10080,Mottalib2018_data.GroupB_FastingGlucose(1:2), "r --o ",'MarkerSize',15,'LineWidth',4)
plot(Mottalib2018_data.timeStartEnd(1:2)/10080,Mottalib2018_data.GroupC_FastingGlucose(1:2), "k --o ",'MarkerSize',15,'LineWidth',4)
ylim([7.5 11])
hold off

% FastingI
subplot(2,3,4)
hold on
a = gca;
set(a,'xtick',[0,10,20],'ytick',[120,140,160],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Fasting Insulin (pM)'});
xlabel('Time [weeks]')
plot(Mottalib2018_data.timeStartEnd(1:2)/10080,Mottalib2018_data.GroupA_FastingI(1:2), "b --o ",'MarkerSize',15,'LineWidth',4)
plot(Mottalib2018_data.timeStartEnd(1:2)/10080,Mottalib2018_data.GroupB_FastingI(1:2), "r --o ",'MarkerSize',15,'LineWidth',4)
plot(Mottalib2018_data.timeStartEnd(1:2)/10080,Mottalib2018_data.GroupC_FastingI(1:2), "k --o ",'MarkerSize',15,'LineWidth',4)
ylim([115 160])
hold off

% BodyFat
subplot(2,3,5)
hold on
a = gca;
set(a,'xtick',[0,10,20],'ytick',[38,40,42,44],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Body fat (%)'});
xlabel('Time (weeks)')
plot(Mottalib2018_data.timeStartEnd(1:2)/10080,Mottalib2018_data.GroupA_BodyFat(1:2), "b --o ",'MarkerSize',15,'LineWidth',4)
plot(Mottalib2018_data.timeStartEnd(1:2)/10080,Mottalib2018_data.GroupB_BodyFat(1:2), "r --o ",'MarkerSize',15,'LineWidth',4)
plot(Mottalib2018_data.timeStartEnd(1:2)/10080,Mottalib2018_data.GroupC_BodyFat(1:2), "k --o ",'MarkerSize',15,'LineWidth',4)
ylim([38 44])
hold off

% Daily Intake food
subplot(2,3,6)
hold on
a = gca;
set(a,'ytick',[0,1000,2000],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
X = categorical({'Baseline A','During study A','Baseline B' , 'During study B','Baseline C','During study C'});
X = reordercats(X,{'Baseline A','During study A','Baseline B' , 'During study B','Baseline C','During study C'});
bar(X,[Mottalib2018_data.GroupA_DailyKcal(1),Mottalib2018_data.GroupA_DailyKcal(2),-10, -10, -10,-10]);
bar(X,[-10,-10,Mottalib2018_data.GroupB_DailyKcal(1),Mottalib2018_data.GroupB_DailyKcal(2), -10,-10]);
bar(X,[-10,-10,-10, -10, Mottalib2018_data.GroupC_DailyKcal(1),Mottalib2018_data.GroupC_DailyKcal(2)]);
ylim([0 2500])
ylabel('Daily Intake food (kcal/day)','FontSize', 25,'FontSmoothing','on')
hold off

%% Lindstrom2012_data
figure('Name', "Lindstrom et al. 2012", 'units', 'normalized', 'outerposition', [0 0 1 1])
 

% Weight
subplot(2,2,1)
hold on
a = gca;
set(a,'xtick',[0,1,2,3,4],'ytick',[60,70,80],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Weight(kg)'});
xlabel('Time (weeks)')
errorbar(Lindstrom2012_data.time/10080,Lindstrom2012_data.weight,Lindstrom2012_data.weightSEM, "-. ",'MarkerSize',15,'LineWidth',4);
xlim([-0.5 4.5])
ylim([60 80])
hold off

% insulin
subplot(2,2,2)
hold on
a = gca;
set(a,'xtick',[0,1,2,3,4],'ytick',[0,40,80],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Fasting Insulin (pM)'});
xlabel('Time (weeks)')
errorbar(Lindstrom2012_data.insulinTime(1:3)/10080,Lindstrom2012_data.insulin(1:3),Lindstrom2012_data.insulinSEM(1:3), "-. ",'MarkerSize',15,'LineWidth',4);
xlim([-0.5 4.5])
ylim([0 80])
hold off

% bodyfat
subplot(2,2,3)
hold on
a = gca;
set(a,'xtick',[0,1,2,3,4],'ytick',[15,20,25],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Bodyfat(%)'});
xlabel('Time (weeks)')
errorbar(Lindstrom2012_data.bodyfatTime(1:2)/10080,Lindstrom2012_data.bodyfat(1:2),Lindstrom2012_data.bodyfatSEM(1:2), "-. ",'MarkerSize',15,'LineWidth',4);
xlim([-0.5 4.5])
ylim([15 26])
hold off

% Dailyintake  
subplot(2,2,4)
a = gca;
set(a,'ytick',[0,2000,4000,6000],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
hold on
X = categorical({'Baseline','During study'});
X = reordercats(X,{'Baseline','During study'});
bar(X,Lindstrom2012_data.Dailyintake(1:2));
ylabel('Daily Intake food (kcal/day)','FontSize', 25,'FontSmoothing','on')
hold off

%% Siero2008_data
figure('Name', "Siero et al. 2008", 'units', 'normalized', 'outerposition', [0 0 1 1])
 
% Weight
subplot(2,3,1)
hold on
a = gca;
set(a,'xtick',[0,4,8,12,16],'ytick',[50,70,90],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Weight(kg)'});
xlabel('Time (weeks)')
plot((Siero2008_data.timestart)/10080,Siero2008_data.p1BW, "--o ",'MarkerSize',15,'LineWidth',4)
plot((Siero2008_data.timestart)/10080,Siero2008_data.p2BW, "--o ",'MarkerSize',15,'LineWidth',4)
plot((Siero2008_data.timestart)/10080,Siero2008_data.p3BW, "--o ",'MarkerSize',15,'LineWidth',4)
plot((Siero2008_data.timestart)/10080,Siero2008_data.p4BW, "--o ",'MarkerSize',15,'LineWidth',4)
plot((Siero2008_data.timestart)/10080,Siero2008_data.p5BW, "--o ",'MarkerSize',15,'LineWidth',4)
legend('P1','P2','P3','P4','P5')
xlim([0 18])
hold off

% bodyfat
subplot(2,3,2)
hold on
a = gca;
set(a,'xtick',[0,4,8,12,16],'ytick',[65,70,75,80],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Total body fat (kg)'});
xlabel('Time (weeks)')
errorbar(Siero2008_data.timestart/10080,Siero2008_data.BW,Siero2008_data.BWSEM, "-. ",'MarkerSize',15,'LineWidth',4);
xlim([0 18])
hold off

% BMR
subplot(2,3,3)
hold on
a = gca;
set(a,'xtick',[0,4,8,12,16],'ytick',[1500,2000,2500,3000,3500],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Basal Metabolic Rate (kcal/day)'});
xlabel('Time (weeks)')
plot((Siero2008_data.timestart)/10080,Siero2008_data.p1BMR, "--o ",'MarkerSize',15,'LineWidth',4)
plot((Siero2008_data.timestart)/10080,Siero2008_data.p2BMR, "--o ",'MarkerSize',15,'LineWidth',4)
plot((Siero2008_data.timestart)/10080,Siero2008_data.p3BMR, "--o ",'MarkerSize',15,'LineWidth',4)
plot((Siero2008_data.timestart)/10080,Siero2008_data.p4BMR, "--o ",'MarkerSize',15,'LineWidth',4)
plot((Siero2008_data.timestart)/10080,Siero2008_data.p5BMR, "--o ",'MarkerSize',15,'LineWidth',4)
legend('P1','P2','P3','P4','P5')
xlim([0 18])
hold off

% bodyfat
subplot(2,3,4)
hold on
a = gca;
set(a,'xtick',[0,4,8,12,16],'ytick',[10,15,20],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Total Fat weight (kg)'});
xlabel('Time (weeks)')
errorbar(Siero2008_data.timestartFat/10080,Siero2008_data.BodyFat,Siero2008_data.BodyFatSEM, "-. ",'MarkerSize',15,'LineWidth',4);
xlim([0 18])
hold off

% kcal
subplot(2,3,5)
hold on
a = gca;
set(a,'xtick',[0,4,8,12,16],'ytick',[1000,2500,4000,5500],'FontSize', 25,'fontname','Arial','FontSmoothing','on')
ylabel({'Calorie intake (kcal/day)'});
xlabel('Time (weeks)')
plot((Siero2008_data.timestart)/10080,Siero2008_data.Dailyintake, "--o ",'MarkerSize',15,'LineWidth',4)
xlim([0 18])
hold off

% Person Specific coverancies
subplot(2,3,6)
hold on
ylabel({'Age(years)'});
xlabel('Height (m)')
title('Person Specific coverancies')
for i=1:5
plot(Siero2008_data.height(i),Siero2008_data.age(i), "x ",'MarkerSize',25,'LineWidth',5)
end
hold off

%% Silfvergren2022

figure('Name', "Silfvergren et al. 2022 - Weight", 'units', 'normalized', 'outerposition', [0 0 1 1])

hold on
set(gca,'ytick',[75,80,85,90,95],'xtick',[0,4,8,12,16],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Weight'; '(kg)'});
xlabel({'Time (weeks)'});
plot((Silfvergren2022Weight_data.Time)/10080,Silfvergren2022Weight_data.Weight,'Color', BlueColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize+40,'LineWidth',LineWidthValue)
ylim([75 95])
xlim([0 16])
hold off
grid on

%%
figure('Name', "Silfvergren et al. 2022 - food", 'units', 'normalized', 'outerposition', [0 0 1 1])

subplot(1,2,1)
hold on
grid on
set(gca,'ytick',[0,10,20,30,40,50,60,70,80,90,100],'xtick',[0,4,8,12,16],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Lean body mass';'(%)'});
xlabel({'Time (weeks)'});
plot((Silfvergren2022Weight_data.Time)/10080,Silfvergren2022Weight_data.LeanmassProcentage*100,'Color', BlueColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
ylim([60 100])
xlim([0 16])
hold off


subplot(1,2,2)
hold on
grid on
set(gca,'ytick',[0,5,10,15,20,30,40,50,60,70,80,90,100],'xtick',[0,4,8,12,16],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Total body fat';'(%)'});
xlabel({'Time (weeks)'});
plot((Silfvergren2022Weight_data.Time)/10080,Silfvergren2022Weight_data.FatProcentage*100,'Color', BlueColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
ylim([0 20])
xlim([0 16])
hold off

%%
figure('Name', "Silfvergren et al. 2022 - food", 'units', 'normalized', 'outerposition', [0 0 1 1])


CaloriesFromFood = [(Silfvergren2022Food_data.Lipids(1:16)*9)'; (Silfvergren2022Food_data.Starch(1:16)*4)'; (Silfvergren2022Food_data.Sugar(1:16)*4)'; (Silfvergren2022Food_data.Fiber(1:16)*4)'; (Silfvergren2022Food_data.Protein(1:16)*4)'];

hold on
set(gca,'ytick',[0,1000,2000,3000],'xtick',[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Calorie intake' ; '(kcal/day)'});
xlabel('Week')
ba = bar(Silfvergren2022Food_data.Week(1:16), CaloriesFromFood,'stacked', 'FaceColor','flat');
ba(1).CData = YellowColor;
ba(2).CData = RedColor;
ba(3).CData = OrangeColor;
ba(4).CData = [0.1 0.1 0.1];
ba(5).CData = GreenColor;
legend('Fat', 'Starch', 'Sugar', 'Fiber', 'Protein','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize)
xlim([0 17])
ylim([0 3100])
hold off
grid on

figure('Name', "Silfvergren et al. 2022 - food", 'units', 'normalized', 'outerposition', [0 0 1 1])

X = categorical({'Energy expendature'});
X = reordercats(X,{'Energy expendature'});

hold on
set(gca,'ytick',[0,200,400,600,800],'xtick',[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Physical activity' ; '(kcal/day)'});
xlabel('Week')
bar(Silfvergren2022Food_data.Week(1:16),Silfvergren2022Food_data.EnergyExpendature(1:16),'FaceColor',BlueColor)
ylim([0 800])
hold off
grid on


%% Silfvergren2022 Glucose

figure('Name', "Silfvergren et al. 2022", 'units', 'normalized', 'outerposition', [0 0 1 1])
 
hold on
grid on
set(gca,'ytick',[0,2,4,6,8,10],'xtick',[0,4,8,12,16],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Glucose';'(mM)'});
xlabel({'Time (weeks)'});
plot((Silfvergren2022CGM_data.time_minutesGlucose)/10080,Silfvergren2022CGM_data.Glucose,'Color', BlueColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize,'LineWidth',LineWidthValue)
plot((Silfvergren2022GlucoseBloodSensor_data.glucoseTime)/10080,Silfvergren2022GlucoseBloodSensor_data.bloodvalue,'Color', RedColor,'LineStyle', '--', 'Marker','.','MarkerSize',MarkerSize-2,'LineWidth',LineWidthValue-2)
legend('Continous glucose sensor','Basal plasma glucose','Location',legendLocation,'Orientation',legendOrientation,'FontSize',LegendSize+3)
ylim([0 10])
xlim([0 16])
hold off
%%
figure('Name', "Silfvergren et al. 2022", 'units', 'normalized', 'outerposition', [0 0 1 1])

hold on
grid on
set(gca,'ytick',[0,0.5,1,1.5,2],'xtick',[0,4,8,12,16],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Ketones'; '(mM)'});
xlabel({'Time (weeks)'});
plot((Silfvergren2022Ketones_data.Time_minutesKetones)/10080,Silfvergren2022Ketones_data.Ketones,'Color', BlueColor,'LineStyle', 'n', 'Marker','.','MarkerSize',MarkerSize+20,'LineWidth',LineWidthValue)
ylim([0 2])
xlim([0 16])
hold off
grid on


%% Silfvergren2022 Overview
figure('Name', "Silfvergren et al. 2022 - food", 'units', 'normalized', 'outerposition', [0 0 1 1])

CaloriesFromFood = [Silfvergren2022Food_data.LipidsTotal(1), Silfvergren2022Food_data.StarchTotal(1), Silfvergren2022Food_data.SugarTotal(1), Silfvergren2022Food_data.FiberTotal(1), Silfvergren2022Food_data.ProteinTotal(1)];
CaloriesFromFoodSEM = [Silfvergren2022Food_data.LipidsTotalSEM(1), Silfvergren2022Food_data.StarchTotalSEM(1), Silfvergren2022Food_data.SugarTotalSEM(1), Silfvergren2022Food_data.FiberTotalSEM(1), Silfvergren2022Food_data.ProteinTotalSEM(1)];

X = categorical({'Fat','Starch','Sugar','Fiber','Protein'});
X = reordercats(X,{'Fat','Starch','Sugar','Fiber','Protein'});

subplot(1,1,1)
hold on
set(gca,'ytick',[0,50,100,150],'xtick',[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Food intake'; '(g/day)'});
bar(X(1),CaloriesFromFood(1),'FaceColor',YellowColor);
bar(X(2),CaloriesFromFood(2),'FaceColor',RedColor);
bar(X(3),CaloriesFromFood(3),'FaceColor',OrangeColor);
bar(X(4),CaloriesFromFood(4),'FaceColor',[0.3 0.3 0.3]);
bar(X(5),CaloriesFromFood(5),'FaceColor',GreenColor);
errorbar(X,CaloriesFromFood,CaloriesFromFoodSEM,'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',4,'LineWidth',4)
ylim([0 150])
hold off
grid on

figure('Name', "Silfvergren et al. 2022 - food", 'units', 'normalized', 'outerposition', [0 0 1 1])


X = categorical({'Food intake', 'Physical activity'});
X = reordercats(X,{'Food intake', 'Physical activity'});

subplot(1,2,1)
hold on
set(gca,'ytick',[0,800,1600,2400],'xtick',[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Calories' ; '(kcal/day)'});
bar(X,[Silfvergren2022Food_data.MealTotal(1),Silfvergren2022Food_data.EnergyExpendature1(1)])
errorbar(X(1),Silfvergren2022Food_data.MealTotal(1),Silfvergren2022Food_data.MealTotalSEM(1),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',4,'LineWidth',4)
errorbar(X(2),Silfvergren2022Food_data.EnergyExpendature1(1),Silfvergren2022Food_data.EnergyExpendatureSEM(1),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',4,'LineWidth',4)
ylim([0 2400])
hold off
grid on

X = categorical({'Total', 'Lean mass' , 'Fat mass'});
X = reordercats(X,{'Total', 'Lean mass' , 'Fat mass'});

Weightloss = [Silfvergren2022Weight_data.AverageWeightchangeDay(1); Silfvergren2022Weight_data.AverageLeanmasschangeDay(1); Silfvergren2022Weight_data.AverageFatchangeDay(1)];

subplot(1,2,2)
hold on
set(gca,'ytick',[-0.1,0,0.1],'xtick',[],'FontSize', FontSize,'fontname','Arial','FontSmoothing','on')
ylabel({'Mass loss'; '(kg)'});
bar(X(1),Weightloss(1),'FaceColor',BlueColor);
bar(X(2),Weightloss(2),'FaceColor',MagentaColor);
bar(X(3),Weightloss(3),'FaceColor',YellowColor);
errorbar(X(1),Silfvergren2022Weight_data.AverageWeightchangeDay(1),Silfvergren2022Weight_data.AverageWeightchangeDaySEM(1),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',4,'LineWidth',4)
errorbar(X(2),Silfvergren2022Weight_data.AverageLeanmasschangeDay(1),Silfvergren2022Weight_data.AverageLeanmasschangeDaySEM(1),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',4,'LineWidth',4)
errorbar(X(3),Silfvergren2022Weight_data.AverageFatchangeDay(1),Silfvergren2022Weight_data.AverageFatchangeDay(1),'Color', BlackColor,'LineStyle', 'n', 'Marker','.','MarkerSize',4,'LineWidth',4)
ylim([-0.2 0.1])
hold off
grid on