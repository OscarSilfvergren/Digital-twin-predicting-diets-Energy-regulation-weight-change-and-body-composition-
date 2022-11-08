clear
clc
close all
format longG

addpath('Costfunctions')
addpath('Data')
addpath('IQMtools')
addpath('Models')
addpath('Optimisation')
addpath('Options')
addpath('Parameters')
addpath('Plotfunctions')
addpath('GeneralFunctions')

% Data
load('ParametersInformation');
load('Bray2012_data');
load('Cortisol_data');
load('DallaMan2007_data');
load('Firth1986_data');
load('Krssak2004_dataDiabetes');
load('Krssak2004Healthy_data');
load('Lerche2009_data');
load('Lindstrom2012_data');
load('Magnusson1992_data');
load('Mottalib2018_data');
load('Rothman1991_data.mat')
load('Siero2008_data');
load('Silfvergren2021_data');
load('Taylor1996_data');
load('Silfvergren2022GlucoseBloodSensor_data');
load('Silfvergren2022Food_data');
load('Silfvergren2022Weight_data');
load('Silfvergren2022Ketones_data');
load('Silfvergren2022_meals');
load('Silfvergren2022CGM_data');

%options2 = optimoptions(@particleswarm,'PlotFcn','pswplotbestf');
optionsParalell = optimoptions(@particleswarm,'InitialSwarmSpan',2000,'PlotFcn','pswplotbestf','MaxStallIterations',200,'FunctionTolerance',1.0e-6,'UseParallel', true,'UseVectorized', false);
options = optimoptions(@particleswarm,'InitialSwarmSpan',2000,'PlotFcn','pswplotbestf','MaxStallIterations',200,'FunctionTolerance',1.0e-6);
i = 1;

if ~exist('IQMsimulate','file')
    fprintf('\n\nFor these scripts to work, the IQM tools toolbox have to be compiled.\n')
    disp('To do this, a valid C-compiler is necessary')
    disp('If the compilation of the toolbox does not, it is likely that a valid C-compiler is missing.')
    choice = input('press enter to continue and compile the toolbox.');
    run('./IQMtools/installIQMtoolsInitial.m')
end

BlueColor    = [0, 0.5, 0.9];
OrangeColor  = [0.8500, 0.3250, 0.0980];
BlackColor   = [0 0 0]; 
GreyColor    = [0.5 0.5 0.5]; 
RedColor     = [0.6350, 0.0780, 0.1840];
GreenColor   = [0.4660, 0.6740, 0.1880];
YellowColor  = [0.9290, 0.6940, 0.1250];
MagentaColor = [0.4940 0.1840 0.5560];
%%
marker      = 'o';
MarkerSize  = 50;
MarkerSizeErrorbars = 80;
LineWidthValue   = 5;
LineWidthMeal    = 20;
legendLocation = 'best';
legendOrientation = 'horizontal';
LegendSize = 30;
FontSize = 45;
CapSizeValue = 20;
