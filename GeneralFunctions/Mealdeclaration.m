function [param] = Mealdeclaration(CarbohydratesAmount,ProteinAmount, LipidsAmount,MealLength,param,pNames)

param(ismember(pNames,"carb_flow"))    = CarbohydratesAmount/MealLength;
param(ismember(pNames,"protein_flow")) = ProteinAmount/MealLength;
param(ismember(pNames,"lipids_flow"))  = LipidsAmount/MealLength;

end

