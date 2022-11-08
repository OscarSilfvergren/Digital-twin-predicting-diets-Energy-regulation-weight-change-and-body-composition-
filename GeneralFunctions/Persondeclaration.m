function [param, initialvalues] = Persondeclaration(Height,male_boolean,female_boolean,LeanBodyWeight_start,Bodyfat_start,initialvalues,sNames,param,pNames)

param(ismember(pNames,"Height"))          = Height;
param(ismember(pNames,"male_boolean"))    = male_boolean;
param(ismember(pNames,"female_boolean"))  = female_boolean;

initialvalues(ismember(sNames,"LeanBodyWeight"))       = LeanBodyWeight_start;
initialvalues(ismember(sNames,"TGA_Adipocyte"))         = Bodyfat_start*1000;

end

