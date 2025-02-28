% % =================================================================================================================
% Function to transform WindIO to DeSiO-input for modeling wind energy
% converter.
% 
% Author: Christian Hente
% Date: 05.05.2022
% updated: 27.10.2022
% =================================================================================================================
function param_study_vel()
% =================================================================================================================
clc; clear all; close all; format short;

my   = [3];
mx   = [10];
vinf = [300];
for i = 1:length(vinf)
    for j = 1:length(my)
        for k = 1:length(mx)
            fun_DeSiO_create_input_Wing(mx(k),my(j),vinf(i));
        end
    end
end

return