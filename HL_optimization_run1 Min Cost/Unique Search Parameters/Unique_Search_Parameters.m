clear
clc
close all

Parameter = load('Layer_EffTEffB_Cost.txt');

unique_parameter = unique(Parameter,'rows');

writematrix(unique_parameter,'unique_parameter.txt','Delimiter','\t');