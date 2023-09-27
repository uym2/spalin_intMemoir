clearvars
close all
clc

f_name = 'p2c.txt';
fid = fopen("sample_names.txt");
names = textscan(fid, '%s');
name_list = string(zeros(2064));

y =plot2tree(f_name);

%slope = round(compute_slope(f_name));
%theta = round(compute_theta(f_name));

%writematrix([slope theta],"slope_and_theta.csv")
%histogram(theta)


% s9
%histogram(slope(1926:2064),max(slope)-min(slope))
% s10
%histogram(slope(1:88),max(slope)-min(slope))
% s8
%histogram(slope(1799:1925),max(slope)-min(slope))

%histogram(slope,max(slope)-min(slope))

% slope = slope(~isinf(slope) & ~isnan(slope));
% var(slope(1:27))
% var(slope(28:65))
% var(slope(66:115))
% var(slope)