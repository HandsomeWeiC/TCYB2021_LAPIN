%% input
tic
clc 
clear all
%% data input
addpath('fun')

X = cell2mat(struct2cell(load('X_high_dimen.mat'))); % ���ݼ�
truth=zeros(300,1);
truth(1:100)=1;
truth(101:200)=2;
truth(201:300)=3;
k = max(truth);  %�������
c = 10;  %ϵ������һ������ֻ��c������ֵ  ��������
lambda = 100;  %���� 
gamma = 0.05;  % �ֵ�ѡȡ����
len = length(gamma);

[y2,idx,Z,A] = LAPIN(X,k,c,lambda,gamma);
result1 = ClusteringMeasure(y2, truth)
result_idx1 = ClusteringMeasure(idx, truth)


% save result_toy_LAPIN_0812_1.mat Result Result_idx