%%%%%this code is used for simulating cells under different parameters and
%%%%%calculating pulse coincidence rate based on Model.m
addpath('path of Model.m'); %add the path of Model.m
result_folder_0='path for storing results'; %path for storing results

cellnum_0=3; %cell number
cellnumforplot=3; %how many cells to plot

mkdir(result_folder_0);
cd(result_folder_0);

a2_0=1; %basal rate of nuclear importing of Msn2 without regulation by PKA
a4_0=1200; %basal rate of nuclear importing of Msn4 without regulation by PKA and X
n2_0=4; %Hill coefficient for PKA regulating Msn2
n4_0=4; %Hill coefficient for PKA regulating Msn4
nx_0=4; %Hill coefficient for X regulating Msn4
b2_0=2; %the rate of nuclear exporting of Msn2
b4_0=2; %the rate of nuclear exporting of Msn4

M2tot_0=200; %total concentration of Msn2
M20_0=0; %initial concentration of Msn2
M4tot_0=200; %total concentration of Msn4
M40_0=0; %initial concentration of Msn4

tau_0=0.001; %step size
timelimit_0=300; %duration of simulation
PKAdecrease_0=4; %level of weak down-regulation of PKA
Xdecrease_0=4; %level of weak down-regulation of X

PKAOscillationPercentage_0=25; %frequency of periodic fluctuations of PKA and X

x_0=6; %X level
threshold_PKA_0=0.25; %proportion of up-regulation of PKA
threshold_X_0=0.25; %proportion of up-regulation of X

PKA_vector=6:0.25:10; %PKA level
k4_k2_vector=0.2:0.1:1; %Kd(Msn4)/Kd(Msn2)

parameter_matrix=[];
for PKA_i=1:length(PKA_vector)
    for k4_k2_i=1:length(k4_k2_vector)
        parameter_matrix=[parameter_matrix;PKA_vector(PKA_i),k4_k2_vector(k4_k2_i)];
    end
end

size_parameter_matrix=size(parameter_matrix);

%%%%%parallel running if cell number is large
% ClusterNumber=84;
% myCluster = parcluster('local');
% myCluster.NumWorkers = ClusterNumber;
% saveProfile(myCluster);
% myCluster ;
% parpool(myCluster,ClusterNumber);

% parfor ii=1:size_parameter_matrix(1)

%%%%%simulate under each parameter and store results in respective folder
for ii=1:size_parameter_matrix(1)
    PKA_0=parameter_matrix(ii,1);
    k4_k2_0=parameter_matrix(ii,2);
    
    result_folder_ii=[result_folder_0,'/','PKA ',num2str(PKA_0),' k4_k2 ',num2str(k4_k2_0)];
    
    if ~exist([result_folder_ii '/data_result.csv'],'file')
        [Msn2TotalFreq,Msn2TotalFreqSe,Msn4TotalFreq,Msn4TotalFreqSe,...
            Msn2OnlyFreq,Msn2OnlyFreqSe,Msn4OnlyFreq,Msn4OnlyFreqSe,OverlapFreq,OverlapFreqSe,...
            OverlapFraction,OverlapFractionErrorbar]=Model(a2_0,a4_0,n2_0,n4_0,nx_0,...
            b2_0,b4_0,...
            k4_k2_0,...
            PKA_0,PKAdecrease_0,x_0,Xdecrease_0,...
            M2tot_0,M4tot_0,M20_0,M40_0,tau_0,timelimit_0,...
            PKAOscillationPercentage_0,...
            cellnum_0,result_folder_ii,cellnumforplot,threshold_PKA_0,threshold_X_0);
        
        str={'x','PKA','k4_k2','threshold_PKA','threshold_X',...
            'Msn2TotalFreq','Msn2TotalFreqSe','Msn4TotalFreq','Msn4TotalFreqSe',...
            'Msn2OnlyFreq','Msn2OnlyFreqSe','Msn4OnlyFreq','Msn4OnlyFreqSe','OverlapFreq','OverlapFreqSe',...
            'OverlapFraction','OverlapFractionErrorbar'};
        data={x_0,PKA_0,k4_k2_0,threshold_PKA_0,threshold_X_0,...
            Msn2TotalFreq,Msn2TotalFreqSe,Msn4TotalFreq,Msn4TotalFreqSe,...
            Msn2OnlyFreq,Msn2OnlyFreqSe,Msn4OnlyFreq,Msn4OnlyFreqSe,OverlapFreq,OverlapFreqSe,...
            OverlapFraction,OverlapFractionErrorbar};
        data_merge=table();
        data_merge{1,:}=str;
        data_merge{2,:}=data;
        writetable(data_merge,[result_folder_ii,'/result.csv']);
        
    end
end

