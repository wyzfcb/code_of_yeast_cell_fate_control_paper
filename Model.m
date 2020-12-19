function [Msn2TotalFreq,Msn2TotalFreqSe,Msn4TotalFreq,Msn4TotalFreqSe,...
    Msn2OnlyFreq,Msn2OnlyFreqSe,Msn4OnlyFreq,Msn4OnlyFreqSe,OverlapFreq,OverlapFreqSe,...
    OverlapFraction,OverlapFractionErrorbar]=Model(a2_0,a4_0,n2_0,n4_0,nx_0,...
    b2_0,b4_0,...
    k4_k2_0,...
    maxPKA_0,PKAdecrease_0,maxX_0,Xdecrease_0,...
    M2tot_0,M4tot_0,M20_0,M40_0,tau_0,timelimit_0,...
    PKAOscillationPercentage_0,...
    cellnum_0,result_folder,cellnumforplot,threshold_PKA,threshold_X)

%This function is used for stochastic simulations of Msn2 and Msn4
%dynamics based on 'tau-leap' method

%If you have any questions, please contact with Yan WU (wy-fcb@pku.edu.cn)

mkdir(result_folder);
cd(result_folder);

a2=a2_0; %basal rate of nuclear importing of Msn2 without regulation by PKA
a4=a4_0; %basal rate of nuclear importing of Msn4 without regulation by PKA and X

n2=n2_0; %Hill coefficient for PKA regulating Msn2
n4=n4_0; %Hill coefficient for PKA regulating Msn4
nx=nx_0; %Hill coefficient for X regulating Msn4

b2=b2_0;b4=b4_0; %the rate of nuclear exporting of Msn2 (Msn4)
k4_k2=k4_k2_0; %Kd(Msn4)/Kd(Msn2)

M2tot=M2tot_0; %total concentration of Msn2
M4tot=M4tot_0; %total concentration of Msn4
M20=M20_0; %initial concentration of Msn2
M40=M40_0; %initial concentration of Msn4
tau=tau_0; %step size
timelimit=timelimit_0; %duration of simulation
 
maxPKA=maxPKA_0; %PKA level
PKAdecrease=PKAdecrease_0; %level of weak down-regulation of PKA
maxX=maxX_0; %X level
Xdecrease=Xdecrease_0; %level of weak down-regulation of X

PKAOscillationPercentage=PKAOscillationPercentage_0; %frequency of periodic fluctuations of PKA and X
cellnum=cellnum_0; %cell number for simulation

minPKA=0; %minimum level after strong down-regulation of PKA
minPKA_2=maxPKA-PKAdecrease; %minimum level after weak down-regulation of PKA
minX=0; %minimum level after strong down-regulation of X
minX_2=maxX-Xdecrease; %minimum level after weak down-regulation of X

tdur = 0:tau:10*timelimit;

for cellindex=1:cellnum
    %%%%%simulate trajectories of PKA and X
    y0 = (square(tdur,PKAOscillationPercentage)+1)/2;
    PKA_0=[];X_0=[];
    for i=1:length(y0)
        if i==1 && y0(i)==0
            PKA_0(i)=maxPKA;
            X_0(i)=maxX;
        elseif i==1 && y0(i)==1
            if rand(1,1)>threshold_PKA
                PKA_0(i)=maxPKA*(1+rand(1,1)/10);
            else
                if rand(1,1)>0.5
                    rand_i1=rand(1,1);
                    PKA_0(i)=maxPKA*rand_i1+minPKA*(1-rand_i1);
                else
                    rand_i1=rand(1,1);
                    PKA_0(i)=maxPKA*rand_i1+minPKA_2*(1-rand_i1);
                end
            end
            
            if rand(1,1)>threshold_X
                X_0(i)=maxX*(1+rand(1,1)/10);
            else
                if rand(1,1)>0.5
                    rand_i1=rand(1,1);
                    X_0(i)=maxX*rand_i1+minX*(1-rand_i1);
                else
                    rand_i1=rand(1,1);
                    X_0(i)=maxX*rand_i1+minX_2*(1-rand_i1);
                end
            end
            
        elseif i>1 && y0(i)>y0(i-1)
            if rand(1,1)>threshold_PKA
                PKA_0(i)=maxPKA*(1+rand(1,1)/10);
            else
                if rand(1,1)>0.5
                    rand_i1=rand(1,1);
                    PKA_0(i)=maxPKA*rand_i1+minPKA*(1-rand_i1);
                 else
                    rand_i1=rand(1,1);
                    PKA_0(i)=maxPKA*rand_i1+minPKA_2*(1-rand_i1);
                 end
            end
            
            if rand(1,1)>threshold_X
                X_0(i)=maxX*(1+rand(1,1)/10);
            else
                if rand(1,1)>0.5
                    rand_i1=rand(1,1);
                    X_0(i)=maxX*rand_i1+minX*(1-rand_i1);
                else
                    rand_i1=rand(1,1);
                    X_0(i)=maxX*rand_i1+minX_2*(1-rand_i1);
                end
            end
        elseif i>1 && y0(i)==y0(i-1)
            PKA_0(i)= PKA_0(i-1);
            X_0(i)= X_0(i-1);

        elseif i>1 && y0(i)<y0(i-1)
            PKA_0(i)= maxPKA;
            X_0(i)= maxX;
        end
    end
    
    index_vector=randperm(numel(0:tau:6*timelimit));
    
    PKA_trajectory=PKA_0(index_vector(1):(index_vector(1)+numel(0:tau:timelimit)-1));
    PKA=  PKA_trajectory(1);
    
    X_trajectory=X_0(index_vector(2):(index_vector(2)+numel(0:tau:timelimit)-1));
    X=  X_trajectory(1);
    
    %%%%%simulate trajectories of Msn2 and Msn4 from trajectories of
    %%%%%PKA and X based on 'tau-leap' method
    x=zeros(1,5);
    x(1,1)=PKA;
    x(1,2)=M20;
    x(1,3)=M40;
    x(1,4)=X;
    x(1,5)=0;
    
    counter=1;
    t=0;
    tcounter=0;
    
    tindex=1;
    while t < timelimit
              
        if counter==1
            counter2=2;
        elseif counter==2
            counter2=1;
        end
              
        a(1)=a2*(M2tot-x(counter,2))/(1+(x(counter,1)*k4_k2)^n2);
        a(2)=a4*(M4tot-x(counter,3))/((1+(x(counter,1))^n4)*(1+(x(counter,4))^nx));
          
        a(3)=b2*x(counter,2);
        a(4)=b4*x(counter,3);
            
        delta_a=[];
        for i=1:length(a)
            delta_a(i)=random('Poisson',a(i)*tau);
        end
        
        x(counter2,:)=x(counter,:);
        
        counter=counter2;
        x(counter2,1)= PKA_trajectory(tindex);
        x(counter2,2)=max(x(counter,2)+delta_a(1)-delta_a(3),0);
        x(counter2,3)=max(x(counter,3)+delta_a(2)-delta_a(4),0);
        x(counter2,4)=X_trajectory(tindex);        
        
        t=t+tau;
        x(counter2,5)=t;
        
        tindex=tindex+1;
        
        if (t-0.5*tcounter)>=0
            tcounter=tcounter+1;
            result(tcounter,:)=x(counter2,:); 
        end
        
    end
    
    time=result(:,5);
    PKAn=result(:,1);
    Msn2n=result(:,2);
    Msn4n=result(:,3);
    Xn=result(:,4);
    
    finalResults{cellindex}.time=time;
    finalResults{cellindex}.PKA=PKAn; %trajectory of PKA level
    finalResults{cellindex}.nuclearMsn2=Msn2n; %trajectory of nuclear localization level of Msn2
    finalResults{cellindex}.nuclearMsn4=Msn4n; %trajectory of nuclear localization level of Msn4
    finalResults{cellindex}.X=Xn; %trajectory of X level

    %     finalResults{cellindex}.b2=b2;
    %     finalResults{cellindex}.b4=b4;
    %     finalResults{cellindex}.k21=k21;
    %     finalResults{cellindex}.k41=k41;
    %     finalResults{cellindex}.n2=n2;
    %     finalResults{cellindex}.n4=n4;
    %
    %      finalResults{cellindex}.a2=a2;
    %     finalResults{cellindex}.a4=a4;
    %
    %     finalResults{cellindex}.M2tot=M2tot;
    %     finalResults{cellindex}.M20=M20;
    %     finalResults{cellindex}.M4tot=M4tot;
    %     finalResults{cellindex}.M40=M40;
    %
    %     finalResults{cellindex}.tau=tau;
    %     finalResults{cellindex}.timelimit=timelimit;
    %     finalResults{cellindex}.maxPKA=maxPKA;
    %     finalResults{cellindex}.minPKA=minPKA;
    %
    %     finalResults{cellindex}.tdur=tdur;
    %     finalResults{cellindex}.PKAOscillationPercentage=PKAOscillationPercentage;
    
    %     save([result_folder,'\','finalResults','.mat'],'finalResults');
end

%%%%%identify pulses of Msn2 and Msn4 for each cell
cellnum=length(finalResults);
for cellindex=1:cellnum
    data_cellindex=finalResults{cellindex};
    framenum_index=1:length(data_cellindex.time);
    Msn2_cellindex=data_cellindex.nuclearMsn2;
    Msn4_cellindex=data_cellindex.nuclearMsn4;
    
    cellTrace{1} = Msn2_cellindex;
    cellTrace{2} = Msn4_cellindex;
    thre_for_idxpulse_Msn2=25+median(Msn2_cellindex); %threshold for identification of Msn2 pulse
    thre_for_idxpulse_Msn4=25+median(Msn4_cellindex); %threshold for identification of Msn4 pulse
    slopeThreshold=0.005;
    %%%%%identify and measure peaks.
    [peakstemp{1:2}] =deal([]);
    [peaks2{1:2}] =deal([]);
    for i = 1 %find and measure peaks of Msn2
        ampThresh = thre_for_idxpulse_Msn2;
        peakstemp{i} = findpeaks_v20190613(framenum_index,cellTrace{i},slopeThreshold,ampThresh); %find peaks
        x = framenum_index; y = cellTrace{i};
        pcount = 1;
        for idxpk = 1:size(peakstemp{i},1) %measure peaks and eliminate shoulder peaks.
            if ~isshoulderpeak(peakstemp{i},idxpk,x,y,ampThresh) && peakstemp{i}(idxpk,2)>3 && peakstemp{i}(idxpk,2)<length(x)-3
                [leftWidth, rightWidth, slmintegraltotal] = measurePeak_J(x,y,peakstemp{i},idxpk,ampThresh);
                if ~isempty(leftWidth)
                    peaks2{i}(pcount,:) = [pcount peakstemp{i}(idxpk,2) peakstemp{i}(idxpk,3) ...
                        peakstemp{i}(idxpk,2)-leftWidth peakstemp{i}(idxpk,2)+rightWidth slmintegraltotal];
                    %peak_count peak_center peak_height peak_leftx
                    %peak_rightx peak_integral
                    pcount = pcount + 1;
                end
            end
        end
        
    end
    
    for i = 2 %find and measure peaks of Msn4
        ampThresh = thre_for_idxpulse_Msn4;
        peakstemp{i} = findpeaks_v20190613(framenum_index,cellTrace{i},slopeThreshold,ampThresh); %find peaks
        x = framenum_index; y = cellTrace{i};
        pcount = 1;
        for idxpk = 1:size(peakstemp{i},1) %measure peaks and eliminate shoulder peaks.
            if ~isshoulderpeak(peakstemp{i},idxpk,x,y,ampThresh) && peakstemp{i}(idxpk,2)>3 && peakstemp{i}(idxpk,2)<length(x)-3
                [leftWidth, rightWidth, slmintegraltotal] = measurePeak_J(x,y,peakstemp{i},idxpk,ampThresh);
                if ~isempty(leftWidth)
                    peaks2{i}(pcount,:) = [pcount peakstemp{i}(idxpk,2) peakstemp{i}(idxpk,3) ...
                        peakstemp{i}(idxpk,2)-leftWidth peakstemp{i}(idxpk,2)+rightWidth slmintegraltotal];
                    %peak_count peak_center peak_height peak_leftx
                    %peak_rightx peak_integral
                    pcount = pcount + 1;
                end
            end
        end
    end
    
    %identify if each Msn2 pulse overlaps with any Msn4 pulses
    if size(peaks2{1},1)~=0
        peaks2{1}(:,7)=0;
        for index=1:size(peaks2{1},1)
            Msn2PulsePeak=peaks2{1,1}(index,2);
            if size(peaks2{2},1)==0
                peaks2{1}(index,7)=0;
            else
                for index2=1:size(peaks2{2},1)
                    if Msn2PulsePeak>=peaks2{1,2}(index2,4) && Msn2PulsePeak<=peaks2{1,2}(index2,5)
                        peaks2{1}(index,7)=1;
                        break;
                    end
                end
            end
        end
    end
    
    %identify if each Msn4 pulse overlaps with any Msn2 pulses
    if size(peaks2{2},1)~=0
        peaks2{2}(:,7)=0;
        for index=1:size(peaks2{2},1)
            Msn4PulsePeak=peaks2{1,2}(index,2);
            if size(peaks2{1},1)==0
                peaks2{2}(index,7)=0;
            else
                for index2=1:size(peaks2{1},1)
                    if Msn4PulsePeak>=peaks2{1,1}(index2,4) && Msn4PulsePeak<=peaks2{1,1}(index2,5)
                        peaks2{2}(index,7)=1;
                        break;
                    end
                end
            end
        end
    end
    
    finalResults{cellindex}.Msn2Pulses = peaks2{1};
    finalResults{cellindex}.Msn4Pulses = peaks2{2};
    finalResults{cellindex}.Msn2PulseFreq=60*size(peaks2{1},1)/timelimit; %Msn2 pulse frequency
    finalResults{cellindex}.Msn4PulseFreq=60*size(peaks2{2},1)/timelimit; %Msn4 pulse frequency
    finalResults{cellindex}.Msn2PulseNum=size(peaks2{1},1); %Msn2 pulse number
    finalResults{cellindex}.Msn4PulseNum=size(peaks2{2},1); %Msn4 pulse number
    if size(peaks2{1},1)==0
        finalResults{cellindex}.OverlapPulseNum=0; %coincident pulse number
    else
        finalResults{cellindex}.OverlapPulseNum=sum(peaks2{1}(:,7));
    end
    finalResults{cellindex}.Msn2OnlyPulseNum = finalResults{cellindex}.Msn2PulseNum-finalResults{cellindex}.OverlapPulseNum; %Msn2-only pulse number
    finalResults{cellindex}.Msn4OnlyPulseNum = finalResults{cellindex}.Msn4PulseNum-finalResults{cellindex}.OverlapPulseNum; %Msn4-only pulse number
    
    finalResults{cellindex}.Msn2OnlyPulseFreq=60*finalResults{cellindex}.Msn2OnlyPulseNum/timelimit; %Msn2-only pulse frequency
    finalResults{cellindex}.Msn4OnlyPulseFreq=60*finalResults{cellindex}.Msn4OnlyPulseNum/timelimit; %Msn4-only pulse frequency
    finalResults{cellindex}.OverlapPulseFreq=60*finalResults{cellindex}.OverlapPulseNum/timelimit; %coincident pulse frequency
      
end
%save([resultFolder '\analyzedresults.mat'],'finalResults');

Msn2PulseFreqVector=[];Msn4PulseFreqVector=[];
Msn2OnlyPulseFreqVector=[];Msn4OnlyPulseFreqVector=[];
OverlapPulseFreqVector=[];

for cellnum_i=1:cellnum
    Msn2PulseFreqVector(cellnum_i)=finalResults{cellnum_i}.Msn2PulseFreq;
    Msn4PulseFreqVector(cellnum_i)=finalResults{cellnum_i}.Msn4PulseFreq;
    Msn2OnlyPulseFreqVector(cellnum_i)=finalResults{cellnum_i}.Msn2OnlyPulseFreq;
    Msn4OnlyPulseFreqVector(cellnum_i)=finalResults{cellnum_i}.Msn4OnlyPulseFreq;
    OverlapPulseFreqVector(cellnum_i)=finalResults{cellnum_i}.OverlapPulseFreq; 
end
Msn2TotalFreq=mean(Msn2PulseFreqVector); %average of Msn2 pulse frequency among cells
Msn2TotalFreqSe=sqrt(var(Msn2PulseFreqVector))/sqrt(length(Msn2PulseFreqVector)); %standard error of Msn2 pulse frequency among cells

Msn4TotalFreq=mean(Msn4PulseFreqVector); %average of Msn4 pulse frequency among cells
Msn4TotalFreqSe=sqrt(var(Msn4PulseFreqVector))/sqrt(length(Msn4PulseFreqVector)); %standard error of Msn4 pulse frequency among cells

Msn2OnlyFreq=mean(Msn2OnlyPulseFreqVector); %average of Msn2-only pulse frequency among cells
Msn2OnlyFreqSe=sqrt(var(Msn2OnlyPulseFreqVector))/sqrt(length(Msn2OnlyPulseFreqVector)); %standard error of Msn2-only pulse frequency among cells

Msn4OnlyFreq=mean(Msn4OnlyPulseFreqVector); %average of Msn4-only pulse frequency among cells
Msn4OnlyFreqSe=sqrt(var(Msn4OnlyPulseFreqVector))/sqrt(length(Msn4OnlyPulseFreqVector)); %standard error of Msn4-only pulse frequency among cells

OverlapFreq=mean(OverlapPulseFreqVector); %average of coincident pulse frequency among cells
OverlapFreqSe=sqrt(var(OverlapPulseFreqVector))/sqrt(length(OverlapPulseFreqVector)); %standard error of coincident pulse frequency among cells

if cellnumforplot>0 %plot cells if cellnumforplot>0
    for cellnum_i_plot=1:cellnumforplot
        
        time_cellnum_i=finalResults{cellnum_i_plot}.time;
        PKAn_cellnum_i=finalResults{cellnum_i_plot}.PKA;
        Msn2n_cellnum_i=finalResults{cellnum_i_plot}.nuclearMsn2;
        Msn4n_cellnum_i=finalResults{cellnum_i_plot}.nuclearMsn4;
        Xn_cellnum_i=finalResults{cellnum_i_plot}.X;

        figure
        subplot(3,1,1);
        plot(time_cellnum_i,Msn2n_cellnum_i,'-black',time_cellnum_i,Msn4n_cellnum_i,'-r','LineWidth',2);
        xlim([0,timelimit])
        ylim([0,max(M2tot,M4tot)])
        legend('Msn2','Msn4');
        
        subplot(3,1,2);
        plot(time_cellnum_i,PKAn_cellnum_i,'-black','LineWidth',2);
        xlim([0,timelimit])
        ylim([minPKA-1,maxPKA+1])
        legend('PKA');
        
        subplot(3,1,3);
        plot(time_cellnum_i,Xn_cellnum_i,'-black','LineWidth',2);
        xlim([0,timelimit])
        ylim([minX-1,maxX+1])
        legend('X');
        
        saveas(gca,[result_folder,'\','cell',num2str(cellnum_i_plot),'.pdf']);
        close;
        
    end
end

%%%%%merge Msn2 (and Msn4) pulses of all cells
Msn2PulseMat=[];
Msn4PulseMat=[];
for cellnum_i=1:cellnum
    Msn2PulseMat=[Msn2PulseMat;finalResults{cellnum_i}.Msn2Pulses];
    Msn4PulseMat=[Msn4PulseMat;finalResults{cellnum_i}.Msn4Pulses];
end

%%%%%calculate pulse coincidence rate based on bootstrap
if size(Msn2PulseMat,1)==0 || size(Msn4PulseMat,1)==0
    OverlapFraction=0;
    OverlapFractionErrorbar=0;
else
    bootstrapNum=1000;
    Msn2PulseChosenNum=ceil(0.5*(size(Msn2PulseMat,1)+size(Msn4PulseMat,1)));
    Msn4PulseChosenNum=ceil(0.5*(size(Msn2PulseMat,1)+size(Msn4PulseMat,1)));
    OverlapFractionVector=[];
    for boot_i=1:bootstrapNum
        Msn2PulseOverlap=Msn2PulseMat(:,7);
        Msn4PulseOverlap=Msn4PulseMat(:,7);
        Msn2PulseOverlapChosen=Msn2PulseOverlap(randi(length(Msn2PulseOverlap),1,Msn2PulseChosenNum));
        Msn4PulseOverlapChosen=Msn4PulseOverlap(randi(length(Msn4PulseOverlap),1,Msn4PulseChosenNum));
        
        OverlapFractionVector(boot_i)=(sum(Msn2PulseOverlapChosen)+sum(Msn4PulseOverlapChosen))/(length(Msn2PulseOverlapChosen)+length(Msn4PulseOverlapChosen));
    end
    OverlapFraction=mean(OverlapFractionVector); %average pulse coincidence rate
    OverlapFractionErrorbar=errorbar(OverlapFractionVector); %95% confidence interval of pulse coincidence rate
    
end

function err=errorbar(data)
meanDATA=mean(data);
delta =data-meanDATA;
d =prctile(delta,[2.5,100-2.5]);
err=0.5*(d(2)-d(1));

%%%%%the following functions are used for identifying and measuring pulses
%%%%%from previous works

function P = findpeaks_v20190613(x,y,SlopeThreshold,AmpThreshold,pdist)
warning off;
smoothwidth = 1;
peakgroup = 3; % the number points  around the top part of the peak that are taken for measurement.
vectorlength=length(y); %frame numbers
d = smooth(deriv(y),smoothwidth); %smooth the first deriv. with a a moving average filter of size 3.
peak=1;
n=round(peakgroup/2+1);
AmpTest=AmpThreshold;
P = [];
for j=smoothwidth:length(y)-smoothwidth,
    if sign(d(j)) > sign (d(j+1)), % Detects zero-crossing
        if d(j)-d(j+1) > SlopeThreshold*y(j), % if slope of derivative is larger than SlopeThreshold
            if y(j) > AmpTest   % if height of peak is larger than AmpThreshold and pair distance below 6.
                PeakX = x(j);
                PeakY = y(j);
                P(peak,:) = [round(peak) PeakX PeakY];% MeasuredWidth  1.0646.*PeakY*MeasuredWidth]; %#ok<AGROW>
                peak=peak+1;
            end
        end
    end
end
%--------------------------------------------------------------------------
function d=deriv(a)
% First derivative of vector
n=length(a);
d=zeros(size(a));
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
    d(j)=a(j)-a(j-1);
end


%--------------------------------------------------------------------------
% if ~isshoulderpeak(peakstemp{i},idxpk,x,y,ampThresh)
function speak = isshoulderpeak(peaks,idxpk,x,y,ampThresh)
%determine if a given peak is a shoulder peak.

speak = 0;
idx_pkpos = find(x==peaks(idxpk,2));
%if it is a left-shoulder peak
for i = idxpk:size(peaks,1)
    if peaks(i,2) <= peaks(idxpk,2)+5 && peaks(i,3) > peaks(idxpk,3)
        heightThresh = 0.5*peaks(idxpk,3); %min(0.5*peaks(i,3),ampThresh);
        idx_peaki = find(x==peaks(i,2));
        if min(y(idx_pkpos:idx_peaki))>heightThresh %if all points between these two peaks are above the half-max of the test peak.
            speak = 1;
        end
    end
end

%if it is a right-shoulder peak
for  i = 1:idxpk
    if peaks(i,2) >= peaks(idxpk,2)-5 && peaks(i,3) > peaks(idxpk,3)
        heightThresh = 0.5*peaks(idxpk,3); %min(0.5*peaks(i,3),ampThresh);
        idx_peaki = find(x==peaks(i,2));
        if min(y(idx_peaki:idx_pkpos))>heightThresh %if all points between these two peaks are above the half-max of the test peak.
            speak = 1;
        end
    end
end


%--------------------------------------------------------------------------
function [leftWidth, rightWidth, slmintegraltotal] = measurePeak_J(x,y,P,idxpk,ampThresh)
%measure peak width by fitting each half of the peak with spline.
warning off;
leftWidth = []; rightWidth = []; slmintegraltotal = [];
xleft = P(idxpk,2)-3;  %idxpk means the pulse number, 2 means cols(3 frame)
idx_xleft = find(x==xleft);
while isempty(idx_xleft)
    xleft = xleft -1;
    idx_xleft = find(x==xleft);
end
idx_pkpos = find(x==P(idxpk,2));  %peak time
lefthalfx = x(idx_xleft:idx_pkpos);  %from 3 frame before peak to peak
lefthalfy = y(idx_xleft:idx_pkpos);  %amp from 3 frame before peak to peak
%           if  P(idxpk,3)>150 && P(idxpk,3)<200 % P(idxpk,2)>150 && P(idxpk,2)<180  && P(idxpk,3)>10 && P(idxpk,3)<20
%                pause(1)
%            end
heightThresh = min(0.5*lefthalfy(end),ampThresh);
%heigth treshold for resetting a peak is either 50% of the peak
%height or the threshold for defining a peak.
if max(lefthalfy)>lefthalfy(end) && min(lefthalfy(find(lefthalfy==max(lefthalfy)):end))>heightThresh  %shoulder peak
    return;
end

if min(lefthalfy)>heightThresh %4 points are not enough to reach threshold height.
    flag = 0;
    while idx_xleft>1 && flag==0 %idx_xleft is the first frame before peak
        idx_xleft = idx_xleft -1; %search to the last frame
        if y(idx_xleft)< heightThresh
            flag = 1;
        elseif max(y(idx_xleft:idx_pkpos))>lefthalfy(end) %current peak is a shoulder peak.
            return;
            %elseif idx_xleft == 1 && y(idx_xleft)> heightThresh
            %    return;
        end
    end
else
    %find the first left point closest to the peak that is below heightThresh
    idxl = find(x==max(lefthalfx(lefthalfy<=heightThresh)));
    idx_xleft = idxl;
end
lefthalfx = x(idx_xleft:idx_pkpos);
lefthalfy = y(idx_xleft:idx_pkpos);

xright= P(idxpk,2)+3;
idx_xright = find(x==xright);
while isempty(idx_xright)
    xright = xright+1;
    idx_xright = find(x==xright);
end
righthalfx = x(idx_pkpos:idx_xright); %determine right half of the peak
righthalfy = y(idx_pkpos:idx_xright);
if max(righthalfy)>righthalfy(1) && min(righthalfy(1:find(righthalfy==max(righthalfy))))>heightThresh
    return;
end

if min(righthalfy)>heightThresh %4 points are not enough to reach half height.
    flag = 0;
    while idx_xright<length(x) && flag==0
        idx_xright = idx_xright + 1;
        if y(idx_xright)< heightThresh
            flag = 1;
        elseif max(y(idx_pkpos:idx_xright))>righthalfy(1) %current peak is a shoulder peak.
            return;
            % elseif idx_xright == length(x) && y(idx_xright)> heightThresh
            %     return;
        end
    end
    %  elseif max(righthalfy(find(righthalfy==min(righthalfy)):end))>0.5*righthalfy(1)
    %xright = P(idxpk,2)+find(righthalfy==min(righthalfy))-1;
    %     idx_xright = idx_pkpos + find(righthalfy==min(righthalfy))-1;
else
    %find the first point right to the peak that is below heightThresh
    idxr = find(x==min(righthalfx(righthalfy<=heightThresh)));
    idx_xright = idxr;
end
righthalfx = x(idx_pkpos:idx_xright);
righthalfy = y(idx_pkpos:idx_xright);

%%%%%%%%up is all about find the left/right closest frame of peak

%         if  P(idxpk,3)>70 && P(idxpk,3)<400 % P(idxpk,2)>150 && P(idxpk,2)<180  && P(idxpk,3)>10 && P(idxpk,3)<20
%                 pause(1)
%         end

slmleft = slmengine(lefthalfx,lefthalfy,'knots',4,'plot','off','increasing','on');
%change off to on
%leftWidth = lefthalfx(end)-slmeval(0.5*slmeval(lefthalfx(end),slmleft),slmleft,-1);
slmlefteval = slmeval(heightThresh,slmleft,-1);
if isnan(slmlefteval)
    slmlefteval = lefthalfx(1);
    slmintegralleft = 0;
else
    slmintegralleft = slmpar(slmleft,'in');  %'in' refers to the integral
end
leftWidth = lefthalfx(end)-slmlefteval;

if idx_xleft == 1
    leftWidth = lefthalfx(end)-1;
end

slmright = slmengine(righthalfx,righthalfy,'knots',4,'plot','off','decreasing','on');
%rightWidth = slmeval(0.5*slmeval(righthalfx(1),slmright),slmright,-1)-righthalfx(1);
slmrighteval = slmeval(heightThresh,slmright,-1);
if isnan(slmrighteval)
    slmrighteval = righthalfx(end);
    slmintegralright = 0;
else
    slmintegralright = slmpar(slmright,'in');  %'in' refers to the integral
end
slmintegraltotal = slmintegralright +slmintegralleft; %the total integral
rightWidth = slmrighteval-righthalfx(1);
if idx_xright == length(x)
    rightWidth = length(x)-righthalfx(1);
end


function [finalMean, stdSample, inlierIdx, outlierIdx] = robustMean(data,dim,k,fit,logical_idx)
%ROBUSTMEAN calculates mean and standard deviation discarding outliers
%
% SYNOPSIS [finalMean, stdSample, inlierIdx, outlierIdx] = robustMean(data,dim,k,fit)
%
% INPUT    data : input data
%          dim  : (opt) dimension along which the mean is taken {1}
%          k    : (opt) #of sigmas at which to place cut-off {3}
%          fit  : (opt) whether or not to use fitting to robustly estimate
%                  the mean from the data that includes outliers.
%                  0 (default): mean is approximated by median(data)
%                  1 : mean is approximated by
%                      fminsearch(@(x)(median(abs(data-x))),median(data))
%                      This option is only available for scalar data
%          logical_idx : (opt) logical indicating whether to return binary
%                        indicies for inlierIdx and outlierIdx
%                        false (default): actual idx, but slow due to find
%                        true           : logical idx, faster
%
%
% OUTPUT   finalMean : robust mean
%          stdSample : std of the data (divide by sqrt(n) to get std of the
%                      mean)
%          inlierIdx : index into data with the inliers
%          outlierIdx: index into data with the outliers
%
% REMARKS  NaN or Inf will be counted as neither in- nor outlier
%          The code is based on (linear)LeastMedianSquares. It could be changed to
%          include weights
%
% REFERENCES
%
% 1. PJ Rousseeuw and AM Leroy. Robust Regression and Outlier Detection.
%    New York: John Wiley & Sons, 1987. ISBN: 978-0-471-48855-2.
% 2. PJ Rousseeuw. Least median of squares regression. J. Am. Stat. Ass.
%    79, 871-880 (1984). DOI: 10.1080/01621459.1984.10477105.
% 3. G Danuser and M Stricker. Parameteric Model Fitting: From Inlier
%    Characterization to Outlier Detection. IEEE Trans Pattern Anal. Mach.
%    Intell. Vol 20, No. 2, March 1998. DOI: 10.1109/34.667884.
% 4. K Jaqaman and G Danuser. Linking Data to models: data regression.
%    Nature Rev Mol Cell Bio. 7, 813-819 (Nov 2006). DOI: 10.1038/nrm2030.
%
% c: jonas, 04/04
% Mark Kittisopikul, November 2015
%
% See also mad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (C) 2017, Danuser Lab - UTSouthwestern
%
% This file is part of u-track.
%
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
%
%


% test input

if isempty(data)
    error('Please supply non-empty data to robustMean')
end
if nargin<2 || isempty(dim)
    % make sure that the dimensinon is correct if there's a vector
    if any(size(data)==1) && ismatrix(data)
        dim = find(size(data)>1);
    else
        dim = 1;
    end
end
if nargin < 3 || isempty(k)
    % Cut-off is roughly at 3 sigma (References 1-3) by default
    k = 3;
end
% mkitti: was bug? only four parameters possible
% if nargin < 5 || isempty(fit)
if nargin < 4 || isempty(fit)
    fit = 0;
end
if nargin < 5 || isempty(logical_idx)
    logical_idx = false;
end

if fit == 1
    % check for vector
    if sum(size(data)>1)>1
        error('fitting is currently only supported for 1D data')
    end
end

insufficientData = true;
if(numel(data) >= 4)
    if(all(isfinite(data(1:4))) ...
            || all(isfinite(data(end-3:end))) ...
            || all(isfinite(data(((1:4)+floor(end/2)-2))))   )
        % Quickly check if first, last, or middle four are finite
        insufficientData = false;
    else
        finiteMap = isfinite(data);
        % Only need to find four
        finiteCount = numel(find(finiteMap,4));
        insufficientData = finiteCount < 4;
    end
end

% if sum(isfinite(data(:))) < 4
if insufficientData
    warning('ROBUSTMEAN:INSUFFICIENTDATA',...
        'Less than 4 data points!')
    finalMean = nanmean(data,dim);
    stdSample = NaN(size(finalMean));
    inlierIdx = find(isfinite(data));
    outlierIdx = [];
    return
end


%========================
% LEAST MEDIAN SQUARES
%========================

% Scale factor that relates Median absolute deviation (MAD) to standard deviation
% See mad(X,1)
% mad2stdSq=1.4826^2; %see same publications

% mad2stdSq = 1/norminv(3/4).^2
% mad2stdSq = 2.198109338317732

% backwards compatible constant
% sprintf('%0.9g',1.4826^2)
mad2stdSq = 2.19810276;

% calc median - reduce dimension dim to length 1
if fit
    % minimize the median deviation from the mean
    medianData = fminsearch(@(x)(median(abs(data-x))),median(data));
else
    medianData = nanmedian(data,dim);
end

% calculate statistics
% res2 = (data-repmat(medianData,blowUpDataSize)).^2;
res2 = bsxfun(@minus,data,medianData).^2;

medRes2 = max(nanmedian(res2,dim),eps);

%testvalue to calculate weights
% testValue=res2./repmat(mad2stdSq*medRes2,blowUpDataSize);
testValue = bsxfun(@rdivide,res2,mad2stdSq*medRes2);

% outlierIdx = testValue > k^2;
% Note: NaNs will always be false in comparison
inlierIdx = testValue <= k^2;
% Old outliers:
% outlierIdx = find(testValue>k^2); % Does not include NaNs
% New outliers:
outlierIdx = ~inlierIdx; % Also includes NaNs

% Prior to Nov 2015, there used to be an if/else statement here depending if
% vector or higher dimensional input was given for data.
nInliers = sum(inlierIdx,dim);

% calculate std of the sample;
if nargout > 1
    
    %% Obsolete code block with bug
    % put NaN wherever there are not enough data points to calculate a
    % standard deviation
    %         goodIdx = sum(isfinite(res2),dim) > 4;
    %mkitti, Oct 29 2015
    % I believe the following commented out lines constitute a bug.
    % goodIdx does not correctly index res2 in the expected manner.
    % Therefore the second output of robustMean.m when supplied with a
    % multidimensional input is invalid.
    %         stdSample = NaN(size(goodIdx));
    %         stdSample(goodIdx)=sqrt(nansum(res2(goodIdx),dim)./(nInliers(goodIdx)-4));
    
    %% outlierIdx should send NaN to zeros also so nansum not needed
    res2(outlierIdx) = 0;
    stdSample = sqrt(sum(res2,dim)./(nInliers-4));
    stdSample(nInliers <= 4) = NaN;
end

%====END LMS=========

%======
% MEAN
%======

data(outlierIdx) = 0;
finalMean = sum(data,dim)./nInliers;

if(nargout > 2 && ~logical_idx)
    % For backwards compatability only
    inlierIdx = find(inlierIdx);
end
if(nargout > 3)
    % For backwards compatability only, exclude NaNs
    % Above, NaNs are included as outliers
    % Next two line should be equivalent
    %      outlierIdx = testValue > k^2;
    outlierIdx(outlierIdx) = ~isnan(testValue(outlierIdx));
    if(~logical_idx)
        outlierIdx = find(outlierIdx);
    end
end

function [slm,xp,yp] = slmengine(x,y,varargin)
warning off;
% slmengine: estimates a spline function from data plus a fit prescription
% usage 1: slm = slmengine(x,y);
% usage 2: slm = slmengine(x,y,prescription);
% usage 3: slm = slmengine(x,y,prop1,val1,prop2,val2,...);
% usage 4: slm = slmengine(x,y,prescription,prop1,val1,prop2,val2,...);
% usage 5: [slm,xp,yp] = slmengine(x,y,prescription,prop1,val1,prop2,val2,...);
%
% Note: slmengine is the command driven tool for fitting a model
%   to data. It uses either a prescription structure (as supplied by
%   slmset) or sets of property/value pairs. Those pairs are defined
%   in the help to slmset. slmengine is also used by slmfit, the gui
%   tool for fitting a curve to data.
%
% Note: The optimization toolbox (lsqlin) is generally required for
% most fits. Some models (those where the break points are also
% estimated) require fmincon.
%
% Note: Some prescriptive parameters or combinations will be
%   inappropriate for the lower degree models.
%
% arguments: (input)
%  x,y     - vectors of data used to fit the model. Any
%            NaN or inf elements (in x or y) will be excluded
%            from the fit.
%
%  prescription - structure generated by slmset to control the
%            fitting process.
%
%            Alternatively, one can simply supply property/value
%            pairs directly to slmengine.
%
% Arguments: (output)
%  slm     - model as a shape-language-model structure,
%            normally constructed by slmfit or slmengine. slm
%            may also be returned in other forms, either a pp
%            structure tht ppval can use, or as a simple array
%            of coefficients.
%
%            Upon return, slm will contain a strucure field 'stats'
%            that holds fields which describe the quality of fit:
%
%            slm.stats.TotalDoF = Total degrees of freedom in the model
%
%            slm.stats.NetDoF = Net degrees of freedom in the spline model.
%               Thus, NetDoF reflects any equality constraints in the model
%
%            slm.stats.R2 = Traditional R-squared coefficient
%
%            slm.stats.R2Adj = Adjusted R^2, accounting for changing
%               degrees of freedom. For a definition, look here:
%
%               http://en.wikipedia.org/wiki/Coefficient_of_determination
%
%            slm.stats.RMSE = Root-Mean-Squared-Error
%
%            slm.stats.ErrorRange = 1x2 vector that contains the minimmum
%               and maximum errors, defined as (Yhat - Y)
%
%            slm.stats.Quartiles = 1x2 vector that contains the 25% and 75%
%               error quartiles. This gives a trimmed measure of the
%               error magnitudes, less any outliers.
%
%  xp, yp  - predicted points along the fitted curve. The number of these
%            points is defined by the predictions property. (see slmset)
%            xp will be a list of equally spaced points.
%
% Example:
%  (See the SLM_tutorial demo for other examples.)
%  x = rand(1,50);
%  y = sin(x*2*pi) + randn(size(x))/10;
%
%  slm = slmengine(x,y,'knots',0:.2:1,'plot','on','concavedown','on','minvalue',0)
%  slm =
%            form: 'slm'
%          degree: 3
%           knots: [6x1 double]
%            coef: [6x2 double]
%    prescription: [1x1 struct]
%               x: [50x1 double]
%               y: [50x1 double]
%
% Returned in a form that ppval or fnval can use:
%
%  slm = slmengine(x,y,'knots',0:.2:1,'plot','on','concavedown','on','minvalue',0,'res','pp')
% slm =
%            form: 'pp'
%          breaks: [0 0.2 0.4 0.6 0.8 1]
%           coefs: [5x4 double]
%          pieces: 5
%           order: 4
%             dim: 1
%    prescription: [1x1 struct]
%
%
% See also: spap2, ppval, fnval
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 1/2/08


% which form does the shape prescription take? Or do we
% just use our defaults for the prescription?
if nargin<2
    error('SLMENGINE:improperdata','Must supply two vectors to fit a model: x & y')
elseif nargin==2
    prescription = slmset;
elseif (nargin>=3)
    % prescription supplied directly or as property/value pairs,
    % or as a prescription, modified by a few set of pairs.
    % slmset will resolve all of these cases.
    prescription = slmset(varargin{:});
end

% check the data for size, turning it into column vectors
x = x(:);
y = y(:);
n = length(x);
if n~=length(y)
    error('SLMENGINE:inconsistentdata','x and y must be the same size')
end

% were there any NaN or inf elements in the data?
k = isnan(x) | isnan(y) | isinf(x) | isinf(y);
if any(k)
    % drop them from the analysis
    x(k) = [];
    y(k) = [];
    
    % also drop corresponding weights if they were supplied
    if ~isempty(prescription.Weights)
        prescription.Weights(k) = [];
    end
    n = length(x);
end

% if weights or errorbars were set, verify the sizes of these
% parameters, compared to the number of data points.
if ~isempty(prescription.Weights)
    prescription.Weights = prescription.Weights(:);
    if n~=length(prescription.Weights)
        error('SLMENGINE:inconsistentweights','Weights vector must be the same length as # of data points')
    end
end
if ~isempty(prescription.ErrorBar)
    EB = prescription.ErrorBar;
    if length(EB) == 1
        prescription.ErrorBar = repmat(EB,n,2);
    elseif (size(EB,1) == 1)
        prescription.ErrorBar = repmat(EB',1,2);
    elseif (size(EB,2) == 1)
        prescription.ErrorBar = repmat(EB,1,2);
    end
end

if prescription.Verbosity > 1
    disp('=========================================')
    disp('Model Prescription')
    disp(prescription)
end

% is this a free knot problem?
if strcmp(prescription.InteriorKnots,'free')
    % the alternative is normally 'fixed', but its free now
    
    % get the starting values for the knots
    % knots vector, dx
    if length(prescription.Knots)==1
        [knots,nk] = chooseknots(prescription.Knots,n,x);
    else
        % we should check that the knots contain the data
        knots = sort(prescription.Knots(:));
        nk= length(knots);
        if (knots(1)>min(x)) || (knots(end)<max(x))
            error('SLMENGINE:inadequateknots',['Knots do not contain the data. Data range: ',num2str([min(x),max(x)])])
        end
    end
    % there must be at least 3 knots
    if nk<3
        error('SLMENGINE:inadequateknots','Free knot estimation requires at least 3 knots')
    end
    
    % set the constraints for fmincon
    % 1. no two knots may lie too close to each other
    tol = min(0.001,0.1/(nk-1));
    mindelta = tol*(knots(end) - knots(1));
    A = zeros(2,nk-2);
    A(1,1) = -1;
    A(2,nk-2) = 1;
    b = [-knots(1) - mindelta ; knots(end) - mindelta];
    A = [A;full(spdiags(repmat([1 -1],nk-3,1),[0 1],nk-3,nk-2))];
    b = [b;repmat(-mindelta,nk-3,1)];
    
    % set up the optimization parameters (fmincon)
    fminconoptions = optimset('fmincon');
    fminconoptions.LargeScale = 'off';
    fminconoptions.Algorithm = 'active-set';
    if prescription.Verbosity>1
        fminconoptions.Display = 'iter';
    else
        fminconoptions.Display = 'off';
    end
    
    % change the prescription for the subsequent internal calls
    prescrip = prescription;
    prescrip.Knots = knots;
    prescrip.InteriorKnots = 'fixed';
    prescrip.Plot = 'off';
    
    % call the optimizer
    intknots = knots(2:(end-1));
    intknots = fmincon(@free_knot_obj,intknots,A,b, ...
        [],[],[],[],[],fminconoptions,x,y,prescrip);
    
    % fit the curve one last time with the final knots
    prescrip.Knots(2:(end-1)) = intknots;
    slm = slmengine(x,y,prescrip);
    
else
    % the knots are fixed, just estimate the model
    
    % do the appropriate fit. break them down into special
    % cases, purely for simplicity of the code.
    switch prescription.Degree
        case {0 'constant'}
            slm = slmengine_constant(x,y,prescription);
        case {1 'linear'}
            slm = slmengine_linear(x,y,prescription);
        case {3 'cubic'}
            slm = slmengine_cubic(x,y,prescription);
    end
    
end

% do we need to report anything? I.e., what was the
% prescription.Verbosity setting?
if prescription.Verbosity>0
    disp('=========================================')
    disp('MODEL STATISTICS REPORT')
    disp(['Number of data points:      ',num2str(n)])
    disp(['Total degrees of freedom:   ',num2str(slm.stats.TotalDoF)])
    disp(['Net degrees of freedom:     ',num2str(slm.stats.NetDoF)])
    disp(['R-squared:                  ',num2str(slm.stats.R2)])
    disp(['Adjusted R-squared:         ',num2str(slm.stats.R2Adj)])
    disp(['RMSE:                       ',num2str(slm.stats.RMSE)])
    disp(['Range of prediction errors: ',num2str(slm.stats.ErrorRange)])
    disp(['Error quartiles (25%, 75%): ',num2str(slm.stats.Quartiles)])
    disp('=========================================')
end

% append the shape prescription to the final structure
% this can be used as a template for fitting other functions,
% but also as documentation for the curve fit. This is why
% it is attached to the model structure, rather than as a
% second return argument.
slm.prescription = prescription;
% also attach the data used to build the model, purely
% for documentation purposes.
slm.x = x;
slm.y = y;

% do we need to plot the curve? I.e., what was the
% prescription.Plot setting?
if strcmp(prescription.Plot,'on')
    % plot the curve
    plotslm(slm)
end

% do we generate points on the curve?
if isempty(prescription.Predictions)
    % none were asked for
    xp = [];
    yp = [];
else
    % generate the required points on the curve
    xp = linspace(min(x),max(x),prescription.Predictions);
    yp = slmeval(xp,slm,0);
end

% do they want the result in a 'pp' form?
if strcmpi(prescription.Result,'pp') && ~strcmpi(slm.form,'pp')
    % convert to pp form
    stats = slm.stats;
    slm = slm2pp(slm);
    
    % make sure the stats field gets copied into the pp form
    slm.stats = stats;
end



% ========================================================
% =========== slmengines for each model degree ============
% ========================================================

% ========================================================
% ============== piecewise constant model ================
% ========================================================
function slm = slmengine_constant(x,y,prescription)
% fits a piecewise constant shape prescriptive model

% check for inappropriate properties for a piecewise
% constant model
property_check(prescription,'constant')

% simple things about data first...
% slmengine has already made it a column vector,
% and ensured compatibility in length between
% x and y
nx = length(x);
% knots vector, dx
if length(prescription.Knots)==1
    [knots,nk] = chooseknots(prescription.Knots,nx,x);
else
    % we should check that the knots contain the data
    knots = sort(prescription.Knots(:));
    nk= length(knots);
    if (knots(1)>min(x)) || (knots(end)<max(x))
        error('SLMENGINE:inadequateknots',['Knots do not contain the data. Data range: ',num2str([min(x),max(x)])])
    end
end
dx = diff(knots);
if any(dx==0)
    error('SLMENGINE:indistinctknots','Knots must be distinct.')
end

% number of coefficients to estimate
% a piecewise constant function has one coefficient
% at each knot, but the last knot is "not".
nc = nk-1;

% create empty design, equality, inequality
% constraints arrays. also rhs vectors for each
Mdes = zeros(0,nc);
Meq = Mdes;
Mineq = Mdes;
Mreg = Mdes;

rhseq = [];
rhsineq = [];
rhsreg = [];

% -------------------------------------
% build design matrix -
% first, bin the data - histc wll do it
[junk,xbin] = histc(x,knots); %#ok
% any point which falls at the top end, is said to
% be in the last bin.
xbin(xbin==nk)=nk-1;
% the design matrix is easy to build - one call to sparse
Mdes = accumarray([(1:nx)',xbin],1,[nx,nc]);
rhs = y;

% have we been given a constraint on the sum of the residuals?
if ~isempty(prescription.SumResiduals)
    Meq = [Meq;sum(Mdes,1)];
    rhseq = [rhseq;prescription.SumResiduals + sum(rhs)];
end

% apply weights
W = prescription.Weights;
if ~isempty(W)
    W = W(:);
    if length(W)~=nx
        error('SLMENGINE:inconsistentweights','Weight vector is not the same length as data')
    end
    
    % scale weight vector
    W = sqrt(nx)*W./norm(W);
    
    W = spdiags(W,0,nx,nx);
    Mdes = W*Mdes;
    rhs = W*rhs;
end

% -------------------------------------
% end conditions do not apply for a piecewise
% constant function

% -------------------------------------
% build regularizer
if nc>2
    % use a simple second order finite difference
    dx1 = dx(1:(nc-2));
    dx2 = dx(2:(nc-1));
    fda = [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), ...
        -2./(dx2.*(dx1+dx2))];
    
    Mreg = zeros(nc-2,nc);
    for i=1:(nc-2);
        Mreg(i,i+[0 1 2]) = fda(i,:);
    end
    rhsreg = zeros(nc-2,1);
end
% scale the regularizer before we apply the
% regularization parameter.
Mreg = Mreg/norm(Mreg,1);

% -------------------------------------
% single point equality constraints
% left hand side
if ~isempty(prescription.LeftValue)
    M = zeros(1,nc);
    M(1) = 1;
    Meq = [Meq;M];
    rhseq = [rhseq;prescription.LeftValue];
end

% right hand side
if ~isempty(prescription.RightValue)
    M = zeros(1,nc);
    M(end) = 1;
    Meq = [Meq;M];
    rhseq = [rhseq;prescription.RightValue];
end

% force curve through an x-y pair
xy = prescription.XY;
if ~isempty(xy)
    n = size(xy,1);
    if any(xy(:,1)<knots(1)) || any(xy(:,1)>knots(end))
        error('SLMENGINE:improperconstraint','XY pairs to force the curve through must lie inside the knots')
    end
    
    [junk,ind] = histc(xy(:,1),knots); %#ok
    ind(ind==(nc+1))=nc;
    M = sparse((1:n)',ind,1,n,nc);
    Meq = [Meq;M];
    rhseq = [rhseq;xy(:,2)];
end

% -------------------------------------
% Integral equality constraint
if ~isempty(prescription.Integral)
    % Rectangle rule. The last knot point
    % adds no contribution to the area, but
    % there are only nk-1 knots anyway.
    M = dx(:)';
    Meq = [Meq;M];
    rhseq = [rhseq;prescription.Integral];
end

% -------------------------------------
% single point inequality constraints
% left hand side minimum
if ~isempty(prescription.LeftMinValue)
    M = zeros(1,nc);
    M(1) = -1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;-prescription.LeftMinValue];
end
% left hand side maximum
if ~isempty(prescription.LeftMaxValue)
    M = zeros(1,nc);
    M(1) = 1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;prescription.LeftMinValue];
end

% right hand side min
if ~isempty(prescription.RightMinValue)
    M = zeros(1,nc);
    M(end) = -1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;-prescription.RightMinValue];
end

% right hand side max
if ~isempty(prescription.RightMaxValue)
    M = zeros(1,nc);
    M(end) = 1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;prescription.RightMaxValue];
end

% -------------------------------------
% Error bar inequality constraints
if ~isempty(prescription.ErrorBar)
    % lower bounds
    Mineq = [Mineq;-Mdes];
    rhsineq = [rhsineq;-(y-prescription.ErrorBar(:,1))];
    
    % upper bounds
    Mineq = [Mineq;Mdes];
    rhsineq = [rhsineq;y+prescription.ErrorBar(:,2)];
end

% -------------------------------------
% constant region(s)
if ~isempty(prescription.ConstantRegion)
    % there was at least one constant regions specified
    cr = prescription.ConstantRegion;
    M = zeros(0,nc);
    n = 0;
    for i=1:size(cr,1)
        % enforce constancy between the given range limits
        for j=1:(nc-1)
            if (knots(j+1)>=cr(i,1)) && (knots(j+1)<cr(i,2))
                n=n+1;
                M(n,j+[0 1]) = [-1 1];
            end
        end
    end
    
    Meq = [Meq;M];
    rhseq = [rhseq;zeros(size(M,1),1)];
    
end

% -------------------------------------
% overall inequalities
% global min value
if ~isempty(prescription.MinValue)
    M = -eye(nc,nc);
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;repmat(-prescription.MinValue,nc,1)];
end

% global max value
if ~isempty(prescription.MaxValue)
    M = eye(nc,nc);
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;repmat(prescription.MaxValue,nc,1)];
end

% SimplePeak and SimpleValley are really just composite
% properties. A peak at x == a is equivalent to a monotone
% increasing function for x<=a, and a monotone decreasing
% function for x>=a. Likewise a valley is just the opposite.
%
% So specifying a peak or a valley will cause any monotonicity
% specs to be ignored. slmset has already checked for this,
% and turned off any other monotonicity based properties.
sp = prescription.SimplePeak;
if isnumeric(sp) && ~isempty(sp)
    prescription.Increasing = [knots(1),sp];
    prescription.Decreasing = [sp,knots(end)];
end
sv = prescription.SimpleValley;
if isnumeric(sv) && ~isempty(sv)
    prescription.Decreasing = [knots(1), sv];
    prescription.Increasing = [sv,knots(end)];
end

% monotonicity?
% increasing regions
incR = prescription.Increasing;
L=0;
if ischar(incR)
    if strcmp(incR,'on')
        L=L+1;
        mono(L).knotlist = 1:nc;
        mono(L).direction = 1;
        mono(L).range = [];
    end
elseif ~isempty(incR)
    for i=1:size(incR,1)
        L=L+1;
        mono(L).knotlist = []; %#ok
        mono(L).direction = 1; %#ok
        mono(L).range = sort(incR(i,:)); %#ok
    end
end
% decreasing regions
decR = prescription.Decreasing;
if ischar(decR)
    if strcmp(decR,'on')
        L=L+1;
        mono(L).knotlist = 1:nc;
        mono(L).direction = -1;
        mono(L).range = [];
    end
elseif ~isempty(decR)
    for i=1:size(decR,1)
        L=L+1;
        mono(L).knotlist = [];
        mono(L).direction = -1;
        mono(L).range = sort(decR(i,:));
    end
end
if L>0
    % there were at least some monotone regions specified
    M = zeros(0,nc);
    n = 0;
    for i=1:L
        if isempty(mono(L).range)
            % the entire range was specified to be monotone
            for j=1:(nc-1)
                n=n+1;
                M(n,j+[0 1]) = [1 -1]*mono(i).direction;
            end
        else
            % only enforce monotonicity between the given range limits
            for j=1:(nc-1)
                if (knots(j)<mono(i).range(2)) && ...
                        (knots(j+1)>=mono(i).range(1))
                    n=n+1;
                    M(n,j+[0 1]) = [1 -1]*mono(i).direction;
                end
            end
            
        end
    end
    
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;zeros(size(M,1),1)];
    
end

% -------------------------------------
% Use envelope inequalities?
switch prescription.Envelope
    case 'off'
        % nothing to do in this case
    case 'supremum'
        % yhat - y >= 0
        Mineq = [Mineq;-Mdes];
        rhsineq = [rhsineq;-rhs];
    case 'infimum'
        % yhat - y <= 0
        Mineq = [Mineq;Mdes];
        rhsineq = [rhsineq;rhs];
end

% -------------------------------------
% scale equalities for unit absolute row sum
if ~isempty(Meq)
    rs = diag(1./sum(abs(Meq),2));
    Meq = rs*Meq;
    rhseq = rs*rhseq;
end
% scale inequalities for unit absolute row sum
if ~isempty(Mineq)
    rs = diag(1./sum(abs(Mineq),2));
    Mineq = rs*Mineq;
    rhsineq = rs*rhsineq;
end

% -------------------------------------
% now worry about the regularization. There are three
% possible cases.
% 1. We have a given regularization parameter
% 2. We have a given rmse that we wish to match
% 3. We must use cross validation to choose the parameter
RP = prescription.Regularization;
if (isnumeric(RP) && (RP>=0)) || ((ischar(RP)) && (strcmpi(RP,'smoothest')))
    % solve the problem using the given regularization parameter
    coef = solve_slm_system(RP,Mdes,rhs,Mreg,rhsreg, ...
        Meq,rhseq,Mineq,rhsineq,prescription);
elseif isnumeric(RP) && (RP<0)
    % we must match abs(RP) as the rmse.
    aim_rmse = abs(RP);
    fminbndoptions = optimset('fminbnd');
    RP = fminbnd(@match_rmse,-6,6,fminbndoptions, ...
        Mdes,rhs,Mreg,rhsreg,Meq,rhseq,Mineq,rhsineq, ...
        aim_rmse,prescription);
    
    % we logged the parameter in the optimization. undo that for
    % the final call
    RP = 10^RP;
    
    % do one final call to get the final coefficients
    coef = solve_slm_system(RP,Mdes,rhs,Mreg,rhsreg, ...
        Meq,rhseq,Mineq,rhsineq,prescription);
    
elseif ischar(RP)
    % its a cross validation problem to solve. drop out each point
    % in turn from the model, then use fminbnd to minimize the
    % predicted sum of squares at the missing points.
    fminbndoptions = optimset('fminbnd');
    % fminbndoptions.Display = 'iter';
    RP = fminbnd(@min_cv,-6,6,fminbndoptions, ...
        Mdes,rhs,Mreg,rhsreg,Meq,rhseq,Mineq,rhsineq,prescription);
    
    % we logged the parameter in the optimization. undo that for
    % the final call
    RP = 10^RP;
    
    % do one final call to get the final coefficients
    coef = solve_slm_system(RP,Mdes,rhs,Mreg,rhsreg, ...
        Meq,rhseq,Mineq,rhsineq,prescription);
end

% -------------------------------------
% unpack coefficients into the result structure
slm.form = 'slm';
slm.degree = 0;
slm.knots = knots;
slm.coef = coef;

% degrees of freedom available
slmstats.TotalDoF = nk - 1;
slmstats.NetDoF = slmstats.TotalDoF - size(Meq,1);
% this function does all of the stats, stuffing into slmstats
slm.stats = modelstatistics(slmstats,Mdes,y,coef);

% ========================================================
% =============== piecewise linear model =================
% ========================================================
function slm = slmengine_linear(x,y,prescription)
% fits a piecewise linear shape prescriptive model

% check for inappropriate properties for a piecewise
% linear model
property_check(prescription,'linear')

% simple things about data first...
% slmengine has already made it a column vector,
% and ensured compatibility in length between
% x and y
nx = length(x);

% knots vector, dx
if length(prescription.Knots)==1
    % just given a number of knots
    [knots,nk] = chooseknots(prescription.Knots,nx,x);
else
    % we should check that the knots contain the data
    knots = sort(prescription.Knots(:));
    nk= length(knots);
    if (knots(1)>min(x)) || (knots(end)<max(x))
        error('SLMENGINE:inadequateknots',['Knots do not contain the data. Data range: ',num2str([min(x),max(x)])])
    end
end
dx = diff(knots);
if any(dx==0)
    error('SLMENGINE:indistinctknots','Knots must be distinct.')
end




% number of coefficients to estimate.
% a piecewise linear Hermite has one coefficient
% at each knot.
nc = nk;

% create empty design, equality, inequality
% constraints arrays. also rhs vectors for each
Mdes = zeros(0,nc);
Meq = Mdes;
Mineq = Mdes;
Mreg = Mdes;

rhseq = [];
rhsineq = [];
rhsreg = [];

% -------------------------------------
% build design matrix -
% first, bin the data - histc wll do it
[junk,xbin] = histc(x,knots); %#ok
% any point which falls at the top end, is said to
% be in the last bin.
xbin(xbin==nk)=nk-1;

% build design matrix
t = (x - knots(xbin))./dx(xbin);
Mdes = accumarray([repmat((1:nx)',2,1), ...
    [xbin;xbin+1]],[1-t;t],[nx,nc]);
rhs = y;

% have we been given a constraint on the sum of the residuals?
if ~isempty(prescription.SumResiduals)
    Meq = [Meq;sum(Mdes,1)];
    rhseq = [rhseq;prescription.SumResiduals + sum(rhs)];
end

% apply weights
W = prescription.Weights;
if ~isempty(W)
    W = W(:);
    if length(W)~=nx
        error('SLMENGINE:inconsistentweights','Weight vector is not the same length as data')
    end
    
    % scale weight vector
    W = sqrt(nx)*W./norm(W);
    
    W = spdiags(W,0,nx,nx);
    Mdes = W*Mdes;
    rhs = W*rhs;
end

% -------------------------------------
% build regularizer, only bother if at least 3 knots
if nc>2
    % use a simple second order finite difference
    dx1 = dx(1:(nc-2));
    dx2 = dx(2:(nc-1));
    fda = [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), ...
        -2./(dx2.*(dx1+dx2))];
    
    Mreg = zeros(nc-2,nc);
    for i=1:(nc-2);
        Mreg(i,i+[0 1 2]) = fda(i,:);
    end
    rhsreg = zeros(nc-2,1);
end
% scale the regularizer before we apply the
% regularization parameter.
Mreg = Mreg/norm(Mreg,1);

% -------------------------------------
% end conditions only apply for a piecewise linear
% function, IF the conditions are periodicity
if strcmp(prescription.EndConditions,'periodic')
    % set the first and last function values equal
    M = zeros(1,nc);
    M([1,end]) = [-1 1];
    Meq = [Meq;M];
    rhseq = [rhseq;0];
end

% -------------------------------------
% single point equality constraints
% left hand side
if ~isempty(prescription.LeftValue)
    M = zeros(1,nc);
    M(1) = 1;
    Meq = [Meq;M];
    rhseq = [rhseq;prescription.LeftValue];
end

% right hand side
if ~isempty(prescription.RightValue)
    M = zeros(1,nc);
    M(end) = 1;
    Meq = [Meq;M];
    rhseq = [rhseq;prescription.RightValue];
end

% left end slope
if ~isempty(prescription.LeftSlope)
    M = zeros(1,nc);
    M(1:2) = [-1 1]/dx(1);
    Meq = [Meq;M];
    rhseq = [rhseq;prescription.LeftSlope];
end

% Right end slope
if ~isempty(prescription.RightSlope)
    M = zeros(1,nc);
    M([nc-1, nc]) = [-1 1]/dx(end);
    Meq = [Meq;M];
    rhseq = [rhseq;prescription.RightSlope];
end

% force curve through an x-y pair
xy = prescription.XY;
if ~isempty(xy)
    n = size(xy,1);
    if any(xy(:,1)<knots(1)) || any(xy(:,1)>knots(end))
        error('SLMENGINE:improperconstraint','XY pairs to force the curve through must lie inside the knots')
    end
    
    [junk,ind] = histc(xy(:,1),knots); %#ok
    ind(ind==(nc+1))=nc;
    t = (xy(:,1) - knots(ind))./dx(ind);
    
    M = sparse(repmat((1:n)',1,2),[ind,ind+1],[1-t,t],n,nc);
    Meq = [Meq;M];
    rhseq = [rhseq;xy(:,2)];
end

% force slope at a point, or set of points
xyp = prescription.XYP;
if ~isempty(xyp)
    n = size(xyp,1);
    if any(xyp(:,1)<knots(1)) || any(xyp(:,1)>knots(end))
        error('SLMENGINE:improperconstraint','X-Y'' pairs to enforce the slope at must lie inside the knots')
    end
    
    [junk,ind] = histc(xyp(:,1),knots); %#ok
    ind(ind==(nc+1))=nc;
    
    M = sparse(repmat((1:n)',1,2),[ind,ind+1],(1./dx(ind))*[-1,1],n,nc);
    Meq = [Meq;M];
    rhseq = [rhseq;xyp(:,2)];
    
end

% -------------------------------------
% constant region(s)
if ~isempty(prescription.ConstantRegion)
    % there was at least one constant region specified
    cr = prescription.ConstantRegion;
    M = zeros(0,nc);
    n = 0;
    for i=1:size(cr,1)
        % enforce constancy between the given range limits
        for j=1:(nc-1)
            if (knots(j+1)>=cr(i,1)) && (knots(j)<cr(i,2))
                n=n+1;
                M(n,j+[0 1]) = [-1 1];
            end
        end
    end
    
    Meq = [Meq;M];
    rhseq = [rhseq;zeros(size(M,1),1)];
    
end

% -------------------------------------
% linear region(s)
% a linear region imply means that across any knot
% inside that region, the slopes must be constant
if ~isempty(prescription.LinearRegion) && (nc>2)
    % there was at least one linear region specified
    lr = prescription.LinearRegion;
    M = zeros(0,nc);
    n = 0;
    for i=1:size(lr,1)
        % enforce constancy of slope between the given range limits
        for j=1:(nc-2)
            if (knots(j+1)>=lr(i,1)) && (knots(j+1)<lr(i,2))
                n=n+1;
                M(n,j+[0 1 2]) = [-1 1 0]/dx(j) - [0 -1 1]/dx(j+1);
            end
        end
    end
    
    Meq = [Meq;M];
    rhseq = [rhseq;zeros(size(M,1),1)];
    
end

% -------------------------------------
% Integral equality constraint
if ~isempty(prescription.Integral)
    % Trapezoidal rule over the support.
    M = ([dx(:)',0] + [0,dx(:)'])/2;
    Meq = [Meq;M];
    rhseq = [rhseq;prescription.Integral];
end

% -------------------------------------
% single point inequality constraints
% left hand side minimum
if ~isempty(prescription.LeftMinValue)
    M = zeros(1,nc);
    M(1) = -1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;-prescription.LeftMinValue];
end
% left hand side maximum
if ~isempty(prescription.LeftMaxValue)
    M = zeros(1,nc);
    M(1) = 1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;prescription.LeftMinValue];
end

% right hand side min
if ~isempty(prescription.RightMinValue)
    M = zeros(1,nc);
    M(end) = -1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;-prescription.RightMinValue];
end

% right hand side max
if ~isempty(prescription.RightMaxValue)
    M = zeros(1,nc);
    M(end) = 1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;prescription.RightMaxValue];
end

% left end slope min
if ~isempty(prescription.LeftMinSlope)
    M = zeros(1,nc);
    M(1:2) = [1 -1]/dx(1);
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;-prescription.LeftMinSlope];
end

% left end slope max
if ~isempty(prescription.LeftMaxSlope)
    M = zeros(1,nc);
    M(1:2) = [-1 1]/dx(1);
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;prescription.LeftMaxSlope];
end

% right end slope min
if ~isempty(prescription.RightMinSlope)
    M = zeros(1,nc);
    M([nc-1,nc]) = [1 -1]/dx(end);
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;-prescription.RightMinSlope];
end

% right end slope max
if ~isempty(prescription.RightMaxSlope)
    M = zeros(1,nc);
    M([nc-1,nc]) = [-1 1]/dx(end);
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;prescription.RightMaxSlope];
end

% -------------------------------------
% overall inequalities
% global min value
if ~isempty(prescription.MinValue)
    M = -eye(nc,nc);
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;repmat(-prescription.MinValue,nc,1)];
end

% global max value
if ~isempty(prescription.MaxValue)
    M = eye(nc,nc);
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;repmat(prescription.MaxValue,nc,1)];
end

% global min slope
if ~isempty(prescription.MinSlope)
    M = full(spdiags((1./dx)*[1 -1],[0 1],nc-1,nc));
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;repmat(-prescription.MinSlope,nc-1,1)];
end

% global max slope
if ~isempty(prescription.MaxSlope)
    M = full(spdiags((1./dx)*[-1 1],[0 1],nc-1,nc));
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;repmat(prescription.MaxSlope,nc-1,1)];
end

% -------------------------------------
% SimplePeak and SimpleValley are really just composite
% properties. A peak at x == a is equivalent to a monotone
% increasing function for x<=a, and a monotone decreasing
% function for x>=a. Likewise a valley is just the opposite.
%
% So specifying a peak or a valley will cause any monotonicity
% specs to be ignored. slmset has already checked for this,
% and turned off any other monotonicity based properties.
sp = prescription.SimplePeak;
if isnumeric(sp) && ~isempty(sp)
    prescription.Increasing = [knots(1),sp];
    prescription.Decreasing = [sp,knots(end)];
end
sv = prescription.SimpleValley;
if isnumeric(sv) && ~isempty(sv)
    prescription.Decreasing = [knots(1), sv];
    prescription.Increasing = [sv,knots(end)];
end

% -------------------------------------
% monotonicity?
% increasing regions
incR = prescription.Increasing;
L=0;
if ischar(incR)
    if strcmp(incR,'on')
        L=L+1;
        mono(L).knotlist = 'all';
        mono(L).direction = 1;
        mono(L).range = [];
    end
elseif ~isempty(incR)
    for i=1:size(incR,1)
        L=L+1;
        mono(L).knotlist = []; %#ok
        mono(L).direction = 1; %#ok
        mono(L).range = sort(incR(i,:)); %#ok
    end
end
% decreasing regions
decR = prescription.Decreasing;
if ischar(decR)
    if strcmp(decR,'on')
        L=L+1;
        mono(L).knotlist = 'all';
        mono(L).direction = -1;
        mono(L).range = [];
    end
elseif ~isempty(decR)
    for i=1:size(decR,1)
        L=L+1;
        mono(L).knotlist = [];
        mono(L).direction = -1;
        mono(L).range = sort(decR(i,:));
    end
end
if L>0
    % there were at least some monotone regions specified
    M = zeros(0,nc);
    n = 0;
    for i=1:L
        if isempty(mono(L).range)
            % the entire range was specified to be monotone
            for j=1:(nc-1)
                n=n+1;
                M(n,j+[0 1]) = [1 -1]*mono(i).direction;
            end
        else
            % only enforce monotonicity between the given range limits
            for j=1:(nc-1)
                if (knots(j)<mono(i).range(2)) && ...
                        (knots(j+1)>=mono(i).range(1))
                    n=n+1;
                    M(n,j+[0 1]) = [1 -1]*mono(i).direction;
                end
            end
            
        end
    end
    
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;zeros(size(M,1),1)];
end

% -------------------------------------
% PositiveInflection and NegativeInflection are really just
% composite properties. A point of inflection at x == a is
% equivalent to a change in the sign of the curvature at
% x == a.
pinf = prescription.PositiveInflection;
if isnumeric(pinf) && ~isempty(pinf)
    prescription.ConcaveDown = [knots(1),pinf];
    prescription.ConcaveUp = [pinf,knots(end)];
end
ninf = prescription.NegativeInflection;
if isnumeric(ninf) && ~isempty(ninf)
    prescription.ConcaveUp = [knots(1),ninf];
    prescription.ConcaveDown = [ninf,knots(end)];
end

% -------------------------------------
% concave up or down?
% regions of postive curvature
cuR = prescription.ConcaveUp;
L=0;
if ischar(cuR)
    if strcmp(cuR,'on')
        L=L+1;
        curv(L).knotlist = 'all';
        curv(L).direction = 1;
        curv(L).range = [];
    end
elseif ~isempty(cuR)
    for i=1:size(cuR,1)
        L=L+1;
        curv(L).knotlist = []; %#ok
        curv(L).direction = 1; %#ok
        curv(L).range = sort(cuR(i,:)); %#ok
    end
end
% decreasing regions
cdR = prescription.ConcaveDown;
if ischar(cdR)
    if strcmp(cdR,'on')
        L=L+1;
        curv(L).knotlist = 'all';
        curv(L).direction = -1;
        curv(L).range = [];
    end
elseif ~isempty(cdR)
    for i=1:size(cdR,1)
        L=L+1;
        curv(L).knotlist = [];
        curv(L).direction = -1;
        curv(L).range = sort(cdR(i,:));
    end
end
if L>0
    % there were at least some regions with specified curvature
    M = zeros(0,nc);
    n = 0;
    for i=1:L
        if isempty(curv(L).range)
            % the entire range was specified to be
            % curved in some direction
            for j=2:(nc-1)
                n=n+1;
                M(n,j+[-1 0 1]) = ([-1 1 0]/dx(j-1) - [0 -1 1]/dx(j))*curv(i).direction;
            end
        else
            % only enforce curvature between the given range limits
            for j=2:(nc-1)
                if (knots(j)<curv(i).range(2)) && ...
                        (knots(j)>=curv(i).range(1))
                    n=n+1;
                    M(n,j+[-1 0 1]) = ([-1 1 0]/dx(j-1) - [0 -1 1]/dx(j))*curv(i).direction;
                end
            end
            
        end
    end
    
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;zeros(size(M,1),1)];
end

% -------------------------------------
% Error bar inequality constraints
if ~isempty(prescription.ErrorBar)
    % lower bounds
    Mineq = [Mineq;-Mdes];
    rhsineq = [rhsineq;-(y-prescription.ErrorBar(:,1))];
    
    % upper bounds
    Mineq = [Mineq;Mdes];
    rhsineq = [rhsineq;y+prescription.ErrorBar(:,2)];
end

% -------------------------------------
% Use envelope inequalities?
switch prescription.Envelope
    case 'off'
        % nothing to do in this case
    case 'supremum'
        % yhat - y >= 0
        Mineq = [Mineq;-Mdes];
        rhsineq = [rhsineq;-rhs];
    case 'infimum'
        % yhat - y <= 0
        Mineq = [Mineq;Mdes];
        rhsineq = [rhsineq;rhs];
end

% -------------------------------------
% scale equalities for unit absolute row sum
if ~isempty(Meq)
    rs = diag(1./sum(abs(Meq),2));
    Meq = rs*Meq;
    rhseq = rs*rhseq;
end
% scale inequalities for unit absolute row sum
if ~isempty(Mineq)
    rs = diag(1./sum(abs(Mineq),2));
    Mineq = rs*Mineq;
    rhsineq = rs*rhsineq;
end

% -------------------------------------
% now worry about the regularization. There are four
% possible cases.
% 1. We have a given regularization parameter
% 2. We have a given rmse that we wish to match
% 3. We must use cross validation to choose the parameter
RP = prescription.Regularization;
if (isnumeric(RP) && (RP>=0)) || ((ischar(RP)) && (strcmpi(RP,'smoothest')))
    % solve the problem using the given regularization parameter
    coef = solve_slm_system(RP,Mdes,rhs,Mreg,rhsreg, ...
        Meq,rhseq,Mineq,rhsineq,prescription);
elseif isnumeric(RP) && (RP<0)
    % we must match abs(RP) as the rmse.
    aim_rmse = abs(RP);
    fminbndoptions = optimset('fminbnd');
    RP = fminbnd(@match_rmse,-6,6,fminbndoptions, ...
        Mdes,rhs,Mreg,rhsreg,Meq,rhseq,Mineq,rhsineq, ...
        aim_rmse,prescription);
    
    % we logged the parameter in the optimization. undo that for
    % the final call
    RP = 10^RP;
    
    % do one final call to get the final coefficients
    coef = solve_slm_system(RP,Mdes,rhs,Mreg,rhsreg, ...
        Meq,rhseq,Mineq,rhsineq,prescription);
    
elseif ischar(RP)
    % its a cross validation problem to solve. drop out each point
    % in turn from the model, then use fminbnd to minimize the
    % predicted sum of squares at the missing points.
    fminbndoptions = optimset('fminbnd');
    % fminbndoptions.Display = 'iter';
    RP = fminbnd(@min_cv,-6,6,fminbndoptions, ...
        Mdes,rhs,Mreg,rhsreg,Meq,rhseq,Mineq,rhsineq,prescription);
    
    % we logged the parameter in the optimization. undo that for
    % the final call
    RP = 10^RP;
    
    % do one final call to get the final coefficients
    coef = solve_slm_system(RP,Mdes,rhs,Mreg,rhsreg, ...
        Meq,rhseq,Mineq,rhsineq,prescription);
end

% -------------------------------------
% unpack coefficients into the result structure
slm.form = 'slm';
slm.degree = 1;
slm.knots = knots;
slm.coef = coef;

% degrees of freedom available
slmstats.TotalDoF = nk;
slmstats.NetDoF = slmstats.TotalDoF - size(Meq,1);
% this function does all of the stats, stuffing into slmstats
slm.stats = modelstatistics(slmstats,Mdes,y,coef);

% ========================================================
% =============== piecewise cubic model ==================
% ========================================================
function slm = slmengine_cubic(x,y,prescription)
% fits a piecewise cubic shape prescriptive model

% check for inappropriate properties for a cubic model
property_check(prescription,'cubic')

% simple things about data first...
% slmengine has already made it a column vector,
% and ensured compatibility in length between
% x and y
nx = length(x);

% knots vector, dx
if length(prescription.Knots)==1
    % just given a number of knots
    [knots,nk] = chooseknots(prescription.Knots,nx,x);
else
    % we should check that the knots contain the data
    knots = sort(prescription.Knots(:));
    nk= length(knots);
    if (knots(1)>min(x)) || (knots(end)<max(x))
        error('SLMENGINE:inadequateknots',['Knots do not contain the data. Data range: ',num2str([min(x),max(x)])])
    end
end
dx = diff(knots);
if any(dx==0)
    error('SLMENGINE:indistinctknots','Knots must be distinct.')
end

% number of coefficients to estimate.
% a piecewise cubic Hermite has two coefficients
% at each knot.
nc = 2*nk;

% create empty inequality constraints arrays. also
% a rhs vector
Mineq = zeros(0,nc);
rhsineq = [];
Meq = zeros(0,nc);
rhseq = [];

% -------------------------------------
% build design matrix -
% first, bin the data - histc wll do it best
[junk,xbin] = histc(x,knots); %#ok
% any point which falls at the top end, is said to
% be in the last bin.
xbin(xbin==nk)=nk-1;

% build design matrix
t = (x - knots(xbin))./dx(xbin);
t2 = t.^2;
t3 = t.^3;
s2 = (1-t).^2;
s3 = (1-t).^3;

vals = [3*s2-2*s3 ; 3*t2-2*t3 ; ...
    -(s3-s2).*dx(xbin) ; (t3-t2).*dx(xbin)];

% the coefficients will be stored in two blocks,
% first nk function values, then nk derivatives.
Mdes = accumarray([repmat((1:nx)',4,1), ...
    [xbin;xbin+1;nk+xbin;nk+xbin+1]],vals,[nx,nc]);
rhs = y;

% have we been given a constraint on the sum of the residuals?
if ~isempty(prescription.SumResiduals)
    Meq = [Meq;sum(Mdes,1)];
    rhseq = [rhseq;prescription.SumResiduals + sum(rhs)];
end

% apply weights
W = prescription.Weights;
if ~isempty(W)
    W = W(:);
    if length(W)~=nx
        error('SLMENGINE:inconsistentweights','Weight vector is not the same length as data')
    end
    
    % scale weight vector
    W = sqrt(nx)*W./norm(W);
    
    W = spdiags(W,0,nx,nx);
    Mdes = W*Mdes;
    rhs = W*rhs;
end

% -------------------------------------
% Regularizer
% We are integrating the piecewise linear f''(x), as
% a quadratic form in terms of the (unknown) second
% derivatives at the knots.
Mreg=zeros(nk,nk);
Mreg(1,1:2)=[dx(1)/3 , dx(1)/6];
Mreg(nk,nk+[-1 0])=[dx(end)/6 , dx(end)/3];
for i=2:(nk-1)
    Mreg(i,i+[-1 0 1])=[dx(i-1)/6 , (dx(i-1)+dx(i))/3 , dx(i)/6];
end
% do a matrix square root. cholesky is simplest & ok since regmat is
% positive definite. this way we can write the quadratic form as:
%    s'*r'*r*s,
% where s is the vector of second derivatives at the knots.
Mreg=chol(Mreg);

% next, write the second derivatives as a function of the
% function values and first derivatives.   s = [sf,sd]*[f;d]
sf = zeros(nk,nk);
sd = zeros(nk,nk);
for i = 1:(nk-1)
    sf(i,i+[0 1]) = [-1 1].*(6/(dx(i)^2));
    sd(i,i+[0 1]) = [-4 -2]./dx(i);
end
sf(nk,nk+[-1 0]) = [1 -1].*(6/(dx(end)^2));
sd(nk,nk+[-1 0]) = [2 4]./dx(end);
Mreg = Mreg*[sf,sd];

% scale the regularizer before we apply the
% regularization parameter.
Mreg = Mreg/norm(Mreg,1);
rhsreg = zeros(nk,1);

% -------------------------------------
% C2 continuity across knots
if strcmp(prescription.C2,'on')
    MC2 = zeros(nk-2,nc);
    for i = 1:(nk-2)
        MC2(i,[i,i+1]) = [6 -6]./(dx(i).^2);
        MC2(i,[i+1,i+2]) = MC2(i,[i+1,i+2]) + [6 -6]./(dx(i+1).^2);
        
        MC2(i,nk+[i,i+1]) = [2 4]./dx(i);
        MC2(i,nk+[i+1,i+2]) = MC2(i,nk+[i+1,i+2]) + [4 2]./dx(i+1);
    end
    
    Meq = [Meq;MC2];
    rhseq = [rhseq;zeros(nk-2,1)];
end

% -------------------------------------
% end conditions (do nothing here if they are 'estimate')
if strcmp(prescription.EndConditions,'periodic')
    % set the first and last function values equal
    M = zeros(3,nc);
    M(1,[1,nk]) = [-1 1];
    % as well as the first derivatives
    M(2,[nk+1,2*nk]) = [-1 1];
    % and the second derivatives
    M(3,[1 2]) = [-6 6]/(dx(1).^2);
    M(3,[nk+1,nk+2]) = [-4 -2]/dx(1);
    M(3,[nk-1,nk]) = [-6 6]/(dx(end).^2);
    M(3,[2*nk-1,2*nk]) = [-2 -4]/dx(end);
    
    Meq = [Meq;M];
    rhseq = [rhseq;zeros(3,1)];
elseif strcmp(prescription.EndConditions,'not-a-knot') && nk==2
    % third derivative is continuous across 2nd knot
    M = zeros(1,nc);
    
    M(1,[1 2]) = [12 -12]/(dx(1).^3);
    M(1,[2 3]) = M(1,2:3) - [12 -12]/(dx(2).^3);
    
    M(1,nk+[1 2]) = [6 6]/(dx(1).^2);
    M(1,nk+[2 3]) = M(1,nk+[2 3]) - [6 6]/(dx(2).^2);
    
    Meq = [Meq;M];
    rhseq = [rhseq;0];
elseif strcmp(prescription.EndConditions,'notaknot') && nk>2
    % third derivative is continuous at 2nd & penultimate knots
    M = zeros(2,nc);
    
    M(1,[1 2]) = [12 -12]/(dx(1).^3);
    M(1,[2 3]) = M(1,2:3) - [12 -12]/(dx(2).^3);
    M(1,nk+[1 2]) = [6 6]/(dx(1).^2);
    M(1,nk+[2 3]) = M(1,nk+[2 3]) - [6 6]/(dx(2).^2);
    
    M(2,nk + [-2 -1]) = [12 -12]/(dx(end - 1).^3);
    M(2,nk + [-1 0]) = M(2,nk + [-1 0]) - [12 -12]/(dx(end).^3);
    M(2,2*nk+[-2 -1]) = [6 6]/(dx(end-1).^2);
    M(2,2*nk+[-1 0]) = M(2,2*nk+[-1 0]) - [6 6]/(dx(end).^2);
    
    Meq = [Meq;M];
    rhseq = [rhseq;0;0];
elseif strcmp(prescription.EndConditions,'natural')
    % second derivative at each end is zero
    M = zeros(2,nc);
    
    M(1,[1 2]) = [-6 6]/(dx(1).^2);
    M(1,nk+[1 2]) = [-4 -2]/dx(1);
    
    M(2,nk+[-1 0]) = [6 -6]/(dx(end).^2);
    M(2,2*nk+[-1 0]) = [2 4]/dx(end);
    
    Meq = [Meq;M];
    rhseq = [rhseq;0;0];
end

% -------------------------------------
% single point equality constraints at a knot
% left hand side
if ~isempty(prescription.LeftValue)
    M = zeros(1,nc);
    M(1,1) = 1;
    Meq = [Meq;M];
    rhseq = [rhseq;prescription.LeftValue];
end

% right hand side
if ~isempty(prescription.RightValue)
    M = zeros(1,nc);
    M(1,nk) = 1;
    Meq = [Meq;M];
    rhseq = [rhseq;prescription.RightValue];
end

% left end slope
if ~isempty(prescription.LeftSlope)
    M = zeros(1,nc);
    M(1,nk+1) = 1;
    Meq = [Meq;M];
    rhseq = [rhseq;prescription.LeftSlope];
end

% Right end slope
if ~isempty(prescription.RightSlope)
    M = zeros(1,nc);
    M(1,2*nk) = 1;
    Meq = [Meq;M];
    rhseq = [rhseq;prescription.RightSlope];
end

% -------------------------------------
% force curve through an x-y pair or pairs
xy = prescription.XY;
if ~isempty(xy)
    n = size(xy,1);
    if any(xy(:,1)<knots(1)) || any(xy(:,1)>knots(end))
        error('SLMENGINE:improperconstraint','XY pairs to force the curve through must lie inside the knots')
    end
    
    [junk,ind] = histc(xy(:,1),knots); %#ok
    ind(ind==(nk))=nk-1;
    
    t = (xy(:,1) - knots(ind))./dx(ind);
    t2 = t.^2;
    t3 = t.^3;
    s2 = (1-t).^2;
    s3 = (1-t).^3;
    
    vals = [3*s2-2*s3 ; 3*t2-2*t3 ; ...
        -(s3-s2).*dx(ind) ; (t3-t2).*dx(ind)];
    
    M = accumarray([repmat((1:n)',4,1), ...
        [ind;ind+1;nk+ind;nk+ind+1]],vals,[n,nc]);
    
    Meq = [Meq;M];
    rhseq = [rhseq;xy(:,2)];
end

% -------------------------------------
% force slope at a point, or a set of points
xyp = prescription.XYP;
if ~isempty(xyp)
    n = size(xyp,1);
    if any(xyp(:,1)<knots(1)) || any(xyp(:,1)>knots(end))
        error('SLMENGINE:improperconstraint','X-Y'' pairs to enforce the slope at must lie inside the knots')
    end
    
    [junk,ind] = histc(xyp(:,1),knots); %#ok
    ind(ind==(nk))=nk-1;
    
    t = (xyp(:,1) - knots(ind))./dx(ind);
    t2 = t.^2;
    s2 = (1-t).^2;
    s = (1-t);
    
    vals = [(-6*s+6*s2)./dx(ind) ; ...
        (6*t-6*t2)./dx(ind) ; (3*s2-2*s) ; (3*t2-2*t)];
    
    M = accumarray([repmat((1:n)',4,1), ...
        [ind;ind+1;nk+ind;nk+ind+1]],vals,[n,nc]);
    
    Meq = [Meq;M];
    rhseq = [rhseq;xyp(:,2)];
end

% -------------------------------------
% force second derivative at a point, or set of points
xypp = prescription.XYPP;
if ~isempty(xypp)
    n = size(xypp,1);
    if any(xypp(:,1)<knots(1)) || any(xypp(:,1)>knots(end))
        error('SLMENGINE:improperconstraint','X-Y'' pairs to enforce y" at must lie inside the knots')
    end
    
    [junk,ind] = histc(xypp(:,1),knots); %#ok
    ind(ind==(nk))=nk-1;
    
    t = (xypp(:,1) - knots(ind))./dx(ind);
    s = (1-t);
    
    vals = [(6 - 12*s)./dx(ind).^2; ...
        (6 - 12*t)./dx(ind).^2 ; ...
        -(6*s - 2)./dx(ind) ; ...
        (6*t - 2)./dx(ind)];
    
    M = accumarray([repmat((1:n)',4,1), ...
        [ind;ind+1;nk+ind;nk+ind+1]],vals,[n,nc]);
    
    Meq = [Meq;M];
    rhseq = [rhseq;xypp(:,2)];
end

% -------------------------------------
% force third derivative at a point, or set of points
xyppp = prescription.XYPPP;
if ~isempty(xyppp)
    n = size(xyppp,1);
    if any(xyppp(:,1)<knots(1)) || any(xyppp(:,1)>knots(end))
        error('SLMENGINE:improperconstraint','X-Y'''''' pairs to enforce y'''''' at must lie inside the knots')
    end
    
    [junk,ind] = histc(xyppp(:,1),knots); %#ok
    ind(ind==(nk))=nk-1;
    
    vals = [12./dx(ind).^3; -12./dx(ind).^3 ; ...
        6./dx(ind).^2 ;6./dx(ind).^2];
    
    M = accumarray([repmat((1:n)',4,1), ...
        [ind;ind+1;nk+ind;nk+ind+1]],vals,[n,nc]);
    
    Meq = [Meq;M];
    rhseq = [rhseq;xyppp(:,2)];
end

% -------------------------------------
% constant region(s)
if ~isempty(prescription.ConstantRegion)
    % there was at least one constant region specified
    cr = prescription.ConstantRegion;
    M = zeros(0,nc);
    n = 0;
    for i=1:size(cr,1)
        % enforce constancy between the given range limits
        flag = true;
        for j=1:(nk-1)
            if (knots(j+1)>=cr(i,1)) && (knots(j)<cr(i,2))
                % f(j) == f(j+1)
                n=n+1;
                M(n,j+[0 1]) = [-1 1];
                
                % also zero f'' on this interval
                if flag
                    % only set this constraint for the first knot
                    % interval in a constant region
                    n=n+1;
                    M(n,j+[0 1 nk nk+1]) = [-6/(dx(j).^2) ; ...
                        6/(dx(j).^2) ; -4/dx(j) ; -2/dx(j)];
                    
                    flag = false;
                end
                
                n=n+1;
                M(n,j+[0 1 nk nk+1]) = [-6/(dx(j).^2) ; ...
                    6/(dx(j).^2) ; 2/dx(j) ; 4/dx(j)];
                
            end
        end
    end
    
    Meq = [Meq;M];
    rhseq = [rhseq;zeros(size(M,1),1)];
    
end

% -------------------------------------
% linear region(s)
% a linear region simply means that across any knot
% inside that region, the slopes must be constant
if ~isempty(prescription.LinearRegion) && (nc>2)
    % there was at least one linear region specified
    lr = prescription.LinearRegion;
    M = zeros(0,nc);
    n = 0;
    for i=1:size(lr,1)
        % enforce constancy of slope between the given range limits
        flag = true;
        for j=1:(nk-1)
            if (knots(j+1)>=lr(i,1)) && (knots(j)<lr(i,2))
                % f'(j) == f'(j+1)
                n=n+1;
                M(n,j+nk+[0 1]) = [-1 1];
                
                % also zero f'' at each end of this interval
                if flag
                    % only set this constraint for the first knot
                    % interval in a linear region
                    n=n+1;
                    M(n,j+[0 1 nk nk+1]) = [-6/(dx(j).^2) ; ...
                        6/(dx(j).^2) ; -4/dx(j) ; -2/dx(j)];
                    
                    flag = false;
                end
                n=n+1;
                M(n,j+[0 1 nk nk+1]) = [6/(dx(j).^2) ; ...
                    -6/(dx(j).^2) ; 2/dx(j) ; 4/dx(j)];
                
            end
        end
    end
    
    Meq = [Meq;M];
    rhseq = [rhseq;zeros(size(M,1),1)];
end

% -------------------------------------
% Integral equality constraint
if ~isempty(prescription.Integral)
    % Simpson's rule over the support will be exact,
    % evaluating at each midpoint. Don't forget that
    % the effective "stepsize" is actually dx(i)/2
    M = zeros(1,nc);
    for i = 1:(nk-1)
        M(1,i+[0 1]) = M(1,i+[0 1]) + dx(i)/6;
        % the midpoint
        M(1,i+[0 1 nk nk+1]) = M(1,i+[0 1 nk nk+1]) + ...
            [1/2 , 1/2 , (1/8).*dx(i) , (-1/8).*dx(i)]*4*dx(i)/6;
    end
    
    Meq = [Meq;M];
    rhseq = [rhseq;prescription.Integral];
end

% -------------------------------------
% single point inequality constraints
% left hand side minimum
if ~isempty(prescription.LeftMinValue)
    M = zeros(1,nc);
    M(1,1) = -1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;-prescription.LeftMinValue];
end

% left hand side maximum
if ~isempty(prescription.LeftMaxValue)
    M = zeros(1,nc);
    M(1,1) = 1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;prescription.LeftMaxValue];
end

% right hand side min
if ~isempty(prescription.RightMinValue)
    M = zeros(1,nc);
    M(1,nk) = -1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;-prescription.RightMinValue];
end

% right hand side max
if ~isempty(prescription.RightMaxValue)
    M = zeros(1,nc);
    M(1,nk) = 1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;prescription.RightMaxValue];
end

% left end slope min
if ~isempty(prescription.LeftMinSlope)
    M = zeros(1,nc);
    M(1,nk+1) = -1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;-prescription.LeftMinSlope];
end

% left end slope max
if ~isempty(prescription.LeftMaxSlope)
    M = zeros(1,nc);
    M(1,nk+1) = 1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;prescription.LeftMaxSlope];
end

% right end slope min
if ~isempty(prescription.RightMinSlope)
    M = zeros(1,nc);
    M(1,nc) = -1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;-prescription.RightMinSlope];
end

% right end slope max
if ~isempty(prescription.RightMaxSlope)
    M = zeros(1,nc);
    M(1,nc) = 1;
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;prescription.RightMaxSlope];
end

% -------------------------------------
% overall inequalities
% global min & max value
if ~isempty(prescription.MinValue) || ~isempty(prescription.MaxValue)
    % sample points for necessary conditions are chebychev
    % nodes. No good reason here, except that they are
    % as good a set of points as any other.
    tm = [0.017037; 0.066987; 0.14645; 0.25; 0.37059; ...
        0.5; 0.62941; 0.75; 0.85355; 0.93301; 0.98296];
    nsamp = length(tm);
    
    ntot = nk+(nk-1)*nsamp;
    Mmin = zeros(ntot,nc);
    Mmax = Mmin;
    
    % constrain values at knots
    if ~isempty(prescription.MinValue)
        Mmin(1:nk,1:nk) = -eye(nk,nk);
    end
    if ~isempty(prescription.MaxValue)
        Mmax(1:nk,1:nk) = eye(nk,nk);
    end
    
    % and intermediate sample points
    t2 = tm.^2;
    t3 = tm.^3;
    s2 = (1-tm).^2;
    s3 = (1-tm).^3;
    vals = [3*s2-2*s3 , 3*t2-2*t3 , -s3+s2 , t3-t2];
    
    for j = 1:(nk-1)
        if ~isempty(prescription.MinValue)
            Mmin((1:nsamp) + (j-1)*nsamp + nk , j+[0 1 nk nk+1]) = ...
                -vals*diag([1 1 dx(j) dx(j)]);
        end
        if ~isempty(prescription.MaxValue)
            Mmax((1:nsamp) + (j-1)*nsamp + nk , j+[0 1 nk nk+1]) = ...
                vals*diag([1 1 dx(j) dx(j)]);
        end
    end
    
    if ~isempty(prescription.MinValue)
        Mineq = [Mineq;Mmin];
        rhsineq = [rhsineq;repmat(-prescription.MinValue,ntot,1)];
    end
    if ~isempty(prescription.MaxValue)
        Mineq = [Mineq;Mmax];
        rhsineq = [rhsineq;repmat(prescription.MaxValue,ntot,1)];
    end
end

% global min or max slope
if ~isempty(prescription.MinSlope) || ~isempty(prescription.MaxSlope)
    % its a transformed monotonicity problem. See the
    % explanation provided below for sufficient conditions
    % for monotonicity. The trick here is to build them
    % for a minimum slope. Do that by effectively building
    % an offset into the slope. I.e., assume that theta is the
    % minimum slope allowed. Implicitly subtract the line
    % theta*x from the curve, then use the standard linear
    % system derived below for monotonicity. This will yield
    % a system of linear constraints that are sufficient for
    % a minimum (or maximum) slope of theta.
    
    %
    MS = zeros(5*(nk - 1) + nk,nc);
    rhsMS = zeros(5*(nk-1) + nk,1);
    
    % constrain the derivatives at each knot
    for n = 1:nk
        MS(n,nk + n) = -1;
    end
    rhsMS(1:nk) = -1;
    
    n = nk + 1;
    for j = 1:(nk - 1)
        % The constraints applied for an increasing function
        % are (in a form that lsqlin will like), i.e., A*x <= b
        %
        %    -delta + theta  <= 0
        %
        % (etc.)
        MS(n,j + [0 1]) = [1 -1];
        rhsMS(n) = -dx(j);
        
        MS(n + 1,j + [0, 1, nk,nk + 1]) = [3, -3, [-1, 1]*dx(j)];
        rhsMS(n + 1) = -3*dx(j);
        
        MS(n + 2,j + [0, 1, nk,nk + 1]) = [3, -3, [1, -1]*dx(j)];
        rhsMS(n + 2) = -3*dx(j);
        
        MS(n + 3,j + [0, 1, nk,nk + 1]) = [9, -9, [1, 2]*dx(j)];
        rhsMS(n + 3) = -6*dx(j);
        
        MS(n + 4,j + [0, 1, nk,nk + 1]) = [9, -9, [2, 1]*dx(j)];
        rhsMS(n + 4) = -6*dx(j);
        
        n = n + 5;
        
    end
    
    % fold these inequalities into any others
    if ~isempty(prescription.MinSlope)
        theta = prescription.MinSlope;
        
        Mineq = [Mineq;MS];
        rhsineq = [rhsineq;theta*rhsMS];
    end
    
    if ~isempty(prescription.MaxSlope)
        % the max slope case is as if we have negated
        % the min slope case.
        theta = prescription.MaxSlope;
        
        Mineq = [Mineq;-MS];
        rhsineq = [rhsineq;-theta*rhsMS];
    end
    
end

% -------------------------------------
% SimplePeak and SimpleValley are really just composite
% properties. A peak at x == a is equivalent to a monotone
% increasing function for x<=a, and a monotone decreasing
% function for x>=a. Likewise a valley is just the opposite.
%
% It is best if a peak or valley occurs at a knot.
%
% So specifying a peak or a valley will cause any monotonicity
% specs to be ignored. slmset has already checked for this,
% and turned off any other monotonicity based properties.
sp = prescription.SimplePeak;
if isnumeric(sp) && ~isempty(sp)
    prescription.Increasing = [knots(1),sp];
    prescription.Decreasing = [sp,knots(end)];
end
sv = prescription.SimpleValley;
if isnumeric(sv) && ~isempty(sv)
    prescription.Decreasing = [knots(1), sv];
    prescription.Increasing = [sv,knots(end)];
end

% -------------------------------------
% monotonicity?
% increasing regions
totalmonotoneintervals = 0;
incR = prescription.Increasing;
L=0;
if ischar(incR)
    if strcmp(incR,'on')
        L=L+1;
        mono(L).knotlist = [1,nk];
        mono(L).direction = 1;
        totalmonotoneintervals = totalmonotoneintervals + nk - 1;
    end
elseif ~isempty(incR)
    for i=1:size(incR,1)
        L=L+1;
        mono(L).knotlist = [max(1,sum(knots <= incR(i,1))), ...
            min(nk,1 + sum(knots < incR(i,2)))]; %#ok
        mono(L).direction = 1; %#ok
        totalmonotoneintervals = totalmonotoneintervals + diff(mono(L).knotlist);
    end
end

% decreasing regions
decR = prescription.Decreasing;
if ischar(decR)
    if strcmp(decR,'on')
        L=L+1;
        mono(L).knotlist = [1,nk];
        mono(L).direction = -1;
        totalmonotoneintervals = totalmonotoneintervals + nk - 1;
    end
elseif ~isempty(decR)
    for i=1:size(decR,1)
        L=L+1;
        mono(L).knotlist = [max(1,sum(knots <= decR(i,1))), ...
            min(nk,1 + sum(knots < decR(i,2)))];
        mono(L).direction = -1;
        totalmonotoneintervals = totalmonotoneintervals + diff(mono(L).knotlist);
    end
end

if L>0
    % there were at least some monotone regions specified
    M = zeros(7*totalmonotoneintervals,nc);
    n = 0;
    for i=1:L
        for j = mono(i).knotlist(1):(mono(i).knotlist(2) - 1)
            % the function must be monotone between
            % knots j and j + 1. The direction over
            % that interval is specified. The constraint
            % system used comes from Fritsch & Carlson, see here:
            %
            % http://en.wikipedia.org/wiki/Monotone_cubic_interpolation
            %
            % Define delta = (y(i+1) - y(i))/(x(i+1) - x(i))
            % Thus delta is the secant slope of the curve across
            % a knot interval. Further, define alpha and beta as
            % the ratio of the derivative at each end of an
            % interval to the secant slope.
            %
            %  alpha = d(i)/delta
            %  beta = d(i+1)/delta
            %
            % Then we have an elliptically bounded region in the
            % first quadrant that defines the set of monotone cubic
            % segments. We cannot define that elliptical region
            % using a set of linear constraints. However, by use
            % of a system of 7 linear constraints, we can form a
            % set of sufficient conditions such that the curve
            % will be monotone. There will be some few cubic
            % segments that are actually monotone, yet lie outside
            % of the linear system formed. This is acceptable,
            % as our linear approximation here is a sufficient
            % one for monotonicity, although not a necessary one.
            % It merely says that the spline may be slightly over
            % constrained, i.e., slightly less flexible than is
            % absolutely necessary. (So?)
            %
            % The 7 constraints applied for an increasing function
            % are (in a form that lsqlin will like):
            %
            %    -delta          <= 0
            %    -alpha          <= 0
            %    -beta           <= 0
            %    -alpha + beta   <= 3
            %     alpha - beta   <= 3
            %     alpha + 2*beta <= 9
            %   2*alpha + beta   <= 9
            %
            % Multiply these inequalities by (y(i+1) - y(i)) to
            % put them into a linear form.
            M(n + 1,j+[0 1]) = [1 -1]*mono(i).direction;
            M(n + 2,nk + j) = -mono(i).direction;
            M(n + 3,nk + j + 1) = -mono(i).direction;
            
            M(n + 4,j + [0, 1, nk,nk + 1]) = mono(i).direction*[3, -3, [-1, 1]*dx(j)];
            M(n + 5,j + [0, 1, nk,nk + 1]) = mono(i).direction*[3, -3, [1, -1]*dx(j)];
            M(n + 6,j + [0, 1, nk,nk + 1]) = mono(i).direction*[9, -9, [1, 2]*dx(j)];
            M(n + 7,j + [0, 1, nk,nk + 1]) = mono(i).direction*[9, -9, [2, 1]*dx(j)];
            
            n = n + 7;
        end
    end
    
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;zeros(size(M,1),1)];
end

% -------------------------------------
% PositiveInflection and NegativeInflection are really just
% composite properties. A point of inflection at x == a is
% equivalent to a change in the sign of the curvature at
% x == a.
pinf = prescription.PositiveInflection;
if isnumeric(pinf) && ~isempty(pinf)
    prescription.ConcaveDown = [knots(1),pinf];
    prescription.ConcaveUp = [pinf,knots(end)];
end
ninf = prescription.NegativeInflection;
if isnumeric(ninf) && ~isempty(ninf)
    prescription.ConcaveUp = [knots(1),ninf];
    prescription.ConcaveDown = [ninf,knots(end)];
end

% -------------------------------------
% concave up or down?
% regions of positive curvature
cuR = prescription.ConcaveUp;
L=0;
if ischar(cuR)
    if strcmp(cuR,'on')
        L=L+1;
        curv(L).knotlist = 'all';
        curv(L).direction = 1;
        curv(L).range = [];
    end
elseif ~isempty(cuR)
    for i=1:size(cuR,1)
        L=L+1;
        curv(L).knotlist = []; %#ok
        curv(L).direction = 1; %#ok
        curv(L).range = sort(cuR(i,:)); %#ok
    end
end
% negative curvature regions
cdR = prescription.ConcaveDown;
if ischar(cdR)
    if strcmp(cdR,'on')
        L=L+1;
        curv(L).knotlist = 'all';
        curv(L).direction = -1;
        curv(L).range = [];
    end
elseif ~isempty(cdR)
    for i=1:size(cdR,1)
        L=L+1;
        curv(L).knotlist = [];
        curv(L).direction = -1;
        curv(L).range = sort(cdR(i,:));
    end
end
if L>0
    % there were at least some regions with specified curvature
    M = zeros(0,nc);
    n = 0;
    for i=1:L
        if isempty(curv(L).range)
            % the entire domain was specified to be
            % curved in some direction
            for j=1:(nk-1)
                n=n+1;
                M(n,j+[0 1]) = curv(i).direction*[6 -6]/(dx(j).^2);
                M(n,nk+j+[0 1]) = curv(i).direction*[4 2]/dx(j);
            end
            n=n+1;
            M(n,nk+[-1 0]) = curv(i).direction*[-6 6]/(dx(end).^2);
            M(n,2*nk+[-1 0]) = curv(i).direction*[-2 -4]/dx(end);
        else
            % only enforce curvature between the given range limits
            % do each knot first.
            for j=1:(nk-1)
                if (knots(j)<curv(i).range(2)) && ...
                        (knots(j)>=curv(i).range(1))
                    
                    n=n+1;
                    M(n,j+[0 1]) = curv(i).direction*[6 -6]/(dx(j).^2);
                    M(n,nk+j+[0 1]) = curv(i).direction*[4 2]/dx(j);
                end
            end
            
            % also constrain at the endpoints of the range
            curv(i).range = max(min(curv(i).range(:),knots(end)),knots(1));
            [junk,ind] = histc(curv(i).range,knots); %#ok
            ind(ind==(nk))=nk-1;
            
            t = (curv(i).range - knots(ind))./dx(ind);
            s = 1-t;
            
            for j = 1:numel(ind)
                M(n+j,ind(j)+[0 1 nk nk+1]) = -curv(i).direction* ...
                    [(6 - 12*s(j))./(dx(ind(j)).^2), (6 - 12*t(j))./(dx(ind(j)).^2) , ...
                    (2 - 6*s(j))./dx(ind(j)), (6*t(j) - 2)./dx(ind(j))];
            end
            
            n = n + numel(ind);
        end
    end
    
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;zeros(size(M,1),1)];
end

% -------------------------------------
% Jerk function increasing or decreasing?
% Just a constraint on the sign of the third derivative.
switch prescription.Jerk
    case 'positive'
        fpppsign = 1;
    case 'negative'
        fpppsign = -1;
    otherwise
        fpppsign = 0;
end
if fpppsign
    % constrain the third derivative sign
    M = zeros(nk-1,nc);
    for i=1:(nk-1)
        % one constraint on every interval
        M(i,i+[0 1]) = fpppsign*[-12 12]/(dx(i).^3);
        M(i,nk+i+[0 1]) = fpppsign*[-6 -6]/dx(i).^2;
    end
    
    Mineq = [Mineq;M];
    rhsineq = [rhsineq;zeros(size(M,1),1)];
end

% -------------------------------------
% Error bar inequality constraints
if ~isempty(prescription.ErrorBar)
    % lower bounds
    Mineq = [Mineq;-Mdes];
    rhsineq = [rhsineq;-(y-prescription.ErrorBar(:,1))];
    
    % upper bounds
    Mineq = [Mineq;Mdes];
    rhsineq = [rhsineq;y+prescription.ErrorBar(:,2)];
end

% -------------------------------------
% Use envelope inequalities?
switch prescription.Envelope
    case 'off'
        % nothing to do in this case
    case 'supremum'
        % yhat - y >= 0
        Mineq = [Mineq;-Mdes];
        rhsineq = [rhsineq;-rhs];
    case 'infimum'
        % yhat - y <= 0
        Mineq = [Mineq;Mdes];
        rhsineq = [rhsineq;rhs];
end

% -------------------------------------
% scale equalities for unit absolute row sum
if ~isempty(Meq)
    rs = diag(1./sum(abs(Meq),2));
    Meq = rs*Meq;
    rhseq = rs*rhseq;
end
% scale inequalities for unit absolute row sum
if ~isempty(Mineq)
    rs = diag(1./sum(abs(Mineq),2));
    Mineq = rs*Mineq;
    rhsineq = rs*rhsineq;
end

% -------------------------------------
% now worry about the regularization. There are three
% possible cases.
% 1. We have a given regularization parameter
% 2. We have a given rmse that we wish to match
% 3. We must use cross validation to choose the parameter
RP = prescription.Regularization;
if (isnumeric(RP) && (RP>=0)) || ((ischar(RP)) && (strcmpi(RP,'smoothest')))
    % solve the problem using the given regularization parameter
    coef = solve_slm_system(RP,Mdes,rhs,Mreg,rhsreg, ...
        Meq,rhseq,Mineq,rhsineq,prescription);
elseif isnumeric(RP) && (RP<0)
    % we must match abs(RP) as the rmse.
    aim_rmse = abs(RP);
    fminbndoptions = optimset('fminbnd');
    RP = fminbnd(@match_rmse,-6,6,fminbndoptions, ...
        Mdes,rhs,Mreg,rhsreg,Meq,rhseq,Mineq,rhsineq, ...
        aim_rmse,prescription);
    
    % we logged the parameter in the optimization. undo that for
    % the final call
    RP = 10^RP;
    
    % do one final call to get the final coefficients
    coef = solve_slm_system(RP,Mdes,rhs,Mreg,rhsreg, ...
        Meq,rhseq,Mineq,rhsineq,prescription);
    
elseif ischar(RP)
    % its a cross validation problem to solve. drop out each point
    % in turn from the model, then use fminbnd to minimize the
    % predicted sum of squares at the missing points.
    fminbndoptions = optimset('fminbnd');
    % fminbndoptions.Display = 'iter';
    RP = fminbnd(@min_cv,-6,6,fminbndoptions, ...
        Mdes,rhs,Mreg,rhsreg,Meq,rhseq,Mineq,rhsineq,prescription);
    
    % we logged the parameter in the optimization. undo that for
    % the final call
    RP = 10^RP;
    
    % do one final call to get the final coefficients
    coef = solve_slm_system(RP,Mdes,rhs,Mreg,rhsreg, ...
        Meq,rhseq,Mineq,rhsineq,prescription);
end

% -------------------------------------
% unpack coefficients into the result structure
slm.form = 'slm';
slm.degree = 3;
slm.knots = knots;
slm.coef = reshape(coef,nk,2);

% generate model statistics
slmstats.TotalDoF = 2*nk;
slmstats.NetDoF = slmstats.TotalDoF - size(Meq,1);
% this function does all of the stats, stuffing into slmstats
slm.stats = modelstatistics(slmstats,Mdes,y,coef);

% ========================================================
% ========= basic linear system solve ====================
% ========================================================
function [coef,lambda] = solve_slm_system(RP,Mdes,rhs,Mreg,rhsreg,Meq,rhseq,Mineq,rhsineq,prescription)
% solves the final linear system of equations for the model coefficients

if prescription.Verbosity > 1
    % Linear solve parameters
    disp('=========================================')
    disp('LINEAR SYSTEM SOLVER')
    disp(['Design matrix shape:    ',num2str(size(Mdes))])
    disp(['Regularizer shape:      ',num2str(size(Mreg))])
    disp(['Equality constraints:   ',num2str(size(Meq))])
    disp(['Inequality constraints: ',num2str(size(Mineq))])
end

% combine design matrix and regularizer
if strcmpi(prescription.Regularization,'smoothest')
    % find the smoothest possible solution
    Mdes = [0.000001*Mdes;Mreg];
    rhs = [0.000001*rhs;rhsreg];
else
    Mdes = [Mdes;RP*Mreg];
    rhs = [rhs;RP*rhsreg];
end

if prescription.Verbosity > 1
    disp(' ')
    disp(['Condition number of the regression: ',num2str(cond(Mdes))])
    disp(' ')
end

% -------------------------------------
% solve the linear system
if isempty(Mineq) && isempty(Meq)
    coef = Mdes\rhs;
    lambda.eqlin=[];
    lambda.ineqlin=[];
elseif isempty(Mineq)
    % with no inequality constraints, lse is faster than
    % is lsqlin. This also allows the use of slm when
    % the optimization toolbox is not present if there
    % are no inequality constraints.
    coef = lse(Mdes,rhs,full(Meq),rhseq);
    
else
    % use lsqlin. first, set the options
    options = optimset('lsqlin');
    if prescription.Verbosity > 1
        options.Display = 'final';
    else
        options.Display = 'off';
    end
    % the Largescale solver will not allow general constraints,
    % either equality or inequality
    options.LargeScale='off';
    
    % and solve
    [coef,junk,junk,exitflag,junk,lambda] = ...
        lsqlin(Mdes,rhs,Mineq,rhsineq,Meq,rhseq,[],[],[],options); %#ok
    
    % was there a feasible solution?
    if exitflag == -2
        coef = nan(size(coef));
    end
    
end

% ========================================================
% ========== chose RP to hit aim rmse ====================
% ========================================================
function delta_rmse = match_rmse(RP,Mdes,rhs,Mreg,rhsreg,Meq,rhseq,Mineq,rhsineq,aim_rmse,prescription)
% returns delta rmse for a given regularization parameter

% the regularization parameter has been logged.
% exponentiate for the solve
RP = 10^RP;

% solve the system at the given RP
coef = solve_slm_system(RP,Mdes,rhs,Mreg,rhsreg, ...
    Meq,rhseq,Mineq,rhsineq,prescription);

% compute rmse
rmse = sqrt(mean((rhs - Mdes*coef).^2));

% delta rmse error from aim
delta_rmse = abs(aim_rmse - rmse);


% ========================================================
% ======= chose RP to minimize cv prediction error =======
% ========================================================
function cvpred = min_cv(RP,Mdes,rhs,Mreg,rhsreg,Meq,rhseq,Mineq,rhsineq,prescription)
% returns delta rmse for a given regularization parameter

% the regularization parameter has been logged.
% exponentiate for the solve
RP = 10^RP;

% solve the system at the given RP, dropping out one point
% at a time from the data.
n = size(Mdes,1);
uselist = true(n,1);
cvpred = 0;
for i=1:n
    ul = uselist;
    ul(i) = false;
    coef = solve_slm_system(RP,Mdes(ul,:),rhs(ul),Mreg, ...
        rhsreg,Meq,rhseq,Mineq,rhsineq,prescription);
    pred = Mdes(i,:)*coef;
    cvpred = cvpred + (pred - rhs(i)).^2;
end


% ========================================================
% ========= choose knot placement for free knots =========
% ========================================================
function rss = free_knot_obj(intknots,x,y,prescription)
% returns residual sum of squares for fmincon to work on

% insert the interior knots from the optimizer
prescription.Knots(2:(end-1)) = intknots;
slm = slmengine(x,y,prescription);

% get predictions, form residuals, square, and sum
switch prescription.Result
    case 'slm'
        rss = sum((slmeval(x,slm,0) - y).^2);
    case 'pp'
        rss = sum((ppval(slm,x) - y).^2);
end

% ========================================================
% ========= set scaling as necessary =========
% ========================================================
function [ytrans,prescription] = scaleproblem(x,y,prescription)
% chooses an appropriate scaling for the problem

if strcmp(prescription.Scaling,'on')
    % scale y so that the minimum value is 1, and the maximum value 2.
    
    
    
    
    
    
    
    
else
    % no scaling is done
    prescription.yshift = 0;
    prescription.yscale = 1;
    ytrans = y;
end


% ========================================================
% ========= chooseknots =========
% ========================================================
function [knots,nk] = chooseknots(K,n,x)
% choose every K'th data point as a knot

if K > 0
    % just given a number of knots, equally spaced
    nk = K;
    knots = linspace(min(x),max(x),nk)';
    knots(end)=max(x); % just to make sure
else
    % we need every abs(K)'th knot
    
    if (K == 0) || (K >= n)
        error('SLMENGINE:knots','Every K''th data point was specified, but K was too large or zero')
    end
    
    Kind = 1:abs(K):n;
    % in case K was not an integer, or did not
    % go to the very end.
    if (n - Kind(end)) < (K/2)
        % expand the last knot interval
        Kind(end) = n;
    else
        % append one extra knot for the last data point
        Kind(end + 1) = n;
    end
    % rounding makes them into indices
    Kind = round(Kind);
    % in case of replicate data points
    x = sort(x);
    knots = unique(x(Kind));
    nk = numel(knots);
end



% ========================================================
% ========= model statistics =========
% ========================================================
function slmstats = modelstatistics(slmstats,Mdes,y,coef)
% generate model statistics, stuffing them into slmstats

% residuals, as yhat - y
resids = Mdes*coef - y;

% RMSE: Root Mean Squared Error
slmstats.RMSE = sqrt(mean(resids.^2));

% R-squared
slmstats.R2 = 1 - sum(resids.^2)./sum((y - mean(y)).^2);

% adjusted R^2
ndata = numel(y);
slmstats.R2Adj = 1 - (1-slmstats.R2)*(ndata - 1)./(ndata - slmstats.NetDoF);

% range of the errors, min to max, as yhat - y
slmstats.ErrorRange = [min(resids),max(resids)];

% compute the 25% and 75% points (quartiles) of the residuals
% (This is consistent with prctile, from the stats TB.)
resids = sort(resids.');
ind = 0.5 + ndata*[0.25 0.75];
f = ind - floor(ind);
ind = min(ndata - 1,max(1,floor(ind)));
slmstats.Quartiles = resids(ind).*(1-f) + resids(ind+1).*f;


% ========================================================
% ========== check for inappropriate properties ==========
% ========================================================
function property_check(prescription,model_degree)
% issues warning messages for inappropriate properties
% for the given model degree.

switch model_degree
    case {3 'cubic'}
        % no properties flagged for cubic models
        
    case {1 'linear'}
        % only (some) curvature properties are flagged
        % for the linear model
        
        if ~isempty(prescription.Jerk)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a linear model: Jerk')
        end
        
        if ~isempty(prescription.XYPP)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a linear model: XYPP')
        end
        
        if ~isempty(prescription.XYPPP)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a linear model: XYPPP')
        end
        
    case {0 'constant'}
        % curvature properties are flagged for
        % the constant model, as well as most slope
        % related properties
        
        if ~isempty(prescription.Jerk)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: Jerk')
        end
        
        if ~isempty(prescription.LeftMaxSlope)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: LeftMaxSlope')
        end
        
        if ~isempty(prescription.LeftMinSlope)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: LeftMinSlope')
        end
        
        if ~isempty(prescription.LeftSlope)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: LeftSlope')
        end
        
        if ~isempty(prescription.LinearRegion)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: LinearRegion')
        end
        
        if ~isempty(prescription.MaxSlope)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: MaxSlope')
        end
        
        if ~isempty(prescription.MinSlope)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: MinSlope')
        end
        
        if ~isempty(prescription.NegativeInflection)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: NegativeInflection')
        end
        
        if ~isempty(prescription.PositiveInflection)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: PositiveInflection')
        end
        
        if ~isempty(prescription.RightMaxSlope)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: RightMaxSlope')
        end
        
        if ~isempty(prescription.RightMinSlope)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: RightMinSlope')
        end
        
        if ~isempty(prescription.RightSlope)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: RightSlope')
        end
        
        if ~isempty(prescription.XYP)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: XYP')
        end
        
        if ~isempty(prescription.XYPP)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: XYPP')
        end
        
        if ~isempty(prescription.XYPPP)
            warning('SLMENGINE:ignoredconstraint', ...
                'Property ignored for a constant model: XYPPP')
        end
        
end




