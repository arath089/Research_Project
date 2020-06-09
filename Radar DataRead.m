fs=17;
TWO_PERSON=false;
BREATHING_BELTS=false;
EMPTY_ROOM_AVAILABLE =false;
TWO_RADARS=false;


distanceStart = 0.25; %meters
distanceEnd = 9.25; %meters
timeStart = 0*1000; %ms
timeEnd = 1000000*1000; %ms
radarResolution = 0.0522; %3.90625/1000; %radar resolution in meters
num_fast_chan=180;
distances = (0:num_fast_chan-2)*radarResolution;
breathingRateRef=24;

D=2.87; %distance between 2 radars
SHOW_RAW_DATA_IN_ALL_BINS = true;
SHOW_ONLY_MAX_VAR_BIN = true;


REMOVE_CLOSE_BINS = true;
REMOVAL_THRESHOLD = ceil(0.75/radarResolution); %Remove 0.75m in front of radar
FOLDER_ALL = '/home/user1/Documents/AES/Data Collection/20190410';%'/home/user1/Documents/AES/Data Collection/18-03-19/twoperson/nobacksupport/'; %'/home/user1/Documents/AES/Data Collection/28022019/twoperson/sitting';%'C:\Users\Miodrag Bolic\Dropbox\Backup\2018 2019\AES\Data Collection\chairtiledtoright'; 
if isunix()
    FOLDER_ALL = strcat(FOLDER_ALL,'/'); %For Linux systems
else
    FOLDER_ALL = strcat(FOLDER_ALL,'\'); %For Windows systems
end

%%
%%% please use correct format when specifying folder location according to
%%% your operating system (Windows uses '/')    



UWB_FILE_NAME ='20190410-001639_radar_data_.csv'; %'radar1/0degree_2p_sit_rd1.csv';%'\rawdata_rdr1\LAB-5130-ABi-SIT-SIDE-BOOK_20190103-104525_1.csv'; %File name of radar data file 'LAB-5130-ABi-SIT-SIDE-BOOK_20181220-221049.csv';
%UWB_FILE_NAME2 = 'radar1/sitting_fcrd2 (2).csv';

%load radar data
radarRawData1a = csvread(strcat(FOLDER_ALL,UWB_FILE_NAME));
%radarRawData1b = csvread(strcat(FOLDER_ALL,UWB_FILE_NAME2));
radarRawData1=radarRawData1a;%[radarRawData1a; radarRawData1b];
radarDataTimes1 = radarRawData1(:,1);radarStartIndex1 = 1;  
radarRawData1=radarRawData1(radarStartIndex1:end,2:num_fast_chan);
  

radarDataNormalized1 =  abs(radarRawData1); %unwrap(angle(radarRawData(radarStartIndex:end,2:end))); %
radarDataNormalized1(:,1:REMOVAL_THRESHOLD)=0;

%% Radar 2


if TWO_RADARS==true
    UWB_FILE_NAME ='sit_front_30beats_1.652m_sec2.csv'; %File name of radar data file 'LAB-5130-ABi-SIT-SIDE-BOOK_20181220-221049.csv';
    %UWB_FILE_NAME2 ='radar2/sitting_fcrd2 (2).csv';
    radarResolution = 0.0522; %3.90625/1000; %radar resolution in meters

    %load radar data
    radarRawData2 = csvread(strcat(FOLDER_ALL,UWB_FILE_NAME));
    %radarRawData2b = csvread(strcat(FOLDER_ALL,UWB_FILE_NAME2));
    %radarRawData2=[radarRawData2a; radarRawData2b];
    radarDataTimes2 = radarRawData2(:,1);
    radarStartIndex2 = 1;    
    radarRawData2=radarRawData2(radarStartIndex2:end,2:num_fast_chan);
    radarDataNormalized2 =  abs(radarRawData2); %unwrap(angle(radarRawData(radarStartIndex:end,2:end))); %
    radarDataNormalized2(:,1:REMOVAL_THRESHOLD)=0;
end

%%
%Reading Breathing Belt data
if BREATHING_BELTS  
fsb=10;
% BR belt 1
UWB_FILE_NAME ='breathingbelt/zixiang/2per_sit_close.csv';
%UWB_FILE_NAME2 ='breathingbelt/abi/sitting_fcrd2 (2).csv';
A = csvread(strcat(FOLDER_ALL,UWB_FILE_NAME));
%Ab = csvread(strcat(FOLDER_ALL,UWB_FILE_NAME2));
%A=[Aa;Ab];
%A=csvread('/home/user1/Documents/AES/Data Collection/sample1202/twoperson/breathingbelt/sitting-15mapart/abi_sitting(1).csv');%'C:\Users\Miodrag Bolic\Dropbox\Backup\2018 2019\AES\Data Collection\chairtiledtoright\breathingbelt\GDX-RB 0K1000Z9\BREATHING_BELT_GDX-RB 0K1000Z9_20190103-104519.csv');
BR_belt1 = resample(A(:,2),length(radarDataNormalized1(:,1)),length(A(:,2)));

if TWO_PERSON == true
% BR belt 2
UWB_FILE_NAME ='breathingbelt/shan/2per_sit_close.csv';
%UWB_FILE_NAME2 ='breathingbelt/shan/sitting_fcrd2 (2).csv';
A1 = csvread(strcat(FOLDER_ALL,UWB_FILE_NAME));
%Ab = csvread(strcat(FOLDER_ALL,UWB_FILE_NAME2));
%A1=[Aa;Ab];
%A1=csvread('/home/user1/Documents/AES/Data Collection/sample1202/twoperson/breathingbelt/sitting-15mapart/shan_sitting (1).csv');
BR_belt2 = resample(A1(:,2),length(radarDataNormalized2(:,1)),length(A1(:,2)));
end
end
%% Empty room
if EMPTY_ROOM_AVAILABLE ==true
%UWB_FILE_NAME = '\empty_test\empty_radar1.csv';%'\rawdata_rdr2\LAB-5130-ABi-SIT-SIDE-BOOK_20190103-104525_1.csv'; %File name of radar data file 'LAB-5130-ABi-SIT-SIDE-BOOK_20181220-221049.csv';
%load radar data
%radarRawEmpty1 = csvread(strcat(FOLDER_ALL,UWB_FILE_NAME));
radarRawEmpty1 = csvread('/home/user1/Documents/AES/Data Collection/2019_04_06/empty_room1.csv');
radarDataTimesEmpty1 = radarRawEmpty1(:,1);
radarStartIndexEmpty1 = 1;    

radarDataEmptyNormalized1 =  abs(radarRawEmpty1(radarStartIndexEmpty1:end,2:num_fast_chan)); %unwrap(angle(radarRawData(radarStartIndex:end,2:end))); 
radarDataEmptyNormalized1(:,1:REMOVAL_THRESHOLD)=0;

x=mean(radarDataEmptyNormalized1);
y=repmat(x,length(radarDataNormalized1),1);
radarDataNormalized1=radarDataNormalized1-y;

empty1=radarRawEmpty1(radarStartIndexEmpty1:end,2:num_fast_chan);
x1=mean(empty1);
radarRawData1=radarRawData1-x1;

%UWB_FILE_NAME = '\empty_test\empty_radar2.csv';%'\rawdata_rdr2\LAB-5130-ABi-SIT-SIDE-BOOK_20190103-104525_1.csv'; %File name of radar data file 'LAB-5130-ABi-SIT-SIDE-BOOK_20181220-221049.csv';
%load radar data
%radarRawEmpty2 = csvread(strcat(FOLDER_ALL,UWB_FILE_NAME));

if TWO_RADARS==true
radarRawEmpty2 = csvread('/home/user1/Documents/AES/Data Collection/2703/1person/radar2/empty_rd2.csv');%('/home/user1/Documents/AES/Data Collection/28022019/empty/radar2/LAB-5130-empty-SIT-SIDE-BOOK_20190228-164528_.csv');

radarDataTimesEmpty2 = radarRawEmpty2(:,1);
radarStartIndexEmpty2 = 1;    

radarDataEmptyNormalized2 =  abs(radarRawEmpty2(radarStartIndexEmpty2:end,2:num_fast_chan)); %unwrap(angle(radarRawData(radarStartIndex:end,2:end))); 
radarDataEmptyNormalized2(:,1:REMOVAL_THRESHOLD)=0;



x=mean(radarDataEmptyNormalized2);
y=repmat(x,length(radarDataNormalized2),1);
radarDataNormalized2=radarDataNormalized2-y;

empty2=radarRawEmpty2(radarStartIndexEmpty2:end,2:num_fast_chan);
x2=mean(empty2);
radarRawData2=radarRawData2-x2;
%radarDataNormalized1_1=analyze_UWB_empty_room(radarDataNormalized1,radarDataEmptyNormalized1);
end
end

%% Detection of presence
Segment_length_pres = 2*fs;
Segment_Num_pres = fix(length(radarDataNormalized1(:,1))/Segment_length_pres);

for j = 1: Segment_Num_pres
   j 
   %[range_out1{j}, MTI1(j)] = analyze_UWB_MTI(radarDataNormalized1((j-1)*Segment_length_pres+1:Segment_ofset(j)+Segment_length-1_pres,:));
   %[range_out2{j}, MTI2(j)] = analyze_UWB_MTI(radarDataNormalized2((j-1)*Segment_length_pres+1:Segment_ofset(j)+Segment_length-1_pres,:));
  % [range_out1(j), MTI1(j)] = analyze_UWB_MTI(radarDataNormalized1((j-1)*Segment_length_pres+1:Segment_ofset(j)+Segment_length-1_pres,:));
  % [range_out2(j), MTI2(j)] = analyze_UWB_MTI(radarDataNormalized2((j-1)*Segment_length_pres+1:Segment_ofset(j)+Segment_length-1_pres,:));
%    figure, for i=1:10
%     hold on, plot(1:15,range_out1(:,i),'o')
%     end

end
% figure, plot(range_out1)
% hold on
% plot(range_out2)
% Movement indicator need to be defined

%% Breathing estimation
Segment_length = 15*fs;
Segment_overlap=0;
Segment_Num =fix(1+(fix(length(radarDataNormalized1(:,1))/Segment_length)-1)/(1-Segment_overlap));
Segment_ofset=1+(0:(Segment_Num-1))*round(Segment_length*(1-Segment_overlap));

%Segment_Num = fix(length(radarDataNormalized1(:,1))/Segment_length);
WT2_S2N = zeros(Segment_Num,num_fast_chan);
%load('NoiseMean.mat')

for j = 1: Segment_Num
   j 
   clear accept_range1 accept_range2 m signalInPersonBin1 signalInPersonBin2 breathingRateVector2 breathingRateVector1 breathing_rate1 breathing_rate2
if BREATHING_BELTS   
   br_index=find(A(:,1)>=radarDataTimes1(Segment_ofset(j)) & A(:,1)<=radarDataTimes1(Segment_ofset(j)+Segment_length-1));
   BR_belt11 = resample(A(br_index,2),length(radarDataNormalized1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,:)),length(br_index));
   if A(br_index(end),1)<radarDataTimes1(Segment_ofset(j)+Segment_length-1)-2
       break
   end
   [breathingRateBelt] = analyze_UWB_breathing( BR_belt11);
   if std(breathingRateBelt)>2
       [b11, c2]=wt_decompose(BR_belt11);
       [breathingRateBelt] = analyze_UWB_breathing( b11');
   end
   TF = find(abs(breathingRateBelt-mean(medfilt1(breathingRateBelt,1)))<1.2*std(medfilt1(breathingRateBelt,1))); %isoutlier(breathingRateBelt, 'quartiles');
   breathingRateRef=mean(breathingRateBelt((TF)));
   if TWO_PERSON==true
        BR_belt22 = resample(A1(br_index,2),length(radarDataNormalized1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,:)),length(br_index));
       [breathingRateBelt2] = analyze_UWB_breathing( BR_belt22);
        if std(breathingRateBelt2)>2
            [b22, c2]=wt_decompose(BR_belt22);
            [breathingRateBelt2] = analyze_UWB_breathing( b22');
        end
        TF2 = find(abs(breathingRateBelt2-mean(medfilt1(breathingRateBelt2,1)))<1.2*std(medfilt1(breathingRateBelt2,1))); %isoutlier(breathingRateBelt, 'quartiles');
        breathingRateRef2=mean(breathingRateBelt2((TF2)));       
   end
end  
   [SVDFiltered1] = analyze_UWB_denoise(radarDataNormalized1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,:));
if TWO_RADARS==true
   [SVDFiltered2] = analyze_UWB_denoise(radarDataNormalized2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,:));
end
  %SVDFiltered1=radarDataNormalized1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,:);
  %SVDFiltered2=radarDataNormalized2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,:);
  for i=1:178
      c=corrcoef(unwrap(angle(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,i)))',(abs(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,i)))');
      r1(i)=abs(c(1,2));
     % c=corrcoef(unwrap(angle(radarRawData2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,i)))',(abs(radarRawData2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,i)))');
     % r2(i)=abs(c(1,2));
  end      
[ranges1, num_of_people1] = analyze_breathing_2people(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,:)');
    if TWO_RADARS==true
    [ranges2, num_of_people2] = analyze_breathing_2people(radarRawData2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,:)');
end
  % [ranges1, num_of_people1] = analyze_UWB_range_2people(unwrap(angle(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,:)))');%(SVDFiltered1);
   %[ranges2, num_of_people2] = analyze_UWB_range_2people(unwrap(angle(radarRawData2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,:)))');%(SVDFiltered2);%
   for(k=1:num_of_people1)
%        if (bandpower(imag(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges1(k)))-mean(imag(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges1(k)))),fs,[0.1 0.6]) > bandpower(real(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges1(k)+1))-mean(real(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges1(k)+1))),fs,[0.1 0.6]))
%            signalInPersonBin1(:,k)= smooth(imag(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges1(k))),3);
%        else
%            signalInPersonBin1(:,k)= smooth(real(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges1(k))),3);
%        end
        signalInPersonBin1(:,k)=unwrap(angle(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges1(k))));
       [fitresult, gof] = createFit1(1:length(signalInPersonBin1(:,k)), signalInPersonBin1(:,k)');
       [b1(k,:), c1(k,:)]=wt_decompose(signalInPersonBin1(:,k)-fitresult(1:length(signalInPersonBin1(:,k))));
       [breathingRateVector1(k,:)] = analyze_UWB_breathing( b1(k,:)' );%analyze_UWB_breathing(signalInPersonBin1(:,k)-mean(signalInPersonBin1(:,k)));%
       TF1 = find(abs(breathingRateVector1(k,:)-mean(medfilt1(breathingRateVector1(k,:),1)))<1.2*std(medfilt1(breathingRateVector1(k,:),1)))
       var1=var(breathingRateVector1(k,TF1));
       if (c1(k,2)>40 & var1< 6 & c1(k,1)*distances(ranges1(k))^0.8>0.4*10^-9 & max(r1(ranges1(k)-1:ranges1(k)+1))>0.5) 
           accept_range1(k)=1
       else 
           accept_range1(k)=0;
       end
       breathing_rate1{j}(k)=mean(breathingRateVector1(k,TF1));
   end
   if num_of_people1>0
    ranges1(find(accept_range1==0))=[];
    breathing_rate1{j}(find(accept_range1==0))=[];
    num_of_people1=size(ranges1,2);
   end
   
if TWO_RADARS==true
    for(k=1:num_of_people2)
%        if (bandpower(imag(radarRawData2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges2(k)))-mean(imag(radarRawData2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges2(k)))),fs,[0.1 0.6]) > bandpower(real(radarRawData2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges2(k)+1))-mean(real(radarRawData2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges2(k)+1))),fs,[0.1 0.6]))
%            signalInPersonBin2(:,k)= smooth(imag(radarRawData2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges2(k))),3);
%        else
%            signalInPersonBin2(:,k)= smooth(real(radarRawData2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges2(k))),3);
%        end
signalInPersonBin2(:,k)=unwrap(angle(radarRawData2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges2(k))));
       [fitresult, gof] = createFit1(1:length(signalInPersonBin2(:,k)), signalInPersonBin2(:,k));
       [b2(k,:), c2(k,:)]=wt_decompose(signalInPersonBin2(:,k)-fitresult(1:length(signalInPersonBin2(:,k))));
        [breathingRateVector2(k,:)] = analyze_UWB_breathing( b2(k,:)'  ); %analyze_UWB_breathing(signalInPersonBin2(:,k)-mean(signalInPersonBin2(:,k)));
        TF2 = find(abs(breathingRateVector2(k,:)-mean(medfilt1(breathingRateVector2(k,:),1)))<1.2*std(medfilt1(breathingRateVector2(k,:),1)))
        var2=var(breathingRateVector2(k,TF2));
        if (c2(k,2)>40 & var2< 6 & c2(k,1)*distances(ranges2(k))^0.8>0.4*10^-9& max(r2(ranges2(k)-1:ranges2(k)+1))>0.5) 
           accept_range2(k)=1
        else 
           accept_range2(k)=0;
        end       
        breathing_rate2{j}(k)=mean(breathingRateVector2(k,TF2));
   end
   if num_of_people2>0
    ranges2(find(accept_range2==0))=[];
    num_of_people2=size(ranges2,2); 
    breathing_rate2{j}(find(accept_range2==0))=[];
   end
%    dist_x(j)=(distances(ranges2(j))^2-distances(ranges1(j))^2+D^2)/(2*D);
%    dist_y(j)=sqrt(distances(ranges2(j))^2-dist_x(j)^2);
   
   
   %TF1 = isoutlier(breathingRateVector1, 'quartiles');
   %var1=var(breathingRateVector1(not(TF1)));
        

   
%    TF2 = isoutlier(breathingRateVector2, 'quartiles');
%    var2=var(breathingRateVector2(not(TF2)));

%    if bandpower(signalInPersonBin1,fs,[0.05 0.6])> bandpower(signalInPersonBin2,fs,[0.05 0.6])
%        breathing_rate(j)=mean(breathingRateVector1((TF1)));
%    else
%        breathing_rate(j)=mean(breathingRateVector2((TF2)));
%    end
%    error(j)=breathingRateRef-breathing_rate(j);
   %presence(radarDataNormalized1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,:), radarDataNormalized2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,:));

    num_detections=0;
       for k1=1:num_of_people1
            for k2=1:num_of_people2
                %[psi,xval] = mexihat(-3,3,60);
                %r1 = xcor_window(psi,unwrap(angle(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges1(k1)))),0,1);
                %r2 = xcor_window(psi,unwrap(angle(radarRawData2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges2(k2)))),0,1);
                %c=corrcoef(r2,r1);
                c_ph=corrcoef(smooth(unwrap(angle(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges1(k1))))),smooth(unwrap(angle(radarRawData2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges2(k2))))));
                m_ph(k1,k2)=abs(c_ph(1,2));
                c_abs=corrcoef(smooth((abs(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges1(k1))))),smooth((abs(radarRawData2(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges2(k2))))));
                m_abs(k1,k2)=abs(c_abs(1,2));
                if (m_ph(k1,k2) > 0.45 | m_abs(k1,k2) > 0.45) & abs(breathing_rate2{j}(k2)-breathing_rate1{j}(k1))<2
                    %tracks{j}=[];
                    if num_detections == 0
                        detections{j}=[ranges1(k1) ranges2(k2)];
                        num_detections=num_detections+1;
                        dist_x=(distances(ranges2(k2))^2-distances(ranges1(k1))^2+D^2)/(2*D);
                        dist_y=sqrt(distances(ranges2(k2))^2-dist_x^2);
                        figure
                        %(2.1^2-2.5^2+2.8^2)/(2*2.8);sqrt(2.1^2-1.0714^2)
                        % (1.5^2-2.7^2+3.25^2)/(2*3.25); sqrt(1.5^2-0.85^2)
                        plot(0.85, 1.24, 'o')
                        text(0.85+0.1, 1.24-0.1, num2str(breathingRateRef));
                        if TWO_PERSON==true
                            plot(2.0286, 1.6270, 'o')
                            text(2.0286+0.1, 1.6270-0.1, num2str(breathingRateRef2));
                        end
                        xlim([0, 5])
                        ylim([0, 5])
                        hold on
                        plot(0, 0, 'x', D, 0, 'x')
                    else
                        detections{j}=[detections{j}; [ranges1(k1) ranges2(k2)]];
                        num_detections=num_detections+1;
                    end
                    dist_x=(distances(ranges2(k2))^2-distances(ranges1(k1))^2+D^2)/(2*D);
                    dist_y=sqrt(distances(ranges2(k2))^2-dist_x^2);
                    plot(dist_x, dist_y, '*')
                    text(dist_x+0.1, dist_y-0.1, num2str(breathing_rate2{j}(k2)));
                    %tracks{j}=[track(j} [dist_x,dist_y]];
                    
                end
                
            end
       end
       if j>1
           if num_detections > 0
               if size(tracks{j-1}, 1) > 0
                   clear delta costMatrix;
                   for i = size(tracks{j-1}, 1):-1:1
                       delta = detections{j} - tracks{j-1}(i, :);
                        costMatrix(i, :) = sqrt(sum(delta .^ 2, 2));
                   end
                   costofnonassignment = 10;
                   [assignments, unassignedTracks, unassignedDetections] = assignmunkres(costMatrix,costofnonassignment);
                   tracks{j}=detections{j}(assignments(:,2), :);
                   unassigned_Detections{j}=detections{j}(unassignedDetections,:);
                   if length(tracks{j})==0
                       if length(unassigned_Detections{j})>0
                           tracks{j}=unassigned_Detections{j};
                       else
                           tracks{j}=tracks{j-1};
                       end
                   end
               else
                   tracks{j}=detections{j};
               end
           else
               tracks{j}=tracks{j-1};
           end
       else
           if num_detections > 0
                tracks{j}=detections{j};
           else
               tracks{j}=[];
           end
       end
       if TWO_PERSON==false 
           if j>1 & num_detections>0
                error(j)=breathingRateRef-breathing_rate1{j}(assignments(:,1));
           end
       end
end
       if TWO_PERSON==false 
                error(j)=breathingRateRef-breathing_rate1{j};
       end
end
%plot(0:7.5:7.5*length(error)-1,error,'o-')
%xlabel('Time (sec)');ylabel('Error in breaths per minute');
%% 
if TWO_RADARS==true
time_i=0:7.5:7.5*length(error)-1;
for i=1:length(error)
    figure(5)
    dist_x=(distances(tracks{i}(1))^2-distances(tracks{i}(2))^2+D^2)/(2*D);
    dist_y=sqrt(distances(tracks{i}(2))^2-dist_x^2);
    plot(time_i(i), sqrt(dist_x^2+ dist_y^2), '*-')
    text(time_i(i)+0.1, sqrt(dist_x^2+ dist_y^2)-(-1)^i*0.01, num2str(error(i)));
    if length(unassigned_Detections{i})>0
        dist_x=(distances(unassigned_Detections{i}(1))^2-distances(unassigned_Detections{i}(2))^2+D^2)/(2*D);
        dist_y=sqrt(distances(unassigned_Detections{i}(2))^2-dist_x^2);
        plot(time_i(i), sqrt(dist_x^2+ dist_y^2), 'o-')
    end
    pause(0.5)
    hold on
    F(i)=getframe;
end

plot(0:1/fs:15-1/fs,BR_belt11/10)
hold on
plot(0:1/fs:15-1/fs,unwrap(angle(radarRawData1(Segment_ofset(j):Segment_ofset(j)+Segment_length-1,ranges1(k)))))
xlabel('Time (sec)');ylabel('Amplitude of Breathing belt (blue) and radar (red)');
end
%movie(F,1,1)
%  v = VideoWriter('newfile.avi','Uncompressed AVI');
% v.FrameRate=1;
% open(v)
% writeVideo(v,F)
% close(v)

% m = ssest(signalInPersonBin1,40);
% for i=1:120
% [y,fit,x0] = compare(SVDFiltered1(i,:)',m);
% fit1(i)=fit;
% end
%% Classification of activity


