clc;clear all;

load Total_Rate_HtoD_1_1.mat; % Total_Rate_DtoH_3_3.mat;1001*100
min_Data=min(min(Total_Rate));
max_Data=max(max(Total_Rate));
mean_Data=sum(sum(Total_Rate))/(1001*200);

TOTAL_aa=[min_Data max_Data mean_Data];

[x,y]=sort(Total_Rate,1);%��������x�ǰ��������еľ���y��x�и�λ��Ԫ����ԭ�����е��кţ�Ҳ��λ��
tt=y(1,:);%ȡ��һ��Ԫ��
[x1,y1]=sort(tt,2);
p_Value=x1/1001;

TOTAL_bb=[y1;p_Value];

T_size=size(y1);
M=importdata('new11.xlsx');
%M.Sheet1(:,2)
for i=1:T_size(2)
    ss=y1(i);
    Name_pathway(i)=M.Sheet1(ss,2);
end
    



% %combine the two data ensembles
% 
% load Total_Rate_DtoH_3_3.mat; % Total_Rate_DtoH_3_3.mat
% sum_Data=Total_Rate;
% 
% load Total_Rate_HtoD_3_3.mat;
% sum_Data=(Total_Rate+sum_Data)/2;
% 
% min_Data=min(min(sum_Data));
% max_Data=max(max(sum_Data));
% mean_Data=sum(sum(sum_Data))/(1001*200);
% TOTAL_aa=[min_Data max_Data mean_Data];
% 
% 
% [x,y]=sort(sum_Data,1);
% tt=y(1,:);
% [x1,y1]=sort(tt,2);
% p_Value=x1/1001;
% 
% TOTAL_bb=[y1;p_Value];
% 
% T_size=size(y1);
% M=importdata('new11.xlsx');
% %M.Sheet1(:,2)
% for i=1:T_size(2)
%     ss=y1(i);
%     Name_pathway(i)=M.Sheet1(ss,2);
% end
    








