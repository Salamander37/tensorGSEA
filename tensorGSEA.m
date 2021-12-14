% %%%%%%% collect data for analysis
% %%%%%%% 2015.5.14
% clear all;
% clc;
% 
% %load  Total_Rate_DtoH_2_2.mat; %temp use
% 
% %% read data
% GG=importdata( 'GSE13270.geneSymbol.ND.txt');  %value data 'GG.textdata' for name of Gene
% M=importdata('new11.xlsx');  %Gene functions which contain a set of Gene 'M.Sheet1'
% 
% total_data=GG.data;  %double
% data_size=size(GG.data);
% gene_type=data_size(1);
% sample_number=data_size(2);
% 
% % data normalization
% for i=1:gene_type
%     ss=total_data(i,:);
%     norm_ss=mapminmax(ss,0,1);
%     normed_total_data(i,:)=norm_ss;
% end
% 
% Search_size=size(M.Sheet1);
% 
% 
% 
% tic
% 
% %exact the corresponding data from sheet data
% for fun_Num=1:Search_size(1)
%     clear NumS;
%     
%     for j=3:Search_size(2)
%         G_Name=M.Sheet1{fun_Num,j};
%         count=0;
%         count=count+~isempty(G_Name);
%         if count~=0
%             ss=strmatch(G_Name,GG.textdata,'exact');
%             NumS(j-2)= ss;        % the corresponding numbers
%         else     %NumS(i,j-2)=0;
%         end
%     end
% 
% 
% %% data writing out
% Total_size=size(NumS);
% %GeneFun_Num=Total_size(1);
% Gene_Num=Total_size(2);
% total_tensor=[];
% 
% %gene function part
% %for j=1:GeneFun_Num
% for ge=1:Gene_Num
%         total_tensor(ge,:)=normed_total_data(NumS(ge),:);
% end
% 
% total_tensor_Size=size(total_tensor);
% 
% total_tensor2=reshape(total_tensor,[total_tensor_Size(1),10,5]);
% 
% Training_tensor=total_tensor2(:,1:5,:);
% Testing_tensor=total_tensor2(:,6:10,:);
% 
% dim=size(Training_tensor);
% Training_tensor=double(Training_tensor);
% Training_tensor2=tensor(Training_tensor,[dim(1)  dim(2) dim(3)]);
% 
% G=tucker_als(Training_tensor2,[dim(1) 1 1]);
% 
% Testing_tensor2=tensor(Testing_tensor,[dim(1)  dim(2) dim(3)]);
% Recons_Test_cor_tensor=ttm(Testing_tensor2,{G.u{1}',G.u{2}',G.u{3}'},[1,2,3]);
% 
% Recons_Test_tensor=ttm(Recons_Test_cor_tensor,{G.u{1},G.u{2},G.u{3}},[1,2,3]);
% 
% Testing_tensor2=double(Testing_tensor2);
% Recons_Test_tensor=double(Recons_Test_tensor);
% 
% Ori_Vec=reshape(Testing_tensor2,dim(1)*dim(2)*dim(3),1);
% Rec_Vec=reshape(Recons_Test_tensor,dim(1)*dim(2)*dim(3),1);
% 
% 
% a=Ori_Vec'*Rec_Vec;
% b=Ori_Vec'*Ori_Vec;
% c=Rec_Vec'*Rec_Vec;
% 
% Genefun_NorCorRate=a/sqrt(b*c);
% 
% % other parts
% other_tensor=normed_total_data;
% other_tensor(find(other_tensor(NumS)),:)=[];
% other_gene_number=gene_type-Gene_Num;
% 
% 
% for iteration_num=1:1000 
%     aaa=[1:other_gene_number];
%     RandomNum=randperm(numel(aaa));
%     randomtensor=other_tensor(RandomNum(1:Gene_Num),:);
%     
%     total_tensor_Size=size(total_tensor);
%     
%     total_tensor2=reshape(randomtensor,[total_tensor_Size(1),10,5]);
%     
%     Training_tensor=total_tensor2(:,1:5,:);
%     Testing_tensor=total_tensor2(:,6:10,:);
%     
%     dim=size(Training_tensor);
%     Training_tensor=double(Training_tensor);
%     Training_tensor2=tensor(Training_tensor,[dim(1)  dim(2) dim(3)]);
%     
%     G=tucker_als(Training_tensor2,[dim(1) 1 1]);
%     
%     Testing_tensor2=tensor(Testing_tensor,[dim(1)  dim(2) dim(3)]);
%     Recons_Test_cor_tensor=ttm(Testing_tensor2,{G.u{1}',G.u{2}',G.u{3}'},[1,2,3]);
%     
%     Recons_Test_tensor=ttm(Recons_Test_cor_tensor,{G.u{1},G.u{2},G.u{3}},[1,2,3]);
%     
%     Testing_tensor2=double(Testing_tensor2);
%     Recons_Test_tensor=double(Recons_Test_tensor);
%     
%     Ori_Vec=reshape(Testing_tensor2,dim(1)*dim(2)*dim(3),1);
%     Rec_Vec=reshape(Recons_Test_tensor,dim(1)*dim(2)*dim(3),1);
%     
%     
%     a=Ori_Vec'*Rec_Vec;
%     b=Ori_Vec'*Ori_Vec;
%     c=Rec_Vec'*Rec_Vec;
%     
%     Other_NorCorRate(iteration_num)=a/sqrt(b*c);
% end
% 
% Total_sub_Rate=[Genefun_NorCorRate Other_NorCorRate];
% 
% Total_Rate(:,fun_Num)=Total_sub_Rate;
% 
% 
% end
% 
% save('Total_Rate_DtoH_1_1.mat','Total_Rate');
% 
% toc
    
%%%%%%% 2015.5.14
%%
clear all;
clc;

%% read data
% GG=importdata('GSE13268.geneSymbol.ND.txt');
% M=importdata('new11.xlsx');
GG=importdata( 'GSE13270.geneSymbol.ND.txt');  %value data 'GG.textdata' for name of Gene
M=importdata('new11.xlsx');  %Gene functions which contain a set of Gene 'M.Sheet1'

% load Total_Rate_HtoD_1_1.mat;%%Total_Rate

total_data=GG.data;  %double
data_size=size(GG.data);%25345*50
gene_type=data_size(1);
sample_number=data_size(2);

% data normalization
for i=1:gene_type
    ss=total_data(i,:);
    norm_ss=mapminmax(ss,0,1);
    normed_total_data(i,:)=norm_ss;
end

Search_size=size(M.Sheet1);

tic
%exact the corresponding data from sheet data
for fun_Num=198:200 %Search_size(1), the number of pathways
% for fun_Num = 154
    %fun_Num
    for j=3:Search_size(2)
        G_Name=M.Sheet1{fun_Num,j};
        count=0;
        count=count+~isempty(G_Name);
        if count~=0
            ss=strmatch(G_Name,GG.textdata,'exact');
            %textdata中基因的编号和data中似乎是错开一位的
%             NumS(j-2)= ss;        % the corresponding numbers.NumS中为对应的基因序号
            NumS(j-2) = ss-1;
        else     %NumS(i,j-2)=0;
        end
    end
    
    %% data writing out
    Total_size=size(NumS);
    %GeneFun_Num=Total_size(1);
    Gene_Num=Total_size(2);
    total_tensor=[];
    
    %gene function part
    %for j=1:GeneFun_Num
    %for one function
    for ge=1:Gene_Num
        total_tensor(ge,:)=normed_total_data(NumS(ge),:);
    end
    
    total_tensor_Size=size(total_tensor);%genenumber*sample(50)
    total_tensor2=reshape(total_tensor,[total_tensor_Size(1),10,5]);%gene,sample,period
    Testing_tensor=total_tensor2(:,1:5,:);%GK, disease state
    Training_tensor=total_tensor2(:,6:10,:);%WKY, normal state
    
    dim=size(Training_tensor);%gene*5*5
    Training_tensor=double(Training_tensor);
    Training_tensor2=tensor(Training_tensor,[dim(1)  dim(2) dim(3)]);
    
    G=tucker_als(Training_tensor2,[dim(1) 4 4]);%%Compress
    %computes the best rank(R1,R2,..,Rn) approximation of tensor X
    
    Testing_tensor2=tensor(Testing_tensor,[dim(1)  dim(2) dim(3)]);
    Recons_Test_cor_tensor=ttm(Testing_tensor2,{G.u{1}',G.u{2}',G.u{3}'},[1,2,3]);
    
    Recons_Test_tensor=ttm(Recons_Test_cor_tensor,{G.u{1},G.u{2},G.u{3}},[1,2,3]);%reconstruct
    
    Testing_tensor2=double(Testing_tensor2);
    Recons_Test_tensor=double(Recons_Test_tensor);
    
    Ori_Vec=reshape(Testing_tensor2,dim(1)*dim(2)*dim(3),1);
    Rec_Vec=reshape(Recons_Test_tensor,dim(1)*dim(2)*dim(3),1);
    
    a=Ori_Vec'*Rec_Vec;
    b=Ori_Vec'*Ori_Vec;
    c=Rec_Vec'*Rec_Vec;
    %calculate NC
    Genefun_NorCorRate=a/sqrt(b*c);
    
    % other parts
    other_tensor=normed_total_data;
    other_tensor(find(other_tensor(NumS)),:)=[];%删除原pathway中的基因
    other_gene_number = gene_type-Gene_Num;
    
    for iteration_num=1:1000
        aaa=[1:other_gene_number];
        RandomNum = randperm(numel(aaa));
        randomtensor = other_tensor(RandomNum(1:Gene_Num),:);        
        total_tensor_Size = size(total_tensor);
        
        total_tensor2=reshape(randomtensor,[total_tensor_Size(1),10,5]);        
        Testing_tensor=total_tensor2(:,1:5,:);
        Training_tensor=total_tensor2(:,6:10,:);
        
        dim=size(Training_tensor);
        Training_tensor=double(Training_tensor);
        Training_tensor2=tensor(Training_tensor,[dim(1)  dim(2) dim(3)]);
        
        G=tucker_als(Training_tensor2,[dim(1) 4 4]);        
        Testing_tensor2=tensor(Testing_tensor,[dim(1)  dim(2) dim(3)]);
        Recons_Test_cor_tensor=ttm(Testing_tensor2,{G.u{1}',G.u{2}',G.u{3}'},[1,2,3]);        
        Recons_Test_tensor=ttm(Recons_Test_cor_tensor,{G.u{1},G.u{2},G.u{3}},[1,2,3]);        
        Testing_tensor2=double(Testing_tensor2);
        Recons_Test_tensor=double(Recons_Test_tensor);
        
        Ori_Vec=reshape(Testing_tensor2,dim(1)*dim(2)*dim(3),1);
        Rec_Vec=reshape(Recons_Test_tensor,dim(1)*dim(2)*dim(3),1);
        
        a=Ori_Vec'*Rec_Vec;
        b=Ori_Vec'*Ori_Vec;
        c=Rec_Vec'*Rec_Vec;
        
        Other_NorCorRate(iteration_num)=a/sqrt(b*c);
    end    
    Total_sub_Rate=[Genefun_NorCorRate Other_NorCorRate];
    Total_Rate(:,fun_Num)=Total_sub_Rate;   
end

toc
save('Total_Rate_13268_20190911_444_all.mat','Total_Rate');
%  save('Total_Rate_HtoD_1_120190910_444_all.mat', 'Total_Rate');
%  save('Total_Rate_HtoD_1_1_new.mat', 'Total_Rate');
    
    
 
 
    
    
    
    
    
    
  