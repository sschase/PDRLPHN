clear 
clc
%==============parameters==================%
F = 1;
alpha=0.2;
scale_d=0.10;
scale_c=0.10;
Savepath   = ['./result/cv/' num2str(F) '/Sd2c_1.mat'];
%==========================================%

%====================================input data==================================%
%没有做过归一，只做过角0处理的相似性矩阵
load('./input/sim.mat');
simc_g   =simc_g_corr;      clear simc_g_corr;
simc_exp =simc_exp_corr;    clear simc_exp_corr;
simd_g   =simd_g_corr;      clear simd_g_corr;
simd_chem=simd_chem_corr;    clear simd_chem_corr;
load('./input/sensi'); %1表示存在，0表示存在以外
fprintf('data loaded!\n');
%================================================================================%


%================================transNet=sensi============================% 
[m,n]=size(sensi);
Score_d2c_sen=zeros(m,n);
Miss_d2c_sen=zeros(m,n);

%--------------missing value-----------------
% drug-->cline 
for j=1:94
    %种子d 
    seed_d=zeros(n,1);
    seed_d(j)=1;
    % missing
    A=sensi;
    miss=find(A(:,j)==0);   
    seed_c=A(:,j);
    fprintf('miss D-->C j=%d\n',j);
    %--------------sim matrix processing---------------%
    simc_cdp=GaussianKernel(A');
    simd_cdp=GaussianKernel(A );
    simc=normdiag(simc_cdp)+normdiag(simc_exp)+normdiag(simc_g);simc=normdiag(simc);
    simd=normdiag(simd_cdp)+normdiag(simd_chem)+normdiag(simd_g);simd=normdiag(simd);
    simc=sparse_sim(simc,scale_c);simc=colnorm(simc);
    simd=sparse_sim(simd,scale_d);simd=colnorm(simd);
    %收敛
    %--------------LP--------------%
    Score=LPMIHN(seed_d,seed_c,simd,simc,A,alpha);
    Miss_d2c_sen(miss,j)=Score(miss);
end
% save (Savepath_M, 'Miss_d2c_sen'); 
% --------------------------LOOCV -----------------------%
% Drug-Cline
for j=1:94  %对一个特定药物dj
    %种子d 
    dj=j;
    seed_d=zeros(n,1);
    seed_d(dj)=1;
    
    A=sensi;%A是cline*drug(行*列)的元素为{0,1}的转移矩阵
    index_i=find(A(:,dj)~=0);   
    
    for i=1:length(index_i)
        ci=index_i(i);
        A=sensi;
        A(ci,dj)=0 ;
        seed_c=A(:,dj);
        %构造转移矩阵
        fprintf('sensi D-->C   c=%d    d=%d\n',ci,dj);           
        %--------------sim matrix processing---------------%
        simc_cdp=GaussianKernel(A');
        simd_cdp=GaussianKernel(A );
        simc=normdiag(simc_cdp)+normdiag(simc_exp)+normdiag(simc_g);simc=normdiag(simc);
        simd=normdiag(simd_cdp)+normdiag(simd_chem)+normdiag(simd_g);simd=normdiag(simd);
        simc=sparse_sim(simc,scale_c);simc=colnorm(simc);
        simd=sparse_sim(simd,scale_d);simd=colnorm(simd);
        %收敛
        %--------------LP--------------%
        Score=LPMIHN(seed_d,seed_c,simd,simc,A,alpha);%LP中自带归一化
        Score_d2c_sen(ci,dj)=Score(ci);        
    end
end
% save (Savepath_S, 'Score_d2c_sen');  
Sd2c_1=Miss_d2c_sen+Score_d2c_sen;
save (Savepath,'Sd2c_1');