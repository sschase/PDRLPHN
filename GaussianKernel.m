%GaussianKernel
% 列向量间的SWscore,用作相似性
function [sim]=GaussianKernel(Net)

[~,n]=size(Net);

gamma = n/sum(diag(Net'*Net));

sim=zeros(n,n);
for i=1:n
    for j=1:n
        vi=Net(:,i);
        vj=Net(:,j);
        d=vi-vj;
        e=d'*d;
        if e~=0
            sim(j,i)=exp( -gamma* e ); %列归一 
        else
            sim(j,i)=0;
        end
    end
end
sim=sim-diag(diag(sim));
end