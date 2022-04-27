function [vector2]=LPMIHN(seed1, seed2, Net1,Net2, transNet,alpha)
% seed1,2为列向量
% Net1,2做列归一
% transNet维度=N2*N1，行归一
seed1=colnorm(seed1);
seed2=colnorm(seed2);
Net1=Net1-diag(diag(Net1));
Net2=Net2-diag(diag(Net2));
Net1=colnorm(Net1);
Net2=colnorm(Net2);
transNet=(colnorm(transNet'))';

vector1=RW(seed1,Net1,alpha);
vector1=colnorm(vector1);

seed2=(1-2*alpha)/alpha*seed2+alpha/(1-alpha)*transNet*vector1;
vector2=RW(seed2,Net2,alpha);
%vector2=colnorm(vector2);
end