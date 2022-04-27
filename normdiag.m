%norm_sim
function [out]= normdiag( sim )
[n,m]=size(sim);
out=zeros(m,m);
sum_row=sum(sim,2);
sum_row(sum_row==0)=1;

for i=1:n
    for j=1:m
out(i,j)=sim(i,j)/sqrt( sum_row(i)*sum_row(j)   );
    end
end
out=out-diag(diag(out));
end

