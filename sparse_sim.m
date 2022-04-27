% sparse_sim for every node
function [out]=sparse_sim(sim, scale)
if scale==1
    out=sim;
    return;
end
out=zeros(size(sim));
for i=1:length(sim)
    [index_n1]= get_n1(sim,i, scale);
    out(index_n1,i)=sim(index_n1,i);
    out(i,index_n1)=sim(i,index_n1);
end

end