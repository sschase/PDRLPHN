function [index_n1]=get_n1 (sim,  point, scale )

simvec=sim(:,point);
[~,index0] = sort(simvec,'descend');
index_n1 =index0( 1:ceil(scale*length(simvec)) );     

end

