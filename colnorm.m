function [ normalizedNetwork ] = colnorm( Net )

deg = sum(Net,1);
zeros = (deg == 0);
deg(zeros) = 1;
deg = 1./deg;
normalizedNetwork = bsxfun(@times, Net, deg);
%normalizedNetwork = sparse(normalizedNetwork);
end