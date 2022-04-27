function[p]=RW(p0,M,alpha)
backprobability=1-alpha;
p1=p0;
num=1;
p=(1-backprobability)*M*p1+backprobability*p0;
while max(p-p1)>1e-10
    num=num+1;
    p1=p;
    p=(1-backprobability)*M*p1+backprobability*p0;
end
end