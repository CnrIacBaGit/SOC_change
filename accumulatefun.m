function [b, accTSMD]=accumulatefun(temp,rain,evap,s_cover,clay,d)

n=length(temp);

% the rate modifying factor for temperature
a = 47.91./ ( 1+ exp( 106.06 ./ (temp+18.27) ) );

% the rate modifying factor for moisture
maxTSMD = -(20 + 1.3*clay-0.01*clay^2)*d/23;
vec=rain-evap;
accTSMD=zeros(n,1);
accTSMDf=zeros(n,1);
i=1;
while vec(i)>0
    i=i+1;
end
if i>1
for j=i:n
    accTSMD(j)= min( max( accTSMD(j-1)+vec(j) , maxTSMD ), 0);
end
end

b=ones(n,1);
ind= accTSMD<=0.444*maxTSMD;
b(ind) = 0.2+ (1-0.2)* (maxTSMD-accTSMD(ind))./(maxTSMD-0.444*maxTSMD);

c=ones(n,1);
c(s_cover==1)=0.6;


end