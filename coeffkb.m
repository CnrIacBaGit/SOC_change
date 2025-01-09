function [kb, Acc, M, Mb] = coeffkb(cly,d,rain,pet,Nm,Ny,T) 

M = -(20+1.3*cly-.01*cly^2)*d/23;
Mb = .444*M;
vec=rain-pet;
Acc=zeros(Nm,1);
%accTSMDf=zeros(n,1);
for yr=0:Ny-1
    i=yr*T+1;
    while vec(i)>0
        i=i+1;
    end
    for j=i:(yr+1)*T
        Acc(j)= min( max( Acc(j-1)+vec(j) , M ), 0);
    end
end
kb = ones(Nm,1);
ind = Acc<=Mb;
kb(ind) = 0.2 + 0.8 * (M-Acc(ind))/(M-Mb);