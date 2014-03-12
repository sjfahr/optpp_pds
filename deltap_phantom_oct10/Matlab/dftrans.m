function [forward,inverse]=dftrans(smallf)

t=numel(smallf);
v=t;

forward=zeros(t,t);
for jj=1:v
    for ii=1:t
        forward(jj,ii)=(1/t)*smallf(ii)*exp(-2*pi*1i*(ii-1)*(jj-1)/t);
    end
end
forward=sum(forward,2);

inverse=zeros(t,t);
for ii=1:t
    for jj=1:v
        inverse(ii,jj)=smallf(jj)*exp(2*pi*1i*(ii-1)*(jj-1)/t);
    end
end
inverse=sum(inverse,2);