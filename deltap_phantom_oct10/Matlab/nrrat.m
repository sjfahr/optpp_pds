function [avgval stderr avgerr] = nrrat(number1,number2)
%input in order of increasing temperature

[timex,timey1,sigmat1]=nrtime(number1);
[timex,timey2,sigmat2]=nrtime(number2);
clf

rat=wiener2(timey2,[3,1])./wiener2(timey1,[3,1]);

min1=find(timey1>5,1,'first');
%min1=7
max1=find(timex==210);

%sigmat1=.5*ones(numel(timex),1)
%sigmat2=sigmat1


raterror=rat.*sqrt((sigmat1./timey1).^2+(sigmat2./timey2).^2);

errorbar(timex(min1:max1),rat(min1:max1),raterror(min1:max1),'or')
xlabel('Time(s)')
ylabel('\Delta T_R/\Delta T_S (^oC)')

avgval=mean(rat(min1:max1));
stderr=std(rat(min1:max1));
avgerr=sqrt(sum(raterror(min1:max1).^2))./numel(raterror(min1:max1));
