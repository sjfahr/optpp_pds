function [zprofx,zprofy]= nrprof(number)
%nrprof plots the profile of heating along the cylinder axis at the maximum
%temperature. averages over 2 columns on either side of maximum column

pixsp=.3906;
n=2;

nums=num2str(number); %imports the tmaps specified by number.
load('/export/home/warp/cmaclell/Nanorods/20101019/Allruns2.mat',['tmaps' nums],['dat2s' nums]);
tmap=eval(['tmaps' nums]);
dat2=eval(['dat2s' nums]);

[maxpix,tmax,sigmat,sigmaamp,int,center] = nrchar(number);

iii=maxpix(1);
jjj=maxpix(2);
kkk=maxpix(3);

y = tmap(:,jjj,kkk);
%toprow=min(find(zprofy));
toprow=int(1);
x = ([0:255]);
y2=zeros(256,1);
for m=jjj-n:jjj+n;
    y2=y2+tmap(:,m,kkk)/(2*n+1);
end
zprofx=(x-toprow)*pixsp;
zprofy=y2;
plot(zprofx,zprofy);
xlabel('Distance from agar/agar-NP interface (mm)')
ylabel('\Delta T ^oC')
hline=refline(0,5)
set(hline,'Color','r')