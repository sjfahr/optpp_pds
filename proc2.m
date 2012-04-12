
if ~exist('ExactDataPyruvate','var')
    %load('pyruvate.mat','ExactDataPyruvate','x')
    %ExactData=ExactDataPyruvate;
    load('lactate.mat','ExactDataLactate','x')
    ExactData=ExactDataLactate;
end
[mat1 mat2 nparm]=size(ExactData);
dmask3d=zeros(mat1,mat2,nparm);
dmask3d(2:end-1,2:end-1,:)=1;
dmask3d=(dmask3d==1);
dmask2d=squeeze(dmask3d(:,:,1));
ampmsk=dmask3d;ampmsk(:,:,2:end)=0;
s1pmsk=dmask3d;s1pmsk(:,:,3:end)=0;s1pmsk(:,:,1)=0;
s2pmsk=dmask3d;s2pmsk(:,:,1:2)=0;s2pmsk(:,:,4)=0;
delmsk=dmask3d;delmsk(:,:,1:3)=0;

FitData=zeros(size(ExactData));
%FitData(2:7,2:7,:)=reshape(x,6,6,4);
FitData(dmask3d)=x; %=reshape(x,mat1-2,mat2-2,nparm);
tt=[10:90];
%tmpk=zeros(256,256);
LoResFitD=zeros(mat1,mat2);
LoResTru=zeros(mat1,mat2);
HiResSiz=128;
tmpk=zeros(HiResSiz);
projsiz=length(radon(tmpk,0));
maxlo=0;
maxhi=0;
%set figures:
fg1=figure(1);set(gcf,'Name','Lo-Res Fit');colormap(hot)
fg2=figure(2);set(gcf,'Name','Hi-Res Fit');colormap(hot)
fg3=figure(3);set(gcf,'Name','Lo-Res Truth');colormap(hot)
fg4=figure(4);set(gcf,'Name','Hi-Res Truth');colormap(hot)
fg5=figure(5);set(gcf,'Name','Projections');colormap(hot)

Projs=zeros(projsiz,tt(end)-tt(1)+1);
NoLoc=zeros(1,length(tt));

LoResheader.SliceThickness = .5;
LoResheader.PixelSpacing = [1.;1.] ;
HiResheader.SliceThickness = .5;
HiResheader.PixelSpacing = [.125;.125] ;
projheader.SliceThickness = .5;
projheader.PixelSpacing = [1.;3.] ;

fps=3;
for tc=1:length(tt)
    t=tt(tc);
    for ti=0:fps-1;
        
        LoResFitD(dmask2d)=FitData(ampmsk).*gampdf(t+ti/fps-FitData(delmsk),FitData(s1pmsk),FitData(s2pmsk));
        tmp=fftshift(fft2(fftshift(LoResFitD)));
        tmpk(HiResSiz/2-mat1/2+1:HiResSiz/2+mat1/2,HiResSiz/2-mat2/2+1:HiResSiz/2+mat2/2)=tmp;
        HiResFitD=fftshift(ifft2(fftshift(tmpk)))*HiResSiz*HiResSiz/mat1/mat2;
        figure(fg1)
        imagesc([-1 1],[-1 1],abs(LoResFitD),[0 0.7])
        writeVTK(LoResFitD,'loresfit',tc, 'image',LoResheader);
        maxlo=max([maxlo max(max(abs(LoResFitD)))]);
        figure(fg2)
        imagesc([-1 1],[-1 1],abs(HiResFitD),[0 0.8])
        writeVTK(HiResFitD,'hiresfit',tc, 'image',HiResheader);
        maxhi=max([maxhi max(max(abs(HiResFitD)))]);
        title(sprintf('T = %d',t))    
        LoResTru(dmask2d)=ExactData(ampmsk).*gampdf(t+ti/fps-ExactData(delmsk),ExactData(s1pmsk),ExactData(s2pmsk));
        tmp=fftshift(fft2(fftshift(LoResTru)));
        tmpk(HiResSiz/2-mat1/2+1:HiResSiz/2+mat1/2,HiResSiz/2-mat2/2+1:HiResSiz/2+mat2/2)=tmp;
        HiResTru=fftshift(ifft2(fftshift(tmpk)))*HiResSiz*HiResSiz/mat1/mat2;
        figure(fg3)
        imagesc([-1 1],[-1 1],abs(LoResTru),[0 0.7])
        writeVTK(LoResTru,'lorestru',tc, 'image',LoResheader);
        figure(fg4)
        imagesc([-1 1],[-1 1],abs(HiResTru),[0 0.8])
        writeVTK(HiResTru,'hirestru',tc, 'image',HiResheader);
        title(sprintf('T = %d',t))
        if ti==0
            prangle=t*111.23*pi/180;
            lh=line([-cos(prangle) cos(prangle)],[-sin(prangle) sin(prangle)],'color',[1 1 1]);
            Projs(:,tc)=radon(abs(HiResTru),t*111.23);
            figure(fg5);
            imagesc(Projs)
            writeVTK(Projs,'projection',tc,'projection',projheader);
            NoLoc1(tc)=sum(sum(abs(HiResTru)));
            NoLoc2(tc)=sum(Projs(:,tc));
            NoLoc3(tc)=sum(sum(LoResTru))*HiResSiz*HiResSiz/mat1/mat2;
        end
    
        %pause
    end
end

figure(6)
for ii=2:7
    for jj=2:7
        plot(tt,ExactData(ii,jj,1).*gampdf(tt-ExactData(ii,jj,4),ExactData(ii,jj,2),ExactData(ii,jj,3)))
        hold on
    end
end
hold off
figure(7)
plot(tt,NoLoc1,'-.',tt,NoLoc2,'--',tt,NoLoc3,':');
legend('sum(sum(data))','sum(radon(data))','sum(sum(loresdata))',0)


