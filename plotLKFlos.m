function plotLKFlos(input_h5_file,inputlkf_h5_file,Tri,Xkp,lambda,k)
% plotLKFlos: showing the LOS results from linear Kalman filter 
%   example
%   >> plotLKFlos('modelinput.h5','input_for_lkf_Su.h5',TChaman2Su,Xkp2,100(lambda),22);
%
%   Masato Furuya (c) 2021-
imdates=h5read('TS_GEOCml2clip.UTM/cum_filt.h5','/imdates');
%N=length(imdates);Yr=zeros(N,1);
Yr=datenum(num2str(imdates(k)),'yyyymmdd');
YrEvt1=datenum('20160513','yyyymmdd');
YrEvt2=datenum('20160710','yyyymmdd');
% Time-difference vector in "days"
dt1=Yr-YrEvt1;
dt2=Yr-YrEvt2;
% Number of triangular meshes
m=size(Tri,1);
Ref=h5read(input_h5_file,'/Ref');
Rnf=h5read(input_h5_file,'/Rnf');
%
[los, Base, Blosam, qt, ~, ~, ~] = DataErrorCov(input_h5_file,inputlkf_h5_file,Tri,lambda,k);
figure('position', [300, 500, 134, 400]);colormap(gca,'Parula');                      %  ;
%figure('position', [300, 500, 400, 134]);colormap(gca,'Parula');                      %  ;
subaxis(3,1,1,'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);surf(Ref,Rnf,los);view(0,90);shading flat;axis tight;caxis([-2.5 2.5]);
if k<25
   % title(strcat('Obs (',num2str(dt1),'d)'))
   text(232200,3420540,2.1,'Obs','FontSize',14)
   text(250000,3378000,2.1,strcat(num2str(dt1),'d'),'FontSize',14)
else
   % title(strcat('Obs (',num2str(dt1),'d/',num2str(dt2),'d)'));
   text(232200,3420540,2.1,'Obs','FontSize',14)
   text(250000,3378000,2.1,strcat(num2str(dt1),'/',num2str(dt2),'d'),'FontSize',14)
end
%subplot(131);xlabel('UTM East (m)');ylabel('UTM North (m)') % colorbar;  title(strcat(num2str(imdates(k)),' (Obs, cm)'));
%                                          lsq                                             log10(Xkp(1:1/2*length(Xkp)))
%set(gca,'xtick',[])
set(gca,'xticklabel',[])
%set(gca,'ytick',[])
set(gca,'yticklabel',[])
ical2=subsample_back(qt,Blosam,Base*Xkp(k).smt(1:2*m));
subaxis(3,1,2,'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);surf(Ref,Rnf,ical2);view(0,90);shading flat;axis tight;caxis([-2.5 2.5]);
if k<25
   %title(strcat('Cal (',num2str(dt1),'d)'))
   text(232200,3420540,2.1,'Cal','FontSize',14)
   text(250000,3378000,2.1,strcat(num2str(dt1),'d'),'FontSize',14)
else
   % title(strcat('Cal (',num2str(dt1),'d/',num2str(dt2),'d)'));
   text(232200,3420540,2.1,'Cal','FontSize',14)
   text(250000,3378000,2.1,strcat(num2str(dt1),'/',num2str(dt2),'d'),'FontSize',14)
end
%subplot(132);xlabel('UTM East (m)');%ylabel('UTM North (m)') colorbar;
%
%set(gca,'xtick',[])
set(gca,'xticklabel',[])
%set(gca,'ytick',[])
set(gca,'yticklabel',[])
subaxis(3,1,3,'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);surf(Ref,Rnf,los-ical2);view(0,90);shading flat;axis tight;caxis([-1.5 1.5]);
%title('Misfit')
text(232200,3420540,2.1,'Misfit','FontSize',14)
if k<25
   %title(strcat('Cal (',num2str(dt1),'d)'))
   text(250000,3378000,2.1,strcat(num2str(dt1),'d'),'FontSize',14)
else
   % title(strcat('Cal (',num2str(dt1),'d/',num2str(dt2),'d)'));
   text(250000,3378000,2.1,strcat(num2str(dt1),'/',num2str(dt2),'d'),'FontSize',14)
end
%subplot(133);xlabel('UTM East (m)');%ylabel('UTM North (m)')   colorbar;  title(strcat(num2str(imdates(k)),' (Misfit, cm)'))
%set(gca,'xtick',[])
set(gca,'xticklabel',[])
%set(gca,'ytick',[])
set(gca,'yticklabel',[])
saveas(gcf,strcat('~/Desktop/Chaman-S1/Chaman2Su/',num2str(imdates(k)),'loscal.png'));


