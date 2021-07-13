function plotLKFslip(Tri,Xkp,k)
% plotLKFlos: showing the LOS results from linear Kalman filter 
%   example
%   >> plotLKFslip(TChaman2Su,Xkp2,22);
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
%Ref=h5read(input_h5_file,'/Ref');
%Rnf=h5read(input_h5_file,'/Rnf');
%
%[los, Base, Blosam, qt, ~, ~, ~] = DataErrorCov(input_h5_file,inputlkf_h5_file,Tri,lambda,k);
figure('position', [300, 500, 320, 700]);colormap(gca,'jet');                      %  ;
%
subaxis(3,1,1,'Spacing', 0.01, 'Padding', 0.02,'MarginTop',0.01,'MarginLeft', 0.15,'MarginRight',0.05,'MarginBottom',0.02);
patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',Xkp(k).smt(1:m),'FaceColor','flat');view(-91,42);...
    caxis([-0.05 0.35]);colormap(gca,'jet');grid on; % xlabel('UTM East (m)');ylabel('UTM North (m)');zlabel('Depth (m)');
%title('Slip estimates (m)') % colorbar; caxis([-0.05 0.35]);
set(gca,'xticklabel',{'250km','255km','260km','265km'},'XGrid','on','Color','none')
set(gca,'yticklabel',{'3380km','3390km','3400km','3410km','3420km','3430km'},'YGrid','on','Color','none')
set(gca,'zticklabel',{'20km','15km','10km','5km','0km'},'ZGrid','on','Color','none')
if k<25
   %title(strcat('Cal (',num2str(dt1),'d)'))
   %text(265000,3420540,2.1,'Cal','FontSize',14)
   text(265000,3395000,-3000,strcat(num2str(dt1),'d'),'FontSize',12)
else
   % title(strcat('Cal (',num2str(dt1),'d/',num2str(dt2),'d)'));
   %text(232200,3420540,2.1,'Cal','FontSize',14)
   text(265000,3395000,-3000,strcat(num2str(dt1),'/',num2str(dt2),'d'),'FontSize',12)
end
%   %
subaxis(3,1,2,'Spacing', 0.01, 'Padding', 0.02,'MarginLeft', 0.15,'MarginRight',0.05,'MarginBottom',0.02);
patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',Xkp(k).smt(2*m+1:3*m),'FaceColor','flat');view(-91,42);...
    caxis([-0.5*10^-4 6*10^-4]);colormap(gca,'jet');grid on; %colorbar; xlabel('UTM East (m)');ylabel('UTM North (m)');zlabel('Depth (m)');
%set(gca,'xticklabel',[],'XGrid','on','Color','none')
set(gca,'xticklabel',{'250km','255km','260km','265km'},'XGrid','on','Color','none')
set(gca,'yticklabel',{'3380km','3390km','3400km','3410km','3420km','3430km'},'YGrid','on')
set(gca,'zticklabel',{'20km','15km','10km','5km','0km'},'ZGrid','on')
if k<25
   %title(strcat('Cal (',num2str(dt1),'d)'))
   %text(265000,3420540,2.1,'Cal','FontSize',14)
   text(265000,3395000,-3000,strcat(num2str(dt1),'d'),'FontSize',12)
else
   % title(strcat('Cal (',num2str(dt1),'d/',num2str(dt2),'d)'));
   %text(232200,3420540,2.1,'Cal','FontSize',14)
   text(265000,3395000,-3000,strcat(num2str(dt1),'/',num2str(dt2),'d'),'FontSize',12)
end
%
subaxis(3,1,3,'Spacing', 0.01, 'Padding', 0.02,'MarginLeft', 0.15,'MarginRight',0.05,'MarginBottom',0.02);
dsigma=dsigma_tde2(Tri,Xkp(k).smt(1:2*m));
patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',dsigma/10^6,...
        'FaceColor','flat');view(-91,42);caxis([0 10]);colormap(gca,'jet');grid on
set(gca,'xticklabel',{'250km','255km','260km','265km'},'XGrid','on','Color','none')
set(gca,'yticklabel',{'3380km','3390km','3400km','3410km','3420km','3430km'},'YGrid','on')
set(gca,'zticklabel',{'20km','15km','10km','5km','0km'},'ZGrid','on')
%title('Slip-rate estimates (m/day)') % caxis([0 5*10^-4]);
%
if k<25
   %title(strcat('Cal (',num2str(dt1),'d)'))
   %text(265000,3420540,2.1,'Cal','FontSize',14)
   text(265000,3395000,-3000,strcat(num2str(dt1),'d'),'FontSize',12)
else
   % title(strcat('Cal (',num2str(dt1),'d/',num2str(dt2),'d)'));
   %text(232200,3420540,2.1,'Cal','FontSize',14)
   text(265000,3395000,-3000,strcat(num2str(dt1),'/',num2str(dt2),'d'),'FontSize',12)
end
saveas(gcf,strcat('~/Desktop/Chaman-S1/Chaman2Su_slip/',num2str(imdates(k)),'slipvel.png'));

