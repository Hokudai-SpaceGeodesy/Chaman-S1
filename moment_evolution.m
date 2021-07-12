function [Mom1,Mom2]=moment_evolution(cum_filt_h5_file,Xkp,nslimit1,nslimit2,bottom,Tri,N1,N2)
%function [Mom1,Mom2]=moment_evolution(cum_filt_h5_file,output_h5_file,Tri,N1,N2)
% Compute and plot how Moment evolved
%

% Define Year. Month. Day of the First and 2nd EQK
yr0=2016;mo0=5;dd0=13;mo1=7;dd1=10;
% Co-seismic Moment by USGS
Mom_Evt1_Coseismic=3.1125e+17;
Mom_Evt2_Coseismic=1.4125e+16;
imdates=h5read(cum_filt_h5_file,'/imdates');
Mom1=[];Mom2=[];
nparam=size(Tri,1)*2;
for i=N1:N2
    %yr=floor(imdates(i)/10000);mmdd=imdates(i)-yr*10000;
    %mo=floor(mmdd/100);dd=mmdd-mo*100;
    obdate=num2str(imdates(i));
    Yr=2016+((datenum(obdate,'yyyymmdd')-datenum(2016,1,1))/365.25);
    %[yr mo dd]
%    p1=h5read(output_h5_file,strcat('/',num2str(i),'/p1tde_evt1'));
    ptmp=Xkp(i).smt(1:nparam);
    %ptmp=Xkp(i).lsq(1:nparam);
    p1=del_spurious(ptmp,Tri,nslimit1,bottom);
%   subplot(211);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',p1(1:1/2*length(p1),1),...
%   'FaceColor','flat');view(-91,42);colormap(gca,'jet');colorbar
    [~,Mom]=momentmag_tde2(Tri,p1);
%   ddif=datenum(double(yr),double(mo),double(dd))-datenum(yr0,mo0,dd0);
    ddif=datenum(obdate,'yyyymmdd')-datenum(yr0,mo0,dd0);
    %Yr=double(yr)+(double(mo-1)*30+double(dd))/365.25;
    Mom1=[Mom1;Yr ddif Mom];
%    p2=h5read(output_h5_file,strcat('/',num2str(i),'/p1tde_evt2'));
    p2=del_spurious(ptmp,Tri,nslimit2,bottom);
%    subplot(212);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',p2(1:1/2*length(p2),1),...
%    'FaceColor','flat');view(-91,42);colormap(gca,'jet');colorbar
    [~,Mom]=momentmag_tde2(Tri,p2);
    Mom2=[Mom2;Yr ddif Mom];
end
%figure(2);
plot(Mom1(:,1),Mom1(:,3),'r*','LineWidth',2,'MarkerSize',6);
%hold on;plot(Mom2(4:39,1),Mom2(4:39,3),'ro','LineWidth',2,'MarkerSize',8);
hold on;plot(Mom2(:,1),Mom2(:,3),'b*','LineWidth',2,'MarkerSize',6);
%hold on;plot(yr0+((mo0-1)*30+dd0)/365.25,Mom_Evt1_Coseismic,'bx','LineWidth',2,'MarkerSize',8);
%hold on;plot(yr0+((mo1-1)*30+dd1)/365.25,Mom_Evt2_Coseismic,'rx','LineWidth',2,'MarkerSize',8);
%hold on;plot(Mom2(1:3,1),Mom2(1:3,3),'r+','LineWidth',2,'MarkerSize',6); 
grid on;

%plot(Mom1(:,1),Mom1(:,2),'bo');
%hold on;plot(Mom2(:,1),Mom2(:,2),'ro');

%hold on;plot(datenum(yr,mo,dd)-datenum(yr0,mo0,dd0),Mom,'ro')
%hold on;plot(datenum(yr,mo,dd)-datenum(yr0,mo0,dd0),Mom,'ro')