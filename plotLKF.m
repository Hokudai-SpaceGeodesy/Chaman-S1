function plotLKF(input_h5_file,inputlkf_h5_file,Tri,Xkp,Okp,lambda,Ni,Ne)
% plotLKF: showing the results of linear Kalman filter 
%   example
%   >> plotLKF('modelinput.h5','input_for_lkf_Su.h5',TChaman2Su,Xkp2,Okp2,0.1,22,35);
%
%   Masato Furuya (c) 2021-
imdates=h5read('TS_GEOCml2clip.UTM/cum_filt.h5','/imdates');
N=length(imdates);Yr=zeros(N,1);
for i=1:N
  Yr(i)=datenum(num2str(imdates(i)),'yyyymmdd');
end
% Time-difference vector in "days"
dt=diff(Yr);
% Number of triangular meshes
m=size(Tri,1);
Ref=h5read(input_h5_file,'/Ref');
Rnf=h5read(input_h5_file,'/Rnf');
% Minus2L = 0;
% logVk = 0;
for k=Ni:Ne %N
    [los, Base, Blosam, qt, ~, ~, ~] = DataErrorCov(input_h5_file,inputlkf_h5_file,Tri,lambda,k);
    %
    figure('position', [300, 500, 1500, 1100]);colormap(gca,'Parula');
    subplot(331);surf(Ref,Rnf,los);view(0,90);shading flat;axis tight;caxis([-2.5 2.5]);colorbar;title(strcat(num2str(imdates(k)),' (Obs, cm)'));
    subplot(331);xlabel('UTM East (m)');ylabel('UTM North (m)')
    %
    %ical=subsample_back(qt,Blosam,Base*Xkp(k).upd(1:2*m));
    %subplot(321);surf(Ref,Rnf,ical);view(0,90);shading flat;axis tight;caxis([-2.5 2.5]);colorbar;title(strcat(num2str(imdates(k)),' (Cal)'))
    %subplot(321);xlabel('UTM East (m)');ylabel('UTM North (m)')
    %                                          lsq                                             log10(Xkp(1:1/2*length(Xkp)))
    ical2=subsample_back(qt,Blosam,Base*Xkp(k).smt(1:2*m));
    subplot(332);surf(Ref,Rnf,ical2);view(0,90);shading flat;axis tight;caxis([-2.5 2.5]);colorbar;title(strcat(num2str(imdates(k)),' (Cal, cm)'))
    subplot(332);xlabel('UTM East (m)');%ylabel('UTM North (m)')
    %
    subplot(333);surf(Ref,Rnf,los-ical2);view(0,90);shading flat;axis tight;caxis([-1.5 1.5]);colorbar;title(strcat(num2str(imdates(k)),' (Misfit, cm)'))
    subplot(333);xlabel('UTM East (m)');%ylabel('UTM North (m)')
    %                                                                                       log10(Xkp(1:1/2*length(Xkp)))
    %subplot(323);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',Xkp(k).upd(1:m),...
    %'FaceColor','flat');view(-91,42);caxis([-0.05 0.35]);colorbar;colormap(gca,'jet');grid on;xlabel('UTM East (m)');ylabel('UTM North (m)');zlabel('Depth (m)');
    %subplot(323);title('Slip distribution estimates')
    %
    %subplot(324);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',Xkp(k).upd(2*m+1:3*m),...
    %'FaceColor','flat');view(-91,42);colorbar;colormap(gca,'jet');grid on;xlabel('UTM East (m)');ylabel('UTM North (m)');zlabel('Depth (m)');
    %subplot(324);title('Slip-rate distribution estimates')
    %
    subplot(334);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',Xkp(k).smt(1:m),...
    'FaceColor','flat');view(-91,42);caxis([-0.05 0.35]);colorbar;colormap(gca,'jet');grid on;xlabel('UTM East (m)');ylabel('UTM North (m)');zlabel('Depth (m)');
    subplot(334);title('Slip estimates (m)') % caxis([-0.05 0.35]);
%   %
    subplot(335);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',Xkp(k).smt(2*m+1:3*m),...
    'FaceColor','flat');view(-91,42);caxis([-0.5*10^-4 6*10^-4]);colorbar;colormap(gca,'jet');grid on;xlabel('UTM East (m)');ylabel('UTM North (m)');zlabel('Depth (m)');
    subplot(335);title('Slip-rate estimates (m/day)') % caxis([0 5*10^-4]);
    %
    %dsigma0=dsigma_tde2(Tri,X0);
    %dsigma1=dsigma_tde2(Tri,X1);
    dsigma=dsigma_tde2(Tri,Xkp(k).smt(1:2*m));
    subplot(336);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',dsigma/10^6,...
        'FaceColor','flat');view(-91,42);caxis([0 10]);colorbar;colormap(gca,'jet');grid on;ylabel('UTM North (m)');
    subplot(336);title('Static stress drop estimates (MPa)')  % 
    %
    subplot(337);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',2*diag(Okp(k).smt(1:m,1:m)),...
    'FaceColor','flat');view(-91,42);caxis([0 0.03]);colorbar;colormap(gca,'jet');grid on;xlabel('UTM East (m)');ylabel('UTM North (m)');zlabel('Depth (m)');
    subplot(337);title('2\sigma slip error (m)')  %
    %subplot(337);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',Xkp(k).lsq(1:m),...
    %'FaceColor','flat');view(-91,42);caxis([-0.05 0.35]);colorbar;colormap(gca,'jet');grid on;xlabel('UTM East (m)');ylabel('UTM North (m)');zlabel('Depth (m)');
    %subplot(337);title('LSQ Slip estimates (m)') % caxis([-0.05 0.35]);
    %
    subplot(338);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',2*diag(Okp(k).smt(2*m+1:3*m,2*m+1:3*m)),...
    'FaceColor','flat');view(-91,42);colorbar;colormap(gca,'jet');grid on;xlabel('UTM East (m)');ylabel('UTM North (m)');zlabel('Depth (m)');
    subplot(338);title('2\sigma slip-rate error (m/day)')  %caxis([0 10^-6]);
    %ical3=subsample_back(qt,Blosam,Base*Xkp(k).lsq(1:2*m));
    %subplot(338);surf(Ref,Rnf,ical3);view(0,90);shading flat;axis tight;caxis([-2.5 2.5]);colorbar;title(strcat(num2str(imdates(k)),' (Cal, cm)'))
    %subplot(338);xlabel('UTM East (m)');%ylabel('UTM North (m)')
    %
    dsigmaerr=dsigma_tde2(Tri,diag(Okp(k).smt(1:2*m,1:2*m)));
    subplot(339);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',2*dsigmaerr/10^6,...
        'FaceColor','flat');view(-91,42);caxis([0 1.0]);colorbar;colormap(gca,'jet');grid on;ylabel('UTM North (m)');
    subplot(339);title('2\sigma stress drop error (MPa)')  % 
    %subplot(339);surf(Ref,Rnf,los-ical3);view(0,90);shading flat;axis tight;caxis([-1.5 1.5]);colorbar;title(strcat(num2str(imdates(k)),' (LSQ Misfit, cm)'))
    %subplot(339);xlabel('UTM East (m)');%ylabel('UTM North (m)')
   % hp3 = get(subplot(2,3,6),'Position');
   % c=colorbar('Position', [hp3(1)+hp3(3)+0.01  hp3(2)  0.01  hp3(3)]);
   % c.Label.String = 'MPa';
   % c.Label.FontName='Helvetica';
   % c.Label.FontSize = 12;
   % Minus2L = (dk - Hk * Xkp(k).pre)'*((Rk + Hk * Okp(k).pre * Hk')^(-1))*(dk - Hk * Xkp(k).pre);
   % logVk = logdet(Rk + Hk * Okp(k).pre * Hk');
    %[k Minus2L logVk]
    refreshdata;drawnow
%    subplot(221);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',Xkp(k).upd(1:m),...
%    'FaceColor','flat');view(-91,42);caxis([-0.1 0.5]);colorbar;colormap(gca,'jet');grid on;xlabel('UTM East (m)');ylabel('UTM North (m)');zlabel('Depth (m)');
%    subplot(221);title('Slip distribution estimates')
%    %
%    subplot(222);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',Xkp(k).upd(2*m+1:3*m),...
%    'FaceColor','flat');view(-91,42);colorbar;colormap(gca,'jet');grid on;xlabel('UTM East (m)');ylabel('UTM North (m)');zlabel('Depth (m)');
%    subplot(222);title('Slip-rate distribution estimates')
%    %
%    subplot(223);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',Xkp(k).smt(1:m),...
%    'FaceColor','flat');view(-91,42);caxis([-0.1 0.5]);colorbar;colormap(gca,'jet');grid on;xlabel('UTM East (m)');ylabel('UTM North (m)');zlabel('Depth (m)');
%    subplot(221);title('Slip distribution estimates')
%    %
%    subplot(224);patch('Faces',Tri.ConnectivityList,'Vertices',Tri.Points,'FaceVertexCData',Xkp(k).smt(2*m+1:3*m),...
%    'FaceColor','flat');view(-91,42);colorbar;colormap(gca,'jet');grid on;xlabel('UTM East (m)');ylabel('UTM North (m)');zlabel('Depth (m)');
%    subplot(224);title('Slip-rate distribution estimates')
%    refreshdata;drawnow

end

