function [los, Base, Blosam, qt, dk, Hk, Rk, pk] = DataErrorCov(input_h5_file,output_h5_file,Tri,lambda,N)
%function [los, Base, Blosam, qt, dk, Hk, Rk, pk] = DataErrorCov(input_h5_file,output_h5_file,Tri,lambda,N)
% Data error covariance matrix
%   詳細説明をここに記述
Ref=h5read(input_h5_file,'/Ref');
Rnf=h5read(input_h5_file,'/Rnf');
filename1=['/',num2str(N),'/los'];
filename2=['/',num2str(N),'/Basea_tde'];
filename3=['/',num2str(N),'/qt'];
los=h5read(output_h5_file,filename1);
Base=h5read(output_h5_file,filename2);
qt=h5read(output_h5_file,filename3);
%
Lt=smoothing_tde2(Tri);
m=size(Tri,1);
[BRef,BRnf,Blosam,~]=subsample(Ref,Rnf,los,qt);
% "Warning" will be issued.
warning off
%  For epoch-by-epoch Least-Squares
p1=(Base'*Base+3000*Lt'*Lt)\(Base'*Blosam(~isnan(Blosam)));
%[p1,~]=lsqnonneg(Base'*Base+7000*Lt'*Lt,Base'*Blosam(~isnan(Blosam)));
% Residual (orig)
resid0 = Blosam(~isnan(Blosam))-Base*p1;
resid = filloutliers(resid0,'linear');
%iresid=subsample_back(qt,Blosam,resid);
%surf(Ref,Rnf,iresid);shading flat;view(0,90);axis tight;colorbar
sigma = std(resid); %
% "pdist" is necessary.
% X=distmat([BRef(~isnan(Blosam)) BRnf(~isnan(Blosam))]);
% "pdist" and "squareform" from Statistics and ML Toolbox
X=squareform(pdist([BRef(~isnan(Blosam)) BRnf(~isnan(Blosam))]));
% Number of Data
n=length(X);
% "ErrorAnalys_Variogram" tells the small "correlation" length of ~500m (or
% bel％ow) but varies from positive to negative...
% -- Constant Lc ---  Both 500m and 2500m were not good.
Cov = sigma^2*exp(-X/1500);
%Cov = exp(-X/1500);
% -- Variable Lc ---  No Good, either.
%if N >= 22 && N <= 24
%    Cov = sigma^2*exp(-X/500); 
%elseif N >=29 && N<= 32
%    Cov = sigma^2*exp(-X/2500);
%else
%    Cov = sigma^2*exp(-X/1500);  %sigma^2*exp(-X/1200);  % 1200  O(sigma^2) ~ 0.1
%end
%  For input to LKF
dk = [Blosam(~isnan(Blosam));zeros(2*m,1);zeros(2*m,1)];
Rk = [Cov zeros(n,2*m) zeros(n,2*m);zeros(2*m,n) lambda^2 * eye(2*m) zeros(2*m);zeros(2*m,n) zeros(2*m) lambda^2 * eye(2*m)];
Hk = [Base zeros(n,2*m);Lt'*Lt zeros(2*m);zeros(2*m) Lt'*Lt];  % 
%  For epoch-by-epoch Least-Squares
pk = lsqnonneg(Base'*Cov^(-1)*Base + lambda*Lt'*Lt,Base'*Cov^(-1)*Blosam(~isnan(Blosam))); %% Seems best
% pk = lsqnonneg(Base'*Base + lambda*Lt'*Lt,Base'*Blosam(~isnan(Blosam)));
% pk = (Base'*Cov^(-1)*Base + lambda*Lt'*Lt)\(Base'*Cov^(-1)*Blosam(~isnan(Blosam)));

%Hk = [Base zeros(n,2*m);lambda*Lt'*Lt zeros(2*m);zeros(2*m) lambda*Lt'*Lt];
% variogram computation based on m-file by Wolfgang Schwanghart
%d = variogram([BRef(~isnan(Blosam)) BRnf(~isnan(Blosam))],resid,'plotit',false,'nrbins',100);
% Covariance(r) = sigma^2 - d.val(semivariogram);
%figure('position', [300, 500, 900, 300])
%subplot(121);plot(d.distance,sigma^2-d.val,'o-')
%subplot(122);surf(Ref,Rnf,iresid);shading flat;view(0,90);axis tight;colorbar
% Fitting "d.distance vs d.val" with C(r)=Sigma*exp(-d/L) to derive Sigma and L
% Fitting "d.distance vs d.val" with C(r)=Sigma*exp(-d/L) to derive Sigma and L was IMPOSSIBLE
% because of the negative values. 
end

