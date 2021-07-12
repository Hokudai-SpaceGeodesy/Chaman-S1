function [Xkp,Okp,Minus2L,Smisfit,Smooth,Snorm] = LKF(input_h5_file,inputlkf_h5_file,Tri,lambda,alpha,Ni,Ne,X0,X1)
%        [Xkp,Okp,Minus2L,Smooth,Snorm]
% "Time-dependent inversion" by Linear Kalman Filter
%   Example  
%   >> [Xkp2,Okp2,Minus2Lx,Smisfitx,Smoothx,Snormx]=LKF('modelinput.h5','input_for_lkf_Su.h5',TChaman2Su,10^(-2),10^(-5),22,74,X0b+X0c,X1b);
%   where 'modelinput.h5' and 'input_for_lkf_Su.h5' are input hdf files.
%         TChaman2Su is the triangular class. Lambda and Alpha are spatial and temporal damping, respectively.
%         We need to assign the initial epoch Ni and final epoch Ne in the "imdates".
%         X0 and X1 are the initial state vectors for Event 1 and Event 2,
%         respectively.
%   Masato Furuya (c) 2021-      
%
imdates=h5read('TS_GEOCml2clip.UTM/cum_filt.h5','/imdates');
N=length(imdates);Yr=zeros(N,1);
for i=1:N
  % --SILLY--
  %  yr=floor(imdates(i)/10000);mmdd=imdates(i)-yr*10000;
  %  mo=floor(mmdd/100);dd=mmdd-mo*100;
  %  Yr(i)=double(yr)+(double(mo-1)*30+double(dd))/365.25;
  Yr(i)=datenum(num2str(imdates(i)),'yyyymmdd');
end
% Time-difference vector in "days"
dt=diff(Yr);
% Number of triangular meshes
m=size(Tri,1);
% Initializing State vector and Prediction-error variance-covariane matrix
k=Ni-1;
% Xkp(k).upd = [zeros(2*m,1);zeros(2*m,1)]; % when started from no-deformation
Xkp(k).upd = X0; % <-- Co-seismic slip and slip-rate distributions (m) are necessary to set.  !!!!
% "How is X0 derived?"
%  p22tde=h5read('output_file.h5','/22/p1tde'); 
%  p22tde=del_spurious(p22tde,TChaman2S,nlimit_evt1,-9000);
%  X0=[p22tde;zeros(1572,1)];
% Initial "slip-rate" is also important. (Mar 20, 2021)
%  X0=[p22tde;0.0008*p22tde(1:786);zeros(786,1)];
% Revised Initial State (Apr 29)
% [p22nng,~]=lsqnonneg(Base22'*Base22+3000*Lt2'*Lt2,Base22'*Blosam22(~isnan(Blosam22)));
% p22Su=del_spurious(p22nng,TChaman2Su,[3392000.0 3407000.0],-10000);
% X0b=[p22Su;0.003*p22Su(1:519);zeros(519,1)]; <-- This includes only Event 1 but there occured a triggered creep event near the location of Event 2.
%   In "p22nng" there's no signal around the triggered slip event.
% X0c=[0.7*p23Su;0.002*p23Su(1:519);zeros(519,1)]; <-- This accounts for the triggered creep event. The factor 0.5 is "conservative".  
% 
S0 = diag(0.005*[ones(1,m) zeros(1,m)]);
Okp(k).upd = [S0 zeros(2*m);zeros(2*m) 5*10^(-8)*eye(2*m)]; % Uncertainies of 0.5cm strike-slip and 10^-8 m/d slip-rate  !!!! 10^-7 gives a growing moment
%Okp(k).upd = [0.005*eye(2*m) zeros(2*m);zeros(2*m) 10^(-8)*eye(2*m)]; % Uncertainies of 0.5cm strike-slip and 10^-7 m/d slip-rate  !!!!
%%lambda = 1e+4;alpha = 2e+4;
Ref=h5read(input_h5_file,'/Ref');
Rnf=h5read(input_h5_file,'/Rnf');
% Likelihood & Trade-off
Minus2L=0;Smisfit=0;Smooth=0;Snorm=0;  
logVk=0;Nd = 0;
%
Lt=smoothing_tde2(Tri);
% How is X1 derived?   "Xkp(k).pre" should be used for Event 1 but we should provide adequate co-seismic offset for Event 2.  
%  p25tde=h5read('output_file.h5','/25/p1tde'); 
%  p25tde=del_spurious(p25tde,TChaman2S,[3412000 3430000],-9000);
%  (xxx) X1=[p25tde;zeros(1572,1)];
%  X1=[p25tde;0.001*p25tde(1:786);zeros(786,1)];  %Initial "slip-rate" is also important. 0.003 is too high (Mar 20, 2021)
%Ie2 = find(X1 ~= 0);
% ---! But X1 should not be used as it is.
%  [p25nng,~]=lsqnonneg(Base25'*Base25+3000*Lt2'*Lt2,Base25'*Blosam25(~isnan(Blosam25)));
%  p25Su=del_spurious(p25nng,TChaman2Su,[3412000.0 3430000.0],-10000); 
%  X1b=[p25Su;0.002*p25Su(1:519);zeros(519,1)]; %% <--- This includes the triggered slip by the E1 before the onset of E2. 
%
%  Filtering
for k=Ni:Ne
    % Transition matrix formation
    Tk = [eye(2*m) dt(k)*eye(2*m);zeros(2*m) eye(2*m)];
    % Process noise variance-covariance matrix
    Qk = alpha^2*[dt(k)^3*eye(2*m)/3 dt(k)^2*eye(2*m)/2;dt(k)^2*eye(2*m)/2 dt(k)*eye(2*m)]; 
    % Prediction
    Xkp(k).pre = Tk * Xkp(k-1).upd;
    Okp(k).pre = Tk * Okp(k-1).upd * Tk' + Qk;
    if k == 25
      Xkp(k).pre = Xkp(k).pre + X1; % Add pre(?)and Co-seismic for the 2nd event  
    end
    %
    [~, Base, Blosam, ~, dk, Hk, Rk, pk] = DataErrorCov(input_h5_file,inputlkf_h5_file,Tri,lambda,k);
    % Epoch-by-epoch lsq. solution
    Xkp(k).lsq = pk;  
    % Update with Kalman gain, "Kak", and "inovation", resulting in unconstrained estimates.
    Kak = Okp(k).pre * Hk' * (Rk + Hk * Okp(k).pre * Hk')^(-1);% Linear Kalman Filter
    Xkp(k).upd = Xkp(k).pre + Kak * (dk - Hk * Xkp(k).pre);
    Okp(k).upd = Okp(k).pre - Kak * Hk * Okp(k).pre;
    %%  For likelihood computation : Does NOT work nicely (visually/intuitively)
    Nd = Nd + length(dk);
    Minus2L = Minus2L + (dk - Hk * Xkp(k).pre)'*((Rk + Hk * Okp(k).pre * Hk')^(-1))*(dk - Hk * Xkp(k).pre); 
    % --Determinant computation needs to be carefully done. The "logdet" is helpful. --
    % Without Rk, logdet is negative. When lambda is ~O(2) or greater, logVk becomes positive.
    logVk = logVk + logdet(Rk + Hk * Okp(k).pre * Hk');
%   Optimization
    resid = Blosam(~isnan(Blosam)) - Base*Xkp(k).upd(1:2*m); n=length(resid);
    Smisfit = Smisfit + resid'*Rk(1:n,1:n)^(-1)*resid; %resid'*resid; %norm(resid);
    Smooth = Smooth + (diag([Lt zeros(2*m);zeros(2*m) Lt]).*Xkp(k).upd)'*(diag([Lt zeros(2*m);zeros(2*m) Lt]).*Xkp(k).upd); %norm(diag([Lt zeros(2*m);zeros(2*m) Lt]).*Xkp(k).smt); 
    Snorm = Snorm + Xkp(k).upd'*Xkp(k).upd; %norm(Xkp(k).smt);
end
%  Loglikelihood (??)
    %sigmafactor = sqrt(Minus2L/Nd);
    %Minus2L = Nd - Nd*log(Nd) + logVk + Nd * log(Minus2L); % Minus2L; %;
    %[Minus2L logVk]
    Minus2L = logVk + Minus2L;
%  Last Epoch
    Xkp(Ne).smt = Xkp(Ne).upd;
    Okp(Ne).smt = Okp(Ne).upd;
%  Smoothing
for k=Ne-1:-1:Ni
    % Transition matrix formation
    Tk = [eye(2*m) dt(k)*eye(2*m);zeros(2*m) eye(2*m)];
    % Smoothing Matrix
    Sk = Okp(k).upd * Tk' * Okp(k).pre^(-1);
    %
    Xkp(k).smt = Xkp(k).upd + Sk * ( Xkp(k+1).upd - Xkp(k+1).pre );
    Okp(k).smt = Okp(k).upd + Sk * ( Okp(k+1).upd - Okp(k+1).pre ) * Sk';
end
%

    
    
