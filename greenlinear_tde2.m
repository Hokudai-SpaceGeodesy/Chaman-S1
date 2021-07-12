function Basis = greenlinear_tde2( BRef,BRnf,Blosam,TRI,AorD )
%GREENLINEAR_TDE Generate Green function for triangular dislocation 
%   element
%   syntax: >> Basis = greenlinear_tde2( BRef,BRnf,Blosam,TRI,AorD )
%            
%   Here, TRI is generated from >> TRI=gmsh2triangulation('xxx.msh')
%   Within this is called "insarlos_green_tde", which further requires
%   "CaldTriDisps" that was developed by Mead (2007).

% Discard NaN location. Note that BRef and BRnf are in UTM coordinates.
BRef=BRef(~isnan(Blosam));
BRnf=BRnf(~isnan(Blosam));

% TRI.Points are vertices. 
% TRI.Points(:,1)=x-component 
% TRI.Points(:,2)=y-component
% TRI.Points(:,3)=z-component <-- negative downward.
% To use Mead(2007)'s code (CalcTriDisps), however, Vz should be positive downward.
% Thus minus-sign is multiplied to "vz" below.

% Number of TDE equals to the length of "TRI.Triangulation".
%inside=inOutStatus(TRI); % TRI is "constrained"
%tri=TRI.Triangulation(inside,:);
nnl=size(TRI,1);
% Two Bases because of strike and dip slip components.
Bases=repmat(BRef,[1 nnl]);
Based=repmat(BRef,[1 nnl]);

Ia=TRI.ConnectivityList; % indices to i-th triangle
for i=1:nnl
    % For each triangular element, generate impulse response at
    % each observation points and stored in the raw.
    % The factor 10^5 is to convert into centimeter-unit.
    % [ia,ib,ic]=TRI.Triangulation(i,:); 
    vx=[TRI.Points(Ia(i,1),1) TRI.Points(Ia(i,2),1) TRI.Points(Ia(i,3),1)];%
    vy=[TRI.Points(Ia(i,1),2) TRI.Points(Ia(i,2),2) TRI.Points(Ia(i,3),2)];%
    vz=[TRI.Points(Ia(i,1),3) TRI.Points(Ia(i,2),3) TRI.Points(Ia(i,3),3)];%
    % A part of Base for Strike slip. Users must decide Asc or Dsc.
    Bases(:,i)=10^5*insarlos_green_tde(BRef,BRnf,vx,vy,-vz,'S',AorD);
    % A part of Base for Dip slip. Users must decide Asc or Dsc.
    Based(:,i)=10^5*insarlos_green_tde(BRef,BRnf,vx,vy,-vz,'D',AorD);
end

Basis=[Bases Based];

function los=insarlos_green_tde(BRef,BRnf,vertx,verty,vertz,SorD,AorD)

ndata=length(BRef);
los=zeros(length(BRef),1);

if SorD == 'S'
    % Us is 3d-displacement caused by unitary (1m=0.001km) 
    % LEFT-LATERAL (minus-sign) strike slip dislocation
    [Us] = CalcTriDisps(BRef, BRnf, zeros(ndata,1), vertx, ...
    verty, vertz, 0.25, -0.001, 0., 0.);
    % Sign-convention is matched to be the same as "d_insar_qt"
    Us.z=-Us.z;
 if AorD == 'A'
    los=0.6178*Us.x+0.1391*Us.y-0.774*Us.z;
 elseif AorD == 'D' 
    % Below is a special version for azimuth offset. 
    los=-0.174*Us.x+0.985*Us.y-0.00*Us.z;
 elseif AorD == 'E' 
    %Below is the data from ESA
    los=-0.6406*Us.x+0.1182*Us.y-0.7587*Us.z;
 end
 %
%else
%if SorD == 'D'
%    % Ud is 3d-displacement caused by unitary (1m=0.001km)
%    % UP-DIP (minus-sign) slip dislocation
%    [Ud] = CalcTriDisps(BRef, BRnf, zeros(ndata,1), vertx, ...
%    verty, vertz, 0.25, 0., 0., +0.001);
%    % Sign-convention (minus) is matched to be the same as "d_insar_qt"
%    Ud.z=-Ud.z;
% if AorD == 'A'
%    los=0.6178*Ud.x+0.1391*Ud.y-0.774*Ud.z; 
% elseif AorD == 'D'
%    % Below is a special version for azimuth offset. 
%    los=-0.174*Ud.x+0.985*Ud.y-0.00*Ud.z;
% elseif AorD == 'E' 
%    %Below is the data from ESA
%    los=-0.6406*Ud.x+0.1182*Ud.y-0.7587*Ud.z;
% end
%end
end
end

end
