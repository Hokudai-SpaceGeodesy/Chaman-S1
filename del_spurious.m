function ptde=del_spurious(ptde,Tri,nslimits,bottom)
% del_spurious: 
%    Delete spurious values from ptde, setting both NS-limits and bottom depth.
%  Usage: 
%  >>ptde27_evt1=del_spurious(ptde27,TChaman2,[3.39e+6 3.41e+6],-9000);
ptdeorg=ptde;
%
ulimit=1e-8;
%Ibig=find(ptde>ulimit);
Ibig=find(abs(ptde)>ulimit);
if nslimits(2)<nslimits(1)
    return
end
%
IbigDeep=find(Tri.Points(Tri.ConnectivityList(Ibig,1),3)<bottom & ...
    Tri.Points(Tri.ConnectivityList(Ibig,2),3)<bottom & ...
    Tri.Points(Tri.ConnectivityList(Ibig,3),3)<bottom);
ptde(Ibig(IbigDeep))=0;
%
IbigDeep=find(Tri.Points(Tri.ConnectivityList(Ibig,1),2)>nslimits(2) & ...
    Tri.Points(Tri.ConnectivityList(Ibig,2),2)>nslimits(2) ...
    & Tri.Points(Tri.ConnectivityList(Ibig,3),2)>nslimits(2));
ptde(Ibig(IbigDeep))=0;
%
IbigDeep=find(Tri.Points(Tri.ConnectivityList(Ibig,1),2)<nslimits(1) & ...
    Tri.Points(Tri.ConnectivityList(Ibig,2),2)<nslimits(1) & ...
    Tri.Points(Tri.ConnectivityList(Ibig,3),2)<nslimits(1));
ptde(Ibig(IbigDeep))=0;
%
% For "unconstrained" solution
ptde(ptde<1e-5)=0;
