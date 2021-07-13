function d = smoothing_tde2( TRI )
% SMOOTHING_TDE provides with Laplacian operator to smooth fault slip 
%   distribution 
%   Input TRI should be generated from gmsh2triangulation (instead of
%   gmsh2trirep)
nall=size(TRI,1);
d=zeros(nall);
iv=[1:nall];

for i=[iv]
    % For the i-th TDE, ...
    ic=incenter(TRI,i); % Center coordinate of the i-th TDE 
    n1=neighbors(TRI,i); % Tell us Neighbors of the i-th TDE, including NaN  
    icn=incenter(TRI,n1(~isnan(n1))'); % Center coordinates of the neighbors 
    %
    hi=dist3d(ic,icn);
    L=sum(hi); %%% 
    d(i,i)=(-2/L)*sum(1./hi);
    %
    for j=[ n1(~isnan(n1)) ]
    %    jnum= iv==j; % (Fixed by matlab)
    %    d(inum,jnum)=(2/L)/dist3d(ic,incenters(TRI,j),Tv,V);    
        d(i,j)=(2/L)/dist3d(ic,incenter(TRI,j));    
    end
end
d=[d zeros(nall);zeros(nall) d];
scal=(max(max(d))-min(min(d)))/2;
d=d/scal;
end
