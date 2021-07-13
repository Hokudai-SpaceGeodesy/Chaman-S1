function [BRef,BRnf,Bblck,Bwght]=subsample(Ref,Rnf,blck,S)
% SUBSAMPLE: >> [BRef,BRnf,Bblck,Bwght]=subsample(Ref,Rnf,blck,S)
%  S is a sparse matrix derived by "qtdecomp".
[i,j,v]=find(S);
BRef=repmat(1,[length(i) 1]);
BRnf=repmat(1,[length(i) 1]);
Bblck=repmat(1,[length(i) 1]);
Bwght=repmat(1,[length(i) 1]);
for k=1:length(i)
  BRef(k,1)=mean(mean( Ref(i(k):i(k)+v(k)-1,j(k):j(k)+v(k)-1) ));
  BRnf(k,1)=mean(mean( Rnf(i(k):i(k)+v(k)-1,j(k):j(k)+v(k)-1) ));
  Bblck(k,1)=blck(i(k),j(k));
  Bwght(k,1)=S(i(k),j(k));
end
