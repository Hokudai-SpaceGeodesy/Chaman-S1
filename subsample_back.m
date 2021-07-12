function dblck_cal=subsample_back(S,Bdnov05m,Rdat_cal)
% SUBSAMPLE_BACK
%
blck_cal=repmat(Bdnov05m,[1 1]);
blck_cal(~isnan(blck_cal))=Rdat_cal;
%
dblck_cal=full(S);
[i,j,v]=find(S);
for k=1:length(i)
  dblck_cal(i(k):i(k)+v(k)-1,j(k):j(k)+v(k)-1)=blck_cal(k);
end
