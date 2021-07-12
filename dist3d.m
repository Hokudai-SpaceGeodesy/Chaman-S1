function h=dist3d(ic,icn)
%function h=dist3d(ic,icn,Tv,V)
n=size(icn,1);
h=[];
%IC3=vert3dplane([V(Tv(1),:);V(Tv(2),:);V(Tv(3),:)],ic(:,1),ic(:,2));
for i=1:n
 %   ICN3=vert3dplane([V(Tv(1),:);V(Tv(2),:);V(Tv(3),:)],icn(i,1),icn(i,2));
    h=[h;sqrt((icn(i,1)-ic(:,1))^2+(icn(i,2)-ic(:,2))^2+(icn(i,3)-ic(:,3))^2)];
end
end

