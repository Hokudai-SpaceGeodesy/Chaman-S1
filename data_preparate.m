function [Ref,Rnf,cliparea,Ref0,Rnf0]=data_preparate(cum_filt_h5_file,modelinput_h5file,N,cornerE,cornerN,eu,nu,refc,refr)
% 
% Usage: e.g.
% >>[Ref,Rnf,cliparea]=data_preparate(cum,imdates,23,226192.166,3433135.453,100,100);
%　
cum=h5read(cum_filt_h5_file,'/cum');
[mrow,ncol]=size(cum(:,:,N)');
% UTM coordinate: Corner values are read from UTM.dem_par
%Ed=[213729.143:60:213729.143+60*(ncol-1)];
Ed=[cornerE:eu:cornerE+eu*(ncol-1)];
%Nd=[3346879.239-60*(nrow-1):60:3346879.239];
Nd=[cornerN:-nu:cornerN-nu*(mrow-1)];
[Ref,Rnf]=meshgrid(Ed,Nd);
Ref0=Ed(refc);
Rnf0=Nd(refr);
% First Check the entire region.
%figure(1);surf(Ref,Rnf,cum(:,:,N)');shading flat;view(0,90);colorbar;caxis([-25 25])
figure(1);surf(-0.1*cum(:,:,N)');shading flat;axis ij;view(0,90);colorbar;caxis([-2.5 2.5])
%
% Then, decide to clip out the following 512x512 area.
%
disp('Need to clip out the modeling area, e.g., 512x512 pixels');
x=input('Enter [rmin rmax cmix cmax] ');
%figure(1);surf(Ref(81:592,31:542),Rnf(81:592,31:542),cum(31:542,81:592,31)');shading flat;view(0,90);colorbar;caxis([-25 25]);axis tight
figure(2);surf(Ref(x(1):x(2),x(3):x(4)),Rnf(x(1):x(2),x(3):x(4)),-0.1*cum(x(3):x(4),x(1):x(2),N)');shading flat;view(0,90);colorbar;caxis([-2.5 2.5]);axis tight
Ref=Ref(x(1):x(2),x(3):x(4));
Rnf=Rnf(x(1):x(2),x(3):x(4));
cliparea=x;
%
h5create(modelinput_h5file,'/Ref',size(Ref))
h5write(modelinput_h5file,'/Ref',Ref)
h5create(modelinput_h5file,'/Rnf',size(Rnf))
h5write(modelinput_h5file,'/Rnf',Rnf)
h5create(modelinput_h5file,'/cliparea',[1 4])
h5write(modelinput_h5file,'/cliparea',cliparea)
%
h5create(modelinput_h5file,'/Ref0',size(Ref0))
h5write(modelinput_h5file,'/Ref0',Ref0)
h5create(modelinput_h5file,'/Rnf0',size(Rnf0))
h5write(modelinput_h5file,'/Rnf0',Rnf0)


