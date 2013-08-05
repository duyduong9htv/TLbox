function data=transect_effdepavg(x1,y1,x2,y2,dr)
% function data=transect_effdepavg(x1,y1,x2,y2,dr)
%addpath /home/Newmaine2006/Fish_Inversion_Codes/Source


  fid=fopen('bathymetry.arr','rb','ieee-le');
  Xsize = fread(fid,1,'int32');
  Ysize = fread(fid,1,'int32');
  Imagedata = zeros(Xsize,Ysize);
  
  Imagedata=fread(fid,[Xsize,Ysize],'float32');
  fclose(fid);
 
  Imagedata2=rot90(Imagedata); 
  clear Imagedata;
 
 bathymetry
 
  datax=[grid_xmin:grid_inc:grid_xmax];
datay=[grid_ymax:-grid_inc:grid_ymin];
 
 len=sqrt((x1-x2)^2+(y1-y2)^2);
 
 rtemp=[0:dr:round(len)];
 
 phii=atan2((y2-y1),(x2-x1));
 
 xtemp=x1+rtemp*cos(phii);
 ytemp=y1+rtemp*sin(phii);
 
 bathytemp=interp2(datax,datay,Imagedata2,xtemp,ytemp);
 
data=[rtemp.' bathytemp.'];

%hh=find(data(:,2)<0);

%if hh(end)<length(data(:,2));

%data=data(1:hh(end),:);
%end;  
  
%uu=find(data(:,2)<-300);

%hh=isempty(uu);

%if hh==0
%data=data(1:uu(1)-1,:);
%end;

data=[data(:,1) abs(data(:,2))];  
 
 
 
