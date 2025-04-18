clc
clear

clc
clear

files={
'FRF-met_wind_derived_201210.nc';	
'FRF-met_barometer_derivedBarom_201210.nc';				
'FRF-ocean_currents_awac-6m_201210.nc'		;	
'FRF-ocean_currents_awac-8m_201210.nc'		;		
'FRF-ocean_waterlevel_eopNoaaTide_201210.nc'	;	
'FRF-ocean_waves_adop-2m_201210.nc';
'FRF-ocean_waves_awac-6m_201210.nc';
'FRF-ocean_waves_awac-8m_201210.nc';
%'FRF-ocean_waves_8m-array_201210.nc'			;
'FRF-ocean_waves_waverider-17m_201210.nc';
'FRF-ocean_waves_waverider-26m_201210.nc'};


%%
for i=3:4 ;

lon(i)=ncread(files{i},'longitude');
lat(i)=ncread(files{i},'latitude');



time{i}=ncread(files{i},'time');
%station_name{i}=ncread(files{i},'station_name');

v{i}=ncread(files{i},'currentNorth');
u{i}=ncread(files{i},'currentEast');

us{i}=ncread(files{i},'currentSpeed');
ud{i}=ncread(files{i},'currentDirection');


h{i}=ncread(files{i},'depth');


end

%lon=-abs(lon);


% for i=5:8;
% station_name{i}=ncread(files{i},'station_name');
% h{i}=ncread(files{i},'waterLevel');
% %station_name{i}=ncread(files{i},'station_name');
% end


lx=lon(3:4);
ly=lat(3:4);

for i=1:2
tn{i}=time{i+2};
un{i}=u{i+2};
vn{i}=v{i+2};
hn{i}=h{i+2};

usn{i}=us{i+2};
udn{i}=ud{i+2};

end


clear time u v h



for i=1:2
time{i}=tn{i}/3600./24.+datenum([1970 1 1])-datenum([2012 1 1])+1;
end

clear tn

sch='../wwm/schout_1_LON.nc';

sw2d='../wwm/schout_1_LON.nc';
sw3d='../wwm/schout_1_VOR.nc';



x=ncread(sw3d,'SCHISM_hgrid_node_x');
y=ncread(sw3d,'SCHISM_hgrid_node_y');
sig=ncread(sw3d,'sigma');

t1=ncread(sch,'time');
t2=ncread(sw2d,'time');
t3=ncread(sw3d,'time');

dep=ncread(sw3d,'depth');

xy_all=[x, y];


hvel3=ncread(sw3d,'hvel');
hvel2=ncread(sw2d,'hvel');
hvel1=ncread(sch,'hvel');




id=load('line.txt');

x1=x(id);
y1=y(id);

x2 = x1(end:-1:1);
y2 = y1(end:-1:1);





obs_w=[x2, y2];


idx = knnsearch(xy_all, obs_w);


[n1 n2 n3 n4]=size(hvel1);

dis(1)=0;
dd(1)=0;
for i=2:length(x2)
    dis(i)=sqrt((x2(i)-x2(i-1))^2+(y2(i)-y2(i-1))^2);
    dd(i)=dis(i)+dd(i-1);
end
    
dd=dd*110;

h=dep(idx);


%plot(dd,-h)




%w1='../ww3.201210.nc'
%w2= '../ali/duck_case_1/ww3.201210.nc';



%for i=1:2
    
u_3d=squeeze(hvel3(1,:,idx,:));
v_3d=squeeze(hvel3(2,:,idx,:));

u_2d=squeeze(hvel2(1,:,idx,:));
v_2d=squeeze(hvel2(2,:,idx,:));

u_1=squeeze(hvel1(1,:,idx,:));
v_1=squeeze(hvel1(2,:,idx,:));

%end


a=(90-71.2)*pi/180;

u3d_c=u_3d*cos(a)+v_3d*sin(a);
u3d_a=v_3d*cos(a)-u_3d*sin(a);


u2d_c=u_2d*cos(a)+v_2d*sin(a);
u2d_a=v_2d*cos(a)-u_2d*sin(a);


u1d_c=u_1*cos(a)+v_1*sin(a);
u1d_a=v_1*cos(a)-u_1*sin(a);




%sin(30/180.*pi)



clear hvel*


t1=t1/3600./24.+301.;
t2=t2/3600./24.+301.;
t3=t3/3600./24.+301.;

close all;




for i=1:length(h)
    z(i,:)=h(i)*sig;
xx(i,1:31)=dd(i);
end


clim_a = [-1.5 0.5];
clim_c = [-0.5 0.5];





k=5

 
 
 

for k=1:length(t2)
    
    figure
    
    set(gcf, 'PaperUnits', 'inches', 'PaperOrientation', 'landscape', 'PaperPosition', [0 0 11 8.5]);

    subplot(2,2,1)
    
    
      colormap(jet);

pcolor(xx',z',u2d_a(:,:,k));


shading flat
   caxis(clim_a);

xlim([0,6.5])

ylim([-20,0])

text(0.1,-18,'Along-shore 2D Coupling','FontSize',12);

 colorbar;
 
 title(datestr(t2(k)));
    
 
     subplot(2,2,2)
    
    
      colormap(jet);

pcolor(xx',z',u2d_c(:,:,k));


shading flat
   caxis(clim_c);

xlim([0,6.5])

ylim([-20,0])
text(0.1,-18,'Cross-shore 2D Coupling','FontSize',12);
 colorbar;
 
     subplot(2,2,3)
    
    
      colormap(jet);

pcolor(xx',z',u3d_a(:,:,k));


shading flat
   caxis(clim_a);

xlim([0,6.5])

ylim([-20,0])

 colorbar;
 text(0.1,-18,'Along-shore 3D Coupling','FontSize',12);   
 
     subplot(2,2,4)
    
    
      colormap(jet);

pcolor(xx',z',u3d_c(:,:,k));


shading flat
   caxis(clim_c);

xlim([0,6.5])

ylim([-20,0])

 colorbar;
 
 
 text(0.1,-18,'Cross-shore 3D Coupling','FontSize',12);
 
    
    
    output_name= [sprintf('current/fig%2d',k) 'pic.png'];
          print('-dpng','-r140', output_name);
    close all
end








