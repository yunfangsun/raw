%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This scripts generates the WW3 mesh file from an ADCIRC   %
% fort.14 fil. The user has the option to include the open  %
% bounday nodes or not :                                    %
% meshADC: the ADCIRC fort.14 file                          %
% meshWW3: the resulting WW3 mesh file (gmsh 2.* format     %
% include_boundaries: 1 (include open boundaries)           %
%                     0 (do not include open boundaries)    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

meshADC='./hgrid_WW3_final.gr3';
meshWW3='alaska_noobc.msh';

%meshADC='/scratch2/STI/coastal_temp/save/Panagiotis.Velissariou/ECGC_2023/ECG120.14';
%meshWW3='ECG120.14.msh';

include_boundaries = 0;

%% Read the ADCIRC fort.14 file
[t,p,b,op,bd,title] = readfort14(meshADC, include_boundaries);

tri = t;
x = p(:,1);
y = p(:,2);
h = b;

%%%%%%%%%% Open boundary points %%%%%%%%%%
%% op.nbdv is a sparse matrix, thus we need to call the "full" function
if ~isempty(op) && op.nope > 0
  OB_ID = zeros(op.neta, 1);
  idx2 = 1;
  for nb = 1 : op.nope
    idx1 = idx2;
    idx2 = idx1 + op.nvdll(1);
    OB_ID(idx1:idx2-1) = full(op.nbdv(1:op.nvdll(nb), nb));
  end
else
  OB_ID = zeros(0, 1);
end

WW3_mesh_write(tri, x, y, h, OB_ID, meshWW3, 0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LOCAL FUNCTIONS LOCAL FUNCTIONS LOCAL FUNCTIONS  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EToV,VX,B,opedat,boudat,title] = readfort14( finame, read_bou )
% Read ADCIRC grid file fort.14

disp('read adcirc fort.14') ;
tic
if ( nargin == 0 )
    finputname = 'fort.14'
else
    finputname = finame ;
end

% open file
fid = fopen(finputname) ;

% read first title line
agrid = fgetl(fid) ;
disp(agrid) ;
title = agrid ;

% read number of nodes and elements
msgline = fgetl(fid) ;
N = sscanf(msgline,'%d %d %*s') ;

% read the node numbers, position and depth
Val = fscanf(fid,'%d %g %g %g \n', [4 N(2)])' ;

% read the triangulation
idx = fscanf(fid,'%d %d %d %d %d \n', [5 N(1)])' ;

% sort read items into vectors/matrices
EToV = idx(:,3:5) ;
% treat non-sequential ordering
VX = NaN(max(EToV(:)),2);
B = NaN(max(EToV(:)),1);
VX(Val(:,1),:) = Val(:,2:3) ;
B(Val(:,1)) = Val(:,4) ;

if read_bou
    % Read in boundary
    % Open boundary
    msgline = fgetl(fid) ;
    nope = sscanf(msgline,'%d %*s') ;
    
    msgline = fgetl(fid) ;
    neta = sscanf(msgline,'%d %*s') ;
    
    nvdll = zeros(1,nope) ;
    ibtypee = zeros(1,nope) ;
    nbdv = sparse(neta,nope) ;
    % tic
    for i = 1: nope
        msgline = fgetl(fid) ;
        
        [varg] = sscanf(msgline,'%d %*s \n') ;
        nvdll(i) = varg ;
        ibtypee(i) = 0 ;
        
        %
        % for k = 1: nvdll(i)
        %
        %    % % nbdv(k,i) = fscanf(fid,'%d \n')
        %    % % fscanf(fid,'%d \n')
        %    msgline = fgetl(fid) ;
        %
        %    nbdv(k,i) = str2num(msgline) ;
        % end
        %
        % Nov 25, 2012, improve reading efficiency
        nbdv(1:nvdll(i),i) = fscanf(fid,'%d \n', nvdll(i) ) ;
    end
    % toc
    
    % ocean boundary
    opedat.nope = nope ;
    opedat.neta = neta ;
    opedat.nvdll = nvdll ;
    opedat.ibtypee = ibtypee ;
    opedat.nbdv = nbdv(1:max(nvdll),:);
    
    % land boundary
    msgline = fgetl(fid) ;
    nbou = sscanf(msgline,'%d %*s') ;
    
    msgline = fgetl(fid) ;
    nvel = sscanf(msgline,'%d %*s') ;
    
    nvell = zeros(1,nbou) ;
    ibtype = zeros(1,nbou) ;
    nbvv = sparse(nvel,nbou) ;
    ibconn = sparse(nvel,nbou) ;
    barinht = sparse(nvel,nbou) ;
    barincfsb = sparse(nvel,nbou) ;
    barincfsp = sparse(nvel,nbou) ;
    
    % tic
    for i = 1: nbou
        msgline = fgetl(fid) ;
        
        [varg,~] = sscanf(msgline,'%d %d %*s \n') ;
        nvell(i) = varg(1) ;
        ibtype(i) = varg(2) ;
        
        switch ( ibtype(i) )
            case {0,1,2,10,11,12,20,21,22,30,60,61,101,52}
                % Nov 15, 2012, improve reading efficiency
                nbvv(1:nvell(i),i) = fscanf(fid,'%d \n', nvell(i) ) ;
            case  {3, 13, 23}
                disp('3 13 23')
                val = fscanf(fid,'%g %g %g \n', [3 nvell(i)] )  ;
                nbvv(1:nvell(i),i) = val(1,:) ;
            case  {4, 24}
                %disp('4 24')
                val = fscanf(fid,'%g %g %g %g %g \n', [5 nvell(i)] )  ;
                nbvv(1:nvell(i),i) = val(1,:) ;
                ibconn(1:nvell(i),i) = val(2,:) ;
                barinht(1:nvell(i),i) = val(3,:) ;
                barincfsb(1:nvell(i),i) = val(4,:) ;
                barincfsp(1:nvell(i),i) = val(5,:) ;
            case  {5, 25}
                %disp('5 25')
                val = fscanf(fid,'%g % g %g %g %g %g %g %g \n', [8 nvell(i)] ) ;
                nbvv(1:nvell(i),i) = val(1,:) ;
                %otherwise
                %    msgline = fgetl(fid) ;
            case  94
                val = fscanf(fid,'%d %d \n', [2 nvell(i)] ) ;
                nbvv(1:nvell(i),1:2) = val' ;
        end
    end
    % toc
    
    % land boundary
    boudat.nbou = nbou ;
    boudat.nvel = nvel ;
    boudat.nvell = nvell ;
    boudat.ibtype = ibtype ;
    boudat.nbvv = nbvv(1:max(nvell),:);
    
    if ( sum(ibtype == 24) > 0 ||  sum(ibtype == 4) > 0 )
        boudat.ibconn = ibconn(1:max(nvell),:);
        boudat.barinht = barinht(1:max(nvell),:);
        boudat.barincfsb = barincfsb(1:max(nvell),:);
        boudat.barincfsp = barincfsp(1:max(nvell),:);
    end
else
    opedat = []; 
    boudat  = []; 
end

fclose(fid) ;

return

end %% function close


function [filename]=WW3_mesh_write(tri,x,y,h,OB_ID,filename,plott)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function writes WW3 grid including nodes             %
% (longitude,latitude,depth), open bounday nodes            %
% and element connections (triangles)                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ali Abdolali August 2018 ali.abdolali@noaa.gov       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tri: triangle (nelem,3);
% x: longitude (nnode,1);
% y: laritude (nnode,1);
% h: depth (nnode,1);
% OB_ID: Open boundary nodes ID
% filename: mesh name
% Plot it, if plott=1

node(:,1)=(1:length(x));
node(:,2)=x;
node(:,3)=y;
node(:,4)=h;

fileID = fopen(filename,'w');
fprintf(fileID,'%s\n', '$MeshFormat');
fprintf(fileID,'%s\n', '2 0 8');
fprintf(fileID,'%s\n', '$EndMeshFormat');
fprintf(fileID,'%s\n', '$Nodes');
fprintf(fileID,'%d\n', length(node(:,1)));

for i=1:length(node(:,1))
    fprintf(fileID,['%d %s %10.10f %s %10.10f %s %10.10f\n'], node(i,1),'', node(i,2),'',node(i,3),'',node(i,4));
end
fprintf(fileID,'%s\n', '$EndNodes');
fprintf(fileID,'%s\n', '$Elements');
fprintf(fileID,'%d\n', length(tri(:,1))+length(OB_ID));
m=0;
for i=1:length(OB_ID)
    m=m+1;
    fprintf(fileID,['%d %s %d %s %d %s %d %s %d %s %d\n'], m,'',15,'',2,'',0,'',0,'',OB_ID(i));
end

for i=1:length(tri(:,1))
    m=m+1;
    fprintf(fileID,['%d %s %d %s %d %s %d %s %d %s %d %s %d %s %d %s %d\n'], m,'',2,'',3,'',0,'',i,'',0,'',tri(i,1),'',tri(i,2),'',tri(i,3));
end
fprintf(fileID,'%s', '$EndElements');
fclose(fileID);

if plott==1
width=880;  % Width of figure for movie [pixels]
height=700;  % Height of figure of movie [pixels]
left=700;     % Left margin between figure and screen edge [pixels]
bottom=200;  % Bottom margin between figure and screen edge [pixels]
 


    figure
set(gcf,'Position', [left bottom width height])



trisurf(tri,node(:,2),node(:,3),node(:,4));
shading interp

view(2);
axis equal
hold on
p1=scatter(node(OB_ID,2),node(OB_ID,3),'xk')
hCbar=colorbar
colormap(jet)
legend([p1],'Open Boundary Nodes')
xlim([min(node(:,2))-(max(node(:,2))-min(node(:,2)))/10 max(node(:,2))+(max(node(:,2))-min(node(:,2)))/10])
ylim([min(node(:,3))-(max(node(:,3))-min(node(:,3)))/10 max(node(:,3))+(max(node(:,3))-min(node(:,3)))/10])

if min(node(:,4))~=max(node(:,4))
caxis([min(node(:,4)) max(node(:,4))])
end

hold on
xlabel('Longitude ^{\circ} ');
ylabel('Latitude ^{\circ} ');
axis equal
box on
grid off
axis on
title(['WW3 Grid - depth [m]'],'fontsize',15)

print(gcf,'-dpng',[filename,'.png'],'-r900');

end

return

end %% function close
