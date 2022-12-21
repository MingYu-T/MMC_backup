%% Transfer variable to Topology to APDL command %%
clear;clc;close all
% - load var 
var = load('var_opt.txt');

% - Initialize parameters
numCPoint = 8;
numSpline = 12;
Var_num = (numCPoint-2)*2+numCPoint;
ord = 3;
E0=1e8;
nu=0.3;
DW = 100e-3;
DH = 50e-3;
thickness = 1e-3;
nelx = 100;
nely = 50;
EW=DW/nelx; EH=DH/nely; 
force = 5;
Kin = 500;
Kout = 100;
alpha=1e-3;              % parameter alpha in the Heaviside function
epsilon=0.1*min(EW,EH);  % regularization parameter epsilon in the Heaviside function

c_ends = zeros(numSpline,4);
for i = 1:4
    c_ends(i,1) = 0;
    c_ends(i,2) = 0;
    c_ends(i,3) = DW;
    c_ends(i,4) = DH;
end
for i = 5:8
    c_ends(i,1) = 0;
    c_ends(i,2) = DH;
    c_ends(i,3) = DW;
    c_ends(i,4) = DH;
end
for i = 9:12
    c_ends(i,1) = 0;
    c_ends(i,2) = 0;
    c_ends(i,3) = 0;
    c_ends(i,4) = DH;
end

% - Initialize Domain
[i0,j0] = meshgrid(0:nelx,0:nely);
gobalcoords0x=i0*EW;
gobalcoords0y=EH*nely-j0*EH;
gobalcoords0xy=[gobalcoords0x(:) gobalcoords0y(:)]';
gobalcoords0=gobalcoords0xy(:);
[il,jl] = meshgrid(0,nely);
loadnid = il*(nely+1)+(nely+1-jl); 
loaddof = 2*loadnid(:)-1; 
[ilo,jlo] = meshgrid(nelx,nely);
outnid = ilo*(nely+1)+(nely+1-jlo);  
outdof = 2*outnid(:)-1; 
[iif,jf] = meshgrid(0,0:nely/5);
fixednid = iif*(nely+1)+(nely+1-jf);
fixeddof = [2*fixednid(:); 2*fixednid(:)-1];
[iif,jf] = meshgrid(0:nelx,nely);
fixednid = iif*(nely+1)+(nely+1-jf);
fixeddof = [2*fixednid(:); fixeddof];
nele = nelx*nely;
ndof = 2*(nelx+1)*(nely+1);
F = sparse(loaddof,1,force,ndof,1);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);                               
edofVec = 2*nodeids(:)+1;                                                              
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*(nely+1) 2*(nely+1)+1 2*(nely+1)-2 2*(nely+1)-1 -2 -1],nele,1);
edotMat=0.5*edofMat(:,[2,4,6,8]);

% - Forming Phi^s
Phi=cell(numSpline,1);
for i = 1:numSpline
    Phi{i} = TDF(var(Var_num*i-Var_num+1:Var_num*i),numCPoint,c_ends(i,:),DH,DW,nelx,nely);
end
Phi_s = Phi{1};
for i = 2:numSpline
    Phi_s = max(Phi_s,Phi{i});
end

 % - Heaviside Function
 H = Heaviside(Phi_s,alpha,nelx,nely,epsilon);
 den=sum( H(edotMat), 2 ) / 4;
 fval = sum(den)*(EW*EH)/(DW*DH) ;
 disp(fval)
 % - Element Density
 temp1 = zeros(nele,1); 
 for i = 1:nele
    temp1(i) = ( H(edotMat(i,1))^2+H(edotMat(i,2))^2+H(edotMat(i,3))^2+H(edotMat(i,4))^2 )/4;  % ersatz material model
 end
 xba = reshape(temp1,nely,nelx);
 xba = flip(xba);
 xba = xba(:);

% - Commands generation
generatecommand(F,nele,ndof,thickness,gobalcoords0,E0,nu,edotMat,fixeddof,loaddof,Kin,Kout,outdof,xba);

%Plot components
[X,Y] = meshgrid(0:DW/nelx:DW,0:DH/nely:DH); 
contourf(X,Y,Phi_s,[0,0])
axis equal;axis([0 DW 0 DH]);

disp('Variable is converted to APDL command.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  generatecommand(F,nele,ndof,thickness,gobalcoords0,E0,nu,edotMat,fixeddof,loaddof,Kin,Kout,outdof,xba)    
fid = fopen('D:\APDL Working Space\command_test.txt','w');
ndot=ndof/2;
fprintf(fid,'/PREP7\n');
%define node
for i=1:ndot
    fprintf(fid,'N,%d,%G,%G,%G\n',i,gobalcoords0(2*i-1),gobalcoords0(2*i),0);
end
%define element type and the thickness
fprintf(fid,'et,1,plane182\nKEYOPT,1,1,0\nKEYOPT,1,3,3\nkeyopt,1,6,0\nR,1,%G, \nMPTEMP,1,0 \n',thickness);
% calculate the element parameters
MU0=E0*xba./(2*(1+nu));
K0=E0*xba./(2*(1-2*nu));
d=2./K0;
snele=0;
for i=1:nele
    if xba(i)>1e-6
    snele = snele+1;
    %define original elements
    fprintf(fid,'TB,HYPE,%d,1,2,NEO\nTBTEMP,0\n',snele);
    fprintf(fid,'TBDATA,,%G,%G,,,,\n',MU0(i),d(i));
    fprintf(fid,'type,1\nmat,%d\nreal,1\nesys,0\n',snele);
    fprintf(fid,'e,%d,%d,%d,%d\n',edotMat(i,1),edotMat(i,2),edotMat(i,3),edotMat(i,4));
    end
        
    
end
%%%%%define node of the spring
coords0w=reshape(gobalcoords0,2,[]);
lx=max(coords0w(1,:))-min(coords0w(1,:));
fprintf(fid,'N,%d,%f,%f,%f\n',ndot+1,gobalcoords0(loaddof)+2*lx,gobalcoords0(loaddof+1),0);
fprintf(fid,'d,%d,ux,0\n',ndot+1);
fprintf(fid,'d,%d,uy,0\n',ndot+1);
fprintf(fid,'N,%d,%f,%f,%f\n',ndot+2,gobalcoords0(outdof)+2*lx,gobalcoords0(outdof+1),0);
fprintf(fid,'d,%d,ux,0\n',ndot+2);
fprintf(fid,'d,%d,uy,0\n',ndot+2);
%%%%%define the spring
fprintf(fid,'ET,2,LINK180\nKEYOPT,2,2,0\nR,2,1, ,0\n');
fprintf(fid,'MPDATA,ex,%d,,%.15f\nmpdata,prxy,%d,,%f\ntype,2\nmat,%d\nreal,2\nesys,0\ne,%d,%d\n',2*nele+1,Kin*2*lx,2*nele+1,0.3,2*nele+1,(loaddof+1)/2,ndot+1);
fprintf(fid,'MPDATA,ex,%d,,%.15f\nmpdata,prxy,%d,,%f\ntype,2\nmat,%d\nreal,2\nesys,0\ne,%d,%d\n',2*nele+2,Kout*2*lx,2*nele+2,0.3,2*nele+2,(outdof+1)/2,ndot+2);
%apply the displacement
nfix=size(fixeddof,1);
for i=1:nfix
    if mod(fixeddof(i),2)==1
        fprintf(fid,'d,%d,ux,0\n',(fixeddof(i)+1)/2);
    else       
        fprintf(fid,'d,%d,uy,0\n',fixeddof(i)/2);
    end
end
%apply the external load
nload=size(loaddof,1);
for i=1:nload
    if mod(loaddof(i),2)==1
    fprintf(fid,'F,%d,fx,%G\n',(loaddof(i)+1)/2,full(F(loaddof(i))));
    else
    fprintf(fid,'F,%d,fy,%G\n',loaddof(i)/2,full(F(loaddof(i))));
    end
end
%solve 
fprintf(fid,'finish\n/sol\nANTYPE,0\nNLGEOM,1\nNSUBST,5,0,0\n');
fprintf(fid,'CNVTOL,U,-1, \n'); 
fprintf(fid,'CNVTOL,F,%G,0.001,2, ,  \n',sum(abs(full(F))));
fprintf(fid,'OUTRES,ERASE\nOUTRES,ALL,ALL\n/status,solu\nsolve\n');
%post processing---the nodal displacement 
fprintf(fid,'/POST1\n');
fprintf(fid,'SET, , ,1, ,1, , \n');
fprintf(fid,'nsel,all\n');
fprintf(fid,'*dim,UARRAY,array,%d,2\n',ndot);
fprintf(fid,'*vget,UARRAY(1,1),node,1,u,x\n');
fprintf(fid,'*vget,UARRAY(1,2),node,1,u,y\n');
fprintf(fid,'*cfopen,NODALDISPLACEMENT,txt\n');
fprintf(fid,'*vwrite,UARRAY(1,1),UARRAY(1,2)\n');
fprintf(fid,'(E13.5,5X,E13.5)\n');
fprintf(fid,'finish\n');
fclose(fid);
end