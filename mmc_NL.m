% MMC-Based Method (B-Spline variable width) (test ver)
% Nonlinear FEM by APDL
% Case: Displacement Inverter
clear; clc; close all; 
%% Initialize - Material Parameters
E0 = 1e8;         % Young's Modulus
nu = 0.3;         % Poisson Ratio
c2=5e-4;
strain0=0.6;

%% Initialize - Boundary Condition and Mesh
volfrac=0.2;
[nele,nelx,nely,DW,DH,EW,EH,thickness,gobalcoords0,edotMat,fixeddof,loaddof,F,ndof,Kin,Kout,outdof]=initialize2Dinverter;
[x,y] = meshgrid(0:EW:DW,0:EH:DH);
DR = 4e-3;  % maximum spline half width

%% Initialize - Component Geometry
N = 12;                 % number of spline
Ncp= 8;                 % number of control points
N2 = 4;                 % number of spline of each group   
Var_num = (Ncp-2)*2+Ncp; % number of design variables for each component
Phi=cell(N,1);          % Topology Description Function
Bx = zeros(N,Ncp);
By = zeros(N,Ncp);
Br = zeros(N,Ncp);
for i = 1:N
    s_index = rem(i,N2);
switch i
    case {1,2,3,4}
        Bx(i,:) = 0:DW/(Ncp-1):DW;
        if s_index == 0
            By(i,:) = [0,0.15,0.05,0.02,0.12,0.07,0.2,1]*DH;
        elseif s_index == 1
            By(i,:) = [0,0.35,0.25,0.22,0.3,0.23,0.4,1]*DH;
        elseif s_index == 2
            By(i,:) = [0,0.55,0.65,0.57,0.7,0.53,0.65,1]*DH;
        else
            By(i,:) = [0,0.75,0.85,0.77,0.9,0.83,0.92,1]*DH;
        end
    case {5,6,7,8}
        Bx(i,:) = 0:DW/(Ncp-1):DW;
        if s_index == 0
            By(i,:) = [1,0.9,0.8,0.82,0.95,0.91,0.8,1]*DH;
        elseif s_index == 1
            By(i,:) = [1,0.7,0.6,0.62,0.75,0.71,0.6,1]*DH;
        elseif s_index == 2
            By(i,:) = [1,0.4,0.3,0.32,0.45,0.41,0.3,1]*DH;
        else
            By(i,:) = [1,0.13,0.02,0.05,0.15,0.11,0.13,1]*DH;
        end
    case {9,10,11,12}
        By(i,:) = 0:DH/(Ncp-1):DH;
        if s_index == 0
            Bx(i,:) = [0,0.83,0.93,0.75,0.97,0.92,0.8,0]*DW;
        elseif s_index == 1
            Bx(i,:) = [0,0.63,0.73,0.55,0.77,0.72,0.6,0]*DW;
        elseif s_index == 2
            Bx(i,:) = [0,0.35,0.3,0.4,0.45,0.42,0.2,0]*DW;
        else
            Bx(i,:) = [0,0.07,0.12,0.04,0.13,0.17,0.1,0]*DW;
        end
end
   Br(i,:) = 1.3e-3*ones(1,Ncp);
end
variable = zeros(Var_num,N);
c_ends = zeros(N,4);
for i = 1:N
    variable(:,i)=[Bx(i,2:end-1)';By(i,2:end-1)';Br(i,1:end)'];
    c_ends(i,1) = Bx(i,1);
    c_ends(i,2) = By(i,1);
    c_ends(i,3) = Bx(i,end);
    c_ends(i,4) = By(i,end);
end % order of var: [spline1(Bx2...Bxn By2...Byn-1 Br1...Brn);spline2( ... )... ;splinen( ... )]

%% Initialize - MMA Parameters  
changeobj = 1;
xy00 = variable(:);
xval = xy00;
xold1 = xy00;
xold2 = xy00;
xmin = [0*ones(Ncp-2,1);0*ones(Ncp-2,1);0.45e-3*ones(Ncp,1)];
xmin = repmat(xmin,N,1);  % minimum of variable:[ BX ; By ; Br ];
xmax = [DW*ones(Ncp-2,1);DH*ones(Ncp-2,1);DR*ones(Ncp,1)];
xmax = repmat(xmax,N,1);  % maximum of variable:[ BX ; By ; Br ];
low = xmin;
upp = xmax;
m = 1;  %number of constraint
nn=Var_num*N;
c=10000*ones(m,1);
d=zeros(m,1);
a0=1;
a=zeros(m,1);
df0dxba = zeros(nele,1);
dxbadx = zeros(nele,nn);
df0dx = zeros(nn,1);
dfdx = zeros(nn,1);

%% Initialize - Optimization Parameters
loop = 0;
maxiter = 1500; 
obj = zeros(1,maxiter);
looplist = linspace(1,maxiter,maxiter);      
alpha=1e-3;              % parameter alpha in the Heaviside function
epsilon=0.1*min(EW,EH);  % regularization parameter epsilon in the Heaviside function
penal=1;                 % penalty for element density
tt=zeros(4,100);
xy_opt = xy00; obj_opt = 100;


%%　Main Loop
while changeobj>1e-4 && loop < maxiter
    tic
    loop = loop+1;
    lp = num2str(loop);

    % (1) Forming Phi^s
    for i = 1:N
        Phi{i} = TDF(xy00(Var_num*i-Var_num+1:Var_num*i),Ncp,c_ends(i,:),DH,DW,nelx,nely);
    end
    Phi_s = Phi{1};
    for i = 2:N
        Phi_s = max(Phi_s,Phi{i});
    end

    % (2) Heaviside Function
    H = Heaviside(Phi_s,alpha,nelx,nely,epsilon);
    
    % (3) Element Density
    temp1 = zeros(nele,1); 
    for i = 1:nele
        temp1(i) = ( H(edotMat(i,1))^2+H(edotMat(i,2))^2+H(edotMat(i,3))^2+H(edotMat(i,4))^2 )/4;  % ersatz material model
    end
    xba = reshape(temp1,nely,nelx);
    xba = flip(xba);
    xba = xba(:);

    % (4) Generate APDL Command 
    generatecommand(F,nele,ndof,thickness,gobalcoords0,E0,nu,xba,penal,edotMat,fixeddof,loaddof,c2,Kin,Kout,outdof);
    tt(1,loop)=toc;

    % (5) Run ANSYS APDL
    [sta, cmd] = dos('apdl.bat','-echo');
    tt(2,loop)=toc;

    % (6) Read APDL Output
    % (i) Objective value
    U=load('NODALDISPLACEMENT.txt');
    delete('NODALDISPLACEMENT.txt');
    strain_energy=load('strainenergy.txt');
    strain_energy99=load('strainenergy99.txt');
    U=U';U=U(:);
    obj(loop)=U(outdof);
    % (ii) Element logarithmic equivilant strain to update c2
    strain=load('eqvstrain.txt');
    strain=reshape(strain,2,[]);
    maxvms=max(strain');svon=maxvms(2);
    if svon<=strain0
        c2=max((svon/strain0)^0.5,0.8)*c2;
    else
        c2=(svon/strain0)^2*c2;        
    end
    c2=max(c2,1e-6);

    % (7) Sensitivity Analysis
    % (i) df0dxba: nele*1
    dstrain_energy=reshape(strain_energy99-strain_energy,2,[]);     
    for i=1:nele
        df0dxba(i)=-1/(F(loaddof)/100)*(dstrain_energy(1,i)*penal/xba(i)-...
            dstrain_energy(2,i)*penal*xba(i)^(penal-1)/(1-xba(i)^penal+1e-10));
    end
    df0dxba = reshape(df0dxba,nely,nelx);
    df0dxba = flip(df0dxba);
    df0dxba = df0dxba(:);

    % (ii) dHdx: nnode*1
    diffH=cell(N,1);
    deltaxy=2*min(EW,EH);
    deltar=DR/10;
    for j=1:N
        for ii=1:Var_num
            xy001=xy00;
            if ii>(Ncp-2)*2
                delta = deltar; 
            else
                delta = deltaxy;
            end
            xy001(ii+(j-1)*Var_num)=xy00(ii+(j-1)*Var_num)+delta;
            PhiD1=TDF(xy001(Var_num*j-Var_num+1:Var_num*j),Ncp,c_ends(j,:),DH,DW,nelx,nely);
            Phi_s1=PhiD1;
            for ik=1:j-1
                Phi_s1=max(Phi_s1,Phi{ik});
            end
            for ik=j+1:N
                Phi_s1=max(Phi_s1,Phi{ik});
            end
            xy002=xy00;
            xy002(ii+(j-1)*Var_num)=xy00(ii+(j-1)*Var_num)-delta;
            PhiD2=TDF(xy002(Var_num*j-Var_num+1:Var_num*j),Ncp,c_ends(j,:),DH,DW,nelx,nely);
            Phi_s2=PhiD2;
            for ik=1:j-1
                Phi_s2=max(Phi_s2,Phi{ik});
            end
            for ik=j+1:N
                Phi_s2=max(Phi_s2,Phi{ik});
            end
            HD1=Heaviside(Phi_s1,alpha,nelx,nely,epsilon);
            HD2=Heaviside(Phi_s2,alpha,nelx,nely,epsilon);
            diffH{j}(:,ii)=(HD1-HD2)/(2*delta);
        end
    end
    % (iii) df0dx
    temp2 = zeros((nelx+1)*(nely+1),nn);
    for j=1:N
        for ii=1:Var_num
            temp2(:,ii+(j-1)*Var_num) =(1/2)* H'.*diffH{j}(:,ii);
        end
    end
    for j = 1:nele
            dxbadx(j,:) = temp2(edotMat(j,1),:)+temp2(edotMat(j,2),:)+temp2(edotMat(j,3),:)+temp2(edotMat(j,4),:);
    end
    df0dx = (df0dxba'*dxbadx)';

    % (iv) dfdx
    for j = 1:nn
        k = 1+floor((j-1)/Var_num) ; s = 1+rem((j-1),Var_num) ;
        dHdx = diffH{k}(:,s) ;
        dden = sum(dHdx(edotMat),2)/4;
        dfdx(j) = sum(dden)*(EW*EH);
    end
    tt(3,loop)=toc;

    % (8) MMA
    den=sum( H(edotMat), 2 ) / 4;
    f0val = obj(loop);
    fval = (sum(den)*(EW*EH)-volfrac*(DW*DH)) ;

    % scaling factors
    f0_scale = 10000;
    f_scale = 10000;
    var_scale = 10000;
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,ss,low,upp] = ...
        mmasub(m,nn,loop,var_scale*xval,var_scale*xmin,var_scale*xmax,var_scale*xold1,var_scale*xold2, ...
        f0_scale*f0val+95,f0_scale*df0dx,f_scale*fval,f_scale*dfdx,var_scale*low,var_scale*upp,a0,a,c,d); 
    xold2 = xold1;
    xold1 = xval;
    low = low/var_scale;
    upp = upp/var_scale;
    changeobj=max(abs(xval-xmma/var_scale));
    xval = xmma/var_scale;
    xy00 = xval;

    % (9) loop information
    disp([' It.: ' sprintf('%4i\t',loop) ' Obj.: ' sprintf('%f\t',obj(loop)) ' Vol.: ' ...
        sprintf('%6.4f\t',(fval+volfrac*(DW*DH))/(DW*DH)) 'ch.:' sprintf('%6.4f\t',changeobj)]);

    % (10) Plots
    figure(1)
    contourf(x, y,Phi_s,[0,0]);
    axis equal;axis([0 DW 0 DH]); 
    txt = ['iter ',lp];
    title(txt)
%     saveas(gcf,lp,'jpeg')
    figure (2)
    plot(looplist(1:loop),obj(1:loop))
    xlabel('loop'), ylabel('obj')
    pause(1e-6)

    % (11) Solution optimal
    if fval<=0 && loop>20 && obj(loop)<obj_opt
        xy_opt=xy00; obj_opt=obj(loop) ;
    end
    writematrix(xy_opt,'var_opt.txt')
    tt(4,loop)=toc;

end


%% Fuction
function [nele,nelx,nely,DW,DH,EW,EH,thickness,gobalcoords0,edotMat,fixeddof,loaddof,F,ndof,Kin,Kout,outdof]=initialize2Dinverter
DW = 100e-3; DH = 50e-3;
nelx = 100; nely = 50;
EW = DW/nelx; EH = DH/nely; thickness = 1e-3;
[i0,j0] = meshgrid(0:nelx,0:nely);
gobalcoords0x = i0*EW;
gobalcoords0y = EH*nely-j0*EH;
gobalcoords0xy = [gobalcoords0x(:) gobalcoords0y(:)]';
gobalcoords0 = gobalcoords0xy(:);
[il,jl] = meshgrid(0,nely);
loadnid = il*(nely+1)+(nely+1-jl); 
loaddof = 2*loadnid(:)-1; 
force = 5;
Kin = 500;
[ilo,jlo] = meshgrid(nelx,nely);
outnid = ilo*(nely+1)+(nely+1-jlo);  
outdof = 2*outnid(:)-1; 
Kout = 100;
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
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*(nely+1) 2*(nely+1)+1 2*(nely+1)-2 2*(nely+1)-1 -2 -1],nele,1);%%%%b=[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],超神节点提取符写法，阐述其相对单元的1节点x方向自由度的位置
edotMat=0.5*edofMat(:,[2,4,6,8]);
end

function  generatecommand(F,nele,ndof,thickness,gobalcoords0,E0,nu,xba,penal,edotMat,fixeddof,loaddof,c2rate,Kin,Kout,outdof)   
fid = fopen('command.txt','w');
ndot=ndof/2;
fprintf(fid,'/PREP7\n');
% define node
for i=1:ndot
fprintf(fid,'N,%d,%G,%G,%G\n',i,gobalcoords0(2*i-1),gobalcoords0(2*i),0);
end
% define element type and the thickness
fprintf(fid,'et,1,plane182\nKEYOPT,1,1,0\nKEYOPT,1,3,3\nkeyopt,1,6,0\nR,1,%G, \nMPTEMP,1,0 \n',thickness);
% calculate the element parameters
MU0=E0*xba.^penal/(2*(1+nu));
K0=E0*xba.^penal/(2*(1-2*nu));
d=2./K0;
c1=(1-xba.^penal)*E0*1e-9/6;
c2=(1-xba.^penal)*E0*c2rate;

for i=1:nele   
        %define original elements
        fprintf(fid,'TB,HYPE,%d,1,2,NEO\nTBTEMP,0\n',2*i-1);
        fprintf(fid,'TBDATA,,%G,%G,,,,\n',MU0(i),d(i));
        fprintf(fid,'type,1\nmat,%d\nreal,1\nesys,0\n',2*i-1);
        fprintf(fid,'e,%d,%d,%d,%d\n',edotMat(i,1),edotMat(i,2),edotMat(i,3),edotMat(i,4));

        %define the additive hyperelastic elements  
        fprintf(fid,'TB,HYPE,%d,1,2,YEOH\nTBTEMP,0\n',2*i);
        fprintf(fid,'TBDATA,,%G,%G,1e-9,1e6,,\n',c1(i),c2(i));
        fprintf(fid,'type,1\nmat,%d\nreal,1\nesys,0\n',2*i);
        fprintf(fid,'e,%d,%d,%d,%d\n',edotMat(i,1),edotMat(i,2),edotMat(i,3),edotMat(i,4));
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
%apply the additional external force
fprintf(fid,'F,%d,fx,%f\n',(outdof+1)/2,full(F(loaddof)/100));
fprintf(fid,'NSUBST,5,0,0\n');
%relax convergence criteria
fprintf(fid,'CNVTOL,F,%G,0.001,2, ,  \n',1e4*sum(abs(full(F))));
fprintf(fid,'solve\n');
%dummy load case
fprintf(fid,'F,%d,fx,%f\n',loaddof,0);
fprintf(fid,'NSUBST,5,0,0\n');
fprintf(fid,'CNVTOL,F,%G,0.001,2, ,  \n',1e4*sum(abs(full(F))));
fprintf(fid,'solve\nfinish\n');
%post processing---write the element strain energy 
fprintf(fid,'/POST1\n');
fprintf(fid,'SET, , ,1, ,1, , \n');
fprintf(fid,'esel,all\n');
fprintf(fid,'*dim,STEN,array,%d,1\n',2*nele);
fprintf(fid,'ETABLE,ea,SENE\n');
fprintf(fid,'*vget,STEN,elem,1,etab,ea\n');
fprintf(fid,'*cfopen,strainenergy,txt\n');
fprintf(fid,'*vwrite,STEN(1,1)\n');
fprintf(fid,'(E13.5)\n');
%post processing---the element logarithmic equivilant strain
fprintf(fid,'*dim,strain1,array,%d,1\n',2*nele);
fprintf(fid,'ETABLE,ab,EPEL,EQV\n');
fprintf(fid,'*vget,strain1,elem,1,etab,ab\n');
fprintf(fid,'*cfopen,eqvstrain,txt\n');
fprintf(fid,'*vwrite,strain1(1,1)\n');
fprintf(fid,'(E13.5)\n');
%post processing---the nodal displacement
fprintf(fid,'nsel,all\n');
fprintf(fid,'*dim,UARRAY,array,%d,2\n',ndot);
fprintf(fid,'*vget,UARRAY(1,1),node,1,u,x\n');
fprintf(fid,'*vget,UARRAY(1,2),node,1,u,y\n');
fprintf(fid,'*cfopen,NODALDISPLACEMENT,txt\n');
fprintf(fid,'*vwrite,UARRAY(1,1),UARRAY(1,2)\n');
fprintf(fid,'(E13.5,5X,E13.5)\n');
%post processing---write the element strain energy 
fprintf(fid,'SET, , ,1, ,2, , \n');
fprintf(fid,'*dim,STEN99,array,%d,1\n',2*nele);
fprintf(fid,'ETABLE,ea99,SENE\n');
fprintf(fid,'*vget,STEN99,elem,1,etab,ea99\n');
fprintf(fid,'*cfopen,strainenergy99,txt\n');
fprintf(fid,'*vwrite,STEN99(1,1)\n');
fprintf(fid,'(E13.5)\n');
% dummy load data
fprintf(fid,'SET, , ,1, ,3, , \n');
fprintf(fid,'*dim,STEN33,array,%d,1\n',2*nele);
fprintf(fid,'ETABLE,ea33,SENE\n');
fprintf(fid,'*vget,STEN33,elem,1,etab,ea33\n');
fprintf(fid,'*cfopen,strainenergy33,txt\n');
fprintf(fid,'*vwrite,STEN33(1,1)\n');
fprintf(fid,'(E13.5)\n');
fprintf(fid,'finish\n');
fclose(fid);
end