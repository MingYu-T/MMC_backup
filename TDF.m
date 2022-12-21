%% Forming Topology Description Function (Phi_s for whole structure)
%% B-spline with variable width
function Phi = TDF(var,Ncp,c_ends,DH,DW,nelx,nely)
% - Arrange variables
    ord=3;
    Bx = zeros(1,Ncp);
    By = zeros(1,Ncp);

    Bx(2:end-1) = var( 1:Ncp-2 );
    By(2:end-1) = var( (1+Ncp-2):(Ncp-2)*2 );
    Br = var(1+(Ncp-2)*2:end)';
    Bx(1) = c_ends(1);
    By(1) = c_ends(2);
    Bx(end) = c_ends(3);
    By(end) = c_ends(4);

% eliminate small components
    for i = 1:numel(Br)
        if Br(i)==0.45e-3
            Br(i) = 0;
        end
    end
    
% - Parameters of Splines
    n = Ncp-1;                    % n+1 is the number of control points
    k = ord;                              % order of the polynomial, i.e. k=degree+1
    [t,Range] = KnotVector(k,n); % knot vector
    SampleSize = 500;
    % Sample Size is the number of divisions between the knot values
    

% - Parameters of Level Set Function 
    [X,Y] = meshgrid(0:DW/nelx:DW,0:DH/nely:DH);
    [ysize,xsize] = size(X);
    Phi = zeros(ysize,xsize);

% - Base Function
    x = zeros(1,(n+k)*SampleSize+1);
    x(((k-1)*SampleSize+1):((n+1)*SampleSize+1)) = 0:1/SampleSize:Range;
    x(((n+1)*SampleSize+2):end) = Range;
    x = x./Range; % normalize

    N1 = FirstOrderBSplineFunction(k,t,x,SampleSize);
    N = GeneralOrderBsplineFunction(k,t,x,SampleSize,N1);

% - B-spline points
    Px = Bx*N(1:n+1,(k-1)*SampleSize+1:(n+1)*SampleSize);
    Py = By*N(1:n+1,(k-1)*SampleSize+1:(n+1)*SampleSize);
    r  = Br*N(1:n+1,(k-1)*SampleSize+1:(n+1)*SampleSize);

% - Level Set Function
    for i = 1:ysize
        for j = 1:xsize
            dmin = 100;
            rn = 0;
            for k = 1:SampleSize*(n-1)
                d = sqrt((Px(k)-X(i,j))^2+(Py(k)-Y(i,j))^2);
                if d < dmin
                    dmin = d;
                    rn = k;
                end
            end
            Phi(i,j) = r(rn)-dmin;   % Phi of the splines
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Knot Vector
function [t,Range] = KnotVector(k,n)
Range = n-k+2;
for i = k+1:n+1
    t(i) = i-k;
end
t(n+2:n+k+1) = Range;
t = t./Range;   % normalize
end

% First Order Base Function
function N1 = FirstOrderBSplineFunction(k,t,x,SampleSize)
m1 = length(t); m2 = length(x);
N1 = zeros(m1-1,m2);
k1 = k;
k2 = m1-k;
for i = k1:k2
    j1 = (i-1)*SampleSize+1;
    j2 = j1+SampleSize-1;
    N1(i,j1:j2) = 1;
end
end

% General Order Base Function
function N = GeneralOrderBsplineFunction(k,t,x,SampleSize,N1)
m1 = length(t); 
N = N1;
for l = 2:k
    d = l;
    k1 = k-(d-1);
    k2 = m1-k;

    i = k1;
    j1 = (i-1)*SampleSize+1;
    j2 = j1+d*SampleSize-1;
    N(i,j1:j2) = (t(i+d)-x(j1:j2))./(t(i+d)-t(i+1)).*N(i+1,j1:j2);
    for i = k1+1:k2-1
        j1 = (i-1)*SampleSize+1;
        j2 = j1+d*SampleSize-1;
        N(i,j1:j2) = (x(j1:j2)-t(i))./(t(i+d-1)-t(i)).*N(i,j1:j2)+(t(i+d)-x(j1:j2))./(t(i+d)-t(i+1)).*N(i+1,j1:j2);
    end
    i = k2;
    j1 = (i-1)*SampleSize+1;
    j2 = j1+d*SampleSize-1;
    N(i,j1:j2) = (x(j1:j2)-t(i))./(t(i+d-1)-t(i)).*N(i,j1:j2);
end
end