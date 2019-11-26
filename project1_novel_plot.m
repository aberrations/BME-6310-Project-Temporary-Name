function [yout,tout,params,F1,F2,F3] = project1(params);

% % % prompt1 = 'Do you want to include stress (Y/N) \n';
INPUT1 = 'Y';
% % % prompt2 = 'Do you want to include cis or trans; the default is cis-trans (C/T/B) \n';
INPUT2 = 'T';
% % % prompt3 = 'Do you want to include Delta-Knockdown (Y/N) \n';
betaD = 50;

Tmax=10; tspan=[0 Tmax]; % set time for simulation

params=defaultparams; % get the default parameters if none provided

% get the connectivity matrix
syms r
R = 6; % set number of rings
k = 1 + double(symsum(6*r, r, 1, R)); % number of cells
[HM] = hexmatrix(R);
params.connectivity=getconnectivityM(HM,k,R);

% setting the initial conditions + noise
y0=getIC(params,k,betaD);
HM = boundary(HM,R);
BASE = 5; %Base for exponential relationship (larger base means more exponential)
[S,HM] = stressfunc(HM,R,BASE);

for j=1:10
    betaD = j*5;
    % run simulation with lateral inhibition
    [tout,yout] = ode15s(@li,tspan,y0,[],params,S,INPUT1,INPUT2,betaD);
    % show time traces of two cells with lateral inhibition
    %new matrix of morphogenic impact
    ymorphD(:,j) = yout(1:130,1);
    ymorphR(:,j) = yout(1:130,k+1);
    ymorphN(:,j) = yout(1:130,2*k+1);
end
beta_mat = [5:5:50];
plot2cells(tout,ymorphD,ymorphR,ymorphN,k,beta_mat)

% show lattice simulation
% % % F1=movielattice(tout,yout,R,HM,k,INPUT1);
% % % F2=movielattice2(tout,yout,R,HM,k,INPUT1);
% % % F3=movielattice3(tout,yout,R,HM,k,INPUT1);

function dy = li(t,y,params,S,INPUT1,INPUT2,betaD) 

nu=params.nu;
%betaD=params.betaD;
betaN=params.betaN;
betaR=params.betaR;
m=params.m;
h=params.h;
M=params.connectivity;
k=length(M);
mu=params.mu;
kc=params.kc;
kt=params.kt;

D = y(1:k); % levels of Delta in cells 1 to k
R = y(k+1:2*k); % levels of Repressor in cells 1 to k
N = y(2*k+1:3*k); % levels of Repressor in cells 1 to k
Dneighbor=M*y(1:k);% average Delta level in the neighboring cells
Nneighbor=M*y(2*k+1:3*k); % average Notch level in the neighboring cells

if INPUT2 == 'T'
    kt = 10;
    kc = 0;
end
if INPUT2 == 'C'
    kt = 0;
    kc = 10;
end
if INPUT2 == 'B'
    kt = 10;
    kc = 10;
end

if INPUT1 == 'Y'
    dN = mu * (betaN - kt.*N.*Dneighbor-kc.*N.*D-N) + S/max(N); % differential equation describing Notch levels
    dD = nu * (betaD.*1./(1 + R.^h)-kt.*D.*Nneighbor-kc.*N.*D-D + S/max(D)); % differential equation describing Delta levels
    dR = betaR.*(N.*Dneighbor).^m./(1 + (N.*Dneighbor).^m)-R; % differential equation describing repressor levels
end

if INPUT1 == 'N'
    dN = mu * (betaN - kt.*N.*Dneighbor-kc.*N.*D-N); % differential equation describing Notch levels
    dD = nu * (betaD.*1./(1 + R.^h)-kt.*D.*Nneighbor-kc.*N.*D-D); % differential equation describing Delta levels
    dR = betaR.*(N.*Dneighbor).^m./(1 + (N.*Dneighbor).^m)-R; % differential equation describing repressor levels
end

dy = [dD;dR;dN]; 
function params=defaultparams

params.nu=1;
%params.betaD=50; %Delta Production (default value: 50)
params.betaN=1;
params.betaR=200;
params.m=1;
params.h=1;
params.sigma=0.2;
params.mu=1;
params.kc=10;
params.kt=1;
function [HM] = hexmatrix(R) %R is number of rings
%Creating array
HM=zeros(2*R+1,2*R+1); % connectivity matrix 
og=R+1; % identifying origin, so origin is (og,og)

%Identifying where cells exist in array as a non-zero value 
    for i=1:R
        ii = i-1;
        j = 2*(i-1); % even (0,2,4,...)
        l = 2*(i-1)+1; % odd (1,3,5,...)
        if (j<og) && (l<og) 
            HM([og-l:og-j og+j:og+l],[1+ii:og-1 og+1:(2*R+1)-ii]) = 1; % majority of cells
        end
        if (j<og)
            HM([og+j og-j],[og:og+(i-1) og-(i-1):og]) = 1; % middle even top/bottom cells
            HM([og-j og+j], og) = 1; % alternating cells at x=og
        end
        if (l<og)
            HM([og+l og-l],[(og+1):og+(i-1) og-(i-1):(og-1)]) = 1; % middle odd top/bottom cells
        end
        HM(og,og) = 1; % cell at origin 
    end
    seq = 1;
    for n=1:(2*R+1)
        for m=1:(2*R+1)
            if HM(n,m)==1
                HM(n,m)=seq;
                seq = seq + 1;
            end
        end
    end
function M = getconnectivityM(HM,k,R)
M = zeros(double(k),double(k));
w=1/6;
for i=1:(2*R+1)
    for j=1:(2*R+1)
        if HM(i,j)~=0   
            s = HM(i,j);
            out = [];
            if i>1
                out(1) = HM(i-1,j); %up
            end
            if j>1
                out(3) = HM(i,j-1); %left
            end
            if i>1 && j>1
                out(5) = HM(i-1,j-1); %up-left
            end
            if i<2*R+1
                out(2) = HM(i+1,j); %down
            end
            if j<2*R+1
                out(4) = HM(i,j+1); %right
            end
            if i<2*R+1 && j>1
                out(7) = HM(i+1,j-1); %down-left
            end
            if i>1 && j<2*R+1
                out(6) = HM(i-1,j+1); %up-right
            end
            if i<2*R+1 && j<2*R+1
            out(8) = HM(i+1,j+1); %down-right
            end
            out(out==0) = [];
            for r=1:length(out)
               M(s,out(r))=w;
            end
        end
    end
end 
function y0=getIC(params,k,betaD)

U=rand(k,1) - 1/2; % a uniform random distribution
epsilon=1e-5;   % multiplicative factor of Delta initial condition
D0=epsilon*betaD.*(1 + params.sigma*U); % initial Delta levels 
R0=zeros(k,1);  % initial repressor levels
N0=params.betaN.*ones(k,1);  % initial Notch levels are betaN
y0=[D0;R0;N0];  % vector of initial conditions

function plot2cells(tout,ymorphD,ymorphR,ymorphN,k,beta_mat);
% % % figure
%%%clf

% % % figure(1); surf(beta_mat,tout(1:130),ymorphD)
% % % title(['Delta Signaling']); ylabel('time [a.u]'); xlabel('Morphogen concentration [a.u]')
% % % figure(2); surf(beta_mat,tout(1:130),ymorphR)
% % % title(['Repressor Signaling']); ylabel('time [a.u]'); xlabel('Morphogen concentration [a.u]')
% % % figure(3); surf(beta_mat,tout(1:130),ymorphN)
% % % title(['Notch Signaling']); ylabel('time [a.u]'); xlabel('Morphogen concentration [a.u]')
% Fit new surface to given data
mnum = length(beta_mat); tnum = length(tout(1:130));
xm = []; xt = []; xD = []; xR = []; xN = [];
for i=1:mnum
    newxm = beta_mat(i)*ones(1,tnum); newxt = tout(1:130)';
    newxD = ymorphD(:,i)'; newxR = ymorphR(:,i)'; newxN = ymorphN(:,i)';
    xm = [xm, newxm]; xt = [xt, newxt];
    xD = [xD, newxD]; xR = [xR, newxR]; xN = [xN, newxN];
end
%'linearinterp'%'cubicinterp')
Dfit = fit([xm',xt'],xD','cubicinterp');
Rfit = fit([xm',xt'],xR','cubicinterp');
Nfit = fit([xm',xt'],xN','cubicinterp');
figure(1); plot(Dfit)
title(['Delta Signaling']); ylabel('time [a.u]'); xlabel('Morphogen concentration [a.u]')
figure(2); plot(Rfit)
title(['Repressor Signaling']); ylabel('time [a.u]'); xlabel('Morphogen concentration [a.u]')
figure(3); plot(Nfit)
title(['Notch Signaling']); ylabel('time [a.u]'); xlabel('Morphogen concentration [a.u]')

function plotHexagon(q0,p0,c,R)

% This function plots a hexagon centered at hex lattice coordinates p,q

s32 = sqrt(3)/4;

q = q0*3/4;
p = p0*2*s32;

if R/2 ~= round(R/2) %If odd
    if q0/2 == round(q0/2)
       if (p0<=R+1)
           p = p-s32;
       end
       if (p0>R+1)
           p = p-s32;
       end
    else
       if (p0>R+1)
           p = p-s32-s32;
       end
    end
end

if R/2 == round(R/2) %If even
    if q0/2 == round(q0/2)
       if (p0<=R+1)
           p = p+s32;
       end
       if (p0>R+1)
           p = p-s32;
       end
    else
       if (p0>R+1)
           p = p;
       end
    end
end

x(1) = p ; x(2) = p+s32; x(3) = p+s32; x(4) = p; x(5) = p-s32; x(6) = p-s32;

y(1) = q-.5; y(2) = q-.25; y(3) = q+.25; y(4) = q+.5; y(5) = q+.25; y(6) = q-.25;

c=min(c,ones(1,3));

patch(x,y,c,'linewidth',2);
function [x,y] = barrier(i,j,R,HM)
global m n v
m = 0; n = 0; v = 0; 

if i==1 || i==2*R+1
    y = HM(i,j,2);
    else 
        if i<R+1 % go up
            for m=1:i-1;
                if HM(i-m,j,2) ~= 0
                    y = HM(i-m,j,2);
                end
            end
            if HM(i,j,2)~=0 %if its a boundary cell
                y = HM(i,j,2);
            end
        end

        if i>R+1 % go down
            for m=1:(2*R+1 - i);
                if HM(i+m,j,2) ~= 0
                    y = HM(i+m,j,2);
                end
            end
            if HM(i,j,2)~=0 % if its a boundary cell
                y = HM(i,j,2);
            end
        end
        if i==R+1
            y = 0;
        end
end

if j==1 || j==2*R+1
    x = HM(i,j,3);
    else 
        if j>R+1 % go right
            for n=1:(2*R+1 - j);
                if HM(i,j+n,3) ~= 0
                    x = HM(i,j+n,3);
                end
            end
            if HM(i,j,3)~=0 % if its a boundary cell
                x = HM(i,j,3);
            end
        end

        if j<R+1 % go left
            for n=1:j-1; % the max distance to move left
                if HM(i,j-n,3) ~= 0 % if a cell is found to the left
                    x = HM(i,j-n,3); % make x equal the x position of the barrier to the left
                end
            end
            if HM(i,j,3)~=0 %if its a boundary cell
                x = HM(i,j,3);
            end
        end
        if j==R+1
            x = 0;
        end
end
function [S,HM] = stressfunc(HM,R,BASE)
NMX = numel(HM(:,:,1)); %Number of elements in matrix
HM(:,:,4) = 0;
for i = 1:(2*R+1)
    for j = 1:(2*R+1)
        if HM(i,j,1)~=0 
            [x,y] = barrier(i,j,R,HM); %Furthest x and y positions for the current cell at (i,j)
            x = abs((R+1) - x); %Distance from furthest x, for normalization
            y = abs((R+1) - y); %Distance from furthest y, for normalization
            if i==R+1
                HM(i,j,4) = (BASE^(abs((j-(R+1)))/x))/ (BASE);
            else
                if j==R+1
                    HM(i,j,4) = (BASE^(abs((i-(R+1)))/y))/ (BASE);
                else
                HM(i,j,4) = (BASE^(abs((i-(R+1)))/y) + BASE^(abs((j-(R+1)))/x))/ (2*BASE);
                end
            end
        end
    end
end

S = reshape(HM(:,:,4)',NMX,1);
S = S(S~=0);
function HM = boundary(HM,R)
%Create boundary identifier in HM 
HM(:,:,2) = 0; %Corresponds to row 
HM(:,:,3) = 0; %Corresponds to column
%Go Left-Right
for i=1:(2*R+1)
    t = 0;
    for j=1:(2*R+1)
        if HM(i,j,1)~=0
            t = t+1;
            if t == 1
                HM(i,j,2) = i;
                HM(i,j,3) = j;
            end
        end
    end
end

%Go Right-Left
for i=(2*R+1):-1:1
    t = 0;
    for j=(2*R+1):-1:1
        if HM(i,j,1)~=0
            t = t+1;
            if t == 1
                HM(i,j,2) = i;
                HM(i,j,3) = j;
            end
        end
    end
end
%Go Down-Up
for j=(2*R+1):-1:1
    t = 0;
    for i=(2*R+1):-1:1
        if HM(i,j,1)~=0
            t = t+1;
            if t == 1
                HM(i,j,2) = i;
                HM(i,j,3) = j;
            end
        end
    end
end

%Go Up-Down
for j=1:(2*R+1)
    t = 0;
    for i=1:(2*R+1)
        if HM(i,j,1)~=0
            t = t+1;
            if t == 1
                HM(i,j,2) = i;
                HM(i,j,3) = j;
            end
        end
    end
end
function F=movielattice(tout,yout,R,HM,k,INPUT1)
count = k;
figure('Name','Delta Levels')
sy1 = sort(yout(end,1:k));
Norm = sy1(round(length(sy1)*0.95)); % find the Delta level in the high Delta cells
frameind=0;
%HM(:,R+1) = 1;
for tind = 1:5:length(tout)   % shows every 5th frame
    clf;
    for i = 1:(2*R+1)
        for j = 1:(2*R+1)
            if HM(i,j,1)~=0 
                [x,y] = barrier(i,j,R,HM);
                ind = HM(i,j,1);
                if (Norm == 0)
                    Norm = 1;
                end
                mycolor = min([yout(tind,ind)/Norm,1]); % defining the normalized color of cell
              
                %The further away from R+1 (the center), the larger mycolor
                %should be (or the closer to 1 it should be)
                if INPUT1 == 'Y'
                    if i==R+1 && j==1 || j==R+1 && i==1 || i==R+1 && j==(2*R+1) || j==R+1 && i==(2*R+1) 
                    mycolor = mycolor * ( abs(((R+1)-i)/(R+1 - y)) + abs(((R+1)-j)/(R+1 - x)) );
                else
                    mycolor = mycolor * ( abs(((R+1)-i)/(R+1 - y)) + abs(((R+1)-j)/(R+1 - x)) ) / 2;
                    end
                end
                plotHexagon(i,j,[1-mycolor,1-mycolor,1],R);
            end
        end
    end
    axis image; axis on; box off; 
    
    frameind=frameind+1;
    F(frameind) = getframe; % generates a movie variable
end
%movie2avi(F,'movielattice'); % save movie in avi format

cellrad = []; 
for k = 1:(2*R+1)
    for l = 1:(2*R+1)
        if HM(k,l,1) ~= 0 
            cellrad(end+1) = HM(k,l,4); 
        end
    end
end

L = length(sy1);
deltaxls = yout(1:end,1:L);
reprexls = yout(1:end,L+1:2*L);
notchxls = yout(1:end,2*L+1:3*L);

if exist('Concentrations.xls', 'file') 
delete('Concentrations.xls'); 
end

xlswrite('Concentrations.xls',{'Delta'},'Sheet1','A1');
xlswrite('Concentrations.xls',{'Cell#:'},'Sheet1','A2');
xlswrite('Concentrations.xls',{'Radius:'},'Sheet1','A3');
xlswrite('Concentrations.xls',{'Time'},'Sheet1','A4');
xlswrite('Concentrations.xls',1:count,'Sheet1','B2');
xlswrite('Concentrations.xls',tout,'Sheet1','A5');
xlswrite('Concentrations.xls',deltaxls,'Sheet1','B5'); 
xlswrite('Concentrations.xls',cellrad,'Sheet1','B3');

xlswrite('Concentrations.xls',{'Repressor'},'Sheet2','A1');
xlswrite('Concentrations.xls',{'Cell#:'},'Sheet2','A2');
xlswrite('Concentrations.xls',{'Radius:'},'Sheet2','A3');
xlswrite('Concentrations.xls',{'Time'},'Sheet2','A4');
xlswrite('Concentrations.xls',1:count,'Sheet2','B2');
xlswrite('Concentrations.xls',tout,'Sheet2','A5');
xlswrite('Concentrations.xls',reprexls,'Sheet2','B5');     
xlswrite('Concentrations.xls',cellrad,'Sheet2','B3');

xlswrite('Concentrations.xls',{'Notch'},'Sheet3','A1');
xlswrite('Concentrations.xls',{'Cell#:'},'Sheet3','A2');
xlswrite('Concentrations.xls',{'Radius:'},'Sheet3','A3');
xlswrite('Concentrations.xls',{'Time'},'Sheet3','A4');
xlswrite('Concentrations.xls',1:count,'Sheet3','B2');
xlswrite('Concentrations.xls',tout,'Sheet3','A5');
xlswrite('Concentrations.xls',notchxls,'Sheet3','B5');
xlswrite('Concentrations.xls',cellrad,'Sheet3','B3');
function F=movielattice2(tout,yout,R,HM,k,INPUT1)
figure('Name','Notch Levels')
sy1 = sort(yout(end,(2*k+1):end));
Norm = sy1(round(length(sy1)*0.95)); % find the Delta level in the high Delta cells
frameind=0;
%HM(:,R+1) = 1;
for tind = 1:5:length(tout)   % shows every 5th frame
    clf;
    for i = 1:(2*R+1)
        for j = 1:(2*R+1)
            if HM(i,j,1)~=0 
                [x,y] = barrier(i,j,R,HM);
                ind = HM(i,j,1);
                if (Norm == 0)
                    Norm = 1;
                end
                mycolor = min([yout(tind,2*k+ind)/Norm,1]); % defining the normalized color of cell
              
                %The further away from R+1 (the center), the larger mycolor
                %should be (or the closer to 1 it should be)
                if INPUT1 == 'Y'
                    if i==R+1 && j==1 || j==R+1 && i==1 || i==R+1 && j==(2*R+1) || j==R+1 && i==(2*R+1) 
                    mycolor = mycolor * ( abs(((R+1)-i)/(R+1 - y)) + abs(((R+1)-j)/(R+1 - x)) );
                else
                    mycolor = mycolor * ( abs(((R+1)-i)/(R+1 - y)) + abs(((R+1)-j)/(R+1 - x)) ) / 2;
                    end
                end
                plotHexagon(i,j,[1,1-mycolor,1-mycolor],R);
            end
        end
    end
    axis image; axis on; box off; 
    
    frameind=frameind+1;
    F(frameind) = getframe; % generates a movie variable
end
function F=movielattice3(tout,yout,R,HM,k,INPUT1)
figure('Name','Notch Cleavage (Repressor) Levels')
sy1 = sort(yout(end,(k+1):k*2));
Norm = sy1(round(length(sy1)*0.95)); % find the Delta level in the high Delta cells
frameind=0;
%HM(:,R+1) = 1;
for tind = 1:5:length(tout)   % shows every 5th frame
    clf;
    for i = 1:(2*R+1)
        for j = 1:(2*R+1)
            if HM(i,j,1)~=0 
                [x,y] = barrier(i,j,R,HM);
                ind = HM(i,j,1);
                if (Norm == 0)
                    Norm = 1;
                end
                mycolor = min([yout(tind,k+ind)/Norm,1]); % defining the normalized color of cell
              
                %The further away from R+1 (the center), the larger mycolor
                %should be (or the closer to 1 it should be)
                if INPUT1 == 'Y'
                    if i==R+1 && j==1 || j==R+1 && i==1 || i==R+1 && j==(2*R+1) || j==R+1 && i==(2*R+1) 
                    mycolor = mycolor * ( abs(((R+1)-i)/(R+1 - y)) + abs(((R+1)-j)/(R+1 - x)) );
                else
                    mycolor = mycolor * ( abs(((R+1)-i)/(R+1 - y)) + abs(((R+1)-j)/(R+1 - x)) ) / 2;
                    end
                end
                plotHexagon(i,j,[1-mycolor,1,1-mycolor],R);
            end
        end
    end
    axis image; axis on; box off; 
    
    frameind=frameind+1;
    F(frameind) = getframe; % generates a movie variable
end
