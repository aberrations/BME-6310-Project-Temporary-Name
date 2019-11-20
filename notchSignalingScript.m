%% Set parameters
nu = 1;
mu = 1;
betaN = 1; % maximal production rate of Notch
betaR = 200; % maximal production rate of Repressor
betaD = 50; % maximal production rate of Delta
m = 1; % cooperativity of repressor activation
h = 1; % cooperativity of Delta inhibition
sigma = 0.2; % for initial Delta levels
kt = 10; % constant for trans-activation (what we simulate here)
kc = 0; % constant for cis-activation (not simulating this here)

%% Create data structures
R = 6; % number of rings
[ hexmat, numCell ] = hexMatrix(R);
connmat = connectivityMatrix(hexmat, numCell, R);
k = length(connmat);

%% Overlay stress function on the hexagonal lattice
hexmat = createBoundaryLayers(hexmat, R);

% b = 0, 0.5, 5
[ stressmat0, hexmat0 ] = stress(hexmat, R, 0);
[ stressmat05, hexmat05 ] = stress(hexmat, R, 0.5);
[ stressmat5, hexmat5 ] = stress(hexmat, R, 5);

%% Initial conditions 
uniDist = rand(numCell, 1) - 1/2; % a uniform random distribution
epsilon = 1e-5;   % multiplicative factor of Delta initial condition
D0 = epsilon * betaD .* (1 + sigma*uniDist); % initial Delta levels 
R0 = zeros(numCell, 1);  % initial repressor levels
N0 = betaN .* ones(numCell, 1);  % initial Notch levels are betaN
y0 = [D0; R0; N0];  % vector of initial conditions

%% Simulations
tspan = [0 100];
[tout, yout0] = ode15s(@notchEq, tspan, y0, [], stressmat05, connmat);
%yout1 = interp1(tout, yout0, 0:1:100);
F = movielattice(tout, yout0, R, hexmat5, k, 0.5);

%% System of ODEs
function dy = notchEq(t, y, S, connmat)
    nu = 1;
    mu = 1;
    betaN = 1; % maximal production rate of Notch
    betaR = 200; % maximal production rate of Repressor
    betaD = 50; % maximal production rate of Delta
    m = 1; % cooperativity of repressor activation
    h = 1; % cooperativity of Delta inhibition
    sigma = 0.2; % for initial Delta levels
    kt = 10; % constant for trans-activation (what we simulate here)
    kc = 0; % constant for cis-activation (not simulating this here)
    k = length(connmat);
    
    D = y(1:k); % levels of Delta in cells 1 to k
    R = y(k+1:2*k); % levels of Repressor in cells 1 to k
    N = y(2*k+1:3*k); % levels of Repressor in cells 1 to k
    Dneighbor = connmat * y(1:k); % average Delta level in the neighboring cells
    Nneighbor = connmat * y(2*k+1:3*k); % average Notch level in the neighboring cells

    dN = mu * (betaN - kt.*N.*Dneighbor-kc.*N.*D-N) + S/max(N); % differential equation describing Notch levels
    dD = nu * (betaD.*1./(1 + R.^h)-kt.*D.*Nneighbor-kc.*N.*D-D + S/max(D)); % differential equation describing Delta levels
    dR = betaR.*(N.*Dneighbor).^m./(1 + (N.*Dneighbor).^m)-R; % differential equation describing repressor levels
    dy = [dD; dR; dN]; 
end

%% Functions for overlaying the stress function
% They created multiple dimensions in hexmat instead of creating different
% matrices to maintain the same variable, I guess. Second dimension in
% hexmat corresponds to the row of the boundary and third dimension of
% hexmat corresponds to the column of the boundary.
function [ hexmat ] = createBoundaryLayers(hexmat, R)
    hexmat(:, :, 2) = 0; % for the row
    hexmat(:, :, 3) = 0; % for the column
    for i=1:(2*R+1)
        for j=1:(2*R+1)
            if hexmat(i, j, 1) ~= 0 % part of the hex lattice
                if i == 1 || i == 2*R+1 % top or bottom boundary
                    hexmat(i, j, 2) = i;
                    hexmat(i, j, 3) = j;
                else
                    if j > 1 && j < 2*R+1
                        if (j ~= R && j ~= R+2) && ...
                           ((hexmat(i, j-1, 1) == 0 && hexmat(i, j+1, 1) ~= 0) || ...
                           (hexmat(i, j-1, 1) ~= 0 && hexmat(i, j+1, 1) == 0))
                            % left or right boundary, skipping the 0 in the
                            % middle (by defn of hexmat)
                            hexmat(i, j, 2) = i;
                            hexmat(i, j, 3) = j;
                        end
                    elseif j == 1 || j == 2*R+1
                        % middle row, end points
                        hexmat(i, j, 2) = i;
                        hexmat(i, j, 3) = j;
                    end
                end
            end
        end
    end                  
end

function [ stressmat, hexmat ] = stress(hexmat, R, b) 
    % b = measure of steepness of stress gradient
    numEle = numel(hexmat(:, :, 1)); % total num of elements
    hexmat(:, :, 4) = 0; % another layer for assigning stress
    
    for i = 1:(2*R+1)
        for j = 1:(2*R+1)
            if hexmat(i,j,1) ~= 0 
                [x,y] = farpoint(i, j, R, hexmat); 
                % Furthest x and y positions for the current cell at (i,j)
                % Normalize along radius (R+1, R+1) = origin
                x = abs((R+1) - x); %Distance from furthest x, for normalization
                y = abs((R+1) - y); %Distance from furthest y, for normalization
                if i == R+1
                    hexmat(i, j, 4) = (b^(abs((j - (R+1)))/x))/b;
                else
                    if j == R+1
                        hexmat(i, j, 4) = (b^(abs((i - (R+1)))/y))/b;
                    else
                        hexmat(i, j, 4) = (b^(abs((i-(R+1)))/y) + ...
                                            b^(abs((j-(R+1)))/x))/(2*b);
                    end
                end
            end
        end
    end

    stressmat = reshape(hexmat(:, :, 4)', numEle, 1);
    stressmat = stressmat(stressmat ~= 0);
    if b == 0
        stressmat = zeros(size(stressmat));
    end
end

function [ x, y ] = farpoint(i, j, R, hexmat) % find furthest point
    
    % Find y of far point
    if i == 1 || i == 2*R+1 % top or bottom boundary
        y = hexmat(i, j , 2); % col = row of this element (opp end)
    else
        if i < R+1 % top half traverse up till boundary
            for m=1:i-1
                if hexmat(i - m, j, 2) ~= 0 % go till boundary
                    y = hexmat(i - m, j, 2);
                end
            end
            if hexmat(i, j, 2) ~= 0 % boundary cell
                y = hexmat(i, j, 2); % opposite boundary
            end
        end
        
        if i > R+1 % bottom half, traverse down till boundary
            for m=1:(2*R+1 - i)
                if hexmat(i + m, j, 2) ~= 0
                    y = hexmat(i + m, j, 2);
                end
            end
            if hexmat(i, j, 2) ~= 0
                y = hexmat(i, j, 2);
            end
        end
        
        if i == R+1
            y = 0;
        end
    end
    
    % Find x of far point, same logic as y
    if j == 1 || j == 2*R+1
        x = hexmat(i, j, 3);
    else
        if j > R+1
            for n=1:(2*R+1 - j)
                if hexmat(i, j+n, 3) ~= 0
                    x = hexmat(i, j+n, 3);
                end
            end
            if hexmat(i, j, 3) ~= 0
                x = hexmat(i, j, 3);
            end
        end
        
        if j < R+1
            for n=1:j-1
                if hexmat(i, j-n, 3) ~= 0
                    x = hexmat(i, j-n, 3);
                end
            end
            if hexmat(i, j, 3) ~= 0
                x = hexmat(i, j, 3);
            end
        end
        if j == R+1
            x = 0;
        end
    end
end

%% Data structures needed for hexagonal layout

% Create a hexagonal matrix "array"
% In a 2D array, put non-zero values (they used sequential numbers) to show
% where the hexagonal matrix values are.
function [ hexmat, numCell ] = hexMatrix(R)
    % R is the number of rings. At the widest part of the hex matrix, the
    % number of elements are R+R+origin.
    hexmat = ones(2*R+1, 2*R+1);
    origin = hexmat(R+1, R+1); % center element
    
    % In each odd row, they have all numbers. In each even row, the middle
    % element is zero. (to account for that one less number and maintain
    % zeros at the end)
    % Populating from 2R+1th row
    for i=1:R
        if mod(i,2) == 1 % even multiple of R
            numEle = (2*R+1) - i;
            numZero = (2*R+1) - numEle; % number of zeros in the mid+sides
            hexmat(R+1-i, R+1) = 0; numZero = numZero - 1; % zero in the mid
            for j=1:(numZero/2) % populate zeros from both sides
                hexmat(R+1-i, j) = 0;
                hexmat(R+1-i, end-j+1) = 0;
            end
        else % odd multiple of R
            numEle = (2*R+1) - i;
            numZero = (2*R+1) - numEle; % number of zeros in the mid+sides
            for j=1:(numZero/2) % populate zeros from both sides
                hexmat(R+1-i, j) = 0;
                hexmat(R+1-i, end-j+1) = 0;
            end
        end
    end
    hexmat(R+2:end, :) = hexmat(R:-1:1, :); %since they are symmetric
    
    % assign numbers to the elements
    numCell = 1;
    for i=1:2*R+1
        for j=1:2*R+1
            if hexmat(i, j) == 1
                hexmat(i, j) = numCell;
                numCell = numCell+1;
            end
        end
    end
    numCell = numCell-1;
end

% Create a connectivity matrix with size = number of cells in hexmat
% For each connection, they assigned a weightage? 1/6
function [ connmat ] = connectivityMatrix(hexmat, numCell, R)
    connmat = zeros(numCell, numCell);
    weight = 1/6;
    
    % Traverse through each cell in hexmat and look at surrounding cells to
    % identify neighboring cells
    for i=1:2*R+1
        for j=1:2*R+1
            cellNum = hexmat(i, j);
            conn = [];
            if cellNum ~= 0
                if i > 1
                    conn(1) = hexmat(i-1, j); % up
                end
                if i < 2*R+1
                    conn(2) = hexmat(i+1, j); % down
                end
                if j > 1
                    conn(3) = hexmat(i, j-1); % left
                end
                if j < 2*R+1
                    conn(4) = hexmat(i, j+1); % right
                end
                if i > 1 && j > 1
                    conn(5) = hexmat(i-1, j-1); % up-left
                end
                if i < 2*R+1 && j > 1
                    conn(6) = hexmat(i+1, j-1); % down-left
                end
                if i > 1 && j < 2*R+1
                    conn(7) = hexmat(i-1, j+1); % up-right
                end
                if i < 2*R+1 && j < 2*R+1
                    conn(8) = hexmat(i+1, j+1); % down-right
                end 
            end
            for x=1:length(conn)
                if conn(x) ~= 0
                    connmat(cellNum, conn(x)) = weight;
                end
            end
        end
    end
end

%% Images
function v=movielattice(tout, yout, R, hexmat, k, b)
    v = VideoWriter('movielattice5.avi'); open(v);
    figure('Name', 'Delta Levels');
    sy1 = sort(yout(end,1:k));
    Norm = sy1(round(length(sy1)*0.95)); % find the Delta level in the high Delta cells
    for tind = 1:5:length(tout)   % shows every 5th frame
        clf;
        for i = 1:(2*R+1)
            for j = 1:(2*R+1)
                if hexmat(i, j, 1) ~= 0 
                    [x,y] = farpoint(i,j,R,hexmat);
                    ind = hexmat(i, j, 1);
                    if (Norm == 0)
                        Norm = 1;
                    end
                    mycolor = min([yout(tind,ind)/Norm,1]); % defining the normalized color of cell

                    %The further away from R+1 (the center), the larger mycolor
                    %should be (or the closer to 1 it should be)
                    if b ~= 0
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
        hold on
        drawnow;
        writeVideo (v, getframe(gcf)); 
        hold off
    end
    close(v);
end

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
end