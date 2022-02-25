% This file implement the simulation for assignment 2 Q1 (b)
% V = V0 at x = 0, x = L
% V = 0 at y = 0, y = W
% dV/dy = 0
% Assume deltaX = deltaY = 1

% Clear all
clearvars
clearvars -global
close all
format shorte

% Make plot pretier
set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultaxesfontsize',20)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'DefaultLineLineWidth', 2);

% Define the dimension
W = 1;  % Width
L = 3/2 * W;  % Height
deltaXY = 0.05;  % Assume deltaX = deltaY

% Declare the V0
V0 = 1;

% Simulation settings for solving analytical series
simSteps = 100; % Declare the simulation steps
pauseTime = 0.01; % Time paused per simulation step in second

% Calculate the dimension of solution matrix
nx = L/deltaXY;
ny = W/deltaXY;  

% GV = F
G = zeros(nx*ny, nx*ny);  % Declare the G matrix
F = zeros(nx*ny, 1);  % Declare the F matrix
% Hold the solution matrix for Finite Difference Method (FDM)
matrixSolFDM = zeros(nx, ny);  
% Hold the solution matrix for analytical series
matrixSolAnalytic = zeros(nx, ny);  

% Constructor the F and G matrix
for ix = 1:nx  % Loop through the width
    for iy = 1:ny  % Loop through the height
        % Calculate the mapping index
        n = mappingEq(ix, iy, ny);
        % Check for boundary condition
        if ix == 1 || ix == nx  % at x = 0 or x = L
            % Setup F matrix
            F(n, 1) = V0;  % V = V0 at x = 0, x = L
            % Setup G matrix
            G(n, n) = 1;
        elseif iy == 1 || iy == ny  % at y = 0 or y = W
            % setup F matrix
            F(n, 1) = 0;  % V = 0 at y = 0, y = W
            % Setup G matrix
            G(n, n) = 1;
        else
            % Setup F matrix
            F(n, 1) = 0;
            % Setup G matrix
            % (V(i+1,j) + V(i, j+1) - 4V(i,j) + V(i-1,j) + V(i, j-1))/delta^2 = 0
            G(n, n) = -4/deltaXY^2;
            % Calculate mapping index
            nxp = mappingEq(ix+1, iy, ny);  % index for V(i+1,j)
            nxm = mappingEq(ix-1, iy, ny);  % index for V(i-1,j)
            nyp = mappingEq(ix, iy+1, ny);  % index for V(i, j+1)
            nym = mappingEq(ix, iy-1, ny);  % index for V(i, j-1)
            % Set the G matrix
            G(n, nxp) = 1/deltaXY^2;
            G(n, nxm) = 1/deltaXY^2;
            G(n, nyp) = 1/deltaXY^2;
            G(n, nym) = 1/deltaXY^2;
        end
    end
end

% Plot the G matrix
figure(1)
spy(G);
title("G matrix")

% Solve for V from GV = F
V = G\F;

% Map back to the 2D region
for iMap = 1:nx*ny
    % Calculate the index for the 2D region
    ix = ceil(iMap/ny);
    iy = mod(iMap, ny);
    if iy == 0
        iy = ny;
    end
    % Assign the value
    matrixSolFDM(ix, iy) = V(iMap);
end

% Plot the solution from the Finite Difference Method
figure(2)
[X,Y] = meshgrid(linspace(0,L,nx), linspace(0,W,ny));
matrixSolFDM = matrixSolFDM';
surf(X,Y,matrixSolFDM)
xlabel("X axis - Length")
ylabel("Y axis - Width")
zlabel("Z axis - Voltage")
title("Finite Difference Method Solution")

% % Calculate using the analytical series
% figure(3)
% for iSim = 1:simSteps
%     % Traverse through the V function
%     
% 
%     % Pause some time
%     pause(pauseTime)
% end


% Calculate using the analytical series
figure(3)
% Calculate the a and b for the analytical series
a = W;
b = L/2;
matrixSolAnalytic = matrixSolAnalytic';
Xshifted = X-b;
% Loop for simulation for analytical solution
for iSim = 1:simSteps
    % Calculate the n in the sum: n = 1,3,5,7...
    n = iSim*2-1;
    % Calculate the sum
    matrixSolAnalytic = matrixSolAnalytic + 4*V0/pi* (1/n)*(cosh(n*pi*Xshifted/a)/cosh(n*pi*b/a)).*sin(n*pi*Y/a);

    % Plot the evolving surface
    figure(3)
    surf(X,Y,matrixSolAnalytic)
    xlabel("X axis - Length")
    ylabel("Y axis - Width")
    zlabel("Z axis - Voltage")
    title("Analytical Series Solution")

    % Pause some time
    pause(pauseTime)
end

% Calculate the difference between FDM and analytical solution
figure(4)
matrixSolDiff = matrixSolFDM - matrixSolAnalytic;
surf(X,Y, matrixSolDiff)
xlabel("X axis - Length")
ylabel("Y axis - width")
zlabel("Z axis - Voltage")
title("Solution Difference (FDM minus Analytical Series)")
