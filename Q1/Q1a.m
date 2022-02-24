% This file implement the simulation for assignment 2 Q1 (a)
% V = V0 at x = 0 and V = 0 at x = L
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
deltaXY = 0.1;  % Assume deltaX = deltaY

% Declare the V0
V0 = 1;

% Calculate the dimension of solution matrix
nx = L/deltaXY;
ny = W/deltaXY;  

% GV = F
G = zeros(nx*ny, nx*ny);  % Declare the G matrix
F = zeros(nx*ny, 1);  % Declare the F matrix
matrixSol = zeros(nx, ny);  % Hold the solution matrix

% Constructor the F and G matrix
for ix = 1:nx  % Loop through the width
    for iy = 1:ny  % Loop through the height
        % Calculate the mapping index
        n = mappingEq(ix, iy, ny);
        % Check for boundary condition
        if ix == 1 || ix == nx
            % Setup F matrix
            if ix == 1
                F(n, 1) = V0;  % V = V0 at x = 0
            else
                F(n, 1) = 0;  % V = 0 at x = L
            end
            % Setup G matrix
            G(n, n) = 1;
        else
            % Setup F matrix
            F(n, 1) = 0;
            % Setup G matrix
            % Because dV/dy = 0 => (V(i+1,j) - 2V(i,j) + V(i-1,j))/deltaX^2 = 0
            G(n, n) = -2/deltaXY^2;
            % Calculate mapping index for V(i+1,j) and V(i-1,j)
            nxp = mappingEq(ix+1, iy, ny);  % index for V(i+1,j)
            nxm = mappingEq(ix-1, iy, ny);  % index for V(i-1,j)
            G(n, nxp) = 1/deltaXY^2;
            G(n, nxm) = 1/deltaXY^2;
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
    matrixSol(ix, iy) = V(iMap);
end

% Plot the solution
figure(2)
[X,Y] = meshgrid(linspace(0,L,nx), linspace(0,W,ny));
matrixSol = matrixSol';
surf(X,Y,matrixSol)
xlabel("X axis - Length")
ylabel("Y axis - Width")
zlabel("Z axis - Voltage")

% View on X-Z plane
view(0,0)



