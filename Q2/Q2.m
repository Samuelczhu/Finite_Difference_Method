% This file implement the simulation for assignment 2 Q2 
% using divergence of sigma times gradient of V equals 0
% Assume boundary condition:
% V = V0 at x = 0 and V = 0 at x = L
% dV/dy = 0 at y = 0, y = W
% Assume deltaX = deltaY

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
L = 1.5;  % Height
Wb = 0.4*W;  % Width of the box
Lb = 0.3*L;  % Length of the box
deltaXY = 0.02;  % Assume deltaX = deltaY

% Declare the V0
V0 = 1;

% Calculate the dimension of solution matrix
nx = L/deltaXY;
ny = W/deltaXY;  
[X,Y] = meshgrid(linspace(0,L,nx), linspace(0,W,ny));

% Declare the matrix for conductivity: Sigma(y,x)
matrixSigma = ones(ny, nx);  % Dimension: ny times nx
xIndexBox = ceil((L-Lb)/(2*deltaXY));  % Find the starting x index for the box
LbIndexRange = ceil(Lb/deltaXY);  % Index range for the length of the box
WbIndexRange = ceil(Wb/deltaXY);  % Index range for the width of the box
% Assign the region for the box
matrixSigma(1:WbIndexRange, xIndexBox:xIndexBox+LbIndexRange) = 10^-2;
matrixSigma(ny-WbIndexRange:ny, xIndexBox:xIndexBox+LbIndexRange) = 10^-2;

% Plot the region for conductivity
figure(1)
surf(X,Y, matrixSigma)
title("Plot of Conductivity Sigma(x,y)")
xlabel("X axis - Length")
ylabel("Y axis - Width")
zlabel("Z axis - Conductivity")
view(0,90)  % View from top

% Declare the matrix for voltage V(y,x)
matrixV = zeros(ny, nx);  % Dimension: ny times nx

% Declare the G matrix and F vector: GV = F
G = zeros(nx*ny, nx*ny);  
F = zeros(nx*ny, 1);

% Construct the G matrix and F vector
for ix = 1:nx
    for iy = 1:ny
        % Calculate the index
        n = mappingEq(ix, iy, ny);
        % Check for the boundary
        if ix==1 || ix==nx || iy ==1 || iy==ny
            G(n,n) = 1;
            % Boundary condition for x
            if ix == 1
                F(n,1) = V0;  % V = V0 at x = 0 
            elseif ix == nx
                F(n,1) = 0;  % and V = 0 at x = L
            elseif iy == 1
                nyp = mappingEq(ix, iy+1, ny);  % dV/dy=0 at y=0
                G(n,nyp) = -1;
            elseif iy == ny
                nym = mappingEq(ix, iy-1, ny);  % dV/dy=0 at y=W
                G(n, nym) = -1;
            end
        else
            % Calculate the sigma
            sigmaxp = (matrixSigma(iy,ix) + matrixSigma(iy,ix+1))/2;
            sigmaxm = (matrixSigma(iy,ix) + matrixSigma(iy, ix-1))/2;
            sigmayp = (matrixSigma(iy,ix) + matrixSigma(iy+1, ix))/2;
            sigmaym = (matrixSigma(iy,ix) + matrixSigma(iy-1, ix))/2;     

            % Calculate mapping index
            nxp = mappingEq(ix+1, iy, ny);  % index for V(i+1,j)
            nxm = mappingEq(ix-1, iy, ny);  % index for V(i-1,j)
            nyp = mappingEq(ix, iy+1, ny);  % index for V(i,j+1)
            nym = mappingEq(ix, iy-1, ny);  % index for V(i,j-1)

            % Setup the G matrix
            G(n,n) = -(sigmaxp+sigmaxm+sigmayp+sigmaym)/deltaXY^2;
            G(n, nxp) = sigmaxp/deltaXY^2;
            G(n, nxm) = sigmaxm/deltaXY^2;
            G(n, nyp) = sigmayp/deltaXY^2;
            G(n, nym) = sigmaym/deltaXY^2;
        end
    end
end

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
    matrixV(iy, ix) = V(iMap);
end

% Plot the solution for V from the Finite Difference Method
figure(2)
surf(X,Y,matrixV)
xlabel("X axis - Length")
ylabel("Y axis - Width")
zlabel("Z axis - Voltage")
title("3D Plot of Voltage V(x,y)")
view(45,45)  % 3-D View 
figure(3)
surf(X,Y,matrixV)
xlabel("X axis - Length")
ylabel("Y axis - Width")
zlabel("Z axis - Voltage")
title("2D Plot of Voltage V(x,y)")
view(0,90)  % View from top

% Solve the electric field
[Ex, Ey] = gradient(-matrixV);
Ex = Ex/deltaXY;
Ey = Ey/deltaXY;
% Plot the Ex field
figure(4)
surf(X,Y,Ex)
xlabel("X axis - Length")
ylabel("Y axis - Width")
zlabel("Z axis - Ex Field")
title("Plot of Ex(x,y)")
view(0,90)  % View from top
% Plot the Ey field
figure(5)
surf(X,Y,Ey)
xlabel("X axis - Length")
ylabel("Y axis - Width")
zlabel("Z axis - Ey Field")
title("Plot of Ey(x,y)")
view(0,90)  % View from top

% Plot the electric field
figure(6)
quiver(X, Y, Ex, Ey);
title("Plot of Electric Field E(x,y)")
xlabel("X axis - Length")
ylabel("Y axis - Width")

% Solve for current density field J = sigma*E
Jx = matrixSigma .* Ex;
Jy = matrixSigma .* Ey;
% Plot the current density
figure(7)
quiver(X,Y, Jx, Jy)
title("Plot of Current Density J(x,y)")
xlabel("X axis - Length")
ylabel("Y axis - Width")



% Create graph of current vs mesh size
vectorDeltaXY = [0.1, 0.05, 0.02, 0.01];
vectorCurrent = zeros(length(vectorDeltaXY), 1);
vectorMeshSize = zeros(length(vectorDeltaXY), 1);  % Hold the value of nx as mesh size
for indexMesh = 1:length(vectorDeltaXY)
    deltaXY = vectorDeltaXY(indexMesh);
    % Calculate the dimension of solution matrix
    nx = ceil(L/deltaXY);
    ny = ceil(W/deltaXY);
    [X,Y] = meshgrid(linspace(0,L,nx), linspace(0,W,ny));
    % Reconstruct the sigma matrix
    matrixSigma = ones(ny, nx);  % Dimension: ny times nx
    xIndexBox = ceil((L-Lb)/(2*deltaXY));  % Find the starting x index for the box
    LbIndexRange = ceil(Lb/deltaXY);  % Index range for the length of the box
    WbIndexRange = ceil(Wb/deltaXY);  % Index range for the width of the box
    % Assign the region for the box
    matrixSigma(1:WbIndexRange, xIndexBox:xIndexBox+LbIndexRange) = 10^-2;
    matrixSigma(ny-WbIndexRange:ny, xIndexBox:xIndexBox+LbIndexRange) = 10^-2;

    % Construct the G matrix and F vector
    G = zeros(nx*ny, nx*ny);  
    F = zeros(nx*ny, 1);
    matrixV = zeros(ny, nx);
    for ix = 1:nx
        for iy = 1:ny
            % Calculate the index
            n = mappingEq(ix, iy, ny);
            % Check for the boundary
            if ix==1 || ix==nx || iy ==1 || iy==ny
                G(n,n) = 1;
                % Boundary condition for x
                if ix == 1
                    F(n,1) = V0;  % V = V0 at x = 0 
                elseif ix == nx
                    F(n,1) = 0;  % and V = 0 at x = L
                elseif iy == 1
                    nyp = mappingEq(ix, iy+1, ny);  % dV/dy=0 at y=0
                    G(n,nyp) = -1;
                elseif iy == ny
                    nym = mappingEq(ix, iy-1, ny);  % dV/dy=0 at y=W
                    G(n, nym) = -1;
                end
            else
                % Calculate the sigma
                sigmaxp = (matrixSigma(iy,ix) + matrixSigma(iy,ix+1))/2;
                sigmaxm = (matrixSigma(iy,ix) + matrixSigma(iy, ix-1))/2;
                sigmayp = (matrixSigma(iy,ix) + matrixSigma(iy+1, ix))/2;
                sigmaym = (matrixSigma(iy,ix) + matrixSigma(iy-1, ix))/2;     

                % Calculate mapping index
                nxp = mappingEq(ix+1, iy, ny);  % index for V(i+1,j)
                nxm = mappingEq(ix-1, iy, ny);  % index for V(i-1,j)
                nyp = mappingEq(ix, iy+1, ny);  % index for V(i,j+1)
                nym = mappingEq(ix, iy-1, ny);  % index for V(i,j-1)

                % Setup the G matrix
                G(n,n) = -(sigmaxp+sigmaxm+sigmayp+sigmaym)/deltaXY^2;
                G(n, nxp) = sigmaxp/deltaXY^2;
                G(n, nxm) = sigmaxm/deltaXY^2;
                G(n, nyp) = sigmayp/deltaXY^2;
                G(n, nym) = sigmaym/deltaXY^2;
            end
        end
    end

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
        matrixV(iy, ix) = V(iMap);
    end

    % Solve the electric field
    [Ex, Ey] = gradient(-matrixV);
    Ex = Ex/deltaXY;
    Ey = Ey/deltaXY;
    % Solve for current density
    Jx = matrixSigma .* Ex;
    % Solve for current
    vectorCurrent(indexMesh) = sum(Jx(:,2)) * deltaXY;
    vectorMeshSize(indexMesh) = nx;  % Save the mesh size
end
% Plot for current vs. mesh size (nx) 
figure(8)
plot(vectorMeshSize, vectorCurrent)
title("Current vs. Mesh Size")
xlabel("Mesh Size (nx)")
ylabel("Current")
grid on





% Create graph of current vs. various bottle necks (different Wb with fixed Lb)
vectorWb = [0.1, 0.2, 0.3, 0.4, 0.45, 0.5]*W;  % Width of the box
Lb = 0.3*L;  % Length of the box
vectorCurrent = zeros(length(vectorWb), 1);
deltaXY = 0.02; 
% Calculate the dimension of solution matrix
nx = L/deltaXY;
ny = W/deltaXY;
[X,Y] = meshgrid(linspace(0,L,nx), linspace(0,W,ny));
for indexWb = 1:length(vectorWb)
    Wb = vectorWb(indexWb);
    % Reconstruct the sigma matrix
    matrixSigma = ones(ny, nx);  % Dimension: ny times nx
    xIndexBox = ceil((L-Lb)/(2*deltaXY));  % Find the starting x index for the box
    LbIndexRange = ceil(Lb/deltaXY);  % Index range for the length of the box
    WbIndexRange = ceil(Wb/deltaXY);  % Index range for the width of the box
    % Assign the region for the box
    matrixSigma(1:WbIndexRange, xIndexBox:xIndexBox+LbIndexRange) = 10^-2;
    matrixSigma(ny-WbIndexRange:ny, xIndexBox:xIndexBox+LbIndexRange) = 10^-2;    

    % Construct the G matrix and F vector
    G = zeros(nx*ny, nx*ny);  
    F = zeros(nx*ny, 1);
    matrixV = zeros(ny, nx);
    for ix = 1:nx
        for iy = 1:ny
            % Calculate the index
            n = mappingEq(ix, iy, ny);
            % Check for the boundary
            if ix==1 || ix==nx || iy ==1 || iy==ny
                G(n,n) = 1;
                % Boundary condition for x
                if ix == 1
                    F(n,1) = V0;  % V = V0 at x = 0 
                elseif ix == nx
                    F(n,1) = 0;  % and V = 0 at x = L
                elseif iy == 1
                    nyp = mappingEq(ix, iy+1, ny);  % dV/dy=0 at y=0
                    G(n,nyp) = -1;
                elseif iy == ny
                    nym = mappingEq(ix, iy-1, ny);  % dV/dy=0 at y=W
                    G(n, nym) = -1;
                end
            else
                % Calculate the sigma
                sigmaxp = (matrixSigma(iy,ix) + matrixSigma(iy,ix+1))/2;
                sigmaxm = (matrixSigma(iy,ix) + matrixSigma(iy, ix-1))/2;
                sigmayp = (matrixSigma(iy,ix) + matrixSigma(iy+1, ix))/2;
                sigmaym = (matrixSigma(iy,ix) + matrixSigma(iy-1, ix))/2;     

                % Calculate mapping index
                nxp = mappingEq(ix+1, iy, ny);  % index for V(i+1,j)
                nxm = mappingEq(ix-1, iy, ny);  % index for V(i-1,j)
                nyp = mappingEq(ix, iy+1, ny);  % index for V(i,j+1)
                nym = mappingEq(ix, iy-1, ny);  % index for V(i,j-1)

                % Setup the G matrix
                G(n,n) = -(sigmaxp+sigmaxm+sigmayp+sigmaym)/deltaXY^2;
                G(n, nxp) = sigmaxp/deltaXY^2;
                G(n, nxm) = sigmaxm/deltaXY^2;
                G(n, nyp) = sigmayp/deltaXY^2;
                G(n, nym) = sigmaym/deltaXY^2;
            end
        end
    end

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
        matrixV(iy, ix) = V(iMap);
    end

    % Solve the electric field
    [Ex, Ey] = gradient(-matrixV);
    Ex = Ex/deltaXY;
    Ey = Ey/deltaXY;
    % Solve for current density
    Jx = matrixSigma .* Ex;
    % Solve for current
    vectorCurrent(indexWb) = sum(Jx(:,2)) * deltaXY;
end
% Plot the current vs various bottle-neck (different Wb with fixed Lb)
figure(9)
plot(vectorWb, vectorCurrent)
title("Current vs various bottle-neck (different Wb with fixed Lb)")
xlabel("Width of the box (Wb)")
ylabel("Current")
grid on








% Create graph of current vs. various bottle necks (different Lb with fixed Wb)
vectorLb = [0.1, 0.2, 0.3, 0.4, 0.45, 0.5]*L;  % Length of the box
Wb = 0.4*W;  % Width of the box
vectorCurrent = zeros(length(vectorLb), 1);
deltaXY = 0.02; 
% Calculate the dimension of solution matrix
nx = L/deltaXY;
ny = W/deltaXY;
[X,Y] = meshgrid(linspace(0,L,nx), linspace(0,W,ny));
for indexLb = 1:length(vectorLb)
    Lb = vectorLb(indexLb);
    % Reconstruct the sigma matrix
    matrixSigma = ones(ny, nx);  % Dimension: ny times nx
    xIndexBox = ceil((L-Lb)/(2*deltaXY));  % Find the starting x index for the box
    LbIndexRange = ceil(Lb/deltaXY);  % Index range for the length of the box
    WbIndexRange = ceil(Wb/deltaXY);  % Index range for the width of the box
    % Assign the region for the box
    matrixSigma(1:WbIndexRange, xIndexBox:xIndexBox+LbIndexRange) = 10^-2;
    matrixSigma(ny-WbIndexRange:ny, xIndexBox:xIndexBox+LbIndexRange) = 10^-2;    

    % Construct the G matrix and F vector
    G = zeros(nx*ny, nx*ny);  
    F = zeros(nx*ny, 1);
    matrixV = zeros(ny, nx);
    for ix = 1:nx
        for iy = 1:ny
            % Calculate the index
            n = mappingEq(ix, iy, ny);
            % Check for the boundary
            if ix==1 || ix==nx || iy ==1 || iy==ny
                G(n,n) = 1;
                % Boundary condition for x
                if ix == 1
                    F(n,1) = V0;  % V = V0 at x = 0 
                elseif ix == nx
                    F(n,1) = 0;  % and V = 0 at x = L
                elseif iy == 1
                    nyp = mappingEq(ix, iy+1, ny);  % dV/dy=0 at y=0
                    G(n,nyp) = -1;
                elseif iy == ny
                    nym = mappingEq(ix, iy-1, ny);  % dV/dy=0 at y=W
                    G(n, nym) = -1;
                end
            else
                % Calculate the sigma
                sigmaxp = (matrixSigma(iy,ix) + matrixSigma(iy,ix+1))/2;
                sigmaxm = (matrixSigma(iy,ix) + matrixSigma(iy, ix-1))/2;
                sigmayp = (matrixSigma(iy,ix) + matrixSigma(iy+1, ix))/2;
                sigmaym = (matrixSigma(iy,ix) + matrixSigma(iy-1, ix))/2;     

                % Calculate mapping index
                nxp = mappingEq(ix+1, iy, ny);  % index for V(i+1,j)
                nxm = mappingEq(ix-1, iy, ny);  % index for V(i-1,j)
                nyp = mappingEq(ix, iy+1, ny);  % index for V(i,j+1)
                nym = mappingEq(ix, iy-1, ny);  % index for V(i,j-1)

                % Setup the G matrix
                G(n,n) = -(sigmaxp+sigmaxm+sigmayp+sigmaym)/deltaXY^2;
                G(n, nxp) = sigmaxp/deltaXY^2;
                G(n, nxm) = sigmaxm/deltaXY^2;
                G(n, nyp) = sigmayp/deltaXY^2;
                G(n, nym) = sigmaym/deltaXY^2;
            end
        end
    end

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
        matrixV(iy, ix) = V(iMap);
    end

    % Solve the electric field
    [Ex, Ey] = gradient(-matrixV);
    Ex = Ex/deltaXY;
    Ey = Ey/deltaXY;
    % Solve for current density
    Jx = matrixSigma .* Ex;
    % Solve for current
    vectorCurrent(indexLb) = sum(Jx(:,2)) * deltaXY;
end
% Plot the current vs various bottle-neck (different Lb with fixed Wb)
figure(10)
plot(vectorLb, vectorCurrent)
title("Current vs various bottle-neck (different Lb with fixed Wb)")
xlabel("Length of the box (Lb)")
ylabel("Current")
grid on




% Graph of current vs sigma inside the box
vectorSigma = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001];  % Length of the box
Wb = 0.4*W;  % Width of the box
Lb = 0.3*L;  % Length of the box
vectorCurrent = zeros(length(vectorSigma), 1);
deltaXY = 0.02; 
% Calculate the dimension of solution matrix
nx = L/deltaXY;
ny = W/deltaXY;
[X,Y] = meshgrid(linspace(0,L,nx), linspace(0,W,ny));
for indexSigma = 1:length(vectorSigma)
    sigmaInBox = vectorSigma(indexSigma);
    % Reconstruct the sigma matrix
    matrixSigma = ones(ny, nx);  % Dimension: ny times nx
    xIndexBox = ceil((L-Lb)/(2*deltaXY));  % Find the starting x index for the box
    LbIndexRange = ceil(Lb/deltaXY);  % Index range for the length of the box
    WbIndexRange = ceil(Wb/deltaXY);  % Index range for the width of the box
    % Assign the region for the box
    matrixSigma(1:WbIndexRange, xIndexBox:xIndexBox+LbIndexRange) = sigmaInBox;
    matrixSigma(ny-WbIndexRange:ny, xIndexBox:xIndexBox+LbIndexRange) = sigmaInBox;    

    % Construct the G matrix and F vector
    G = zeros(nx*ny, nx*ny);  
    F = zeros(nx*ny, 1);
    matrixV = zeros(ny, nx);
    for ix = 1:nx
        for iy = 1:ny
            % Calculate the index
            n = mappingEq(ix, iy, ny);
            % Check for the boundary
            if ix==1 || ix==nx || iy ==1 || iy==ny
                G(n,n) = 1;
                % Boundary condition for x
                if ix == 1
                    F(n,1) = V0;  % V = V0 at x = 0 
                elseif ix == nx
                    F(n,1) = 0;  % and V = 0 at x = L
                elseif iy == 1
                    nyp = mappingEq(ix, iy+1, ny);  % dV/dy=0 at y=0
                    G(n,nyp) = -1;
                elseif iy == ny
                    nym = mappingEq(ix, iy-1, ny);  % dV/dy=0 at y=W
                    G(n, nym) = -1;
                end
            else
                % Calculate the sigma
                sigmaxp = (matrixSigma(iy,ix) + matrixSigma(iy,ix+1))/2;
                sigmaxm = (matrixSigma(iy,ix) + matrixSigma(iy, ix-1))/2;
                sigmayp = (matrixSigma(iy,ix) + matrixSigma(iy+1, ix))/2;
                sigmaym = (matrixSigma(iy,ix) + matrixSigma(iy-1, ix))/2;     

                % Calculate mapping index
                nxp = mappingEq(ix+1, iy, ny);  % index for V(i+1,j)
                nxm = mappingEq(ix-1, iy, ny);  % index for V(i-1,j)
                nyp = mappingEq(ix, iy+1, ny);  % index for V(i,j+1)
                nym = mappingEq(ix, iy-1, ny);  % index for V(i,j-1)

                % Setup the G matrix
                G(n,n) = -(sigmaxp+sigmaxm+sigmayp+sigmaym)/deltaXY^2;
                G(n, nxp) = sigmaxp/deltaXY^2;
                G(n, nxm) = sigmaxm/deltaXY^2;
                G(n, nyp) = sigmayp/deltaXY^2;
                G(n, nym) = sigmaym/deltaXY^2;
            end
        end
    end

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
        matrixV(iy, ix) = V(iMap);
    end

    % Solve the electric field
    [Ex, Ey] = gradient(-matrixV);
    Ex = Ex/deltaXY;
    Ey = Ey/deltaXY;
    % Solve for current density
    Jx = matrixSigma .* Ex;
    % Solve for current
    vectorCurrent(indexSigma) = sum(Jx(:,2)) * deltaXY;
end
% Plot of current vs sigma
figure(11)
plot(vectorSigma, vectorCurrent)
title("Current versus Conductivity (Sigma)")
xlabel("Conductivity (Sigma)")
ylabel("Current")
grid on

