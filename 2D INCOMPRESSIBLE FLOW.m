% simulateFlow.m
% Simulates 2D incompressible flow with particle advection
% Uses Jacobi for pressure projection and RK4 for advection of both
% velocity field (semi-Lagrangian) and Lagrangian tracer particles.

clear; clc; close all;
%% 1. Simulation parameters
s  = 200;                          % grid size (rows)
ar = 2;                            % aspect ratio (cols = ar * rows)
J  = [0 1 0; 1 0 1; 0 1 0]/4;      % Jacobi stencil
dt = 1;                            % time step (for clarity)
% create coordinate grid
[X, Y] = meshgrid(1:s*ar, 1:s);
% initialize fields
p  = zeros(s, s*ar);    % pressure
vx = zeros(s, s*ar);    % x-velocity
vy = zeros(s, s*ar);    % y-velocity

%% 2. initialize tracer particles - making a larger, more visible block
[px, py] = meshgrid(10:20, 90:110);
px = px(:);  py = py(:);
px0 = px;    py0 = py;    % keep for inflow

%% 3. set up figure
f = figure('Name','2D Flow & Particles','NumberTitle','off', 'Position', [100, 100, 800, 400]);
axis([0, s*ar, 0, s]);

%% 4. main loop
frameCount = 0;
while ishandle(f) && frameCount < 500  % Add a maximum frame count for safety
    frameCount = frameCount + 1;
    
    % 4.1 inject a stronger horizontal jet for better visibility
    vx(90:110, 10:20) = 1.0;  % Increased velocity for better movement
    
    % 4.2 build RHS (-divergence of provisional velocity)
    [dvx_dx, ~] = gradient(vx);
    [~, dvy_dy] = gradient(vy);
    divv = dvx_dx + dvy_dy;
    rhs = -divv;
    
    % 4.3 solve Poisson via Jacobi (50 iterations - reduced for speed)
    nIter = 50;
    for k = 1:nIter
        p = conv2(p, J, 'same') + rhs/2;
    end
    
    % 4.4 pressure-correction: project onto ∇·v=0
    [dpdx, dpdy] = gradient(p);
    vx(2:end-1,2:end-1) = vx(2:end-1,2:end-1) - dpdx(2:end-1,2:end-1);
    vy(2:end-1,2:end-1) = vy(2:end-1,2:end-1) - dpdy(2:end-1,2:end-1);
    
    % 4.5 advect the velocity field backward (semi-Lagrangian)
    % RK4 for velocity field
    k1x = interp2(vx, X, Y, 'linear', 0);
    k1y = interp2(vy, X, Y, 'linear', 0);
    k2x = interp2(vx, X - 0.5*dt*k1x, Y - 0.5*dt*k1y, 'linear', 0);
    k2y = interp2(vy, X - 0.5*dt*k1x, Y - 0.5*dt*k1y, 'linear', 0);
    k3x = interp2(vx, X - 0.5*dt*k2x, Y - 0.5*dt*k2y, 'linear', 0);
    k3y = interp2(vy, X - 0.5*dt*k2x, Y - 0.5*dt*k2y, 'linear', 0);
    k4x = interp2(vx, X - dt*k3x, Y - dt*k3y, 'linear', 0);
    k4y = interp2(vy, X - dt*k3x, Y - dt*k3y, 'linear', 0);
    PX = X - dt*(k1x + 2*k2x + 2*k3x + k4x)/6;
    PY = Y - dt*(k1y + 2*k2y + 2*k3y + k4y)/6;
    vx = interp2(vx, PX, PY, 'linear', 0);
    vy = interp2(vy, PX, PY, 'linear', 0);
    
    % 4.6 advect tracer particles forward
    k1x = interp2(vx, px, py, 'linear', 0);
    k1y = interp2(vy, px, py, 'linear', 0);
    k2x = interp2(vx, px + 0.5*dt*k1x, py + 0.5*dt*k1y, 'linear', 0);
    k2y = interp2(vy, px + 0.5*dt*k1x, py + 0.5*dt*k1y, 'linear', 0);
    k3x = interp2(vx, px + 0.5*dt*k2x, py + 0.5*dt*k2y, 'linear', 0);
    k3y = interp2(vy, px + 0.5*dt*k2x, py + 0.5*dt*k2y, 'linear', 0);
    k4x = interp2(vx, px + dt*k3x, py + dt*k3y, 'linear', 0);
    k4y = interp2(vy, px + dt*k3x, py + dt*k3y, 'linear', 0);
    px = px + dt*(k1x + 2*k2x + 2*k3x + k4x)/6;
    py = py + dt*(k1y + 2*k2y + 2*k3y + k4y)/6;
    
    % 4.7 re-inject inflow particles every 10 frames to maintain density
    if mod(frameCount, 10) == 0
        px = [px; px0];
        py = [py; py0];
    end
    
    % 4.8 plot with increased point size and color for better visibility
    cla;  % Clear previous particles
    scatter(px, py, 10, 'filled', 'MarkerFaceColor', 'b');
    
    % Optional: add velocity field visualization (uncomment to see)
    % hold on;
    % [xq, yq] = meshgrid(20:20:s*ar, 20:20:s);
    % vxq = interp2(vx, xq, yq);
    % vyq = interp2(vy, xq, yq);
    % quiver(xq, yq, vxq, vyq, 'r');
    % hold off;
    
    axis equal tight;
    axis([0, s*ar, 0, s]);
    title(['2D Incompressible Flow with Particle Advection (Frame: ' num2str(frameCount) ')']);
    drawnow;
    
    % Check if there are any particles still in frame
    visibleParticles = sum((px > 0 & px < s*ar & py > 0 & py < s));
    if visibleParticles < 10 && frameCount > 50
        disp(['Warning: Only ' num2str(visibleParticles) ' particles visible. Adding more.']);
        px = [px; px0];
        py = [py; py0];
    end
end