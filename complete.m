clear all; close all; clc;

%% --- MAIN SETUP ---
% Define Input Angle (Yellow Link, Loop 1)
theta_in_deg = 19.94;
theta_in = deg2rad(theta_in_deg); % Renamed to 'theta_in' to match the solver

% Create Figure
figure('Color','w','Position',[50 50 1200 600]); 

% --- PLOT CASE 1: OPEN CONFIGURATION ---
subplot(1,2,1); hold on; axis equal; grid on; box on;
title('Case 1: Open Configuration');
xlabel('X (mm)'); ylabel('Y (mm)');
plot_full_mechanism(theta_in, 1); 

% --- PLOT CASE 2: CROSSED CONFIGURATION ---
subplot(1,2,2); hold on; axis equal; grid on; box on;
title('Case 2: Crossed Configuration');
xlabel('X (mm)'); ylabel('Y (mm)');
plot_full_mechanism(theta_in, 2); 

%% --- FUNCTION: Calculate and Plot ---
function plot_full_mechanism(theta_in, config_mode)
    
    %% 1. LOOP 1 CALCULATION (Green-Yellow-Grey)
    d1 = 210; a1 = 180; b1 = 180; c1 = 118; 
    
    % Solver: Input Yellow (b1) -> Find Green (a1)
    Px = b1 * cos(theta_in) - d1;
    Py = b1 * sin(theta_in);
    P_mag_sq = Px^2 + Py^2;
    
    A1 = 2*a1*Py; 
    B1 = 2*a1*Px; 
    C1 = a1^2 + P_mag_sq - c1^2;
    
    [th_Green_1, th_Green_2, ok1] = solve_angle_all(A1, B1, C1);
    
    if ~ok1
        text(0,0,'Loop 1 Failed (No Assembly)'); 
        return; 
    end
    
    % Select Solution based on Config Mode or Geometry
    % Case 1 logic: Green usually below ground for this mechanism style
    if sin(th_Green_1) < 0
        th_Green = th_Green_1;
    else
        th_Green = th_Green_2;
    end
    
    % Output Loop 1 (Grey)
    % Calculate vector explicitly to get the angle
    R_Grey_Vec = a1*exp(1j*th_Green) + b1*exp(1j*theta_in) - d1;
    th_Grey_L1 = angle(R_Grey_Vec);
    
    %% 2. LOOP 2 CALCULATION (Brown-Blue-LightBlue)
    d2 = 210; a2 = 118; b2 = 210; c2 = 236; 
    
    % Input: Brown follows Grey L1
    theta_Brown = th_Grey_L1; 
    
    % Solver: Input Brown (a2) -> Find LightBlue (c2)
    % Using Vector Loop Method (Robust)
    X2 = d2 - a2 * cos(theta_Brown);
    Y2 = -a2 * sin(theta_Brown);
    
    A2 = 2 * c2 * Y2;
    B2 = 2 * c2 * X2;
    C2 = X2^2 + Y2^2 + c2^2 - b2^2;
    
    [th_LB_1, th_LB_2, ok2] = solve_angle_all(A2, B2, C2);
    
    if ~ok2
        text(0,0,'Loop 2 Failed'); 
        return; 
    end
    
    % Selection Logic for Loop 2
    if config_mode == 1
        % Generally "Open" or specific geometric constraint
        % Heuristic: Pick the one that puts the joint 'lower' if needed
        if sin(th_LB_1) < sin(th_LB_2), th_LightBlue = th_LB_1; else, th_LightBlue = th_LB_2; end
    else
        % Crossed
        if sin(th_LB_1) > sin(th_LB_2), th_LightBlue = th_LB_1; else, th_LightBlue = th_LB_2; end
    end
    
    % Calculate Coupler (Blue)
    R_Blue_Vec = d2 + c2*exp(1j*th_LightBlue) - a2*exp(1j*theta_Brown);
    th_Blue = angle(R_Blue_Vec);

    %% 3. LOOP 3 CALCULATION (LightBlue-Red-Grey)
    d3 = 210; a3 = 118; b3 = 210; c3 = 118; 
    
    % Input: Light Blue from Loop 2
    theta_LB_In = th_LightBlue;
    
    % Solver: Input LightBlue (a3) -> Find Grey (c3)
    X3 = d3 - a3 * cos(theta_LB_In);
    Y3 = -a3 * sin(theta_LB_In);
    
    A3 = 2 * c3 * Y3;
    B3 = 2 * c3 * X3;
    C3 = X3^2 + Y3^2 + c3^2 - b3^2;
    
    [th_Grey_L3_1, th_Grey_L3_2, ok3] = solve_angle_all(A3, B3, C3);
    
    if ~ok3
        text(0,0,'Loop 3 Failed'); 
        return; 
    end
    
    % Selection Logic for Loop 3
    % We need to check which solution aligns best with the structure
    % Usually we define this by where the Red link ends up
    th_Grey_L3 = th_Grey_L3_1; % Default
    th_Red     = angle(d3 + c3*exp(1j*th_Grey_L3_1) - a3*exp(1j*theta_LB_In));
    
    % Simple heuristic: If config 2, swap
    if config_mode == 2
         th_Grey_L3 = th_Grey_L3_2;
         th_Red     = angle(d3 + c3*exp(1j*th_Grey_L3_2) - a3*exp(1j*theta_LB_In));
    end

    %% --- PLOTTING ---
    % Define Vectors
    R_Green   = a1 * exp(1j*th_Green);
    R_Yellow  = b1 * exp(1j*theta_in);
    
    R_Brown   = a2 * exp(1j*theta_Brown);
    R_Blue    = b2 * exp(1j*th_Blue);
    R_LB_L2   = c2 * exp(1j*th_LightBlue); 
    
    R_LB_L3   = a3 * exp(1j*theta_LB_In); 
    R_Red     = b3 * exp(1j*th_Red);
    R_Grey_L3 = c3 * exp(1j*th_Grey_L3);
    
    % Draw Ground
    plot([0 d1], [0 0], 'Color', [1 0.07 0.57], 'LineWidth', 2, 'LineStyle', '--');
    text(0, -20, 'O2'); text(d1, -20, 'O4');
    
    % Draw Loop 1
    quiver(0,0, real(R_Green), imag(R_Green), 0, 'g', 'LineWidth', 3, 'MaxHeadSize',0.4); 
    quiver(real(R_Green), imag(R_Green), real(R_Yellow), imag(R_Yellow), 0, 'y', 'LineWidth', 3, 'MaxHeadSize',0.4); 
    
    % Draw Loop 2
    quiver(d1, 0, real(R_Brown), imag(R_Brown), 0, 'Color', [0.6 0.3 0], 'LineWidth', 3, 'MaxHeadSize',0.4); 
    J_Brown = d1 + R_Brown;
    quiver(real(J_Brown), imag(J_Brown), real(R_Blue), imag(R_Blue), 0, 'b', 'LineWidth', 3, 'MaxHeadSize',0.4); 
    quiver(0, 0, real(R_LB_L2), imag(R_LB_L2), 0, 'c', 'LineWidth', 3, 'MaxHeadSize',0.4); 
    
    % Draw Loop 3
    J_LB_Mid = R_LB_L3; % Start of Red
    quiver(real(J_LB_Mid), imag(J_LB_Mid), real(R_Red), imag(R_Red), 0, 'r', 'LineWidth', 3, 'MaxHeadSize',0.4); 
    quiver(d1, 0, real(R_Grey_L3), imag(R_Grey_L3), 0, 'Color', [0.3 0.3 0.3], 'LineWidth', 3, 'MaxHeadSize',0.4); 
    
    % Joints
    plot(real(R_Green), imag(R_Green), 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
    plot(real(J_Brown), imag(J_Brown), 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
    plot(real(J_LB_Mid), imag(J_LB_Mid), 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
    
    xlim([-100 400]); ylim([-350 250]);
end

%% --- HELPER: Solve Angle ---
function [th1, th2, status] = solve_angle_all(A, B, C)
    % Solves A*sin(th) + B*cos(th) + C = 0
    % Using t = tan(th/2) substitution
    
    Aq = C - B; 
    Bq = 2 * A; 
    Cq = C + B;
    
    det_val = Bq^2 - 4*Aq*Cq;
    
    if det_val < 0
        th1 = 0; th2 = 0; status = 0; return; 
    end
    
    t1 = (-Bq + sqrt(det_val)) / (2*Aq);
    t2 = (-Bq - sqrt(det_val)) / (2*Aq);
    
    th1 = 2*atan(t1);
    th2 = 2*atan(t2);
    status = 1;
end
