clear all; close all; clc;

%% --- MAIN CONFIGURATION ---
% Input Angle (Yellow Link, Loop 1)
theta_in_deg = 19.94;
theta_in = deg2rad(theta_in_deg); 

% Figure Setup
figure('Color','w','Position',[50 50 1200 600]); 

%% --- PLOT 1: CASE 1 (Target: Below Ground / Open) ---
subplot(1,2,1); hold on; axis equal; grid on; box on;
title('Case 1: Target Config (All Below Ground)');
xlabel('X (mm)'); ylabel('Y (mm)');
plot_mechanism(theta_in, 1); % Mode 1 = Try to pick 'Below Ground' solution

%% --- PLOT 2: CASE 2 (Alternative / Crossed) ---
subplot(1,2,2); hold on; axis equal; grid on; box on;
title('Case 2: Alternative Config (Crossed)');
xlabel('X (mm)'); ylabel('Y (mm)');
plot_mechanism(theta_in, 2); % Mode 2 = Pick alternative roots


%% =========================================================================
%% FUNCTION: SOLVE AND PLOT MECHANISM
%% =========================================================================
function plot_mechanism(theta_in, mode)

    %% --- 1. SOLVE LOOP 1 (Green-Yellow-Grey) ---
    d1 = 210; a1 = 180; b1 = 180; c1 = 118; 
    
    % Solver: Input Yellow (b1) -> Find Green (a1)
    % Equation: a1*exp(j*th2) - c1*exp(j*th4) = d1 - b1*exp(j*th3)
    Px = b1 * cos(theta_in) - d1;
    Py = b1 * sin(theta_in);
    P_mag_sq = Px^2 + Py^2;
    
    A1 = 2*a1*Py; 
    B1 = 2*a1*Px; 
    C1 = a1^2 + P_mag_sq - c1^2;
    
    [th_Green_1, th_Green_2, ok1] = solve_angle_eq(A1, B1, C1);
    if ~ok1, text(0,0,'Loop 1 Failed'); return; end
    
    % Selection L1: Always pick Green Below Ground (Assumption for both cases)
    if sin(th_Green_1) < 0
        th_Green = th_Green_1;
    else
        th_Green = th_Green_2;
    end
    
    % Find Grey L1 (Output)
    R_Grey_Vec_L1 = a1*exp(1j*th_Green) + b1*exp(1j*theta_in) - d1;
    th_Grey_L1 = angle(R_Grey_Vec_L1);

    
    %% --- 2. SOLVE LOOP 2 (Brown-Blue-LightBlue) ---
    d2 = 210; a2 = 118; b2 = 210; 
    c2 = 118; % *** Radius of Light Blue (236/2) ***
    
    % Input: Brown follows Grey L1
    theta_Brown = th_Grey_L1; 
    
    % Solver: Input Brown (a2) -> Find LightBlue (c2)
    X2 = d2 - a2 * cos(theta_Brown);
    Y2 = -a2 * sin(theta_Brown);
    
    A2 = 2 * c2 * Y2;
    B2 = 2 * c2 * X2;
    C2 = X2^2 + Y2^2 + c2^2 - b2^2;
    
    [th_LB_1, th_LB_2, ok2] = solve_angle_eq(A2, B2, C2);
    if ~ok2, text(0,0,'Loop 2 Failed'); return; end
    
    % Selection L2:
    if mode == 1
        % Case 1: Pick 'Below Ground' (y < 0)
        if sin(th_LB_1) < 0 && sin(th_LB_2) < 0
             % Both down, pick geometric open (usually smaller angle diff)
             % Or just picking the first valid negative
             if sin(th_LB_1) < sin(th_LB_2), th_LightBlue = th_LB_1; else, th_LightBlue = th_LB_2; end
        elseif sin(th_LB_1) < 0
             th_LightBlue = th_LB_1;
        else
             th_LightBlue = th_LB_2;
        end
    else
        % Case 2: Pick the OTHER one (Crossed)
        % First find the 'Below' one, then pick the other
        if sin(th_LB_1) < 0 && sin(th_LB_2) > 0
             th_LightBlue = th_LB_2; % Pick the upper one
        elseif sin(th_LB_2) < 0
             th_LightBlue = th_LB_1;
        else
             th_LightBlue = th_LB_1; % Fallback
        end
    end
    
    % Find Blue (Coupler)
    R_Blue_Vec = d2 + c2*exp(1j*th_LightBlue) - a2*exp(1j*theta_Brown);
    th_Blue = angle(R_Blue_Vec);

    
    %% --- 3. SOLVE LOOP 3 (LightBlue-Red-Grey) ---
    d3 = 210; 
    a3 = 118; % *** Radius of Light Blue (236/2) *** b3 = 210; c3 = 118; 
    
    % Input: Light Blue from Loop 2 (SAME LINK)
    theta_LB_In = th_LightBlue;
    
    % Solver: Input LightBlue (a3) -> Find Grey (c3)
    X3 = d3 - a3 * cos(theta_LB_In);
    Y3 = -a3 * sin(theta_LB_In);
    
    A3 = 2 * c3 * Y3;
    B3 = 2 * c3 * X3;
    C3 = X3^2 + Y3^2 + c3^2 - b3^2;
    
    [th_Grey_L3_1, th_Grey_L3_2, ok3] = solve_angle_eq(A3, B3, C3);
    if ~ok3, text(0,0,'Loop 3 Failed'); return; end
    
    % Selection L3:
    if mode == 1
        % Case 1: Pick 'Below Ground'
        if sin(th_Grey_L3_1) < 0, th_Grey_L3 = th_Grey_L3_1; else, th_Grey_L3 = th_Grey_L3_2; end
    else
        % Case 2: Pick Alternative
        if sin(th_Grey_L3_1) >= 0, th_Grey_L3 = th_Grey_L3_1; else, th_Grey_L3 = th_Grey_L3_2; end
    end
    
    % Find Red (Coupler)
    R_Red_Vec = d3 + c3*exp(1j*th_Grey_L3) - a3*exp(1j*theta_LB_In);
    th_Red = angle(R_Red_Vec);
    
    %% --- 4. PLOTTING ---
    % Define Vectors
    R_Ground = d1 * exp(1j*0);
    
    % Loop 1
    R_Green   = a1 * exp(1j*th_Green);
    R_Yellow  = b1 * exp(1j*theta_in);
    
    % Loop 2
    R_Brown   = a2 * exp(1j*theta_Brown);
    R_Blue    = b2 * exp(1j*th_Blue);
    R_LB_L2   = c2 * exp(1j*th_LightBlue); 
    
    % Loop 3
    % Light Blue L3 is same vector as L2 (same angle, same length 118)
    R_LB_L3   = a3 * exp(1j*theta_LB_In); 
    R_Red     = b3 * exp(1j*th_Red);
    R_Grey_L3 = c3 * exp(1j*th_Grey_L3);
    
    % --- Draw Geometry ---
    % Ground (Pink)
    plot([0 d1], [0 0], 'Color', [1 0.07 0.57], 'LineWidth', 2, 'LineStyle', '--');
    text(0, -20, 'O2'); text(d1, -20, 'O4');
    
    % Loop 1
    quiver(0,0, real(R_Green), imag(R_Green), 0, 'g', 'LineWidth', 3, 'MaxHeadSize',0.4); 
    quiver(real(R_Green), imag(R_Green), real(R_Yellow), imag(R_Yellow), 0, 'y', 'LineWidth', 3, 'MaxHeadSize',0.4); 
    
    % Loop 2
    quiver(d1, 0, real(R_Brown), imag(R_Brown), 0, 'Color', [0.6 0.3 0], 'LineWidth', 3, 'MaxHeadSize',0.4); 
    J_Brown = d1 + R_Brown;
    quiver(real(J_Brown), imag(J_Brown), real(R_Blue), imag(R_Blue), 0, 'b', 'LineWidth', 3, 'MaxHeadSize',0.4); 
    
    % Light Blue (Unified Link)
    % วาดเส้นเดียวจาก O2 ไปปลาย (เพราะ Loop 2 และ Loop 3 ใช้เวกเตอร์เดียวกัน)
    quiver(0, 0, real(R_LB_L2), imag(R_LB_L2), 0, 'c', 'LineWidth', 4, 'MaxHeadSize',0.4); 
    
    % Loop 3
    % Red starts from tip of Light Blue
    J_LB_Tip = R_LB_L3; 
    quiver(real(J_LB_Tip), imag(J_LB_Tip), real(R_Red), imag(R_Red), 0, 'r', 'LineWidth', 3, 'MaxHeadSize',0.4); 
    quiver(d1, 0, real(R_Grey_L3), imag(R_Grey_L3), 0, 'Color', [0.3 0.3 0.3], 'LineWidth', 3, 'MaxHeadSize',0.4); 
    
    % Joints
    plot(real(R_Green), imag(R_Green), 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
    plot(real(J_Brown), imag(J_Brown), 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
    plot(real(J_LB_Tip), imag(J_LB_Tip), 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
    
    % Limits
    xlim([-100 400]); ylim([-350 200]);
end

%% --- HELPER: SOLVE ANGLE EQ ---
function [th1, th2, status] = solve_angle_eq(A, B, C)
    % Solves A sin(th) + B cos(th) + C = 0
    Aq = C - B; Bq = 2*A; Cq = C + B;
    det = Bq^2 - 4*Aq*Cq;
    
    if det < 0, th1=0; th2=0; status=0; return; end
    
    t1 = (-Bq + sqrt(det))/(2*Aq);
    t2 = (-Bq - sqrt(det))/(2*Aq);
    
    th1 = 2*atan(t1);
    th2 = 2*atan(t2);
    status = 1;
end
