clear all; close all; clc;

%% --- MAIN SETUP ---
% กำหนดค่ามุม Input หลัก (ก้านสีเหลือง Loop 1)
theta_Yellow_Input_Deg = 19.94;
theta_Yellow_Input_Rad = deg2rad(theta_Yellow_Input_Deg);

% สร้าง Figure รวม
figure('Color','w','Position',[50 50 1200 600]); 

% --- PLOT CASE 1: OPEN CONFIGURATION (ตรงตามรูปวาดของคุณ) ---
subplot(1,2,1); hold on; axis equal; grid on; box on;
title('Complete Mechanism: Case 1 (Open Config - Matches Sketch)');
xlabel('X (mm)'); ylabel('Y (mm)');
plot_full_mechanism(theta_Yellow_Input_Rad, 1); % เรียกฟังก์ชันวาด Case 1

% --- PLOT CASE 2: CROSSED CONFIGURATION (เผื่อไว้ดูเปรียบเทียบ) ---
subplot(1,2,2); hold on; axis equal; grid on; box on;
title('Complete Mechanism: Case 2 (Crossed Config)');
xlabel('X (mm)'); ylabel('Y (mm)');
plot_full_mechanism(theta_Yellow_Input_Rad, 2); % เรียกฟังก์ชันวาด Case 2


%% --- FUNCTION: Calculate and Plot Everything Together ---
function plot_full_mechanism(theta_in, config_mode)
    
    %% 1. LOOP 1 CALCULATION (Green-Yellow-Grey)
    d1 = 210; a1 = 180; b1 = 180; c1 = 118; % Ground, Green, Yellow, Grey
    
    % Solver (Yellow Input)
    Px = b1 * cos(theta_in) - d1;
    Py = b1 * sin(theta_in);
    P_mag_sq = Px^2 + Py^2;
    
    A1 = 2*a1*Py; B1 = 2*a1*Px; C1 = a1^2 + P_mag_sq - c1^2;
    [th_Green, ok1] = solve_angle(A1, B1, C1, config_mode);
    if ~ok1, text(0,0,'Loop 1 Failed'); return; end
    
    % Output Loop 1
    th_Grey_L1 = angle(a1*exp(1j*th_Green) + b1*exp(1j*theta_in) - d1);

    %% 2. LOOP 2 CALCULATION (Brown-Blue-LightBlue)
    d2 = 210; a2 = 118; b2 = 210; c2 = 168; % Ground, Brown, Blue, LightBlue
    
    % Input Mapping: Brown รับมุมมาจาก Grey ของ Loop 1
    theta_Brown = th_Grey_L1; 
    
    % Solver (Brown Input)
    K1=d2/a2; K2=d2/c2; K3=(a2^2-b2^2+c2^2+d2^2)/(2*a2*c2);
    A2 = cos(theta_Brown)-K1-K2*cos(theta_Brown)+K3;
    B2 = -2*sin(theta_Brown);
    C2 = K1-(K2+1)*cos(theta_Brown)+K3;
    
    [th_LightBlue, ok2] = solve_angle(A2, B2, C2, config_mode);
    if ~ok2, text(0,0,'Loop 2 Failed'); return; end
    
    % Calculate Coupler (Blue)
    th_Blue = angle(d2 + c2*exp(1j*th_LightBlue) - a2*exp(1j*theta_Brown));

    %% 3. LOOP 3 CALCULATION (LightBlue-Red-Grey)
    d3 = 210; a3 = 118; b3 = 210; c3 = 118; % Ground, LightBlue, Red, Grey
    
    % Input Mapping: Light Blue จาก Loop 2 (แต่ยาว 118 ใน Loop นี้)
    % *** สำคัญ: ใช้มุม th_LightBlue เดียวกัน เพราะคือก้านเดียวกัน ***
    theta_LB_In = th_LightBlue;
    
    % Solver (LightBlue Input)
    K1=d3/a3; K2=d3/c3; K3=(a3^2-b3^2+c3^2+d3^2)/(2*a3*c3);
    A3 = cos(theta_LB_In)-K1-K2*cos(theta_LB_In)+K3;
    B3 = -2*sin(theta_LB_In);
    C3 = K1-(K2+1)*cos(theta_LB_In)+K3;
    
    [th_Grey_L3, ok3] = solve_angle(A3, B3, C3, config_mode);
    if ~ok3, text(0,0,'Loop 3 Failed'); return; end
    
    % Calculate Coupler (Red)
    th_Red = angle(d3 + c3*exp(1j*th_Grey_L3) - a3*exp(1j*theta_LB_In));

    %% --- PLOTTING ALL LINKS ---
    
    % Define Vectors
    R_Ground = d1 * exp(1j*0);
    
    % Loop 1 Vectors
    R_Green = a1 * exp(1j*th_Green);
    R_Yellow = b1 * exp(1j*theta_in);
    R_Grey_L1 = c1 * exp(1j*th_Grey_L1);
    
    % Loop 2 Vectors (Brown starts at O4)
    R_Brown = a2 * exp(1j*theta_Brown);
    R_Blue  = b2 * exp(1j*th_Blue);
    R_LB_L2 = c2 * exp(1j*th_LightBlue);
    
    % Loop 3 Vectors (LightBlue starts at O2)
    % *** ใช้มุม th_LightBlue เดียวกัน ***
    R_LB_L3 = a3 * exp(1j*theta_LB_In); 
    R_Red   = b3 * exp(1j*th_Red);
    R_Grey_L3 = c3 * exp(1j*th_Grey_L3);
    
    % Coordinates
    O2 = 0;
    O4 = d1;
    
    % --- DRAWING ---
    
    % 1. Ground (Pink)
    plot([0 d1], [0 0], 'Color', [1 0.07 0.57], 'LineWidth', 4, 'LineStyle', '--');
    text(0, -20, 'O2'); text(d1, -20, 'O4');
    
    % 2. Loop 1 (Green-Yellow-Grey)
    quiver(0,0, real(R_Green), imag(R_Green), 0, 'g', 'LineWidth', 3, 'MaxHeadSize',0.4); % Green
    quiver(real(R_Green), imag(R_Green), real(R_Yellow), imag(R_Yellow), 0, 'y', 'LineWidth', 3, 'MaxHeadSize',0.4); % Yellow
    % quiver(d1, 0, real(R_Grey_L1), imag(R_Grey_L1), 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', ':'); % Grey (L1 Ref - ซ่อนไว้เพื่อให้รูปดูสะอาดตา)
    
    % 3. Loop 2 (Brown-Blue-LightBlue)
    % Brown (Connects to Grey L1 at O4)
    quiver(d1, 0, real(R_Brown), imag(R_Brown), 0, 'Color', [0.6 0.3 0], 'LineWidth', 3, 'MaxHeadSize',0.4); % Brown
    % Blue (Connects Brown to LightBlue)
    J_Brown = d1 + R_Brown;
    quiver(real(J_Brown), imag(J_Brown), real(R_Blue), imag(R_Blue), 0, 'b', 'LineWidth', 3, 'MaxHeadSize',0.4); % Blue
    % Light Blue (Output of L2 - ส่วนที่ยาวกว่า)
    % *** พล็อตส่วนนี้ก่อน ***
    quiver(0, 0, real(R_LB_L2), imag(R_LB_L2), 0, 'c', 'LineWidth', 3, 'MaxHeadSize',0.4); % Light Blue (Long)
    
    % 4. Loop 3 (LightBlue-Red-Grey)
    % Light Blue (Input of L3 - ส่วนที่สั้นกว่า) 
    % *** พล็อตทับลงไป มันจะซ้อนกันพอดีเพราะมุมเดียวกัน ***
    quiver(0, 0, real(R_LB_L3), imag(R_LB_L3), 0, 'c', 'LineWidth', 3, 'MaxHeadSize',0.4); % Light Blue (Short)
    
    % Red (Connects LB to Grey)
    J_LB_Short = R_LB_L3;
    quiver(real(J_LB_Short), imag(J_LB_Short), real(R_Red), imag(R_Red), 0, 'r', 'LineWidth', 3, 'MaxHeadSize',0.4); % Red
    % Grey (Output of L3)
    quiver(d1, 0, real(R_Grey_L3), imag(R_Grey_L3), 0, 'Color', [0.3 0.3 0.3], 'LineWidth', 3, 'MaxHeadSize',0.4); % Grey (Final)
    
    % Joints
    plot(real(R_Green), imag(R_Green), 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
    plot(real(J_Brown), imag(J_Brown), 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
    plot(real(J_LB_Short), imag(J_LB_Short), 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
    
    % Adjust View
    ylim([-200 350]); xlim([-100 350]); % ปรับแกนให้เห็นด้านบนมากขึ้นตามรูปวาด
end

%% --- HELPER: Solve Angle ---
function [theta, status] = solve_angle(A, B, C, mode)
    det = B^2 + A^2 - C^2;
    if det < 0, theta=0; status=0; return; end
    
    % Quadratic form logic
    Aq = C - B; Bq = 2*A; Cq = C + B;
    dq = Bq^2 - 4*Aq*Cq;
    
    if dq < 0, theta=0; status=0; return; end
    
    t1 = (-Bq + sqrt(dq))/(2*Aq);
    t2 = (-Bq - sqrt(dq))/(2*Aq);
    
    if mode == 1, theta = 2*atan(t1);
    else, theta = 2*atan(t2); end
    status = 1;
end
