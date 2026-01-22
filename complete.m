clear all; close all; clc;

%% --- MAIN SETUP ---
% กำหนดค่ามุม Input หลัก (ก้านสีเหลือง Loop 1)
theta_Yellow_Input_Deg = 19.94;
theta_Yellow_Input_Rad = deg2rad(theta_Yellow_Input_Deg);

% สร้าง Figure รวม
figure('Color','w','Position',[50 50 1200 600]); 

% --- PLOT CASE 1: OPEN CONFIGURATION (Main Case) ---
subplot(1,2,1); hold on; axis equal; grid on; box on;
title('Complete Mechanism: Case 1 (All Below Ground Logic)');
xlabel('X (mm)'); ylabel('Y (mm)');
plot_full_mechanism(theta_Yellow_Input_Rad, 1); 

% --- PLOT CASE 2: CROSSED CONFIGURATION (Alternative) ---
subplot(1,2,2); hold on; axis equal; grid on; box on;
title('Complete Mechanism: Case 2');
xlabel('X (mm)'); ylabel('Y (mm)');
plot_full_mechanism(theta_Yellow_Input_Rad, 2); 


%% --- FUNCTION: Calculate and Plot Everything Together ---
function plot_full_mechanism(theta_in, config_mode)
    
    %% 1. LOOP 1 CALCULATION (Green-Yellow-Grey)
    d1 = 210; a1 = 180; b1 = 180; c1 = 118; % Ground, Green, Yellow, Grey
    
    % Solver (Yellow Input) -> Find Green & Grey
    Px = b1 * cos(theta_in) - d1;
    Py = b1 * sin(theta_in);
    P_mag_sq = Px^2 + Py^2;
    
    % สมการ A sin(t) + B cos(t) + C = 0 สำหรับ Green (a1)
    A1 = 2*a1*Py; B1 = 2*a1*Px; C1 = a1^2 + P_mag_sq - c1^2;
    
    [th_Green_1, th_Green_2, ok1] = solve_angle_all(A1, B1, C1);
    if ~ok1, text(0,0,'Loop 1 Failed'); return; end
    
    % เลือก Solution ที่ Green อยู่ใต้ Ground (y < 0)
    if sin(th_Green_1) < 0
        th_Green = th_Green_1;
    else
        th_Green = th_Green_2;
    end
    
    % หา Theta Grey (Output L1)
    th_Grey_L1 = angle(a1*exp(1j*th_Green) + b1*exp(1j*theta_in) - d1);

    %% 2. LOOP 2 CALCULATION (Brown-Blue-LightBlue)
    d2 = 210; a2 = 118; b2 = 210; c2 = 236; % Ground, Brown, Blue, LightBlue (Total=236)
    
    % Input Mapping: Brown รับมุมจาก Grey L1
    % เช็คเงื่อนไข: Brown ต้องอยู่ใต้ Ground
    % ถ้า Grey L1 ชี้ขึ้น (y>0) แต่เราต้องการ Brown ชี้ลง เราอาจต้องใช้มุมตรงข้าม หรือ Grey L1 มันชี้ลงอยู่แล้ว
    % สมมติว่าเป็นก้านเดียวกัน Rigidly connected
    theta_Brown = th_Grey_L1; 
    
    % Solver (Brown Input) -> Find Blue & LightBlue
    K1=d2/a2; K2=d2/c2; K3=(a2^2-b2^2+c2^2+d2^2)/(2*a2*c2);
    A2 = cos(theta_Brown)-K1-K2*cos(theta_Brown)+K3;
    B2 = -2*sin(theta_Brown);
    C2 = K1-(K2+1)*cos(theta_Brown)+K3;
    
    [th_LB_1, th_LB_2, ok2] = solve_angle_all(A2, B2, C2);
    if ~ok2, text(0,0,'Loop 2 Failed'); return; end
    
    % เลือก Solution ที่ Light Blue อยู่ใต้ Ground (y < 0)
    if sin(th_LB_1) < -0.01 % Tolerance
        th_LightBlue = th_LB_1;
    elseif sin(th_LB_2) < -0.01
        th_LightBlue = th_LB_2;
    else
        % ถ้าไม่ลงทั้งคู่ ให้เลือกตาม config_mode หรือตัวที่ต่ำกว่า
        if sin(th_LB_1) < sin(th_LB_2), th_LightBlue = th_LB_1; else, th_LightBlue = th_LB_2; end
    end
    
    % Calculate Coupler (Blue)
    th_Blue = angle(d2 + c2*exp(1j*th_LightBlue) - a2*exp(1j*theta_Brown));

    %% 3. LOOP 3 CALCULATION (LightBlue-Red-Grey)
    d3 = 210; a3 = 118; b3 = 210; c3 = 118; % Ground, LightBlue (Half), Red, Grey
    
    % Input Mapping: Light Blue จาก Loop 2 (มุมเดิม แต่ยาว 118)
    theta_LB_In = th_LightBlue;
    
    % Solver (LightBlue Input) -> Find Red & Grey
    K1=d3/a3; K2=d3/c3; K3=(a3^2-b3^2+c3^2+d3^2)/(2*a3*c3);
    A3 = cos(theta_LB_In)-K1-K2*cos(theta_LB_In)+K3;
    B3 = -2*sin(theta_LB_In);
    C3 = K1-(K2+1)*cos(theta_LB_In)+K3;
    
    [th_Grey_L3_1, th_Grey_L3_2, ok3] = solve_angle_all(A3, B3, C3);
    if ~ok3, text(0,0,'Loop 3 Failed'); return; end
    
    % เลือก Solution ที่ Grey อยู่ใต้ Ground (y < 0)
    % และเช็ค Red ให้อยู่ใต้ Ground ด้วย
    
    % Case A: Check th_Grey_L3_1
    th_Red_A = angle(d3 + c3*exp(1j*th_Grey_L3_1) - a3*exp(1j*theta_LB_In));
    y_Red_A = imag(a3*exp(1j*theta_LB_In) + b3*exp(1j*th_Red_A)); % ปลาย Red
    
    % Case B: Check th_Grey_L3_2
    th_Red_B = angle(d3 + c3*exp(1j*th_Grey_L3_2) - a3*exp(1j*theta_LB_In));
    y_Red_B = imag(a3*exp(1j*theta_LB_In) + b3*exp(1j*th_Red_B)); % ปลาย Red
    
    % Selection Logic
    if (sin(th_Grey_L3_1) < 0 && y_Red_A < 0)
         th_Grey_L3 = th_Grey_L3_1; th_Red = th_Red_A;
    elseif (sin(th_Grey_L3_2) < 0 && y_Red_B < 0)
         th_Grey_L3 = th_Grey_L3_2; th_Red = th_Red_B;
    else
         % ถ้าไม่มีอันไหนลงหมด ให้เลือกตาม config_mode
         if config_mode == 1
             th_Grey_L3 = th_Grey_L3_1; th_Red = th_Red_A;
         else
             th_Grey_L3 = th_Grey_L3_2; th_Red = th_Red_B;
         end
    end

    %% --- PLOTTING ALL LINKS ---
    
    % Define Vectors
    R_Ground = d1 * exp(1j*0);
    
    % Loop 1
    R_Green   = a1 * exp(1j*th_Green);
    R_Yellow  = b1 * exp(1j*theta_in);
    R_Grey_L1 = c1 * exp(1j*th_Grey_L1);
    
    % Loop 2 (Brown starts at O4)
    R_Brown = a2 * exp(1j*theta_Brown);
    R_Blue  = b2 * exp(1j*th_Blue);
    R_LB_L2 = c2 * exp(1j*th_LightBlue); % ยาว 236
    
    % Loop 3 (LightBlue starts at O2)
    R_LB_L3   = a3 * exp(1j*theta_LB_In); % ยาว 118 (ซ้อนกับ L2)
    R_Red     = b3 * exp(1j*th_Red);
    R_Grey_L3 = c3 * exp(1j*th_Grey_L3);
    
    % --- DRAWING ---
    
    % 1. Ground (Pink)
    plot([0 d1], [0 0], 'Color', [1 0.07 0.57], 'LineWidth', 4, 'LineStyle', '--');
    text(0, -20, 'O2'); text(d1, -20, 'O4');
    
    % 2. Loop 1 (Green-Yellow-Grey)
    quiver(0,0, real(R_Green), imag(R_Green), 0, 'g', 'LineWidth', 3, 'MaxHeadSize',0.4, 'DisplayName','Green'); 
    quiver(real(R_Green), imag(R_Green), real(R_Yellow), imag(R_Yellow), 0, 'y', 'LineWidth', 3, 'MaxHeadSize',0.4, 'DisplayName','Yellow'); 
    % Grey L1 (ซ่อนหรือแสดงเป็น Ghost)
    % quiver(d1, 0, real(R_Grey_L1), imag(R_Grey_L1), 0, 'k', 'LineWidth', 1, 'LineStyle', ':'); 
    
    % 3. Loop 2 (Brown-Blue-LightBlue)
    % Brown (Connects to Grey L1 at O4 -> Rigid/Same link?)
    quiver(d1, 0, real(R_Brown), imag(R_Brown), 0, 'Color', [0.6 0.3 0], 'LineWidth', 3, 'MaxHeadSize',0.4, 'DisplayName','Brown'); 
    
    % Blue (Connects Brown to LightBlue Tip)
    J_Brown = d1 + R_Brown;
    quiver(real(J_Brown), imag(J_Brown), real(R_Blue), imag(R_Blue), 0, 'b', 'LineWidth', 3, 'MaxHeadSize',0.4, 'DisplayName','Blue'); 
    
    % Light Blue (Output of L2 - ยาว 236)
    quiver(0, 0, real(R_LB_L2), imag(R_LB_L2), 0, 'c', 'LineWidth', 3, 'MaxHeadSize',0.4, 'DisplayName','Light Blue (236)'); 
    
    % 4. Loop 3 (LightBlue-Red-Grey)
    % Light Blue (Input of L3 - ยาว 118 - ซ้อนทับ)
    % ไม่ต้องวาดซ้ำก็ได้ หรือวาดทับ
    
    % Red (Connects LB Midpoint to Grey)
    J_LB_Mid = R_LB_L3;
    quiver(real(J_LB_Mid), imag(J_LB_Mid), real(R_Red), imag(R_Red), 0, 'r', 'LineWidth', 3, 'MaxHeadSize',0.4, 'DisplayName','Red'); 
    
    % Grey (Output of L3)
    quiver(d1, 0, real(R_Grey_L3), imag(R_Grey_L3), 0, 'Color', [0.3 0.3 0.3], 'LineWidth', 3, 'MaxHeadSize',0.4, 'DisplayName','Grey'); 
    
    % Joints
    plot(real(R_Green), imag(R_Green), 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
    plot(real(J_Brown), imag(J_Brown), 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
    plot(real(J_LB_Mid), imag(J_LB_Mid), 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
    plot(real(J_Brown+R_Blue), imag(J_Brown+R_Blue), 'ko', 'MarkerFaceColor','w', 'MarkerSize',6); % ปลาย Blue
    
    % Adjust View
    ylim([-350 250]); xlim([-100 400]);
    legend('Location','NorthEast');
end

%% --- HELPER: Solve Quadratic for Angle ---
function [th1, th2, status] = solve_angle_all(A, B, C)
    det = B^2 + A^2 - C^2;
    if det < 0, th1=0; th2=0; status=0; return; end
    
    A_q = C - B; B_q = 2 * A; C_q = C + B;
    det_q = B_q^2 - 4*A_q*C_q;
    
    if det_q < 0, th1=0; th2=0; status=0; return; end
    
    t1 = (-B_q + sqrt(det_q)) / (2*A_q);
    t2 = (-B_q - sqrt(det_q)) / (2*A_q);
    
    th1 = 2*atan(t1);
    th2 = 2*atan(t2);
    status = 1;
end
