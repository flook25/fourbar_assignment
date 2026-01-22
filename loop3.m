clear all; close all; clc;

%% --- STEP 1: SOLVE LOOP 2 (To get Light Blue Angle) ---
% Parameters Loop 2
L_Ground_L2    = 210; 
L_Brown_Input  = 118; % Input (a)
L_Blue_Coupler = 210; % Coupler (b)
L_LightBlue_Out= 168; % Output (c) - Length for Loop 2

% Input Angle for Loop 2 (Brown Link at O4)
theta_brown_deg = -81.65;
theta2_L2 = deg2rad(theta_brown_deg);

% Norton Parameters for Loop 2
d = L_Ground_L2; a = L_Brown_Input; b = L_Blue_Coupler; c = L_LightBlue_Out;
K1 = d/a; K2 = d/c; K3 = (a^2 - b^2 + c^2 + d^2)/(2*a*c);

A = cos(theta2_L2) - K1 - K2*cos(theta2_L2) + K3;
B = -2*sin(theta2_L2);
C = K1 - (K2+1)*cos(theta2_L2) + K3;

det = B^2 - 4*A*C;
if det < 0, error('Loop 2 Assembly failed'); end

% Calculate Theta Output (Light Blue) - Using Case 2 (Crossed/Open depending on config)
% จากผลลัพธ์ก่อนหน้า เราใช้ t4_2
t4_2 = 2*atan2(-B - sqrt(det), 2*A);
theta_LightBlue = t4_2; % นี่คือมุม Input ที่จะส่งไป Loop 3

%% --- STEP 2: SOLVE LOOP 3 (Using Light Blue as Input) ---
fprintf('----------------------------------------\n');
fprintf('STARTING LOOP 3 CALCULATION\n');
fprintf('----------------------------------------\n');

% 1. SYSTEM PARAMETERS (Loop 3) - ใช้ข้อมูลใหม่ที่คุณให้มา
L_Ground_L3 = 210; % d
L_LightBlue_In = 118; % a (Link 2 at O2) - Note: Length is 118 here
L_Red_Coupler  = 210; % b (Link 3)
L_Grey_Out     = 118; % c (Link 4 at O4)

% กำหนดตัวแปรสำหรับสูตร Norton Loop 3
d = L_Ground_L3;
a = L_LightBlue_In; % Input (Light Blue)
b = L_Red_Coupler;  % Coupler (Red)
c = L_Grey_Out;     % Output (Grey)

% *** CRITICAL: INPUT THETA ***
% ใช้มุม Light Blue ที่ได้จาก Loop 2
theta2_L3 = theta_LightBlue; 
fprintf('Input Theta for Loop 3 (Light Blue): %.2f deg\n', rad2deg(theta2_L3));

% 2. SOLVER (Forward Kinematics for Loop 3)
K1 = d/a;
K2 = d/c;
K3 = (a^2 - b^2 + c^2 + d^2)/(2*a*c);
K4 = d/b;
K5 = (c^2 - d^2 - a^2 - b^2)/(2*a*b);

% Find Theta 4 (Output - Grey)
A_eq = cos(theta2_L3) - K1 - K2*cos(theta2_L3) + K3;
B_eq = -2*sin(theta2_L3);
C_eq = K1 - (K2+1)*cos(theta2_L3) + K3;

det_L3 = B_eq^2 - 4*A_eq*C_eq;
if det_L3 < 0, error('Loop 3 Assembly failed'); end

% Calculate Theta 4 (Grey) - 2 Solutions
t4_sol1 = 2*atan2(-B_eq + sqrt(det_L3), 2*A_eq);
t4_sol2 = 2*atan2(-B_eq - sqrt(det_L3), 2*A_eq);

% Find Theta 3 (Coupler - Red)
D_eq = cos(theta2_L3) - K1 + K4*cos(theta2_L3) + K5;
E_eq = -2*sin(theta2_L3);
F_eq = K1 + (K4-1)*cos(theta2_L3) + K5;

det3_L3 = E_eq^2 - 4*D_eq*F_eq;
t3_sol1 = 2*atan2(-E_eq + sqrt(det3_L3), 2*D_eq);
t3_sol2 = 2*atan2(-E_eq - sqrt(det3_L3), 2*D_eq);

% เลือก Case ที่ต้องการแสดง (เลือก Case 1 สำหรับ Loop นี้ หรือลองสลับดูรูป)
theta_Grey = t4_sol1;
theta_Red  = t3_sol1;

%% 3. PLOTTING LOOP 3
% คำนวณพิกัดเวกเตอร์
R_Ground    = d * exp(1j * 0);
R_LightBlue = a * exp(1j * theta2_L3); % Input (Light Blue @ 118mm)
R_Red       = b * exp(1j * theta_Red);
R_Grey      = c * exp(1j * theta_Grey);

% Coordinates
O2 = 0;
O4 = R_Ground;
J_LB_Red = O2 + R_LightBlue;   % จุดต่อ LightBlue-Red
J_Red_Grey = O4 + R_Grey;      % จุดต่อ Red-Grey

figure('Color','w','Position',[100 100 800 600]); 
hold on; axis equal; grid on; box on;
title(['Loop 3 Analysis: Light Blue Input \theta = ' num2str(rad2deg(theta2_L3), '%.2f') '^\circ']);
xlabel('X (mm)'); ylabel('Y (mm)');

% 1. Ground (Pink)
plot([0 real(R_Ground)], [0 imag(R_Ground)], 'Color', [1 0.07 0.57], 'LineWidth', 3, 'LineStyle', '--');
text(0, -20, 'O2'); text(real(O4), -20, 'O4');

% 2. Light Blue Link (O2 -> Joint)
quiver(0, 0, real(R_LightBlue), imag(R_LightBlue), 0, 'c', 'LineWidth', 5, 'MaxHeadSize', 0.5, 'DisplayName', 'Light Blue (118)');

% 3. Red Link (Coupler)
% วาดจากปลาย Light Blue ไปหาปลาย Grey
quiver(real(J_LB_Red), imag(J_LB_Red), real(R_Red), imag(R_Red), 0, 'r', 'LineWidth', 5, 'MaxHeadSize', 0.5, 'DisplayName', 'Red');

% 4. Grey Link (O4 -> Joint)
% วาดจาก O4 ขึ้นไป
quiver(real(O4), imag(O4), real(R_Grey), imag(R_Grey), 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 5, 'MaxHeadSize', 0.5, 'DisplayName', 'Grey');

% Joints
plot(real(J_LB_Red), imag(J_LB_Red), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
plot(real(J_Red_Grey), imag(J_Red_Grey), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);

legend('Location','NorthEast');
xlim([-50 350]); ylim([-250 250]);

% แสดงผลลัพธ์
fprintf('RESULTS FOR LOOP 3\n');
fprintf('========================================\n');
fprintf('Theta Light Blue (Input): %8.2f deg\n', rad2deg(theta2_L3));
fprintf('Theta Red (Coupler)     : %8.2f deg\n', rad2deg(theta_Red));
fprintf('Theta Grey (Output)     : %8.2f deg\n', rad2deg(theta_Grey));
fprintf('========================================\n');
