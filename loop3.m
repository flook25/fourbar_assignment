clear all; close all; clc;

%% --- STEP 1: PRE-CALCULATE INPUT FROM LOOP 2 ---
% คำนวณหา Theta Light Blue ที่ถูกต้องจาก Loop 2 (กรณี Brown ชี้ลง)
L_Ground_L2 = 210; L_Brown = 118; L_Blue = 210; L_LightBlue_L2 = 168;
theta_Brown = deg2rad(-81.65); % Brown Input

% Solver Loop 2 (Simplified for getting input)
d = L_Ground_L2; a = L_Brown; b = L_Blue; c = L_LightBlue_L2;
K1=d/a; K2=d/c; K3=(a^2-b^2+c^2+d^2)/(2*a*c);
A=cos(theta_Brown)-K1-K2*cos(theta_Brown)+K3; B=-2*sin(theta_Brown); C=K1-(K2+1)*cos(theta_Brown)+K3;
det = B^2 - 4*A*C;
% เลือกคำตอบที่ทำให้ Light Blue ชี้ลง (Below Ground)
t4_1 = 2*atan2(-B+sqrt(det), 2*A);
t4_2 = 2*atan2(-B-sqrt(det), 2*A);
if sin(t4_1) < 0, theta_LB_In = t4_1; else, theta_LB_In = t4_2; end

fprintf('Input Theta from Loop 2 (Light Blue): %.2f deg\n', rad2deg(theta_LB_In));

%% --- STEP 2: SOLVE LOOP 3 (Target) ---
% SYSTEM PARAMETERS (Loop 3)
L_Ground    = 210; % d
L_LightBlue = 118; % a (Input for L3, Length is 118 here)
L_Red       = 210; % b (Coupler)
L_Grey      = 118; % c (Output)

% Map to Norton Variables
d = L_Ground;
a = L_LightBlue; % Input
b = L_Red;       % Coupler
c = L_Grey;      % Output

% Input Angle (From Step 1)
theta2 = theta_LB_In; 

% 1. Calculate Constants
K1 = d/a;
K2 = d/c;
K3 = (a^2 - b^2 + c^2 + d^2)/(2*a*c);
K4 = d/b;
K5 = (c^2 - d^2 - a^2 - b^2)/(2*a*b);

% 2. Calculate Theta 4 (Grey - Output)
A = cos(theta2) - K1 - K2*cos(theta2) + K3;
B = -2*sin(theta2);
C = K1 - (K2+1)*cos(theta2) + K3;

det = B^2 - 4*A*C;
if det < 0, error('Assembly failed'); end

t4_1 = 2*atan2(-B + sqrt(det), 2*A);
t4_2 = 2*atan2(-B - sqrt(det), 2*A);

% 3. Calculate Theta 3 (Red - Coupler)
D = cos(theta2) - K1 + K4*cos(theta2) + K5;
E = -2*sin(theta2);
F = K1 + (K4-1)*cos(theta2) + K5;

det3 = E^2 - 4*D*F;
t3_1 = 2*atan2(-E + sqrt(det3), 2*D);
t3_2 = 2*atan2(-E - sqrt(det3), 2*D);

% --- SELECTION LOGIC: ALL LINKS BELOW GROUND ---
% เราต้องการให้จุดต่อ (Joints) มีค่า Y < 0
% Joint 1 (LB-Red): y = a*sin(theta2) (รู้อยู่แล้วว่าเป็นลบ)
% Joint 2 (Red-Grey): y = c*sin(theta4) (ต้องเช็คค่านี้)

% Check Case 1
y_grey_1 = c * sin(t4_1);
% Check Case 2
y_grey_2 = c * sin(t4_2);

if y_grey_1 < 0
    theta_Grey = t4_1;
    % คู่กับ Red ที่เหมาะสม (ปกติ Parallelogram จะคู่กับ Red ~ 0 หรือ 180)
    % เช็ค Vector Loop: R_LB + R_Red - R_Grey - R_Ground = 0
    % ลองคู่กับ t3_1
    err1 = abs(a*exp(1j*theta2) + b*exp(1j*t3_1) - c*exp(1j*theta_Grey) - d);
    if err1 < 1e-4
        theta_Red = t3_1;
    else
        theta_Red = t3_2;
    end
else
    theta_Grey = t4_2;
    err1 = abs(a*exp(1j*theta2) + b*exp(1j*t3_1) - c*exp(1j*theta_Grey) - d);
    if err1 < 1e-4
        theta_Red = t3_1;
    else
        theta_Red = t3_2;
    end
end

%% 3. PLOTTING
R_Ground    = d * exp(1j * 0);
R_LightBlue = a * exp(1j * theta2);
R_Red       = b * exp(1j * theta_Red);
R_Grey      = c * exp(1j * theta_Grey);

% Coordinates
O2 = 0;
O4 = R_Ground;
J_LB_Red = O2 + R_LightBlue;   % Joint 1
J_Red_Grey = O4 + R_Grey;      % Joint 2 (Note: Grey vector starts at O4 in Norton derivation)

figure('Color','w','Position',[100 100 800 600]); 
hold on; axis equal; grid on; box on;
title(['Loop 3: All Links Below Ground (Input \theta_{LB} = ' num2str(rad2deg(theta2), '%.2f') '^\circ)']);
xlabel('X (mm)'); ylabel('Y (mm)');

% 1. Ground (Pink)
plot([0 real(R_Ground)], [0 imag(R_Ground)], 'Color', [1 0.07 0.57], 'LineWidth', 3, 'LineStyle', '--');
text(0, 15, 'O2'); text(real(O4), 15, 'O4');

% 2. Light Blue Link (O2 -> Joint 1)
quiver(0, 0, real(R_LightBlue), imag(R_LightBlue), 0, 'c', 'LineWidth', 5, 'MaxHeadSize', 0.5, 'DisplayName', 'Light Blue');

% 3. Red Link (Joint 1 -> Joint 2)
% สร้าง Vector จากจุด J_LB_Red ไป J_Red_Grey
vec_Red = J_Red_Grey - J_LB_Red;
quiver(real(J_LB_Red), imag(J_LB_Red), real(vec_Red), imag(vec_Red), 0, 'r', 'LineWidth', 5, 'MaxHeadSize', 0.5, 'DisplayName', 'Red');

% 4. Grey Link (O4 -> Joint 2)
quiver(real(O4), imag(O4), real(R_Grey), imag(R_Grey), 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 5, 'MaxHeadSize', 0.5, 'DisplayName', 'Grey');

% Joints
plot(real(J_LB_Red), imag(J_LB_Red), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
plot(real(J_Red_Grey), imag(J_Red_Grey), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);

legend('Location','NorthEast');
xlim([-50 300]); ylim([-250 100]); % ปรับแกน Y ให้เน้นด้านล่าง

% แสดงผลลัพธ์
fprintf('========================================\n');
fprintf('RESULTS LOOP 3 (All Below Ground)\n');
fprintf('========================================\n');
fprintf('Theta Light Blue (Input): %8.2f deg\n', rad2deg(theta2));
fprintf('Theta Red (Coupler)     : %8.2f deg\n', rad2deg(theta_Red));
fprintf('Theta Grey (Output)     : %8.2f deg\n', rad2deg(theta_Grey));
fprintf('========================================\n');
