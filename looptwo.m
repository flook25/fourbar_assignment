clear all; close all; clc;

%% 1. SYSTEM PARAMETERS (Loop 2)
% ตามข้อมูลที่คุณให้มา:
L1 = 210; % d (Ground)
L2 = 168; % a (Light Blue Link - Input)
L3 = 210; % b (Blue Link - Coupler)
L4 = 118; % c (Brown Link - Output)

% กำหนดตัวแปรสำหรับเข้าสูตร
d = L1;
a = L2;
b = L3;
c = L4;

% --- INPUT FROM PREVIOUS STEP ---
theta2_deg = 90;      % Input Angle (Green/Light Blue)
theta2 = deg2rad(theta2_deg); 

%% 2. SOLVER (Forward Kinematics: Known Theta2 -> Find Theta3, Theta4)
% คำนวณค่าคงที่ K
K1 = d/a;
K2 = d/c;
K3 = (a^2 - b^2 + c^2 + d^2) / (2*a*c);

% คำนวณสัมประสิทธิ์ A, B, C
A_eq = cos(theta2) - K1 - K2*cos(theta2) + K3;
B_eq = -2 * sin(theta2);
C_eq = K1 - (K2+1)*cos(theta2) + K3;

% ตรวจสอบ Discriminant
det = B_eq^2 - 4*A_eq*C_eq;
if det < 0
    error('Assembly failed: No solution for this input angle.');
end

% --- คำนวณหา Theta4 (Output - Brown) ---
% Case 1: Open Circuit (มักใช้เครื่องหมาย + หรือ - ขึ้นอยู่กับ configuration)
t4_sol1 = 2 * atan2( (-B_eq + sqrt(det)), (2*A_eq) );
% Case 2: Crossed Circuit
t4_sol2 = 2 * atan2( (-B_eq - sqrt(det)), (2*A_eq) );

% --- คำนวณหา Theta3 (Coupler - Blue) ---
% ฟังก์ชันหา Theta3 จาก Vector Loop
get_theta3 = @(t4) atan2( (c*sin(t4) - a*sin(theta2)), (c*cos(t4) + d - a*cos(theta2)) );

t3_sol1 = get_theta3(t4_sol1);
t3_sol2 = get_theta3(t4_sol2);

% --- เลือก Case ที่ต้องการแสดงผล ---
% ลองเช็คดูว่า Case ไหนให้มุม Blue ใกล้ 90 หรือตามรูปโจทย์
% ปกติถ้าแขนไขว้กันจะเป็น Crossed (Sol 2)
% เดี๋ยวจะ Plot ทั้ง 2 แบบให้เลือกครับ

%% 3. PLOTTING FUNCTION
figure('Color','w','Position',[100 100 1000 500]); 

% --- PLOT CASE 1 (Open) ---
subplot(1,2,1); hold on; axis equal; grid on;
title('Case 1: Open Configuration');
plot_mechanism(a, b, c, d, theta2, t3_sol1, t4_sol1);
xlim([-100 350]); ylim([-150 250]);

% --- PLOT CASE 2 (Crossed) ---
subplot(1,2,2); hold on; axis equal; grid on;
title('Case 2: Crossed Configuration');
plot_mechanism(a, b, c, d, theta2, t3_sol2, t4_sol2);
xlim([-100 350]); ylim([-150 250]);

% แสดงผลลัพธ์เป็นตัวเลข
fprintf('========================================\n');
fprintf('RESULTS FOR LOOP 2 (Input = %.2f deg)\n', theta2_deg);
fprintf('========================================\n');
fprintf('CASE 1 (Open):\n');
fprintf('  Theta 3 (Blue) : %8.2f deg\n', rad2deg(t3_sol1));
fprintf('  Theta 4 (Brown): %8.2f deg\n', rad2deg(t4_sol1));
fprintf('----------------------------------------\n');
fprintf('CASE 2 (Crossed):\n');
fprintf('  Theta 3 (Blue) : %8.2f deg\n', rad2deg(t3_sol2));
fprintf('  Theta 4 (Brown): %8.2f deg\n', rad2deg(t4_sol2));
fprintf('========================================\n');


%% --- INTERNAL FUNCTION FOR PLOTTING ---
function plot_mechanism(a, b, c, d, th2, th3, th4)
    % Position Vectors
    R_Input  = a * exp(1j * th2);    % Light Blue
    R_Coupler= b * exp(1j * th3);    % Blue
    R_Output = c * exp(1j * th4);    % Brown
    R_Ground = d * exp(1j * 0);      % Ground

    % Coordinates
    O2 = 0;
    O4 = d;
    J1 = R_Input;             % จุดต่อ Input-Coupler
    J2 = R_Input + R_Coupler; % จุดต่อ Coupler-Output

    % Plot Vectors (Quiver Style)
    % Ground (Black)
    plot([0 real(R_Ground)], [0 imag(R_Ground)], 'k--', 'LineWidth', 2);
    text(0, -20, 'O2'); text(real(R_Ground), -20, 'O4');

    % Input Link (Light Blue - Cyan)
    quiver(0, 0, real(R_Input), imag(R_Input), 0, 'c', 'LineWidth', 4, 'MaxHeadSize', 0.5);
    
    % Coupler Link (Blue - Dark Blue)
    quiver(real(J1), imag(J1), real(R_Coupler), imag(R_Coupler), 0, 'b', 'LineWidth', 4, 'MaxHeadSize', 0.5);
    
    % Output Link (Brown) - Plot from O4
    % Vector Loop: R_Output goes from O4 to J2
    quiver(real(R_Ground), imag(R_Ground), real(R_Output), imag(R_Output), 0, 'Color', [0.6 0.3 0], 'LineWidth', 4, 'MaxHeadSize', 0.5);

    % Joints
    plot(real(J1), imag(J1), 'ko', 'MarkerFaceColor', 'w');
    plot(real(J2), imag(J2), 'ko', 'MarkerFaceColor', 'w');
end
