clear all; close all; clc;

%% 1. SYSTEM PARAMETERS (Loop 2)
% กำหนดค่าความยาวตามที่คุณระบุ
L1 = 210; % d (Ground)
L2 = 168; % a (Light Blue Link - Input)
L3 = 210; % b (Blue Link - Coupler)
L4 = 118; % c (Brown Link - Output)

a = L2;
b = L3;
c = L4;
d = L1;

% Input Angle (จาก Loop ที่แล้ว)
q2d = -30.69; 
q2 = deg2rad(q2d); % radian angle

%% 2. CALCULATION (Professor's Pattern with K1-K5)
% คำนวณค่า K ทั้ง 5 ตัว
K1 = d/a;
K2 = d/c;
K3 = (a^2 - b^2 + c^2 + d^2) / (2*a*c);
K4 = d/b;
K5 = (c^2 - d^2 - a^2 - b^2) / (2*a*b);

% คำนวณสัมประสิทธิ์ A-F
A = cos(q2) - K1 - K2*cos(q2) + K3;
B = -2*sin(q2);
C = K1 - (K2+1)*cos(q2) + K3;

D = cos(q2) - K1 + K4*cos(q2) + K5;
E = -2*sin(q2);
F = K1 + (K4-1)*cos(q2) + K5;

% คำนวณมุม q4 (Output Link)
q41 = 2*atan((-B + sqrt(B^2 - 4*A*C)) / (2*A)); % Solution 1
q42 = 2*atan((-B - sqrt(B^2 - 4*A*C)) / (2*A)); % Solution 2

% คำนวณมุม q3 (Coupler Link)
q31 = 2*atan((-E + sqrt(E^2 - 4*D*F)) / (2*D)); % Solution 1
q32 = 2*atan((-E - sqrt(E^2 - 4*D*F)) / (2*D)); % Solution 2

% แปลงเป็นองศาเพื่อดูผลลัพธ์
q41d = rad2deg(q41);
q42d = rad2deg(q42);
q31d = rad2deg(q31);
q32d = rad2deg(q32);

fprintf('--- Results ---\n');
fprintf('Input q2 = %.2f deg\n', q2d);
fprintf('Set 1: q3 = %.2f, q4 = %.2f\n', q31d, q41d);
fprintf('Set 2: q3 = %.2f, q4 = %.2f\n', q32d, q42d);

%% 3. VECTOR ANALYSIS & PLOTTING
% สร้าง Vector ตามแพทเทิร์นอาจารย์
RA = a*exp(1j*q2);

% --- Case 1: (Crossed?) ---
% หมายเหตุ: ต้องจับคู่ q3 กับ q4 ให้ถูกคู่ (ปกติ 1 คู่ 1, 2 คู่ 2 แต่อาจสลับกันได้ ต้องเช็ค)
% ลองจับคู่ Set 1 (q31, q41)
RBA1 = b*exp(1j*q31); 
RB1_Left = RA + RBA1;     % หาจุด B จากทางซ้าย (Input -> Coupler)

RO4O2 = d*exp(1j*0);
RBO41 = c*exp(1j*q41);
RB1_Right = RO4O2 + RBO41; % หาจุด B จากทางขวา (Ground -> Output)

% --- Case 2: (Open?) ---
% ลองจับคู่ Set 2 (q32, q42)
RBA2 = b*exp(1j*q32); 
RB2_Left = RA + RBA2;

RBO42 = c*exp(1j*q42);
RB2_Right = RO4O2 + RBO42;

% แยก Component เพื่อ Plot (ใช้ Case 2: Open เป็นตัวอย่างหลักตามโค้ดอาจารย์)
% *คุณเลือกเปลี่ยนตัวแปรตรงนี้ได้ว่าจะเอา set 1 หรือ 2 มาพล็อต*
% สมมติเลือก Set 2 (Open) มาโชว์
RAx = real(RA); RAy = imag(RA);
RBA2x = real(RBA2); RBA2y = imag(RBA2);
RB2x = real(RB2_Left); RB2y = imag(RB2_Left);

RO4O2x = real(RO4O2); RO4O2y = imag(RO4O2);
RBO42x = real(RBO42); RBO42y = imag(RBO42);

%% 4. PLOT (Quiver Style)
figure('Color','w'); hold on; axis equal; grid on;
title(['Loop 2 Analysis: Input \theta_2 = ' num2str(q2d) '^\circ']);

% วาด Input (Red)
quiver(0, 0, RAx, RAy, 0, 'r', 'MaxHeadSize', 0.5, 'LineWidth', 2, 'DisplayName', 'Input (a)');

% วาด Coupler (Blue) - ต่อจากหัวลูกศร Input
quiver(RAx, RAy, RBA2x, RBA2y, 0, 'b', 'MaxHeadSize', 0.5, 'LineWidth', 2, 'DisplayName', 'Coupler (b)');

% วาด Vector ผลลัพธ์จุด B (Green) - ลากจากจุดกำเนิดไปจุด B
quiver(0, 0, RB2x, RB2y, 0, 'g', 'MaxHeadSize', 0.5, 'LineWidth', 2, 'DisplayName', 'Position B');

% วาด Ground (Black)
quiver(0, 0, RO4O2x, RO4O2y, 0, 'k', 'MaxHeadSize', 0.5, 'LineWidth', 2, 'DisplayName', 'Ground (d)');

% วาด Output (Black/Brown) - ลากจาก O4 ไปจุด B
quiver(RO4O2x, RO4O2y, RBO42x, RBO42y, 0, 'Color', [0.6 0.4 0.2], 'MaxHeadSize', 0.5, 'LineWidth', 2, 'DisplayName', 'Output (c)');

legend('Location', 'best');
xlabel('X'); ylabel('Y');

% เช็คความถูกต้อง (Vector Loop Closure Error)
Error_Dist = abs(RB2_Left - RB2_Right);
fprintf('Loop Closure Error: %.4e\n', Error_Dist);
if Error_Dist > 1e-4
    warning('Warning: The loop does not close properly. Try swapping q3/q4 pairs.');
end
