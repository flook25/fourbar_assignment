clear all; close all; clc;

%% 1. SYSTEM PARAMETERS (Loop 3: Light Blue -> Blue -> Brown)
L1 = 210; % d (Ground)
L2 = 168; % a (Light Blue Link - Input)
L3 = 210; % b (Blue Link - Coupler)
L4 = 118; % c (Brown Link - Output)

a = L2;
b = L3;
c = L4;
d = L1;

% Input Angle (จาก Loop ก่อนหน้า)
theta2_deg = -30.69;
q2 = deg2rad(theta2_deg); % q2 in radians

%% 2. CALCULATION (Professor's Pattern with K1-K5)
% คำนวณค่า K constants
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

% คำนวณมุม q4 (Output) - มี 2 คำตอบ
q41 = 2*atan((-B + sqrt(B^2 - 4*A*C)) / (2*A));
q42 = 2*atan((-B - sqrt(B^2 - 4*A*C)) / (2*A));

% คำนวณมุม q3 (Coupler) - มี 2 คำตอบ
q31 = 2*atan((-E + sqrt(E^2 - 4*D*F)) / (2*D));
q32 = 2*atan((-E - sqrt(E^2 - 4*D*F)) / (2*D));

% แปลงเป็น Degree เพื่อดูผลลัพธ์
q41d = rad2deg(q41);
q42d = rad2deg(q42);
q31d = rad2deg(q31);
q32d = rad2deg(q32);

fprintf('q2 input = %.2f\n', theta2_deg);
fprintf('Solution 1: q3 = %.2f, q4 = %.2f\n', q31d, q41d);
fprintf('Solution 2: q3 = %.2f, q4 = %.2f\n', q32d, q42d);

%% 3. VECTOR DEFINITION (Exp style)
RA = a*exp(1j*q2);

% --- Set 1 (Crossed Circuit?) ---
RBA1 = b*exp(1j*q31);
RBO41 = c*exp(1j*q41);
RB1 = RA + RBA1; % Resultant Vector B from Path 1 (Input + Coupler)

% --- Set 2 (Open Circuit?) ---
RBA2 = b*exp(1j*q32);
RBO42 = c*exp(1j*q42);
RB2 = RA + RBA2; % Resultant Vector B from Path 1 (Input + Coupler)

RO4O2 = d*exp(1j*0); % Ground Vector

% แยก Component เพื่อ Plot
RAx = real(RA); RAy = imag(RA);
RO4O2x = real(RO4O2); RO4O2y = imag(RO4O2);

% Set 1 Components
RBA1x = real(RBA1); RBA1y = imag(RBA1);
RBO41x = real(RBO41); RBO41y = imag(RBO41);
RB1x = real(RB1); RB1y = imag(RB1);

% Set 2 Components
RBA2x = real(RBA2); RBA2y = imag(RBA2);
RBO42x = real(RBO42); RBO42y = imag(RBO42);
RB2x = real(RB2); RB2y = imag(RB2);

%% 4. PLOTTING (Quiver Style like Professor)
figure('Color','w','Position',[100 100 1000 500]);

% --- Plot Solution 1 ---
subplot(1,2,1); hold on; axis equal; grid on;
title(['Solution 1 (q3=' num2str(q31d,'%.1f') ')']);
xlabel('X'); ylabel('Y');

% 1. Input Vector (Red) - RA
quiver(0, 0, RAx, RAy, 0, 'r', 'MaxHeadSize', 0.5, 'LineWidth', 2);

% 2. Coupler Vector (Blue) - RBA (ต่อจาก RA)
quiver(RAx, RAy, RBA1x, RBA1y, 0, 'b', 'MaxHeadSize', 0.5, 'LineWidth', 2);

% 3. Resultant Vector B (Green) - RB (จากจุดกำเนิดไปหาจุดต่อ)
quiver(0, 0, RB1x, RB1y, 0, 'g', 'MaxHeadSize', 0.5, 'LineWidth', 2);

% 4. Ground Vector (Black) - RO4O2
quiver(0, 0, RO4O2x, RO4O2y, 0, 'k', 'MaxHeadSize', 0.5, 'LineWidth', 2);

% 5. Output Vector (Black/Grey) - RBO4 (ต่อจาก Ground)
quiver(RO4O2x, RO4O2y, RBO41x, RBO41y, 0, 'k', 'MaxHeadSize', 0.5, 'LineWidth', 2);


% --- Plot Solution 2 ---
subplot(1,2,2); hold on; axis equal; grid on;
title(['Solution 2 (q3=' num2str(q32d,'%.1f') ')']);
xlabel('X'); ylabel('Y');

% 1. Input Vector (Red) - RA
quiver(0, 0, RAx, RAy, 0, 'r', 'MaxHeadSize', 0.5, 'LineWidth', 2);

% 2. Coupler Vector (Blue) - RBA (ต่อจาก RA)
quiver(RAx, RAy, RBA2x, RBA2y, 0, 'b', 'MaxHeadSize', 0.5, 'LineWidth', 2);

% 3. Resultant Vector B (Green) - RB (จากจุดกำเนิดไปหาจุดต่อ)
quiver(0, 0, RB2x, RB2y, 0, 'g', 'MaxHeadSize', 0.5, 'LineWidth', 2);

% 4. Ground Vector (Black) - RO4O2
quiver(0, 0, RO4O2x, RO4O2y, 0, 'k', 'MaxHeadSize', 0.5, 'LineWidth', 2);

% 5. Output Vector (Black/Grey) - RBO4 (ต่อจาก Ground)
quiver(RO4O2x, RO4O2y, RBO42x, RBO42y, 0, 'k', 'MaxHeadSize', 0.5, 'LineWidth', 2);
