clear all; close all; clc;

%% 1. SYSTEM PARAMETERS (Loop 3: Light Blue -> Blue -> Brown)
% กำหนดค่าความยาว (L) ตามที่คุณระบุ
L1 = 210; % d (Ground)
L2 = 168; % a (Light Blue - Input)
L3 = 210; % b (Blue - Coupler)
L4 = 118; % c (Brown - Output)

% Mapping ตัวแปรเข้าสูตรอาจารย์
a = L2;
b = L3;
c = L4;
d = L1;

% Input Angle (จาก Loop ก่อนหน้า)
theta2_deg = -30.69; 
q2 = deg2rad(theta2_deg); % q2 in radians

%% 2. ANALYTIC SOLVER (Lecture Pattern with K1-K5)
% คำนวณค่า K constants
K1 = d/a;
K2 = d/c;
K3 = (a^2 - b^2 + c^2 + d^2) / (2*a*c);
K4 = d/b;
K5 = (c^2 - d^2 - a^2 - b^2) / (2*a*b);

% คำนวณสัมประสิทธิ์ A, B, C สำหรับหา q4 (Output)
A = cos(q2) - K1 - K2*cos(q2) + K3;
B = -2*sin(q2);
C = K1 - (K2+1)*cos(q2) + K3;

% คำนวณสัมประสิทธิ์ D, E, F สำหรับหา q3 (Coupler)
D = cos(q2) - K1 + K4*cos(q2) + K5;
E = -2*sin(q2);
F = K1 + (K4-1)*cos(q2) + K5;

% --- Solve for q4 (Output Angles) ---
% Solution 1 & 2
q4_1 = 2*atan((-B + sqrt(B^2 - 4*A*C)) / (2*A));
q4_2 = 2*atan((-B - sqrt(B^2 - 4*A*C)) / (2*A));

% --- Solve for q3 (Coupler Angles) ---
% Solution 1 & 2
q3_1 = 2*atan((-E + sqrt(E^2 - 4*D*F)) / (2*D));
q3_2 = 2*atan((-E - sqrt(E^2 - 4*D*F)) / (2*D));

% แปลงเป็นองศา (และปรับให้เป็นค่าบวก 0-360 เพื่อความสวยงามแบบ code ตัวอย่าง)
q4_1d = mod(rad2deg(q4_1), 360);
q4_2d = mod(rad2deg(q4_2), 360);
q3_1d = mod(rad2deg(q3_1), 360);
q3_2d = mod(rad2deg(q3_2), 360);

%% 3. POSITION VECTORS CALCULATIONS
% Vector Input (a) - ใช้ร่วมกันทั้ง 2 cases
RA = a * exp(1j * q2);

% --- Case 1 (Set 1) ---
% จับคู่ q3_1 กับ q4_1 (ตรวจสอบโดยธรรมชาติของสูตรมักจะคู่กัน)
RBA1 = b * exp(1j * q3_1);      % Vector Coupler (Blue)
RBO4_1 = c * exp(1j * q4_1);    % Vector Output (Brown)
RB1 = RA + RBA1;                % ตำแหน่งจุด B จากฝั่งซ้าย

% --- Case 2 (Set 2) ---
% จับคู่ q3_2 กับ q4_2
RBA2 = b * exp(1j * q3_2);      % Vector Coupler (Blue)
RBO4_2 = c * exp(1j * q4_2);    % Vector Output (Brown)
RB2 = RA + RBA2;                % ตำแหน่งจุด B จากฝั่งซ้าย

% Ground Vector
RO4O2 = d * exp(1j * 0);

% แยก Component เพื่อ Plot (Real/Imag)
% Case 1
RAx = real(RA);         RAy = imag(RA);
RBA1x = real(RBA1);     RBA1y = imag(RBA1);
RB1x = real(RB1);       RB1y = imag(RB1);
RBO4_1x = real(RBO4_1); RBO4_1y = imag(RBO4_1);

% Case 2
RBA2x = real(RBA2);     RBA2y = imag(RBA2);
RB2x = real(RB2);       RB2y = imag(RB2);
RBO4_2x = real(RBO4_2); RBO4_2y = imag(RBO4_2);

RO4O2x = real(RO4O2);   RO4O2y = imag(RO4O2);

%% 4. PLOTTING (Side-by-Side Comparison)
figure('Color','w','Position',[50 50 1200 500]);

% --- Subplot 1: Solution Set 1 ---
subplot(1,2,1); hold on; axis equal; grid on;
title(['CASE 1: \theta_3=' num2str(q3_1d,'%.1f') '\circ, \theta_4=' num2str(q4_1d,'%.1f') '\circ']);
xlabel('X (mm)'); ylabel('Y (mm)');

% Plot Vectors (Quiver)
quiver(0, 0, RAx, RAy, 0, 'c', 'LineWidth', 3, 'MaxHeadSize', 0.5);          % Input (Light Blue/Cyan)
quiver(RAx, RAy, RBA1x, RBA1y, 0, 'b', 'LineWidth', 3, 'MaxHeadSize', 0.5);  % Coupler (Blue)
quiver(0, 0, RO4O2x, RO4O2y, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5);    % Ground
quiver(RO4O2x, RO4O2y, RBO4_1x, RBO4_1y, 0, 'Color', [0.6 0.3 0], 'LineWidth', 3, 'MaxHeadSize', 0.5); % Output (Brown)

% Joints
plot(0,0,'ko'); plot(RO4O2x,RO4O2y,'ko');
plot(RAx,RAy,'ko','MarkerFaceColor','w');
plot(RB1x,RB1y,'ko','MarkerFaceColor','w');

% --- Subplot 2: Solution Set 2 ---
subplot(1,2,2); hold on; axis equal; grid on;
title(['CASE 2: \theta_3=' num2str(q3_2d,'%.1f') '\circ, \theta_4=' num2str(q4_2d,'%.1f') '\circ']);
xlabel('X (mm)'); ylabel('Y (mm)');

% Plot Vectors (Quiver)
quiver(0, 0, RAx, RAy, 0, 'c', 'LineWidth', 3, 'MaxHeadSize', 0.5);          % Input (Light Blue/Cyan)
quiver(RAx, RAy, RBA2x, RBA2y, 0, 'b', 'LineWidth', 3, 'MaxHeadSize', 0.5);  % Coupler (Blue)
quiver(0, 0, RO4O2x, RO4O2y, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5);    % Ground
quiver(RO4O2x, RO4O2y, RBO4_2x, RBO4_2y, 0, 'Color', [0.6 0.3 0], 'LineWidth', 3, 'MaxHeadSize', 0.5); % Output (Brown)

% Joints
plot(0,0,'ko'); plot(RO4O2x,RO4O2y,'ko');
plot(RAx,RAy,'ko','MarkerFaceColor','w');
plot(RB2x,RB2y,'ko','MarkerFaceColor','w');

% Add Labels to Case 2 for clarity
text(RAx/2, RAy/2+20, 'Light Blue', 'Color','c', 'FontWeight','bold');
text(RAx + RBA2x/2, RAy + RBA2y/2 + 20, 'Blue', 'Color','b', 'FontWeight','bold');
text(RO4O2x + RBO4_2x/2, RO4O2y + RBO4_2y/2 - 20, 'Brown', 'Color',[0.6 0.3 0], 'FontWeight','bold');

% --- Output Results to Command Window ---
fprintf('========================================\n');
fprintf('RESULTS (Lecture Method K1-K5)\n');
fprintf('========================================\n');
fprintf('Set 1:\n  Theta3 (Blue) = %.2f deg\n  Theta4 (Brown) = %.2f deg\n', q3_1d, q4_1d);
fprintf('Set 2:\n  Theta3 (Blue) = %.2f deg\n  Theta4 (Brown) = %.2f deg\n', q3_2d, q4_2d);
fprintf('========================================\n');
