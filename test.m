clear all; close all; clc;

%% 1. SYSTEM PARAMETERS (Loop 1: Green -> Yellow -> Grey)
L1 = 210; % d (Ground)
L2 = 180; % a (Green Link - Crank)
L3 = 180; % b (Yellow Link - Coupler/Input Known)
L4 = 118; % c (Grey Link - Output)

% Parameter mapping
a = L2;
b = L3;
c = L4;
d = L1;

% โจทย์กำหนดมุม Link 3 (Yellow) มาให้
theta3_deg = 19.94;
q3 = deg2rad(theta3_deg); 

%% 2. SOLVER (Given q3, Find q2 and q4)
% เราจะใช้สูตรที่ประยุกต์มาจาก Vector Loop เพื่อหา q2 ก่อน
% สมการ: P*cos(q2) + Q*sin(q2) + R = 0

Kx = b*cos(q3) - d;
Ky = b*sin(q3);

P = 2*a*Kx;
Q = 2*a*Ky;
R = a^2 + Kx^2 + Ky^2 - c^2;

% แก้สมการ Quadratic หา t = tan(q2/2)
% (R - P)t^2 + 2Qt + (R + P) = 0
A_quad = R - P;
B_quad = 2*Q;
C_quad = R + P;

det = B_quad^2 - 4*A_quad*C_quad;
if det < 0
    error('No solution');
end

% หา q2 (Green) - มี 2 คำตอบ
t1 = (-B_quad + sqrt(det)) / (2*A_quad);
t2 = (-B_quad - sqrt(det)) / (2*A_quad);

q2_1 = 2*atan(t1);
q2_2 = 2*atan(t2);

% หา q4 (Grey) จาก q2 ที่ได้
% c*cos(q4) = a*cos(q2) + Kx
% c*sin(q4) = a*sin(q2) + Ky

% Set 1
val_cos1 = (a*cos(q2_1) + Kx) / c;
val_sin1 = (a*sin(q2_1) + Ky) / c;
q4_1 = atan2(val_sin1, val_cos1);

% Set 2
val_cos2 = (a*cos(q2_2) + Kx) / c;
val_sin2 = (a*sin(q2_2) + Ky) / c;
q4_2 = atan2(val_sin2, val_cos2);

% แปลงเป็น Degree
q2_1d = rad2deg(q2_1); q4_1d = rad2deg(q4_1);
q2_2d = rad2deg(q2_2); q4_2d = rad2deg(q4_2);

fprintf('=== Results ===\n');
fprintf('Known Yellow(q3) = %.2f\n', theta3_deg);
fprintf('Set 1: Green(q2)=%.2f, Grey(q4)=%.2f\n', q2_1d, q4_1d);
fprintf('Set 2: Green(q2)=%.2f, Grey(q4)=%.2f\n', q2_2d, q4_2d);

%% 3. VECTOR DEFINITION (Exp style like Professor)
% สร้าง Vector จากผลลัพธ์ที่ได้
% เราจะ Plot ทั้ง 2 sets เพื่อเทียบดูว่าอันไหน Grey ตั้งฉาก

% --- Set 1 ---
RA_1 = a*exp(1j*q2_1);      % Green
RBA_1 = b*exp(1j*q3);       % Yellow (Known q3)
RB_1 = RA_1 + RBA_1;        % Resultant B (Loop Left side)

RO4O2 = d*exp(1j*0);        % Ground
RBO4_1 = c*exp(1j*q4_1);    % Grey
RB_Right_1 = RO4O2 + RBO4_1; % Resultant B (Loop Right side)

% --- Set 2 ---
RA_2 = a*exp(1j*q2_2);      % Green
RBA_2 = b*exp(1j*q3);       % Yellow (Known q3)
RB_2 = RA_2 + RBA_2;

RBO4_2 = c*exp(1j*q4_2);    % Grey

% แยก Component
RO4O2x = real(RO4O2); RO4O2y = imag(RO4O2);

% Set 1 Components
RA1x = real(RA_1); RA1y = imag(RA_1);
RBA1x = real(RBA_1); RBA1y = imag(RBA_1);
RB1x = real(RB_1); RB1y = imag(RB_1);
RBO41x = real(RBO4_1); RBO41y = imag(RBO4_1);

% Set 2 Components
RA2x = real(RA_2); RA2y = imag(RA_2);
RBA2x = real(RBA_2); RBA2y = imag(RBA_2);
RB2x = real(RB_2); RB2y = imag(RB_2);
RBO42x = real(RBO4_2); RBO42y = imag(RBO4_2);

%% 4. PLOTTING (Professor's Quiver Style)
figure('Color','w','Position',[100 100 1000 500]);

% --- Plot Set 1 ---
subplot(1,2,1); hold on; axis equal; grid on;
title(['Set 1: Grey Angle = ' num2str(q4_1d,'%.1f')]);
xlabel('X'); ylabel('Y');

% Loop Left Side (Green + Yellow)
quiver(0, 0, RA1x, RA1y, 0, 'g', 'LineWidth', 3, 'MaxHeadSize', 0.5); % Green
quiver(RA1x, RA1y, RBA1x, RBA1y, 0, 'y', 'LineWidth', 3, 'MaxHeadSize', 0.5); % Yellow

% Loop Right Side (Ground + Grey)
quiver(0, 0, RO4O2x, RO4O2y, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Ground
quiver(RO4O2x, RO4O2y, RBO41x, RBO41y, 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 3, 'MaxHeadSize', 0.5); % Grey

% Resultant Check (Green line to B)
plot([0 RB1x], [0 RB1y], 'g:'); 
plot(RB1x, RB1y, 'ko', 'MarkerFaceColor', 'w');


% --- Plot Set 2 ---
subplot(1,2,2); hold on; axis equal; grid on;
title(['Set 2: Grey Angle = ' num2str(q4_2d,'%.1f')]);
xlabel('X'); ylabel('Y');

% Loop Left Side (Green + Yellow)
quiver(0, 0, RA2x, RA2y, 0, 'g', 'LineWidth', 3, 'MaxHeadSize', 0.5); % Green
quiver(RA2x, RA2y, RBA2x, RBA2y, 0, 'y', 'LineWidth', 3, 'MaxHeadSize', 0.5); % Yellow

% Loop Right Side (Ground + Grey)
quiver(0, 0, RO4O2x, RO4O2y, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Ground
quiver(RO4O2x, RO4O2y, RBO42x, RBO42y, 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 3, 'MaxHeadSize', 0.5); % Grey

% Resultant Check
plot(RB2x, RB2y, 'ko', 'MarkerFaceColor', 'w');
