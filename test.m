clear all; close all; clc;

%% 1. SYSTEM PARAMETERS (Loop 1: Yellow -> Brown -> Grey)
L1 = 210; % d (Ground O2-O4)
L2 = 118; % a (Yellow Link - Input)
L3 = 180; % b (Brown Link - Coupler)
L4 = 236; % c (Grey Link - Output) Note: Length to pivot O4

a = L2;
b = L3;
c = L4;
d = L1;

% Input Angle (Yellow Link)
theta2_deg = 19.94; 
q2 = deg2rad(theta2_deg); 

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

% --- Solve for q4 (Output - Grey) ---
q41 = 2*atan((-B + sqrt(B^2 - 4*A*C)) / (2*A));
q42 = 2*atan((-B - sqrt(B^2 - 4*A*C)) / (2*A));

% --- Solve for q3 (Coupler - Brown) ---
q31 = 2*atan((-E + sqrt(E^2 - 4*D*F)) / (2*D));
q32 = 2*atan((-E - sqrt(E^2 - 4*D*F)) / (2*D));

% แปลงเป็น Degree
q41d = mod(rad2deg(q41), 360);
q42d = mod(rad2deg(q42), 360);
q31d = mod(rad2deg(q31), 360);
q32d = mod(rad2deg(q32), 360);

fprintf('=== Loop 1 Results ===\n');
fprintf('Input Yellow = %.2f deg\n', theta2_deg);
fprintf('Set 1: Brown(q3)=%.2f, Grey(q4)=%.2f\n', q31d, q41d);
fprintf('Set 2: Brown(q3)=%.2f, Grey(q4)=%.2f\n', q32d, q42d);

%% 3. VECTOR DEFINITION (Exp style)
RA = a*exp(1j*q2); % Yellow Vector

% --- Set 1 (Closed/Open?) ---
RBA1 = b*exp(1j*q31);   % Brown
RBO41 = c*exp(1j*q41);  % Grey
RB1 = RA + RBA1; 

% --- Set 2 (Closed/Open?) ---
RBA2 = b*exp(1j*q32);   % Brown
RBO42 = c*exp(1j*q42);  % Grey
RB2 = RA + RBA2;

RO4O2 = d*exp(1j*0); % Ground

% แยก Component
RAx = real(RA); RAy = imag(RA);
RO4O2x = real(RO4O2); RO4O2y = imag(RO4O2);

% Set 1
RBA1x = real(RBA1); RBA1y = imag(RBA1);
RBO41x = real(RBO41); RBO41y = imag(RBO41);
RB1x = real(RB1); RB1y = imag(RB1);

% Set 2
RBA2x = real(RBA2); RBA2y = imag(RBA2);
RBO42x = real(RBO42); RBO42y = imag(RBO42);
RB2x = real(RB2); RB2y = imag(RB2);

%% 4. PLOTTING (Side-by-Side)
figure('Color','w','Position',[100 100 1000 500]);

% --- Plot Solution 1 ---
subplot(1,2,1); hold on; axis equal; grid on;
title(['Solution 1: Grey Angle = ' num2str(q41d,'%.1f') ' deg']);
xlabel('X'); ylabel('Y');
% Draw Vectors
quiver(0, 0, RAx, RAy, 0, 'y', 'LineWidth', 3, 'MaxHeadSize', 0.5); % Yellow
quiver(RAx, RAy, RBA1x, RBA1y, 0, 'Color', [0.6 0.3 0], 'LineWidth', 3, 'MaxHeadSize', 0.5); % Brown
quiver(0, 0, RO4O2x, RO4O2y, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Ground
quiver(RO4O2x, RO4O2y, RBO41x, RBO41y, 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 3, 'MaxHeadSize', 0.5); % Grey
% Joints
plot(0,0,'ko'); plot(RO4O2x,RO4O2y,'ko');
plot(RAx,RAy,'ko','MarkerFaceColor','w');
plot(RB1x,RB1y,'ko','MarkerFaceColor','w');


% --- Plot Solution 2 ---
subplot(1,2,2); hold on; axis equal; grid on;
title(['Solution 2: Grey Angle = ' num2str(q42d,'%.1f') ' deg']);
xlabel('X'); ylabel('Y');
% Draw Vectors
quiver(0, 0, RAx, RAy, 0, 'y', 'LineWidth', 3, 'MaxHeadSize', 0.5); % Yellow
quiver(RAx, RAy, RBA2x, RBA2y, 0, 'Color', [0.6 0.3 0], 'LineWidth', 3, 'MaxHeadSize', 0.5); % Brown
quiver(0, 0, RO4O2x, RO4O2y, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Ground
quiver(RO4O2x, RO4O2y, RBO42x, RBO42y, 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 3, 'MaxHeadSize', 0.5); % Grey
% Joints
plot(0,0,'ko'); plot(RO4O2x,RO4O2y,'ko');
plot(RAx,RAy,'ko','MarkerFaceColor','w');
plot(RB2x,RB2y,'ko','MarkerFaceColor','w');
