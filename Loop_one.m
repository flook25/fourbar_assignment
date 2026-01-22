clear all; close all; clc;

%% 1. SYSTEM PARAMETERS (กำหนดตัวแปรตามโจทย์)
L1 = 210; % d (Ground)
L2 = 180; % a (Green Link)
L3 = 180; % b (Yellow Coupler Link)
L4 = 118; % c (Grey Link)

a = L2;
b = L3;
c = L4;
d = L1;

% --- INPUT SPECIFICATION ---
% โจทย์กำหนดมุมของ Link 3 (Yellow) เป็น Input
theta3_deg = 19.94;
theta3 = deg2rad(theta3_deg); % q3 in radians

%% 2. ANALYTIC SOLVER (ประยุกต์ Concept สูตรอาจารย์หา theta2)
% จากสมการ Vector Loop: a*e^(j*q2) + b*e^(j*q3) - c*e^(j*q4) - d = 0
% เราทราบ q3 ต้องการหา q2 (Green) และ q4 (Grey)
% จัดรูปสมการให้อยู่ในรูปแบบ A*sin(q2) + B*cos(q2) + C = 0 เพื่อเข้าสูตร

% คำนวณสัมประสิทธิ์ A, B, C สำหรับกรณี Coupler Input (q3 known)
% (หมายเหตุ: สูตรนี้ डिrive มาจาก Vector Loop เดียวกับของอาจารย์)
K_A = 2*a*b*sin(theta3);
K_B = 2*a*(b*cos(theta3) - d);
K_C = a^2 + b^2 + d^2 - c^2 - 2*b*d*cos(theta3);

% แปลงสัมประสิทธิ์เพื่อเข้าสูตร Quadratic (t = tan(q2/2))
% เปรียบเทียบกับรูปแบบของอาจารย์: A_quad*t^2 + B_quad*t + C_quad = 0
A_quad = K_C - K_B;
B_quad = 2 * K_A;
C_quad = K_C + K_B;

% ตรวจสอบ Discriminant
det = B_quad^2 - 4*A_quad*C_quad;
if det < 0
    error('Assembly failed: No real solution for given input angle.');
end

% คำนวณมุม q2 (Green Link) ทั้ง 2 คำตอบ (Open/Crossed)
t_1 = (-B_quad + sqrt(det)) / (2*A_quad);
t_2 = (-B_quad - sqrt(det)) / (2*A_quad);

q2_sol1 = 2*atan(t_1);
q2_sol2 = 2*atan(t_2);

% --- เลือกคำตอบ (Selection) ---
% เงื่อนไข: Green Link ต้องอยู่ใต้ Ground (มุมเป็นลบ หรือ sin(q2) < 0)
if sin(q2_sol1) < 0
    q2 = q2_sol1;
else
    q2 = q2_sol2;
end
q2_deg = rad2deg(q2);

% --- คำนวณ q4 (Grey Link) ---
% เมื่อได้ q2 แล้ว สามารถหา q4 ได้จากสมการ Vector Loop
% R4 = R2 + R3 - R1
R2_vec = a * exp(1j * q2);
R3_vec = b * exp(1j * theta3);
R1_vec = d * exp(1j * 0);
R4_vec = R2_vec + R3_vec - R1_vec;

q4 = angle(R4_vec);
q4_deg = rad2deg(q4);

%% 3. POSITION VECTORS FOR PLOTTING (ใช้รูปแบบ exp แบบอาจารย์)
% สร้าง Vector ของแต่ละก้าน
R_Green = a * exp(1j * q2);       % Link 2 (Green)
R_Yellow = b * exp(1j * theta3);  % Link 3 (Yellow)
R_Grey = c * exp(1j * q4);        % Link 4 (Grey)
R_Ground = d * exp(1j * 0);       % Link 1 (Ground)

% แยก Component Real/Imaginary สำหรับ quiver
R_Green_x = real(R_Green); R_Green_y = imag(R_Green);
R_Yellow_x = real(R_Yellow); R_Yellow_y = imag(R_Yellow);
R_Grey_x = real(R_Grey); R_Grey_y = imag(R_Grey);
R_Ground_x = real(R_Ground); R_Ground_y = imag(R_Ground);

% คำนวณพิกัดจุดปลาย (Absolute Position)
% O2 = (0,0)
Pos_J1 = R_Green;            % ปลาย Green
Pos_J2 = R_Green + R_Yellow; % ปลาย Yellow (ควรเท่ากับ R_Ground + R_Grey)

%% 4. DISPLAY RESULTS
fprintf('========================================\n');
fprintf('RESULTS (ANALYTIC METHOD - LECTURE STYLE)\n');
fprintf('========================================\n');
fprintf('Input Theta3 (Yellow): %8.2f deg\n', theta3_deg);
fprintf('----------------------------------------\n');
fprintf('Calculated Theta2 (Green): %8.2f deg\n', q2_deg);
fprintf('Calculated Theta4 (Grey):  %8.2f deg\n', q4_deg);
fprintf('========================================\n');

%% 5. PLOTTING (ใช้ quiver ตามสไตล์อาจารย์)
figure('Color','w','Position',[100 100 700 500]); hold on; axis equal; grid on; box on;
title(['Position Analysis (Lecture Concept): Input \theta_3 = ' num2str(theta3_deg) '^\circ']);
xlabel('X (mm)'); ylabel('Y (mm)');

% Plot using quiver (Vector style)
% Ground (O2 -> O4) - สีชมพูเข้ม (Dark Pink/Magenta)
quiver(0, 0, R_Ground_x, R_Ground_y, 0, 'Color', [0.8, 0, 0.5], 'MaxHeadSize', 0.5, 'LineWidth', 3, 'DisplayName', 'Ground');

% Link 2 (Green) - O2 -> J1
quiver(0, 0, R_Green_x, R_Green_y, 0, 'g', 'MaxHeadSize', 0.5, 'LineWidth', 3, 'DisplayName', 'Green');

% Link 3 (Yellow) - J1 -> J2
quiver(R_Green_x, R_Green_y, R_Yellow_x, R_Yellow_y, 0, 'y', 'MaxHeadSize', 0.5, 'LineWidth', 3, 'DisplayName', 'Yellow');

% Link 4 (Grey) - O4 -> J2 (วาดจาก O4 ขึ้นไปชน J2 เพื่อให้หัวลูกศรชนกันแบบ Vector Loop ปิด)
% หรือวาดแบบอาจารย์คือจาก O4 ไปปลาย
quiver(R_Ground_x, R_Ground_y, R_Grey_x, R_Grey_y, 0, 'Color', [0.5 0.5 0.5], 'MaxHeadSize', 0.5, 'LineWidth', 3, 'DisplayName', 'Grey');

% Add Labels
text(0, -10, 'O2', 'FontSize', 12, 'FontWeight', 'bold');
text(R_Ground_x, -10, 'O4', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast');

% ปรับแกนให้สวยงาม
xlim([-50 300]); ylim([-200 150]);
