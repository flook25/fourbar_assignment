clear all; close all; clc;

%% --- MAIN SETUP ---
% กำหนดค่ามุม Input หลัก (ก้านสีเหลือง Loop 1)
theta_Yellow_Input_Deg = 19.94;
theta_Yellow_Input_Rad = deg2rad(theta_Yellow_Input_Deg);

% เรียกฟังก์ชันคำนวณสำหรับทั้ง 2 กรณี
% Case 1: Open Configuration (ชุดคำตอบที่ 1 ของแต่ละลูป)
analyze_and_plot_mechanism(theta_Yellow_Input_Rad, 1);

% Case 2: Crossed Configuration (ชุดคำตอบที่ 2 ของแต่ละลูป)
analyze_and_plot_mechanism(theta_Yellow_Input_Rad, 2);

%% --- FUNCTION: Analyze and Plot Mechanism ---
function analyze_and_plot_mechanism(theta_in, config_mode)
    % config_mode: 1 = Open (Solution 1), 2 = Crossed (Solution 2)
    
    fprintf('\n========================================\n');
    if config_mode == 1
        fprintf('ANALYSIS: CASE 1 (OPEN CONFIGURATION)\n');
        fig_name = 'Case 1: Open Configuration';
    else
        fprintf('ANALYSIS: CASE 2 (CROSSED CONFIGURATION)\n');
        fig_name = 'Case 2: Crossed Configuration';
    end
    fprintf('========================================\n');

    %% --- LOOP 1 ANALYSIS (Green-Yellow-Grey) ---
    % Parameters
    d1 = 210; % Ground
    a1 = 180; % Green (Link 2 @ O2)
    b1 = 180; % Yellow (Link 3 Coupler - INPUT)
    c1 = 118; % Grey (Link 4 @ O4)
    
    % Solver (Coupler Input -> Find Green & Grey)
    % Vector Loop: R_Green + R_Yellow - R_Grey - R_Ground = 0
    % แก้หา Theta_Green (t2) และ Theta_Grey (t4)
    
    % P = b*e^j*t3 - d
    Px = b1 * cos(theta_in) - d1;
    Py = b1 * sin(theta_in);
    P_mag_sq = Px^2 + Py^2;
    
    % รูปแบบสมการสำหรับหา Green (a): A*sin(t2) + B*cos(t2) + C = 0
    A_L1 = 2 * a1 * Py;
    B_L1 = 2 * a1 * Px;
    C_L1 = a1^2 + P_mag_sq - c1^2;
    
    [th_Green, status1] = solve_quadratic_angle(A_L1, B_L1, C_L1, config_mode);
    
    if status1 == 0, error('Loop 1 Assembly Failed'); end
    
    % หา Theta_Grey (Output Loop 1)
    % R_Grey = R_Green + R_Yellow - R_Ground
    R_Grey_Vec = a1*exp(1j*th_Green) + b1*exp(1j*theta_in) - d1;
    th_Grey = angle(R_Grey_Vec);
    
    fprintf('LOOP 1 RESULTS:\n');
    fprintf('  Theta Yellow (Input): %8.2f deg\n', rad2deg(theta_in));
    fprintf('  Theta Green         : %8.2f deg\n', rad2deg(th_Green));
    fprintf('  Theta Grey (Output) : %8.2f deg\n', rad2deg(th_Grey));

    %% --- LOOP 2 ANALYSIS (Brown-Blue-LightBlue) ---
    % Parameters
    d2 = 210; % Ground
    a2 = 118; % Brown (Link 2 @ O4 - INPUT from Loop 1)
    b2 = 210; % Blue (Coupler)
    c2 = 168; % Light Blue (Link 4 @ O2 - Output)
    
    % Mapping Input: สมมติว่า Brown เชื่อมกับ Grey (Loop 1)
    % ถ้า config_mode == 1 (Open), Grey ชี้ขึ้น, Brown อาจต้องชี้ลงตามโจทย์?
    % แต่เพื่อให้เป็นระบบ เราจะใช้ค่าที่ส่งต่อมาโดยตรง
    % *** ปรับแก้: ถ้าต้องการให้ Brown ชี้ลง (Below Ground) ใน Case ที่ถูกต้อง
    % เราอาจต้องกลับเครื่องหมาย หรือเลือก Solution ที่ถูกต้องจาก Loop 1
    % ในที่นี้ใช้ค่าดิบจาก Loop 1 ส่งไปก่อน
    theta_Brown = th_Grey; 
    
    % Solver (Crank Input -> Find Blue & Light Blue)
    % Standard Norton for Crank Input (a2)
    K1 = d2/a2; K2 = d2/c2; K3 = (a2^2 - b2^2 + c2^2 + d2^2)/(2*a2*c2);
    
    A_L2 = cos(theta_Brown) - K1 - K2*cos(theta_Brown) + K3;
    B_L2 = -2*sin(theta_Brown);
    C_L2 = K1 - (K2+1)*cos(theta_Brown) + K3;
    
    [th_LightBlue_L2, status2] = solve_quadratic_angle(A_L2, B_L2, C_L2, config_mode);
    if status2 == 0, error('Loop 2 Assembly Failed'); end
    
    % Find Coupler (Blue) - ใช้ Vector Loop Check
    % R_Blue = R_Ground + R_LightBlue - R_Brown (Ground Direction O4->O2 Reversed in Norton?)
    % Standard: R2(Brown) + R3(Blue) = R1(Ground) + R4(LightBlue) -> R3 = d + R4 - R2
    R_Blue_Vec = d2 + c2*exp(1j*th_LightBlue_L2) - a2*exp(1j*theta_Brown);
    th_Blue = angle(R_Blue_Vec);
    
    fprintf('LOOP 2 RESULTS:\n');
    fprintf('  Theta Brown (Input) : %8.2f deg\n', rad2deg(theta_Brown));
    fprintf('  Theta Blue          : %8.2f deg\n', rad2deg(th_Blue));
    fprintf('  Theta LightBlue(Out): %8.2f deg\n', rad2deg(th_LightBlue_L2));

    %% --- LOOP 3 ANALYSIS (LightBlue-Red-Grey) ---
    % Parameters
    d3 = 210; % Ground
    a3 = 118; % Light Blue (Link 2 @ O2 - INPUT from Loop 2)
    b3 = 210; % Red (Coupler)
    c3 = 118; % Grey (Link 4 @ O4 - Output)
    
    % Mapping Input: Light Blue จาก Loop 2
    theta_LightBlue_In = th_LightBlue_L2;
    
    % Solver (Crank Input -> Find Red & Grey)
    K1 = d3/a3; K2 = d3/c3; K3 = (a3^2 - b3^2 + c3^2 + d3^2)/(2*a3*c3);
    
    A_L3 = cos(theta_LightBlue_In) - K1 - K2*cos(theta_LightBlue_In) + K3;
    B_L3 = -2*sin(theta_LightBlue_In);
    C_L3 = K1 - (K2+1)*cos(theta_LightBlue_In) + K3;
    
    [th_Grey_L3, status3] = solve_quadratic_angle(A_L3, B_L3, C_L3, config_mode);
    if status3 == 0, error('Loop 3 Assembly Failed'); end
    
    % Find Coupler (Red)
    R_Red_Vec = d3 + c3*exp(1j*th_Grey_L3) - a3*exp(1j*theta_LightBlue_In);
    th_Red = angle(R_Red_Vec);
    
    fprintf('LOOP 3 RESULTS:\n');
    fprintf('  Theta LightBlue(In) : %8.2f deg\n', rad2deg(theta_LightBlue_In));
    fprintf('  Theta Red           : %8.2f deg\n', rad2deg(th_Red));
    fprintf('  Theta Grey (Output) : %8.2f deg\n', rad2deg(th_Grey_L3));

    %% --- PLOTTING ---
    figure('Color','w','Position',[100+(config_mode-1)*50 100 800 600]); 
    hold on; axis equal; grid on; box on;
    title(fig_name);
    xlabel('X (mm)'); ylabel('Y (mm)');
    
    % Vector Definitions
    % Loop 1
    R1_G  = d1 * exp(1j*0);
    R1_Gr = a1 * exp(1j*th_Green);
    R1_Y  = b1 * exp(1j*theta_in);
    R1_Gy = c1 * exp(1j*th_Grey);
    
    % Loop 2 (O4 is Origin relative to vector math, but physically at d1)
    R2_Br = a2 * exp(1j*theta_Brown);
    R2_Bl = b2 * exp(1j*th_Blue);
    R2_LB = c2 * exp(1j*th_LightBlue_L2);
    
    % Loop 3
    R3_LB = a3 * exp(1j*theta_LightBlue_In);
    R3_Rd = b3 * exp(1j*th_Red);
    R3_Gy = c3 * exp(1j*th_Grey_L3);
    
    % Coordinates
    O2 = 0;
    O4 = d1;
    
    % Plot Ground (Pink)
    plot([0 d1], [0 0], 'Color', [1 0.07 0.57], 'LineWidth', 3, 'LineStyle', '--');
    text(0, -15, 'O2'); text(d1, -15, 'O4');
    
    % Plot Loop 1 Links
    quiver(0,0, real(R1_Gr), imag(R1_Gr), 0, 'g', 'LineWidth', 4); % Green
    quiver(real(R1_Gr), imag(R1_Gr), real(R1_Y), imag(R1_Y), 0, 'y', 'LineWidth', 4); % Yellow
    quiver(d1, 0, real(R1_Gy), imag(R1_Gy), 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'LineStyle', ':'); % Grey L1 (Ghost)
    
    % Plot Loop 2 Links
    % Brown starts at O4
    quiver(d1, 0, real(R2_Br), imag(R2_Br), 0, 'Color', [0.6 0.3 0], 'LineWidth', 4); % Brown
    % Blue connects Brown to LightBlue
    J_Brown = d1 + R2_Br;
    quiver(real(J_Brown), imag(J_Brown), real(R2_Bl), imag(R2_Bl), 0, 'b', 'LineWidth', 4); % Blue
    % Light Blue L2 (Output)
    quiver(0, 0, real(R2_LB), imag(R2_LB), 0, 'c', 'LineWidth', 2, 'LineStyle', ':'); % LightBlue L2 (Ghost)
    
    % Plot Loop 3 Links
    % Light Blue L3 (Input) - Should overlap with L2 (same angle, diff length)
    quiver(0, 0, real(R3_LB), imag(R3_LB), 0, 'c', 'LineWidth', 4); % Light Blue
    % Red connects LB to Grey
    J_LB = R3_LB;
    quiver(real(J_LB), imag(J_LB), real(R3_Rd), imag(R3_Rd), 0, 'r', 'LineWidth', 4); % Red
    % Grey L3 (Output)
    quiver(d1, 0, real(R3_Gy), imag(R3_Gy), 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 4); % Grey
    
    % Joints
    plot(real(R1_Gr), imag(R1_Gr), 'ko', 'MarkerFaceColor','w');
    plot(real(R1_Gr+R1_Y), imag(R1_Gr+R1_Y), 'ko', 'MarkerFaceColor','w');
    plot(real(d1+R2_Br), imag(d1+R2_Br), 'ko', 'MarkerFaceColor','w');
    
    xlim([-100 350]); ylim([-250 250]);
end

%% --- HELPER: Solve Quadratic for Angle (A sin + B cos + C = 0) ---
function [theta, status] = solve_quadratic_angle(A, B, C, mode)
    % Solves for theta using tangent half-angle substitution or direct atan2
    % Returns solution based on 'mode' (1 or 2)
    
    det = B^2 + A^2 - C^2; % Note: For form Asin + Bcos + C = 0
    % Or use Quadratic form for t = tan(th/2): (C-B)t^2 + 2At + (C+B) = 0
    
    A_q = C - B;
    B_q = 2 * A;
    C_q = C + B;
    
    det_q = B_q^2 - 4*A_q*C_q;
    
    if det_q < 0
        theta = 0; status = 0; return;
    end
    
    t1 = (-B_q + sqrt(det_q)) / (2*A_q);
    t2 = (-B_q - sqrt(det_q)) / (2*A_q);
    
    th1 = 2*atan(t1);
    th2 = 2*atan(t2);
    
    if mode == 1
        theta = th1;
    else
        theta = th2;
    end
    status = 1;
end
