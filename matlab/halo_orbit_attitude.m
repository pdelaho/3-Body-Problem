%% Circular Restricted 3 Body Problem Library for a 3D spacecraft (no torques)
% Propagation of the motion of the center of mass + attitude of the
% spacecraft about it, using the Modified Rodrigues Parameters or the
% quaternion representation

% Data provided for the Earth-Moon system, especially about L2

% Pauline de la HOGUE MORAN 06.25.2024

%% Data for the CR3BP considering the Earth-Moon system
% All quantities have been adimensionned using r12 for the distances and TU
% for the time quantities

r12 = 389703; % km, distance between primary attractors
mu = 1.215058560962404e-2; % no unit, mass parameter of the system
TU = 382981; % s, inverse of the relative angular frequency between the
% two primary attractors

% Using the same inertia tensor as in the paper from Calaon and Schaub 2022
inertia = [6.67 0     0;
           0    41.87 0;
           0    0     41.87] * 1e-3;

% Assuming no disturbances for now
torques = [0 0 0];

%% Computing intial conditions for the Halo orbit
% Initial guess using https://ssd.jpl.nasa.gov/tools/periodic_orbits.html

x0  = 1.0836947694764694;
y0  = 0.0;
z0  = 0.06378527023677243;
vx0 = 0.0;
vy0 = 0.2783751354704641;
vz0 = 0.0;
% L1x = 0.83691513;
L2x = 1.15568217;

initial_position = [x0 y0 z0]'; % nd
initial_velocity = [vx0 vy0 vz0]'; % nd
initial_STM = [1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1]';
% Initial angular velocity about the center of mass in principal axes
initial_omega = [0 0 1e-6]'; % rad/s
% Initial attitude aligned with the inertial reference frame
initial_MRP = [0 0 0]';
initial_quat = [0 0 0 1]';
initial_conditions_MRP = [initial_position; initial_velocity; 
                          initial_omega; initial_MRP];
initial_conditions_quat = [initial_position; initial_velocity;
                           initial_omega; initial_quat];

period = 3.315677359772673;
t_sim = [0 period]; % plotting for just one orbit
[t_orbit, y_orbit] = ode113(@(t, y) halo_propagator_quaternion(t, y, mu, ...
    inertia, torques), t_sim, initial_conditions_quat);

%% Plot of the initial guess for the position of the center of mass
figure(1)
plot3(y_orbit(:, 1) * r12 * 1e-3 - (1 - mu) * r12 * 1e-3, y_orbit(:, 2) * r12 * 1e-3, ...
    y_orbit(:, 3) * r12 * 1e-3)
% Add the Earth and the Moon in the visualization
% hold on
% plot3(-mu*r12*1e-3, 0, 0,'ro')
hold on
plot3(0, 0, 0, 'bo')
hold on
plot3(L2x * r12 * 1e-3 - (1 - mu) * r12 * 1e-3, 0, 0, 'go')
axis equal
grid on
xlabel('X axis [nd]')
ylabel('Y axis [nd]')
zlabel('Z axis [nd]')
legend('Halo orbit', 'Moon', 'L2')
hold off

%% Optimization using the single-shooting differenciation correction
adjusted_conditions = [initial_position; initial_velocity; initial_STM; 
                       mu];
tf = period / 2;
tol = 1e-5;

for i=1:1000
    i
%     adjusted_conditions
    [t_temp, y_temp] = ode113(@halo_propagator_point_mass_with_STM, ...
        [0 tf], adjusted_conditions);
    f = [y_temp(end, 2) y_temp(end, 4) y_temp(end, 6)]';
    norm(f)
    if norm(f)<tol
        adjusted_conditions(1) = y_temp(1, 1);
        adjusted_conditions(5) = y_temp(1, 5);
        break
    else
        state_end = halo_propagator_point_mass_with_STM(t_temp(end), ...
            y_temp(end, :));
        % Computing the Jacobian
        df = [state_end(13) state_end(17) state_end(2);
              state_end(25) state_end(29) state_end(4);
              state_end(37) state_end(41) state_end(6)];
        new_x = single_shooting([adjusted_conditions(1); ...
            adjusted_conditions(5); tf], f, df);
        adjusted_conditions(1) = new_x(1);
        adjusted_conditions(5) = new_x(2);
        tf = new_x(3);
    end
end

% new period is 2*tf and changed only the x position and y velocity
initial_conditions_MRP(1)  = adjusted_conditions(1);
initial_conditions_MRP(5)  = adjusted_conditions(5);
initial_conditions_quat(1) = adjusted_conditions(1);
initial_conditions_quat(5) = adjusted_conditions(5);

%% Plotting the Halo orbit with the adjusted initial conditions
[t_orbit,y_orbit] = ode113(@(t,y) halo_propagator_quaternion(t, y, mu, ...
    inertia, torques), [0 4 * tf], initial_conditions_quat);

figure(2)
plot3(y_orbit(:, 1) - (1 - mu), y_orbit(:, 2), y_orbit(:, 3))
% Add the Earth and the Moon in the visualization
% hold on
% plot3(-mu*r12*1e-3, 0, 0,'ro')
hold on
plot3(0, 0, 0, 'bo')
hold on
plot3(L2x - (1 - mu), 0, 0, 'go')
axis equal
grid on
xlabel('X axis [nd]')
ylabel('Y axis [nd]')
zlabel('Z axis [nd]')
legend('Halo orbit', 'Moon', 'L2')
hold off

%% Position and attitude propagation functions

function statedot = halo_propagator_MRP(t, state, mu, I, T)
% Computing the derivative of the state vector

% Inputs

% t - time
% state = [x y z xdot ydot zdot omegax omegay omegaz sigma1 sigma2 sigma3] 
% - state vector [position, velocity, angular velocity about principal
% axes, Modified Rodrigues Parameters (MRP)]
% mu - mass parameter of the 3-body problem
% I - inertia tensor of the spacecraft in principal axes
% T - external torques expressed in principal axes

% Output

% statedot - derivative of the state vector at time t given the state
% vector

    x  = state(1);
    y  = state(2);
    z  = state(3);
    vx = state(4);
    vy = state(5);
%     mu = state(7);
    rB1 = sqrt((x     + mu)^2 + y^2 + z^2); % distance to primary attractor
    rB2 = sqrt((x - 1 + mu)^2 + y^2 + z^2); % distance to secondary attractor
    omega = [state(7)  state(8)  state(9)]'; % angular velocity
    sigma = [state(10) state(11) state(12)]; % MRP

    if norm(sigma) > 1 % want to work with the MRP within the unitary ball
        sigma = - sigma / norm(sigma)^2;
    end

    statedot = zeros(12, 1); % 7
    statedot(1:3) = state(4:6);
    statedot(4) = x + 2 * vy - (1 - mu) * (x + mu) / (rB1^3) - mu * (x - 1 + mu) / (rB2^3);
    statedot(5) = y - 2 * vx - (1 - mu) * y / (rB1^3) - mu * y / (rB2^3);
    statedot(6) = - (1 - mu) * z / (rB1^3) - mu * z / (rB2^3);
    
    % Derivative of the angular velocity using Euler equations
    statedot(7) = (T(1) + (I(2, 2) - I(3, 3)) * omega(2) * omega(3)) / I(1,1);
    statedot(8) = (T(2) + (I(3, 3) - I(1, 1)) * omega(1) * omega(3)) / I(2,2);
    statedot(9) = (T(3) + (I(1, 1) - I(2, 2)) * omega(1) * omega(2)) / I(3,3);

    % Derivative of the MRP using the paper by Calaon and Schaub from 2022
    scross = [0         - sigma(3)  sigma(2);
               sigma(3)  0          - sigma(1);
              - sigma(2) sigma(1)   0];
    B = (1 - norm(sigma)^2) * eye(3) + 2 * scross + 2 * sigma * sigma';
    statedot(10:12) = (1 / 4) * B * omega;
    
end

function statedot = halo_propagator_quaternion(t, state, mu, I, T)
% Computing the derivative of the state vector

% Inputs

% t - time
% state = [x y z xdot ydot zdot omegax omegay omegaz q1 q2 q3 q4] - state 
% vector [position, velocity, angular velocity about principal axes,
% quaternion (scalar last)]
% mu - mass parameter of the 3-body problem
% I - inertia tensor of the spacecraft in principal axes
% T - external torques expressed in principal axes

% Output

% statedot - derivative of the state vector at time t given the state
% vector

    x = state(1);
    y = state(2);
    z = state(3);
    vx = state(4);
    vy = state(5);
%     mu = state(7);
    rB1 = sqrt((x+mu)^2 + y^2 + z^2); % distance to primary attractor
    rB2 = sqrt((x-1+mu)^2 + y^2 + z^2); % distance to secondary attractor
    omega = [state(7) state(8) state(9)]';
    quaternion = [state(10) state(11) state(12) state(13)]';
    quaternion = quaternion / norm(quaternion); % normalize the quaternion
    % to ensure that its norm stays equal to 1

    statedot = zeros(13,1);
    statedot(1:3) = state(4:6);
    statedot(4) = x + 2*vy - (1-mu)*(x+mu)/(rB1^3) - mu*(x-1+mu)/(rB2^3);
    statedot(5) = y - 2*vx - (1-mu)*y/(rB1^3) - mu*y/(rB2^3);
    statedot(6) = -(1-mu)*z/(rB1^3) - mu*z/(rB2^3);

    % Derivative of the angular velocity using Euler equations
    statedot(7) = (T(1) + (I(2,2) - I(3,3))*omega(2)*omega(3))/I(1,1);
    statedot(8) = (T(2) + (I(3,3) - I(1,1))*omega(1)*omega(3))/I(2,2);
    statedot(9) = (T(3) + (I(1,1) - I(2,2))*omega(1)*omega(2))/I(3,3);

    % Detivative of the quaternion
    Omega = [0 omega(3) -omega(2) omega(1);
             -omega(3) 0 omega(1) omega(2);
             omega(2) -omega(1) 0 omega(3);
             -omega(1) -omega(2) -omega(3) 0];
    statedot(10:13) = (1/2) * Omega * quaternion;
end

function statedot = halo_propagator_point_mass_with_STM(t,state)
% Computing the derivative of the state vector, including the state
% transition matrix

% Inputs

% t - time
% state = [x y z xdot ydot zdot a11 a12 a13 a14 a15 a16 a21 a22 a23 a24 a25
% a26 a31 a32 a33 a34 a35 a36 a41 a42 a43 a44 a45 a46 a51 a52 a53 a54 a56
% a61 a62 a63 a64 a65 a66 mu] - state vector [position, velocity, state
% transition matrix, mass parameter]

% Output

% statedot - derivative of the state vector at time t given the state
% vector

    x = state(1);
    y = state(2);
    z = state(3);
    vx = state(4);
    vy = state(5);
    mu = state(43);
    rB1 = sqrt((x+mu)^2 + y^2 + z^2); % distance to primary attractor
    rB2 = sqrt((x-1+mu)^2 + y^2 + z^2); % distance to secondary attractor

    statedot = zeros(43,1); % 3 position + 3 velocity + 6x6 state
    % transition matrix + 1 constant parameter mu
    statedot(1:3) = state(4:6);
    statedot(4) = x + 2*vy - (1-mu)*(x+mu)/(rB1^3) - mu*(x-1+mu)/(rB2^3);
    statedot(5) = y - 2*vx - (1-mu)*y/(rB1^3) - mu*y/(rB2^3);
    statedot(6) = -(1-mu)*z/(rB1^3) - mu*z/(rB2^3);

    % Computing the second derivative of the pseudo-potential
    % U = (1-mu)/rB1 + mu/rB2 + (1/2)*(x^2 + y^2)
    dUdxx = 1 - (1-mu)/(rB1^3) + 3*(1-mu)*(x+mu)^2/(rB1^5) - mu/(rB2^3) ...
        + 3*mu*(x-1+mu)^2/(rB2^5);
    dUdyy = 1 - (1-mu)/(rB1^3) + 3*(1-mu)*y^2/(rB1^5) - mu/(rB2^3) ...
        + 3*mu*y^2/(rB2^5);
    dUdzz = -(1-mu)/(rB1^3) + 3*(1-mu)*z^2/(rB1^5) - mu/(rB2^3) ...
        + 3*mu*z^2/(rB2^5);
    dUdxy = 3*(1-mu)*(x+mu)*y/(rB1^5) + 3*mu*(x-1+mu)*y/(rB2^5);
    dUdxz = 3*(1-mu)*(x+mu)*z/(rB1^5) + 3*mu*(x-1+mu)*z/(rB2^5);
    dUdyz = 3*(1-mu)*y*z/(rB1^5) + 3*mu*y*z/(rB2^5);

    % Computing the matrix A such that [xdot; xddot] = A*[x; xdot]
    % Linearization of the equations of motion
    A = [0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1;
         dUdxx dUdxy dUdxz 0 2 0;
         dUdxy dUdyy dUdyz -2 0 0;
         dUdxz dUdyz dUdzz 0 0 0];

    STM = [state(7) state(8) state(9) state(10) state(11) state(12);
           state(13) state(14) state(15) state(16) state(17) state(18);
           state(19) state(20) state(21) state(22) state(23) state(24);
           state(25) state(26) state(27) state(28) state(29) state(30);
           state(31) state(32) state(33) state(34) state(35) state(36);
           state(37) state(38) state(39) state(40) state(41) state(42)];

    % Computing the derivative of the state transition matrix
    dSTMdt = A*STM;
    statedot(7:12) = dSTMdt(1,:);
    statedot(13:18) = dSTMdt(2,:);
    statedot(19:24) = dSTMdt(3,:);
    statedot(25:30) = dSTMdt(4,:);
    statedot(31:36) = dSTMdt(5,:);
    statedot(37:42) = dSTMdt(6,:);
end