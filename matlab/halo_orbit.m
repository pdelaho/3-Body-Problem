%% Circular Restricted 3 Body-Problem Library for point mass
% Data provided for the Earth-Moon system

% Pauline de la HOGUE MORAN 06.24.2024

%% Data for the circular restricted 3-body problem (CR3BP) Earth-Moon
% All quantities have been adimensionned using r12 for the distances and TU
% for the time quantities

r12 = 389703; % km, distance between primary attractors
mu = 1.215058560962404e-2; % no unit, mass parameter of the system
TU = 382981; % s, inverse of the relative angular frequency between the
% two primary attractors

%% Computing intial conditions for the Halo orbit
% Initial guess using https://ssd.jpl.nasa.gov/tools/periodic_orbits.html
x0  = 1.0225293572044083E+0; % nd, for Lyapounov orbit about L1: 8.2967381582787081E-1
y0  = 8.3946315691394390E-27; % nd, for Lyapounov orbit about L1: 4.5691881617111996E-29
z0  = 1.8244450659922309E-1; % nd, for Lyapounov orbit about L1: -2.4095847440443644E-32
vx0 = 5.0465868764760129E-14; % nd, for Lyapounov orbit about L1: 2.7691850370932105E-16
vy0 = -1.0435995011600951E-1; % nd, for Lyapounov orbit about L1: 6.4159717067070993E-2
vz0 = -8.0801218229562511E-13; % nd, for Lyapounov orbit about L1: 4.2674206516771668E-32
L1x = 0.83691513; % nd, along the x direction
L2x = 1.15568217;  % nd, along the x direction

initial_position = [x0 y0 z0]';
initial_velocity = [vx0 vy0 vz0]';
initial_STM = [1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1]';
initial_conditions = [initial_position; initial_velocity; initial_STM; mu];

period = 1.5178035526683518E+0; % in TU, for Lyapounov orbit about L1: 2.7041588513971861E+0
t_sim = [0 period];
[t_orbit, y_orbit] = ode113(@halo_propagator_point_mass_with_STM, t_sim, ...
    initial_conditions);

%% Plot of the initial guess
figure(1)
plot3(y_orbit(:, 1) * r12 * 1e-3 - (1 - mu) * r12 * 1e-3, y_orbit(:, 2) * r12 * 1e-3, ...
    y_orbit(:, 3) * r12 * 1e-3)
% Add the Earth and the Moon in the visualization
% hold on
% plot3(-mu*r12*1e-3, 0, 0,'ro')
hold on
plot3(0, 0, 0, 'mo')
hold on
plot3(L2x * r12 * 1e-3 - (1 - mu) * r12 * 1e-3, 0, 0, 'go')
axis equal
grid on
xlabel('X axis *1e3[km]')
ylabel('Y axis *1e3[km]')
zlabel('Z axis *1e3[km]')
legend('Halo orbit', 'Moon', 'L2')
hold off

%% Single-shooting differenciation correction to correct initial guess
% Could make a different function for this part

[adjusted_conditions, tf] = opti(initial_conditions, period, 1000, 1e-5);

%% Plotting the Halo orbit with the adjusted initial conditions
% adjusted_conditions,tf
[t_orbit,y_orbit] = ode113(@halo_propagator_point_mass_with_STM, ...
    [0 2 * tf], adjusted_conditions);

figure(2)
plot3(y_orbit(:, 1) * r12 * 1e-3 - (1 - mu) * r12 * 1e-3, y_orbit(:, 2) * r12 * 1e-3, ...
    y_orbit(:, 3) * r12 *1e-3)
% Add the Earth and the Moon in the visualization
% hold on
% plot3(-mu*r12*1e-3, 0, 0,'ro')
hold on
plot3(0, 0, 0, 'mo')
hold on
plot3(L2x * r12 * 1e-3 - (1 - mu) * r12 * 1e-3, 0, 0, 'go')
axis equal
grid on
xlabel('X axis *1e3[km]')
ylabel('Y axis *1e3[km]')
zlabel('Z axis *1e3[km]')
legend('Halo orbit', 'Moon', 'L2')
hold off

%% Computing and plotting the manifolds trajectories

% Get the monodromy matrix = state transition matrix after one period
M = [y_orbit(end,  7) y_orbit(end,  8) y_orbit(end,  9) y_orbit(end, 10) y_orbit(end, 11) y_orbit(end, 12);
     y_orbit(end, 13) y_orbit(end, 14) y_orbit(end, 15) y_orbit(end, 16) y_orbit(end, 17) y_orbit(end, 18);
     y_orbit(end, 19) y_orbit(end, 20) y_orbit(end, 21) y_orbit(end, 22) y_orbit(end, 23) y_orbit(end, 24);
     y_orbit(end, 25) y_orbit(end, 26) y_orbit(end, 27) y_orbit(end, 28) y_orbit(end, 29) y_orbit(end, 30);
     y_orbit(end, 31) y_orbit(end, 32) y_orbit(end, 33) y_orbit(end, 34) y_orbit(end, 35) y_orbit(end, 36);
     y_orbit(end, 37) y_orbit(end, 38) y_orbit(end, 39) y_orbit(end, 40) y_orbit(end, 41) y_orbit(end, 42)];

% Getting the stable and unstable eigenvectors of the monodromy matrix
[vectors,values] = eig(M);

% The eigen values > 1 are stable and < 1 are unstable
vs = [];
vu = [];
lambdas = [];
lambdau = [];
for i=1:length(values(1, :))
    values(i, i)
    if (imag(values(i, i)) == 0) && (real(values(i, i)) > 10)
        % unstable eigenvalue
        lambdau = [lambdau values(i, i)];
        vu = [vu vectors(:, i)];
    end
    if (imag(values(i, i)) == 0) && (real(values(i, i)) < 0.1)
        % stable eigenvalue
        lambdas = [lambdas values(i, i)];
        vs = [vs vectors(:, i)];
    end
end

% Create different starting points for the manifold trajectories
n = 50; % number of manifold trajectories to be plotted
t_starts = linspace(0, period, n+1);
eps = 1e-4; % magnitude of the perutrbation to get manifold trajectories
% Using the function to scale the perturbation
eps_stable = scale_perturbation(mu, r12, adjusted_conditions(1:6), 2 * tf, ...
    vs(:, 1), vu(:, 1), M, true)
eps_unstable = scale_perturbation(mu, r12, adjusted_conditions(1:6), 2 * tf, ...
    vs(:, 1), vu(:, 1), M, false)
% eps = eps_stable
t_manifold = 4; % how long we want to simulate the manifolds

figure(3) % starting to plot oustide
for i=2:length(t_starts)
    % Get the point of the orbit at which the manifold starts
    [t_temp, y_temp] = ode113(@halo_propagator_point_mass_with_STM, ...
        [0 t_starts(i)], adjusted_conditions);
    STM = [y_temp(end,  7) y_temp(end,  8) y_temp(end,  9) y_temp(end, 10) y_temp(end, 11) y_temp(end, 12);
           y_temp(end, 13) y_temp(end, 14) y_temp(end, 15) y_temp(end, 16) y_temp(end, 17) y_temp(end, 18);
           y_temp(end, 19) y_temp(end, 20) y_temp(end, 21) y_temp(end, 22) y_temp(end, 23) y_temp(end, 24);
           y_temp(end, 25) y_temp(end, 26) y_temp(end, 27) y_temp(end, 28) y_temp(end, 29) y_temp(end, 30);
           y_temp(end, 31) y_temp(end, 32) y_temp(end, 33) y_temp(end, 34) y_temp(end, 35) y_temp(end, 36);
           y_temp(end, 37) y_temp(end, 38) y_temp(end, 39) y_temp(end, 40) y_temp(end, 41) y_temp(end, 42)];
    % Get the eigenvectors at the same time thanks to the state transition
    % matrix
    vs_t = STM * vs;
    vu_t = STM * vu;
    
    % Perturb the initial conditions to get the manifold trajectories
    manifold_conds = [];
    manifold_condu = [];
    for j=1:length(vs_t(1, :))
        cond_plus = y_temp(end, 1:6)' + eps * vs_t(:, j) / norm(vs_t(:, j));
%         cond_plus = adjusted_conditions(1:6) + eps * vs_t(:,j)/norm(vs_t(:,j));
        cond_minus = y_temp(end, 1:6)' - eps * vs_t(:, j) / norm(vs_t(:, j));
        manifold_conds = [manifold_conds cond_plus cond_minus]; % cond_minus];
    end
    for j=1:length(vu_t(1, :))
        cond_plus = y_temp(end, 1:6)' + eps * vu_t(:, j) / norm(vu_t(:, j));
%         cond_plus = adjusted_conditions(1:6) + eps * vu_t(:,j)/norm(vu_t(:,j));
        cond_minus = y_temp(end, 1:6)' - eps * vu_t(:, j) / norm(vu_t(:, j));
        manifold_condu = [manifold_condu cond_plus cond_minus]; % cond_minus];
    end

    % Now we propagate the manifold trajectories
%     manifold_condu
    for j=1:length(manifold_conds(1, :))
        [t_prop, y_prop] = ode113(@halo_propagator_point_mass, ...
            [0 - t_manifold], [manifold_conds(:, j); mu]);
        plot3(y_prop(:, 1) - (1 - mu), y_prop(:, 2), y_prop(:, 3), 'b') % 'DisplayName','Stable Manifold'
        hold on
    end
    for j=1:length(manifold_condu(1, :))
        [t_prop, y_prop] = ode113(@halo_propagator_point_mass, ...
            [0 t_manifold], [manifold_condu(:, j); mu]);
        plot3(y_prop(:, 1) - (1 - mu), y_prop(:, 2), y_prop(:, 3), 'k') % 'DisplayName','Unstable Manifold'
        hold on
    end
end
plot3(y_orbit(:, 1) - (1 - mu), y_orbit(:, 2), ...
    y_orbit(:, 3),'r','DisplayName','Original halo orbit', 'LineWidth',1)
hold on
plot3(0, 0, 0, 'mo','DisplayName','Moon', 'LineWidth',1)
hold on
plot3(L1x - (1 - mu), 0, 0, 'go','DisplayName','L2','LineWidth',1)
hold on
plot3(L2x - (1 - mu), 0, 0, 'go','DisplayName','L2','LineWidth',1)
% hold on
% plot3(-1, 0, 0, 'bo','DisplayName','Earth')
axis equal
grid on
title('Manifolds, stable in blue, unstable in black, Moon in magenta, L2 in cyan, L1 in green')
xlabel('X axis [nd]')
ylabel('Y axis [nd]')
zlabel('Z axis [nd]')
% legend
hold off

%% Position propagation functions

function statedot = halo_propagator_point_mass(t, state)
% Computing the derivative of the state vector

% Inputs

% t - time
% state = [x y z xdot ydot zdot mu] - state vector [position, velocity, 
% mass parameter]

% Output

% statedot - derivative of the state vector at time t given the state
% vector

    x   = state(1);
    y   = state(2);
    z   = state(3);
    vx  = state(4);
    vy  = state(5);
    mu  = state(7);
    rB1 = sqrt((x     + mu)^2 + y^2 + z^2); % distance to primary attractor
    rB2 = sqrt((x - 1 + mu)^2 + y^2 + z^2); % distance to secondary attractor

    statedot = zeros(7, 1); % 7
    statedot(1:3) = state(4:6);
    statedot(4) = x + 2 * vy - (1 - mu) * (x + mu) / (rB1^3) - mu * (x - 1 + mu) / (rB2^3);
    statedot(5) = y - 2 * vx - (1 - mu) * y / (rB1^3) - mu * y / (rB2^3);
    statedot(6) = - (1 - mu) * z / (rB1^3) - mu * z / (rB2^3);
end

function statedot = halo_propagator_point_mass_with_STM(t, state)
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

    x   = state(1);
    y   = state(2);
    z   = state(3);
    vx  = state(4);
    vy  = state(5);
    mu  = state(43);
    rB1 = sqrt((x     + mu)^2 + y^2 + z^2); % distance to primary attractor
    rB2 = sqrt((x - 1 + mu)^2 + y^2 + z^2); % distance to secondary attractor

    statedot = zeros(43, 1); % 3 position + 3 velocity + 6x6 state
    % transition matrix + 1 constant parameter mu
    statedot(1:3) = state(4:6);
    statedot(4) = x + 2 * vy - (1 - mu) * (x + mu) / (rB1^3) - mu * (x - 1 + mu) / (rB2^3);
    statedot(5) = y - 2 * vx - (1 - mu) * y / (rB1^3) - mu * y / (rB2^3);
    statedot(6) = - (1 - mu) * z / (rB1^3) - mu * z / (rB2^3);

    % Computing the second derivative of the pseudo-potential
    % U = (1-mu)/rB1 + mu/rB2 + (1/2)*(x^2 + y^2)
    dUdxx = 1 - (1 - mu) / (rB1^3) + 3 * (1 - mu) * (x + mu)^2 / (rB1^5) - mu / (rB2^3) ...
        + 3 * mu * (x - 1 + mu)^2 / (rB2^5);
    dUdyy = 1 - (1 - mu) / (rB1^3) + 3 * (1 - mu) * y^2 / (rB1^5) - mu / (rB2^3) ...
        + 3 * mu * y^2 / (rB2^5);
    dUdzz = - (1 - mu) / (rB1^3) + 3 * (1 - mu) * z^2 / (rB1^5) - mu / (rB2^3) ...
        + 3 * mu * z^2 / (rB2^5);
    dUdxy = 3 * (1 - mu) * (x + mu) * y / (rB1^5) + 3 * mu * (x - 1 + mu) * y / (rB2^5);
    dUdxz = 3 * (1 - mu) * (x + mu) * z / (rB1^5) + 3 * mu * (x - 1 + mu) * z / (rB2^5);
    dUdyz = 3 * (1 - mu) * y * z / (rB1^5) + 3 * mu * y * z / (rB2^5);

    % Computing the matrix A such that [xdot; xddot] = A*[x; xdot]
    % Linearization of the equations of motion
    A = [0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1;
         dUdxx dUdxy dUdxz 0 2 0;
         dUdxy dUdyy dUdyz -2 0 0;
         dUdxz dUdyz dUdzz 0 0 0];

    STM = [state(7)  state(8)  state(9)  state(10) state(11) state(12);
           state(13) state(14) state(15) state(16) state(17) state(18);
           state(19) state(20) state(21) state(22) state(23) state(24);
           state(25) state(26) state(27) state(28) state(29) state(30);
           state(31) state(32) state(33) state(34) state(35) state(36);
           state(37) state(38) state(39) state(40) state(41) state(42)];

    % Computing the derivative of the state transition matrix
    dSTMdt = A * STM;
    statedot(7:12) = dSTMdt(1,:);
    statedot(13:18) = dSTMdt(2,:);
    statedot(19:24) = dSTMdt(3,:);
    statedot(25:30) = dSTMdt(4,:);
    statedot(31:36) = dSTMdt(5,:);
    statedot(37:42) = dSTMdt(6,:);
end

function new_initial_state = single_shooting(init_state, res,jac)
    new_initial_state = init_state - pinv(jac) * res;
end

function [new_initial_conditions, tf] = opti(initial_conditions, period, imax, tol)
    adjusted_conditions = initial_conditions;
    tf = period / 2;

    for i=1:imax
        i
    %     adjusted_conditions
        [t_temp, y_temp] = ode113(@halo_propagator_point_mass_with_STM, ...
            [0 tf], adjusted_conditions);
        f = [y_temp(end, 2) y_temp(end, 4) y_temp(end, 6)]';
        norm(f)
        if norm(f) < tol
            adjusted_conditions(1) = y_temp(1, 1);
            adjusted_conditions(5) = y_temp(1, 5);
            break
        else
            % if we wanna keep the starting z altitude constant but we have
            % to change the period of the orbit
            state_end = halo_propagator_point_mass_with_STM(t_temp(end), ...
                y_temp(end, :));
            % Computing the Jacobian
    %         df = [state_end(13) state_end(17) state_end(2);
    %               state_end(25) state_end(29) state_end(4);
    %               state_end(37) state_end(41) state_end(6)];
            df = [y_temp(end, 13) y_temp(end, 17) state_end(2);
                  y_temp(end, 25) y_temp(end, 29) state_end(4);
                  y_temp(end, 37) y_temp(end, 41) state_end(6)];
            new_x = single_shooting([adjusted_conditions(1); ...
                adjusted_conditions(5); tf], f, df);
            adjusted_conditions(1) = new_x(1);
            adjusted_conditions(5) = new_x(2);
            tf = new_x(3);
        end
    end

    new_initial_conditions = adjusted_conditions;
end

function eps = scale_perturbation(mu, r12, init_state, period, vs0, vu0, monodromy, stable)
    % Setting tolerances
    relativeTolerance = 0.1;
    absoluteTolerance_km = 100;
    absoluteTolerance = absoluteTolerance_km / r12;

    % Assigning tolerances
    relTol = relativeTolerance;
    absTol = absoluteTolerance;

    % Perturbations to try
    perturbations_list_km = [0.01 0.05 0.1 0.5 1 5 10 50 100 500 1000];
    perturbations_list = perturbations_list_km / r12;

    if stable
        ef = monodromy \ vs0;
    else
        ef = monodromy \ vu0;
    end

    % Compute position error
    efPosition = norm(ef(1:3));

    eps = perturbations_list(1);

    for i=1:length(perturbations_list)
        epsilon = perturbations_list(i);
        predictedError = epsilon * efPosition;

        % Compute error by propagating 1 period
        if stable
            x0 = init_state + epsilon * vs0;
            [~, manifold] = ode113(@halo_propagator_point_mass,[0 - period], ...
                [x0; mu]); 
        else
            x0 = init_state + epsilon * vu0;
            [~, manifold] = ode113(@halo_propagator_point_mass,[0 period], ...
                [x0; mu]); 
        end

        % Actual position error
        positionError = norm(manifold(end, 1:3) - init_state(1:3));

        % Break loop if tolerances are broken
        if abs((positionError - predictedError) / positionError) > relTol && abs(positionError - predictedError) > absTol
            break
        end

        eps = epsilon;
    end

end