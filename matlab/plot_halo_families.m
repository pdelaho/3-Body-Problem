%% Plotting the halo families for the CR3BP
% Data provided for the Earth-Moon system
% L2 Halo (Northern) family

% Pauline de la HOGUE MORAN 07.01.2024

%% Data for the circular restricted 3-body problem (CR3BP) Earth-Moon
% All quantities have been adimensionned using r12 for the distances and TU
% for the time quantities

r12 = 389703; % km, distance between primary attractors
mu = 1.215058560962404e-2; % no unit, mass parameter of the system
TU = 382981; % s, inverse of the relative angular frequency between the
% two primary attractors

%% Computing intial conditions for the Halo orbit
% Initial guess using https://ssd.jpl.nasa.gov/tools/periodic_orbits.html
x0 = [1.1808985497899205E+0; 1.1807299126562403E+0; 1.1802432007231030E+0; 
    1.1794814105552625E+0; 1.1784496916179399E+0; 1.1772049627468157E+0; 
    1.1757357974855753E+0; 1.1740583104016364E+0; 1.1721212352395958E+0;
    1.1699952083306175E+0; 1.1677006591947434E+0; 1.1652591082218164E+0;
    1.1626929203920950E+0; 1.1599142517039056E+0; 1.1571604941196401E+0;
    1.1542263493262155E+0; 1.1513800635090152E+0; 1.1485259549360594E+0;
    1.1455413442894928E+0; 1.1427323976670050E+0; 1.1398248385779477E+0;
    1.1369842439274436E+0; 1.1340654944747190E+0; 1.1310749646265861E+0;
    1.1281911619625644E+0; 1.1252564593293650E+0; 1.1222780089124029E+0;
    1.1192632353862111E+0; 1.1162196069793926E+0; 1.1131543850560104E+0;
    1.1098928777210029E+0; 1.1068039033622386E+0; 1.1035300322591290E+0;
    1.1002579217768833E+0; 1.0969915644292789E+0; 1.0935531551495765E+0;
    1.0901265457836595E+0; 1.0867128634138150E+0; 1.0831337559276397E+0]; % nd

y0 = [-2.5444988241150091E-26; -3.9224900862468671E-27; -4.0131763697449768E-27;
    -3.3912217928592086E-27; -5.4426459074481046E-27; 3.4536760079285946E-27;
    -2.4315077582384360E-28; 1.4616208552339063E-28; 4.3976196849120970E-29;
    -2.0139114798129488E-28; 4.0672599863558720E-28; 8.3768387931258446E-28;
    -1.6721129610932788E-28; -2.3173380562008736E-28; -6.0743358593137044E-28;
    -6.2944595408249574E-28; 7.1390708216286326E-29; 3.1565990810495496E-28;
    1.6192114933067259E-27; -9.4723070232695753E-28; 2.3639470684970793E-27;
    2.4486207109097478E-27; -2.1292148592828467E-27; -9.8852194218501502E-28;
    7.6541498074543159E-27; 3.4582554314059270E-27; 4.9312297519955555E-28;
    5.9854005070830337E-28; 1.0397089411165433E-27; 2.5445211165619221E-27;
    -3.2370214994299834E-27; -3.7763309589608048E-27; 1.4131620744948173E-27;
    1.6426007027565626E-27; -1.3906893609387349E-27; -3.0099620762482007E-27;
    -3.3929555608982004E-28; 8.0304267092243950E-30; 1.1706458956133718E-26]; % nd

z0 = [1.0295054075242347E-4; 1.3121069431778331E-2; 2.5617923919654416E-2;
    3.7178973271809496E-2; 4.8140435768697247E-2; 5.8225138268011067E-2;
    6.7801780194122280E-2; 7.6909584828747435E-2; 8.5877005926413993E-2;
    9.4422999379980818E-2; 1.0256492290418533E-1; 1.1031384321787237E-1;
    1.1767600543672109E-1; 1.2493104137908298E-1; 1.3152090959188770E-1;
    1.3799262593245051E-1; 1.4380485557811437E-1; 1.4922846889506389E-1;
    1.5451145148373718E-1; 1.5915200933380999E-1; 1.6364154884722404E-1;
    1.6773758461592672E-1; 1.7166247666605308E-1; 1.7539783791199678E-1;
    1.7873519017004103E-1; 1.8187257007554328E-1; 1.8479677771071973E-1;
    1.8749672012259902E-1; 1.8996401464715618E-1; 1.9219343403238076E-1;
    1.9429272244184692E-1; 1.9603030754249040E-1; 1.9761601341425358E-1;
    1.9894909370174629E-1; 2.0004092398988046E-1; 2.0094617631480924E-1;
    2.0161377613193826E-1; 2.0206059408320026E-1; 2.0231055618100155E-1]; % nd

vx0 = [3.3765359485568778E-15; 2.5207120235612143E-15; 2.6505334124465135E-15;
    3.5657215709679162E-15; -2.5260992021190302E-15; 9.2876289344915468E-16;
    -2.1060536016712307E-15; -2.7746128158444920E-15; -2.8399759059065597E-15;
    -4.9656594086928137E-16; -1.0417705214131128E-15; -2.1459660712392417E-15;
    6.6655743243173158E-16; 1.2307479474947080E-15; 1.4776297513298329E-15;
    4.1794997459487238E-15; 3.0283461616922221E-15; 4.2917354658873129E-15;
    3.4242613219894274E-15; 5.9154586036713890E-15; 5.3692753492503482E-15;
    4.1935762479173162E-15; 3.9613389302666476E-15; 3.4952440956136861E-15;
    -1.2630506284953841E-15; -5.5891287188789060E-16; -3.1687920458283919E-16;
    5.2868361571034389E-16; -2.3824535810123310E-15; -9.2916156853283232E-16;
    4.8664558894610642E-16; 1.1792722431053536E-15; 8.1496282395178728E-16;
    7.4879463226110969E-16; 4.9047891808896165E-15; 3.4271448531393373E-15;
    1.1386658764907653E-14; 1.1192016952016542E-14; 8.3292712433950044E-15]; % nd

vy0 = [-1.5585631393981156E-1; -1.5684752073236072E-1; -1.5954905069606692E-1;
    -1.6338139359639692E-1; -1.6797405042301869E-1; -1.7282597331742994E-1;
    -1.7782655254524235E-1; -1.8280860876068225E-1; -1.8782217382326816E-1;
    -1.9261435418133890E-1; -1.9712508121327701E-1; -2.0131471324908751E-1;
    -2.0515853269673287E-1; -2.0877759243379301E-1; -2.1188791433341295E-1;
    -2.1474891071483865E-1; -2.1712922972872784E-1; -2.1916532303628991E-1;
    -2.2095226684799862E-1; -2.2233902828510901E-1; -2.2349331292348754E-1;
    -2.2436040946170996E-1; -2.2499619681915137E-1; -2.2539063500356940E-1;
    -2.2553347911252272E-1; -2.2544675199362010E-1; -2.2512563128868165E-1;
    -2.2456703699414379E-1; -2.2376980465129670E-1; -2.2273474373483881E-1;
    -2.2138261166524692E-1; -2.1986831677048768E-1; -2.1802034988367752E-1;
    -2.1592883120638881E-1; -2.1360265558047956E-1; -2.1090275486068658E-1;
    -2.0796204285184380E-1; -2.0479030592633110E-1; -2.0121168055184560E-1]; % nd

vz0 = [5.5263881873244218E-18; 1.1007889483810504E-15; 5.4577430591805649E-15;
    4.4439926659633611E-15; 9.8199001002925144E-15; 1.1032674450864492E-14;
    1.7594110167841084E-14; 1.4809505158598828E-14; 1.6774526301835190E-16;
    3.8537854741228954E-15; 2.2798663794896066E-15; -4.7168468693615866E-16;
    -2.5760642407696903E-15; -4.7596301953914308E-15; -1.2326017791389893E-14;
    -8.0790091421123428E-15; -1.9547204857788442E-15; -2.8427505695990620E-15;
    -2.2634800484517227E-15; -3.7801420730129741E-15; -4.8011160052121524E-15;
    -4.5209590114403008E-15; -3.2397902355862227E-15; -4.8486629407260712E-15;
    -3.4633466525785658E-15; -1.5886999426493357E-15; -3.4416453607522604E-15;
    -6.7328707463436131E-15; -9.8830793422649480E-15; -9.2655794365851240E-15;
    -4.9325030141434507E-15; 1.3528599320325566E-15; -8.8844976835138455E-16;
    -4.9874747883170076E-15; 6.6909069485186488E-15; 1.0078939073787943E-15;
    -2.8195151432675860E-14; -3.6204794938766712E-14; -3.8933834984161984E-14]; % nd

L1x = 0.83691513; % nd, along the x direction
L2x = 1.15568217;  % nd, along the x direction

% initial_position = [x0 y0 z0]';
% initial_velocity = [vx0 vy0 vz0]';
initial_STM = [1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1]';
% initial_conditions = [initial_position; initial_velocity; initial_STM; mu];

period = [3.4155308065628454E+0; 3.4141232367543708E+0; 3.4101539733562030E+0;
    3.4041666282888854E+0; 3.3963805494683674E+0; 3.3873299978687315E+0; 
    3.3769669271238048E+0; 3.3653986536399993E+0; 3.3522346943474246E+0;
    3.3378824002096996E+0; 3.3223749314308342E+0; 3.3057376675493617E+0;
    3.2879941141468882E+0; 3.2683796978298347E+0; 3.2484349509162689E+0;
    3.2265351027754612E+0; 3.2045676617775203E+0; 3.1817478106944685E+0;
    3.1569602099666483E+0; 3.1327035113287223E+0; 3.1065871267784422E+0;
    3.0800303869394452E+0; 3.0516246962540947E+0; 3.0213062208215136E+0;
    2.9908794086451036E+0; 2.9587012252195288E+0; 2.9247896381352474E+0;
    2.8891920020877224E+0; 2.8519851641669658E+0; 2.8132733875316878E+0;
    2.7707861925069994E+0; 2.7293971867214561E+0; 2.6844146154370518E+0;
    2.6384276046753588E+0; 2.5916248068584653E+0; 2.5415355492607352E+0;
    2.4909336664289770E+0; 2.4399955442394532E+0; 2.3861839336103214E+0]; % in TU
% t_sim = [0 period];
% [t_orbit,y_orbit] = ode113(@halo_propagator_point_mass_with_STM, t_sim, ...
%     initial_conditions);

%% Single-shooting differenciation correction to correct initial guess
% Could make a different function for this part

figure(1)
for i=1:length(x0)
    T = period(i);
    initial_conditions = [x0(i); y0(i); z0(i); vx0(i); vy0(i); vz0(i); 
        initial_STM; mu];
    [adjusted_conditions,tf] = opti(initial_conditions,T,1000,1e-5);
    x0(i) = adjusted_conditions(1);
    y0(i) = adjusted_conditions(2);
    z0(i) = adjusted_conditions(3);
    vx0(i) = adjusted_conditions(4);
    vy0(i) = adjusted_conditions(5);
    vz0(i) = adjusted_conditions(6);
    period(i) = 2*tf;
    [t_orbit,y_orbit] = ode113(@halo_propagator_point_mass_with_STM, [0 period(i)], ...
        adjusted_conditions);
    plot3(y_orbit(:,1)*r12*1e-3-(1-mu)*r12*1e-3, y_orbit(:,2)*r12*1e-3, ...
        y_orbit(:,3)*r12*1e-3, 'b', 'DisplayName', 'Halo orbit')
    hold on
end
hold on
plot3(0, 0, 0, 'm*', 'DisplayName','Moon')
hold on
plot3(L2x*r12*1e-3-(1-mu)*r12*1e-3, 0, 0, 'ro', 'DisplayName', 'L2')
axis equal
grid on
xlabel('X axis *1e3[km]')
ylabel('Y axis *1e3[km]')
zlabel('Z axis *1e3[km]')
legend
hold off

%% Plotting the Halo orbit with the adjusted initial conditions

% figure(2)
% for i=1:length(x0)
%     initial_conditions = [x0(i)]
%     [t_orbit,y_orbit] = ode113(@halo_propagator_point_mass_with_STM, ...
%         [0 2*tf], adjusted_conditions);
%     
%     figure(2)
%     plot3(y_orbit(:,1)*r12*1e-3-(1-mu)*r12*1e-3, y_orbit(:,2)*r12*1e-3, ...
%         y_orbit(:,3)*r12*1e-3)
%     % Add the Earth and the Moon in the visualization
%     % hold on
%     % plot3(-mu*r12*1e-3, 0, 0,'ro')
%     
% end
% hold on
% plot3(0, 0, 0, 'm*')
% hold on
% plot3(L2x*r12*1e-3-(1-mu)*r12*1e-3, 0, 0, 'go')
% axis equal
% grid on
% xlabel('X axis *1e3[km]')
% ylabel('Y axis *1e3[km]')
% zlabel('Z axis *1e3[km]')
% legend('Halo orbit', 'Moon', 'L2')
% hold off
%% Position propagation functions

function statedot = halo_propagator_point_mass(t,state)
% Computing the derivative of the state vector

% Inputs

% t - time
% state = [x y z xdot ydot zdot mu] - state vector [position, velocity, 
% mass parameter]

% Output

% statedot - derivative of the state vector at time t given the state
% vector

    x = state(1);
    y = state(2);
    z = state(3);
    vx = state(4);
    vy = state(5);
    mu = state(7);
    rB1 = sqrt((x+mu)^2 + y^2 + z^2); % distance to primary attractor
    rB2 = sqrt((x-1+mu)^2 + y^2 + z^2); % distance to secondary attractor

    statedot = zeros(7,1); % 7
    statedot(1:3) = state(4:6);
    statedot(4) = x + 2*vy - (1-mu)*(x+mu)/(rB1^3) - mu*(x-1+mu)/(rB2^3);
    statedot(5) = y - 2*vx - (1-mu)*y/(rB1^3) - mu*y/(rB2^3);
    statedot(6) = -(1-mu)*z/(rB1^3) - mu*z/(rB2^3);
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

function new_initial_state = single_shooting(init_state,res,jac)
    new_initial_state = init_state - pinv(jac)*res;
end

function [new_initial_conditions,tf] = opti(initial_conditions,period,imax,tol)
    adjusted_conditions = initial_conditions;
    tf = period/2;

    for i=1:imax
        i
    %     adjusted_conditions
        [t_temp,y_temp] = ode113(@halo_propagator_point_mass_with_STM, ...
            [0 tf], adjusted_conditions);
        f = [y_temp(end,2) y_temp(end,4) y_temp(end,6)]';
        norm(f)
        if norm(f)<tol
            adjusted_conditions(1) = y_temp(1,1);
            adjusted_conditions(5) = y_temp(1,5);
            break
        else
            % if we wanna keep the starting z altitude constant but we have
            % to change the period of the orbit
            state_end = halo_propagator_point_mass_with_STM(t_temp(end), ...
                y_temp(end,:));
            % Computing the Jacobian
    %         df = [state_end(13) state_end(17) state_end(2);
    %               state_end(25) state_end(29) state_end(4);
    %               state_end(37) state_end(41) state_end(6)];
            df = [y_temp(end,13) y_temp(end,17) state_end(2);
                  y_temp(end,25) y_temp(end,29) state_end(4);
                  y_temp(end,37) y_temp(end,41) state_end(6)];
            new_x = single_shooting([adjusted_conditions(1); ...
                adjusted_conditions(5); tf], f, df);
            adjusted_conditions(1) = new_x(1);
            adjusted_conditions(5) = new_x(2);
            tf = new_x(3);
        end
    end

    new_initial_conditions = adjusted_conditions;
end
