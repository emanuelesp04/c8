clear all; close all; clc;

%% Inizializzazioni 

% Tempo di simulazione
T = [0, 1, 1.5, 16, 16.5, 17, 23.1, 23.6, 24.1, 41.6, 42.1, 42.6, 48.7, 49.2, 49.7, 52.2, 52.7, 53.7];
t_fin = T(end);

% Condizioni iniziali
X0=0;      Y0=0;        Z0 = -1500;        % Posizione iniziale
psi0 = 0;  theta0 = 0;  phi0 = 0;          % Angoli di Eulero

% Variabili p,q,r
q_max = convangvel(15.0, 'deg/s', 'rad/s');
p_max = convangvel(27.27, 'deg/s', 'rad/s');
r_max = 0;

vBreakPointsQ = [T; ...
                 0, 0, q_max, q_max, 0, 0, 0, 0, q_max, q_max, 0, 0, 0, 0, q_max, q_max, 0, 0];
vBreakPointsP = [T; ...
                 0, 0, 0, 0, 0, p_max, p_max, 0, 0, 0, 0, p_max, p_max, 0, 0, 0, 0, 0];
p = @(t) interp1(vBreakPointsP(1,:), vBreakPointsP(2,:), t, 'pchip');
q = @(t) interp1(vBreakPointsQ(1,:), vBreakPointsQ(2,:), t, 'pchip');
r = @(t) 0*t;   % stessa dimensione di t (con @(t) 0 il plot di r(t) non compare)


% Variabili u,v,w
u0 = convvel(380.0,'km/h','m/s'); 
v0 = convvel( 0.0,'km/h','m/s');  
w0 = convvel( 0.0,'km/h','m/s'); 

u = @(t) interp1( ...
   [ 0, t_fin/30, t_fin/10, t_fin/5, 0.8*t_fin, 0.9*t_fin, t_fin], ...
   [u0,   u0,   u0,      u0,    u0,   u0,    u0], ...
   t, 'pchip');
v = @(t) 0*t;
w = @(t) 0*t;


N = 120;                   
t = linspace (0, t_fin, N); 

% Plot p(t), q(t), r(t)
figure (1)
plot(t, p(t)*180/pi, '.-b', t, q(t)*180/pi, '.-r', t, r(t)*180/pi, '.-g');
hold on; grid on;
xlabel ('t (s)'); ylabel ('(deg/s)');
legend ('p(t)','q(t)', 'r(t)'); 
title ('Componenti della velocità angolare');

% Plot u(t)
figure (2)
plot(t, u(t), '.-b');
hold on; grid on;
xlabel ('t (s)'); ylabel ('(m/s)');
legend('u(t)');
title ('Componente u della velocità del baricentro');

%% Risoluzione sistema

% Condizioni iniziali quaternioni
Q0 = angle2quat(psi0, theta0, phi0);

% Equazioni cinematiche
dQuatdt = @(t, Q) ... 
0.5*[    0, -p(t), -q(t), -r(t);
      p(t),     0,  r(t), -q(t);
      q(t), -r(t),     0,  p(t);
      r(t),  q(t), -p(t),     0] * Q;

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-6*ones(1,4));

% Integrazione Runge-Kutta
[vTime, vQuat] = ode45(dQuatdt, [0 t_fin], Q0, options);

figure(3)
% Plot quaternioni 
subplot 121 
plot( ...
vTime, vQuat(:,1), '-b',  ... % q0
vTime, vQuat(:,2), '-r', ... % qx
vTime, vQuat(:,3), '-g', ... % qy
vTime, vQuat(:,4), '-k'   ... % qz
);
grid on;
xlabel('t (s)'); 
legend('q_0(t)','q_x(t)','q_y(t)','q_z(t)');
title('Componenti dei quaternioni')
axis ([0 inf -1 1]);

% Plot angoli di Eulero
subplot 122 
[vpsi, vtheta, vphi] = quat2angle(vQuat); 
plot( ...
vTime, convang(vpsi,'rad','deg'), '-b', ...
vTime, convang(vtheta,'rad','deg'), '-r', ...
vTime, convang(vphi,'rad','deg'), '-g' );
grid on;
xlabel('t (s)'); ylabel('(deg)');
legend('\psi(t)','\theta(t)','\phi(t)');
title('Angoli di Eulero');
axis ([0 inf -190 190]);

figure(4)
% Componenti velocità angolari in assi body
subplot 211
plot( ...
vTime, convangvel(p(vTime),'rad/s','deg/s'), '-b', ...
vTime, convangvel(q(vTime),'rad/s','deg/s'), '-r', ...
vTime, convangvel(r(vTime),'rad/s','deg/s'), '-g');
grid on;
xlabel('\it t (s)'); ylabel('(deg/s)');
legend('p(t)','q(t)','r(t)');
title('Componenti velocità angolare in assi body');

% Componenti delle velocità lineari in assi body 
subplot 212  
plot(vTime, u(vTime), '-b');
grid on;
xlabel('\it t (s)'); 
ylabel('(m/s)');
legend('u(t)');
title('Componente velocità lungo l''asse x in assi body');

% Funzione di interpolazione temporale della storia dei quaternioni
Quat = @(t) ...                    
[interp1(vTime,vQuat(:,1),t), ... 
 interp1(vTime,vQuat(:,2),t), ...
 interp1(vTime,vQuat(:,3),t), ... 
 interp1(vTime,vQuat(:,4),t)]; 

% Posizione iniziale
PosE0 = [X0;Y0;Z0];

% Navigation equations
dPosEdt = @(t,PosE) transpose(quat2dcm(Quat(t)))*[u(t);v(t);w(t)]; 

options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3*ones(3,1));

% Integrazione Runge-Kutta
[vTime2, vPosE] = ode45(dPosEdt, vTime, PosE0, options);

% La condizione iniziale PosE0 = [X0;Y0;Z0] e' gia' inclusa nella soluzione:
% NON va sommata di nuovo (altrimenti Z parte da 2*Z0 = -3000)
vXe = vPosE(:,1);
vYe = vPosE(:,2);
vZe = vPosE(:,3);

%% Grafici traiettoria

% Plot delle coordinate del baricentro espresse nel riferimento earth
figure(5)
plot(vTime, vXe, '-b',...
     vTime, vYe, '-r', ...
     vTime, vZe, '-g');
grid on;
xlabel('t (s)'); ylabel('(m)');
legend({'x_{G,E}(t)','y_{G,E}(t)','z_{G,E}(t)'}, 'Location', 'northwest');
title('Coordinate del baricentro negli assi Earth');
xlim([0 t_fin]);

% Plot body and trajectory 
h_fig6 = figure(6);
grid on;
theView = [1,1,0.5];
scaleFactor = 0.008;

% ode45 restituisce istanti NON equispaziati (passi piccoli dove p e q
% cambiano, grandi altrove): campionando ogni <step> punti i velivoli si
% addensano. Si ricampiona su una griglia temporale uniforme: con u costante
% i velivoli risultano equidistanti lungo la traiettoria.
nBody = 18;                          % numero di velivoli da disegnare
step  = 20;                          % punti di traiettoria tra due velivoli
tPlot = linspace(0, t_fin, nBody*step + 1)';

vXp = interp1(vTime, vXe, tPlot, 'spline');
vYp = interp1(vTime, vYe, tPlot, 'spline');
vZp = interp1(vTime, vZe, tPlot, 'spline');
vQuatP = interp1(vTime, vQuat, tPlot);           % interpolazione lineare...
vQuatP = vQuatP ./ vecnorm(vQuatP, 2, 2);        % ...e rinormalizzazione

plotTrajectoryAndBody(h_fig6,vXp,vYp,vZp,vQuatP,scaleFactor,step,theView);

% Titolo ed etichette DOPO la funzione (che fa clf e li cancellerebbe)
figure(h_fig6);
title('Traiettoria Cuban Eight');
xlabel('x_E (m)'); ylabel('y_E (m)'); zlabel('h = -z_E (m)');
hold on;

% Assi NED con origine nel punto iniziale della traiettoria.
% Il grafico usa le coordinate (x, y, h = -z_E), quindi:
%   N = +x_E -> [L 0 0],  E = +y_E -> [0 L 0],  D = +z_E -> [0 0 -L] (quota che scende)
L = 300;                                   % lunghezza delle frecce (m)
O = [vXe(1), vYe(1), -vZe(1)];             % origine (x, y, h)
quiver3(O(1),O(2),O(3), L,0,0, 0, 'r', 'LineWidth',2, 'MaxHeadSize',0.5); % N
quiver3(O(1),O(2),O(3), 0,L,0, 0, 'g', 'LineWidth',2, 'MaxHeadSize',0.5); % E
quiver3(O(1),O(2),O(3), 0,0,-L, 0, 'b', 'LineWidth',2, 'MaxHeadSize',0.5); % D
text(O(1)+1.1*L, O(2),         O(3),         'N', 'Color','r', 'FontWeight','bold');
text(O(1),         O(2)+1.1*L, O(3),         'E', 'Color','g', 'FontWeight','bold');
text(O(1),         O(2),         O(3)-1.1*L, 'D', 'Color','b', 'FontWeight','bold');

xlim([min(vXe)-300, max(vXe)+300]);   % prima era min(vXe)+300: tagliava l'inizio
ylim([-800 800]);
zlim([min([-vZe; O(3)-1.2*L])-100, max(-vZe)+300]);   % quota h = -z_E, freccia D inclusa
daspect([1 1 1]);
view(theView);
