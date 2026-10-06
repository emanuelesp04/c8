function plotTrajectoryAndBody(h_fig,vXe,vYe,vZe,vQuat,scale_factor,step,theView,varargin)
%
%   function plotTrajectoryAndBody(h_fig,vXe,vYe,vZe,vQuat,scale_factor,step,theView,varargin)
%   INPUT:
%   vXe, vYe, vZe       vectors of CG coordinates (m), Earth axes (NED)
%   vQuat               vector of orientation's quaternion (-)
%   scale_factor        normalization factor (scalar)
%                       (related to body aircraft dimension, <1 magnifies the body)
%   step                attitude sampling factor (scalar)
%                       (the body is shown each <step> points along the flight path)
%
%   *******************************
%   Author: Agostino De Marco, Università di Napoli Federico II
%
%   Function inspired by "trajectory3" function on Matlab Central File Exchange: 
%   https://www.mathworks.com/matlabcentral/fileexchange/5656-trajectory-and-attitude-plot-version-3
%

%% Sanity checks
if (length(vXe)~=length(vYe))||(length(vXe)~=length(vZe))||(length(vYe)~=length(vZe))
    disp('  Error:');
    disp('      Wrong dimensions of CG coordinate vectors. Check the size.');
    return;
end
if size(vQuat,1)~=length(vXe)
    disp('  Error:');
    disp('      Size mismatch between quaternion history and CG coordinate vectors.');
    return
end
if step>=length(vXe)
    disp('  Error:');
    disp('      Attitude samplig factor out of range. Reduce step.');
    return
end

if step<1
    step = 1;
end

% Le coordinate x e y vengono disegnate SENZA modifiche (niente "trick"
% vXe = -vXe - min(-vXe), che spostava l'origine in x = max(vXe)).
% Sull'asse verticale si disegna la quota h = -z_E (positiva verso l'alto).
% Il passaggio (x,y,z) -> (x,y,h) e' una riflessione: per non vedere
% l'aereo "specchiato" si inverte anche l'asse Y del grafico (vedi in fondo),
% cosi' il risultato visivo e' una rotazione di 180 gradi attorno a x.
vHe = -vZe;                 % quota (m)
Tph = diag([1 1 -1]);       % da assi Earth (NED) ad assi del grafico (x,y,h)

%% Misc.
movie   = nargout;
cur_dir = pwd;

%% Read aircraft data file 'aircraft.mat' and prepare vertices
load aircraft.mat; % Coords of vertices matche with body-axes definitions

V(:,1) = V(:,1)-round(sum(V(:,1))/size(V,1));
V(:,2) = V(:,2)-round(sum(V(:,2))/size(V,1));
V(:,3) = V(:,3)-round(sum(V(:,3))/size(V,1));

Xb_nose_tip = max(abs(V(:,1)));
V = V./(scale_factor*Xb_nose_tip);

%% How many body visualizations along the flight path
usr_modulo = mod(length(vXe),step);

%% The visualization loop
frame = 0;
for i=1:step:(length(vXe)-usr_modulo)
    
    if movie || (i == 1)
        clf(h_fig);
        plot3(vXe,vYe,vHe);
        grid on;
        hold on;
        light;
    end

    % Transf. matrix from Earth- to body-axes 
    Tbe = quat2dcm(vQuat(i,:));
    % Transf. matrix from Body- to Earth-axes 
    Teb = Tbe';
            
    %% Vertices in Earth-axis coordinates
    Vb = Tph*Teb*V';
    Vb = Vb';
    
    P_on_traj = [vXe(i) vYe(i) vHe(i)];
    X0 = repmat(P_on_traj,size(Vb,1),1);
    Vb = Vb + X0;

    %% plot body-axes
    Xb = transpose( (2/scale_factor)*Tph*Teb*[1;0;0] ); % CG-to-fuselage nose
    Yb = transpose( (2/scale_factor)*Tph*Teb*[0;1;0] ); % CG-to-right wing
    Zb = transpose( (2/scale_factor)*Tph*Teb*[0;0;1] ); % Pilot's head-to-feet direction
    quiver3( ...
        P_on_traj(1),P_on_traj(2),P_on_traj(3), ...
        Xb(1),Xb(2),Xb(3), ...
        'r','linewidth',2.5 ...
    ); hold on
    quiver3( ...
        P_on_traj(1),P_on_traj(2),P_on_traj(3), ...
        Yb(1),Yb(2),Yb(3), ...
        'g','linewidth',2.5 ...
    ); hold on
    quiver3( ...
        P_on_traj(1),P_on_traj(2),P_on_traj(3), ...
        Zb(1),Zb(2),Zb(3), ...
        'b','linewidth',2.5 ...
    ); hold on
    
    %% Display aircraft shape
    p = patch('faces', F, 'vertices' ,Vb);
    set(p, 'facec', [1 0 0]);          
    set(p, 'EdgeColor','none');
    if movie || (i == 1)
        view(theView);
        axis equal;
    end
    
    if movie
        if i == 1
            ax = axis;
        else
            axis(ax);
        end
        lighting phong
        frame = frame + 1;
        M(frame) = getframe;
    end
end % end-of-for

hold on;

%% Plot the flight path
% some useful lines
plot3([min(vXe) max(vXe)]*1.1, [max(vYe) max(vYe)]*1.1, [0 0], 'color', ones(3,1)*.8);
plot3([min(vXe) max(vXe)]*1.1, [min(vYe) min(vYe)]*1.1, [0 0], 'color', ones(3,1)*.8);

plot3([max(vXe) max(vXe)]*1.1, [min(vYe) max(vYe)]*1.1, [0 0], 'color', ones(3,1)*.8);
plot3([min(vXe) min(vXe)]*1.1, [min(vYe) max(vYe)]*1.1, [0 0], 'color', ones(3,1)*.8);

% curves (ones(...) invece di vYe./vYe: con y = 0 dava NaN e la curva spariva)
plot3(vXe, vYe, ones(size(vHe))*min(vHe), '-', 'color', [.5 .5 .5]) % Ground Track
plot3(vXe, ones(size(vYe))*min(vYe), vHe, '-', 'color', [.5 .5 .5]) % Trajectory in a vertical plane
plot3(vXe, vYe, vHe, 'k', 'linewidth', 1.4); % 3D trajectory

lighting phong;
daspect([1 1 1]);

xlabel('x_E'); ylabel('y_E'); zlabel('h = -z_E');

%% Plot Earth axes
quiver3( ...
    0,0,0, ...
    abs(1.1*max(vXe)),0,0, ...
    'r','linewidth',2.5 ...
); hold on;
quiver3( ...
    0,0,0, ...
    0,abs(1.1*max(vYe)),0, ...
    'g','linewidth',2.5 ...
); hold on;
quiver3( ...
    0,0,0, ...
    0,0,-max([abs(min(vHe)),0.18*abs(max(vXe))]), ... % z_E punta verso il basso
    'b','linewidth',2.5 ...
); hold on;

xlim([min([min([-0.05*abs(max(vXe)),1.1*min(vXe)]),0]), max([1.1*max(vXe),1])]);
ylim([min([min([-0.05*abs(max(vYe)),1.1*min(vYe)]),0]), max([1.1*max(vYe),1])]);
zlim([min([min(vHe), 0]), max([1.1*max(vHe), 1])]);

% Asse Y invertito per compensare la riflessione z -> h (vedi sopra)
set(gca,'YDir','reverse');

cd (cur_dir);
end
