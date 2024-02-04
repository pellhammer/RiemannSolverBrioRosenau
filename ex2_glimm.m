% ex2_glimm.m: 
% Minimal example of the random choice method (Glimm scheme) using a background viscosity dependent 
% Riemann solver for the Brio Rosenau model.
%
% This programm uses the matlab function riemannSolverGlimm2x2.m
%
% AUTHOR:
% Valentin Pellhammer
% Department of Mathematics and Statistics,
% University of Konstanz, 78457 Konstanz
% email adress: valentin.pellhammer@uni-konstanz.de
% homepage: http://www.math.uni-konstanz.de/~pellhammer/
%
% Date: April 2023
% Last revision: February 2024

clc;close all;clear all;

% define discretisation parameters
xmin = -1;
xmax = 1;
dt = 0.0005;
N = 400;
tmax = 1;
dx = (xmax-xmin)/(N-1);
x =  linspace(xmin, xmax,N);


% define model 
a=3;
F = @(u,v) [a* u.^2 + v.^2;2*u.*v];

BoundaeryCond = 'transmissive';
% BoundaeryCond = 'periodic';

% Initial data %
% Riemann data 
uL = 0;
vL = -1;
uR = 1/3;
vR = 1/2;
init =  [(x<0) .* uL + (x>=0) .* uR;(x<0) .* vL + (x>=0) .* vR];
% other data
    init = [cos(2.*x);sin(6.*x)];

sol = init;

% initialize plots
subplot(2,1,1)
h1 = plot(x,sol(1,:),'k*');hold on;
h1.MarkerSize = 4;

T = title('Random choice method (press any key)');
axis([xmin xmax,-1 2]);
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$u$','Interpreter','latex','Rotation',0,'FontSize',14);
subplot(2,1,2)
h2 = plot(x,sol(2,:),'k*');hold on;
h2.MarkerSize = 4;

xlabel('$x$','Interpreter','latex','FontSize',14);
axis([xmin xmax,-3 3]);
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$v$','Interpreter','latex','Rotation',0,'FontSize',14);
an = annotation('textbox',[0.8,0.07,0,0],'string','$t=0$','Interpreter','latex','FontSize',14);

pause;
T.String = 'Random choice method';

for i= 0:dt:tmax
  if strcmp(BoundaeryCond,'transmissive')
  %  transmissive boundary conditions
    solm = [sol(:,1),sol(:,1:end-1)];
    solp = [sol(:,2:end),sol(:,end)]; 
  else
%     % periodic boundary conditions
   solm = [sol(:,end),sol(:,1:end-1)];
   solp = [sol(:,2:end),sol(:,1)]; 
  end
           
    theta = ones(1,length(x)) .*rand(1,1);

    [Sm1,Sm2] = riemannSolverGlimm2x2(solm,sol,a,theta.*(dx/dt));
    
    [Sp1,Sp2] = riemannSolverGlimm2x2(sol,solp,a,(theta-1).*(dx/dt));

    sol = (theta<=(1/2)).*[Sm1;Sm2] + (theta>(1/2)) .* [Sp1;Sp2];
    
    % data update
    set(h1,'Ydata',sol(1,:));     
    set(h2,'Ydata',sol(2,:));
    an.String = ['$t=',num2str(i),'$'];
    
    
    
    pause(0.05);
end

