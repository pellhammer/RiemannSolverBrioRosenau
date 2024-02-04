% ex1_godunov.m: 
% Minimal example of the godunov method using a background viscosity dependent 
% Riemann solver for the Brio Rosenau model.
%
% This programm uses the matlab function riemannSolverGodunov2x2.m
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
dt = 0.001;
N = 350;
tmax = 1;
dx = (xmax-xmin)/(N-1);
x =  linspace(xmin, xmax,N);


% define model 
a=3;
F = @(u,v) [a* u.^2 + v.^2;2*u.*v];

% boundary conditions
BoundaeryCond = 'transmissive';
% BoundaeryCond = 'periodic';


% Initial data %
uL = -1/2;
vL = 0;
uR = 0.2;
vR = 0;
init =  [(x<0) .* uL + (x>=0) .* uR;(x<0) .* vL + (x>=0) .* vR];
% other data
init = [cos(6.*x);sin(4.*x)];
init = [3/4*x.*sin(6*x);-sin(7*(x-0.3))];
% init = [-sin(7*(x-0.3));exp(-x.*x)];





sol=init;

%initialize plots
subplot(2,1,1)
h1 = plot(x,sol(1,:),'-k*','LineWidth',1);hold on;
h1.MarkerSize = 4;
grid on;
T = title('Godunov method (press any key)');
axis([xmin xmax,-1 1.3]);
xlabel('x','Interpreter','latex');
ylabel('u','Interpreter','latex','Rotation',0,'FontSize',14);

subplot(2,1,2)
h2 = plot(x,sol(2,:),'-k*','LineWidth',1);hold on;
h2.MarkerSize = 4;
grid on;
xlabel('$x$','Interpreter','latex','FontSize',14);
axis([xmin xmax,-2 2]);
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$v$','Interpreter','latex','Rotation',0,'FontSize',14);
an = annotation('textbox',[0.8,0.07,0,0],'string','$t=0$','Interpreter','latex','FontSize',14);

pause;
T.String = 'Godunov method';
kk = 1;
 

% generate vido
% wobj = VideoWriter('test1.avi');
% wobj.FrameRate = 10;                  % frames per second (video speed)
% open(wobj); 



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
%     
    [Sp1,Sp2] = riemannSolverGodunov2x2(sol,solp,a);

    [Sm1,Sm2] = riemannSolverGodunov2x2(solm,sol,a);

    sol = sol - (dt./dx).*(F(Sp1,Sp2) - F(Sm1,Sm2));
   
    set(h1,'Ydata',sol(1,:));
    set(h2,'Ydata',sol(2,:));
    an.String = ['$t=',num2str(i),'$'];


%M(kk) = getframe; activate when making video


  pause(0.001);
     kk = kk+1;
end



