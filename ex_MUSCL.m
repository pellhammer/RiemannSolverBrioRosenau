% ex1_godunov.m: 
% Minimal example of the MUSCL using a background viscosity dependent 
% Riemann solver for the Brio Rosenau model.
%
% The following flux limiter is used: superbee
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
N = 250;
tmax = 1;
dx = (xmax-xmin)/(N-1);
x =  linspace(xmin, xmax,N);


% define model 
a=3;
F = @(u,v) [a* u.^2 + v.^2;2*u.*v];

% boundary conditions
BoundaeryCond = 'transmissive';
 %BoundaeryCond = 'periodic'; % up to now only for periodic boundary

% Initial data %
% Riemann data (transmissive boundary cond. recomended; cf. main loop)
uL = -1/2;
vL = 0;
uR = 0.2;
vR = 0;
init =  [(x<0) .* uL + (x>=0) .* uR;(x<0) .* vL + (x>=0) .* vR];
% other data
%init = [cos(6.*x);sin(4.*x)];
init = [x.*sin(8*x);-sin(7*(x-0.3))];
% init = [-sin(7*(x-0.3));exp(-x.*x)];


sol=init;

%initialize plots
subplot(2,1,1)
h1 = plot(x,sol(1,:),'-k','LineWidth',3);hold on;
grid on;
T = title('Godunov method (press any key)');
axis([xmin xmax,-1 1]);
xlabel('x','Interpreter','latex');
ylabel('u','Interpreter','latex','Rotation',0,'FontSize',14);

subplot(2,1,2)
h2 = plot(x,sol(2,:),'-k','LineWidth',3);hold on;
grid on;
xlabel('$x$','Interpreter','latex','FontSize',14);
axis([xmin xmax,-2 2]);
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$v$','Interpreter','latex','Rotation',0,'FontSize',14);
an = annotation('textbox',[0.8,0.07,0,0],'string','$t=0$','Interpreter','latex','FontSize',14);

pause;
T.String = 'Godunov method';
kk = 1;


% uncomment to generate video
% wobj = VideoWriter('test1.avi');
% wobj.FrameRate = 10;                  % frames per second (video speed)
% open(wobj); 

beta = 2;

for i= 0:dt:tmax


    DeltaBar = fluxlimiter(sol,beta);

    UL = sol - (1/2)*DeltaBar;
    UR = sol + (1/2)*DeltaBar;


    ULbar = UL + (1/2)*(dt/dx)*(F(UL(1,:),UL(2,:)) - F(UR(1,:),UR(2,:))); 
    URbar = UR + (1/2)*(dt/dx)*(F(UL(1,:),UL(2,:)) - F(UR(1,:),UR(2,:))); 


  if strcmp(BoundaeryCond,'transmissive')
  %  transmissive boundary conditions
    solm = [sol(:,1),sol(:,1:end-1)];
    solp = [sol(:,2:end),sol(:,end)]; 
  else
%     % periodic boundary conditions
   solm = [sol(:,end),sol(:,1:end-1)];
   solp = [sol(:,2:end),sol(:,1)]; 
  end



    RP1L = ULbar;
    RP1R = [URbar(:,2:end),URbar(:,1)]; 
    [Sp1,Sp2] = riemannSolverGodunov2x2(RP1L,RP1R,a);


    RP2L = [ULbar(:,end),ULbar(:,1:end-1)];
    RP2R = URbar;
    [Sm1,Sm2] = riemannSolverGodunov2x2(RP2L,RP2R,a);

    sol = sol - (dt./dx).*(F(Sp1,Sp2) - F(Sm1,Sm2));
   
    set(h1,'Ydata',sol(1,:));
    set(h2,'Ydata',sol(2,:));
    an.String = ['$t=',num2str(i),'$'];


% only for generating video
% M(kk) = getframe


  pause(0.001);
     kk = kk+1;
end



































