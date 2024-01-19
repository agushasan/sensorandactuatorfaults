%% Code by Agus Hasan

clear;
clc;

load xArrayPlanned.mat;
load thetaArray5.mat;
load thetaArray9.mat;
load thetaArrayA5.mat;
load thetaArrayA20.mat;


%% System parameters
tf = 40;
dt = 0.01;
t  = dt:dt:tf;

Af = 1*eye(4);

A = eye(4);
B = [0 0;0 0;dt 0;0 dt];
C = eye(4);

Ab  = [A 0*eye(4);
       Af*dt*C eye(4)-Af*dt];
Bb  = [B;zeros(4,2)];
Cb  = [0 0 0 0 1 0 0 0;
      0 0 0 0 0 1 0 0;
      0 0 0 0 0 0 1 0;
      0 0 0 0 0 0 0 1];

%% Control
u = [1 0.1]';

%% Initialization
x         = [0 0 0 0 0 0 0 0]';
xhat      = [0 0 0 0 0 0 0 0]';
xArray    = [];
xhatArray = [];
thetahat  = [0;0;0;0;0;0];
thetaArray = [];
thetahatArray = [];

Pplus       = 5*eye(rank(Ab));
QF          = 0.001*eye(rank(Ab));
RF          = 0.001*eye(rank(C));
a           = 0.999;
UpsilonPlus = 1*zeros(8,6);
S           = 0.1*eye(6);
Upsilon     = 0;
lambda      = 0.992; 
theta       = [0.0;0.0;0.0;0.0;0.0;0.0];

%% Simulation

for k = 1:(tf/dt)
   
    Phi = [[0 0;0 0;-dt*u(1) 0;0 -dt*u(2)] 0*eye(4);
       zeros(4,2) Af*dt];

    if k>1000
        theta       = [0.12;0.1;0;0;0;0];
    end

    if k>2000
        theta       = [0.2;0.0;0;0;3;0.2];
    end

    if k>3000
        theta       = [0.3;0.45;0;0;4;0.5];
    end

    xArray          = [xArray x];
    xhatArray       = [xhatArray xhat];
    thetaArray      = [thetaArray theta];
    thetahatArray   = [thetahatArray thetahat];
    
    % Simulating the system
    x = Ab*x+dt*[x(3)*cos(x(4));x(3)*sin(x(4));0;0;0;0;0;0]+Bb*u+Phi*theta;
    % Taking measurement
    y = Cb*x;
    
    FX = Ab+dt*[0 0 cos(xhat(4)) -xhat(3)*sin(xhat(4)) 0 0 0 0;
                0 0 sin(xhat(4)) xhat(3)*cos(xhat(4)) 0 0 0 0;
                0 0 0 0 0 0 0 0;
                0 0 0 0 0 0 0 0;
                0 0 0 0 0 0 0 0;
                0 0 0 0 0 0 0 0;
                0 0 0 0 0 0 0 0;
                0 0 0 0 0 0 0 0];

    % Estimation using observer
     Pmin  = FX*Pplus*FX'+QF;
     Sigma = Cb*Pmin*Cb'+RF;
     KF    = Pmin*Cb'*inv(Sigma);
     Pplus = (eye(rank(Ab))-KF*Cb)*Pmin;
     
     ytilde = y-Cb*xhat;
     QF    = a*QF + (1-a)*(KF*(ytilde*ytilde')*KF');    
     RF    = a*RF + (1-a)*(ytilde*ytilde'+Cb*Pmin*Cb');
 
    Upsilon = (eye(rank(Ab))-KF*Cb)*FX*UpsilonPlus+(eye(rank(Ab))-KF*Cb)*Phi;
    Omega   = Cb*FX*UpsilonPlus+Cb*Phi;
    Lambda  = inv(lambda*Sigma+Omega*S*Omega');
    Gamma   = S*Omega'*Lambda;
    S       = (1/lambda)*S-(1/lambda)*S*Omega'*Lambda*Omega*S;
    UpsilonPlus = Upsilon;
    
    thetahat  = thetahat + Gamma*(y-Cb*xhat);
    xhat      = Ab*xhat+dt*[xhat(3)*cos(xhat(4)) xhat(3)*sin(xhat(4)) 0 0 0 0 0 0]'+Bb*u+Phi*thetahat+KF*(y-Cb*xhat)+Upsilon*Gamma*(y-Cb*xhat);
end

%thetaAf=thetahatArray;

% Plotting
figure(1)
subplot(2,2,[1 3])
plot(xArrayPlanned(1,:),xArrayPlanned(2,:),'-b','LineWidth',6)
hold on;
plot(xArray(1,:),xArray(2,:),'-r','LineWidth',6)
hold on;
plot(xhatArray(1,:),xhatArray(2,:),':g','LineWidth',6)
hold on;
plot(xhatArray(1,1),xhatArray(2,1),'or','LineWidth',6);
hold on;
plot(xhatArray(1,end),xhatArray(2,end),'ok','LineWidth',6);
hold on;
plot(xArrayPlanned(1,end),xArrayPlanned(2,end),'ok','LineWidth',6);
set(gca,'color','white','FontSize',14)
grid on;
grid minor;
ylabel('q (m)','FontSize',36)
legend('trajectory without faults','trajectory with faults','estimated trajectory with faults','start','end')
xlabel('p (m)','FontSize',36)
set(gca,'FontSize',36)
subplot(2,2,2)
plot(t,xArrayPlanned(3,:),'-b','LineWidth',6)
hold on;
plot(t,xArray(3,:),'-r','LineWidth',6)
hold on;
plot(t,xhatArray(3,:),':g','LineWidth',6)
set(gca,'color','white','FontSize',14)
grid on;
grid minor;
ylabel('U (m/s)','FontSize',36)
legend('without faults','with faults','estimated with faults')
set(gca,'FontSize',36)
subplot(2,2,4)
plot(t,xArrayPlanned(4,:),'-b','LineWidth',6)
hold on;
plot(t,xArray(4,:),'-r','LineWidth',6)
hold on;
plot(t,xhatArray(4,:),':g','LineWidth',6)
set(gca,'color','white','FontSize',14)
grid on;
grid minor;
ylabel('\chi (Rad)','FontSize',36)
xlabel('time (s)','FontSize',36)
set(gca,'FontSize',36)

figure(2)
subplot(4,2,1)
plot(t,thetaArray(1,:),'-g','LineWidth',6)
hold on;
plot(t,thetahatArray(1,:),':m','LineWidth',6)
hold on;
plot(t,thetaArrayA5(1,:),':c','LineWidth',6)
hold on;
plot(t,thetaArrayA20(1,:),':k','LineWidth',6)
set(gca,'color','white','LineWidth',3,'FontSize',14)
grid on;
grid minor;
ylabel('\theta_1^a','FontSize',36)
%legend('actual','estimated (\lambda=0.992)','estimated (\lambda=0.995)','estimated (\lambda=0.999)')
legend('actual','estimated (A_f=I)','estimated (A_f=5I)','estimated (A_f=20I)')
set(gca,'FontSize',32)
set(gca,'XTick',[])
subplot(4,2,2)
plot(t,thetaArray(1,:)-thetahatArray(1,:),':m','LineWidth',6)
hold on;
plot(t,thetaArray(1,:)-thetaArrayA5(1,:),':c','LineWidth',6)
hold on;
plot(t,thetaArray(1,:)-thetaArrayA20(1,:),':k','LineWidth',6)
set(gca,'color','white','LineWidth',3,'FontSize',14)
grid on;
grid minor;
%legend('error (\lambda=0.992)','error (\lambda=0.995)','error (\lambda=0.999)')
legend('error (A_f=I)','error (A_f=5I)','error (A_f=20I)')
set(gca,'FontSize',32)
set(gca,'XTick',[])
subplot(4,2,3)
plot(t,thetaArray(2,:),'-g','LineWidth',6)
hold on;
plot(t,thetahatArray(2,:),':m','LineWidth',6)
hold on;
plot(t,thetaArrayA5(2,:),':c','LineWidth',6)
hold on;
plot(t,thetaArrayA20(2,:),':k','LineWidth',6)
set(gca,'color','white','LineWidth',3,'FontSize',14)
grid on;
grid minor;
ylabel('\theta_2^a','FontSize',36)
set(gca,'FontSize',32)
set(gca,'XTick',[])
subplot(4,2,4)
plot(t,thetaArray(2,:)-thetahatArray(2,:),':m','LineWidth',6)
hold on;
plot(t,thetaArray(2,:)-thetaArrayA5(2,:),':c','LineWidth',6)
hold on;
plot(t,thetaArray(2,:)-thetaArrayA20(2,:),':k','LineWidth',6)
set(gca,'color','white','LineWidth',3,'FontSize',14)
grid on;
grid minor;
set(gca,'FontSize',32)
set(gca,'XTick',[])
subplot(4,2,5)
plot(t,thetaArray(5,:),'-g','LineWidth',6)
hold on;
plot(t,thetahatArray(5,:),':m','LineWidth',6)
hold on;
plot(t,thetaArrayA5(5,:),':c','LineWidth',6)
hold on;
plot(t,thetaArrayA20(5,:),':k','LineWidth',6)
set(gca,'color','white','LineWidth',3,'FontSize',14)
grid on;
grid minor;
ylabel('\theta_3^s','FontSize',36)
set(gca,'FontSize',32)
set(gca,'XTick',[])
subplot(4,2,6)
plot(t,thetaArray(5,:)-thetahatArray(5,:),':m','LineWidth',6)
hold on;
plot(t,thetaArray(5,:)-thetaArrayA5(5,:),':c','LineWidth',6)
hold on;
plot(t,thetaArray(5,:)-thetaArrayA20(5,:),':k','LineWidth',6)
set(gca,'color','white','LineWidth',3,'FontSize',14)
grid on;
grid minor;
set(gca,'FontSize',32)
set(gca,'XTick',[])
subplot(4,2,7)
plot(t,thetaArray(6,:),'-g','LineWidth',6)
hold on;
plot(t,thetahatArray(6,:),':m','LineWidth',6)
hold on;
plot(t,thetaArrayA5(6,:),':c','LineWidth',6)
hold on;
plot(t,thetaArrayA20(6,:),':k','LineWidth',6)
set(gca,'color','white','LineWidth',3,'FontSize',14)
grid on;
grid minor;
ylabel('\theta_4^s','FontSize',36)
xlabel('time (s)','FontSize',48)
set(gca,'FontSize',32)
subplot(4,2,8)
plot(t,thetaArray(6,:)-thetahatArray(6,:),':m','LineWidth',6)
hold on;
plot(t,thetaArray(6,:)-thetaArrayA5(6,:),':c','LineWidth',6)
hold on;
plot(t,thetaArray(6,:)-thetaArrayA20(6,:),':k','LineWidth',6)
set(gca,'color','white','LineWidth',3,'FontSize',14)
grid on;
grid minor;
xlabel('time (s)','FontSize',48)
set(gca,'FontSize',32)