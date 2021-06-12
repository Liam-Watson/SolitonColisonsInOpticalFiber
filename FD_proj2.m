clear all; %Ensure testing is not effected by past computation
clc;
tiledlayout(2,2)
figCount = 0;

for v = [0.5, 1, 2 , 10]
    figCount = figCount + 1;
L = 100; %Overall simulation interval
N = 2000; %Number of space subdivisions
h = L/N; %Size of mesh spacing
x = [0:h:(L-h)]'; %Space mesh
tau = h^2/3; %Time step

A1 = 1; %Inital amplitude of soliton 1
A2 = 1; %Intial amplitude of soliton 2
S = 0.0; %Nonlinarity saturation
xPos1 = 30; %Starting position of wave peak in x direction
xPos2 = 60; %Starting position of wave peak in x direction
v1 = v;
v2 = -v;
%v1 = 1; %Velocity of soliton 1
%v2 = -1; %%Velocity of soliton 2
psi = A1*sech(A1*(x-xPos1)).*exp(1i*v1*(x-xPos1))...%Soliton 1 Inital
      + A2*sech(A2*(x-xPos2)).*exp(1i*v2*(x-xPos2));%Soliton 2 Inital

%plot(x,abs(psi)) %Plot inital configuration. Removed for testing speed

s = tau/h^2; %code simplifcation
s2 = s/2; %code simplifcation

A = diag((1i-s)*ones(N,1),0)... %Set the FD matrix values on diagonals
    +diag(s2*ones(N-1,1),1)...
    +diag(s2*ones(N-1,1),-1);
A(1,N) = s2; %Set periodic BC elements in matrix A
A(N,1) = s2;
A = inv(A);
rhs = zeros(size(x)); %Create RHS vector
time = 6000/v; %Max time

X = [x]; %Used for 3D plot
Y = [abs(psi)];%Used for 3D plot
Z = [0:time];%Used for 3D plot
for t=1:time %Loop until max time reached   
    rhs(2:N-1) = (1i + s)*psi(2:N-1)-s2*psi(1:N-2)... %Update RHS
        - s2*psi(3:N)...
        - (2*tau*(abs(psi(2:N-1)).^2).*psi(2:N-1))./(1+S*sin(abs(psi(2:N-1)).^2));
        %- (tau*(abs(psi(1:N-2)).^2).*psi(1:N-2))./(1+S*sin(abs(psi(1:N-2))))...
        %- (tau*(abs(psi(3:N)).^2).*psi(3:N))./(1+S*sin(abs(psi(3:N))));
        
    
    rhs(1) = (1i + s)*psi(1) -s2*psi(2) -s2*psi(N)...
        -(2*tau*abs(psi(N)).^2.*psi(N))/(1+S*sin(abs(psi(N)).^2));
        %-(tau*abs(psi(N)).^2.*psi(N))/(1+S*sin(abs(psi(N))))...
        %-(tau*abs(psi(2)).^2.*psi(2))/(1+S*sin(abs(psi(2))));
        
    
    rhs(N) = (1i + s)*psi(N) -s2*psi(1) -s2*psi(N-1)...
            -(2*tau*abs(psi(N)).^2.*psi(N))/(1+S*sin(abs(psi(N)).^2));
            %-(tau*abs(psi(N-1)).^2.*psi(N-1))/(1+S*sin(abs(psi(N-1))))...
            %-(tau*abs(psi(1)).^2.*psi(1))/(1+S*sin(abs(psi(1))));
    %psi = A\rhs;
    psi= A*rhs; %Solve for psi
    Y = [abs(psi) Y]; %%Used for 3D plot
    Const = trapz(x,abs(psi).^2); %The conserved quantity
    
    %plot(x,abs(psi)), drawnow %Plot 2D updated solution for time t
end
%subplot(2,2,figCount);
figure(figCount);
h = surf(X,Z,Y') %Plot 3D, x, t, psi
set(h,'LineStyle','none')
title("|v|=" + v + ",   D=" + compose("%9.7f",Const))
pause;
end