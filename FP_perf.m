clear all; %Ensure testing is not effected by past computation
clc;
format long
tic
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

v1 = 1; %Velocity of soliton 1
v2 = -1; %%Velocity of soliton 2
psi = A1*sech(A1*(x-xPos1)).*exp(1i*v1*(x-xPos1))...%Soliton 1 Inital
      + A2*sech(A2*(x-xPos2)).*exp(1i*v2*(x-xPos2));%Soliton 2 Inital



s = tau/h^2; %code simplifcation
s2 = s/2; %code simplifcation

A = diag((1i-s)*ones(N,1),0)... %Set the FD matrix values on diagonals
    +diag(s2*ones(N-1,1),1)...
    +diag(s2*ones(N-1,1),-1);
A(1,N) = s2; %Set periodic BC elements in matrix A
A(N,1) = s2;
%A = inv(A);
rhs = zeros(size(x)); %Create RHS vector
time = 6000; %Max time

for t=1:time %Loop until max time reached   
    rhs(2:N-1) = (1i + s)*psi(2:N-1)-s2*psi(1:N-2)... %Update RHS
        - s2*psi(3:N)...
        - (2*tau*(abs(psi(2:N-1)).^2).*psi(2:N-1))./(1+S*sin(abs(psi(2:N-1)).^2));
        
    
    rhs(1) = (1i + s)*psi(1) -s2*psi(2) -s2*psi(N)...
        -(2*tau*abs(psi(N)).^2.*psi(N))/(1+S*sin(abs(psi(N)).^2));
        
    
    rhs(N) = (1i + s)*psi(N) -s2*psi(1) -s2*psi(N-1)...
            -(2*tau*abs(psi(N)).^2.*psi(N))/(1+S*sin(abs(psi(N)).^2));

    psi = A\rhs;
    %psi= A*rhs; %Solve for psi
    %Const = trapz(x,abs(psi).^2) %The conserved quantity
    
end
toc
