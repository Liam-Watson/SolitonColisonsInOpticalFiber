clear all; %Ensure testing is not effected by past computation
clc;
format long
tic
L = 100; %Simulation space interval
N = 150; %Number of mesh points
h = L/N; %Size of mesh spacing
x = 0:h:(L-h); %Discretized spacial interval
tau = 8.333333333333335e-04;
%tau = h^2/3;%The method is explicit so there is a stability condition
A1 = 1; %Aplitude of soliton 1
A2 = 1; %Aplitude of soliton 2
xPos1 = 30;
xPos2 = 60;
S = 0.0; %Saturation of nonlinearity
v1 = 1;  %Velocity of soliton 1
v2 = -1; %Velocity of soliton 1
psi = A1*sech(A1*(x-xPos1)).*exp(1i*v1*(x-xPos1))...%Soliton 1 Inital
      + A2*sech(A2*(x-xPos2)).*exp(1i*v2*(x-xPos2));%Soliton 2 Inital
%plot(x,abs(psi))%Plot inital configuration. Removed for testing speed
time=6000; %Max time

vec2 = 4*pi^2/L^2*[0:N/2 - 1 -N/2:-1].^2;
for ti=1:time
    psi = psi.*exp((tau*2.i*abs(psi).^2)./(1+S*sin(abs(psi).^2))); %Update psi
    psi = ifft(fft(psi).*exp(-1i*tau*vec2)); %Update psi with FFT and IFFT approxmiation
    %Const = trapz(x,abs(psi).^2) %The conserved quantity
    %plot(x,abs(psi)), drawnow
end
toc
%nexttitle
%subplot(2,2,figCount);
%[X,Z,Y] = reducem(X,Z,Y,0.000001);
%h = surf(X,Z,Y) %Plot 3D, x, t, psi
%set(h,'LineStyle','none')