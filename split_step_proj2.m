clear all; %Ensure testing is not effected by past computation
clc;
tiledlayout(2,2)
figCount = 0;

%for xx = [5,10,30, 40]
%    figCount = figCount + 1;
L = 100; %Simulation space interval
N = 2000; %Number of mesh points
h = L/N; %Size of mesh spacing
x = 0:h:(L-h); %Discretized spacial interval
%tau = h^2/3;%The method is explicit so there is a stability condition
tau = 8.333333333333335e-04;
A1 = 1; %Aplitude of soliton 1
A2 = 1; %Aplitude of soliton 2
S = 0.5; %Saturation of nonlinearity
%xPos1 = 50 - xx;
%xPos2 = 50 + xx;
xPos1 = 30;  %Starting position of wave peak in x direction
xPos2 = 60; %Starting position of wave peak in x direction
%v1 = v;
%v2 = -v;
v1 = 5;  %Velocity of soliton 1
v2 = -5; %Velocity of soliton 1
psi = A1*sech(A1*(x-xPos1)).*exp(1i*v1*(x-xPos1))...%Soliton 1 Inital
      + A2*sech(A2*(x-xPos2)).*exp(1i*v2*(x-xPos2));%Soliton 2 Inital
%plot(x,abs(psi))%Plot inital configuration. Removed for testing speed
%time=200*xx; %Max time
time = 6000;
%Used for 3D plot
X = [x];
Y = [abs(psi)];
Z = [0:time];
%%%%%%%%%%%%%%%%
vec2 = 4*pi^2/L^2*[0:N/2 - 1 -N/2:-1].^2;
for ti=1:time
    psi = psi.*exp((tau*2.i*abs(psi).^2)./(1+S*sin(abs(psi).^2))); %Update psi
    psi = ifft(fft(psi).*exp(-1i*tau*vec2)); %Update psi with FFT and IFFT approxmiation
    Const = trapz(x,abs(psi).^2); %The conserved quantity
    %X = [X;x];
    Y = [Y;abs(psi)]; %Used for 3D plot
    
    %plot(x,abs(psi)), drawnow %Plot 2D updated solution for time t
end
%nexttitle
%subplot(2,2,figCount);
%[X,Z,Y] = reducem(X,Z,Y,0.000001);
h = surf(X,Z,Y) %Plot 3D, x, t, psi
set(h,'LineStyle','none')
%Const = round(Const,6,'significant');
%title("|x1 - x2|=" + 2*xx + ",   D=" + compose("%9.7f",Const))
%end