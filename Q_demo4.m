%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; 
close all; 
clc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for solving the problem in the interval 0 < x < L
potential =5; 
L = 100; % Interval Length
N = 400; % No of points
x = linspace(0,L,N)'; % Coordinate vector
dx = x(2) - x(1); % Coordinate step
% Parameters for making intial momentum space wavefunction phi(k)
ko = 4; % Peak momentum
a = 10; % Momentum width parameter
dk = 2*pi/L; % Momentum step
km=N*dk; % Momentum limit
k=linspace(0,+km,N)'; % Momentum vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make psi(x,0) from Gaussian kspace wavefunction phi(k) using
% fast fourier transform :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = exp(-a*(k-ko).^2).*exp(-1i*6*k.^2); % unnormalized phi(k)
psi = ifft(phi); % multiplies phi by expikx and integrates vs. x
psi = psi/sqrt(psi'*psi*dx); % normalize the psi(x,0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expectation value of energy; e.g. for the parameters
% chosen above <E> = 2.062.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avgE = phi'*0.5*diag(k.^2,0)*phi*dk/(phi'*phi*dk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE POTENTIAL U(X): Either U = 0 OR
% U = step potential of height avgE that is located at x=L/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch potential 
    case 1 
    U = 0*heaviside(x-(L/2)); % free particle wave packet evolution
    case 2 
    U = avgE*heaviside(x-(L/2)); % scattering off step potential
    case 3 
        U = avgE*100* (heaviside(x+(L/2)) + heaviside(x-(L/2)));
    case 4
        w = L/20; a = 2*w;
    U = avgE* 100*( heaviside(x+w-a) - heaviside(x-w-a) ...
    + heaviside(x+w+a) - heaviside(x-w+a));

    case 5
        U = avgE .* [2*ones(floor(N/4)-L/2, 1); L*ones(floor(N/2), 1); 2*ones(floor(N/4)+L/2,1)];
        
        

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite-difference representation of Laplacian and Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e1 = zeros(N, N); 
e = ones(N,1); Lap = spdiags([e -2*e e],[-1 0 1],e1)/dx^2;
H = -(1/2)*Lap + spdiags(U,0,e1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for computing psi(x,T) at different times 0 < T < TF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NT = 200; % No. of time steps
TF = 29; 
T = linspace(0,TF,NT); % Time vector
dT = T(2)-T(1); % Time step
hbar = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time displacement operator E=exp(-iHdT/hbar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = expm(-1i*full(H)*dT/hbar); % time desplacement operator

Usc = U/max(abs(U)); % rescale U for plotting
%***************************************************************
% Simulate rho(x,T) and plot for each T
%***************************************************************
pdb_state= zeros(N, NT); 
for t = 1:NT % time index for loop
% calculate probability density rho(x,T)
psi = E*psi; % calculate new psi from old psi
rho = conj(psi).*psi; % rho(x,T)
pdb_state(:, t) = rho; 
plot(x,rho,x, Usc, '--k', 'Linewidth', 2.5); % plot rho(x,T) vs. x
axis([0 L 0 0.15]); % set x,y axis parameters for plotting
xlabel('x [m]', 'FontSize', 24);
ylabel('probability density [1/m]','FontSize', 24);
pause(0.05); % pause between each frame displayed
end
% Calculate Reflection probability
R= 0; 
for a=1:N/2
R=R+rho(a);
RR(a) = R; 
end
R=RR*dx;

close; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PLOTTING FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color', 'w', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); 

subplot(131);
set(gca, 'Fontsize', 12); 
plot(x,rho,x, Usc, '--k', 'Linewidth', 2.5);
axis([0 L 0 0.15]); % set x,y axis parameters for plotting
xlabel('x [m]', 'FontSize', 24);
ylabel('probability density [1/m]','FontSize', 24);
lgnd_str = ['T = ',num2str(T(t))];
legend(lgnd_str); 
xlabel('x [m]', 'FontSize', 24);
ylabel('probability density [1/m]','FontSize', 24);
title('Potential vs superposition of states'); 

subplot(132);
set(gca, 'Fontsize', 12); 
plot(1:N/2,R, 'Linewidth', 2.5);
xlabel('x [m]', 'FontSize', 24);
ylabel('Reflectance ','FontSize', 24);
lgnd_str = ['T = ',num2str(T(t))];
legend(lgnd_str); 
xlabel('x [m]', 'FontSize', 24);
ylabel('probability density [1/m]','FontSize', 24);
title('Potential vs superposition of states'); 



subplot(133);
set(gca, 'Fontsize', 12); 
mesh(x* ones(1, NT), ones(N,  1)*T, pdb_state); 
shading interp; 
axis tight; 
xlabel('x [m]', 'FontSize', 24);
ylabel('time [s]','FontSize', 24);
zlabel('probability density [1/m]','FontSize', 24);
title('Potential vs superposition of states'); 
colorbar; 
colormap('jet');
view(0, 90); 
