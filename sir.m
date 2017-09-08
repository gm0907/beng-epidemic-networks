%% Constants And Initial Conditions
t0 = 0; T = 20;         % time

beta = .9;               % contact rate
b = .25;                % birth/death rate
gamma = .75;              % recovery rate
R0 = beta/(b+gamma);    % reproduction number

N = 100;                % population size
y0 = [90 10 0]';         % initial number of individuals [S0 I0 R0]

[t,y] = ode45(@(t,y) solver(t,y,beta,b,gamma,N),[t0 T],y0);

figure(1)
plot(t,y, t,ones(size(t))*b*N/(b+gamma)*(1-1/R0),'--')
xlabel('Time $t$')
ylabel('Number of individuals')
legend({'$I(t)$', '$S(t)$', '$R(t)$', '$Nb \, (1-1/\mathcal{R}_0)/(b+\gamma)$'},'Interpreter','latex')

%% Phase Plot
figure(2)
plot(y(:,2),y(:,1))
xlabel('Number of infectious individuals $I$')
ylabel('Number of susceptible individuals $S$')
grid on

axes('position',[.18 .45 .2 .2])
box on
plot(y(:,2),y(:,1))
title('Zoom')
axis([0 2 15 25])

%% DTMC SIS Model
dt = 0.01; t = 0:dt:T;  % time
steps = T/dt;           % number of time steps

S = zeros(steps+1,1);   % susceptible individuals
I = zeros(steps+1,1);   % infectives

S(1) = y0(1);           % initial number of susceptible individuals
I(1) = y0(2);           % initial number of infectives

for k = 1:steps
    r = rand;
    if r <= beta*I(k)*S(k)/N*dt
        S(k+1) = S(k)-1;
        I(k+1) = I(k)+1;
    elseif (r > beta*I(k)*S(k)/N*dt) && (r <= (beta*I(k)*S(k)/N+gamma*I(k))*dt)
        S(k+1) = S(k);
        I(k+1) = I(k)-1;
    elseif (r > (beta*I(k)*S(k)/N+gamma*I(k))*dt) && (r <= (beta*I(k)*S(k)/N+(gamma+b)*I(k))*dt)
        S(k+1) = S(k)+1;
        I(k+1) = I(k)-1;
    elseif (r > (beta*I(k)*S(k)/N+(gamma+b)*I(k))*dt) && (r <= (beta*I(k)*S(k)/N+gamma*I(k)+b*(N-S(k)))*dt)
    	S(k+1) = S(k)+1;
    	I(k+1) = I(k);
    elseif (r > (beta*I(k)*S(k)/N+gamma*I(k)+b*(N-S(k)))*dt) && (r <= 1)
    	S(k+1) = S(k);
    	I(k+1) = I(k);
    end
end
R = N - I - S;          % removed individuals

figure(3)
plot(t,[S I R], t,ones(size(t))*b*N/(b+gamma)*(1-1/R0),'--')
xlabel('Time $t$')
ylabel('Number of individuals')
legend({'$I(t)$', '$S(t)$', '$R(t)$', '$Nb \, (1-1/\mathcal{R}_0)/(b+\gamma)$'},'Interpreter','latex')

%% Function For The Integration Of The ODE System
function yp = solver(~,y,beta,b,gamma,N)
    S = y(1); I = y(2); R = y(3);
    dS = -beta*S*I/N + b*(I+R);
    dI = beta*S*I/N - (gamma+b)*I;
    dR = gamma*I - b*R;
    yp = [dS; dI; dR];
end