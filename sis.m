%% Constants And Initial Conditions
T = 20; t = 0:.01:T;    % time

beta = 1;               % contact rate
b = .35;                % birth/death rate
gamma = .75;            % recovery rate
R0 = beta/(b+gamma);    % reproduction number

N = 100;                % population size
I0 = 20;                % initial number of infectives

%% Deterministic SIS Model
I = N*(1-1/R0)./(1+((1-1/R0)*N/I0-1)*exp(-(1-1/R0)*beta*t));
S = N - I;

figure(1)
plot(t,I, t,S, t,ones(size(t))*N*(1-1/R0),'--')
xlabel('Time $t$')
ylabel('Number of individuals')
legend({'$I(t)$','$S(t)$','$N \, (1-1/\mathcal{R}_0)$'},'Interpreter','latex')

%% DTMC SIS Model
dt = 0.01; t = 0:dt:T;  % time
steps = T/dt;           % number of time steps

I = zeros(steps+1,1);   % infectives
I(1) = I0;              % initial number of infectives

for k = 1:steps
    r = rand;           % random number
    if r <= beta*I(k)*(N-I(k))/N*dt
        I(k+1) = I(k) + 1;
    elseif (r > beta*I(k)*(N-I(k))/N*dt) && (r <= (beta*I(k)*(N-I(k))/N+(b+gamma)*I(k))*dt)
        I(k+1) = I(k)-1;
    elseif (r > (beta*I(k)*(N-I(k))/N+(b+gamma)*I(k))*dt) && (r <= 1)
        I(k+1) = I(k);
    end
end
S = N - I;              % susceptible individuals

figure(2)
plot(t,I, t,S, t,ones(size(t))*N*(1-1/R0),'--')
xlabel('Time $t$')
ylabel('Number of individuals')
legend({'$I(t)$', '$S(t)$', '$N \, (1-1/\mathcal{R}_0)$'},'Interpreter','latex')