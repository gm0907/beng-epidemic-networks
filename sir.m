function start
model_title = 'SIR Epidemics';

%% Constants and initial conditions
param.gamma = 0.75;
param.b = 0.25;
param.beta = 3;
R0= param.beta/(param.gamma+param.b);
initial.S = 98;
initial.I = 2;
initial.R = 0;
N=initial.S + initial.I + initial.R;

end_time = 20;

%% Deterministic SIR model
function deriv = ode_system (t, x, param)
    S = x(1);
    I = x(2);
    R = x(3);
    dS = -param.beta*S*I/N+param.b*(I+R);
    dI = param.beta*S*I/N-(param.gamma+param.b)*I; dR = param.gamma*I-param.b*R;
    deriv = [dS; dI; dR];
end

% Extract initial values from the 'initial' structure and collect them
% in a column vector for use in 'ode45'.
initial_values = [];
variable_names = fieldnames(initial);
for i=1:length(variable_names)
    initial_values = [initial_values; initial.(variable_names{i})];
end

% integrate the ODE system
[t, y] = ode45(@(t, x) ode_system(t, x, param), ... [0 end_time], ...
               initial_values, ...
               []);

%% DTMC SIR model
dt=0.01;
t1=0:dt:end_time;
n=end_time/dt;
S2=zeros(n+1);
I2=zeros(n+1);

I2(1)=initial.I;
S2(1)=initial.S;

for k=1:n
    u=rand;
    if u<=param.beta*I2(k)*S2(k)/N*dt
        S2(k+1)=S2(k)-1;
        I2(k+1)=I2(k)+1;
    elseif u>param.beta*I2(k)*S2(k)/N*dt && u<=param.beta*I2(k)*S2(k)/N*dt+...
            param.gamma*I2(k)*dt
        S2(k+1)=S2(k);
        I2(k+1)=I2(k)-1;
    elseif u>param.beta*I2(k)*S2(k)/N*dt+param.gamma*I2(k)*dt && u<=param.beta*...
            I2(k)*S2(k)/N*dt+(param.gamma+param.b)*I2(k)*dt
        S2(k+1)=S2(k)+1;
        I2(k+1)=I2(k)-1;
    elseif u>param.beta*I2(k)*S2(k)/N*dt+(param.gamma+param.b)*I2(k)*dt &&...
                u<=param.beta*I2(k)*S2(k)/N*dt+param.gamma*I2(k)*dt+param.b*(N-S2(k))*dt
    	S2(k+1)=S2(k)+1;
    	I2(k+1)=I2(k);
    elseif u>param.beta*I2(k)*S2(k)/N*dt+param.gamma*I2(k)*dt+param.b*...
                    (N-S2(k))*dt && u<=1
    	S2(k+1)=S2(k);
    	I2(k+1)=I2(k);
    end
end
    
R2=N-I2-S2;

%% Plotting
% prepare legend texts
legend_texts = cell(length(variable_names), 1);
for i=1:length(variable_names)
    text = [variable_names{i} '(t)'];
    legend_texts{i} = text;
end

% plot the results
figure(1);
subplot(2,1,1)
plot(t, y);
xlabel('time');
ylabel('number of individuals');
title(model_title);
if R0>=1
    hold on
    plot(t,ones(size(t))*N*param.b/(param.b+param.gamma)*(1-1/R0),'-.c')
elseif R0<=1
    hold on
    plot(t,zeros(size(t)),'-.c')
end
legend(legend_texts, 'N*b/(b+gamma)*(1-1/R0)');
subplot(2,1,2)
plot(t1,S2, t1,I2, t1,R2);
xlabel('time');
ylabel('number of individuals');
if R0>=1
    hold on
    plot(t,ones(size(t))*N*param.b/(param.b+param.gamma)*(1-1/R0),'-.c')
elseif R0<=1
    hold on
    plot(t,zeros(size(t)),'-.c')
end
end