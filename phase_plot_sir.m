%% Constants and initial values
S_0=98;
I_0=2;
N=100;
gamma=1;
R0=5;
b=0.01;
beta=R0*(gamma+b);

%% Equilibrium point (Xeq,Yeq)
Xeq=N/R0;
Yeq=N*b/beta*(R0-1);

%% Phase plane (i,s)
[S,I]=meshgrid(0:100, 0:100);
Sdot = -beta*I.*S/N+b*(N-S);
Idot = beta*S.*I/N-(b+gamma)*I;

% Plotting the vector field
quiver(S,I,Sdot,Idot)
hold on
% Adding trajectories
startx = [80 80];
starty = [50 60];
streamline(S,I,Sdot,Idot,startx,starty)
grid
plot(Xeq,Yeq,'or');
str1='(S*,I*)';
text(Xeq,Yeq+0.090,str1);
xlabel('number of susceptible');
ylabel('number of infected');