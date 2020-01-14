% Author: Jesse Sharp; Last Update: 13/11/2019
%  
% This solves problem with Payoff =
% dt*(a1*trapz(U)+a2*trapz(V.^2)+a3*trapz(y(2,:))), corresponding to  
% Fig. S6 of the Supplementary material by Sharp et al. 

%% Set-up
clear
            
Tfinal = 20;  %Specified final time
dt = 0.001; %Time-step
N = floor(Tfinal/dt+1); %Number of nodes in time discretisation
t_y = linspace(0,Tfinal,N); %Time discretisation for plotting
omega1 = 0.97; %Portion of previous iteration's control maintained when updating control u
omega2 = 0.97; %Portion of previous iteration's control maintained when updating control v
RelTol = 1e-4; %Desired relative tolerance for convergence on u and v
MaxIters = 1000; %Number of iterations to perform before giving up if convergence is not reached
U = zeros(1,N); %Initial guess for the chemotherapy control
V = zeros(1,N); %Initial guess for the stem cell transplant control

%Model parameters
ps = 0.5; %Proliferation of S
pa = 0.43; %Proliferation of A
pl = 0.27; %Proliferation of L
gs = 0.14; %Differentiation of S to A
ga = 0.44; %Differentiation of A to D
gl = 0.05; %Differentiation of L to T
ux = 0.275; %Migration of D into the blood stream
ut = 0.3; %Migration of T into the blood stream
Alpha = 0.015; %Michaelis-Menten kinetic parameter alpha
Gamma = 0.1; %Michaelis-Menten kinetic parameter gamma
kappa = 1; %Rate chemotherapy control kills progenitor blood cells relative to leukaemic cells

%Pay-off weightings
a1 = 1; %Weighting on negative impact of chemotherapy control
a2 = 1; %Weighting on cost of stem cell transplant control
a3 = 1; %Weighting on negative impact of leukaemia

%Minimum (lower) and maximum (upper) admissable controls (bounds). 
Ulower = 0; %Lower bound on chemotherapy control
Uupper = 0.1; %Upper bound on chemotherapy control
Vlower = 0; %Lower bound on stem cell transplant control
Vupper = 0.1; %Upper bound on stem cell transplant control


%Initialisation of the state. 
%Here we use the approximate steady state of the uncontrolled system for the given parameters
A(1) = 0.3255; %Initial progenitor blood cell population
L(1) = 0.3715; %Initial leukaemic stem cell population
y(:,1) = [A(1),L(1)];

%Other set-up
iterations = 0; %initialise iteration count
RelTolStoreU = []; %Initialise vector to monitor convergence of chemotherapy control
RelTolStoreV = []; %Initialise vector to monitor convergence of stem cell transplant control

%State equations
State = @(y,U,V) [gs*(1-gs/ps)+pa*y(1)*(1-(y(1)+y(2)))-ga*y(1)-kappa*U*y(1)+V;
    pl*y(2)*(1-(y(1)+y(2)))-gl*y(2)-Alpha*y(2)/(Gamma+y(2))-U*y(2)]; 

%Costate equations
Costate = @(Lambda,y,U,V) [-Lambda(1)*(pa*(1-y(1)-y(2))-pa*y(1)-ga-kappa*U)+Lambda(2)*pl*y(2); 
    -a3+Lambda(1)*pa*y(1)-Lambda(2)*(pl*(1-y(1)-y(2))-pl*y(2)-gl-Alpha/(Gamma+y(2))+Alpha*y(2)/((Gamma+y(2))^2)-U)];

 %% Forward-backward sweep
while(iterations<MaxIters)
      
uold = U;%Store control from previous iteration
vold = V; %Store control from previous iteration
i = 0; %Initialise loop variable
t = 0; %Initialise time for forward sweep

%Forward sweep using fourth-order Runge-Kutta scheme
while i < N-1
    t = t + dt;
    i=i+1;
    k1 = State(y(:,i),U(i),V(i));
    k2 = State(y(:,i)+dt*k1/2,0.5*(U(i)+U(i+1)),0.5*(V(i)+V(i+1)));
    k3 = State(y(:,i)+dt*k2/2,0.5*(U(i)+U(i+1)),0.5*(V(i)+V(i+1)));
    k4 = State(y(:,i)+dt*k3,U(i+1),V(i+1));
    y(:,i+1) = y(:,i) + (dt/6)*(k1+2*k2+2*k3+k4);
end

Lambda = zeros(2,length(y)); %Initialise Lambda (Costate)
Lambda(:,end) = [0,0]; %Apply transversality conditions to obtain final time condition on Lambda (costate)
t = Tfinal; %Initialise time for backward sweep
i = 0; %Initialise loop variable
j = N; %Initialise loop variable

%Backward sweep using fourth-order Runge-Kutta scheme 
while j > 1 
    t = t - dt;
    i = i+1;
    j = N-i;
    k1 = Costate(Lambda(:,j+1),y(:,j+1),U(j+1),V(j+1));
    k2 = Costate(Lambda(:,j+1)-dt*k1/2,0.5*(y(:,j)+y(:,j+1)),0.5*(U(j)+U(j+1)),0.5*(V(j)+V(j+1)));
    k3 = Costate(Lambda(:,j+1)-dt*k2/2,0.5*(y(:,j)+y(:,j+1)),0.5*(U(j)+U(j+1)),0.5*(V(j)+V(j+1)));
    k4 = Costate(Lambda(:,j+1)-dt*k3,y(:,j),U(j),V(j));
    Lambda(:,j) = Lambda(:,j+1) - (dt/6)*(k1+2*k2+2*k3+k4);
end

%Control updating. Current logic assumes non-negative lower bounds.
Uupdate = max(Uupper*((sign((-y(1,:).*kappa.*Lambda(1,:)-y(2,:).*Lambda(2,:)+a1))-1)/(-2)),...
    Ulower*((sign((-y(1,:).*kappa.*Lambda(1,:)-y(2,:).*Lambda(2,:)+a1))+1)/(2))); %Determine updated bang-bang chemotherapy control
U = omega1*U + (1-omega1)*Uupdate; %Actual updated control after applying relaxation to aid convergence

Vupdate = max(min(Vupper,-0.5*Lambda(1,:)/a2),Vlower); %Determine updated continuous stem cell transplant control
V = omega2*V + (1-omega2)*Vupdate; %Actual updated control after applying relaxation to aid convergence

% %To view the solution as it converges, uncomment this block (line 113 to line 126) to produce interim
% %figures
 
% box on
% line1 = plot(t_y,y(1,:),'LineWidth',2);
% hold on
% line2 = plot(t_y,y(2,:),'LineWidth',2);
% line3 = plot(t_y,U,'k--','LineWidth',2);
% line4 = plot(t_y,V,'m--','LineWidth',2);
% hL = legend([line1,line2,line3,line4],{'A','L','u*','v*'},'Location','northeast');
% ylabel('State','fontsize',18);
% xlabel('Time','fontsize',18);
% axis([0,Tfinal,0,1])
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 18)
% hold off
% drawnow 

% if either control update is identically zero after several iterations,
% update directly and flag as converged
if iterations > 150 && max(Vupdate) == 0
    V = Vupdate;
end

if iterations > 150 && max(Uupdate) == 0
    U = Uupdate;
end

%Check for convergence
RelTolTestU = RelTol*norm(U) - norm(U-uold);
RelTolStoreU = [RelTolStoreU RelTolTestU];

RelTolTestV = RelTol*norm(V) - norm(V-vold);
RelTolStoreV = [RelTolStoreV RelTolTestV];

if RelTolTestU >= 0 && RelTolTestV >= 0
    fprintf('Specified relative tolerance of %g has been met \n\r',RelTol)
    break
end  

%Update and display iteration count
iterations = iterations+1;
fprintf('Iterations completed: %d  \n\r',iterations) 
end

Payoff = dt*(a1*trapz(U)+a2*trapz(V.^2)+a3*trapz(y(2,:))); %Calculates the pay-off 

%% Plot final state with optimal control
colours = [ 

    237/255  177/255  32/255 

    20/255  42/255  140/255

]; %Define colours for plot
 
figure
set(gca, 'ColorOrder', colours);
box on
line1 = plot(t_y,y(1,:),'LineWidth',2);
hold on
line2 = plot(t_y,y(2,:),'LineWidth',2);
line3 = plot(t_y,U,'k--','LineWidth',2);
line4 = plot(t_y,V,'m--','LineWidth',2);
hL = legend([line1,line2,line3,line4],{'A','L','u*','v*'},'Location','northeast');
ylabel('State','fontsize',18);
xlabel('Time','fontsize',18);
axis([0,Tfinal,0,1])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 18)
hold off
