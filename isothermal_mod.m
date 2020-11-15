%% Effect of radius on isothermal heat release

%% General properties
Rho_ce = 3.16; %g/cm3 == density of cement
Rho_w = 1; %g/cm3 == density of water
v_cbw = 0.25; %gmass of chemically bound water/gcement
wg = 0.15; %gmass of physically bound water/gcement
T_initial = 293; %K == initial temp

wc_ratio = 0.5;% Water to cement ration
k_3 = (0.4/wc_ratio);

%% Time discretization
t = 500; %hrs
nt = 10000; %no of time steps
delta_t = t/nt; %time step
r = (t/delta_t)+1; %total no of time steps
t_vec = 0:delta_t:t; %time vector

%% Cement composition
C3S = 54;
C2S = 19;
C3A = 11;
C4AF = 10;
Gypsum = 6;

% C3S = 62;
% C2S = 15;
% C3A = 8;
% C4AF = 9;
% Gypsum = 6;

r0 = 9*10^-4;%radius in cm
% r0 = 6.18*10^-4;

%% Kinetic parameters
beta_1 = 1000;
beta_2 = 1000;
E_R = 5400;
beta_3 = 7500;

% beta_1 = 9;
% beta_2 = 804;
% E_R = 6000;
% beta_3 = 5942;

R3(1,1) = (((6*10^-12)*(C3S + C3A)) + (2*10^-10))*exp(-beta_1*((1/T_initial)-(1/293)));
% S3(1,1) = (0.0003*C3S + 0.0186)*exp(-beta_2*((1/T_initial)-(1/293)));
S3(1,1) = 0.3;
T3(1,1) = 7*r0/((7*10^-10 - (8*10^-12)*C2S)*(exp(-beta_3*((1/T_initial)-(1/293)))));
U3(1,1) = 0.55/(((8*10^-8)*C3S + (10^-6))*exp(-E_R*((1/T_initial)-(1/293))));

% R3(1,1) = (3.5*10^-10)*exp(-beta_1*((1/T_initial)-(1/293)));
% S3(1,1) = (2.82*10^-2)*exp(-beta_2*((1/T_initial)-(1/293)));
% T3(1,1) = r0/((3*10^-11)*(exp(-beta_3*((1/T_initial)-(1/293)))));
% U3(1,1) = 1/((2.5*10^-6)*exp(-E_R*((1/T_initial)-(1/293))));
%% Alpha initialization
alpha_mat = zeros(101, r);
dalphadt = zeros(101, r);
dQ_mat = zeros(101, r);

dalphadt(:, 1) = 0;
alpha_mat(:, 1) = 0;

%% Total heat of hydration
He = (120*C3S + 62*C2S + 100*C4AF + 320*(C3A+Gypsum))*0.01;%Cal/g

%% Experimental results


%% Integrate using explicit method

for idx = 2:r-1

        alpha_mat(1,idx) = alpha_mat(1, idx-1) + dalphadt(1, idx-1)*delta_t;%first boundary isothermal case    
        
        dalphadt(1, idx) = ((3*Rho_w/(r0*Rho_ce*(v_cbw+wg))) * power(1 - k_3 * alpha_mat(1, idx), (2.6-4*wc_ratio))) / ((1 / ((R3(1,1) * power(alpha_mat(1, idx), -1.5)) + S3(1,1) * power(alpha_mat(1, idx), 3))) + T3(1,1) / log(alpha_mat(1, idx)) + U3(1,1) * power(1 - alpha_mat(1, idx), -2/3) - T3(1,1) * power(1 - alpha_mat(1, idx), -1/3) / log(alpha_mat(1, idx)));    
        dQ_mat(1, idx) = He*dalphadt(1, idx);
        dalphadt(1, end) = dalphadt(1, end-1);
        
      
end

%% Plot
% figure(1)
% plot(t_vec(1,:),dalphadt(1,:),'r');
% xlabel('time(hrs)');
% ylabel('dalpha/dt');
% title('dalpha/dt vs time for isothermal condition');
% axis([0 50 0 0.02]);
% hold on;
% figure(2)
% plot(t_vec(1,:), alpha_mat(1,:),'g');
% xlabel('time(hrs)');
% ylabel('alpha');
% title('alpha vs time for isothermal condition');
% hold on;
figure(3)
plot(t_vec(1,:), dQ_mat(1,:), 'b');
xlabel('time(hrs)');
ylabel('Hydration heat (Cal/g/hour)');
title('Rate of isothermal heat release');
axis([0 50 0 3]);
hold on;
% scatter(Time(:,1), dQ_dt(:,1));
% legend('Predicted', 'Experimental');


