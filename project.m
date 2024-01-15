clc
clear all
 
%CONSTANTS
 hct = 0.35;
 kh = 0.2;
 kl = 1-kh;
 e = 1/3;
 fh = 0.85;
 fl = 1-fh;              % eqn 17
 Knr =3/1000;
 Quf_d = 9.1/1000;      % reference - paper no. 26
 Quf_id =0;
 alpha = 1.14/1000;      % reference - paper no. 26
 CO = 5.8;
 Qb_d = 280/1000;
 Qb_id =0;
 Qd= 500/1000;
 
 
 %estimating flowrates by sloving linear equations
 Qbp_d = (1-hct)*Qb_d;               %eqn 6
 Qs_d = (1-hct)*(CO - Qb_d);         %eqn 15
 Qhp_d = fh*Qs_d ;                   %eqn 16
 Qlp_d = fl*Qs_d;                    %equ 16
 Qbp_id = (1-hct)*Qb_id;             %eqn 6
 Qs_id = (1-hct)*(CO - Qb_id);       %eqn 15
 Qhp_id = fh*Qs_id ;                 %eqn 16
 Qlp_id = fl*Qs_id;                  %equ 16
 
 %for estimating value of constant Kd
 KoA = 800;     
 c1 = (-1)*(KoA*(Qb_d -Qd))/(Qb_d*Qd*1000); 
 c2 = exp(c1);
 Kd_d= (c2 - 1)*Qd*Qb_d/(c2*Qd-Qb_d);       % Michael`s Equation
 Kd_id =0;
 
 
 
 
%DEFINING TIME PERIODS OF BOTH THE PERIODS
 tspan_dialysis =(0:20:240);
 tspan_inter_dialysis =(240:20:480);
 
% DATA OF PATIENT 1 to 10;
 
 kip_p =[43 46 51 48 39 28 54 24 40 70]*(10^-3);
 vd_p = [18.47 10.51 14.28 15.03 15.21 16.26 14.39 10.84 13.10 14.12];
 fp_p = [0.41 0.43 0.27 0.38 0.27 0.33 0.46 0.54 0.39 0.42];
 Gb2m_p =[0.134 0.121 0.125 0.091 0.136 0.128 0.122 0.155 0.171 0.126];
 cart_0 =[32 30 28 26 24 22 20 18 17 15];
 
 %-----Analaysis for Ind Patients--------------------
 
 %ANALYSIS FOR PATIENT 1
 
%INITIAL CONDITIONS
Conc_1 = (Qhp_d+Qlp_d+Kd_d)*cart_0(1)/(Qs_d+ Quf_d);
 
x0_init = [vd_p(1)*kh*0.25 vd_p(1)*kh*0.75 vd_p(1)*kl*.25 vd_p(1)*kl*0.75 Conc_1 Conc_1 Conc_1 Conc_1];
 
%DURING DIALYSIS
 [t1_d,x1_d] = ode45(@(t1_d,x1_d)ODE_during_dialysis(x1_d,kip_p(1),vd_p(1),fp_p(1),Gb2m_p(1),CO),tspan_dialysis,x0_init);
 
%INITIAL CONDITIONS for INTER_DIALYSIS
 n = length(t1_d);
 x0_i_init = x1_d(n,:);
 
 
%INTER-DIALYSIS FOR PATIENT 1
 [t1_i,x1_i] = ode45(@(t1_i,x1_i)ODE_during_interDialysis(x1_i,kip_p(1),vd_p(1),fp_p(1),Gb2m_p(1),CO),tspan_inter_dialysis,x0_i_init);
 
 t =[t1_d;t1_i];
 C_d1 =x1_d(:,5:8);
 C_i1 =x1_i(:,5:8);
 C_1 =[C_d1;C_i1];

 subplot(3,2,3)
 plot(t,C_1)
 title('Patient 1')
 xlabel('Time (min)')
 ylabel('Chp/Clp and Chi/Cli B2M Concentration (mg/L)')
 legend('HFR Plasma' , 'HFR Interstitium', 'LFR Plasma', 'LFR Interstitium')
 
 % MAKING Cart FOR PATIENT 1
 Cart_d_1 =[];
 Cart_id_1=[];
     
 % B2M CONC OF PATIENT 1
 Cart_d_1 = ((Qhp_d+kh*Quf_d)*x1_d(:,5)+(Qlp_d+kl*Quf_d)*x1_d(:,7))/(Qhp_d+Qlp_d+Kd_d);
 Cart_id_1 = ((Qhp_id+kh*Quf_id)*x1_i(:,5)+(Qlp_id+kl*Quf_id)*x1_i(:,7))/(Qhp_id+Qlp_id+Kd_id);
 Cart_1 =[Cart_d_1;Cart_id_1];
     
 subplot(3,2,1)
 plot(t,Cart_1)
 title('Patient 1')
 xlabel('Time (min)')
 ylabel('B2M Concentration (mg/L')
     
     
%ANALYSIS FOR PATIENT 10
 
%INITIAL CONDITIONS
Conc_10 = (Qhp_d+Qlp_d+Kd_d)*cart_0(10)/(Qs_d+ Quf_d);
 
x10_init = [vd_p(10)*kh*0.25 vd_p(10)*kh*0.75 vd_p(10)*kl*.25 vd_p(10)*kl*0.75 Conc_10 Conc_10 Conc_10 Conc_10];
 
%DURING FOR PATIENT 10
[t10_d,x10_d] = ode45(@(t10_d,x10_d)ODE_during_dialysis(x10_d,kip_p(10),vd_p(10),fp_p(10),Gb2m_p(10),CO),tspan_dialysis,x10_init);
 
%INITIAL CONDITIONS for INTER_DIALYSIS
n = length(t10_d);
x0_10_init = x10_d(n,:);
 
 
%INTER-DIALYSIS FOR PATIENT 10
[t10_id,x10_id] = ode45(@(t10_id,x10_id)ODE_during_interDialysis(x10_id,kip_p(10),vd_p(10),fp_p(10),Gb2m_p(10),CO),tspan_inter_dialysis,x0_10_init);
 
C_d10 =x10_d(:,5:8);
C_i10 =x10_id(:,5:8);
C_10 =[C_d10;C_i10];
subplot(3,2,4)
plot(t,C_10)
title('Patient 10')
xlabel('Time (min)')
ylabel('Chp/Clp and Chi/Cli B2M Concentration (mg/L)')
legend('HFR Plasma' , 'HFR Interstitium', 'LFR Plasma', 'LFR Interstitium')
 
% MAKING Cart FOR PATIENT 10
Cart_d_10 =[];
Cart_id_10=[];
 
% GRAPH OF PATIENT 10
Cart_d_10 = ((Qhp_d+kh*Quf_d)*x10_d(:,5)+(Qlp_d+kl*Quf_d)*x10_d(:,7))/(Qhp_d+Qlp_d+Kd_d);
Cart_id_10 = ((Qhp_id+kh*Quf_id)*x10_id(:,5)+(Qlp_id+kl*Quf_id)*x10_id(:,7))/(Qhp_id+Qlp_id+Kd_id);
Cart_10 =[Cart_d_10;Cart_id_10];
t =[t10_d;t10_id];
subplot(3,2,2)
plot(t,Cart_10)
title('Patient 10')
xlabel('Time (min)')
ylabel('B2M Concentration (mg/L)')
 
 
% Calculation of decrease in rebound when CO is increasing
co_6 =6;
rebound_c6 =[];
 
for i=1:length(vd_p)
 
%calculation of compartmental volume using vd
% 25% of volume goes to plasma region
vd_pl = 0.25*vd_p(i); 
vd_i = 0.75*vd_p(i);
% kh fraction of plasma or interstial goes to high flow
 
% calculation of compartmental concentrations
Conc6 = (Qhp_d+Qlp_d+Kd_d)*cart_0(i)/(Qs_d+ Quf_d);
x0_d = [vd_pl*kh vd_i*kh vd_pl*kl vd_i*kl Conc6 Conc6 Conc6 Conc6];
[t_p,x_p] = ode45(@(t_p,x_p)ODE_during_dialysis(x_p,kip_p(i),vd_p(i),fp_p(i),Gb2m_p(i),co_6),tspan_dialysis,x0_d);
n_d = length(t_p);
x0_id = x_p(n_d,:);
[t_pd,x_pd] = ode45(@(t_pd,x_pd)ODE_during_interDialysis(x_pd,kip_p(i),vd_p(i),fp_p(i),Gb2m_p(i),co_6),tspan_inter_dialysis,x0_id);
% calculating rebound 
rebound_c6(i)= ((x_pd(5,5)-x_p(5,5))/(x_p(1,5)-x_p(5,5)))*100;          % Equation 24
end
% rebound_c6
 
%for co =12
co_12 =12;
rebound_c12 =[];
 
for i=1:length(vd_p)
%calculation of compartmental volume using vd
% 25% of volume goes to plasma region
vd_pl = .25*vd_p(i); 
vd_i = .75*vd_p(i);
% kh fraction of plasma or interstial goes to high flow
 
% calculation of compartmental concentrations
Conc12 = (Qhp_d+Qlp_d+Kd_d)*cart_0(i)/(Qs_d+ Quf_d);
 
x0_d = [vd_pl*kh vd_i*kh vd_pl*kl vd_i*kl Conc12 Conc12 Conc12 Conc12];
[t_p,x_p] = ode45(@(t_p,x_p)ODE_during_dialysis(x_p,kip_p(i),vd_p(i),fp_p(i),Gb2m_p(i),co_12),tspan_dialysis,x0_d);
n_d = length(t_p);
x0_id = x_p(n_d,:);
[t_pd,x_pd] = ode45(@(t_pd,x_pd)ODE_during_interDialysis(x_pd,kip_p(i),vd_p(i),fp_p(i),Gb2m_p(i),co_12),tspan_inter_dialysis,x0_id);
rebound_c12(i)= ((x_pd(5,5)-x_p(5,5))/(x_p(1,5)-x_p(5,5)))*100;
end
% rebound_c12
rebound_decrease_co = [];
for i=1:length(vd_p)
rebound_decrease_co(i) = rebound_c6(i)-rebound_c12(i);                % Equation 24
end
disp('Decrease in Rebound% due to excercise keeping Kip constant:');
disp(rebound_decrease_co);
 
%ANALYSIS FOR PATIENT 1 when CO changes from 6 to 12
 
%INITIAL CONDITIONS were calculated above
 
%DURING DIALYSIS
[tc1_6,xc1_6] = ode45(@(tc1_6,xc1_6)ODE_during_dialysis(xc1_6,kip_p(1),vd_p(1),fp_p(1),Gb2m_p(1),co_6),tspan_dialysis,x0_init);
 
%INITIAL CONDITIONS for INTER_DIALYSIS
n = length(tc1_6);
xc0_1i_6 = xc1_6(n,:);
%INTER-DIALYSIS FOR PATIENT 1 when CO =6
[tc1_i_6,xc1_i_6] = ode45(@(tc1_i_6,xc1_i_6)ODE_during_interDialysis(xc1_i_6,kip_p(1),vd_p(1),fp_p(1),Gb2m_p(1),co_6),tspan_inter_dialysis,xc0_1i_6);
%plot(t1_i,x1_i(:,5))
 
% MAKING Cart FOR PATIENT 1
Cart_d6_1 =[];
Cart_id6_1=[];
 
% B2M CONC OF PATIENT 1
Cart_d6_1 = ((Qhp_d+kh*Quf_d)*xc1_6(:,5)+(Qlp_d+kl*Quf_d)*xc1_6(:,7))/(Qhp_d+Qlp_d+Kd_d);
Cart_id6_1 = ((Qhp_id+kh*Quf_id)*xc1_i_6(:,5)+(Qlp_id+kl*Quf_id)*xc1_i_6(:,7))/(Qhp_id+Qlp_id+Kd_id);
Cart_1_6 =[Cart_d6_1;Cart_id6_1];
 
 
 
%ANALYSIS FOR PATIENT 1 when CO =12
 
%INITIAL CONDITIONS were calculated above
 
%DURING DIALYSIS
 [tc1_12,xc1_12] = ode45(@(tc1_12,xc1_12)ODE_during_dialysis(xc1_12,kip_p(1),vd_p(1),fp_p(1),Gb2m_p(1),co_12),tspan_dialysis,x0_init);
 
%INITIAL CONDITIONS for INTER_DIALYSIS
 n = length(tc1_12);
 xc0_i1_12 = xc1_12(n,:);
%INTER-DIALYSIS FOR PATIENT 1 when CO =12
 [t1_i_12,x1_i_12] = ode45(@(t1_i_12,x1_i_12)ODE_during_interDialysis(x1_i_12,kip_p(1),vd_p(1),fp_p(1),Gb2m_p(1),co_12),tspan_inter_dialysis,xc0_i1_12);
 %plot(t1_i,x1_i(:,5))
 
% MAKING Cart FOR PATIENT 1
Cart_d12_1 =[];
Cart_id12_1=[];
 
% B2M CONC OF PATIENT 1
Cart_d12_1 = ((Qhp_d+kh*Quf_d)*xc1_12(:,5)+(Qlp_d+kl*Quf_d)*xc1_12(:,7))/(Qhp_d+Qlp_d+Kd_d);
Cart_id12_1 = ((Qhp_id+kh*Quf_id)*x1_i_12(:,5)+(Qlp_id+kl*Quf_id)*x1_i_12(:,7))/(Qhp_id+Qlp_id+Kd_id);
Cart_1_12 =[Cart_d12_1;Cart_id12_1];
subplot(3,2,5)
plot(t,Cart_1_12,t,Cart_1_6,'-.red');
 title('Patient 1')
 xlabel('Time (min)')
 ylabel('B2M Concentration (mg/L)')
 legend('100% increse in CO','Nominal CO')
 
%---------------- For Patient 10
 
 %ANALYSIS FOR PATIENT 10 when CO changes from 6 to 12
 
%INITIAL CONDITIONS were calculated above
 
%DURING DIALYSIS
 [tc10_6,xc10_6] = ode45(@(tc10_6,xc10_6)ODE_during_dialysis(xc10_6,kip_p(10),vd_p(10),fp_p(10),Gb2m_p(10),co_6),tspan_dialysis,x10_init);
 
%INITIAL CONDITIONS for INTER_DIALYSIS
 n = length(tc10_6);
 xc0_10i_6 = xc10_6(n,:);
%INTER-DIALYSIS FOR PATIENT 10 when CO =6
 [tc10_i_6,xc10_i_6] = ode45(@(tc10_i_6,xc10_i_6)ODE_during_interDialysis(xc10_i_6,kip_p(10),vd_p(10),fp_p(10),Gb2m_p(10),co_6),tspan_inter_dialysis,xc0_10i_6);
 %plot(t1_i,x1_i(:,5))
 
% MAKING Cart FOR PATIENT 10
Cart_d6_10 =[];
Cart_id6_10=[];
 
% B2M CONC OF PATIENT 1
Cart_d6_10 = ((Qhp_d+kh*Quf_d)*xc10_6(:,5)+(Qlp_d+kl*Quf_d)*xc10_6(:,7))/(Qhp_d+Qlp_d+Kd_d);
Cart_id6_10 = ((Qhp_id+kh*Quf_id)*xc10_i_6(:,5)+(Qlp_id+kl*Quf_id)*xc10_i_6(:,7))/(Qhp_id+Qlp_id+Kd_id);
Cart_10_6 =[Cart_d6_10;Cart_id6_10];
 
 
 
 %ANALYSIS FOR PATIENT 10 when CO =12
 
%INITIAL CONDITIONS were calculated above
 
%DURING DIALYSIS
 [tc10_12,xc10_12] = ode45(@(tc10_12,xc10_12)ODE_during_dialysis(xc10_12,kip_p(10),vd_p(10),fp_p(10),Gb2m_p(10),co_12),tspan_dialysis,x10_init);
 
%INITIAL CONDITIONS for INTER_DIALYSIS
 n = length(tc10_12);
 xc0_i10_12 = xc10_12(n,:);
%INTER-DIALYSIS FOR PATIENT 10 when CO =12
 [t10_i_12,x10_i_12] = ode45(@(t10_i_12,x10_i_12)ODE_during_interDialysis(x10_i_12,kip_p(1),vd_p(1),fp_p(1),Gb2m_p(1),co_12),tspan_inter_dialysis,xc0_i10_12);
 %plot(t1_i,x1_i(:,5))
 
% MAKING Cart FOR PATIENT 10
Cart_d12_10 =[];
Cart_id12_10=[];
 
% B2M CONC OF PATIENT 10
Cart_d12_10 = ((Qhp_d+kh*Quf_d)*xc10_12(:,5)+(Qlp_d+kl*Quf_d)*xc10_12(:,7))/(Qhp_d+Qlp_d+Kd_d);
Cart_id12_10 = ((Qhp_id+kh*Quf_id)*x10_i_12(:,5)+(Qlp_id+kl*Quf_id)*x10_i_12(:,7))/(Qhp_id+Qlp_id+Kd_id);
Cart_10_12 =[Cart_d12_10;Cart_id12_10];
subplot(3,2,6)
plot(t,Cart_10_12,t,Cart_10_6,'-.red');
 title('Patient 10')
 xlabel('Time (min)')
 ylabel('B2M Concentration (mg/L)')
 legend('100% increse in CO','Nominal CO')
%------------- for change in Kip 
 rebound_for_kip1 =[];
 rebound_for_kip2 =[];
 for i=1:length(kip_p)
 % calculation of compartmental volume using vd
 % 25% of volume goes to plasma region
 vd_pl = 0.25*vd_p(i); 
 vd_i = 0.75*vd_p(i);
 % kh fraction of plasma or interstial goes to high flow
 
 % calculation of compartmental concentrations
 Conc_kip = (Qhp_d+Qlp_d+Kd_d)*cart_0(i)/(Qs_d+ Quf_d);
 x_initial = [vd_pl*kh vd_i*kh vd_pl*kl vd_i*kl Conc_kip Conc_kip Conc_kip Conc_kip];
 [t_k1,x_k1] = ode45(@(t_k1,x_k1)ODE_during_dialysis(x_k1,kip_p(i),vd_p(i),fp_p(i),Gb2m_p(i),CO),tspan_dialysis,x_initial);
 x0_k_id1 = x_k1(length(t_k1),:);
 [t_id_k1,x_id_k1] = ode45(@(t_id_k1,x_id_k1)ODE_during_interDialysis(x_id_k1,kip_p(i),vd_p(i),fp_p(i),Gb2m_p(i),CO),tspan_inter_dialysis,x0_k_id1);
 %calculating rebound
 rebound_for_kip1(i)= (x_id_k1(n,5)-x_k1(n,5))/(x_k1(1,5)-(x_k1(n,5)))*100;         %Equation 24
 
 [t_k2,x_k2] = ode45(@(t_k2,x_k2)ODE_during_dialysis(x_k2,kip_p(i)*1.15,vd_p(i),fp_p(i),Gb2m_p(i),CO),tspan_dialysis,x_initial);
 x0_k_id2 = x_k2(length(t_k2),:);
 [t_id_k2,x_id_k2] = ode45(@(t_id_k2,x_id_k2)ODE_during_interDialysis(x_id_k2,kip_p(i)*1.15,vd_p(i),fp_p(i),Gb2m_p(i),CO),tspan_inter_dialysis,x0_k_id2);
 % calculation rebound
 rebound_for_kip2(i)= (x_id_k2(n,5)-x_k2(n,5))/(x_k2(1,5)-(x_k2(n,5)))*100;         %Equation 24
 end
 rebound_decrease_kip = rebound_for_kip1-rebound_for_kip2;
 disp('Decrease in Rebound% due to excercise keeping CO constant:');
 disp(rebound_decrease_kip);
 
%---------------------Analysis of ind patient ends-----------------
 
 % B2M CONC OF PATIENT 1
 
 MR_240 =[];
 MR_480 =[];
 % ESTIMATION OF REMOVED TOXIN MASS
 for i = 1:length(kip_p)
 %calculation of compartmental volume using vd
 % 25% of volume goes to plasma region
 vd_pl = 0.25*vd_p(i);
 vd_i = 0.75*vd_p(i);
 % kh fraction of plasma or interstial goes to high flow
 
 % calculation of compartmental concentrations
 Conc_k = (Qhp_d+Qlp_d+Kd_d)*cart_0(i)/(Qs_d+ Quf_d);
 x0_d = [vd_pl*kh vd_i*kh vd_pl*kl vd_i*kl Conc_k Conc_k Conc_k Conc_k]; %initial condition
 [t_d,x_d] = ode45(@(t_d,x_d)ODE_during_dialysis(x_d,kip_p(i),vd_p(i),fp_p(i),Gb2m_p(i),CO),tspan_dialysis,x0_d);
 m_0_d = x_d(1,1)*x_d(1,5) + x_d(1,2)*x_d(1,6) + x_d(1,3)*x_d(1,7) + x_d(1,4)*x_d(1,8);
 n_d = length(t_d);
 m_240_d =x_d(n_d,1)*x_d(n_d,5) + x_d(n_d,2)*x_d(n_d,6) + x_d(n_d,3)*x_d(n_d,7) + x_d(n_d,4)*x_d(n_d,8);
 MR_240(i) = (m_0_d - m_240_d);                %Equation 22
 
 x0_id = x_d(n_d,:);
 [t_id,x_id] = ode45(@(t_id,x_id)ODE_during_interDialysis(x_id,kip_p(i),vd_p(i),fp_p(i),Gb2m_p(i),CO),tspan_inter_dialysis,x0_id);
 
 n = length(t_id);
 m_480_id =x_id(n,1)*x_id(n,5) + x_id(n,2)*x_id(n,6) + x_id(n,3)*x_id(n,7) + x_id(n,4)*x_id(n,8);
 MR_480(i) = (m_0_d - m_480_id);               %Equation 23
 end
 disp('Removed Toxin mass:');
 disp(MR_480);
 
%APPLYING MASS BALANCE DURING DIALYSIS AND FORMING DIFFERENTIAL EQUATIONS
function dxdt = ODE_during_dialysis(x,kip,vd,fp,Gb2M,CO)
dxdt =zeros(8,1);
 
 %CONSTANTS
 hct = 0.35;
 kh = 0.2;
 kl = 1-kh;
 e = 1/3;
 fh = 0.85;
 fl = 1-fh;                  % eqn 17
 Knr =3/1000;
 Quf = 9.1/1000;             % Referrence -paper no 26
 alpha = 1.14;               % Referrence -paper no 26
 
 %CO =5.8;
 Qb = 280/1000;
 Qd = 500/1000;
 Qbp = (1-hct)*Qb;          % eqn 6,7
 Qs = (1-hct)*(CO - Qb);    % eqn 15
 Qhp = fh*Qs ;              % eqn 16
 Qlp = fl*Qs;               % equ 16
 KoA = 800;                 % ml/min
 c1 = (-1)*(KoA*(Qb -Qd))/(Qb*Qd*1000);                              %constant for kd
 c2 = exp(c1);
 Kd = (c2 -1)*Qd*Qb/(c2*Qd-Qb);                                      %Michael`s Equation
 Cart = ((Qhp+kh*Quf)*x(5)+(Qlp+kl*Quf)*x(7))/(Qhp+Qlp+Kd);          % Equation 9
 %dVhp /dt
  % index of Vhp =1, Vhi=2, Vlp =3, Vli =4 
 dxdt(1) = -e*Quf*kh*fp;            % Equation 10
 dxdt(2) = -e*Quf*kh*(1-fp);        % Equation 11
 dxdt(3) = -e*Quf*kl*fp;            % Equation 12
 dxdt(4) = -e*Quf*kl*(1-fp);        % Equation 13
 
 
 %dCVdt - after applying chain rule we have got the following eqns
 % index of Chp =5, Chi=6, Clp =7, Cli =8
 dxdt(5) = (kh*kip*(x(6)-x(5)) + kh*Quf*(x(6)-x(5)) + Qhp*(Cart-x(5)) - kh*Knr*x(5) + kh*Gb2M*fp +x(5)*e*Quf*kh*fp )/x(1) ;     % Equation 1
 dxdt(6) = (-kh*kip*(x(6)-x(5)) - kh*Quf*x(6) + kh*Gb2M *(1-fp) + x(6)*e*Quf*kh*(1-fp))/x(2);                                   % Equation 2
 dxdt(7) = (kl*kip*(x(8)-x(7)) + kl*Quf*(x(8)-x(7)) + Qlp*(Cart-x(7)) - kl*Knr*x(7) + kl*Gb2M*fp +x(7)*e*Quf*kl*fp )/x(3) ;     % Equation 3
 dxdt(8) = (-kl*kip*(x(8)-x(7)) - kl*Quf*x(8) + kl*Gb2M *(1-fp) + x(8)*e*Quf*kl*(1-fp))/x(4);                                   % Equation 4
 
 
end
 
 
%APPLYING MASS BALANCE DURING INTER DIALYSIS PERIOD AND FORMING DIFFERENTIAL EQUATIONS
function dxdt = ODE_during_interDialysis(x,kip,vd,fp,Gb2M,co)
 dxdt =zeros(8,1);
 
 %CONSTANTS
 hct = 0.35;
 kh = 0.2;
 kl = 1-kh;
 e = 1/3;
 fh = 0.85;
 fl = 1-fh;          % eqn 17
 Knr =3/1000;        %L/min
 Quf = 0;            % from a referreed paper26
 alpha = 1.14/1000; %L/min % above paper
 
 Qb = 0;
 Qd = 500/1000;          %L/min
 Qbp = (1-hct)*Qb;       %eqn 6,7
 Qs = (1-hct)*(co - Qb); %eqn 15
 Qhp = fh*Qs ;           %eqn 16
 Qlp = fl*Qs;            %equ 16
 Kd = 0;
 Cart = ((Qhp+kh*Quf)*x(5)+(Qlp+kl*Quf)*x(7))/(Qhp+Qlp+Kd);
 %dVhp /dt
 % index of Vhp =1, Vhi=2, Vlp =3, Vli =4 
 dxdt(1) = e*alpha*kh*fp;               % Equation 19
 dxdt(2) = e*alpha*kh*(1-fp);
 dxdt(3) = e*alpha*kl*fp;
 dxdt(4) = e*alpha*kl*(1-fp);
 
 
 %dCVdt - after applying chain rule we have got the following eqns
 % index of chp =5, chi =6, clp =7, cli =8
 dxdt(5) = (kh*kip*(x(6)-x(5)) + e*alpha*kh*x(5)*(1-fp) + Qhp*(Cart-x(5)) - kh*Knr*x(5) + kh*Gb2M*fp -alpha*x(5)*e*kh*fp )/x(1) ;       % Equation 18
 dxdt(6) = (-kh*kip*(x(6)-x(5))+ e*alpha*kh*x(5)*fp + kh*Gb2M *(1-fp) -alpha*x(6)*e*kh*(1-fp))/x(2);
 dxdt(7) = (kl*kip*(x(8)-x(7)) + e*alpha*kh*x(7)*(1-fp) + Qlp*(Cart-x(7)) - kl*Knr*x(7) + kl*Gb2M*fp -alpha*x(7)*e*kl*fp )/x(3) ;
 dxdt(8) = (-kl*kip*(x(8)-x(7))+ e*alpha*kh*x(7)*fp + kl*Gb2M *(1-fp) -alpha*x(8)*e*kl*(1-fp))/x(4);
 
end