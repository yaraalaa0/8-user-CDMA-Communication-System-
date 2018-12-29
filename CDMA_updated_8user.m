SNR = [-25:1:15];
L = length(SNR);
BER_BPSK_Simulation=zeros(1,L);
err =zeros(1,L);
N=500; %number of symbols transmitted
Global_counter = 0;
code1 = randi([0,1],1,32);
code2 = randi([0,1],1,32);
code3 = randi([0,1],1,32);
code4 = randi([0,1],1,32);
code5 = randi([0,1],1,32);
code6 = randi([0,1],1,32);
code7 = randi([0,1],1,32);
code8 = randi([0,1],1,32);
code1 = 2*(code1-0.5);
code2 = 2*(code2-0.5);
code3 = 2*(code3-0.5);
code4 = 2*(code4-0.5);
code5 = 2*(code5-0.5);
code6 = 2*(code6-0.5);
code7 = 2*(code7-0.5);
code8 = 2*(code8-0.5);

for si=1:length(SNR) 
while (err(si)< 600)
Global_counter =Global_counter+1;
trans = zeros(1,32*N);
X_un1 = randi([0,1],1,N);
X_un2 = randi([0,1],1,N);
X_un3 = randi([0,1],1,N);
X_un4 = randi([0,1],1,N);
X_un5 = randi([0,1],1,N);
X_un6 = randi([0,1],1,N);
X_un7 = randi([0,1],1,N);
X_un8 = randi([0,1],1,N);
X1 = 2*(X_un1-0.5);
X2 = 2*(X_un2-0.5);
X3 = 2*(X_un3-0.5);
X4 = 2*(X_un4-0.5);
X5 = 2*(X_un5-0.5);
X6 = 2*(X_un6-0.5);
X7 = 2*(X_un7-0.5);
X8 = 2*(X_un8-0.5);
counter = 0;
for i=1:N
    
    index = (counter*32)+1;
    trans(index:index+31) = (X1(i)*code1) + (X2(i)*code2) + (X3(i)*code3) + (X4(i)*code4) + (X5(i)*code5) + (X6(i)*code6) + (X7(i)*code7) + (X8(i)*code8);
    counter = counter +1;
end


N0 = 1/(10^(SNR(si)/10));
noisevec = sqrt(15*N0)*randn(size(trans));
R = trans + noisevec;

%% Integrating over the code %%

res1 = zeros(N);
res2 = zeros(N);
res3 = zeros(N);
res4 = zeros(N);
res5 = zeros(N);
res6 = zeros(N);
res7 = zeros(N);
res8 = zeros(N);
%Result = zeros(L,N);
counter = 0;
for k=1:N
    index = (counter*32)+1;
    res1(k)= sum(R(index:index+31).*code1)/32;
    res2(k)= sum(R(index:index+31).*code2)/32;
    res3(k)= sum(R(index:index+31).*code3)/32;
    res4(k)= sum(R(index:index+31).*code4)/32;
    res5(k)= sum(R(index:index+31).*code5)/32;
    res6(k)= sum(R(index:index+31).*code6)/32;
    res7(k)= sum(R(index:index+31).*code7)/32;
    res8(k)= sum(R(index:index+31).*code8)/32;
    counter = counter +1;
end

result1 = res1 > 0;
result2 = res2 > 0;
result3 = res3 > 0;
result4 = res4 > 0;
result5 = res5 > 0;
result6 = res6 > 0;
result7 = res7 > 0;
result8 = res8 > 0;
diff1 = zeros(1,N);
diff2 = zeros(1,N);
diff3 = zeros(1,N);
diff4 = zeros(1,N);
diff5 = zeros(1,N);
diff6 = zeros(1,N);
diff7 = zeros(1,N);
diff8 = zeros(1,N);
for i=1:1:N
diff1(i) = result1(i)-X_un1(i);
diff2(i) = result2(i)-X_un2(i);
diff3(i) = result3(i)-X_un3(i);
diff4(i) = result4(i)-X_un4(i);
diff5(i) = result5(i)-X_un5(i);
diff6(i) = result6(i)-X_un6(i);
diff7(i) = result7(i)-X_un7(i);
diff8(i) = result8(i)-X_un8(i);
end
err(si) = err(si) + sum(abs(diff1)) + sum(abs(diff2)) + sum(abs(diff3)) + sum(abs(diff4)) + sum(abs(diff5)) + sum(abs(diff6)) + sum(abs(diff7)) + sum(abs(diff8));


end

%% Calculating BER of BPSK of the simulation %%
BER_BPSK_Simulation(si) = err(si)/(Global_counter*N*8);
Global_counter = 0;
end 
SNR_1u = [-25:1:5];
L_1u = length(SNR_1u);
BER_BPSK_Simulation_1u=zeros(1,L_1u);
err_1u =zeros(1,L_1u);
N_1u=500; %number of symbols transmitted
Global_counter_1u = 0;
code_1u = randi([0,1],1,32);
code_1u = 2*(code_1u-0.5);

for si=1:length(SNR_1u) 
while (err_1u(si)< 200)
Global_counter_1u =Global_counter_1u+1;
trans_1u = zeros(1,32*N_1u);
X_un_1u = randi([0,1],1,N_1u);
X_1u = 2*(X_un_1u-0.5);
counter = 0;
for i=1:N_1u
    index = (counter*32)+1;
    trans_1u(index:index+31)= X_1u(i)*code_1u;
    counter = counter +1;
end


N0 = 1/(10^(SNR_1u(si)/10));
noisevec = sqrt(15*N0)*randn(size(trans_1u));
R_1u = trans_1u + noisevec;

%% Integrating over the code %%

res = zeros(N_1u);
%Result = zeros(L,N);
counter = 0;
for k=1:N_1u
    index = (counter*32)+1;
    res(k)= sum(R_1u(index:index+31).*code_1u)/32;
    counter = counter +1;
end

res2 = res > 0;
diff = zeros(1,N_1u);
for i=1:1:N_1u
diff(i) = res2(i)-X_un_1u(i);
end
err_1u(si) = err_1u(si) + sum(abs(diff));


end

%% Calculating BER of BPSK of the simulation %%
BER_BPSK_Simulation_1u(si) = err_1u(si)/(Global_counter_1u*N_1u);
Global_counter_1u = 0;
end 

theoreticalBER = 0.5*erfc(sqrt(10.^(SNR/10)));
figure(1);

semilogy(SNR_1u,BER_BPSK_Simulation_1u,'or','LineWidth',2)
hold on;
xlabel('SNR (dB)')
ylabel('BER')
title('SNR Vs BER plot for 1 user CDMA system in AWGN Channel')
figure(1);

semilogy(SNR,BER_BPSK_Simulation,'o-g','LineWidth',2)
hold on;
xlabel('SNR (dB)')
ylabel('BER')
title('SNR Vs BER plot for 8 users CDMA system in AWGN Channel')

figure(1);
semilogy(SNR,theoreticalBER,'blad-','LineWidth',2)
%plot(SNR,BER_BPSK_Simulation)

legend('1 user CDMA AWGN Simulated','8 user CDMA AWGN Simulated','CDMA AWGN Theoretical')
axis([-25 15 10^-5 0.5]);
grid on;