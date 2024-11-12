totaltime = 20;
deltat = 0.01;
numsteps = totaltime/deltat;
t = zeros(numsteps, 1);
cIrna = zeros(numsteps,1);
crorna = zeros(numsteps,1);
cIprot = zeros(numsteps,1);
croprot = zeros(numsteps,1);

for (i=1:numsteps-1)
    dcIprot_dt = (50*cIrna(i)) - (1.2*cIprot(i));
    dcIrna_dt = 50*(1-((croprot(i)^2)/(10+croprot(i)^2))) - (1.2*cIrna(i));
    dcroprot_dt = (50*crorna(i)) - (0.8*croprot(i));
    dcrorna_dt = 50*(1-((cIprot(i)^2)/(10+cIprot(i)^2))) - (0.8*crorna(i));
    
    cIrna(i+1) = cIrna(i) + dcIrna_dt*deltat;
    crorna(i+1) = crorna(i) + dcrorna_dt*deltat;
    cIprot(i+1) = cIprot(i) + dcIprot_dt*deltat;
    croprot(i+1) = croprot(i) + dcroprot_dt*deltat;
    
    t(i+1) = t(i) + deltat;
end
figure(1)
plot(t,cIrna, 'DisplayName','cI RNA')
hold on
plot(t, crorna, 'DisplayName', 'cro RNA')
plot(t, cIprot, 'DisplayName', 'cI protein')
plot(t, croprot, 'DisplayName', 'cro protein')
legend()
xlabel('Time (s)')
ylabel('Concentration (molecules/cell)')
title('Concentration of each species vs time. Starting conditions 0.')
hold off

cIrna = zeros(numsteps,1);
crorna = zeros(numsteps,1);
cIprot = zeros(numsteps,1);
croprot = zeros(numsteps,1);

crorna(1) = 20;

for (i=1:numsteps-1)
    dcIprot_dt = (50*cIrna(i)) - (1.2*cIprot(i));
    dcIrna_dt = 50*(1-((croprot(i)^2)/(10+croprot(i)^2))) - (1.2*cIrna(i));
    dcroprot_dt = (50*crorna(i)) - (0.8*croprot(i));
    dcrorna_dt = 50*(1-((cIprot(i)^2)/(10+cIprot(i)^2))) - (0.8*crorna(i));
    
    cIrna(i+1) = cIrna(i) + dcIrna_dt*deltat;
    crorna(i+1) = crorna(i) + dcrorna_dt*deltat;
    cIprot(i+1) = cIprot(i) + dcIprot_dt*deltat;
    croprot(i+1) = croprot(i) + dcroprot_dt*deltat;
    
    t(i+1) = t(i) + deltat;
end
figure(2)
plot(t,cIrna, 'DisplayName','cI RNA')
hold on
plot(t, crorna, 'DisplayName', 'cro RNA')
plot(t, cIprot, 'DisplayName', 'cI protein')
plot(t, croprot, 'DisplayName', 'cro protein')
legend()
xlabel('Time (s)')
ylabel('Concentration (molecules/cell)')
title('Concentration of each species vs time. Starting conditions cro RNA = 20.')
hold off

cIrna = zeros(numsteps,1);
crorna = zeros(numsteps,1);
cIprot = zeros(numsteps,1);
croprot = zeros(numsteps,1);

cIrna(1) = 50;

for (i=1:numsteps-1)
    dcIprot_dt = (50*cIrna(i)) - (1.2*cIprot(i));
    dcIrna_dt = 50*(1-((croprot(i)^2)/(10^2+croprot(i)^2))) - (1.2*cIrna(i));
    dcroprot_dt = (50*crorna(i)) - (0.8*croprot(i));
    dcrorna_dt = 50*(1-((cIprot(i)^2)/(10^2+cIprot(i)^2))) - (0.8*crorna(i));
    
    cIrna(i+1) = cIrna(i) + dcIrna_dt*deltat;
    crorna(i+1) = crorna(i) + dcrorna_dt*deltat;
    cIprot(i+1) = cIprot(i) + dcIprot_dt*deltat;
    croprot(i+1) = croprot(i) + dcroprot_dt*deltat;
    
    t(i+1) = t(i) + deltat;
end
figure(3)
plot(t,cIrna, 'DisplayName','cI RNA')
hold on
plot(t, crorna, 'DisplayName', 'cro RNA')
plot(t, cIprot, 'DisplayName', 'cI protein')
plot(t, croprot, 'DisplayName', 'cro protein')
legend()
xlabel('Time (s)')
ylabel('Concentration (molecules/cell)')
title('Concentration of each species vs time. Starting conditions cI RNA = 50.')
hold off
