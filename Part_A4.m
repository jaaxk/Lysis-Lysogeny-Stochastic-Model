totaltime = 500;
deltat = 0.01;
numsteps = totaltime/deltat;
cIrna = zeros(numsteps,1);
crorna = zeros(numsteps,1);
cIprot = zeros(numsteps,1);
croprot = zeros(numsteps,1);

%Initial conditions of lysogeny steady state:
cIrna(1) = 41.666;
crorna(1) = 0.00020737;
cIprot(1) = 1736.1;
croprot(1) = 0.013;

for (s=4:0.01:5)
    
for (i=1:numsteps-1)
    dcIprot_dt = (50*cIrna(i)) - (1.2*cIprot(i)) - (s*cIprot(i));
    dcIrna_dt = 50*(1-((croprot(i)^2)/(10^2+croprot(i)^2))) - (1.2*cIrna(i));
    dcroprot_dt = (50*crorna(i)) - (0.8*croprot(i));
    dcrorna_dt = 50*(1-((cIprot(i)^2)/(10^2+cIprot(i)^2))) - (0.8*crorna(i));
    
    cIrna(i+1) = cIrna(i) + dcIrna_dt*deltat;
    crorna(i+1) = crorna(i) + dcrorna_dt*deltat;
    cIprot(i+1) = cIprot(i) + dcIprot_dt*deltat;
    croprot(i+1) = croprot(i) + dcroprot_dt*deltat;
end
figure(1)
plot(cIprot,croprot)
hold on
plot(cIprot(end),croprot(end),'r*')
if (croprot(end) > cIprot(end))
    s
    break
end
end

title('System going from lysogeny to lysis under stressful conditions')
xlabel('cI protein concentration')
ylabel('cro protein concentration')
