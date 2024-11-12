totaltime = 20;
deltat = 0.01;
numsteps = totaltime/deltat;

hold on

for (j=0:20)
    for (k=0:20)
        
        cIrna = zeros(numsteps,1);
        crorna = zeros(numsteps,1);
        cIprot = zeros(numsteps,1);
        croprot = zeros(numsteps,1);

        crorna(1) = j;
        cIrna(1) = k;
        
        for (i=1:numsteps-1)
            dcIprot_dt = (50*cIrna(i)) - (1.2*cIprot(i));
            dcIrna_dt = 50*(1-((croprot(i)^2)/(10^2+croprot(i)^2))) - (1.2*cIrna(i));
            dcroprot_dt = (50*crorna(i)) - (0.8*croprot(i));
            dcrorna_dt = 50*(1-((cIprot(i)^2)/(10^2+cIprot(i)^2))) - (0.8*crorna(i));
    
            cIrna(i+1) = cIrna(i) + dcIrna_dt*deltat;
            crorna(i+1) = crorna(i) + dcrorna_dt*deltat;
            cIprot(i+1) = cIprot(i) + dcIprot_dt*deltat;
            croprot(i+1) = croprot(i) + dcroprot_dt*deltat;
    
        end
        croprot(end)
        cIprot(end)
        figure(1)
        plot(cIprot,croprot,'k')
        plot(cIprot(end),croprot(end),'r*')
    end
end
title('Initial RNA conditions from 0 to 20')
xlabel('cI protein concentration')
ylabel('cro protein concentration')
hold off


for (j=0:500:2000)
    for (k=0:500:2000)
        
        cIrna = zeros(numsteps,1);
        crorna = zeros(numsteps,1);
        cIprot = zeros(numsteps,1);
        croprot = zeros(numsteps,1);

        crorna(1) = j;
        cIrna(1) = k;
        
        for (i=1:numsteps-1)
            dcIprot_dt = (50*cIrna(i)) - (1.2*cIprot(i));
            dcIrna_dt = 50*(1-((croprot(i)^2)/(10^2+croprot(i)^2))) - (1.2*cIrna(i));
            dcroprot_dt = (50*crorna(i)) - (0.8*croprot(i));
            dcrorna_dt = 50*(1-((cIprot(i)^2)/(10^2+cIprot(i)^2))) - (0.8*crorna(i));
    
            cIrna(i+1) = cIrna(i) + dcIrna_dt*deltat;
            crorna(i+1) = crorna(i) + dcrorna_dt*deltat;
            cIprot(i+1) = cIprot(i) + dcIprot_dt*deltat;
            croprot(i+1) = croprot(i) + dcroprot_dt*deltat;
    
        end
        figure(2)
        plot(cIprot,croprot, 'k')
        hold on
        plot(cIprot(end),croprot(end),'r*')

    end
end
title('Initial RNA conditions from 0 to 2000')
xlabel('cI protein concentration')
ylabel('cro protein concentration')
hold off


