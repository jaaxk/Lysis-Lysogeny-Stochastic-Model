%first two graphs with 50000 steps:
steps = 50000;

for (j=1:20)
%Define population vector cI RNA, cI prot, cro RNA, cro prot:
x = zeros(steps,4);
%Set initial conditions to lysogeny state
x(1,1) = 41.66;
x(1,2) = 1736.1;

%Define stress constant s based on previously acquired value
s = 6.2;

t = zeros(steps,1);
%Define how each reaction changes population vector:
R1 = [0, 1, 0, 0];
R2 = [0, -1, 0, 0];
R3 = [1, 0, 0, 0];
R4 = [-1, 0, 0, 0];
R5 = [0, 0, 0, 1];
R6 = [0, 0, 0, -1];
R7 = [0, 0, 1, 0];
R8 = [0, 0, -1, 0];

for (i=1:steps)
    v1 = 50*x(i,1);
    v2 = 1.2*x(i,2);
    v3 = 50*(1-((x(i,4)^2)/(10^2+x(i,4)^2)));
    v4 = 1.2*x(i,1);
    v5 = 50*x(i,3);
    v6 = 0.8*x(i,4);
    v7 = 50*(1-((x(i,2)^2)/(10^2+x(i,2)^2)));
    v8 = 0.8*x(i,3);
    %New rate for degradation term
    v9 = s*x(i,2);
    rates=[v1, v2, v3, v4, v5, v6, v7, v8, v9];
    vtot = sum(rates);
    tau = -log(rand())/vtot;
    t(i+1) = t(i) + tau;
    rxn = find(cumsum(rates) > (vtot*rand()) ,1);
    switch rxn
        case 1
            x(i+1,:) = x(i,:) + R1;
        case 2
            x(i+1,:) = x(i,:) + R2;
        case 3
            x(i+1,:) = x(i,:) + R3;
        case 4
            x(i+1,:) = x(i,:) + R4;
        case 5
            x(i+1,:) = x(i,:) + R5;
        case 6
            x(i+1,:) = x(i,:) + R6;
        case 7
            x(i+1,:) = x(i,:) + R7;
        case 8
            x(i+1,:) = x(i,:) + R8;
        case 9
            x(i+1,:) = x(i,:) + R2;
            
    end
    
    
end


hold on
    figure(1)
    plot(t, x(:,1), 'green')
    plot(t, x(:,2), 'blue')
    plot(t, x(:,3), 'yellow')
    plot(t, x(:,4), 'red')
    
    figure(2)
    plot(x(:,2), x(:,4), 'blue')
    plot(x(end,2),x(end,4), 'r*')
end

figure(1)
legend('cI RNA', 'cI protein', 'cro RNA', 'cro protein')
xlabel('Time (s)')
ylabel('Number of molecules')
title('Stochastic model of cro-cI gene network. S=6.2, 50,000 steps.')

figure(2)
xlabel('cI protein')
ylabel('cro protein')
title('cro vs cI protein concentration. S=6.2, 50,000 steps.')

%second two plots with 10000 steps:
steps = 100000;

for (j=1:20)
%Define population vector cI RNA, cI prot, cro RNA, cro prot:
x = zeros(steps,4);
%Set initial conditions to lysogeny state
x(1,1) = 41.66;
x(1,2) = 1736.1;

%Define stress constant s based on previously acquired value
s = 6.2;

t = zeros(steps,1);
%Define how each reaction changes population vector:
R1 = [0, 1, 0, 0];
R2 = [0, -1, 0, 0];
R3 = [1, 0, 0, 0];
R4 = [-1, 0, 0, 0];
R5 = [0, 0, 0, 1];
R6 = [0, 0, 0, -1];
R7 = [0, 0, 1, 0];
R8 = [0, 0, -1, 0];

for (i=1:steps)
    v1 = 50*x(i,1);
    v2 = 1.2*x(i,2);
    v3 = 50*(1-((x(i,4)^2)/(10^2+x(i,4)^2)));
    v4 = 1.2*x(i,1);
    v5 = 50*x(i,3);
    v6 = 0.8*x(i,4);
    v7 = 50*(1-((x(i,2)^2)/(10^2+x(i,2)^2)));
    v8 = 0.8*x(i,3);
    %New rate for degradation term
    v9 = s*x(i,2);
    rates=[v1, v2, v3, v4, v5, v6, v7, v8, v9];
    vtot = sum(rates);
    tau = -log(rand())/vtot;
    t(i+1) = t(i) + tau;
    rxn = find(cumsum(rates) > (vtot*rand()) ,1);
    switch rxn
        case 1
            x(i+1,:) = x(i,:) + R1;
        case 2
            x(i+1,:) = x(i,:) + R2;
        case 3
            x(i+1,:) = x(i,:) + R3;
        case 4
            x(i+1,:) = x(i,:) + R4;
        case 5
            x(i+1,:) = x(i,:) + R5;
        case 6
            x(i+1,:) = x(i,:) + R6;
        case 7
            x(i+1,:) = x(i,:) + R7;
        case 8
            x(i+1,:) = x(i,:) + R8;
        case 9
            x(i+1,:) = x(i,:) + R2;
            
    end
    
    
end


hold on
    figure(3)
    plot(t, x(:,1), 'green')
    hold on
    plot(t, x(:,2), 'blue')
    plot(t, x(:,3), 'yellow')
    plot(t, x(:,4), 'red')
    
    figure(4)
    plot(x(:,2), x(:,4), 'blue')
    plot(x(end,2),x(end,4), 'r*')
end

figure(3)
legend('cI RNA', 'cI protein', 'cro RNA', 'cro protein')
xlabel('Time (s)')
ylabel('Number of molecules')
title('Stochastic model of cro-cI gene network. S=6.2, 100,000 steps.')

figure(4)
xlabel('cI protein')
ylabel('cro protein')
title('cro vs cI protein concentration. S=6.2, 100,000 steps.')




