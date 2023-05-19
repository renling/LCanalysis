function [Pa, PH, PD, PA, PAD] = PoWSlotPdfUB(alpha, Delta, N)

beta = 1 - alpha;       % fraction of adv mining power
g = exp(-alpha*Delta);  % prob of lagger

Pa = beta;              % Pr(A=1, H=0)
PH = zeros(N, 1);       % Pr(H = i)
PD = zeros(N, N);       % PD(i, j) = Pr(D=j | H=i)
PA = cell(N, N);        % distribution of A (0 to N-1) given H=i, D=j
PAtest = zeros(N, N); 

NSteps = 1001;
step = Delta / NSteps;
dt = [step/2:step:Delta]';
PTExp = exppdf(dt, 1/alpha) .* step ./ (1-g);
PTExp(NSteps+1) = 0;

PAtmp0 = zeros(N, 1);
PAtmp1 = zeros(N, 1);

for i = 1:N
    PH(i) = alpha * g * (1-g)^(i-1);    % H is a simple geometric r.v.
end

% If H=1, Depth is 1, A is Poisson on Delta
PD(1, 1) = 1;                          
PA{1, 1} = poisspdf(0:N-1, beta*Delta)';

PTstart = zeros(NSteps+1, 1);
PTjump = zeros(NSteps+1, 1);
PTstart(1) = 1;
PTjump(1) = 1;
for i = 2:N    
    PTstart = conv(PTstart, PTExp);
    PTjump_new = conv(PTjump, PTExp);
    
    PDdiff0 = sum(PTjump_new(1:NSteps));     % Pr of D stays unchanged vs. increments
    if i == 2
        PDdiff0 = 1;
    end
    PD(i, 1:i) = conv(PD(i-1, 1:i-1), [PDdiff0, 1-PDdiff0]);           

    % Tailgater time distribution conditioned on no jump vs. jump
    PTG0 = zeros(NSteps, 1);
    PTG1 = zeros(NSteps, 1);
    for ds = 0:NSteps
        PTG0(1:NSteps-ds) = PTG0(1:NSteps-ds) + PTjump(ds+1) .* PTExp(1:NSteps-ds);
        PTG1(NSteps-ds+1:NSteps) = PTG1(NSteps-ds+1:NSteps) + PTjump(ds+1) .* PTExp(NSteps-ds+1:NSteps);
    end        
    dt = (1:NSteps) .* step - step/2;
    
    % A during tailgater time conditioned on no jump vs. jump
    for a = 0:N-1
        PAtmp0(a+1) = dot(PTG0, poisspdf(a, beta .* dt));
        PAtmp1(a+1) = dot(PTG1, poisspdf(a, beta .* dt));
    end 
    
    for d = 1:i-1
        PAtmp = zeros(2*N-1, 1);
        if d == 1 || d <= i-2
            PAtmp = PAtmp + conv(PA{i-1, d}, PAtmp0);
        end        
        if d > 1
            PAtmp = PAtmp + conv(PA{i-1, d-1}, PAtmp1);
        end
        PA{i, d} = PAtmp(1:N);
        PAtest(i,d) = sum(PAtmp);
    end
           
    PTjump = PTjump_new(1:NSteps+1);
    PTjump(1) = sum(PTjump_new(NSteps+1:end));
    PTjump = PTjump ./ sum(PTjump);
        % If jumped, reset last jump time
        
end

for i = 1:N
    for d = 1:i-1
        PA{i, d} = PA{i, d} ./ PD(i, d);
    end
end
    
% marginal distribution of D, A
PAD = zeros(N, N);       % PAD(i, j) = Pr(A=i-1, D=j-1)
for h = 1:N
    for j = 1:max(1, h-1)
        PAD(:,j+1) = PAD(:,j+1) + PH(h) * PD(h, j) * PA{h, j}; 
    end
end
PAD(2,1) = Pa;

end