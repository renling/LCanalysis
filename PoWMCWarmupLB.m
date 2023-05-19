function St0 = PoWMCWarmupLB(PAD, Alphabet, States)

% PAD(i,j) = Pr(A=i-1, D=j-1)

WarmUpSteps = 200;

Z = States+1;           
St0 = zeros(2*States+1, 1);
St0(Z) = 1;             % This corresponds to state 0

% warm up transition
Q0 = zeros(2*States+1, 2*States+1);
for i = 0:Alphabet-1        % i is #a
    for j = 0:Alphabet-1    % j is depth of h
        for k = 0:States
            if abs(k+i-j) > States
                continue
            end
            dst = max(0, k+i-j);          
            Q0(Z+dst, Z+k) = Q0(Z+dst, Z+k) + PAD(i+1,j+1);
        end        
    end
end

St0 = (Q0^WarmUpSteps) * St0;

end