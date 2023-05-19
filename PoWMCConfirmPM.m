function St2 = PoWMCConfirmPM(K, Pa, PH, PD, PA, PAD, St0, Alphabet, States)

% alpha is fraction of honest mining power
% K is the number of confirmation needed, Bitcoin uses 5
% Pa is the prob of (a=1, h=0)
% PH, PD, PA, PAD: see PoWSlotPdf
% St0 is the state after warmup
% St2 is the state at the time of potential confirm 

Z = States+1; 

% pre-confirmation transition
St1 = zeros(2*States+1, 2*K);    % tracking 2D state (adv lead, # blocks)
St2 = zeros(2*States+1, 1); 
for i = 0:States
    lmin = 2*K - i;
    if lmin > 0
        St1(Z+i, lmin) = St0(Z+i);
    else
        St2(Z+i) = St0(Z+i);
    end
end

St1_n = zeros(2*States+1, 2*K); 
for step = 1:2*K
    St1_n(:) = 0;
                           
    for k = -States:States
        for i = 0:Alphabet-1        % i is #a         
            for j = 0:Alphabet-1    % j is depth of h
                if abs(k+i-j) > States
                    continue
                end  
                dst = k+i-j;
                Prob = PAD(i+1, j+1);
                nBlks = j + i;
                lmin = min(nBlks, 2*K);                  
                for l = nBlks+1:2*K
                    St1_n(Z+dst, l-nBlks) = St1_n(Z+dst, l-nBlks) + St1(Z+k, l) * Prob;
                end                        
                St2(Z+dst) = St2(Z+dst) + sum(St1(Z+k,1:lmin)) * Prob;
            end        
        end
    end    
    St1 = St1_n;
end


end
