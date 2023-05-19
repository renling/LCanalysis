function St2 = PoWMCConfirmUB(K, Pa, PH, PD, PA, St0, Alphabet, States)

% alpha is fraction of honest mining power
% K is the number of confirmation needed, Bitcoin uses 5
% Pa is the prob of (a=1, h=0)
% PH, PD, PA: see PoWSlotPdf
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
    
    % transition from k to k+1 due to (a=1,h=0) 
    for k = -States:States-1
        for l = 2:2*K
            St1_n(Z+k+1,l-1) = St1_n(Z+k+1,l-1) + St1(Z+k,l) * Pa;    
        end
        St2(Z+k+1) = St2(Z+k+1) + St1(Z+k,1) * Pa;
    end  
    
    % transition from 0 to -1 due to (a=0,h=1)
    Ph = PH(1) * PA{1,1}(1);
    for l = 2:2*K
        St1_n(Z-1,l-1) = St1_n(Z-1,l-1) + St1(Z,l) * Ph;    
    end
    St2(Z-1) = St2(Z-1) + St1(Z,1) * Ph;
    
    % general case
    for h = 1:Alphabet
        for j = 1:max(1, h-1)          
            for i = 0:Alphabet-1        % symbol is (H=h, D=j, A=i) 
                for k = -States:States
                    if abs(k+i-j) > States
                        continue
                    end
                    if k == 0 && h == 1 && i == 0   % handled earlier separately
                        continue
                    end                  
                    nBlks = h + i;
                    lmin = min(nBlks, 2*K);                    
                    Prob = PH(h) * PD(h, j) * PA{h, j}(i+1);
                    if k < -i || k > j 
                        dst = k+i-j;
                    else
                        dst = min(0,k) + i;
                    end
                    for l = nBlks+1:2*K
                        St1_n(Z+dst, l-nBlks) = St1_n(Z+dst, l-nBlks) + St1(Z+k, l) * Prob;    
                    end                                            
                    St2(Z+dst) = St2(Z+dst) + sum(St1(Z+k,1:lmin)) * Prob;                             
                end
            end
        end
    end       

    St1 = St1_n;
end

end


