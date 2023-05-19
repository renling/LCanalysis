function St2 = PoSMCConfirmPM(K, Pa, PH, PD, PA, PAD, St0, Alphabet, States)

% alpha is fraction of honest mining power
% K is the number of confirmation needed, Bitcoin uses 5
% Pa is the prob of (a=1, h=0)
% PH, PD, PA, PAD: see PoWSlotPdf
% St0 is the state after warmup
% St2 is the state at the time of potential confirm 

Z = States+1; 

% pre-confirmation transition
St1 = zeros(2*States+1, 2*States+1, 2*K);    % tracking 3D state (adv lead, hidden lead, # blocks)
St2 = zeros(2*States+1, 1); 
for i = 0:States
    lmin = 2*K - i;
    if lmin > 0
        St1(Z+i, Z+i, lmin) = St0(Z+i);
    else
        St2(Z+i) = St0(Z+i);
    end
end

St1_n = zeros(2*States+1, 2*States+1, 2*K);
for step = 1:2*K
    St1_n(:) = 0;
    
    for i = 0:Alphabet-1        % i is #a         
        for j = 0:Alphabet-1    % j is depth of h
            Prob = PAD(i+1, j+1);

            % the positive diagonal case
            for k = 0:States                    
                dst = k+i-j;
                if abs(dst) > States
                    continue
                end
                nBlks = j + i;   
                lmin = min(nBlks, 2*K); 
                dst_x = dst;
                dst_y = max(0, dst); 
                for l = nBlks+1:2*K
                    St1_n(Z+dst_x, Z+dst_y, l-nBlks) = St1_n(Z+dst_x, Z+dst_y, l-nBlks) + St1(Z+k, Z+k, l) * Prob;    
                end                        
                St2(Z+dst_x) = St2(Z+dst_x) + sum(St1(Z+k, Z+k, 1:lmin)) * Prob;                 
            end

            % Then the 2D plane
            for x = -States:-1
                for y = 0:States
                    nBlks = j + i;
                    lmin = min(nBlks, 2*K); 
                    dst_x = x+i-j;
                    dst_y = y+i-j;
                    if dst_x >= 0
                        dst_x = dst_y;
                    end
                    dst_y = max(0, dst_y);  % y is never negative
                    if abs(dst_x) > States || abs(dst_y) > States
                        continue
                    end
                    for l = nBlks+1:2*K
                        St1_n(Z+dst_x, Z+dst_y, l-nBlks) = St1_n(Z+dst_x, Z+dst_y, l-nBlks) + St1(Z+x, Z+y, l) * Prob;    
                    end                        
                    St2(Z+dst_x) = St2(Z+dst_x) + sum(St1(Z+x, Z+y, 1:lmin)) * Prob;  
                end                          
            end
        end
    end       

    St1 = St1_n;
end

end


