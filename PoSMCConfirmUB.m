function St2 = PoSMCConfirmUB(K, Pa, PH, PD, PA, St0, Alphabet, States)

% alpha is fraction of honest mining power
% K is the number of confirmation needed
% Pa is the prob of (a=1, h=0)
% PH, PD, PA: see PoSSlotPdf
% St0 is the state after warmup
% St2 is the state at the time of potential confirm 

Z = States+1; 

% pre-confirmation transition
St1 = zeros(2*States+1, 2*States+1, 2*K);    % tracking 3D state (adv lead, hidden lead, # blocks counting down)
St2 = zeros(2*States+1, 1); 
for i = 0:States
    lmin = K - i;            % is adv's lead (i) already big enough to win? lead = k wins in PoS
    if lmin > 0
        St1(Z+i, Z+i, lmin) = St0(Z+i);
    else
        St2(Z+i) = St0(Z+i);            
    end
end

St1_n = zeros(2*States+1, 2*States+1, 2*K);
for step = 1:2*K
    St1_n(:) = 0;
    
    % transition from k to k+1 due to (a=1,h=0), each adv block counts down by 2
    for k = 0:States-1
        for l = 3:2*K
            St1_n(Z+k+1,Z+k+1,l-2) = St1_n(Z+k+1,Z+k+1,l-2) + St1(Z+k,Z+k,l) * Pa;     
        end
        St2(Z+k+1) = St2(Z+k+1) + sum(St1(Z+k,Z+k,1:2)) * Pa;   % count down to 0, potential commit
    end  
    for y = 0:States-1  % x from -1 to 0, teleport case
        for l = 3:2*K
            St1_n(Z+y+1,Z+y+1,l-2) = St1_n(Z+y+1,Z+y+1,l-2) + St1(Z-1,Z+y,l) * Pa;    
        end
        St2(Z+y+1) = St2(Z+y+1) + sum(St1(Z-1,Z+y,1:2)) * Pa;
    end
    for x = -States:-2 
        for y = 0:States-1
            for l = 3:2*K
                St1_n(Z+x+1,Z+y+1,l-2) = St1_n(Z+x+1,Z+y+1,l-2) + St1(Z+x,Z+y,l) * Pa;    
            end
            St2(Z+x+1) = St2(Z+x+1) + sum(St1(Z+x,Z+y,1:2)) * Pa; 
        end
    end  
    
    % transition from 0 to -1 due to (a=0,h=1), each honest block counts down by 1
    Ph = PH(1) * PA{1,1}(1);
    for l = 2:2*K
        St1_n(Z-1,Z,l-1) = St1_n(Z-1,Z,l-1) + St1(Z,Z,l) * Ph;    
    end
    St2(Z-1) = St2(Z-1) + St1(Z,Z,1) * Ph;
           
    % general case        
    for h = 1:Alphabet
        for j = 1:max(1, h-1)          
            for i = 0:Alphabet-1        % symbol is (H=h, D=j, A=i) 
                Prob = PH(h) * PD(h, j) * PA{h, j}(i+1);
                
                % the positive diagonal case
                for k = 0:States
                    if abs(k+i-j) > States
                        continue
                    end
                    if k == 0 && h == 1 && i == 0
                        continue    % (0,0) with h=1, a=0 is already handled
                    end
                    nBlks = h + 2*i;    % each adv block counts down by 2
                    lmin = min(nBlks, 2*K);                                   
                    if k > j 
                        dst = k+i-j;
                    else
                        dst = i;
                    end
                    for l = nBlks+1:2*K 
                        St1_n(Z+dst, Z+dst, l-nBlks) = St1_n(Z+dst, Z+dst, l-nBlks) + St1(Z+k, Z+k, l) * Prob;    
                    end                        
                    St2(Z+dst) = St2(Z+dst) + sum(St1(Z+k, Z+k, 1:lmin)) * Prob;                  
                end
                
                 % Then the 2D plane
                for x = -States:-1
                    for y = 0:States
                        nBlks = h + 2*i;    % each adv block counts down by 2
                        lmin = min(nBlks, 2*K);  
                        if x < -i
                            dst_x = x+i-j;
                            dst_y = max(0, y+i-j);  % y is never negative                        
                        else
                            dist2z = -x; 
                            telePos = y - x;
                            dst_x = telePos + (i-dist2z) - j;
                            dst_x = max(0, dst_x);
                            dst_y = dst_x;
                        end
                        if min(dst_x,dst_y) < -States || max(dst_x,dst_y) > States
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
    end       

    St1 = St1_n;
end

end


