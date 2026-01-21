function Error = PoSMCConfirmLB(K, Pa, PH, PD, PA, PAD, St0, Alphabet, States)

% alpha is fraction of honest mining power
% K is the number of confirmation needed, Bitcoin uses 5
% Pa is the prob of (a=1, h=0)
% PH, PD, PA, PAD: see PoWSlotPdf
% St0 is the state after warmup
% St2 is the state at the time of potential confirm 

Error = 0;
SimulateSteps = 80;
Z = States+1; 

St = zeros(States+1, 2*States+1, K+1);    % tracking 3D state (L-p, S-p, K-S). See Section 4.7 of CRYPTO23 paper
for i = 0:States
    St(1+i, Z, 1) = St0(Z+i);           % L-p = i, S = p, # blocks = 0 
end

St_n = zeros(2*States+1, 2*States+1, K+1);
for step = 1:SimulateSteps
    St_n(:) = 0;
        
    % transition due to (a=1,h=0) 
    Symb_Prob = Pa;
    for l = 0:States-1
        for s = -States:States-1
            for k = 1:K                
                if s+1 >= 0 && k+1 > K
                    Error = Error + St(1+l,Z+s,k) * Pa;
                else                    
                    St_n(1+l+1,Z+s+1,k+1) = St_n(1+l+1,Z+s+1,k+1) + St(1+l,Z+s,k) * Pa;
                end
            end
        end
    end
              
    % general case        
    for h = 1:Alphabet
        for j = 1:max(1, h-1)          
            for i = 0:Alphabet-1        % symbol is (H=h, D=j, A=i) 
                Prob = PH(h) * PD(h, j) * PA{h, j}(i+1); 
                Symb_Prob = Symb_Prob + Prob;                                
                for l = 0:States                 
                    for s = -States:l     
                        for k = 1:K   
                            if s < 0 % Case 1: s < p
                                dp = j;
                                stmp = s + i - dp;
                                ltmp = max(0, l+i-dp);   % honest on top of p, adv on top of L                                                       
                            else % Case 2: s >= p
                                dp = j;
                                stmp = max(0, s+i-dp);  % honest on top o p, av on top of S      
                                ltmp = l + i - dp;
                            end                                                        
                            lnew = max(ltmp, stmp);
                            snew = min(ltmp, stmp);
                            knew = min(k + snew + dp - s, K+1); 
                            knew = max(knew, 1);
                            if lnew > States || abs(snew) > States
                                continue;
                            end
                            if snew >= 0 && knew > K                           
                                Error = Error + St(1+l,Z+s,k) * Prob;
                            else    
                                St_n(1+lnew,Z+snew,knew) = St_n(1+lnew,Z+snew,knew) + St(1+l,Z+s,k) * Prob;
                            end
                        end
                    end
                end
            end
        end
    end       

    St = St_n;
end

end


