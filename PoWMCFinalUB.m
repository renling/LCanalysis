function Error = PoWMCFinalUB(PAD, St2, Alphabet, States)

% PAD(i,j) = Pr(A=i-1, D=j-1)
% St2 is the state at the time of potential confirm 

SimulateSteps = 200;
Z = States+1; 

% post-confirmation, never transitions from >=0 to <0
Q2 = zeros(2*States+1, 2*States+1);
for i = 0:Alphabet-1
    for j = 0:Alphabet-1
        for k = -States:States
            if abs(k+i-j) > States
                continue
            end              
            if k < -i || k > j  
                dst = k+i-j;
                Q2(Z+dst, Z+k) = Q2(Z+dst, Z+k) + PAD(i+1,j+1);
            else    % k must be small in this case
                dst = min(0,k) + i;
                Q2(Z+dst, Z+k) = Q2(Z+dst, Z+k) + PAD(i+1,j+1);
                    % Upon h=1 and a=0, would have transitioned from 0 to -1
                    % but stayed at 0 because attack already succeeded
            end
        end        
    end
end

St2 = (Q2^SimulateSteps) * St2;
Error = sum(St2(Z:end));

end


