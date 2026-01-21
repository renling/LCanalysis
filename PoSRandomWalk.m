alpha = 0.8;
D = 2 * 1/20;       % Cardano 20s per block

Alphabet = 8;
States = 24;
KK = 20;
% Alphabet is max possible epoch length
%   need to be large enough to ensure numeric precision of P(j, 2) 
%   as well as neligible probability of larger j
% States is the number of states in the Markov chain tracked
% KK is the max number of confirmation to evaluate

[Pa, PH, PD, PA, PAD] = PoSSlotPdf(alpha, D, Alphabet);

% PoS warmup and final stages are the same as PoW
St0 = PoWMCWarmupUB(PAD, Alphabet, States); 
ErrorUB = zeros(KK, 1);
tic
for K = 1:KK
%    St2 = PoSMCConfirmUB(K, Pa, PH, PD, PA, St0, Alphabet, States);
%    ErrorUB(K) = PoWMCFinalUB(PAD, St2, Alphabet, States);
end
toc

St0 = PoWMCWarmupLB(PAD, Alphabet, States);
ErrorLB1 = zeros(KK, 1);
ErrorLB2 = zeros(KK, 1);
tic
for K = 1:KK
    % private mining as lower bound     
    %St2 = PoSMCConfirmPM(K, Pa, PH, PD, PA, PAD, St0, Alphabet, States);
    %ErrorLB1(K) = PoWMCFinalLB(PAD, St2, Alphabet, States);    
    ErrorLB2(K) = PoSMCConfirmLB(K, Pa, PH, PD, PA, PAD, St0, Alphabet, States);
end
toc
ErrorLB = ErrorLB2;

%csvwrite('graph-ADA-2s-0.2-upper.csv', [(1:KK)',ErrorUB])
csvwrite('graph-ADA-2s-0.2-lower.csv', [(1:KK)',ErrorLB])

