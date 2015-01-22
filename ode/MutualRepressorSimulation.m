%=========================================================================%
% Mutual repressor switch
%
% Deterministic Fractional Occupancy Model
% iGEM 2012 Team Slovenia - http://2012.igem.org/Team:Slovenia
%=========================================================================%
function A=MutualRepressorSimulation(n1,n2,n3,n4,leak_TALA_KRAB, leak_TALB_KRAB, signalPattern)
    TIME = 1200;            % Simulation duration.
    SamplesPerStep = 100;
    SAMPLES = TIME*SamplesPerStep;
    dt = 1/SamplesPerStep;  % Time step.

    % Graphs data.
    grafTALA_KRAB = zeros(1,SAMPLES);
    grafTALB_KRAB = zeros(1,SAMPLES);
    grafBFP = zeros(1,SAMPLES);
    grafMCitrin = zeros(1,SAMPLES);
    %grafSignal1 = zeros(1,SAMPLES);
    %grafSignal2 = zeros(1,SAMPLES);
    counter=1;     

    % Initial values/concentrations.
    BFP = 0;
    MCitrin = 0;
    TALA_KRAB = 0;
    TALB_KRAB = 0;
    PIPKRAB = 0;
    EKRAB = 0;

    Signal1 = 0;
    Signal2 = 0;

    % Fractional occupancy parameters - association constants (k1-k4) and exponents (n1-n4).
    k1 = 1;     
    k2 = 1;
    k3 = 1;
    k4 = 1;

    Ka = 3;     % Amount of activator required for 50% activation (minimal promoter).
    Kr = 1;     % Amount of repressor required for 50% repression (constitutive promoter).

    % Production & degradation parameters.
    k_BFP = 10;              % BFP production rate (active promoter).
    kb_BFP = 0.00;           % BFP basal rate (leaking - inactive promoter).
    deg_BFP = 0.1;           % BFP degradation rate.

    k_citrin = 10;              % MCitrin production rate (active pormoter).
    kb_citrin = 0.00;           % MCitrin basal rate (leaking - inactive promoter).
    deg_citrin = 0.1;           % MCitrin degradation rate.

    k1_TALA_KRAB = 10;         % TALA_KRAB production rate from construct 1.
    %kb1_TALA_KRAB = 0.00;      % TALA_KRAB basal rate from construct 1.
    kb1_TALA_KRAB = leak_TALA_KRAB * k1_TALA_KRAB;
    k3_TALA_KRAB = 10;         % TALA_KRAB production rate from construct 3.
    %kb3_TALA_KRAB = 0.00;      % TALA_KRAB basal rate from construct 3.
    kb3_TALA_KRAB = leak_TALA_KRAB * k3_TALA_KRAB;
    deg_TALA_KRAB = 0.1;       % TALA_KRAB degradation rate.

    k2_TALB_KRAB = 10;         % TALB_KRAB production rate from construct 2.
    kb2_TALB_KRAB = 0.00;      % TALB_KRAB basal rate from construct 2.
    kb2_TALB_KRAB = leak_TALB_KRAB * k2_TALB_KRAB;
    k4_TALB_KRAB = 10;         % TALB_KRAB production rate from construct 4.
    kb4_TALB_KRAB = 0.00;      % TALB_KRAB basal rate from construct 4.
    kb4_TALB_KRAB = leak_TALB_KRAB * k4_TALB_KRAB;
    deg_TALB_KRAB = 0.1;       % TALB_KRAB degradation rate.

    k_PIPKRAB = 10;             % PIPKRAB production rate.
    deg_PIPKRAB = 0.1;          % PIPKRAB degradation rate.

    k_EKRAB = 10;              % EKRAB production rate.
    deg_EKRAB = 0.1;           % EKRAB degradation rate.

    deg_PIPKRAB_fast = 500;     % PIPKRAB degradation rate when signal1 is present.
    deg_PIPKRAB_normal = 0.1;   % PIPKRAB degradation rate when signal1 is absent (normal degradation rate).

    deg_EKRAB_fast = 500;      % EKRAB degradation rate when signal2 is present.
    deg_EKRAB_normal = 0.1;    % EKRAB degradation rate when signal2 is absent (normal degradation rate).


    %========== SIMULATION ==========%
    for t=0:dt:(TIME-1)                   
    % State 1.
    if t==0
        Signal1 = signalPattern(1,1);
        Signal2 = signalPattern(2,1);
    end
    if t==300
        Signal1 = signalPattern(1,2);
        Signal2 = signalPattern(2,2);
    end
    if t==600
        Signal1 = signalPattern(1,3);
        Signal2 = signalPattern(2,3);
    end



    % Calculate fractional occupancies
    f1 = Kr/(Kr+k1*TALB_KRAB^n1);            % Construct 1 (BFP).
    f2 = Kr/(Kr+k2*TALA_KRAB^n2);            % Construct 2 (citrin).    
    f3 = Kr/(Kr+k3*PIPKRAB^n3);               % Construct 3.    
    f4 = Kr/(Kr+k4*EKRAB^n4);                 % Construct 4.    

    % Signal1 present -> PIPKRAB cannot bind (cannot repress construct 3).
    % This is modeled as fast PIPKRAB degradation in the presence of signal1.
    if Signal1==1
    deg_PIPKRAB = deg_PIPKRAB_fast;
    else
    deg_PIPKRAB = deg_PIPKRAB_normal;
    end

    % Signal2 present -> EKRAB cannot bind (cannot repress construct 4).
    if Signal2==1
    deg_EKRAB = deg_EKRAB_fast;
    else
    deg_EKRAB = deg_EKRAB_normal;
    end


    % ODEs.
    % Form: dProtein/dt = production_rate * P(active promoter) + basal_rate*P(inactive promoter) - degradation_rate*Protein
    dBFP = dt * (k_BFP*f1 + kb_BFP*(1-f1) - deg_BFP*BFP);
    dMCitrin = dt * (k_citrin*f2 + kb_citrin*(1-f2) - deg_citrin*MCitrin);
    dTALA_KRAB = dt * (k1_TALA_KRAB*f1 + kb1_TALA_KRAB*(1-f1) + k3_TALA_KRAB*f3 + kb3_TALA_KRAB*(1-f3) - deg_TALA_KRAB*TALA_KRAB);
    dTALB_KRAB = dt * (k2_TALB_KRAB*f2 + kb2_TALB_KRAB*(1-f2) + k4_TALB_KRAB*f4 + kb4_TALB_KRAB*(1-f4) - deg_TALB_KRAB*TALB_KRAB);
    dPIPKRAB = dt * (k_PIPKRAB - deg_PIPKRAB*PIPKRAB);
    dEKRAB = dt * (k_EKRAB - deg_EKRAB*EKRAB);

    % Update values.
    BFP = BFP + dBFP;
    MCitrin = MCitrin + dMCitrin;
    TALA_KRAB = TALA_KRAB + dTALA_KRAB;
    TALB_KRAB = TALB_KRAB + dTALB_KRAB;
    PIPKRAB = PIPKRAB + dPIPKRAB;
    EKRAB = EKRAB + dEKRAB;

    if (PIPKRAB < 0) PIPKRAB = 0; end
    if (EKRAB < 0) EKRAB = 0; end


    % Graphs.
    grafBFP(counter) = BFP;
    grafMCitrin(counter) = MCitrin;
    %grafSignal1(counter) = Signal1
    %grafSignal2(counter) = Signal2
    counter = counter+1;    
    end

    time = (1:SAMPLES)/SamplesPerStep;
    A = [time', grafBFP', grafMCitrin'];
endfunction





