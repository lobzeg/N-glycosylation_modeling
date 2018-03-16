function [sx,ssx] = adjustment(a,sx,ssx,v_m,k_m,change_compartment_rate,concentration,div,n_e)

compartment = a(11);
sia = 0;

% Normalize VPTC
for i = 1:4
    for j = 1:n_e
        s(j,i) = sx(j,i) / ssx(i);
    end;
end;

while compartment < 5
    
    % Determine possible reactions and likelihood of each of them
    [e, m_m, m_sum] = available_reactions(a,s,v_m,k_m,change_compartment_rate,concentration,n_e,sia);

    % Time until the next reaction
    if m_m > 0
        time = -log(rand()) / m_m;
    end;
    
    % Take into account time until the next reaction for all possible reactions (i.e. time glycan spent being available for reaction X)
    if e(1) == 1
        sx(1,compartment) = sx(1,compartment) + time;
        sx(2,compartment) = sx(2,compartment) + time;
        sx(3,compartment) = sx(3,compartment) + time;
    end;
    
    if e(2) == 1
        sx(1,compartment) = sx(1,compartment) + time;
        sx(2,compartment) = sx(2,compartment) + time;
        sx(3,compartment) = sx(3,compartment) + time;
    end;
    
    if e(3) == 1
        sx(1,compartment) = sx(1,compartment) + time;
        sx(2,compartment) = sx(2,compartment) + time;
        sx(3,compartment) = sx(3,compartment) + time;
    end;
    
    for j = 4:9
        if e(j) == 1
            sx(j,compartment) = sx(j,compartment) + time;
        end;
    end;
    
    if e(14) == 1
        sx(14,compartment) = sx(14,compartment) + time;
    end;            
    
    if e(10) == 1
        sx(10,compartment) = sx(10,compartment) + time;
        sx(11,compartment) = sx(11,compartment) + time;
        sx(12,compartment) = sx(12,compartment) + time;
        sx(13,compartment) = sx(13,compartment) + time;
    end;                        
    if e(11) == 1
        sx(10,compartment) = sx(10,compartment) + time;
        sx(11,compartment) = sx(11,compartment) + time;
        sx(12,compartment) = sx(12,compartment) + time;
        sx(13,compartment) = sx(13,compartment) + time;
    end;                        
    if e(12) == 1
        sx(10,compartment) = sx(10,compartment) + time;
        sx(11,compartment) = sx(11,compartment) + time;
        sx(12,compartment) = sx(12,compartment) + time;
        sx(13,compartment) = sx(13,compartment) + time;
    end;                        
    if e(13) == 1
        sx(10,compartment) = sx(10,compartment) + time;
        sx(11,compartment) = sx(11,compartment) + time;
        sx(12,compartment) = sx(12,compartment) + time;
        sx(13,compartment) = sx(13,compartment) + time;
    end;
    
    if e(15) == 1
        sx(15,compartment) = sx(15,compartment) + time;
        sx(16,compartment) = sx(16,compartment) + time;
        sx(17,compartment) = sx(17,compartment) + time;
        sx(18,compartment) = sx(18,compartment) + time;
    end;
    if e(16) == 1
        sx(15,compartment) = sx(15,compartment) + time;
        sx(16,compartment) = sx(16,compartment) + time;
        sx(17,compartment) = sx(17,compartment) + time;
        sx(18,compartment) = sx(18,compartment) + time;
    end;
    if e(17) == 1
        sx(15,compartment) = sx(15,compartment) + time;
        sx(16,compartment) = sx(16,compartment) + time;
        sx(17,compartment) = sx(17,compartment) + time;
        sx(18,compartment) = sx(18,compartment) + time;
    end;
    if e(18) == 1
        sx(15,compartment) = sx(15,compartment) + time;
        sx(16,compartment) = sx(16,compartment) + time;
        sx(17,compartment) = sx(17,compartment) + time;
        sx(18,compartment) = sx(18,compartment) + time;
    end;
        
    if e(19) == 1
        sx(19,compartment) = sx(19,compartment) + time;
        sx(20,compartment) = sx(20,compartment) + time;
        sx(21,compartment) = sx(21,compartment) + time;
        sx(22,compartment) = sx(22,compartment) + time;
    end;
    if e(20) == 1
        sx(19,compartment) = sx(19,compartment) + time;
        sx(20,compartment) = sx(20,compartment) + time;
        sx(21,compartment) = sx(21,compartment) + time;
        sx(22,compartment) = sx(22,compartment) + time;
    end;
    if e(21) == 1
        sx(19,compartment) = sx(19,compartment) + time;
        sx(20,compartment) = sx(20,compartment) + time;
        sx(21,compartment) = sx(21,compartment) + time;
        sx(22,compartment) = sx(22,compartment) + time;
    end;
    if e(22) == 1
        sx(19,compartment) = sx(19,compartment) + time;
        sx(20,compartment) = sx(20,compartment) + time;
        sx(21,compartment) = sx(21,compartment) + time;
        sx(22,compartment) = sx(22,compartment) + time;
    end;
    
    % Take into account time until the next reaction as time spent in the system (i.e. time glycan spent in Golgi)
    ssx(compartment) = ssx(compartment) + time;
    
    % Perform a random reaction (from available ones)
    [a,sia] = perform_reaction(m_m,m_sum,a,sia,n_e);
    compartment = a(11);
end;

% Decrease weight of previous adjustment runs
for i = 1:4                    
    for j = 1:n_e
        sx(j,i) = sx(j,i) * div;                        
    end;
    ssx(i) = ssx(i) * div;
end;
end
