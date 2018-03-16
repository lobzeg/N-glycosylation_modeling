function [e, m_m, m_sum] = available_reactions(a,s,v_m,k_m,change_compartment_rate,concentration,n_e,sia)

% e(1)    ManI
% e(2)    ManI
% e(3)    ManI
% e(4)    ManII
% e(5)    GnTI
% e(6)    GnTII
% e(7)    GnTIII
% e(8)    GnTIV
% e(9)    GnTV
% e(10)   GalT1
% e(11)   GalT2
% e(12)   GalT3
% e(13)   GalT4
% e(14)   FucT
% e(15)   SiaT1
% e(16)   SiaT2
% e(17)   SiaT3
% e(18)   SiaT4
% e(19)   GnTE1
% e(20)   GnTE2
% e(21)   GnTE3
% e(22)   GnTE4
% e(23)   compartment change

glucose = a(1);
manose_first_branch = a(2);
manose_second_branch = a(3);
manose_third_branch = a(4);
bisecting_glcnac = a(5);
fucose = a(6);
first_arm = a(7);
second_arm = a(8);
third_arm = a(9);
fourth_arm = a(10);
compartment = a(11);

k_m7 = k_m(7);
k_m8 = k_m(8);
k_m9 = k_m(9);
k_m_gal = k_m(10);
k_m_sia = k_m(15);
k_m_e1 = k_m(19);
k_m_e2 = k_m(21);
t1 = 1 / 36;

% Concentration of enzymes depends on compartment
% t is a variable taking care of it
if compartment == 1
    for k = 1:(n_e - 14)
        t(k) = 3 * t1 / 5;
    end;
    for k = (n_e - 13) : (n_e - 1)
        if k ~= 14
        t(k) = 0;
        end;
    end;
end;
t(14) = 3 * t1 / 5;

if compartment == 2
    for k = 1:(n_e - 14)
        t(k) = 9 * t1 / 5;
    end;
    for k = (n_e - 13) : (n_e - 1)
        if k ~= 14
        t(k) = t1 / 5;
        end;
    end;
end;
t(14) = 9 * t1 / 5;

if compartment == 3
    for k = 1:(n_e - 14)
        t(k) = 6 * t1 / 5;
    end;
    for k = (n_e - 13) : (n_e - 1)
        if k ~= 14
        t(k) = 4 * t1 / 5;
        end;
    end;
    
end;
t(14) = 6 * t1 / 5;

if compartment == 4
    for k = 1:(n_e - 14)
        t(k) = 2 * t1 / 5;
    end;
    for k = (n_e - 13) : (n_e - 1)
        if k ~= 14
        t(k) = 15 * t1 / 5;
        end;
    end;
    
end;
t(14) = 2 * t1 / 5;

for j = 1:n_e
    e(j) = 0;
end;
e(n_e) = 1;

% Determine possible reactions based on the glycan structure
if (manose_first_branch > 0)
    e(1) = 1;
end;
if (manose_second_branch > 1)
    e(2) = 1;
end;
if (manose_third_branch > 1)
    e(3) = 1;
end;

if (manose_second_branch == 1 && manose_third_branch == 1 && first_arm == 1 && bisecting_glcnac == 0)
    e(4) = 1;
end;

if (manose_first_branch == 0 && manose_second_branch == 1 && manose_third_branch == 1 && first_arm == 0)
    e(5) = 1;
end;

if (manose_second_branch == 0 && manose_third_branch == 0 && first_arm == 1 && bisecting_glcnac == 0 && third_arm == 0)
    e(6) = 1;
end;

if (first_arm > 0 && bisecting_glcnac == 0 && first_arm < 2 && second_arm < 2 && third_arm < 2 && fourth_arm < 2 )
    e(7) = 1;
    if third_arm > 0
        k_m(7) = k_m7 * 0.048;
    end;
end;

if (first_arm > 0 && bisecting_glcnac == 0 && second_arm == 0)
    e(8) = 1;
    kef = 1;
    if third_arm == 0
        t(8) = t(8) / 5;
    end;
    if fourth_arm > 1 && third_arm < 2
        kef = kef * 1.5;
    end;
    if fourth_arm > 0
        kef = kef * 0.178;
    end;
    k_m(8) = k_m8 * kef;
end;

if (third_arm == 1 && bisecting_glcnac == 0 && fourth_arm == 0)
    e(9) = 1;
    if second_arm > 0
        k_m(9) = k_m9 * 0.692;
    end;
end;
    
if (rem(first_arm,3) == 1)
    e(10) = 1;
    if bisecting_glcnac > 0 && third_arm > 0
        k_m(10) = k_m_gal * 3.62;
    end;
end;

if (rem(second_arm,3) == 1)
    e(11) = 1;
    if bisecting_glcnac > 0 && third_arm > 0
        k_m(11) = k_m_gal * 3.62;
    end;
end;

if (rem(third_arm,3) == 1)
    e(12) = 1;
    if bisecting_glcnac > 0 && third_arm > 0
        k_m(12) = k_m_gal * 3.62;
    end;
end;

if (rem(fourth_arm,3) == 1)
    e(13) = 1;
    if bisecting_glcnac > 0 && third_arm > 0
        k_m(13) = k_m_gal * 3.62;
    end;
end;
            
if ( fucose == 0 && bisecting_glcnac == 0 && first_arm == 1 && second_arm < 2 && third_arm < 2 && fourth_arm < 2 )
    e(14) = 1;
end;
    
if( rem(first_arm,3) == 2 )
    e(15) = 1;
    if sia == 1;
        k_m(15) = k_m_sia * 5;
    end;
end;
    
if( rem(second_arm,3) == 2 )
    e(16) = 1;
    if sia == 1;
        k_m(16) = k_m_sia * 5;
    end;
end;

if( rem(third_arm,3) == 2 )
    e(17) = 1;
    if sia == 1;
        k_m(17) = k_m_sia * 5;
    end;
end;
    
if( rem(fourth_arm,3) == 2 )
    e(18) = 1;
    if sia == 1;
        k_m(18) = k_m_sia * 5;
    end;
end;
    
if( rem(first_arm,3) == 2 )
    e(19) = 1;
    if second_arm == 0
        k_m(19) = k_m_e1 * 4;
    end;
end;
    
if( rem(second_arm,3) == 2 )
    e(20) = 1;
    if first_arm == 0
        k_m(20) = k_m_e1 * 4;
    end;
end;
    
if( rem(third_arm,3) == 2 )
    e(21) = 1;
    if fourth_arm == 0
        k_m(21) = k_m_e2 * 4;
    end;
end;
    
if( rem(fourth_arm,3) == 2 )
    e(22) = 1;
    if third_arm == 0
        k_m(22) = k_m_e2 * 4;
    end;
end;

% Determine likelihood for each reaction    
m_m = 0;
for j = 1:n_e - 1
    if (e(j) == 1)
        m(j) = v_m(j) * t(j) / (k_m(j) + s(j,compartment) * concentration);                                        
    else
        m(j) = 0;
    end;
    m_m = m_m + m(j);
    m_sum(j) = m_m;
end;
m(n_e) = change_compartment_rate;
m_m = m_m + m(n_e);
m_sum(n_e) = m_m;

end
