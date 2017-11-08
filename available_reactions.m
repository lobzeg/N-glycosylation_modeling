function [e, MM, Msum] = available_reactions(a,s,Vm,Km,ChangeCompartmentRate,concentration,NE,sia)

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

Km7 = Km(7);
Km8 = Km(8);
Km9 = Km(9);
KmGal = Km(10);
KmSia = Km(15);
KmE1 = Km(19);
KmE2 = Km(21);
T1 = 1 / 36;

% Concentration of enzymes depends on compartment
% T is a variable taking care of it
if compartment == 1
    for k = 1:(NE - 14)
        T(k) = 3 * T1 / 5;
    end;
    for k = (NE - 13) : (NE - 1)
        if k ~= 14
        T(k) = 0;
        end;
    end;
end;
T(14) = 3 * T1 / 5;

if compartment == 2
    for k = 1:(NE - 14)
        T(k) = 9 * T1 / 5;
    end;
    for k = (NE - 13) : (NE - 1)
        if k ~= 14
        T(k) = T1 / 5;
        end;
    end;
end;
T(14) = 9 * T1 / 5;

if compartment == 3
    for k = 1:(NE - 14)
        T(k) = 6 * T1 / 5;
    end;
    for k = (NE - 13) : (NE - 1)
        if k ~= 14
        T(k) = 4 * T1 / 5;
        end;
    end;
    
end;
T(14) = 6 * T1 / 5;

if compartment == 4
    for k = 1:(NE - 14)
        T(k) = 2 * T1 / 5;
    end;
    for k = (NE - 13) : (NE - 1)
        if k ~= 14
        T(k) = 15 * T1 / 5;
        end;
    end;
    
end;
T(14) = 2 * T1 / 5;

for j = 1:NE
    e(j) = 0;
end;
e(NE) = 1;

% Determine possible reactions based on glycan structure
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
        Km(7) = Km7 * 0.048;
    end;
end;

if (first_arm > 0 && bisecting_glcnac == 0 && second_arm == 0)
    e(8) = 1;
    kef = 1;
    if third_arm == 0
        T(8) = T(8) / 5;
    end;
    if fourth_arm > 1 && third_arm < 2
        kef = kef * 1.5;
    end;
    if fourth_arm > 0
        kef = kef * 0.178;
    end;
    Km(8) = Km8 * kef;
end;

if (third_arm == 1 && bisecting_glcnac == 0 && fourth_arm == 0)
    e(9) = 1;
    if second_arm > 0
        Km(9) = Km9 * 0.692;
    end;
end;
    
if (rem(first_arm,3) == 1)
    e(10) = 1;
    if bisecting_glcnac > 0 && third_arm > 0
        Km(10) = KmGal * 3.62;
    end;
end;

if (rem(second_arm,3) == 1)
    e(11) = 1;
    if bisecting_glcnac > 0 && third_arm > 0
        Km(11) = KmGal * 3.62;
    end;
end;

if (rem(third_arm,3) == 1)
    e(12) = 1;
    if bisecting_glcnac > 0 && third_arm > 0
        Km(12) = KmGal * 3.62;
    end;
end;

if (rem(fourth_arm,3) == 1)
    e(13) = 1;
    if bisecting_glcnac > 0 && third_arm > 0
        Km(13) = KmGal * 3.62;
    end;
end;
            
if ( fucose == 0 && bisecting_glcnac == 0 && first_arm == 1 && second_arm < 2 && third_arm < 2 && fourth_arm < 2 )
    e(14) = 1;
end;
    
if( rem(first_arm,3) == 2 )
    e(15) = 1;
    if sia == 1;
        Km(15) = KmSia * 5;
    end;
end;
    
if( rem(second_arm,3) == 2 )
    e(16) = 1;
    if sia == 1;
        Km(16) = KmSia * 5;
    end;
end;

if( rem(third_arm,3) == 2 )
    e(17) = 1;
    if sia == 1;
        Km(17) = KmSia * 5;
    end;
end;
    
if( rem(fourth_arm,3) == 2 )
    e(18) = 1;
    if sia == 1;
        Km(18) = KmSia * 5;
    end;
end;
    
if( rem(first_arm,3) == 2 )
    e(19) = 1;
    if second_arm == 0
        Km(19) = KmE1 * 4;
    end;
end;
    
if( rem(second_arm,3) == 2 )
    e(20) = 1;
    if first_arm == 0
        Km(20) = KmE1 * 4;
    end;
end;
    
if( rem(third_arm,3) == 2 )
    e(21) = 1;
    if fourth_arm == 0
        Km(21) = KmE2 * 4;
    end;
end;
    
if( rem(fourth_arm,3) == 2 )
    e(22) = 1;
    if third_arm == 0
        Km(22) = KmE2 * 4;
    end;
end;

% Determine likelihood of each reaction    
MM = 0;
for j = 1:NE - 1
    if (e(j) == 1)
        m(j) = Vm(j) * T(j) / (Km(j) + s(j,compartment) * concentration);                                        
    else
        m(j) = 0;
    end;
    MM = MM + m(j);
    Msum(j) = MM;
end;
m(NE) = ChangeCompartmentRate;
MM = MM + m(NE);
Msum(NE) = MM;

end