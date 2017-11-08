function [a,sia] = perform_reaction(MM,Msum,a,sia,NE)

Y = random('unif',0,MM);
for j = 1:NE
    if Y < Msum(j)
        KK = j;
        Y = MM + 1;
    end;
end;

switch KK
    case 1
        a(2) = a(2) - 1;
    case 2
        a(3) = a(3) - 1;
    case 3
        a(4) = a(4) - 1;
    case 4
        a(3) = a(3) - 1;
        a(4) = a(4) - 1;
    case 5
        a(7) = 1;
    case 6
        a(9) = 1;
    case 7
        a(5) = 1;
    case 8
        a(8) = 1;
    case 9
        a(10) = 1;
    case 10
        a(7) = a(7) + 1;
    case 11
        a(8) = a(8) + 1;
    case 12
        a(9) = a(9) + 1;
    case 13
        a(10) = a(10) + 1;
    case 14
        a(6) = 1;
    case 15;
        a(7) = a(7) + 1;
        sia = 1;
    case 16;
        a(8) = a(8) + 1;
        sia = 1;
    case 17;
        a(9) = a(9) + 1;
        sia = 1;
    case 18;
        a(10) = a(10) + 1;
        sia = 1;
    case 19;
        a(7) = a(7) + 2;
    case 20;
        a(8) = a(8) + 2;
    case 21;
        a(9) = a(9) + 2;
    case 22;
        a(10) = a(10) + 2;
    case 23
        a(11) = a(11) + 1;
end;
end