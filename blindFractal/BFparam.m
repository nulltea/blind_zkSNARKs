//// Set parameters of blind fractal

lpack := Log(2,pack);

if lpack eq 6 then
    b1 := 2^2; // len_z
    w1 := 2^0;
    b2 := 2^2; // len_L
    w2 := 2^1;
elif lpack eq 7 then
    b1 := 2^2; // len_z
    w1 := 2^1;
    b2 := 2^2;  // len_L
    w2 := 2^2;
elif lpack eq 8 then
    b1 := 2^3; // len_z
    w1 := 2^0;
    b2 := 2^3; // len_L
    w2 := 2^1;
elif lpack eq 9 then
    b1 := 2^2; // len_z
    w1 := 2^3;
    b2 := 2^4; // len_L
    w2 := 2^4;
end if;

b3 := b2; // 2*len_z
w3 := w2;
b4 := b2; // len_L
w4 := w2;

