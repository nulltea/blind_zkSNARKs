//// Timings for HE scheme in sec

get_timings := function(curqbits)
    
    if n eq 2^13 then
        print "WARNING: don't have all timings for full dynamic scaling";
        if curqbits le 120 then
            Tadd := 0.00001; 
            Tptct := 0.00002;
            Taut := 0.0004;
            Tctct := 0.00344;
        elif curqbits le 240 then
            Tadd := 0.00004; 
            Tptct := 0.0008;
            Taut := 0.00132;
            Tctct := 0.00728;
        else
            Tadd := 0.00007; 
            Tptct := 0.0018;
            Taut := 0.00354;
            Tctct := 0.0145;
        end if;
    elif n eq 2^14 or n eq 12288 then 
        if curqbits le 120 then
            if curqbits le 60 then print "WARNING: measure smaller mod timings"; end if; // not possible in current GBFV implementation
            Tadd := 0.00002; 
            Tptct := 0.00006;
            Taut := 0.00084;
            Tctct := 0.00558;
        elif curqbits eq 180 then
            Tadd := 0.00004; 
            Tptct := 0.0001;
            Taut := 0.00166;
            Tctct := 0.00944;
        elif curqbits eq 240 then
            Tadd := 0.00004; 
            Tptct := 0.0001;
            Taut := 0.00166;
            Tctct := 0.00946;
        elif curqbits eq 300 then
            Tadd := 0.0001; 
            Tptct := 0.00024;
            Taut := 0.00398;
            Tctct := 0.01812;
        elif curqbits eq 360 then
            Tadd := 0.00012; 
            Tptct := 0.00028;
            Taut := 0.00552;
            Tctct := 0.02306;
        else // curqbits eq 420 bits
            Tadd := 0.00014; 
            Tptct := 0.00034;
            Taut := 0.00724;
            Tctct := 0.0284;
        end if;
    end if;

    return Tadd, Tptct, Taut, Tctct;
end function;

if n eq 2^13 then
    Tenc := 0.007; // TODO: not yet updated
elif n eq 2^14 or n eq 12288 then
    Tenc := 0.013;  // TODO: not yet updated
end if;

//// Noise growth for HE scheme

getHEnoise := function(t)
    Nadd := 0.5; // uses avg case
    Naut := 0; // most of the time does not increase noise
    Nms := 4.6;

    i := Log(2,Degree(t)/7);
    assert m eq (2^11)*3*7;
    if i eq 5 then
        return 8.5, Nadd, Naut, 6.3, 10.2, Nms;
    elif i eq 6 then
        return 8.5, Nadd, Naut, 7.3, 11.1, Nms;
    elif i eq 7 then
        return 8.6, Nadd, Naut, 9.1, 13.0, Nms;
    elif i eq 8 then
        return 8.8, Nadd, Naut, 12.9, 16.9, Nms;
    elif i eq 9 then
        return 8.7, Nadd, Naut, 20.9, 24.7, Nms;
    elif i eq 10 then
        return 8.4, Nadd, Naut, 37, 41, Nms;
    else
        error "No noise estimates for t(x)";
    end if;
end function;

