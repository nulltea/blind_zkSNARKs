//// Some helper functions

// for storing previous cumulative values
pCadd := 0;
pCptct := 0;
pCaut := 0;
pCctct := 0;
pN := N0;

modswitch := procedure(~N, MSeffect)
    if MSeffect gt 0 then tmp := 1; else tmp := 0; end if;
    N := Log(2, 2^(N-MSeffect) + tmp*(2^Nms));
end procedure;

cumulate := procedure(Cadd, Cptct, Caut, Cctct, ~pN, ~pCadd, ~pCptct, ~pCaut, ~pCctct,\
    text, ~T, ~qbits, ~MSeffect, ~Nmax)

    assert Nmax lt qbits;
    
    Tadd, Tptct, Taut, Tctct := get_timings(qbits);

    if print_info then print "INFO: estimating timing with modulus", qbits, "bits"; end if;

    Nnew := Nmax - pN;
    Tnew := (Cadd-pCadd)*Tadd + (Cptct-pCptct)*Tptct + (Caut-pCaut)*Taut + (Cctct-pCctct)*Tctct;
    if print_nice then
        print "Timing:", RealField(5) ! Tnew/3600, "(h), Noise:",Round(Nnew),"for",text;
    else
        print text,"&",Round(Nnew),"&",Ceiling(Cadd-pCadd),"&",Ceiling(Cptct-pCptct),"&",\
            Ceiling(Caut-pCaut),"&",Cctct-pCctct,"\\\\";
        print "\\hline";
    end if;
    
    // checking for possible modswitch
    MSeffect := qbits - Ceiling((qbits-Nmax)/RNSmodbits)*RNSmodbits;
    MSeffect := (MSeffect + Abs(MSeffect))/2; // make zero if negative
    
    // prevent from modswitching to 60 bits in zkdel (since switches to 48bits in PoD)
    if (not vcoed) and qbits - MSeffect eq 60 then
        MSeffect := 0;
    end if;

    qbits -:= MSeffect;
    assert MSeffect ge 0;
    
    // possibly perform the modswitch
    modswitch(~Nmax, MSeffect);
    // don't forget to also modswitch other ciphertexts with less noise

    pN := Nmax;
    T +:= Tnew;
    pCadd := Cadd;
    pCptct := Cptct;
    pCaut := Caut;
    pCctct := Cctct;

end procedure;

// Returns the expected number of partitions hit when sampling k distinct
// elements out of a set of n that is divided into m equal partitions
E_num_hits := function(n, m, k)
    n := Integers() ! n;
    m := Integers() ! m;
    assert n mod m eq 0;
    p := n/m;
    if k eq n then
        return m;
    else
        return m*(1-&*[(n - p - i)/(n - i) : i in [0..k-1]]);
    end if;
end function;

get_estimates := function(len, avghw)
    len := Integers() ! len;
    assert len mod pack eq 0;
    l := len/pack;
    // expected num of nonzero plaintext diagonals per subrow
    t1 := E_num_hits(len*pack, len, avghw*pack);
    // expected num of nonzero diagonal "types" per subcolumn
    t2 := E_num_hits((pack-1)*len, pack-1, avghw*(pack-1)); // -1 because one does not require aut
    
    logpack := Log(2,pack);
    assert 2^Ceiling(logpack) eq pack;
    a := 2^Ceiling(logpack/2);
    b := 2^Floor(logpack/2);
    assert a*b eq pack;
    // expected num of nonzero diagonals in inner part of BSGS
    t3 := E_num_hits(pack*b*(a-1)*l, a-1, avghw*b*(a-1));
    t3 := a-1; // non-sparse BSGS
    // expected num of nonzero diagonals in outer part of BSGS
    t4 := E_num_hits(pack*(b-1)*a*l, l*(b-1), avghw*(b-1)*a);
    t4 := l*(b-1); // non-sparse BSGS
    // SPARSENESS DOES NOT SEEM TO MAKE A DIFFERENCE
    
    return t1, t2, t3, t4;
end function;

// Time/Noise for Matrix-Vector multiply
C_MV := procedure(~Cadd, ~Cptct, ~Caut, len : avghw := len, rep := 1, halfzero := false) 
    t1, t2, t3, t4 := get_estimates(len, avghw);
    l := len/pack;

    if halfzero then hz := 2; else hz := 1; end if;
    
    for i in [1..rep] do
        // Halevi-Shoup method in submatrices
        // can consider each of the l subrows independently
        Cptct +:= l*t1/hz;
        Cadd +:= l*(t1/hz-1);
        Caut +:= (l/hz)*Min(t2, t3+t4);
    end for;
    if t2 gt t3+t4 and print_info then
        print "INFO: using BSGS";
    end if;
end procedure;

N_MV := function(len : avghw := len, halfzero := false) 
    t1, t2, t3, t4 := get_estimates(len, avghw);
    l := len/pack;

    if halfzero then hz := 2; else hz := 1; end if;
    
    // Halevi-Shoup method in submatrices
    if t2 lt t3+t4 then
        return Nptct + Naut + Ceiling(Log(2,t1/hz))*Nadd;
    else
        return Nptct + 2*Naut + Ceiling(Log(2,(t3+1)*(t4+l)/hz))*Nadd;
    end if;
end function;

// Time/Noise for (inverse) NTT
C_NTT := procedure(~Cadd, ~Cptct, ~Caut, len, base : w := 1, halfzero := false)
    assert len mod pack eq 0;
    columns := pack*w;
    rows := len/columns;
    assert base^Ceiling(Log(base,rows)) eq rows;
 
    // row-wise NTT
    temp := (base-1)*Log(base,rows)*rows*w;
    Cadd +:= temp;
    Cptct +:= temp;
    
    // multiplication by twiddles can be squashed in
    
    // column-wise NTT
    if halfzero and w gt 1 then
        // this assumes that the MV-product NTT is performed first
        C_MV(~Cadd, ~Cptct, ~Caut, columns : rep := rows, halfzero := true);
    else
        C_MV(~Cadd, ~Cptct, ~Caut, columns : rep := rows);
    end if;

end procedure;

N_NTT := function(len, base : w := 1)
    assert len mod pack eq 0;
    columns := pack*w;
    rows := len/columns;
    assert base^Ceiling(Log(base,rows)) eq rows;
    
    return N_MV(columns) + Log(base, rows)*((base-1)*Nadd + Nptct);
end function;
