/*
  Exact computational certificates for
  "Quadratic Cremona transformations of P^4".

  All computations are over Q. They supplement, but do not replace, the
  geometric arguments in the paper. Load this file in Magma and call

      RunChecks();

  The procedures print a short confirmation and stop with an assertion if a
  certificate fails.
*/

QQ := Rationals();
print "Loading Cremona_P4_checks.m -- build 2026-07-18";

function NormalizePolynomial(f)
    if f eq 0 then
        return f;
    end if;
    return f / LeadingCoefficient(f);
end function;

function DeleteRowMatrix(M, k)
    R := BaseRing(M);
    entries := [];
    for i in [1..Nrows(M)] do
        if i ne k then
            for j in [1..Ncols(M)] do
                Append(~entries, M[i,j]);
            end for;
        end if;
    end for;
    return Matrix(R, Nrows(M)-1, Ncols(M), entries);
end function;

function MaximalMinorsFiveByFour(M)
    assert Nrows(M) eq 5 and Ncols(M) eq 4;
    return [ Determinant(DeleteRowMatrix(M,i)) : i in [1..5] ];
end function;

function SameRowSpace(A, B)
    if Ncols(A) ne Ncols(B) or Rank(A) ne Rank(B) then
        return false;
    end if;
    R := BaseRing(A);
    entries := [];
    for i in [1..Nrows(A)] do
        for j in [1..Ncols(A)] do
            Append(~entries,A[i,j]);
        end for;
    end for;
    for i in [1..Nrows(B)] do
        for j in [1..Ncols(B)] do
            Append(~entries,B[i,j]);
        end for;
    end for;
    C := Matrix(R,Nrows(A)+Nrows(B),Ncols(A),entries);
    return Rank(C) eq Rank(A);
end function;

function BinaryQuarticInvariants(f, s, r)
    a := MonomialCoefficient(f, s^4);
    b := MonomialCoefficient(f, s^3*r);
    c := MonomialCoefficient(f, s^2*r^2);
    d := MonomialCoefficient(f, s*r^3);
    e := MonomialCoefficient(f, r^4);
    I := 12*a*e - 3*b*d + c^2;
    J := 72*a*c*e + 9*b*c*d - 27*a*d^2 - 27*b^2*e - 2*c^3;
    return I, J;
end function;

/* ----------------------------------------------------------------------- */
/* Plane boundary: Smith form, residual quartics, and cross-ratio.         */
/* ----------------------------------------------------------------------- */

function QuadricPairs()
    pairs := [];
    for i in [1..5] do
        for j in [i..5] do
            Append(~pairs, <i,j>);
        end for;
    end for;
    return pairs;
end function;

function EvaluationRow(p)
    return [ p[pair[1]]*p[pair[2]] : pair in QuadricPairs() ];
end function;

function PolarRow(p, d)
    row := [];
    for pair in QuadricPairs() do
        i := pair[1];
        j := pair[2];
        if i eq j then
            Append(~row, 2*p[i]*d[i]);
        else
            Append(~row, p[i]*d[j] + p[j]*d[i]);
        end if;
    end for;
    return row;
end function;

function PlaneConditionMatrix(c)
    Rt<t> := PolynomialRing(QQ);
    e0 := Vector(Rt,[1,0,0,0,0]);
    e1 := Vector(Rt,[0,1,0,0,0]);
    e2 := Vector(Rt,[0,0,1,0,0]);
    e3 := Vector(Rt,[0,0,0,1,0]);
    e4 := Vector(Rt,[0,0,0,0,1]);

    a := e0 + t*(-2*e0 + e1 + e2 + c*e3 - e4);
    cc := e2 + t*(e0 + 2*e2 + e3 - 3*e4);
    b := e1 + t*(e0 - 3*e1 + 3*e2 - e4);
    d := e2 + t*(e0 - 2*e1 - 2*e2 + 2*e3);
    p := e0 + e1 + e2 + t*(e0 + 3*e1 + e2);
    q := e0 + e1 + e3 + 2*e4
         + t*(-3*e0 + 2*e1 + 3*e2 - 3*e3 - 2*e4);
    delta1 := e0 - e2 + t*(2*e0 + 3*e1 - 2*e2 - 2*e3 + 2*e4);
    delta2 := e1 - e2 + t*(-2*e0 + 3*e1 + e2 + 2*e4);

    rows := [
        EvaluationRow(a), PolarRow(a,cc), EvaluationRow(cc),
        EvaluationRow(b), PolarRow(b,d), EvaluationRow(d),
        EvaluationRow(p), PolarRow(p,delta1), PolarRow(p,delta2),
        EvaluationRow(q)
    ];
    return Matrix(Rt,rows);
end function;

print "PlaneConditionMatrix loaded";

function ExpectedPlaneLimit(c)
    if c eq -3 then
        return Matrix(QQ,[
          [0,0,0,44,0, 0,0,0,24, 0,-84,28, 0,0,-23],
          [0,0,0,0,22, 0,0,0,40, 0,-30,10, 0,0,-31],
          [0,0,0,0,0, 0,0,44,-24, 0,-48,16, 0,0,1],
          [0,0,0,0,0, 0,0,0,0, 0,0,0, 4,0,-1],
          [0,0,0,0,0, 0,0,0,0, 0,0,0, 0,2,-1]
        ]);
    end if;
    assert c eq -1;
    return Matrix(QQ,[
      [0,0,0,68,0, 0,0,0,-40, 0,12,-4, 0,0,3],
      [0,0,0,0,34, 0,0,0,100, 0,-30,10, 0,0,-67],
      [0,0,0,0,0, 0,0,68,24, 0,-48,16, 0,0,-29],
      [0,0,0,0,0, 0,0,0,0, 0,0,0, 4,0,-1],
      [0,0,0,0,0, 0,0,0,0, 0,0,0, 0,2,-1]
    ]);
end function;

function ResidualMatrix(c)
    Rb<r,s> := PolynomialRing(QQ,2);
    if c eq -3 then
        M := Matrix(Rb,[
          [44*s,24*r,-84*s+28*r,-23*r^2],
          [22*r,40*r,-30*s+10*r,-31*r^2],
          [0,44*s-24*r,-48*s+16*r,r^2],
          [0,0,0,4*s^2-r^2],
          [0,0,0,2*s*r-r^2]
        ]);
        g := (r-3*s)*(r-2*s)*(2*r-5*s)*(3*r-s);
    else
        assert c eq -1;
        M := Matrix(Rb,[
          [68*s,-40*r,12*s-4*r,3*r^2],
          [34*r,100*r,-30*s+10*r,-67*r^2],
          [0,68*s+24*r,-48*s+16*r,-29*r^2],
          [0,0,0,4*s^2-r^2],
          [0,0,0,2*s*r-r^2]
        ]);
        g := (r-3*s)*(r-2*s)*(r+5*s)*(2*r-s);
    end if;
    return M, g, r, s;
end function;

procedure PlaneBoundaryCertificate()
    for c in [-3,-1] do
        M := PlaneConditionMatrix(c);
        S, U, V := SmithForm(M);
        assert Rank(M) eq 10;
        assert &and[ S[i,i] ne 0 : i in [1..10] ];

        Kentries := [];
        for i in [1..15] do
            for j in [11..15] do
                Append(~Kentries, V[i,j]);
            end for;
        end for;
        K := Matrix(BaseRing(M),15,5,Kentries);

        K0entries := [];
        for i in [1..15] do
            for j in [1..5] do
                Append(~K0entries, Evaluate(K[i,j],0));
            end for;
        end for;
        K0 := Matrix(QQ,15,5,K0entries);

        E := ExpectedPlaneLimit(c);
        assert Rank(K0) eq 5;
        assert SameRowSpace(Transpose(K0),E);

        R, g, r, s := ResidualMatrix(c);
        minors := MaximalMinorsFiveByFour(R);
        nonzero := [ f : f in minors | f ne 0 ];
        assert #nonzero gt 0;
        gcdmin := GCD(nonzero);
        assert NormalizePolynomial(gcdmin) eq NormalizePolynomial(g);
        quotients := [ f div gcdmin : f in nonzero ];
        assert GCD([g] cat quotients) eq 1;

        I, J := BinaryQuarticInvariants(g,s,r);
        if c eq -3 then
            assert I eq 129 and J eq -1134;
            assert J^2/I^3 eq QQ!47628/79507;
        else
            assert I eq 961 and J eq -35282;
            assert J^2/I^3 eq QQ!1244819524/887503681;
        end if;
    end for;
    print "PlaneBoundaryCertificate: passed";
end procedure;

/* ----------------------------------------------------------------------- */
/* Type I lift: flat-limit Hilbert polynomial and infinitesimal stabilizer. */
/* ----------------------------------------------------------------------- */

function DegreeTwoMonomials(vars)
    mons := [];
    for i in [1..#vars] do
        for j in [i..#vars] do
            Append(~mons, vars[i]*vars[j]);
        end for;
    end for;
    return mons;
end function;

function CoefficientVector(f, mons)
    return Vector(QQ,[ MonomialCoefficient(f,m) : m in mons ]);
end function;

procedure TypeIFlatLimitCertificate()
    S<t,X0,X1,X2,X3,Z> := PolynomialRing(QQ,6,"grevlex");
    IFam := ideal<S |
      Z*X1-t*(X0*X1-2*X2*X3),
      Z*X2+t*(X0*X2-2*X1^2),
      Z*X3+t*(X0*X3-2*X1*X2),
      Z^2-t^2*(X0^2-4*X1*X3),
      X1*X3-X2^2,
      X0*X2-X1^2-X3^2>;
    Sat := Saturation(IFam,ideal<S | t>);

    P<Y0,Y1,Y2,Y3,Y4> := PolynomialRing(QQ,5,"grevlex");
    sp0 := hom<S -> P | [0,Y0,Y1,Y2,Y3,Y4]>;
    I0 := ideal<P | [ sp0(f) : f in Basis(Sat) ]>;
    expected := ideal<P |
      Y4*Y1,Y4*Y2,Y4*Y3,Y4^2,
      Y1*Y3-Y2^2,Y0*Y2-Y1^2-Y3^2>;
    assert I0 eq expected;
    HP0, reg0 := HilbertPolynomial(I0);
    d0 := Parent(HP0).1;
    assert HP0 eq 4*d0+1;

    sp1 := hom<S -> P | [1,Y0,Y1,Y2,Y3,Y4]>;
    I1 := ideal<P | [ sp1(f) : f in Basis(Sat) ]>;
    HP1, reg1 := HilbertPolynomial(I1);
    d1 := Parent(HP1).1;
    assert HP1 eq 4*d1+1;

    Rt<tt> := PolynomialRing(QQ);
    C := Matrix(Rt,[
      [1,0,0,0,1],
      [0,1,0,0,0],
      [0,0,1,0,0],
      [0,0,0,1,0],
      [tt,0,0,0,-tt]
    ]);
    assert Determinant(C) eq -2*tt;
    print "TypeIFlatLimitCertificate: passed";
end procedure;

procedure TypeIStabilizerCertificate()
    P<X0,X1,X2,X3,Z> := PolynomialRing(QQ,5,"grevlex");
    vars := [X0,X1,X2,X3,Z];
    mons := DegreeTwoMonomials(vars);
    W := [
      Z*(X2-X1),
      Z*(X3-2*X1),
      Z*(Z-3*X1),
      3*X1*X3-3*X2^2-X1*Z,
      3*X0*X2-3*X1^2-3*X3^2+4*X1*Z
    ];
    WM := Matrix(QQ,[ Eltseq(CoefficientVector(f,mons)) : f in W ]);
    assert Rank(WM) eq 5;
    ann := Basis(NullspaceOfTranspose(WM));
    assert #ann eq 10;

    columns := [];
    for a in [1..5] do
        for b in [1..5] do
            col := [];
            for f in W do
                action := vars[b]*Derivative(f,a);
                cv := CoefficientVector(action,mons);
                col cat:= [ &+[cv[k]*n[k] : k in [1..15]] : n in ann ];
            end for;
            Append(~columns,col);
        end for;
    end for;

    Lentries := [];
    for i in [1..50] do
        for j in [1..25] do
            Append(~Lentries, columns[j][i]);
        end for;
    end for;
    L := Matrix(QQ,50,25,Lentries);
    assert Rank(L) eq 24;
    print "TypeIStabilizerCertificate: passed";
end procedure;

/* ----------------------------------------------------------------------- */
/* Determinantal and Pfaffian models in (H1) and (H2).                      */
/* ----------------------------------------------------------------------- */

function PrincipalPfaffians(M)
    assert Nrows(M) eq 5 and Ncols(M) eq 5;
    ans := [];
    for k in [1..5] do
        inds := [ i : i in [1..5] | i ne k ];
        entries := [];
        for i in inds do
            for j in inds do
                Append(~entries,M[i,j]);
            end for;
        end for;
        Mk := Matrix(BaseRing(M),4,4,entries);
        Append(~ans,Pfaffian(Mk));
    end for;
    return ans;
end function;

procedure DeterminantalPfaffianCertificates()
    P<x,y,u,v,w> := PolynomialRing(QQ,5,"grevlex");

    N := Matrix(P,[ [u,-v,-y-w,x], [0,u,v,w] ]);
    mins := [];
    for i in [1..4] do
        for j in [i+1..4] do
            Append(~mins,N[1,i]*N[2,j]-N[1,j]*N[2,i]);
        end for;
    end for;
    IC := ideal<P | mins>;
    HPC, regC := HilbertPolynomial(IC);
    dC := Parent(HPC).1;
    assert HPC eq 4*dC+1;
    assert ideal<P | mins> eq ideal<P |
        u^2,u*v,u*w,u*y-v^2+u*w,u*x+v*w,x*v+y*w+w^2>;

    M5 := Matrix(P,[
      [0,0,0,u,v],
      [0,0,u,-w,x],
      [0,-u,0,v,y],
      [-u,w,-v,0,0],
      [-v,-x,-y,0,0]
    ]);
    pf5 := PrincipalPfaffians(M5);
    I5 := ideal<P | pf5>;
    HP5, reg5 := HilbertPolynomial(I5);
    d5 := Parent(HP5).1;
    assert HP5 eq 5*d5;
    assert ideal<P | pf5> eq ideal<P |
        u^2,u*v,u*x+v*w,u*y-v^2,x*v+y*w>;

    Q<x0,x1,x2,uu,vv> := PolynomialRing(QQ,5,"grevlex");
    f := x0*x2-x1^2;
    IY := ideal<Q |
      f,x0*uu+x1*vv,x1*uu+x2*vv,uu^2,uu*vv,vv^2>;
    HPY, regY := HilbertPolynomial(IY);
    dY := Parent(HPY).1;
    assert HPY eq 4*dY+1;

    MH2 := Matrix(Q,[
      [0,0,x1,-x0,vv],
      [0,0,x2,-x1,-uu],
      [-x1,-x2,0,uu,0],
      [x0,x1,-uu,0,0],
      [-vv,uu,0,0,0]
    ]);
    pfH2 := PrincipalPfaffians(MH2);
    IH2 := ideal<Q | pfH2>;
    HPH2, regH2 := HilbertPolynomial(IH2);
    dH2 := Parent(HPH2).1;
    assert HPH2 eq 5*dH2;
    assert ideal<Q | pfH2> eq ideal<Q |
      uu^2,uu*vv,x0*uu+x1*vv,x1*uu+x2*vv,f>;

    print "DeterminantalPfaffianCertificates: passed";
end procedure;

/* ----------------------------------------------------------------------- */
/* Resultants in (Q) and final image relation in (R).                       */
/* ----------------------------------------------------------------------- */

procedure QuadraticBranchCertificates()
    A<alpha,y4,v,E,G> := PolynomialRing(QQ,5);
    U<r> := PolynomialRing(A);
    discarded := r^2+E;
    inverseeq := alpha*r^2-y4*v^2*r+alpha*E-y4*G;
    res := Resultant(discarded,inverseeq);
    assert res eq y4^2*(G^2+v^4*E);
    disc := y4^2*v^4-4*alpha*(alpha*E-y4*G);
    assert disc eq y4^2*v^4+4*alpha*y4*G-4*alpha^2*E;

    C<g0,g1,g2,g3,e0,e1,e2,u0,v0> := PolynomialRing(QQ,9);
    GG := g0*u0^3+g1*u0^2*v0+g2*u0*v0^2+g3*v0^3;
    EE := e0*u0^2+e1*u0*v0+e2*v0^2;
    identity := GG^2+v0^4*EE;

    // If identity vanishes, specialization at v0=0 gives g0=0.
    psiV0 := hom<C -> C | [g0,g1,g2,g3,e0,e1,e2,u0,0]>;
    assert psiV0(identity) eq g0^2*u0^6;

    // After g0=0, factor v0^2 and specialize the quotient at v0=0;
    // this gives g1=0.
    phi0 := hom<C -> C | [0,g1,g2,g3,e0,e1,e2,u0,v0]>;
    H0 := (g1*u0^2+g2*u0*v0+g3*v0^2)^2+v0^2*EE;
    assert phi0(identity) eq v0^2*H0;
    assert psiV0(H0) eq g1^2*u0^4;

    // After g0=g1=0, the identity is exactly
    // v0^4*((g2*u0+g3*v0)^2+EE).
    phi01 := hom<C -> C | [0,0,g2,g3,e0,e1,e2,u0,v0]>;
    assert phi01(identity) eq v0^4*((g2*u0+g3*v0)^2+EE);

    P<x0,x1,u,vv,z> := PolynomialRing(QQ,5,"grevlex");
    assert u*vv^2+vv*(-u*vv) eq 0;
    assert (-u)*(-u*vv)+(-vv)*u^2 eq 0;

    cases := [ <1,0,1>, <0,0,1>, <1,1,0> ];
    for triple in cases do
        r0 := QQ!triple[1];
        r1 := QQ!triple[2];
        r2 := QQ!triple[3];
        q0 := x1*u+r2*z^2;
        q1 := -x0*u+x1*vv-r1*z^2;
        q2 := -x0*vv+r0*z^2;
        q3 := z*u;
        q4 := z*vv;
        assert vv^2*q0-u*vv*q1+u^2*q2 eq
               z^2*(r0*u^2+r1*u*vv+r2*vv^2);
        I := ideal<P | q0,q1,q2,q3,q4>;
        HP, reg := HilbertPolynomial(I);
        d := Parent(HP).1;
        assert HP eq 3*d+4;
    end for;
    assert 0^2-4*1*1 ne 0;
    assert 0^2-4*0*1 eq 0;
    assert 1^2-4*1*0 ne 0;

    print "QuadraticBranchCertificates: passed";
end procedure;

procedure LowRankNormalFormCertificates()
    P<x0,x1,x2,x3,z> := PolynomialRing(QQ,5);
    q0 := x0*x2;
    q1 := x0*x3;
    q2 := x1*x2;
    q3 := x1*x3;
    assert q0*q3-q1*q2 eq 0;
    print "LowRankNormalFormCertificates: passed";
end procedure;

procedure RunChecks()
    PlaneBoundaryCertificate();
    TypeIFlatLimitCertificate();
    TypeIStabilizerCertificate();
    DeterminantalPfaffianCertificates();
    QuadraticBranchCertificates();
    LowRankNormalFormCertificates();
    print "All certificates passed.";
end procedure;

print "Cremona_P4_checks.m loaded successfully.";
print "Call RunChecks();";

/* Uncomment the following line to run all checks immediately after loading. */
/* RunChecks(); */
