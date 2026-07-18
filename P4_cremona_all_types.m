// ------------------------------------------------------------
// Description of the script.
// ------------------------------------------------------------
//
// This script checks explicit quadratic rational maps
//
//        phi : P^4 ---> P^4
//
// defined by five quadratic forms.  For each example appearing
// in the paper, the script computes the projective degrees
//
//        [d0,d1,d2,d3,d4]
//
// by residual intersections.  In particular, if the computed list is
//
//        [1,2,d2,d3,1],
//
// then the map is generically finite of degree one, hence birational,
// and has abbreviated multidegree (2,d2,d3).
//
// The computation of projective degrees is done as follows.  For fixed k,
// we intersect k general members of the given linear system with 4-k
// general linear forms in P^4.  This gives a zero-dimensional scheme
// after removing the base locus by saturation.  Its degree is the k-th
// projective degree.  Since the choices are random, the computation is
// repeated several times; the most frequent value is returned.
//
// The main functions are:
//
//   RandomLinearForm()
//      returns a random non-zero linear form in x0,...,x4.
//
//   RandomLinearCombination(forms)
//      returns a random non-zero linear combination of the given forms.
//
//   SaturateByIdeal(I,J)
//      computes the saturation of I with respect to J by iterated
//      colon ideals.  This removes the fixed base locus from a residual
//      intersection.
//
//   ResidualDegree(forms,k,Ibase)
//      computes one random residual intersection degree contributing
//      to the k-th projective degree.
//
//   ProjectiveDegrees(forms : Trials := 6)
//      computes the full list [d0,d1,d2,d3,d4] of projective degrees.
//      The optional parameter Trials controls how many random residual
//      intersections are used for each k.
//
//   TryInverseMap(forms)
//      asks Magma to test whether the corresponding rational map is
//      invertible and, if possible, computes the inverse.  This step is
//      optional and may be slow for some examples.
//
//   CheckExample(name, forms, expected)
//      prints the dimension and degree of the saturated base scheme,
//      computes the projective degrees, and compares them with the
//      expected list.
//
//   FivePfaffians(...)
//      returns the five 4x4 Pfaffians of a 5x5 skew-symmetric matrix.
//      This is used for the elliptic normal quintic example.
//
//   Example_abc()
//      returns the five quadratic forms defining the example with
//      abbreviated multidegree (a,b,c).  When two different geometric
//      types have the same multidegree, the function name records the
//      type.
//
//   RunAll()
//      runs all checks.  The optional parameters are:
//        Trials     : number of random trials for each projective degree;
//        TryInverse : if true, Magma tries IsInvertible and Inverse;
//        PrintForms : if true, the five quadrics are printed.
//
// The script is written over the finite field GF(32003).  This is enough
// to certify the expected projective degrees for the displayed explicit
// examples.  The same computations can be repeated over another large
// finite field by changing the definition of K.
// ------------------------------------------------------------

SetSeed(12345);

K := GF(32003);
P4<x0,x1,x2,x3,x4> := ProjectiveSpace(K,4);
R := CoordinateRing(P4);
vars := [x0,x1,x2,x3,x4];

// ------------------------------------------------------------
// Random auxiliaries.
// ------------------------------------------------------------

function RandomLinearForm()
    repeat
        f := &+[Random(K)*vars[i] : i in [1..5]];
    until f ne 0;
    return f;
end function;

function RandomLinearCombination(forms)
    repeat
        f := &+[Random(K)*forms[i] : i in [1..#forms]];
    until f ne 0;
    return f;
end function;

function SaturateByIdeal(I,J)
    old := I;
    for t in [1..30] do
        new := ColonIdeal(old,J);
        if new eq old then
            return new;
        end if;
        old := new;
    end for;
    return old;
end function;

function MostFrequent(vals)
    U := Setseq(Seqset(vals));
    best := U[1];
    bestc := 0;
    for v in U do
        c := #[w : w in vals | w eq v];
        if c gt bestc then
            best := v;
            bestc := c;
        end if;
    end for;
    return best;
end function;

function ResidualDegree(forms,k,Ibase)
    eqs := [RandomLinearCombination(forms) : i in [1..k]] cat
           [RandomLinearForm() : i in [1..4-k]];
    I := ideal<R | eqs>;
    Ires := SaturateByIdeal(I,Ibase);
    X := Scheme(P4,Ires);
    if Dimension(X) ne 0 then
        return -1;
    end if;
    return Degree(X);
end function;

function ProjectiveDegrees(forms : Trials := 6)
    Ibase := Saturation(ideal<R | forms>);
    degs := [1];
    for k in [1..4] do
        vals := [];
        for t in [1..Trials] do
            d := ResidualDegree(forms,k,Ibase);
            if d ge 0 then
                Append(~vals,d);
            end if;
        end for;
        if #vals eq 0 then
            Append(~degs,-1);
        else
            if #Seqset(vals) gt 1 then
                print "  warning: non-constant residual degrees for k =", k, vals;
            end if;
            Append(~degs,MostFrequent(vals));
        end if;
    end for;
    return degs;
end function;

procedure TryInverseMap(forms)
    phi := map<P4 -> P4 | forms>;
    print "  IsInvertible:";
    try
        b := IsInvertible(phi);
        print " ", b;
        if b then
            psi := Inverse(phi);
            invforms := DefiningPolynomials(psi);
            print "  inverse degrees:", [Degree(g) : g in invforms];
        end if;
    catch e
        print "  IsInvertible/Inverse failed or did not return:";
        print e;
    end try;
end procedure;

procedure CheckExample(name, forms, expected : Trials := 6, TryInverse := false, PrintForms := false)
    print "============================================================";
    print name;
    print "expected projective degrees:", expected;
    if PrintForms then
        print "forms:";
        print forms;
    end if;
    Ibase := Saturation(ideal<R | forms>);
    B := Scheme(P4,Ibase);
    print "base dimension:", Dimension(B);
    print "base degree:", Degree(B);
    degs := ProjectiveDegrees(forms : Trials := Trials);
    print "computed projective degrees:", degs;
    if degs eq expected then
        print "status: OK";
    else
        print "status: CHECK";
    end if;
    if TryInverse then
        TryInverseMap(forms);
    end if;
end procedure;

// ------------------------------------------------------------
// Pfaffians for the elliptic quintic example.
// ------------------------------------------------------------

function FivePfaffians(a12,a13,a14,a15,a23,a24,a25,a34,a35,a45)
    return [
        a23*a45 - a24*a35 + a25*a34,
        a13*a45 - a14*a35 + a15*a34,
        a12*a45 - a14*a25 + a15*a24,
        a12*a35 - a13*a25 + a15*a23,
        a12*a34 - a13*a24 + a14*a23
    ];
end function;

// ------------------------------------------------------------
// Examples with d2 = 4.
// ------------------------------------------------------------

function Example_248()
    return [
        x2*x3,
        x2*x4,
        x3*x4,
        x1*(x3 + x4 - x2),
        x0*(x1 + x2 + x3 + x4) - x1*x2
    ];
end function;

function Example_247()
    return [
        x1*x2 - x0*x4,
        x1*x3 + 2*x1*x4 - 3*x0*x4,
        4*x2^2 - x1*x3 - 3*x0*x4,
        4*x2*x3 - 3*x1*x3 - x0*x4,
        4*x3^2 - x1*x3 - 3*x0*x4
    ];
end function;

// Reduced-line de Jonquieres type.  This is a second geometric
// family with projective degrees [1,2,4,7,1].  Its positive-dimensional
// base hull is a reduced line; the residual zero-dimensional scheme
// consists of a length-seven multiplicity-nine point and one reduced point.
function Example_247_reduced_line_dejonquieres()
    return [
        3*x0*x1 - 4*x0*x2 + x0*x4 - x1^2 + 2*x1*x2 - x1*x3,
        6*x0*x1 - 3*x0*x2 - 7*x1^2 + 4*x1*x2 + x1*x4 - x2*x3,
        39366*x0^2 + 28039*x0*x1 - 200188*x0*x2 - 738*x0*x3
            + 50047*x0*x4 - 66185*x1^2 + 124781*x1*x2
            - 12413*x1*x4 + 12552*x2^2 - 3138*x2*x4 - 4620*x3^2,
        9414*x0^2 + 59411*x0*x1 - 62252*x0*x2 + 3138*x0*x3
            + 1703*x0*x4 - 98665*x1^2 + 18049*x1*x2
            + 32003*x1*x4 + 14088*x2^2 - 3522*x2*x4 - 4620*x3*x4,
        10566*x0^2 - 53801*x0*x1 + 17972*x0*x2 + 3522*x0*x3
            - 4493*x0*x4 + 24955*x1^2 - 21499*x1*x2
            + 6247*x1*x4 - 38568*x2^2 + 28122*x2*x4 - 4620*x4^2
    ];
end function;

function Example_246_conic()
    return [
        x0*x2 - x1^2,
        x1*x3,
        x1*x4,
        x2*x3 - x3^2,
        x2*x4 - x4^2
    ];
end function;

// Genus -2 ribbon type.  This is a third geometric family with
// projective degrees [1,2,4,6,1].  Its one-dimensional base component
// is a degree-two locally complete-intersection ribbon supported on a line,
// and the residual base scheme consists of three reduced points.
function Example_246_genus_minus_two_ribbon()
    return [
        x0*x2 + x1*x3 + x2*x4 + x3^2,
        x0*x3 + x1*x4 + x2^2 + x4^2,
        2*x2^2 - 7*x2*x3 - x2*x4 + 8*x3^2,
        10*x2^2 - 19*x2*x3 - 5*x2*x4 + 8*x3*x4,
        18*x2^2 - 15*x2*x3 - 25*x2*x4 + 8*x4^2
    ];
end function;

function Example_246_two_skew_lines()
    return [
        x0*x2 + 2*x0*x4 - 2*x2*x4,
        x0*x3 + 3*x0*x4 - 3*x2*x4,
        x1*x2 + 2*x0*x4 - 2*x2*x4,
        x1*x3 + 3*x0*x4 - 3*x2*x4,
        x4*(x1 - x3 - 2*x0 + 2*x2)
    ];
end function;

function Example_245_twisted_cubic()
    return [
        x0*x2 - x1^2,
        x0*x3 - x1*x2,
        x1*x3 - x2^2,
        x0*x4,
        x3*x4
    ];
end function;

function Example_245_conic_plus_line()
    return [
        -4*x0*x2 + x0*x3 + 2*x0*x4,
        -3*x0*x2 + 2*x0*x4 + x1*x2,
        -4*x0*x2 + 2*x0*x4 + x1*x3,
        -4*x0*x2 + 3*x0*x4 + x1*x4,
        8*x0*x2 - 5*x0*x4 + x2*x4 - x3^2
    ];
end function;

function Example_244_nonreduced_line()
    return [
        x0*x3,
        x0*x4,
        x0*x1 - x1^2,
        x0*x2 - x2^2,
        x1*x2
    ];
end function;

function Example_244_rational_quartic()
    return [
        x0*x3 - x1*x2,
        x1*x3 - x2^2,
        x2*x4 - x3^2,
        x0*x2 - x1^2 + x0*x4 - x1*x3,
        x0*x4 - x1*x3 - x1*x4 + x2*x3
    ];
end function;

function Example_243_elliptic_quintic()
    return FivePfaffians(
        x0, x1, x2, x3,
        x2, x3, x4,
        x4, x0 + x1, x1 + x2
    );
end function;

function Example_242()
    return [
        x0*x1,
        x0^2,
        x3^2 - x1*x2,
        x1*x4,
        x0*x3
    ];
end function;

// ------------------------------------------------------------
// Examples with d2 = 3 and d2 = 2.
// ------------------------------------------------------------

function Example_234()
    return [
        x0*x3 - x3*x4 + x1*x4,
        x0*x4 - x3*x4 + x1*x4,
        x1*x3 + x1*x4,
        x2*x3,
        x2*x4
    ];
end function;

function Example_233()
    return [
        x1*x3 + 2*x2*x3 - 3*x4^2,
        2*x2*x3 + x3*x4 - 3*x4^2,
        2*x2*x3 + x0*x4 - 3*x4^2,
        x1*x4 - x4^2,
        3*x2*x3 + x2*x4 - 4*x4^2
    ];
end function;

function Example_232()
    return [
        x0*x1,
        x0*x2,
        x1*x2,
        x2*x3,
        x1*x4
    ];
end function;

function Example_222()
    return [
        x1^2 + x2^2 + x3^2 + x4^2,
        x0*x1,
        x0*x2,
        x0*x3,
        x0*x4
    ];
end function;

// ------------------------------------------------------------
// Run all checks.
// ------------------------------------------------------------

procedure RunAll(: Trials := 6, TryInverse := false, PrintForms := false)
    CheckExample("(2,4,8), punctual", Example_248(), [1,2,4,8,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,4,7), original type", Example_247(), [1,2,4,7,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,4,7), reduced-line de Jonquieres type", Example_247_reduced_line_dejonquieres(), [1,2,4,7,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,4,6), conic type", Example_246_conic(), [1,2,4,6,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,4,6), genus -2 ribbon type", Example_246_genus_minus_two_ribbon(), [1,2,4,6,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,4,6), two skew lines type", Example_246_two_skew_lines(), [1,2,4,6,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,4,5), twisted cubic type", Example_245_twisted_cubic(), [1,2,4,5,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,4,5), conic plus line type", Example_245_conic_plus_line(), [1,2,4,5,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,4,4), non-reduced line type", Example_244_nonreduced_line(), [1,2,4,4,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,4,4), rational quartic type", Example_244_rational_quartic(), [1,2,4,4,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,4,3), elliptic normal quintic", Example_243_elliptic_quintic(), [1,2,4,3,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,4,2)", Example_242(), [1,2,4,2,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,3,4)", Example_234(), [1,2,3,4,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,3,3)", Example_233(), [1,2,3,3,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,3,2)", Example_232(), [1,2,3,2,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
    CheckExample("(2,2,2)", Example_222(), [1,2,2,2,1] : Trials := Trials, TryInverse := TryInverse, PrintForms := PrintForms);
end procedure;

RunAll(: Trials := 6, TryInverse := false, PrintForms := false);
