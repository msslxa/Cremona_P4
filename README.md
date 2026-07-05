// ------------------------------------------------------------
// Description of the Magma script.
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
