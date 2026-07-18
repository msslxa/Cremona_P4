# MAGMA COMPUTATIONS

This repository contains two independent Magma scripts accompanying the paper

Quadratic Cremona transformations of P^4

The two files have complementary purposes:

P4_cremona_all_types.m

computes the projective degrees of the explicit quadratic Cremona transformations appearing in the paper, while

Cremona_P4_checks.m

verifies the exact algebraic calculations used in some of the proofs.

1. P4_cremona_all_types.m

---
This script checks explicit quadratic rational maps

phi : P^4 ---> P^4

defined by five quadratic forms. For each example appearing in the paper, it computes the projective degrees

[d0,d1,d2,d3,d4]

by residual intersections.

For a fixed integer k, the script intersects k general members of the quadratic linear system with 4-k general hyperplanes in P^4. It then removes the base locus by saturation. The degree of the resulting zero-dimensional residual scheme is the k-th projective degree.

Since the intersections use random choices, the computation is repeated several times and the most frequent value is returned.

If the resulting list is

[1,2,d2,d3,1],

then the map is generically finite of degree one, hence birational, and has abbreviated multidegree

(2,d2,d3).

The script contains examples of the following geometric types:

- punctual transformations;
- reduced and non-reduced line types;
- de Jonquieres transformations;
- conic types;
- genus -2 ribbon types;
- two skew lines;
- twisted cubics;
- a conic together with a line;
- rational normal quartics;
- elliptic normal quintics;
- all multidegrees occurring in the paper.

The computations are performed over the finite field

GF(32003).

The random seed is fixed, so the calculations are reproducible. The same computations can be repeated over another sufficiently large finite field by changing the definition of the ground field.

The main functions are the following.

RandomLinearForm()

returns a random nonzero linear form in x0,...,x4.

RandomLinearCombination(forms)

returns a random nonzero linear combination of the forms in the input sequence.

SaturateByIdeal(I,J)

computes the saturation of I with respect to J by iterated colon ideals. This removes the fixed base locus from a residual intersection.

ResidualDegree(forms,k,Ibase)

computes one residual intersection contributing to the k-th projective degree.

ProjectiveDegrees(forms : Trials := 6)

computes the full list of projective degrees. The optional parameter Trials specifies the number of random residual intersections used for each degree.

TryInverseMap(forms)

asks Magma to test whether the corresponding rational map is invertible and, when possible, computes its inverse. This calculation is optional and can be slow.

CheckExample(name,forms,expected)

prints the dimension and degree of the saturated base scheme, computes the projective degrees, and compares them with the expected list.

FivePfaffians(...)

returns the five 4 by 4 Pfaffians of a 5 by 5 skew-symmetric matrix. It is used for the elliptic normal quintic example.

Example_abc()

returns the five quadrics defining an example with abbreviated multidegree (a,b,c). When distinct geometric types have the same multidegree, the function name also records the geometric type.

RunAll()

runs all examples.

The optional parameters of RunAll are:

Trials

the number of random computations performed for each projective degree;

TryInverse

if set to true, Magma calls IsInvertible and Inverse;

PrintForms

if set to true, the five defining quadrics are printed.

To run the script, use

load "P4_cremona_all_types.m";

The file automatically executes

RunAll(
    : Trials := 6,
      TryInverse := false,
      PrintForms := false
);

## 2. Cremona_P4_checks.m

This script contains exact algebraic certificates for calculations used in the proofs of the paper.

All computations are performed over the rational numbers. The script does not use random choices. Each procedure verifies a finite algebraic identity or calculation and stops with an assertion if the expected result fails.

To load the file and run all certificates, use

load "Cremona_P4_checks.m";
RunChecks();

Each procedure prints a confirmation when its checks pass.

## 2.1. PlaneBoundaryCertificate()

This procedure verifies the degeneration used to place systems having a plane in the base locus in the boundary of the component associated with two skew lines.

It checks:

- the Smith form of the matrix of conditions imposed on quadrics;
- the five-dimensional limiting linear system;
- the maximal minors of the residual matrices;
- the two residual binary quartics;
- the binary quartic invariants I and J;
- the two values of J^2/I^3;
- the non-constancy of the cross-ratio in the one-parameter family.

The two certified values are

47628/79507

and

1244819524/887503681.

Since these values are distinct, the cross-ratio is not constant along the family.

## 2.2. TypeIFlatLimitCertificate()

This procedure verifies the flat degeneration used for the type I lift of a quadratic Cremona transformation of P^3.

It checks:

- saturation with respect to the deformation parameter;
- the ideal of the special fiber;
- the Hilbert polynomial 4m+1 of the special fiber;
- the Hilbert polynomial 4m+1 of a general fiber;
- the linear independence of the five parametrizing sections away from the special fiber.

## 2.3. TypeIStabilizerCertificate()

This procedure computes the infinitesimal stabilizer of the explicit type I system.

It verifies that the linear map

gl_5 ---> Hom(W, Sym^2(V^*)/W)

has rank 24. Its kernel therefore consists only of scalar matrices. Hence the corresponding stabilizer in PGL_5 is finite.

## 2.4. DeterminantalPfaffianCertificates()

This procedure checks the determinantal and Pfaffian models occurring in the branches (H1) and (H2).

It verifies:

- the 2 by 2 minors of the determinantal quartic model;
- the Hilbert polynomial 4m+1 of the determinantal quartic;
- the five principal Pfaffians of the elliptic quintic model;
- the Hilbert polynomial 5m of the elliptic quintic;
- the Hilbert polynomial 4m+1 of the double-conic model;
- the five Pfaffians of the model arising from the double conic;
- the equality between the computed ideals and the ideals used in the paper.

## 2.5. QuadraticBranchCertificates()

This procedure verifies the algebraic identities used in the quadratic discriminant branch (Q).

It checks the resultant identity

Res_r(
    r^2+E,
    alpha*r^2-y4*v^2*r+alpha*E-y4*G
)
=
y4^2*(G^2+v^4*E).

It also checks:

- the discriminant of the residual quadratic equation;
- the coefficient comparisons implying the required divisibility of G;
- the successive specializations used to obtain g0=0 and g1=0;
- the factorization obtained after g0=g1=0;
- the syzygies of the normal form;
- the Hilbert polynomial 3m+4 in the reduced, double-root, and degenerate residual cases.

## 2.6. LowRankNormalFormCertificates()

This procedure verifies the final equation used in the low-rank branch (R).

For the normal form

x0*x2,
x0*x3,
x1*x2,
x1*x3,
z^2,

it checks the identity

(x0*x2)*(x1*x3) - (x0*x3)*(x1*x2) = 0.

Therefore the image is contained in the quadric

y0*y3-y1*y2 = 0,

and the corresponding rational map is not dominant.

3. Running the scripts

---
To compute the projective degrees of all explicit examples, use

load "P4_cremona_all_types.m";

To run all exact certificates used in the proofs, use

load "Cremona_P4_checks.m";
RunChecks();

The first script performs randomized residual-intersection computations over GF(32003). The second script performs deterministic exact computations over the rational numbers.
