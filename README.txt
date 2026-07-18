## Magma computations

This repository contains two independent Magma scripts accompanying the paper *Quadratic Cremona transformations of (\mathbb P^4)*.

### `P4_cremona_all_types.m`

This script contains explicit quadratic rational maps

[
\varphi:\mathbb P^4\dashrightarrow\mathbb P^4
]

defined by five quadratic forms. The examples represent the geometric types and multidegrees considered in the paper.

For each example, the script computes the projective degrees

[
[d_0,d_1,d_2,d_3,d_4]
]

by residual intersections. For fixed (k), it intersects (k) general members of the quadratic linear system with (4-k) general hyperplanes and removes the base locus by saturation. The degree of the resulting zero-dimensional residual scheme is the (k)-th projective degree.

The computations are repeated several times with random choices, and the most frequent result is returned. A list

[
[1,2,d_2,d_3,1]
]

shows that the map has degree one and hence is birational, with abbreviated multidegree ((2,d_2,d_3)).

The script includes examples of:

* punctual transformations;
* reduced and non-reduced line types;
* de Jonquières transformations;
* conic and ribbon types;
* two skew lines;
* twisted cubics;
* a conic together with a line;
* rational normal quartics;
* elliptic normal quintics;
* all the multidegrees occurring in the paper.

The computations are performed over the finite field (\operatorname{GF}(32003)). The random seed is fixed, so the calculations are reproducible. The field can be replaced by another sufficiently large finite field.

To run all examples, use

```magma
load "P4_cremona_all_types.m";
```

The file automatically calls

```magma
RunAll(: Trials := 6, TryInverse := false, PrintForms := false);
```

The optional parameters have the following meaning:

* `Trials`: number of random residual intersections used for each projective degree;
* `TryInverse`: asks Magma to test invertibility and, when possible, compute the inverse;
* `PrintForms`: prints the five quadrics defining each example.

### `Cremona_P4_checks.m`

This script contains exact algebraic certificates for the computations used in the proofs of the paper. All calculations are performed over (\mathbb Q).

Unlike `P4_cremona_all_types.m`, this file does not use random choices. Each procedure verifies a finite identity or algebraic calculation and stops with an assertion if the expected result fails.

The script contains the following certificates.

#### `PlaneBoundaryCertificate()`

This procedure verifies the degeneration used to place systems with a plane in the base locus in the boundary of the component associated with two skew lines. It checks:

* the Smith form of the matrix of conditions imposed on quadrics;
* the five-dimensional limiting linear system;
* the maximal minors of the residual matrices;
* the two residual binary quartics;
* the invariants (I), (J), and (J^2/I^3);
* the non-constancy of the cross-ratio in the one-parameter family.

#### `TypeIFlatLimitCertificate()`

This procedure verifies the flat degeneration used for the type I lift of a quadratic Cremona transformation of (\mathbb P^3). It checks:

* saturation with respect to the deformation parameter;
* the ideal of the special fiber;
* the Hilbert polynomial (4m+1) of both the special and general fibers;
* the linear independence of the five parametrizing sections away from the special fiber.

#### `TypeIStabilizerCertificate()`

This procedure computes the infinitesimal stabilizer of the explicit type I system. It verifies that the linear map

[
\mathfrak{gl}_5\longrightarrow
\operatorname{Hom}
\left(
W,\operatorname{Sym}^2(V^*)/W
\right)
]

has rank (24). Thus its kernel consists only of scalar matrices, and the corresponding stabilizer in (\operatorname{PGL}_5) is finite.

#### `DeterminantalPfaffianCertificates()`

This procedure checks the determinantal and Pfaffian models occurring in the branches ((H1)) and ((H2)). It verifies:

* the (2\times2) minors of the determinantal quartic model;
* its Hilbert polynomial (4m+1);
* the five principal Pfaffians of the elliptic quintic model;
* its Hilbert polynomial (5m);
* the Pfaffian model arising from the double-conic configuration;
* the equality between the computed ideals and the ideals used in the paper.

#### `QuadraticBranchCertificates()`

This procedure verifies the algebraic identities used in the quadratic discriminant branch ((Q)). It checks:

* the resultant
  [
  \operatorname{Res}_r
  \bigl(r^2+E,,
  \alpha r^2-y_4v^2r+\alpha E-y_4G
  \bigr)
  ======

  y_4^2(G^2+v^4E);
  ]
* the discriminant of the residual quadratic equation;
* the coefficient comparisons implying the required divisibility of (G);
* the syzygies of the normal form;
* the Hilbert polynomial (3m+4) in the reduced, double-root, and degenerate residual cases.

#### `LowRankNormalFormCertificates()`

This procedure verifies the final equation used in the low-rank branch ((R)). For the normal form

[
x_0x_2,\quad x_0x_3,\quad x_1x_2,\quad x_1x_3,\quad z^2,
]

it checks the identity

[
(x_0x_2)(x_1x_3)-(x_0x_3)(x_1x_2)=0.
]

Hence the image is contained in the quadric

[
y_0y_3-y_1y_2=0
]

and the rational map is not dominant.

To run all exact certificates, use

```magma
load "Cremona_P4_checks.m";
RunChecks();
```

Each procedure prints a confirmation when its checks pass. The program stops with an assertion if any certificate fails.

The two scripts have complementary purposes:

* `P4_cremona_all_types.m` computes the projective degrees of the explicit Cremona transformations over a large finite field;
* `Cremona_P4_checks.m` verifies over (\mathbb Q) the exact algebraic calculations used in the proofs.
