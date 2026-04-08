# Border Waring Rank Computations

Macaulay2 code for computing canonical forms, apolar ideals, and multiplication
operators for local Artinian Gorenstein algebras from Casnati's classification.

Companion code for the paper *Debordering Concise Forms* by Aviv, Rafael, and Amir.

## Requirements

- [Macaulay2](https://macaulay2.com/) (tested with version 1.25)
- Macaulay2 packages: `SpechtModule`, `CorrespondenceScrolls`

## Usage

Edit the parameters at the top of `Main.m2`:

```m2
idealDegree = 5;    -- e from Casnati's classification
numVars = 1;         -- n (number of variables)
idealNumber = 1;     -- ideal type (1-26)
```

Run:

```bash
M2 --script Main.m2
```

Output:
- `output/example.tex` -- LaTeX for the paper (form, apolar ideal, etc.)
- `output/draft.tex` -- Multiplication operators (not included in paper)

## Algebra Types

The 26 ideal types correspond to Casnati's classification of local Gorenstein
algebra isomorphism types of multiplicity at most 9.

## File Structure

- `Main.m2` -- driver script, produces LaTeX output
- `HelperFunctions.m2` -- apolar pairing, Hankel operators, multiplication operators
- `GenerateGorensteinAlgebra.m2` -- algebra generation from Casnati classification
- `ExtractBorderForms.m2` -- canonical form extraction
- `Parameters.m2` -- prime characteristic p
- `Ideals/` -- one file per Casnati algebra type (1-26)
