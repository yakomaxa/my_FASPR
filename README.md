K.S.: Substantial modifications (topN rotamer and Chi2 limitation) were introduced and are under validation.

# INTRODUCTION

FASPR is a fast and accurate program for protein side-chain packing, which is an important step in conventional, energy-based protein structure prediction and protein design.

# USAGE

<b> $path/FASPR -i input.pdb -o output.pdb [-s sequence.txt] [-n N] [-n CHAIN:POS[ICODE]:N ...] [-m sitemap.txt] [-a CHI2MIN:CHI2MAX[:trp]] </b>

-i: the input backbone in pdb format for packing. There should be no missing main-chain atoms (N, CA, C and O). If side-chain atoms are included in input.pdb, they are ignored by FASPR.

-o: the repacked model output in pdb format by FASPR. The residue positions are kept identical to the input pdb file.

-s (optional): an input sequence to be repacked on the input protein backbone. The sequence should be written as one-single line of one-letter alphabet of amino acid types, e.g., ACDEFGHIKLMNPQRSTVWYYWVTSRQPNMLKIHGFEDCA. Only 20 canonical amino-acid types are allowed in the sequence. When the input sequence is the same as the one extracted from the pdb structure, the amino acid side-chain conformations are repacked only. When the input sequence is different, mutations will be introduced. Thus, <b>FASPR can be used to build mutant models efficiently</b>. If you want to fix the conformation of some residues during packing, you can specify them in lower-case letters while keeping other residues in upper-case, e.g., acdefghiklmnpqrstvwyYWVTSRQPNMLKIHGFEDCA.

-n N (optional): restrict packing to at most the N highest-probability rotamers per residue, ranked by their Dunbrack backbone-dependent probability for that residue's (phi, psi) context. N must be a positive integer. If a residue has fewer than N candidate rotamers, all of them are kept. Omitting -n preserves the default FASPR behavior (no restriction). This can be used to trade some accuracy for faster packing on large structures.

-n CHAIN:POS[ICODE]:N (optional, repeatable): set a rotamer-count limit for one specific residue instead of the whole structure, e.g. -n A:15:3 limits chain A residue 15 to its 3 highest-probability rotamers. Use an insertion code directly after the position if needed (e.g. A:15A:3). Use "_" as CHAIN to mean a blank/space chain id. This form can be repeated to set limits for several residues, and can be combined with a single plain -n N to also set a default limit for every residue not explicitly listed. A site that doesn't match any residue in the input structure produces a warning and is ignored (packing still proceeds). This inline form cannot be combined with -m (see below) in the same run.

-m sitemap.txt (optional): an alternative to repeated -n CHAIN:POS:N flags for setting many per-residue rotamer limits at once. The file has one override per line: `CHAIN POS[ICODE] N`, e.g.:

```
# chain  resi  N
A        15    3
A        22A   1
```

Blank lines and lines starting with # are ignored. As with the inline form, -m can be combined with a single plain -n N for the default limit on unlisted residues, but -m and inline -n CHAIN:POS:N overrides cannot be used together in the same run (choose one mechanism).

Note: a per-site -n/-m override has no effect on a residue that is fixed via lower-case in -s (its conformation is never touched by packing) or on ALA/GLY (no side-chain rotamers to restrict). FASPR still matches the override to the residue and prints a warning explaining why it was a no-op, rather than silently ignoring it.

-a CHI2MIN:CHI2MAX[:trp] (optional): FASPR's equivalent of Rosetta's `LimitAromaChi2` task operation. By default FASPR (like Rosetta's Dunbrack energy term) will happily use PHE, TYR and HIS rotamers whose chi2 is far from 90 degrees, which Rosetta's own documentation notes are rarely physically realistic. `-a` makes the rotamer selection INCLUSIVE instead: only PHE/TYR/HIS rotamers whose chi2 falls in `[CHI2MIN,CHI2MAX]` are kept as candidates; everything else is discarded before packing. The usual default window is `70:110` (e.g. `-a 70:110`). TRP is excluded from the restriction by default (its energy landscape is more permissive for this), add a third field of exactly `trp` (e.g. `-a 70:110:trp`) to also restrict TRP. CHI2MIN/CHI2MAX may be real numbers, and must satisfy `0 <= CHI2MIN < CHI2MAX <= 180`. Before comparing, a rotamer's chi2 is folded into `[0,180)` by adding 180 if it is negative (this folding is only used for the accept/reject test -- it does not alter the rotamer's actual stored chi2 or its built geometry). If, for some residue, no rotamer in its Dunbrack-library candidate set survives the window, FASPR falls back to keeping that residue's full original candidate set rather than leaving it with zero rotamers. Omitting `-a` entirely preserves the default FASPR behavior (no restriction). `-a` composes with `-n`/`-m`: the chi2 filter is applied first, and any top-N restriction is then applied to the surviving, chi2-acceptable rotamers.

Please remember to put the executable FASPR program and the binary Dunbrack rotamer library 'dun2010bbdep.bin' in the same directory. Otherwise, the program will report "error! cannot find rotamer library dun2010bbdep.bin". Also, do not change the name of 'dun2010bbdep.bin' because it is hard-coded in the source code.

Examples:

Use at most the 5 highest-probability rotamers per residue:
<b> ./FASPR -i input.pdb -o output.pdb -n 5 </b>

Limit only chain A residue 15 to 1 rotamer and chain A residue 20 to 2 rotamers, with 10 as the default elsewhere:
<b> ./FASPR -i input.pdb -o output.pdb -n 10 -n A:15:1 -n A:20:2 </b>

Same, using a site-map file instead:
<b> ./FASPR -i input.pdb -o output.pdb -n 10 -m sitemap.txt </b>

Reject unrealistic aromatic chi2 rotamers (PHE/TYR/HIS only, Rosetta LimitAromaChi2 defaults):
<b> ./FASPR -i input.pdb -o output.pdb -a 70:110 </b>

Same, also restricting TRP, combined with a global top-5 rotamer limit:
<b> ./FASPR -i input.pdb -o output.pdb -a 70:110:trp -n 5 </b>

# INSTALLATION
We recommend users to download the FASPR source-code package to your computer and build the FASPR executable on your own. After downloading and unzipping the package, change into the $path/FASPR/ directory and run "<b>g++ -O3 --fast-math -o FASPR src/*.cpp</b>" if you are working on UNIX or Linux. For Mac users, use "-fast-math" or ignore it. If you are working on the Windows system, you need to install the g++ compiler first.

# COPYRIGHT & CONTACT
Copyright (c) Xiaoqiang Huang. FASPR is free to academic users. For suggestions, please contact xiaoqiah@umich.edu or xiaoqiah@outlook.com.

# REFERENCES
1. Huang x, Pearce R, Zhang Y, FASPR: an open-source tool for fast and accurate protein side-chain packing. Bioinformatics (2020) 36: 3758-3765.
2. Huang X, Pearce R, Zhang Y. Toward the Accuracy and Speed of Protein Side-Chain Packing: A Systematic Study on Rotamer Libraries. J. Chem. Inf. Model. 2020; 60:410-420.
