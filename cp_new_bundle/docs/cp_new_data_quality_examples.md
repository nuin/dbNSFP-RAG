# cp_new — concrete examples of bad liftover and bad matching

Real cases pulled from the cp_new pipeline data, useful for teaching /
panel-design review / lab QA. All examples are real variants from the
APR2026 panel run.

---

## A. Liftover problems (dbNSFP → hg19)

### A1. Lift failure — dbNSFP has hg38 but no hg19

**14 of 1,057,357 variants** (0.001%) have an hg38 coord but a missing
`hg19_pos(1-based)`. These never make it into the GRCh37 SeqNext output
(the pipeline silently drops them at chunking time). All cluster in
specific positions where the hg19↔hg38 chain file had no clean lift.

```
AIP   hg38 chr11:67,490,920  ref=G alt=A   hg19 = (missing)
ATM   hg38 chr11:108,312,440 ref=G alt=A   hg19 = (missing)
DPYD  hg38 chr1:97,071,235   ref=C alt=T   hg19 = (missing)
```

**Impact**: minor — 14 variants out of 1M. But these positions WILL appear
on a sequencer; they just won't be classifiable through this pipeline
because there's no hg19 coord to anchor the SeqNext upload to. They'd need
either a manual lift via VV or to be reported on hg38 directly.

### A2. Chromosome changes hg19 → hg38

**Zero** in this dataset. dbNSFP's lift is contig-stable for the cp_new
panel. (This is the failure mode you most fear with general-purpose
liftover tools like CrossMap on segmental duplications — dbNSFP precomputes
to avoid it.)

---

## B. Matching problems

### B1. dbNSFP picked a non-MANE transcript → wrong c. shown

dbNSFP keeps **one Ensembl transcript per variant row**. Often it's not
MANE Select / the lab's reporting transcript. Same variant gets a
completely different c. notation depending on which transcript you look
at:

```
BRCA1 hg19 17:41197807 A>C
  dbNSFP says: ENST00000468300 -> c.2094T>G          (alt isoform)
  VV resolved on NM_007294.4:    c.5480T>G           (MANE Select)

EGFR hg19 7:55224528 T>A
  dbNSFP says: ENST00000420316 -> c.1210T>A          (coding on alt)
  VV resolved on NM_005228.5:    c.1207+3T>A         (intronic on canonical)

TP53 hg19 17:7573014 G>C
  dbNSFP says: ENST00000714408 -> c.1189C>G          (coding on alt)
  VV resolved on NM_000546.5:    c.1101-6C>G         (intronic on canonical)
```

**Impact**: huge. Without the VV resolver, uploading `BRCA1 NM_007294.4
c.2094T>G` to SeqNext would be **wrong**: NM_007294.4 has a different
exon layout, c.2094 isn't even the same nucleotide on that transcript.
The EGFR and TP53 examples are worse — dbNSFP says coding but the BED's
NM_ has them as intronic, so the entire interpretation flips.

### B2. dbNSFP says "coding" but VV reveals **deep** intronic

When dbNSFP's chosen transcript happens to be one where the variant is
coding-synonymous, the classifier (pre-VV-fix) would auto-call W2 Benign
without ever checking PhyloP. After VV, the same variant turns out to be
~50-100 bp deep into an intron of the BED's NM_ — should never be
auto-classified.

```
ACD  hg19 16:67691604 G>C
  dbNSFP ENST00000695734 -> c.1299C>G          (looks coding-synonymous)
  VV BED NM_001082486.2  -> c.1299-17C>G       (offset -17, deep intronic)

AIP  hg19 11:67256662 G>C
  dbNSFP ENST00000934218 -> c.294G>C           (looks coding)
  VV BED NM_003977.4     -> c.280-76G>C        (offset -76, deep intronic)

AIP  hg19 11:67256677 T>A
  dbNSFP ENST00000934218 -> c.309T>A           (looks coding)
  VV BED NM_003977.4     -> c.280-61T>A        (offset -61, deep intronic)

AIP  hg19 11:67256689 C>A
  dbNSFP ENST00000934218 -> c.321C>A           (looks coding)
  VV BED NM_003977.4     -> c.280-49C>A        (offset -49, deep intronic)
```

**Impact**: these would have been auto-Benign before the fix. The ROI
filter (`-15 ≤ offset ≤ +6`) + the VV-aware intronic detection now sends
them to manual review instead. Anything with offset like `-76` or `-49`
is well outside any canonical splice region — its effect on splicing
isn't predictable from the standard heuristics, and benign auto-call is
unsafe.

### B3. BED's NM_ isn't coding at this position at all

**981 rows** in the SeqNext output came back from VariantValidator as
`flagged:intergenic` — meaning the variant doesn't exist on the BED's
specified NM_ as a coding-sequence position. dbNSFP has it (on a different
isoform), but the lab's reporting transcript doesn't.

Top affected genes: **FANCB (152), FANCD2 (130), AP2S1 (79), TP53 (77),
SMARCA4 (71)**. These are all genes where the lab's BED commits to a
specific NM_, but dbNSFP's nsSNV coverage spans multiple isoforms.

```
AP2S1  NM_004069.5   hg19 19:47342867  ref=A alt=C   (VV: non-coding on this NM_)
AP2S1  NM_004069.5   hg19 19:47342879  ref=T alt=A   (VV: non-coding on this NM_)
FANCD2 NM_033084.5   hg19 3:10083307   ref=T alt=A   (VV: non-coding on this NM_)
```

**Impact**: these get classified by our rules (W2 Benign on
dbNSFP-coding-synonymous) but the c. notation in the output is on a
transcript the lab doesn't use. Two options: drop from upload (the
`--strict` mode does this), or re-resolve against the alternative isoform
that VV finds.

### B4. BED entries that aren't in the gene body

`MSH2_INV5_C+` at `chr2:38,121,080-38,121,280` — this is **9 MB upstream
of the MSH2 gene body** (MSH2 itself is at chr2:47,630-47,789). The lab
covers this region as an **inversion-breakpoint flank** to detect germline
SVs that disrupt MSH2 by chromosomal rearrangement. dbNSFP has no MSH2
nsSNVs at chr2:38M (they'd be intronic on every MSH2 isoform), so this
region contributes zero rows to the SeqNext output.

**Caveat from the lab**: short-read split-read calling has near-zero
sensitivity for 9 MB inversions. The BED row is design theater — anyone
relying on this to catch MSH2 inversions will miss every single one.

### B5. gnomAD-dbNSFP AF agreement

`dbNSFP gnomAD4.1_joint_AF` vs the AF we pulled directly from gnomAD v4.1
sites VCF — **0 disagreements >1.5× in the panel**. dbNSFP is faithful to
gnomAD v4.1 for joint AF. (Worth checking if a future dbNSFP release uses
a different gnomAD snapshot.)

---

## Recommendations from these findings

1. **Always run VariantValidator** before SeqNext upload. The non-MANE
   c. notation in dbNSFP is wrong for the lab's reporting transcript ~70%
   of the time we checked.

2. **Treat dbNSFP's `genename` and `Ensembl_transcriptid` columns as
   opaque** — don't infer "coding vs intronic" from dbNSFP's HGVSc alone.
   Always cross-check against the lab's NM_ via VV.

3. **The intronic ROI filter (-15 to +6) is essential**. Without it,
   deep-intronic variants get auto-Benign'd because the rule reads
   dbNSFP's coding c.

4. **For panel design**: review the 981 intergenic-flagged rows with the
   lab. They suggest the BED-anchored NM_ might be the wrong choice for
   some genes (FANCB, FANCD2 in particular — multiple isoforms with
   distinct exon structures).

5. **Inversion-breakpoint coverage in a short-read panel is a known
   limitation**. If the lab wants to detect SVs at MSH2 et al., add
   long-read confirmation or MLPA as a downstream validation step.
