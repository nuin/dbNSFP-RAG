**Subject:** Re: MSH3 review — gap is real, fixing

You're right. The bigger issue is that I had ~1,000 MSH3 high-FAF variants
sitting in a "needs HGVSc resolution" parking lot that never made it into
the upload. Across all 68 genes that's 17,244 high-frequency variants
(many at 70–99% gnomAD frequency, including indels) that should be Benign
FAF >5% but were silently excluded because I couldn't resolve their HGVSc.

That's the root cause for what you and Celia are seeing:
- High-freq dels/dups: parked as "noncoding indel, no HGVSc"
- High-freq intronic: parked as "noncoding, no HGVSc"
- Many syn variants: included via the local enumeration catalog, but indels
  aren't enumerated there either

I'm running MyVariant.info batch resolution on all 17,244 right now. Each
variant gets its proper HGVSc anchored to the project's RefSeq NM_, then I
add them to the upload as Benign FAF >5%. Should take ~15 minutes.

**Re your suggestion to do the high-frequency variants separately** — yes,
that's exactly the right model. The independent gnomAD pull at FAF≥5% is
ALREADY done (17,244 rows). The only missing step was HGVS resolution.
After this run, those become a clean drop-in source for the upload.

**Re missense calls being missed** — separate issue. The project's W3 rule
requires FAF >0.1% AND REVEL <0.29 AND SpliceAI ≤0.1. ClinVar B/LB missense
calls that don't meet all three thresholds (most rare ones don't) end up as
`rule_fails:W3_*` rather than Likely_benign. If you want to relax W3 (e.g.
allow REVEL-only or computational-only LB without the FAF gate), say the
word and I can produce a v3.2.

**Re synonymous catalog** — it's a local enumeration I built (script
`scripts/enumerate_synonymous_catalog.py`) that walks every CDS exon in
hg19 RefSeq, generates each codon, and lists every SNV that would produce
a synonymous substitution per the standard codon table (BioPython). It
was built specifically because dbNSFP excludes synonymous variants and
ClinVar only has the syn variants someone has submitted. Re-verified
yesterday: all 276,885 entries translate to genuine syn (BioPython
re-translation matches). It does NOT include indels — for the 1-2 nt
syn indels you'd expect, they're in gnomAD but not in this catalog.

**Re MSH3/DHFR overlap** — the multi-gene `genename` filter bug WAS real
but I fixed it back in May. dbNSFP rows with `genename = "MSH3;DHFR"` are
now correctly matched on MSH3. So the overlap shouldn't be losing rows
from the MSH3 column. But MyVariant might return DHFR annotations first
for variants in the overlap region — I'm anchoring resolution to MSH3's
NM_002439.5 specifically to avoid that.

Will report back with v3.2 numbers shortly.

Paulo
