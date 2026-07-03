**Subject:** cp_new updates — 4 fixes done + MBD4 & AXIN2 for review

Thanks for the thorough MSH3 review. All four action items are done; MBD4 and
AXIN2 attached for the next round.

**1. Synonymous "PhastCons unmeasured" — resolved (0 remain).**
Rather than removing them, I pulled the actual PhastCons + PhyloP from the UCSC
hg19 100-way tracks (same track dbNSFP bundles, just queried directly), so every
synonymous row now has conservation. No "unmeasured" entries left.
- Synonymous with PhastCons < 1.0 → Benign (BP7 supportable)
- Synonymous with PhastCons ≥ 1.0 → relabeled "Synonymous ... conserved
  (review, not auto-benign)" and KEPT in the list so you can filter/interpret
  them (per Paulo's call). Filter on column D or the PhastCons column.

**2. Verified missense are filtered on SpliceAI — and closed a real gap.**
Yes, splicing was used: W3 requires SpliceAI ≤0.1. But you were right that many
missense had a blank SpliceAI column — those had passed on an *assumed* 0 because
no measured value existed. I ran SpliceAI (masked + unmasked) on all 230 of them.
Now populated, and **7-8 missense that actually have SpliceAI >0.1 were dropped**
(they'd been passing by assumption). Every missense in the list now has real
masked + unmasked SpliceAI values.

**3. 5'UTR variants removed.**
Removed all c.-N (5'UTR) variants (16 from each ROI list). Note: 3'UTR (c.*N)
variants are still in — you only mentioned 5'UTR. Say the word if you want those
gone too (~170 rows).

**4. MBD4 and AXIN2 for review** — attached (v3.4, c.-based gnomAD anchor,
ROI [-20,+10]):

| gene | rows | Benign W2 | Synonymous conserved (review) | Benign W1 | LB W3 (missense) |
|---|---:|---:|---:|---:|---:|
| MBD4  | 283 | 242 | 33  | 4 | 4 |
| AXIN2 | 659 | 403 | 233 | 7 | 16 |

Both have 100% PhastCons, 100% gnomad_id, and every missense carries REVEL +
SpliceAI (masked & unmasked) + p.HGVS.

Reminder on the columns for your lookups:
- gnomad_id (e.g. 5-80768028-G-A) pastes straight into gnomAD v4.1
- g_hg19 / g_hg38 both present (we're an hg19 shop; hg38 is liftover for gnomAD)
- p_hgvs for the protein change

Files: `REVIEW_MBD4_AXIN2/MBD4_v3_4_roi_review.tsv`,
`REVIEW_MBD4_AXIN2/AXIN2_v3_4_roi_review.tsv`

Paulo
