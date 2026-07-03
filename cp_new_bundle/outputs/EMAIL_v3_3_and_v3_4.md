**Subject:** cp_new upload — two versions (v3.3 and v3.4)

Sending two versions of the cp_new upload. Same rules (W1/W2/W3), same gnomAD
anchoring decision — two different ways of checking gnomAD membership:

- **v3.3 (41,514 rows)** — matches variants to gnomAD by **genomic position**
  (hg38 coordinates; liftover used internally for the check only, never for any
  coordinate that reaches the upload).
- **v3.4 (41,450 rows)** — matches by **HGVS c. on the project transcript**
  against gnomAD v4's own HGVSc. **No liftover** for 63 of 67 genes (the c. is
  build-independent, as you pointed out). The 4 genes where our NM_ differs from
  gnomAD's MANE Select — DPYD, GPR161, MITF, SMARCA4 — fall back to
  position-matching so they don't drop.

The two agree to within 64 rows (0.15%), which is the main point: the
gnomAD-observed filtering is robust either way. v3.4 is the cleaner, liftover-free
design.

Both files:
- column D is the classification SeqNext imports
- every variant is observed in gnomAD v4.1
- classified purely by W1/W2/W3 (no ClinVar evidence used)
- 67/68 genes (TERC is a non-coding RNA — no W1/W2/W3 applies)

Files:
- v3.3: `new_genes_68_v3/UPLOAD_v3.tsv`
- v3.4: `new_genes_68_v3_4/UPLOAD_v3_4.tsv`

One open question for the 4 mismatch genes: do you want the project NM_ or MANE
Select? Choosing MANE would make v3.4 100% c.-based (no liftover anywhere).

Paulo
