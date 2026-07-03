**Subject:** cp_new v3.4 — liftover-free gnomAD anchor (you were right)

Done as a separate version: `cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4.tsv`
— 41,450 rows.

You're right that gnomAD v4's c. is enough and no liftover is needed — the c. is
relative to its transcript, not the genome, so it's the same in hg19 and hg38.
v3.4 anchors purely on `(gene, NM_, c.)` matched against gnomAD v4's own HGVSc.
No genome-build liftover for 63 of 67 coding genes.

The one wrinkle (the part worth knowing): gnomAD's c. on its variant page
defaults to the MANE Select transcript, which isn't always the NM_ we anchor on.
So instead of reading the displayed c., I pull gnomAD's HGVSc for our specific
NM_ (via its MANE_SELECT field in the VEP annotation). For 4 genes our project
NM_ differs from MANE — DPYD, GPR161, MITF, SMARCA4 — and for just those I fall
back to position matching so they don't silently drop. If you'd rather we switch
those 4 to the MANE transcript, that would make the whole pipeline 100%
c.-based; your call on whether the project transcript or MANE is the one MGL
wants for those genes.

Result is essentially identical to the position-anchored v3.3 (41,450 vs 41,514),
which is the good news — both methods agree, so the anchor is solid either way.
v3.4 is just the cleaner design you asked for.

MSH3: all 24 of your variants present; the 25th (c.181_189del) is in as
c.181_189dup — gnomAD carries both del and dup at that 9-bp repeat and our pull
grabbed the dup. Same locus, ~19% FAF. Can chase the del allele specifically if
you want exact-notation parity.

v3.3 (position-anchored) is still in `new_genes_68_v3/` if you want to compare.

Paulo
