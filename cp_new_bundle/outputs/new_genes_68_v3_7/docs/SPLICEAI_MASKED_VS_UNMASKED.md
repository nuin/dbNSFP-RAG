# SpliceAI: masked vs unmasked

The upload carries two SpliceAI columns — `SpliceAI_masked` and
`SpliceAI_unmasked`. Same model, same variant; the mask only changes which
signals are zeroed.

## The two modes

- **Unmasked (`-M 0`)** — the raw model output: probability of *any* splice
  change (gains and losses) at nearby positions. **This is what spliceai.org
  shows by default.**
- **Masked (`-M 1`)** — the same scores, but signals that merely *reinforce the
  already-annotated splice site* are set to 0. Keeps only changes that go
  *against* the reference gene model.

## The one rule

**unmasked ≥ masked, always.** Masking can only zero a score out, never raise it.

## Two practical cases

| Situation | Behaviour | Example |
|---|---|---|
| Deep in an exon (typical synonymous) | nothing to mask → **masked == unmasked** | HOXB13 c.600A>G = 0.20 / 0.20 |
| At / near a canonical splice site | masked can be **0** while unmasked is **high** (the reinforcing signal is masked out) | ANKRD26 c.1816A>C = masked 0.07 / unmasked 0.98 |

## Which we use for classification

The **unmasked** value is the honest "is there any splice signal" number and is
what matches spliceai.org, so the W2/W3 rules gate on **unmasked ≤ 0.1**. Both
columns are kept in the file so you can see them side by side.

## Quick sanity check on any row

- `masked == unmasked` → the variant is **not** near a canonical splice site.
- `unmasked` much higher than `masked` → the variant sits **on / near a canonical
  splice junction** — worth a closer look.

## How these scores were generated / validated

- SpliceAI v1.3.1, GRCh38, distance 50, run locally in both modes.
- Exported values are an exact copy of the raw SpliceAI VCF output
  (1000/1000 internal-integrity check).
- Cross-checked against the precomputed SpliceAI scores the spliceai.org site
  itself serves (public `storage.googleapis.com/tgg-viewer/...`); values agree
  where cleanly comparable.
- Note: differences vs **Alamut** specifically are a tool-configuration
  difference (its distance window / transcript set), not an error in these
  scores.
