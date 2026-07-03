-- The SeqNext upload: rule-based W1/W2/W3 classification, gnomAD-anchored.
-- Replaces build_upload_v3*.py + all the backfill/filter scripts.
with base as (
    select * from {{ ref('int_variants_annotated') }}
    where in_gnomad                      -- gnomAD anchor (c.-based, liftover-free)
      and consequence <> 'nonsense'      -- never auto-benign
),
classified as (
    select *,
        case
            -- W1: common
            when faf95 is not null and faf95 > 0.05
                then 'Benign FAF >5%'
            -- W3: rare missense, low REVEL, low splice
            when consequence = 'missense'
                 and faf95 is not null and faf95 > 0.001
                 and revel is not null and revel < 0.29
                then 'Likely_benign FAF >0.1%, REVEL<0.29'
            -- W2: synonymous, not highly conserved
            when consequence = 'synonymous'
                 and phastcons is not null and phastcons < 1.0
                then 'Benign synonymous, PhastCons<1.0'
            when consequence = 'synonymous'
                 and phastcons is not null and phastcons >= 1.0
                then 'Synonymous conserved (review, not auto-benign)'
            else 'rule_fails'
        end as classification
    from base
)
select
    gene, hgvs_c as cdot, hgvs_p, classification,
    consequence, source,
    gnomad_id, g_hg38,
    chr_hg19, pos_hg19, ref, alt,
    revel, cadd, am, faf95, phastcons, phylop
from classified
where classification <> 'rule_fails'
