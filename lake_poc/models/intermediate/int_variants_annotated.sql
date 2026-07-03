-- One row per candidate variant (missense/nonsense from dbNSFP + enumerated
-- synonymous), annotated with the gnomAD c.-anchor (g./p./gnomad_id + observed
-- flag) and FAF. This is the join layer that replaces the pile of Python scripts.
with coding as (
    select
        gene, hgvs_c, hgvs_p, consequence,
        chr_hg19, pos_hg19, ref, alt,
        revel, cadd, am, phastcons, phylop, faf95,
        'dbnsfp' as source
    from {{ ref('stg_dbnsfp') }}
),
syn as (
    select
        gene, hgvs_c, cast(null as varchar) as hgvs_p, 'synonymous' as consequence,
        chr_hg19, pos_hg19, ref, alt,
        cast(null as double) as revel, cast(null as double) as cadd,
        cast(null as varchar) as am,
        cast(null as double) as phastcons, cast(null as double) as phylop,
        cast(null as double) as faf95,
        'synonymous_catalog' as source
    from {{ ref('stg_syn_catalog') }}
),
unioned as (
    select * from coding
    union all
    select * from syn
),
-- gnomAD c.-anchor: match (gene, c.) -> gnomad_id/g./p.  (transcript-agnostic
-- within a gene; the c. is unique to its variant)
anchored as (
    select
        u.*,
        h.gnomad_id, h.g_hg38, h.p_hgvs,
        (h.gnomad_id is not null) as in_gnomad
    from unioned u
    left join (
        select gene, cdot, any_value(gnomad_id) as gnomad_id,
               any_value(g_hg38) as g_hg38, any_value(p_hgvs) as p_hgvs
        from {{ ref('stg_gnomad_hgvsc') }}
        group by gene, cdot
    ) h
      on u.gene = h.gene and u.hgvs_c = h.cdot
),
-- Conservation backfill: prefer dbNSFP's value, fall back to the UCSC track
-- (dbNSFP is sparse for synonymous positions -- this join is what makes the
-- synonymous rows classifiable, the crux of the manual rework this week).
with_cons as (
    select
        a.* exclude (phastcons, phylop),
        coalesce(a.phastcons, c.phastcons) as phastcons,
        coalesce(a.phylop, c.phylop)       as phylop
    from anchored a
    left join {{ ref('stg_conservation') }} c
      on a.chr_hg19 = c.chr_hg19 and a.pos_hg19 = c.pos_hg19
)
select * from with_cons
