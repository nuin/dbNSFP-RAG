-- Non-synonymous coding variants from dbNSFP (missense + nonsense).
select
    gene,
    chr_hg38, cast(pos_hg38 as bigint) as pos_hg38,
    chr_hg19, cast(pos_hg19 as bigint) as pos_hg19,
    ref, alt, aaref, aaalt,
    hgvs_c, hgvs_p,
    try_cast(revel as double)     as revel,
    try_cast(cadd  as double)     as cadd,
    am,
    try_cast(faf95 as double)     as faf95,
    try_cast(phastcons as double) as phastcons,
    try_cast(phylop as double)    as phylop,
    case
        when aaref = aaalt then 'synonymous'
        when aaalt in ('X','*') then 'nonsense'
        else 'missense'
    end as consequence
from {{ source('raw','raw_dbnsfp') }}
where hgvs_c <> ''
