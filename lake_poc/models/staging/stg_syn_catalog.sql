select
    gene, transcript,
    chr_grch37 as chr_hg19, cast(pos_grch37 as bigint) as pos_hg19,
    ref, alt, hgvs_c
from {{ source('raw','raw_syn_catalog') }}
