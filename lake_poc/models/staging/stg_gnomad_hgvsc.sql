select
    gene,
    split_part(nm, '.', 1) as nm_base,
    cdot,
    gnomad_id, g_hg38, p_hgvs
from {{ source('raw','raw_gnomad_hgvsc') }}
