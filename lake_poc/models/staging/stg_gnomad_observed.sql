select
    gene, chr_hg38, cast(pos_hg38 as bigint) as pos_hg38,
    ref, alt,
    try_cast(faf95 as double) as faf95
from {{ source('raw','raw_gnomad_observed') }}
