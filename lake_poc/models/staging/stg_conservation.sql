select chr_hg19, cast(pos_hg19 as bigint) as pos_hg19,
       try_cast(phastcons as double) as phastcons,
       try_cast(phylop as double) as phylop
from {{ source('raw','raw_conservation') }}
