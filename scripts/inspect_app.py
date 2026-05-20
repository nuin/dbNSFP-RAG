#!/usr/bin/env python3
"""Streamlit web UI for the cp_new variant inspector.

Wraps scripts/inspect_variant.py's fetchers and renders the report in a
browser-friendly layout: sidebar for variant input + gene browser, main
panel with collapsible sections per data source.

Run:
  uv pip install streamlit
  uv run streamlit run scripts/inspect_app.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent))
from inspect_variant import (
    Variant, find_in_per_gene, find_in_classifications,
    fetch_myvariant, fetch_variantvalidator, fetch_ensembl_gene,
    FINAL, DATA,
)

st.set_page_config(page_title="cp_new variant inspector", page_icon="🧬", layout="wide")

st.title("cp_new variant inspector")
st.caption(
    "Alamut-style lookup across local cp_new cache + public APIs "
    "(myvariant.info, VariantValidator, Ensembl REST, ClinVar). "
    "Source files at `data/exports/cp_new/`."
)

# ------------------------------------------------------------------------
# Sidebar: variant input + gene browser
# ------------------------------------------------------------------------

with st.sidebar:
    st.subheader("Inspect a variant")
    mode = st.radio("Input form", ["gene + HGVS c.", "chr:pos:ref:alt", "rsID"], index=0)
    gene = transcript = hgvs = coord = rsid = ""
    build = "hg19"
    if mode == "gene + HGVS c.":
        gene = st.text_input("Gene symbol", "GOT2")
        hgvs = st.text_input("HGVS c. notation", "c.816C>T")
        transcript = st.text_input("Transcript (optional NM_)", "")
    elif mode == "chr:pos:ref:alt":
        gene = st.text_input("Gene (helps lookup)", "GOT2")
        coord = st.text_input("Coordinates", "16:58750604:G:A")
        build = st.selectbox("Build", ["hg19", "hg38"], index=0)
    else:
        rsid = st.text_input("dbSNP rsID", "rs1058192")

    no_external = st.checkbox("Local cache only (skip APIs)", value=False)
    submit = st.button("Inspect", type="primary")

    st.divider()
    st.subheader("Browse by gene")
    gene_browse = st.text_input("Gene symbol to browse", "")
    n_show = st.slider("Rows", 5, 100, 20)


# ------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------

def build_variant(gene, hgvs, coord, build, rsid, transcript) -> Variant:
    v = Variant()
    if gene: v.gene = gene.upper()
    if transcript: v.transcript = transcript
    if hgvs: v.hgvs_c = hgvs
    if coord:
        m = re.match(r"^(?:chr)?([\dXYM]+)[:-](\d+)[:-]([ACGTN-]+)[:-]([ACGTN-]+)$", coord, re.I)
        if m:
            chrom, pos, ref, alt = m.groups()
            if build == "hg38":
                v.chr_hg38, v.pos_hg38 = chrom, int(pos)
            else:
                v.chr_hg19, v.pos_hg19 = chrom, int(pos)
            v.ref, v.alt = ref.upper(), alt.upper()
    if rsid: v.rsid = rsid.lower()
    return v


@st.cache_data(show_spinner=False)
def _cached_myvariant(key: str, v_json: str):
    return fetch_myvariant(Variant(**json.loads(v_json)))


@st.cache_data(show_spinner=False)
def _load_final():
    if not FINAL.exists():
        return pd.DataFrame()
    return pd.read_csv(FINAL, sep="\t", dtype=str).fillna("")


# ------------------------------------------------------------------------
# Main: variant report
# ------------------------------------------------------------------------

if submit:
    v = build_variant(gene, hgvs, coord, build, rsid, transcript)

    if v.gene:
        st.markdown(f"### {v.gene} — {v.hgvs_c or v.key()}")
    else:
        st.markdown(f"### {v.key()}")

    col_id, col_class = st.columns([1, 1])

    with col_id:
        with st.expander("Identity", expanded=True):
            per_gene = find_in_per_gene(v) if v.gene else None
            if per_gene:
                st.write({
                    "gene": per_gene.get("genename","").split(";")[0],
                    "Ensembl tx": per_gene.get("Ensembl_transcriptid","").split(";")[0],
                    "HGVSc (dbNSFP)": per_gene.get("HGVSc_snpEff","").split(";")[0],
                    "HGVSp (dbNSFP)": per_gene.get("HGVSp_snpEff","").split(";")[0],
                    "hg19": f"chr{per_gene.get('hg19_chr')}:{per_gene.get('hg19_pos(1-based)')} {per_gene.get('ref')}>{per_gene.get('alt')}",
                    "hg38": f"chr{per_gene.get('#chr')}:{per_gene.get('pos(1-based)')} {per_gene.get('ref')}>{per_gene.get('alt')}",
                })
            else:
                st.info("Not in our cp_new dbNSFP cache (synonymous variants etc.).")

    with col_class:
        with st.expander("Our pipeline classification", expanded=True):
            cls = find_in_classifications(v)
            if cls:
                for c in cls:
                    st.success(f"**{c['classification']}**")
                    st.code(f"source={c.get('source')}  vv_status={c.get('vv_status')}\n"
                            f"transcript={c.get('transcript')}\n"
                            f"hgvs_c={c.get('hgvs_c')}", language="text")
            else:
                st.warning("No matching classification in cp_new_seqnext_FINAL.tsv.")

    with st.expander("Frequencies (gnomAD)"):
        if per_gene:
            st.dataframe(pd.DataFrame([{
                "gnomAD v4.1 joint AF":     per_gene.get("gnomad_v41_af_joint"),
                "v4.1 grpmax FAF95":         per_gene.get("gnomad_v41_faf95_grpmax"),
                "v4.1 grpmax ancestry":      per_gene.get("gnomad_v41_grpmax_anc"),
                "v4.1 nhomalt":              per_gene.get("gnomad_v41_nhomalt"),
                "v4.1 FILTER":               per_gene.get("gnomad_v41_filter"),
                "v2.1.1 controls AF":        per_gene.get("gnomAD2.1.1_exomes_controls_AF"),
                "1000G AF":                  per_gene.get("1000Gp3_AF"),
            }]).T.rename(columns={0:"value"}))
        else:
            st.info("No dbNSFP frequencies for this variant in local cache.")

    with st.expander("Predictors"):
        if per_gene:
            preds = {k: per_gene.get(v) for k, v in [
                ("REVEL","REVEL_score"), ("CADD phred","CADD_phred"),
                ("AlphaMissense pred","AlphaMissense_pred"),
                ("AlphaMissense score","AlphaMissense_score"),
                ("SIFT","SIFT_pred"), ("Polyphen2 HDIV","Polyphen2_HDIV_pred"),
                ("MutationTaster","MutationTaster_pred"), ("BayesDel","BayesDel_addAF_pred"),
                ("PROVEAN","PROVEAN_pred"), ("ClinPred","ClinPred_pred"),
                ("MetaSVM","MetaSVM_pred"),
            ]}
            st.dataframe(pd.DataFrame([preds]).T.rename(columns={0:"value"}))
        else:
            st.info("No dbNSFP predictors.")

    with st.expander("Splice & conservation"):
        if per_gene:
            st.dataframe(pd.DataFrame([{
                "SpliceAI masked DS_MAX":    per_gene.get("spliceai_ds_max_masked"),
                "PhastCons100way":           per_gene.get("phastCons100way_vertebrate"),
                "PhyloP100way":              per_gene.get("phyloP100way_vertebrate"),
                "PhyloP470way mammalian":    per_gene.get("phyloP470way_mammalian"),
                "GERP++ RS":                 per_gene.get("GERP++_RS"),
                "Interpro domain":           per_gene.get("Interpro_domain"),
            }]).T.rename(columns={0:"value"}))
        else:
            st.info("No conservation/splice scores in cache.")

    with st.expander("Clinical assertions (ClinVar)"):
        if per_gene:
            st.write({
                "ClinVar ID":     per_gene.get("clinvar_id"),
                "ClinVar sig":    per_gene.get("clinvar_clnsig"),
                "ClinVar review": per_gene.get("clinvar_review"),
                "ClinVar trait":  per_gene.get("clinvar_trait"),
            })
        if not no_external:
            mv = _cached_myvariant(v.key(), json.dumps({k: getattr(v, k) for k in
                ("gene","transcript","hgvs_c","chr_hg19","pos_hg19","chr_hg38","pos_hg38",
                 "ref","alt","rsid","build_input")}))
            if mv and "clinvar" in mv:
                st.markdown("**myvariant.info ClinVar payload:**")
                st.json(mv["clinvar"], expanded=False)

    with st.expander("Cross-references (myvariant.info aggregator)"):
        if no_external:
            st.info("External APIs disabled.")
        else:
            mv = _cached_myvariant(v.key(), json.dumps({k: getattr(v, k) for k in
                ("gene","transcript","hgvs_c","chr_hg19","pos_hg19","chr_hg38","pos_hg38",
                 "ref","alt","rsid","build_input")}))
            if mv:
                refs = {}
                if "dbsnp" in mv:
                    refs["dbSNP"] = mv["dbsnp"].get("rsid")
                if "cosmic" in mv:
                    refs["COSMIC"] = mv["cosmic"].get("cosmic_id") if isinstance(mv["cosmic"], dict) else "present"
                if "evs" in mv:
                    refs["EVS"] = "present (older Exome Variant Server)"
                st.write(refs)
                with st.expander("Raw myvariant.info JSON"):
                    st.json(mv, expanded=False)
            else:
                st.warning("No myvariant.info entry.")


# ------------------------------------------------------------------------
# Gene browser
# ------------------------------------------------------------------------

if gene_browse:
    st.subheader(f"All classified variants for {gene_browse.upper()}")
    final = _load_final()
    if final.empty:
        st.error(f"Missing {FINAL}")
    else:
        sub = final[final["gene"] == gene_browse.upper()]
        if sub.empty:
            st.warning(f"No classifications for {gene_browse}.")
        else:
            st.write(f"{len(sub):,} rows  (showing first {n_show})")
            display_cols = [c for c in ["gene","transcript","hgvs_c","classification",
                            "FAF95_grpmax","REVEL","SpliceAI_masked","PhastCons100way",
                            "ClinVar_sig","source"] if c in sub.columns]
            st.dataframe(sub[display_cols].head(n_show), height=600)
            st.download_button(
                "Download full TSV",
                sub.to_csv(sep="\t", index=False).encode(),
                file_name=f"{gene_browse}_classified.tsv",
            )
