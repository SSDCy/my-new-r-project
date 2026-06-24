# app.py
# Streamlit Proteomics Platform
# Fixed: ESM2 batch buttons feedback, STRING UniProt fallback always shows mapping.

import streamlit as st
import pandas as pd
import numpy as np
import warnings
import plotly.express as px
import tempfile
import os
import requests
import gzip
import shutil
import re
import json
import time

# ---------- Utility modules ----------
from utils.data_io import load_expression_data, load_sample_info
from utils.preprocessing import apply_missing_filter, impute_missing, run_combat_py
from utils.normalization import total_intensity_normalize
from utils.de_analysis import run_de_analysis
from utils.visualization import (
    plot_volcano, plot_heatmap, plot_venn, plot_upset,
    plot_volcano_single_annotated, build_combined_plot
)
from utils.cd_search import (
    batch_cd_search_simple, batch_cd_search_with_cache,
    parse_cd_tsv, load_cd_cache, get_cached_tasks, CD_CHUNK_SIZE
)
from utils.export_utils import create_excel_report
from utils.data_quality import (
    plot_missing_heatmap, plot_valid_values_per_sample,
    plot_venn_by_group, get_peptide_sequences_table,
    plot_peptide_length_histogram, plot_sample_correlation_heatmap,
    plot_pca_raw_by_group
)
from utils.eggnog_annotation import (
    run_eggnog_annotation_local,
    parse_eggnog_manual_file,
    EMAPPER_PATH,
    load_go_names,
    load_ko_names,
    get_go_name,
    get_ko_name
)
from utils.esm_search import (
    load_esm_model, build_reference_library, get_esm_embedding,
    find_similar_proteins, get_sequence_for_id, is_pretrained_available,
    ESM_AVAILABLE, batch_search_top1, fetch_uniprot_annotations_batch,
    parse_fasta, compute_fasta_hash, is_uniprot_accession, extract_accession
)
from utils.enrichment import get_go_background, enrich_go
from utils.llm_summary import generate_function_summary

from utils.string_ppi import (
    extract_uniprot_from_eggnog,
    run_blast_mapping,
    call_string_api,
    build_ppi_network_html
)

warnings.filterwarnings("ignore", category=RuntimeWarning)

st.set_page_config(page_title="Proteomics Platform", layout="wide")

# ---------- Debug helper ----------
def debug(msg: str):
    timestamp = time.strftime("%H:%M:%S")
    print(f"[{timestamp}][DEBUG] {msg}", flush=True)

# ---------- Helper functions ----------
def parse_manual_uniprot_ids(text: str) -> list:
    ids = [x.strip() for x in text.replace(",", "\n").replace(";", "\n").splitlines() if x.strip()]
    if any(x.startswith(">") for x in ids):
        st.error("FASTA sequences are not allowed here.")
        st.stop()
    if not ids:
        st.warning("No identifiers entered.")
        st.stop()
    return ids

def build_annotations_from_eggnog(eggnog_df):
    if eggnog_df is None or eggnog_df.empty:
        return {}
    id_col = eggnog_df.columns[0]
    annotations = {}
    if 'GOs' not in eggnog_df.columns:
        return {}
    for _, row in eggnog_df.iterrows():
        pid = str(row[id_col])
        gos = str(row.get('GOs', ''))
        go_list = [g.strip() for g in gos.split(',') if g.strip() and g.strip() != '-']
        if go_list:
            annotations[pid] = {'go': go_list}
    return annotations

def enrich_esm_result_with_eggnog(result_df, eggnog_df):
    debug("Enriching ESM2 results with eggNOG annotations...")
    if eggnog_df is None or eggnog_df.empty:
        return result_df
    id_col = eggnog_df.columns[0]
    lookup = {}
    for _, row in eggnog_df.iterrows():
        pid = str(row[id_col])
        gene = str(row.get('Preferred_name', '') or row.get('Description', ''))
        protein_name = str(row.get('Description', ''))
        lookup[pid] = {
            'Gene': gene,
            'Protein_Name': protein_name,
            'GO': str(row.get('GOs', '')),
            'EC': str(row.get('EC', ''))
        }
    for idx, protein_id in enumerate(result_df['Protein_ID']):
        if protein_id in lookup:
            result_df.at[idx, 'Gene'] = lookup[protein_id]['Gene']
            result_df.at[idx, 'Protein_Name'] = lookup[protein_id]['Protein_Name']
            result_df.at[idx, 'GO'] = lookup[protein_id]['GO']
            result_df.at[idx, 'EC'] = lookup[protein_id]['EC']
        else:
            result_df.at[idx, 'Gene'] = ''
            result_df.at[idx, 'Protein_Name'] = ''
            result_df.at[idx, 'GO'] = ''
            result_df.at[idx, 'EC'] = ''
    return result_df

# ---------- Batch search step functions ----------
def init_batch_search_state(session_state, total_sequences: int):
    session_state.batch_search_in_progress = True
    session_state.batch_search_index = 0
    session_state.batch_search_results = []
    session_state.batch_search_total = total_sequences

def run_batch_search_step(query_seqs, model, batch_converter, custom_library, session_state, step_size=100):
    total = session_state.batch_search_total
    start = session_state.batch_search_index
    end = min(start + step_size, total)
    ids = list(custom_library.keys())
    embeddings = np.stack([custom_library[k] for k in ids], axis=0)
    norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
    norms[norms == 0] = 1
    embeddings_norm = embeddings / norms
    for i in range(start, end):
        seq = query_seqs[i]
        try:
            q_emb = get_esm_embedding(seq, model, batch_converter)
        except Exception:
            session_state.batch_search_results.append({'query_index': i, 'best_id': None, 'similarity': 0.0})
            continue
        q_norm = np.linalg.norm(q_emb)
        if q_norm == 0:
            session_state.batch_search_results.append({'query_index': i, 'best_id': None, 'similarity': 0.0})
            continue
        sims = np.dot(embeddings_norm, q_emb / q_norm)
        best_idx = np.argmax(sims)
        session_state.batch_search_results.append({
            'query_index': i, 'best_id': ids[best_idx], 'similarity': float(sims[best_idx])
        })
    session_state.batch_search_index = end
    if end >= total:
        session_state.batch_search_in_progress = False
        return True
    return False

def display_batch_results_and_enrichment(results, target_ids, annotations):
    if results is not None and len(results) > 0:
        st.subheader("Top Match Proteins")
        df = pd.DataFrame({
            "Query No.": [r['query_index']+1 for r in results],
            "Best Match": [r['best_id'] if r['best_id'] else 'Failed' for r in results],
            "Similarity": [round(r['similarity'],4) for r in results]
        })
        st.dataframe(df, width='stretch')
        if all(round(r['similarity'],4) >= 0.9999 for r in results):
            st.warning("All similarities are ~1 because the query sequences are identical to those in the reference library. For GO enrichment, this is expected.")

    if annotations and len(target_ids) > 0:
        valid_targets = [pid for pid in target_ids if pid in annotations]
        if not valid_targets:
            st.warning("None of the selected proteins have GO annotations.")
            return
        if len(valid_targets) < len(target_ids):
            st.info(f"Only {len(valid_targets)} out of {len(target_ids)} proteins have GO annotations.")
        st.subheader("GO Enrichment Analysis")
        bg_list = list(annotations.keys())
        st.info(f"Background proteins ({len(bg_list)}): {', '.join(bg_list[:20])}{'...' if len(bg_list)>20 else ''}")
        if len(valid_targets) == len(annotations):
            st.warning("Target set is identical to the background set. Please select a subset.")
            return
        go_bg = get_go_background(annotations)
        enrich_df = enrich_go(valid_targets, annotations, go_bg)
        if len(enrich_df) > 0:
            fig = px.scatter(enrich_df, x='enrichment_ratio', y='GO', size='count',
                             color='adjusted_p', color_continuous_scale='Reds_r',
                             title="Enriched GO Terms")
            st.plotly_chart(fig, width='stretch')
            st.dataframe(enrich_df)
        else:
            st.info("No significantly enriched GO terms (p<0.05).")
    else:
        st.warning("No annotations available for enrichment.")
        if st.session_state.get("eggnog_result") is not None:
            if st.button("Use eggNOG annotations for enrichment", key="batch_use_eggnog"):
                annots = build_annotations_from_eggnog(st.session_state.eggnog_result)
                st.session_state.custom_annotations = annots
                st.rerun()
        else:
            if not st.session_state.get("use_pretrained", True) and "custom_library" in st.session_state:
                if st.button("Load UniProt annotations now", key="batch_load_annot"):
                    lib_ids = list(st.session_state.custom_library.keys())
                    annots = fetch_uniprot_annotations_batch(lib_ids)
                    st.session_state.custom_annotations = annots
                    st.rerun()

# ---------- Initialize session state ----------
defaults = {
    'raw_expr_df': None, 'expr_df': None, 'lfq_cols': None, 'sample_names': None,
    'sample_info': None, 'processed': None, 'norm_data': None, 'groups': {},
    'comparisons': [], 'de_results': {}, 'page': "Home", 'cleaning_stats': None,
    'cd_result': None, 'eggnog_result': None, 'integrated_data': None,
    'integrated_protein_list': None, 'integrated_go_data': None,
    'integrated_enrich': None, 'integrated_use_api': False,
    'ref_fasta_text': "", 'eggnog_fasta_text': "",
    'batch_final_results': None, 'batch_target_ids': None,
    'last_search_result': None, 'current_summary': None,
    'custom_annotations': None, 'use_pretrained': False,
    'custom_ids': None,
    'query_id_f1': "", 'query_seq_f1': "",
    'batch_text': "",
    'esm_model': None, 'esm_alphabet': None, 'esm_batch_converter': None
}
for key, val in defaults.items():
    if key not in st.session_state:
        st.session_state[key] = val

# ---------- Sidebar navigation ----------
st.sidebar.title("Navigation")
pages = [
    "Home", "Data Upload", "Data Quality", "Preprocessing",
    "Define Groups & Comparisons", "Differential Analysis",
    "Visualization", "CD-Search", "Export", "EggNOG Annotation",
    "ESM2 Similarity Search", "STRING PPI Network",
    "Integrated Analysis"
]
current_page = pages.index(st.session_state.page) if st.session_state.page in pages else 0
page = st.sidebar.radio("Go to", pages, index=current_page)
st.session_state.page = page

# ---------- Page: Home ----------
if page == "Home":
    st.title("Proteomics Differential Analysis Platform")
    st.write("From MaxQuant data to differential expression, annotation, and report.")
    if st.button("Start Analysis -> Upload Data"):
        st.session_state.page = "Data Upload"
        st.rerun()

# ---------- Page: Data Upload ----------
elif page == "Data Upload":
    st.header("Data Upload")
    st.info("Data cleaning (removing reverse hits, contaminants) will be performed automatically during Preprocessing.")
    expr_file = st.file_uploader("Upload proteinGroups.txt", type=['txt'])
    sample_file = st.file_uploader("Upload sample info (CSV or Excel)", type=['csv', 'xlsx'])
    if expr_file and sample_file:
        df_clean, lfq_cols, sample_names, col_type, clean_stats, df_raw = load_expression_data(expr_file, clean=False)
        si = load_sample_info(sample_file)
        if df_clean is not None and si is not None:
            st.session_state.raw_expr_df = df_raw
            st.session_state.expr_df = df_clean
            st.session_state.lfq_cols = lfq_cols
            st.session_state.sample_names = sample_names
            st.session_state.sample_info = si
            st.session_state.cleaning_stats = None
            st.success("Data loaded successfully! (Raw data, cleaning will be done in Preprocessing)")
            st.write(f"**Data dimensions**: {df_clean.shape[0]} rows x {df_clean.shape[1]} cols")
            st.write("First 5 rows:", df_clean.head())
            st.write("Sample info:", si.head())

# ---------- Page: Data Quality ----------
elif page == "Data Quality":
    st.header("Data Quality")
    if st.session_state.expr_df is None:
        st.warning("Please upload data first.")
    else:
        expr_df = st.session_state.expr_df
        lfq_cols = st.session_state.lfq_cols
        sample_info = st.session_state.sample_info
        fig_miss = plot_missing_heatmap(expr_df, lfq_cols)
        if fig_miss: st.pyplot(fig_miss)
        fig_valid = plot_valid_values_per_sample(expr_df, lfq_cols)
        if fig_valid: st.pyplot(fig_valid)
        venn_figs = plot_venn_by_group(expr_df, lfq_cols)
        if venn_figs:
            cols = st.columns(3)
            for i, (gname, fig) in enumerate(venn_figs):
                with cols[i % 3]: st.pyplot(fig)
        else:
            st.info("No Venn diagrams available (need 2-3 samples per group)")
        merged, _ = get_peptide_sequences_table(expr_df)
        if merged is not None:
            fig_len = plot_peptide_length_histogram(merged)
            if fig_len: st.pyplot(fig_len)
        else:
            st.info("No peptide sequence data.")
        fig_cor = plot_sample_correlation_heatmap(expr_df, lfq_cols)
        if fig_cor: st.pyplot(fig_cor)
        fig_pca = plot_pca_raw_by_group(expr_df, lfq_cols, sample_info)
        if fig_pca: st.pyplot(fig_pca)

# ---------- Page: Preprocessing ----------
elif page == "Preprocessing":
    st.header("Preprocessing")
    if st.session_state.expr_df is None:
        st.warning("Please upload data first.")
    else:
        max_missing = st.slider("Max missing fraction", 0.0, 1.0, 0.5)
        impute_method = st.selectbox("Imputation method", ["none","knn","ppca","quantile","minvalue"], index=3)
        quantile_prob = 0.01
        if impute_method == "quantile":
            quantile_prob = st.number_input("Quantile", 0.001, 0.5, 0.01, 0.01)
        lfq_cols = st.session_state.lfq_cols
        baseline_sample = st.selectbox("Baseline sample for normalization", lfq_cols,
                                        index=lfq_cols.index("LFQ intensity WT-1") if "LFQ intensity WT-1" in lfq_cols else 0)
        if st.button("Run Preprocessing"):
            start_time = time.time()
            df = st.session_state.expr_df.copy()

            if st.session_state.cleaning_stats is None:
                debug("Performing data cleaning...")
                reverse_removed, contam_removed, con_removed = 0, 0, 0
                if 'Reverse' in df.columns:
                    mask = df['Reverse'].str.contains(r'\+', na=False)
                    reverse_removed = mask.sum()
                    df = df[~mask]
                if 'Potential contaminant' in df.columns:
                    mask = df['Potential contaminant'].str.contains(r'\+', na=False)
                    contam_removed = mask.sum()
                    df = df[~mask]
                if 'Protein IDs' in df.columns:
                    mask = df['Protein IDs'].str.startswith('CON_', na=False)
                    con_removed = mask.sum()
                    df = df[~mask]
                    if 'Master protein IDs' not in df.columns:
                        df['Master protein IDs'] = df['Protein IDs'].str.split(';').str[0]
                st.session_state.cleaning_stats = {
                    'original': st.session_state.raw_expr_df.shape[0],
                    'reverse_removed': reverse_removed,
                    'contaminant_removed': contam_removed,
                    'con_removed': con_removed,
                    'retained': df.shape[0]
                }
                st.session_state.expr_df = df
                debug(f"Data cleaned: {df.shape[0]} rows retained.")
            else:
                debug("Data already cleaned, skipping cleaning step.")

            if st.session_state.cleaning_stats:
                cs = st.session_state.cleaning_stats
                st.write(f"**Cleaning Summary**: Original {cs['original']}, Removed reverse {cs['reverse_removed']}, contaminants {cs['contaminant_removed']}, CON_ {cs['con_removed']}. Retained {cs['retained']}")

            df = apply_missing_filter(df, lfq_cols, max_missing, st.session_state.sample_info)
            df = impute_missing(df, lfq_cols, impute_method, quantile_prob=quantile_prob, k=10)
            df = total_intensity_normalize(df, lfq_cols, baseline_sample, st.session_state.sample_names)
            st.session_state.processed = df
            elapsed = time.time() - start_time
            st.success(f"Preprocessing completed in {elapsed:.1f} seconds.")
            st.dataframe(df.head())

# ---------- Page: Define Groups & Comparisons ----------
elif page == "Define Groups & Comparisons":
    st.header("Define Groups & Comparisons")
    if st.session_state.sample_info is not None and st.session_state.lfq_cols:
        si = st.session_state.sample_info
        lfq_cols = st.session_state.lfq_cols
        short_from_cols = [c[len('LFQ intensity '):] if c.startswith('LFQ intensity ') else c for c in lfq_cols]
        short_std = [s.replace('.','_').strip() for s in short_from_cols]
        col_std = [c.replace('.','_').strip() for c in lfq_cols]
        info_index_std = [str(idx).replace('.','_').strip() for idx in si.index]
        match_map = {}
        for i, (fc, sh, fstd) in enumerate(zip(lfq_cols, short_std, col_std)):
            if sh in info_index_std:
                match_map[fc] = sh
            elif fstd in info_index_std:
                match_map[fc] = fstd
            else:
                for idxn in info_index_std:
                    if sh in idxn or idxn in sh:
                        match_map[fc] = idxn
                        break
        groups = {}
        used_col = None
        if 'SubGroup' in si.columns:
            for fc, midx in match_map.items():
                sg = si.loc[si.index[info_index_std.index(midx)], 'SubGroup']
                groups.setdefault(sg, []).append(fc)
            used_col = 'SubGroup'
        elif 'Group' in si.columns:
            for fc, midx in match_map.items():
                g = si.loc[si.index[info_index_std.index(midx)], 'Group']
                groups.setdefault(g, []).append(fc)
            used_col = 'Group'
        else:
            prefixes = {}
            for s in short_std:
                parts = s.split('_')
                pref = '_'.join(parts[:-1]) if len(parts)>=2 else s
                prefixes.setdefault(pref, []).append(s)
            for pref, slist in prefixes.items():
                flist = [lfq_cols[short_std.index(s)] for s in slist if s in short_std]
                if flist:
                    groups[pref] = flist
            used_col = 'prefix'
        st.session_state.groups = groups
        st.write(f"### Groups (based on {used_col})")
        for gn, cols in groups.items():
            st.write(f"**{gn}** ({len(cols)} samples): {', '.join(cols)}")
        if groups:
            st.subheader("Add comparison")
            c1, c2 = st.columns(2)
            with c1: treat = st.selectbox("Treatment", list(groups.keys()), key="treat")
            with c2: ctrl = st.selectbox("Control", list(groups.keys()), key="ctrl")
            if st.button("Add comparison"):
                if treat != ctrl:
                    comp_name = f'{treat} vs {ctrl}'
                    if comp_name not in [c['name'] for c in st.session_state.comparisons]:
                        st.session_state.comparisons.append({'treat':treat, 'ctrl':ctrl, 'name':comp_name})
                        st.success(f"Comparison '{comp_name}' added.")
                else:
                    st.error("Treatment and control must be different.")
            st.subheader("Add all vs WT")
            wt_candidates = [g for g in groups if 'wt' in g.lower() or g.lower()=='control']
            if wt_candidates:
                wt_group = wt_candidates[0]
                others = [g for g in groups if g != wt_group]
                if others and st.button(f"Add all vs {wt_group}"):
                    added = 0
                    for tr in others:
                        cname = f'{tr} vs {wt_group}'
                        if cname not in [c['name'] for c in st.session_state.comparisons]:
                            st.session_state.comparisons.append({'treat':tr, 'ctrl':wt_group, 'name':cname})
                            added += 1
                    if added: st.success(f"{added} comparisons added.")
        if st.session_state.comparisons:
            st.subheader("Current comparisons")
            for comp in st.session_state.comparisons:
                st.write(f"- **{comp['name']}**: {comp['treat']} vs {comp['ctrl']}")

# ---------- Page: Differential Analysis ----------
elif page == "Differential Analysis":
    st.header("Differential Analysis")
    if not st.session_state.comparisons:
        st.warning("Please define comparisons first.")
    elif st.session_state.processed is None:
        st.warning("Please run preprocessing first.")
    else:
        method = st.selectbox("Statistical method", ["t-test","wilcoxon","limma"])
        fc_up = st.number_input("FC up >", value=1.2)
        fc_down = st.number_input("FC down <", value=0.84)
        p_cut = st.number_input("P-value threshold", value=0.05)
        col1,col2,col3,col4 = st.columns(4)
        with col1: min_treat_valid = st.number_input("Min valid replicates treatment", min_value=1, value=2)
        with col2: min_ctrl_valid = st.number_input("Min valid replicates control", min_value=1, value=2)
        with col3: min_rep_ttest = st.number_input("Min replicates for t-test", min_value=1, value=2)
        with col4: pass
        col5,col6 = st.columns(2)
        with col5: min_rep_inc = st.number_input("Min replicates for Increase", min_value=1, value=2)
        with col6: min_rep_dec = st.number_input("Min replicates for Decrease", min_value=1, value=2)
        min_unique_pep = st.number_input("Min Unique Peptides", min_value=1, value=2)
        use_pep_filter = st.checkbox("Enable Unique Peptide filter", value=True)
        if st.button("Run all comparisons"):
            start_time = time.time()
            norm_df = st.session_state.processed.copy()
            if use_pep_filter and 'Unique peptides' in norm_df.columns:
                norm_df['Unique peptides'] = pd.to_numeric(norm_df['Unique peptides'], errors='coerce')
                norm_df = norm_df.dropna(subset=['Unique peptides'])
                norm_df = norm_df[norm_df['Unique peptides'] >= min_unique_pep]
            results = {}
            for comp in st.session_state.comparisons:
                res = run_de_analysis(norm_df, st.session_state.groups[comp['treat']],
                                      st.session_state.groups[comp['ctrl']],
                                      fc_up=fc_up, fc_down=fc_down, p_cut=p_cut,
                                      method=method, min_treat_valid=min_treat_valid,
                                      min_ctrl_valid=min_ctrl_valid, min_rep_ttest=min_rep_ttest,
                                      min_rep_inc=min_rep_inc, min_rep_dec=min_rep_dec)
                if res is not None:
                    results[comp['name']] = res
            st.session_state.de_results = results
            elapsed = time.time() - start_time
            st.success(f"Differential analysis completed in {elapsed:.1f} seconds.")
            for name, res in results.items():
                st.write(f"**{name}**: total {len(res)}, Up {(res['regulation']=='Up').sum()}, Down {(res['regulation']=='Down').sum()}")

# ---------- Page: Visualization ----------
elif page == "Visualization":
    st.header("Visualization")
    if not st.session_state.de_results:
        st.warning("Please run differential analysis first.")
    else:
        comp_name = st.selectbox("Select comparison", list(st.session_state.de_results.keys()))
        res = st.session_state.de_results[comp_name]
        fig_vol = plot_volcano(res)
        st.plotly_chart(fig_vol)
        if st.button("Generate Heatmap"):
            start_time = time.time()
            sig = res[res['regulation'].isin(['Up','Down','Increase','Decrease'])]['Master protein IDs']
            if len(sig):
                mat = st.session_state.processed.set_index('Master protein IDs').loc[lambda x: x.index.isin(sig)]
                numeric_cols = mat.select_dtypes(include=[np.number]).columns
                if len(numeric_cols):
                    mat = mat[numeric_cols].apply(pd.to_numeric, errors='coerce')
                    heatmap_fig = plot_heatmap(mat)
                    st.pyplot(heatmap_fig)
                    elapsed = time.time() - start_time
                    st.success(f"Heatmap generated in {elapsed:.1f} seconds.")
                else:
                    st.warning("No numeric data for heatmap.")
            else:
                st.warning("No significant proteins.")

# ---------- Page: CD-Search ----------
elif page == "CD-Search":
    st.header("NCBI CD-Search")
    st.info("Note: The 'Accession' column in CD-Search results refers to NCBI conserved domain IDs, not UniProt protein accessions.")
    def ensure_fasta_header(text):
        lines = text.strip().splitlines()
        if lines and not lines[0].startswith(">"):
            st.info("No FASTA header detected. Adding '>query' automatically.")
            text = ">query\n" + text
        return text

    tab1, tab2, tab3 = st.tabs(["Simple Search", "Batch Search (with cache)", "Load from Cache"])
    with tab1:
        st.subheader("Simple Search (<950 sequences)")
        fasta_text = st.text_area("Paste FASTA sequences", height=200, key="cd_simple_fasta")
        if st.button("Run Search", key="cd_simple_btn"):
            start_time = time.time()
            if fasta_text.strip():
                fasta_text = ensure_fasta_header(fasta_text)
                with st.spinner("Searching..."):
                    result = batch_cd_search_simple(fasta_text)
                    if result is not None:
                        st.session_state.cd_result = result
                        elapsed = time.time() - start_time
                        st.success(f"Search completed in {elapsed:.1f} seconds. {len(result)} records")
                        st.dataframe(result)
                    else:
                        st.error("Search failed.")
            else:
                st.warning("Please enter FASTA sequences.")
    with tab2:
        st.subheader("Batch Search (auto-split + cache + resume)")
        st.markdown(f"Auto splits into chunks of {CD_CHUNK_SIZE}, saves to cache.")
        batch_fasta = st.text_area("Paste FASTA (supports large input)", height=200, key="cd_batch_fasta")
        if batch_fasta.strip():
            from utils.cd_search import parse_fasta_for_cd
            ids, seqs = parse_fasta_for_cd(ensure_fasta_header(batch_fasta))
            n = len(ids)
            batches = (n + CD_CHUNK_SIZE - 1)//CD_CHUNK_SIZE
            st.info(f"Detected **{n}** sequences, will be split into **{batches}** batches.")
        if st.button("Start Batch Search", key="cd_batch_btn", type="primary"):
            start_time = time.time()
            if not batch_fasta.strip():
                st.warning("No sequences provided.")
            else:
                batch_fasta = ensure_fasta_header(batch_fasta)
                progress_bar = st.progress(0)
                status_text = st.empty()
                def update_progress(frac):
                    progress_bar.progress(frac)
                    status_text.text(f"Progress: {int(frac*100)}%")
                with st.spinner("Searching (check terminal for progress)..."):
                    combined, task_dir = batch_cd_search_with_cache(batch_fasta, progress_callback=update_progress)
                progress_bar.progress(1.0)
                if combined is not None and not combined.empty:
                    st.session_state.cd_result = combined
                    elapsed = time.time() - start_time
                    st.success(f"Batch search completed in {elapsed:.1f} seconds. {len(combined)} records")
                    st.dataframe(combined.head(20))
                    csv = combined.to_csv(index=False).encode()
                    st.download_button("Download CSV", csv, file_name="cd_search_results.csv")
                else:
                    st.error("Search failed. Check terminal logs.")
    with tab3:
        st.subheader("Load from Cache")
        cached_tasks = get_cached_tasks()
        if not cached_tasks:
            st.info("No cached tasks found.")
        else:
            selected = st.selectbox("Choose task", cached_tasks)
            if st.button("Load", key="cd_cache_load"):
                start_time = time.time()
                df, name = load_cd_cache(selected)
                if df is not None and not df.empty:
                    st.session_state.cd_result = df
                    elapsed = time.time() - start_time
                    st.success(f"Loaded {len(df)} records in {elapsed:.1f} seconds.")
                    st.dataframe(df.head(20))
    if st.session_state.cd_result is not None:
        st.info(f"CD-Search result loaded: {st.session_state.cd_result.shape[0]} rows")

# ---------- Page: Export ----------
elif page == "Export":
    st.header("Export")
    if st.session_state.processed is None or not st.session_state.de_results:
        st.warning("Preprocessing and differential analysis must be completed first.")
    else:
        st.subheader("Download single plot")
        comp = st.selectbox("Comparison", list(st.session_state.de_results.keys()))
        fmt = st.selectbox("Format", ["svg","tiff"])
        w = st.number_input("Width (inch)", 5,30,10)
        h = st.number_input("Height (inch)", 5,30,8)
        title = st.text_input("Custom title (leave blank for default)")
        if st.button("Download Plot"):
            start_time = time.time()
            res = st.session_state.de_results[comp]
            fig = plot_volcano_single_annotated(res, title=title or comp)
            if fig:
                with tempfile.NamedTemporaryFile(suffix=f".{fmt}", delete=False) as tmp:
                    fig.savefig(tmp.name, dpi=300, bbox_inches='tight')
                    with open(tmp.name,'rb') as f:
                        st.download_button("Download", f, file_name=f"volcano_{comp}.{fmt}")
                elapsed = time.time() - start_time
                st.success(f"Plot generated in {elapsed:.1f} seconds.")
        st.subheader("Excel Report")
        if st.button("Generate Excel Report"):
            start_time = time.time()
            wb = create_excel_report(
                raw_data=st.session_state.raw_expr_df,
                clean_data=st.session_state.expr_df,
                norm_data=st.session_state.processed,
                de_results=st.session_state.de_results,
                comparisons=st.session_state.comparisons,
                groups=st.session_state.groups,
                sample_info=st.session_state.sample_info
            )
            with tempfile.NamedTemporaryFile(delete=False, suffix='.xlsx') as tmp:
                wb.save(tmp.name)
                with open(tmp.name,'rb') as f:
                    st.download_button('Download Report', f, file_name='proteomics_report.xlsx')
            elapsed = time.time() - start_time
            st.success(f"Excel report generated in {elapsed:.1f} seconds.")

# ---------- Page: EggNOG Annotation (REVISED) ----------
elif page == "EggNOG Annotation":
    st.header("EggNOG Annotation")
    st.markdown("Use local eggnog-mapper or upload manual result file.")
    
    # 数据库状态检查
    st.subheader("Database Status")
    data_dir = os.path.join(os.path.dirname(EMAPPER_PATH), "data")
    db_ready = os.path.exists(os.path.join(data_dir, "eggnog.db")) and os.path.exists(os.path.join(data_dir, "eggnog_proteins.dmnd"))
    if db_ready:
        st.success("Database ready")
    else:
        st.warning("Database not fully prepared")
    with st.expander("Download tools", expanded=False):
        c1,c2 = st.columns(2)
        with c1:
            if st.button("Download eggnog.db.gz"):
                start_time = time.time()
                url = "https://downloads.eggnogdb.org/emapper/emapperdb-5.0.2/eggnog.db.gz"
                dest = os.path.join(data_dir, "eggnog.db.gz")
                os.makedirs(data_dir, exist_ok=True)
                with st.spinner("Downloading..."):
                    r = requests.get(url, stream=True)
                    with open(dest, 'wb') as f:
                        for chunk in r.iter_content(1024*1024):
                            f.write(chunk)
                elapsed = time.time() - start_time
                st.success(f"Download finished in {elapsed:.1f} seconds.")
        with c2:
            if st.button("Decompress eggnog.db.gz"):
                start_time = time.time()
                gz = os.path.join(data_dir, "eggnog.db.gz")
                db = os.path.join(data_dir, "eggnog.db")
                if os.path.exists(gz):
                    import gzip
                    with gzip.open(gz,'rb') as src, open(db,'wb') as dst:
                        shutil.copyfileobj(src, dst)
                    elapsed = time.time() - start_time
                    st.success(f"Decompressed in {elapsed:.1f} seconds.")

    # ----- NEW: 如果已有注释结果，直接显示，无需重新运行 -----
    if st.session_state.eggnog_result is not None:
        st.success("Existing annotation result loaded from session.")
        df = st.session_state.eggnog_result
        st.dataframe(df)

        st.subheader("GO & KEGG Summary (from existing result)")
        # 加载映射文件
        go_obo_path = r"D:\myProject\go-basic.obo"
        kegg_file_path = r"D:\myProject\TBtools.KeggBackEnd"
        debug("Loading GO names from obo file...")
        go_name_map = load_go_names(go_obo_path)
        debug("Loading KO names from KEGG backend file...")
        ko_name_map = load_ko_names(kegg_file_path)

        if 'GOs' in df.columns:
            all_go = []
            for gos in df['GOs'].dropna():
                all_go.extend([g.strip() for g in gos.split(',') if g.strip() and g.strip() != '-'])
            if all_go:
                go_counts = pd.Series(all_go).value_counts().head(20)
                go_labels = [get_go_name(go_id, go_name_map) for go_id in go_counts.index]
                debug(f"GO plot labels (first 5): {go_labels[:5]}")
                fig_go = px.bar(x=go_labels, y=go_counts.values,
                                labels={'x':'GO Term','y':'Count'}, title='Top 20 GO Terms')
                st.plotly_chart(fig_go, width='stretch')
            else:
                debug("No valid GO terms found.")
        else:
            debug("'GOs' column not found in annotation result.")

        if 'KEGG_ko' in df.columns:
            all_ko = []
            for kos in df['KEGG_ko'].dropna():
                all_ko.extend([k.strip() for k in kos.split(',') if k.strip() and k.strip() != '-'])
            if all_ko:
                ko_counts = pd.Series(all_ko).value_counts().head(20)
                ko_labels = [get_ko_name(ko_id, ko_name_map) for ko_id in ko_counts.index]
                debug(f"KO plot labels (first 5): {ko_labels[:5]}")
                fig_ko = px.bar(x=ko_labels, y=ko_counts.values,
                                labels={'x':'KEGG KO','y':'Count'}, title='Top 20 KEGG Orthologs')
                st.plotly_chart(fig_ko, width='stretch')
            else:
                debug("No valid KO terms found.")
        else:
            debug("'KEGG_ko' column not found in annotation result.")

        st.info("You can still re-run annotation below if needed.")
        st.markdown("---")

    # 原有标签页：本地运行和手动上传
    tab_local, tab_upload = st.tabs(["Local Run", "Manual Upload"])
    with tab_local:
        st.subheader("Paste FASTA to run local emapper.py")
        if st.session_state.ref_fasta_text:
            st.info("Reference FASTA already loaded.")
        else:
            st.warning("No reference FASTA loaded. Upload here or in ESM2 page.")
            uploaded_ref = st.file_uploader("Upload reference FASTA", type=["fasta","fa","txt"], key="eggnog_ref_upload")
            if uploaded_ref:
                st.session_state.ref_fasta_text = uploaded_ref.getvalue().decode("utf-8")
                st.success("Reference FASTA loaded.")
        with st.expander("Auto-fill from DE proteins", expanded=True):
            if st.session_state.de_results:
                comp = st.selectbox("Comparison", list(st.session_state.de_results.keys()), key="eggnog_de_comp")
                res = st.session_state.de_results[comp]
                reg_filter = st.multiselect("Regulation filter", ['Up','Down','Increase','Decrease'], default=['Up','Down'], key="eggnog_de_reg")
                mask = res['regulation'].isin(reg_filter)
                protein_ids = res[mask]['Master protein IDs'].tolist()
                if protein_ids:
                    selected = st.multiselect("Proteins to annotate", protein_ids, default=protein_ids[:5] if len(protein_ids)>5 else protein_ids)
                    if st.button("Fill FASTA"):
                        if not st.session_state.ref_fasta_text:
                            st.error("No reference FASTA.")
                        else:
                            fasta_dict = parse_fasta(st.session_state.ref_fasta_text)
                            lines = [f">{pid}\n{fasta_dict[pid]}" for pid in selected if pid in fasta_dict]
                            if lines:
                                st.session_state.eggnog_fasta_text = "\n".join(lines)
                                st.success(f"Filled {len(lines)} sequences.")
                else:
                    st.info("No proteins matching filter.")
        st.text_area("FASTA sequences", height=220, key="eggnog_fasta_text")
        if st.button("Run Local Annotation"):
            fasta_content = st.session_state.eggnog_fasta_text
            if not fasta_content.strip():
                st.warning("No sequences entered.")
            else:
                start_time = time.time()
                df = run_eggnog_annotation_local(fasta_content)
                if df is not None:
                    st.session_state.eggnog_result = df
                    elapsed = time.time() - start_time
                    st.success(f"Annotation completed in {elapsed:.1f} seconds.")
                    # 立即显示结果
                    st.dataframe(df)
                    # 加载映射并绘图
                    go_obo_path = r"D:\myProject\go-basic.obo"
                    kegg_file_path = r"D:\myProject\TBtools.KeggBackEnd"
                    debug("Loading GO names from obo file...")
                    go_name_map = load_go_names(go_obo_path)
                    debug("Loading KO names from KEGG backend file...")
                    ko_name_map = load_ko_names(kegg_file_path)
                    st.subheader("GO & KEGG Summary")
                    if 'GOs' in df.columns:
                        all_go = []
                        for gos in df['GOs'].dropna():
                            all_go.extend([g.strip() for g in gos.split(',') if g.strip() and g.strip() != '-'])
                        if all_go:
                            go_counts = pd.Series(all_go).value_counts().head(20)
                            go_labels = [get_go_name(go_id, go_name_map) for go_id in go_counts.index]
                            debug(f"GO plot labels (first 5): {go_labels[:5]}")
                            fig_go = px.bar(x=go_labels, y=go_counts.values,
                                            labels={'x':'GO Term','y':'Count'}, title='Top 20 GO Terms')
                            st.plotly_chart(fig_go, width='stretch')
                    if 'KEGG_ko' in df.columns:
                        all_ko = []
                        for kos in df['KEGG_ko'].dropna():
                            all_ko.extend([k.strip() for k in kos.split(',') if k.strip() and k.strip() != '-'])
                        if all_ko:
                            ko_counts = pd.Series(all_ko).value_counts().head(20)
                            ko_labels = [get_ko_name(ko_id, ko_name_map) for ko_id in ko_counts.index]
                            debug(f"KO plot labels (first 5): {ko_labels[:5]}")
                            fig_ko = px.bar(x=ko_labels, y=ko_counts.values,
                                            labels={'x':'KEGG KO','y':'Count'}, title='Top 20 KEGG Orthologs')
                            st.plotly_chart(fig_ko, width='stretch')
                else:
                    st.error("Annotation failed.")
    with tab_upload:
        st.subheader("Upload eggNOG result file (TSV)")
        upload_file = st.file_uploader("Upload annotation file", type=['tsv','txt','csv'])
        if upload_file:
            content = upload_file.getvalue().decode('utf-8')
            df = parse_eggnog_manual_file(content)
            if df is not None and not df.empty:
                st.session_state.eggnog_result = df
                st.success(f"Parsed {len(df)} records.")
                st.dataframe(df)
                # 加载映射并绘图
                go_obo_path = r"D:\myProject\go-basic.obo"
                kegg_file_path = r"D:\myProject\TBtools.KeggBackEnd"
                debug("Loading GO names from obo file...")
                go_name_map = load_go_names(go_obo_path)
                debug("Loading KO names from KEGG backend file...")
                ko_name_map = load_ko_names(kegg_file_path)
                if 'GOs' in df.columns or 'KEGG_ko' in df.columns:
                    st.subheader("GO & KEGG Summary")
                    if 'GOs' in df.columns:
                        all_go = []
                        for gos in df['GOs'].dropna():
                            all_go.extend([g.strip() for g in gos.split(',') if g.strip() and g.strip() != '-'])
                        if all_go:
                            go_counts = pd.Series(all_go).value_counts().head(20)
                            go_labels = [get_go_name(go_id, go_name_map) for go_id in go_counts.index]
                            fig_go = px.bar(x=go_labels, y=go_counts.values,
                                            labels={'x':'GO Term','y':'Count'}, title='Top 20 GO Terms')
                            st.plotly_chart(fig_go, width='stretch')
                    if 'KEGG_ko' in df.columns:
                        all_ko = []
                        for kos in df['KEGG_ko'].dropna():
                            all_ko.extend([k.strip() for k in kos.split(',') if k.strip() and k.strip() != '-'])
                        if all_ko:
                            ko_counts = pd.Series(all_ko).value_counts().head(20)
                            ko_labels = [get_ko_name(ko_id, ko_name_map) for ko_id in ko_counts.index]
                            fig_ko = px.bar(x=ko_labels, y=ko_counts.values,
                                            labels={'x':'KEGG KO','y':'Count'}, title='Top 20 KEGG Orthologs')
                            st.plotly_chart(fig_ko, width='stretch')

# ---------- Page: ESM2 Similarity Search ----------
elif page == "ESM2 Similarity Search":
    st.header("ESM2 Protein Similarity Search")
    st.caption("Function 1: Single query similarity search & AI summary. Function 2: Batch annotation & enrichment.")
    if not ESM_AVAILABLE:
        st.error("ESM2 not available. Install PyTorch and fair-esm.")
        st.stop()

    if st.session_state.esm_model is None:
        start_time = time.time()
        with st.spinner("Loading ESM2 model (first time only)..."):
            model, alphabet, batch_converter = load_esm_model()
            if model is not None:
                st.session_state.esm_model = model
                st.session_state.esm_alphabet = alphabet
                st.session_state.esm_batch_converter = batch_converter
                elapsed = time.time() - start_time
                st.success(f"ESM2 model loaded in {elapsed:.1f} seconds.")
            else:
                st.error("Model loading failed")
                st.stop()

    model = st.session_state.esm_model
    batch_converter = st.session_state.esm_batch_converter

    st.subheader("Reference Library")
    pretrained_available = is_pretrained_available()
    use_pretrained = st.checkbox("Use pretrained Swiss-Prot library", value=pretrained_available, disabled=not pretrained_available)
    st.session_state.use_pretrained = use_pretrained

    if st.session_state.ref_fasta_text:
        st.info("Reference FASTA already loaded.")
        if not use_pretrained and "custom_library" not in st.session_state:
            if st.button("Build/Update Custom Library from existing FASTA", key="build_from_existing"):
                start_time = time.time()
                with st.spinner("Processing..."):
                    emb_dict, ids = build_reference_library(st.session_state.ref_fasta_text, model, batch_converter, use_cache=True)
                    st.session_state.custom_library = emb_dict
                    st.session_state.custom_ids = ids
                    elapsed = time.time() - start_time
                    st.success(f"Library built/updated in {elapsed:.1f} seconds. {len(emb_dict)} proteins")
    else:
        if not use_pretrained:
            custom_fasta = st.file_uploader("Upload custom reference FASTA", type=["fasta","fa","txt"], key="global_ref_fasta")
            if custom_fasta:
                fasta_text = custom_fasta.getvalue().decode("utf-8")
                st.session_state.ref_fasta_text = fasta_text
                hash_val = compute_fasta_hash(fasta_text)
                if os.path.exists(f"data/custom_cache_{hash_val}.npz"):
                    st.info("Cache detected, will be loaded directly.")
                if st.button("Build/Update Custom Library"):
                    start_time = time.time()
                    with st.spinner("Processing..."):
                        emb_dict, ids = build_reference_library(fasta_text, model, batch_converter, use_cache=True)
                        st.session_state.custom_library = emb_dict
                        st.session_state.custom_ids = ids
                        elapsed = time.time() - start_time
                        st.success(f"Library built in {elapsed:.1f} seconds. {len(emb_dict)} proteins")

    if use_pretrained:
        st.info("Using pretrained Swiss-Prot library.")
        if "custom_library" in st.session_state:
            del st.session_state.custom_library
    elif "custom_library" in st.session_state and st.session_state.custom_library is not None:
        st.info(f"Custom library with {len(st.session_state.custom_library)} proteins.")
    else:
        st.info("No library selected or built yet.")

    st.markdown("---")

    def show_de_protein_selector(label="Select DE proteins", key_suffix=""):
        if not st.session_state.get("de_results"):
            st.warning("No differential analysis results available.")
            return None, []
        comps = list(st.session_state.de_results.keys())
        comp_name = st.selectbox("Comparison", comps, key=f"esm_de_comp_{key_suffix}")
        res = st.session_state.de_results[comp_name]
        reg_filter = st.multiselect("Regulation filter", ['Up','Down','Increase','Decrease'], default=['Up','Down'], key=f"esm_de_reg_{key_suffix}")
        if not reg_filter:
            return comp_name, []
        mask = res['regulation'].isin(reg_filter)
        protein_ids = res[mask]['Master protein IDs'].tolist()
        if not protein_ids:
            st.info("No proteins match the selected filters.")
            return comp_name, []
        selected = st.multiselect(label, protein_ids, key=f"esm_de_sel_{key_suffix}")
        return comp_name, selected

    # Function 1
    with st.expander("Function 1: Single query similarity & AI summary", expanded=True):
        st.markdown("Enter a protein ID or sequence, or select from DE proteins.")
        with st.expander("Or select from differential expression proteins", expanded=False):
            comp_name_f1, selected_f1 = show_de_protein_selector("Select one protein", "f1")
            if selected_f1:
                first_pid = selected_f1[0]
                st.info(f"Will search for: {first_pid}")
                st.session_state.query_id_f1 = first_pid
                if st.session_state.ref_fasta_text:
                    fasta_dict = parse_fasta(st.session_state.ref_fasta_text)
                    seq = fasta_dict.get(first_pid, "")
                    st.session_state.query_seq_f1 = seq if seq else ""
        c1, c2 = st.columns(2)
        with c1:
            query_id = st.text_input("Protein ID", key="query_id_f1")
        with c2:
            query_seq = st.text_area("Or paste protein sequence", key="query_seq_f1", height=120)
        top_n = st.slider("Number of results", 5, 50, 10)
        if st.button("Search Similar Proteins", key="single_search"):
            start_time = time.time()
            final_seq = None
            if query_seq.strip():
                lines = query_seq.strip().splitlines()
                final_seq = "".join([l for l in lines if not l.startswith(">")])
            elif query_id.strip():
                final_seq = get_sequence_for_id(query_id.strip(), fasta_text=st.session_state.ref_fasta_text)
                if final_seq is None:
                    st.error("Could not retrieve sequence.")
                    st.stop()
            else:
                st.error("Provide either ID or sequence.")
                st.stop()
            custom_lib = st.session_state.get("custom_library")
            result_df = find_similar_proteins(final_seq, top_n, model, batch_converter,
                                              use_pretrained=use_pretrained, custom_library=custom_lib)
            if "Error" in result_df.columns:
                st.error(result_df["Error"][0])
            else:
                if st.session_state.get("eggnog_result") is not None:
                    result_df = enrich_esm_result_with_eggnog(result_df, st.session_state.eggnog_result)
                st.session_state.last_search_result = result_df
                elapsed = time.time() - start_time
                st.success(f"Search completed in {elapsed:.1f} seconds. Found {len(result_df)} proteins.")
        if "last_search_result" in st.session_state and st.session_state.last_search_result is not None:
            res = st.session_state.last_search_result
            st.dataframe(res, width='stretch')
            if 'Gene' in res.columns and (res['Gene'] == '').any():
                st.info("Only proteins present in the eggNOG annotation will have Gene/Protein_Name/GO/EC information.")
            st.subheader("AI Functional Summary")
            selected_idx = st.selectbox("Select protein:", range(len(res)),
                                        format_func=lambda i: f"#{i+1}: {res.iloc[i]['Protein_ID']} (sim {res.iloc[i]['Similarity']:.3f})")
            use_api = st.checkbox("Use DeepSeek AI (requires API key)", key="use_api_summary")
            if st.button("Generate Summary", key="gen_summary"):
                start_time = time.time()
                row = res.iloc[selected_idx]
                summary = generate_function_summary(
                    protein_id=row['Protein_ID'], protein_name=row.get('Protein_Name',''),
                    similarity_score=row['Similarity'], similar_protein=res.iloc[0]['Protein_ID'],
                    go_terms=row.get('GO',''), ec_numbers=row.get('EC',''), use_api=use_api
                )
                st.session_state.current_summary = summary
                elapsed = time.time() - start_time
                st.success(f"Summary generated in {elapsed:.1f} seconds.")
            if 'current_summary' in st.session_state:
                st.success(st.session_state.current_summary)

    # Function 2
    with st.expander("Function 2: Batch annotation & enrichment", expanded=False):
        st.markdown("Input multiple sequences or select from DE proteins, or directly run GO enrichment on all DE proteins.")
        st.info("In 'ESM similarity + enrichment' mode, you can directly search selected proteins without pasting sequences manually.")
        analysis_mode = st.radio("Analysis mode", ["ESM similarity + enrichment", "Direct GO enrichment (all DE proteins)"], horizontal=True, key="batch_analysis_mode")

        if "batch_search_in_progress" not in st.session_state:
            st.session_state.batch_search_in_progress = False
            st.session_state.batch_search_index = 0
            st.session_state.batch_search_results = []
            st.session_state.batch_search_total = 0
            st.session_state.batch_search_seqs = []

        if analysis_mode == "Direct GO enrichment (all DE proteins)":
            st.info("This will run enrichment analysis directly on all differentially expressed proteins (using current comparison and regulation filter) without ESM similarity search.")
            comps = list(st.session_state.de_results.keys())
            comp_name = st.selectbox("Comparison", comps, key="batch_dir_comp")
            res = st.session_state.de_results[comp_name]
            reg_filter = st.multiselect("Regulation filter", ['Up','Down','Increase','Decrease'], default=['Up','Down'], key="batch_dir_reg")
            if st.button("Run GO Enrichment on all filtered DE proteins", key="batch_direct_go"):
                if not reg_filter:
                    st.warning("Please select at least one regulation type.")
                else:
                    mask = res['regulation'].isin(reg_filter)
                    target_ids = res[mask]['Master protein IDs'].tolist()
                    if not target_ids:
                        st.warning("No proteins match the filter.")
                    else:
                        annotations = st.session_state.get('custom_annotations')
                        if not annotations and st.session_state.get("eggnog_result") is not None:
                            annotations = build_annotations_from_eggnog(st.session_state.eggnog_result)
                            if annotations:
                                st.session_state.custom_annotations = annotations
                        if not annotations:
                            st.error("No annotations available. Please run eggNOG annotation or load UniProt annotations first.")
                        else:
                            display_batch_results_and_enrichment([], target_ids, annotations)
        else:
            with st.expander("Or select multiple DE proteins", expanded=True):
                comp_name_f2, selected_f2 = show_de_protein_selector("Select proteins for batch", "f2")
                debug(f"ESM2 batch: selected proteins = {selected_f2}")

                col1, col2 = st.columns(2)
                with col1:
                    if st.button("Fill batch sequences", key="fill_batch_btn", help="Extract sequences and fill the text area below."):
                        debug("Fill batch sequences button clicked.")
                        if not st.session_state.ref_fasta_text:
                            st.error("Reference FASTA is required.")
                        elif not selected_f2:
                            st.warning("Please select at least one protein.")
                        else:
                            fasta_dict = parse_fasta(st.session_state.ref_fasta_text)
                            seqs = []
                            for pid in selected_f2:
                                seq = fasta_dict.get(pid, "")
                                if seq:
                                    seqs.append(seq)
                                else:
                                    st.warning(f"Sequence not found for {pid}")
                            if seqs:
                                text_lines = []
                                for i, (pid, s) in enumerate(zip(selected_f2, seqs)):
                                    text_lines.append(f">{pid}")
                                    text_lines.append(s)
                                st.session_state.batch_text = "\n".join(text_lines)
                                st.success(f"Filled {len(seqs)} sequences into the text area below.")
                            else:
                                st.warning("No sequences could be extracted.")
                with col2:
                    if st.button("Run ESM Search on Selected Proteins", key="run_esm_selection", help="Directly run ESM search using the selected proteins' sequences."):
                        debug("Run ESM Search button clicked.")
                        if not st.session_state.ref_fasta_text:
                            st.error("Reference FASTA is required.")
                        elif "custom_library" not in st.session_state or not st.session_state.custom_library:
                            st.error("Build reference library first.")
                        elif not selected_f2:
                            st.warning("Please select at least one protein.")
                        else:
                            fasta_dict = parse_fasta(st.session_state.ref_fasta_text)
                            seqs = []
                            for pid in selected_f2:
                                seq = fasta_dict.get(pid, "")
                                if seq:
                                    seqs.append(seq)
                                else:
                                    st.warning(f"Sequence not found for {pid}")
                            if seqs:
                                st.session_state.batch_final_results = None
                                init_batch_search_state(st.session_state, len(seqs))
                                st.session_state.batch_search_seqs = seqs
                                st.rerun()
                            else:
                                st.warning("No sequences could be extracted.")

            if st.session_state.get("batch_final_results") and not st.session_state.batch_search_in_progress:
                results = st.session_state.batch_final_results
                target_ids = st.session_state.batch_target_ids or []
                annotations = st.session_state.get('custom_annotations')
                if not annotations and st.session_state.get("eggnog_result") is not None:
                    annotations = build_annotations_from_eggnog(st.session_state.eggnog_result)
                    if annotations:
                        st.session_state.custom_annotations = annotations
                display_batch_results_and_enrichment(results, target_ids, annotations or {})

            elif st.session_state.batch_search_in_progress:
                total = st.session_state.batch_search_total
                current = st.session_state.batch_search_index
                st.progress(current/total if total else 0)
                st.write(f"Processing: {current}/{total}")
                lib = st.session_state.get("custom_library")
                if lib is None:
                    st.error("Library lost.")
                    st.session_state.batch_search_in_progress = False
                    st.rerun()
                else:
                    done = run_batch_search_step(st.session_state.batch_search_seqs, model, batch_converter, lib, st.session_state, 100)
                    if done:
                        results = st.session_state.batch_search_results
                        target_ids = [r['best_id'] for r in results if r['best_id']]
                        st.session_state.batch_final_results = results
                        st.session_state.batch_target_ids = target_ids
                        annotations = st.session_state.get('custom_annotations')
                        if not annotations and st.session_state.get("eggnog_result") is not None:
                            annotations = build_annotations_from_eggnog(st.session_state.eggnog_result)
                            if annotations:
                                st.session_state.custom_annotations = annotations
                        display_batch_results_and_enrichment(results, target_ids, annotations or {})
                        st.session_state.batch_search_in_progress = False
                        st.rerun()
                    else:
                        st.rerun()
                st.stop()

            else:
                st.markdown("---")
                st.write("Or paste sequences manually:")
                input_mode = st.radio("Input method", ["Paste sequences", "Upload FASTA"], horizontal=True)
                batch_seqs = []
                if input_mode == "Paste sequences":
                    txt = st.text_area("One sequence per line (FASTA headers allowed)", key="batch_text", height=150)
                    if txt:
                        lines = txt.strip().split('\n')
                        seq = ""
                        for line in lines:
                            if line.startswith('>'):
                                if seq: batch_seqs.append(seq); seq = ""
                            else: seq += line.strip()
                        if seq: batch_seqs.append(seq)
                else:
                    batch_file = st.file_uploader("Upload FASTA", type=["fasta","fa","txt"], key="batch_file")
                    if batch_file:
                        batch_seqs = list(parse_fasta(batch_file.read().decode()).values())
                        st.success(f"Parsed {len(batch_seqs)} sequences")

                if st.button("Load UniProt annotations for current library", disabled=(use_pretrained or "custom_library" not in st.session_state)):
                    start_time = time.time()
                    lib_ids = list(st.session_state.custom_library.keys())
                    if lib_ids:
                        valid = sum(is_uniprot_accession(extract_accession(x)) for x in lib_ids[:20])
                        if valid < 2:
                            st.warning("Library IDs are not standard UniProt.")
                        else:
                            with st.spinner("Fetching annotations..."):
                                annots = fetch_uniprot_annotations_batch(lib_ids)
                                st.session_state.custom_annotations = annots
                                elapsed = time.time() - start_time
                                st.success(f"Fetched {len(annots)} annotations in {elapsed:.1f} seconds.")
                    else:
                        st.warning("Library empty.")

                if st.button("Run Batch Annotation & Enrichment", disabled=(len(batch_seqs)==0)):
                    if "custom_library" not in st.session_state or not st.session_state.custom_library:
                        st.error("Build reference library first.")
                    else:
                        st.session_state.batch_final_results = None
                        init_batch_search_state(st.session_state, len(batch_seqs))
                        st.session_state.batch_search_seqs = batch_seqs
                        st.rerun()

# ---------- Page: STRING PPI Network ----------
elif page == "STRING PPI Network":
    st.header("STRING PPI Network")
    st.info("Query protein-protein interactions via STRING. No local database required.")
    source_option = st.radio("Protein source", ["Manual input / FASTA", "Select from DE results"])
    fasta_text = ""
    fasta_seqs = {}
    original_ids = []
    selected_de_protein = None
    if source_option == "Manual input / FASTA":
        fasta_text = st.text_area("Paste FASTA sequence(s)", height=200, placeholder=">Pt_Chr0100005\nMKT...")
    else:
        if not st.session_state.de_results:
            st.warning("No DE results.")
        else:
            comp = st.selectbox("Comparison", list(st.session_state.de_results.keys()), key="string_de_comp")
            res = st.session_state.de_results[comp]
            reg_filter = st.multiselect("Regulation filter", ['Up','Down','Increase','Decrease'], default=['Up','Down'], key="string_de_reg")
            mask = res['regulation'].isin(reg_filter)
            protein_ids = res[mask]['Master protein IDs'].tolist()
            if protein_ids:
                selected_de_protein = st.selectbox("Protein of interest", protein_ids)
                auto_seq = ""
                if st.session_state.ref_fasta_text:
                    fasta_dict = parse_fasta(st.session_state.ref_fasta_text)
                    if selected_de_protein in fasta_dict:
                        auto_seq = fasta_dict[selected_de_protein]
                if not auto_seq:
                    st.warning("Sequence not found in reference FASTA, please paste manually.")
                else:
                    st.success("Sequence auto-filled from reference.")
                fasta_text = st.text_area("Sequence", value=auto_seq, height=150)
            else:
                st.warning("No proteins matching filter.")
    if fasta_text.strip():
        fasta_seqs = parse_fasta(fasta_text)
        original_ids = list(fasta_seqs.keys())
        if not fasta_seqs:
            fasta_text = ">query\n" + fasta_text.strip()
            fasta_seqs = parse_fasta(fasta_text)
            original_ids = list(fasta_seqs.keys())
        st.success(f"Parsed {len(original_ids)} sequences")
    if source_option == "Select from DE results" and selected_de_protein:
        if selected_de_protein not in original_ids:
            original_ids.append(selected_de_protein)

    st.subheader("Mapping method")
    mapping_method = st.radio("Choose method", ["Manual UniProt ID", "Use eggNOG", "BLAST sequence"], index=2)
    manual_uniprot_text = ""
    if mapping_method == "Manual UniProt ID":
        manual_uniprot_text = st.text_area("UniProt Accession per line", height=150, placeholder="P93004\nQ9S7W5")

    st.subheader("Species & parameters")
    species_presets = {
        "A. thaliana (3702)": 3702, "B. distachyon (15368)": 15368, "B. taurus (9913)": 9913,
        "C. elegans (6239)": 6239, "D. melanogaster (7227)": 7227, "D. rerio (7955)": 7955,
        "E. coli K-12 (83333)": 83333, "G. gallus (9031)": 9031, "H. sapiens (9606)": 9606,
        "M. musculus (10090)": 10090, "O. sativa Japonica (39947)": 39947, "R. norvegicus (10116)": 10116,
        "S. cerevisiae (4932)": 4932, "S. pombe (4896)": 4896, "S. lycopersicum (4081)": 4081,
        "S. tuberosum (4113)": 4113, "T. aestivum (4565)": 4565, "Z. mays (4577)": 4577,
    }
    common_names = {
        3702: "Thale cress", 15368: "Purple false brome", 9913: "Cattle", 6239: "Roundworm",
        7227: "Fruit fly", 7955: "Zebrafish", 83333: "E. coli K-12", 9031: "Chicken",
        9606: "Human", 10090: "House mouse", 39947: "Rice", 10116: "Norway rat",
        4932: "Baker's yeast", 4896: "Fission yeast", 4081: "Tomato", 4113: "Potato",
        4565: "Common wheat", 4577: "Maize",
    }
    auto_species = 3702
    eggnog_df = st.session_state.get("eggnog_result")
    if eggnog_df is not None and not eggnog_df.empty:
        seed_col = None
        for col in eggnog_df.columns:
            if 'seed_ortholog' in col.lower():
                seed_col = col
                break
        if seed_col:
            first_seed = str(eggnog_df[seed_col].iloc[0])
            m = re.match(r'(\d+)\.', first_seed)
            if m:
                auto_species = int(m.group(1))
                species_name = [k for k,v in species_presets.items() if v == auto_species]
                common = common_names.get(auto_species, "")
                if species_name:
                    st.info(f"Detected species from eggNOG: {species_name[0]} ({common}) - ID {auto_species}")
                else:
                    st.info(f"Detected species ID from eggNOG: {auto_species}{' ('+common+')' if common else ''}")
    else:
        st.info("If you don't know the species, common model organisms are listed below. You can also enter any NCBI Taxonomy ID.")
    default_species = auto_species if auto_species in species_presets.values() else 3702
    preset_label = st.selectbox("Preset species", list(species_presets.keys()),
                                index=list(species_presets.values()).index(default_species) if default_species in species_presets.values() else 0)
    species = species_presets[preset_label]
    custom_species = st.text_input("Or enter NCBI Taxonomy ID", value=str(species))
    try:
        species = int(custom_species)
    except:
        st.warning("Invalid species ID, using preset.")
    required_score = st.slider("Required score", 0, 1000, 400)
    if st.button("Fetch STRING Network"):
        start_time = time.time()
        id_mapping = {}
        with st.spinner("Preparing UniProt IDs..."):
            if mapping_method == "Manual UniProt ID":
                if not manual_uniprot_text.strip():
                    st.error("No UniProt IDs.")
                    st.stop()
                manual_ids = parse_manual_uniprot_ids(manual_uniprot_text)
                if not manual_ids: st.stop()
                id_mapping = {uid:uid for uid in manual_ids}
            elif mapping_method == "Use eggNOG":
                eggnog_df = st.session_state.get("eggnog_result")
                if eggnog_df is None:
                    st.error("No eggNOG result.")
                    st.stop()
                query_ids = original_ids if original_ids else [selected_de_protein] if selected_de_protein else []
                id_mapping = extract_uniprot_from_eggnog(eggnog_df, query_ids, fallback_to_seed=True)
                if not id_mapping:
                    st.error("No UniProt extracted. Try BLAST or manual input.")
                else:
                    st.success(f"Extracted {len(id_mapping)} mappings.")
            else:
                if not fasta_text.strip(): st.error("No FASTA input."); st.stop()
                if not fasta_seqs: st.error("Could not parse FASTA."); st.stop()
                id_mapping = run_blast_mapping(fasta_seqs, max_retries=2, max_wait_per_job=180)
                if not id_mapping: st.error("BLAST failed.")
        if not id_mapping: st.stop()
        uniprot_ids = list(dict.fromkeys(id_mapping.values()))
        map_df = pd.DataFrame(list(id_mapping.items()), columns=["Original", "UniProt"])
        valid_pattern = re.compile(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?$')
        map_df["Valid UniProt"] = map_df["UniProt"].apply(lambda x: "Yes" if valid_pattern.match(str(x)) else "No")
        st.dataframe(map_df, width='stretch')
        if all(map_df["Valid UniProt"] == "No"):
            st.warning("None of the extracted IDs are standard UniProt accessions. STRING may not recognize them. Consider using BLAST or manual UniProt input.")
        else:
            st.info("UniProt accessions marked 'Yes' are suitable for direct input in 'Manual UniProt ID'.")
        with st.spinner("Querying STRING..."):
            network_json = call_string_api(uniprot_ids, species=species, required_score=required_score)
        if network_json is None:
            st.error("STRING API call failed.")
        elif len(network_json) == 0:
            st.warning("No interactions found.")
        else:
            elapsed = time.time() - start_time
            st.success(f"Network with {len(network_json)} edges fetched in {elapsed:.1f} seconds.")
            html = build_ppi_network_html(network_json, id_mapping)
            st.components.v1.html(html, height=800, scrolling=True)
            st.download_button("Download JSON", data=json.dumps(network_json, indent=2), file_name="string_network.json")

# ---------- Page: Integrated Analysis ----------
elif page == "Integrated Analysis":
    st.header("Integrated Analysis")
    st.markdown("Automatically uses eggNOG annotations to filter CD-Search, run ESM2 similarity, and query STRING.")

    if st.session_state.eggnog_result is None:
        st.warning("No eggNOG annotation result found. Please run eggNOG annotation first.")
        st.stop()

    eggnog_df = st.session_state.eggnog_result
    egg_ids = eggnog_df[eggnog_df.columns[0]].tolist()

    use_eggnog = st.checkbox("Use eggNOG annotations for integrated analysis", value=True)

    # CD-Search filtering
    st.subheader("CD-Search Results (filtered by eggNOG proteins)")
    if st.session_state.cd_result is not None:
        cd_df = st.session_state.cd_result
        id_col = cd_df.columns[0]
        cd_filtered = cd_df[cd_df[id_col].isin(egg_ids)] if use_eggnog else cd_df
        st.dataframe(cd_filtered.head(20), width='stretch')
        st.info(f"Showing {len(cd_filtered)} out of {len(cd_df)} CD-Search records")
    else:
        st.warning("No CD-Search results loaded.")

    # ESM2 Similarity Search
    st.subheader("ESM2 Similarity Search (using first eggNOG protein)")
    if st.session_state.esm_model is None:
        st.warning("ESM2 model not loaded yet. Please go to ESM2 page once to load the model.")
    else:
        if use_eggnog:
            first_egg_id = egg_ids[0]
            st.info(f"Will search for: {first_egg_id}")
            if st.session_state.ref_fasta_text:
                fasta_dict = parse_fasta(st.session_state.ref_fasta_text)
                seq = fasta_dict.get(first_egg_id, "")
                if not seq:
                    st.error(f"Sequence for {first_egg_id} not found in reference FASTA.")
                else:
                    if st.button("Run ESM2 Similarity Search"):
                        with st.spinner("Searching..."):
                            result_df = find_similar_proteins(
                                seq, top_n=10,
                                model=st.session_state.esm_model,
                                batch_converter=st.session_state.esm_batch_converter,
                                use_pretrained=st.session_state.get("use_pretrained", False),
                                custom_library=st.session_state.get("custom_library")
                            )
                            if "Error" in result_df.columns:
                                st.error(result_df["Error"][0])
                            else:
                                result_df = enrich_esm_result_with_eggnog(result_df, eggnog_df)
                                st.session_state.last_search_result = result_df
                                st.success("Search completed.")
            else:
                st.error("Reference FASTA not loaded. Please upload it in ESM2 or EggNOG page.")
            if st.session_state.get("last_search_result") is not None:
                st.dataframe(st.session_state.last_search_result, width='stretch')
                if 'Gene' in st.session_state.last_search_result.columns and (st.session_state.last_search_result['Gene'] == '').any():
                    st.info("Only proteins present in the eggNOG annotation will have Gene/Protein_Name/GO/EC information.")
        else:
            st.info("Uncheck the toggle to manually control ESM2 from its dedicated page.")

    # STRING PPI Network
    st.subheader("STRING PPI Network (using eggNOG UniProt IDs)")
    if use_eggnog:
        target_ids = egg_ids[:5]
        id_mapping = extract_uniprot_from_eggnog(eggnog_df, target_ids, fallback_to_seed=True)
        if not id_mapping:
            st.warning("No valid UniProt IDs could be extracted from eggNOG. You can try BLAST or manual input on the STRING page.")
        else:
            uniprot_ids = list(dict.fromkeys(id_mapping.values()))
            st.info(f"Extracted {len(uniprot_ids)} UniProt IDs: {', '.join(uniprot_ids[:5])}")
            species = 3702
            seed_col = None
            for col in eggnog_df.columns:
                if 'seed_ortholog' in col.lower():
                    seed_col = col
                    break
            if seed_col:
                first_seed = str(eggnog_df[seed_col].iloc[0])
                m = re.match(r'(\d+)\.', first_seed)
                if m:
                    species = int(m.group(1))
            required_score = st.slider("Required score", 0, 1000, 400, key="integrated_string_score")
            if st.button("Fetch STRING Network"):
                with st.spinner("Querying STRING..."):
                    network_json = call_string_api(uniprot_ids, species=species, required_score=required_score)
                if network_json is None:
                    st.error("STRING API call failed.")
                elif len(network_json) == 0:
                    st.warning("No interactions found.")
                else:
                    html = build_ppi_network_html(network_json, id_mapping)
                    st.components.v1.html(html, height=600, scrolling=True)
                    st.download_button("Download JSON", data=json.dumps(network_json, indent=2), file_name="string_network.json")
    else:
        st.info("Uncheck the toggle to manually use STRING from its dedicated page.")

    # AI Functional Summary
    st.subheader("AI Functional Summary")
    st.markdown("Generate an AI-powered functional summary for a selected eggNOG-annotated protein.")
    protein_list = egg_ids
    selected_protein = st.selectbox("Select a protein from eggNOG:", protein_list)
    use_api = st.checkbox("Use DeepSeek AI (requires API key)", value=False, key="integrated_use_api")
    if st.button("Generate Summary", key="integrated_gen_summary"):
        protein_row = eggnog_df[eggnog_df[eggnog_df.columns[0]] == selected_protein].iloc[0]
        protein_name = str(protein_row.get('Description', '') or protein_row.get('Preferred_name', ''))
        go_terms = str(protein_row.get('GOs', ''))
        ec_numbers = str(protein_row.get('EC', ''))
        start_time = time.time()
        summary = generate_function_summary(
            protein_id=selected_protein,
            protein_name=protein_name,
            similarity_score=None,
            similar_protein="",
            go_terms=go_terms,
            ec_numbers=ec_numbers,
            use_api=use_api
        )
        elapsed = time.time() - start_time
        st.success(f"Summary generated in {elapsed:.1f} seconds.")
        st.write(summary)

    # Export integrated report
    st.subheader("Export Integrated Report")
    if st.button("Generate Integrated Report CSV"):
        export_parts = {}
        if st.session_state.cd_result is not None:
            cd_df = st.session_state.cd_result
            if use_eggnog:
                id_col = cd_df.columns[0]
                cd_sub = cd_df[cd_df[id_col].isin(egg_ids)]
                export_parts["CD_Search"] = cd_sub
            else:
                export_parts["CD_Search"] = cd_df
        if st.session_state.get("last_search_result") is not None:
            export_parts["ESM2_Search"] = st.session_state.last_search_result
        if use_eggnog:
            export_parts["EggNOG"] = eggnog_df
        if export_parts:
            combined_list = []
            for source, df in export_parts.items():
                temp = df.copy()
                temp["Source"] = source
                combined_list.append(temp)
            combined = pd.concat(combined_list, ignore_index=True)
            csv = combined.to_csv(index=False).encode()
            st.download_button("Download CSV", csv, file_name="integrated_analysis.csv")
        else:
            st.info("No data to export.")

# ---------- End of pages ----------
