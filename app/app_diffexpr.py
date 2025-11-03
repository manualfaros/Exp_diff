from __future__ import annotations
import io, re, base64
from typing import Optional, Dict
import numpy as np
import pandas as pd
import altair as alt
import streamlit as st

# ----------------------------
# Configuración de página
# ----------------------------
st.set_page_config(page_title="EDA expresión diferencial", page_icon="🧬", layout="wide")

# ----------------------------
# Tema Altair
# ----------------------------
def make_alt_theme():
    def _theme():
        return {"config": {
            "axis": {"labelFontSize": 12, "titleFontSize": 13},
            "legend": {"labelFontSize": 12, "titleFontSize": 13},
            "view": {"strokeWidth": 0},
            "title": {"fontSize": 16}
        }}
    return _theme

alt.themes.register("clean", make_alt_theme())
alt.themes.enable("clean")

# ----------------------------
# Utilidades
# ----------------------------
def bytes_download_link(data: bytes, filename: str, label: str) -> str:
    b64 = base64.b64encode(data).decode()
    return f'<a download="{filename}" href="data:file/octet-stream;base64,{b64}">{label}</a>'
@st.cache_data(show_spinner=False)
def read_table(file_bytes: bytes, filename: str, sep_choice: str, nrows_preview: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Lee la tabla a partir de los *bytes* del archivo, para evitar problemas
    con el puntero del UploadedFile y el cache.
    Usa:
      - \t por defecto para .tsv
      - autodetección para otros si sep_choice == "auto"
    """
    name = filename.lower()

    # Caso Excel
    if name.endswith((".xlsx", ".xls")):
        df = pd.read_excel(io.BytesIO(file_bytes))
    else:
        # Determinar el separador real
        if sep_choice == "auto":
            if name.endswith(".tsv"):
                # Para TSV, asumimos tabulador
                sep_real = "\t"
            else:
                # Sniff en un buffer separado
                sniff_df = pd.read_csv(io.BytesIO(file_bytes), sep=None, engine="python", nrows=1000)
                sep_real = sniff_df.attrs.get("delimiter", ",")
        elif sep_choice == "\\t":
            sep_real = "\t"
        else:
            sep_real = sep_choice

        # Lectura definitiva en otro buffer limpio
        df = pd.read_csv(io.BytesIO(file_bytes), sep=sep_real)

    preview = df.head(nrows_preview).copy()
    return df, preview


def detectar_contrastes(df: pd.DataFrame) -> Dict[str, Dict[str, Optional[str]]]:
    patrones = {
        "logFC": re.compile(r"(.+?)_logFC$", re.IGNORECASE),
        "AveExpr": re.compile(r"(.+?)_AveExpr$", re.IGNORECASE),
        "adjP": re.compile(r"(.+?)_adj\.P\.Val$", re.IGNORECASE),
        "P": re.compile(r"(.+?)_P\.Value$", re.IGNORECASE),
    }
    info = {}
    for col in df.columns:
        col = str(col)
        for tipo, pat in patrones.items():
            m = pat.match(col)
            if m:
                pref = m.group(1)
                if pref not in info:
                    info[pref] = {"logFC": None, "AveExpr": None, "adjP": None, "P": None}
                info[pref][tipo if tipo != "adjP" else "adjP"] = col
                break
    validos = {}
    for pref, d in info.items():
        if d["logFC"] and (d["adjP"] or d["P"]):
            validos[pref] = d
    return validos

def elegir_columna_gen(df: pd.DataFrame) -> str | None:
    candidatos = ["SYMBOL", "gene_id", "Geneid", "GeneID", "gene", "Gene", "locus_tag"]
    for c in candidatos:
        if c in df.columns:
            return c
    return None

def pretty_contrast_name(pref: str) -> str:
    return re.sub(r"^PRJNA\d+_[^_]+_", "", pref)

def preparar_deg_por_contraste(df: pd.DataFrame, gene_col: str, contraste: str, info_contrastes):
    meta = info_contrastes[contraste]
    logfc_col = meta["logFC"]
    ave_col = meta["AveExpr"]
    padj_col = meta["adjP"] if meta["adjP"] else meta["P"]
    work = df.rename(columns={
        gene_col: "gene",
        logfc_col: "logFC",
        padj_col: "padj"
    }).copy()
    if ave_col:
        work = work.rename(columns={ave_col: "AveExpr"})
    else:
        work["AveExpr"] = np.nan
    work["padj"] = pd.to_numeric(work["padj"], errors="coerce")
    work["logFC"] = pd.to_numeric(work["logFC"], errors="coerce")
    work["neg_log10_padj"] = -np.log10(np.clip(work["padj"], 1e-300, 1.0))
    return work[["gene", "logFC", "AveExpr", "padj", "neg_log10_padj"]]

# ----------------------------
# Sidebar: carga
# ----------------------------
st.sidebar.title("⚙️ Carga de datos")

uploaded = st.sidebar.file_uploader("Sube tabla de expresión diferencial (TSV/CSV/XLSX)", type=["tsv","csv","txt","xlsx","xls"])
sep_choice = st.sidebar.selectbox("Separador", ["auto", ",", ";", "\\t", "|"])
preview_n = st.sidebar.number_input("Filas de preview", 5, 100, 10)

df = None
info_contrastes = {}
gene_col = None
if uploaded is not None:
    file_bytes = uploaded.getvalue()
    df, preview = read_table(file_bytes, uploaded.name, sep_choice, preview_n)
    st.sidebar.success(f"{df.shape[0]} filas × {df.shape[1]} columnas")
    info_contrastes = detectar_contrastes(df)
    st.sidebar.info(f"{len(info_contrastes)} contrastes detectados")
    gene_col = elegir_columna_gen(df)
    if gene_col:
        st.sidebar.write(f"Columna de gen detectada: `{gene_col}`")
    else:
        gene_col = st.sidebar.selectbox("Selecciona la columna de gen", df.columns)

# ----------------------------
# Tabs
# ----------------------------
tabs = st.tabs(["📊 Exploración", "🌋 Volcano / MA", "📑 DEGs", "🔍 Gen"])

# ----------------------------
# Tab 1: Exploración
# ----------------------------
with tabs[0]:
    st.header("📊 Exploración general")
    if df is None:
        st.info("Sube una tabla para comenzar.")
    else:
        st.subheader("Preview")
        st.dataframe(preview, use_container_width=True)
        st.subheader("Contrastes detectados")
        if not info_contrastes:
            st.warning("No se encontraron contrastes con *_logFC / *_adj.P.Val")
        else:
            rows = []
            for pref, meta in info_contrastes.items():
                rows.append({
                    "ID": pref,
                    "Nombre": pretty_contrast_name(pref),
                    "logFC": meta["logFC"],
                    "adj.P.Val": meta["adjP"] or meta["P"]
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True)

# ----------------------------
# Tab 2: Volcano / MA
# ----------------------------
with tabs[1]:
    st.header("🌋 Volcano / MA plot")
    if df is None or not info_contrastes or gene_col is None:
        st.info("Carga datos válidos para ver volcano/MA.")
    else:
        contraste = st.selectbox("Contraste", sorted(info_contrastes.keys()), format_func=pretty_contrast_name)
        deg = preparar_deg_por_contraste(df, gene_col, contraste, info_contrastes)

        col1, col2 = st.columns(2)
        with col1:
            thr_logfc = st.number_input("|logFC| mínimo", 0.0, 10.0, 1.0, 0.1)
        with col2:
            thr_padj = st.number_input("adj.P.Val máximo", 0.0, 1.0, 0.05, 0.01)

        def cat(row):
            if np.isnan(row["logFC"]) or np.isnan(row["padj"]):
                return "No eval"
            if abs(row["logFC"]) >= thr_logfc and row["padj"] <= thr_padj:
                return "Up" if row["logFC"] > 0 else "Down"
            return "No sig"

        deg["cat"] = deg.apply(cat, axis=1)
        volcano = alt.Chart(deg).mark_circle(size=50, opacity=0.7).encode(
            x=alt.X("logFC", title="log2 Fold Change"),
            y=alt.Y("neg_log10_padj", title="-log10(adj.P.Val)"),
            color=alt.Color("cat:N", title="Categoría"),
            tooltip=["gene", "logFC", "padj"]
        ).properties(height=400, title=f"Volcano – {pretty_contrast_name(contraste)}")
        st.altair_chart(volcano, use_container_width=True)

        ma = alt.Chart(deg).mark_circle(size=50, opacity=0.7).encode(
            x=alt.X("AveExpr", title="AveExpr"),
            y=alt.Y("logFC", title="log2 Fold Change"),
            color=alt.Color("cat:N"),
            tooltip=["gene", "AveExpr", "logFC", "padj"]
        ).properties(height=400, title=f"MA plot – {pretty_contrast_name(contraste)}")
        st.altair_chart(ma, use_container_width=True)

# ----------------------------
# Tab 3: DEGs
# ----------------------------
with tabs[2]:
    st.header("📑 Tabla de genes diferencialmente expresados")
    if df is None or not info_contrastes:
        st.info("Sube un archivo con contrastes detectados.")
    else:
        contraste = st.selectbox("Contraste", sorted(info_contrastes.keys()), format_func=pretty_contrast_name, key="contraste_tabla")
        deg = preparar_deg_por_contraste(df, gene_col, contraste, info_contrastes)
        thr_logfc = st.number_input("|logFC| mínimo", 0.0, 10.0, 1.0, 0.1, key="thr1")
        thr_padj = st.number_input("adj.P.Val máximo", 0.0, 1.0, 0.05, 0.01, key="thr2")
        mask = (deg["logFC"].abs() >= thr_logfc) & (deg["padj"] <= thr_padj)
        deg_filt = deg.loc[mask].sort_values("padj")
        st.write(f"Genes significativos: {deg_filt.shape[0]}")
        st.dataframe(deg_filt, use_container_width=True)
        out = io.StringIO(); deg_filt.to_csv(out, index=False)
        st.markdown(bytes_download_link(out.getvalue().encode(), "DEGs_filtrados.csv", "⬇️ Descargar CSV"), unsafe_allow_html=True)

# ----------------------------
# Tab 4: Explorador por gen
# ----------------------------
with tabs[3]:
    st.header("🔍 Explorador por gen")
    if df is None or not info_contrastes or gene_col is None:
        st.info("Carga un dataset con contrastes detectados.")
    else:
        genes = sorted(df[gene_col].dropna().astype(str).unique().tolist())
        gene = st.selectbox("Gen", genes)
        sub = df[df[gene_col].astype(str) == gene].copy()
        rows = []
        for pref, meta in info_contrastes.items():
            logfc_col = meta["logFC"]
            padj_col = meta["adjP"] if meta["adjP"] else meta["P"]
            ave_col = meta["AveExpr"]
            rows.append({
                "Contraste": pretty_contrast_name(pref),
                "logFC": sub[logfc_col].iloc[0] if logfc_col in sub.columns else np.nan,
                "AveExpr": sub[ave_col].iloc[0] if ave_col and ave_col in sub.columns else np.nan,
                "adj.P.Val": sub[padj_col].iloc[0] if padj_col in sub.columns else np.nan
            })
        res = pd.DataFrame(rows)
        st.dataframe(res, use_container_width=True)
        chart = alt.Chart(res).mark_bar().encode(
            x=alt.X("Contraste:N", sort="-y"),
            y=alt.Y("logFC:Q"),
            tooltip=["Contraste", "logFC", "adj.P.Val"]
        ).properties(height=400, title=f"logFC por contraste – {gene}")
        st.altair_chart(chart, use_container_width=True)
