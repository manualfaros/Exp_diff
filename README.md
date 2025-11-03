# 🧬 EDA de Expresión Diferencial (Streamlit)

Esta aplicación permite explorar resultados de **análisis de expresión diferencial** (por ejemplo, de edgeR o limma) de forma interactiva.

Está basada en un flujo clásico de **EDA (Exploratory Data Analysis)** y permite visualizar **Volcano plots**, **MA plots**, tablas de **DEGs** y exploración por gen.

---

## 🚀 Funcionalidades principales

- **Carga de tablas** (`.tsv`, `.csv`, `.xlsx`) con resultados de expresión diferencial.  
- **Detección automática de contrastes** a partir de columnas tipo:
- **Visualizaciones interactivas (Altair):**
- Volcano plot (`logFC` vs `-log10(adj.P.Val)`)
- MA plot (`AveExpr` vs `logFC`)
- **Filtro de DEGs** por umbral de |logFC| y FDR (`adj.P.Val`)
- **Explorador por gen:** muestra logFC, AveExpr y FDR de un gen específico en todos los contrastes.
- **Descarga directa de DEGs filtrados** (`⬇️ Descargar CSV`).

---

## 🛠️ Instalación

### 1. Clonar el repositorio

```bash
git clone https://github.com/<tu_usuario>/<tu_repo>.git
cd <tu_repo>
