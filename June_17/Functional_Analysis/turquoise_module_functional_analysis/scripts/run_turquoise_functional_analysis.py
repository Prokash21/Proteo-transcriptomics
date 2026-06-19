from pathlib import Path

import pandas as pd
import requests
import seaborn as sns
from matplotlib import pyplot as plt


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
INPUT_FILE = ANALYSIS_DIR / "inputs" / "turquoise_module_genes.csv"
PROCESSED_DIR = ANALYSIS_DIR / "results" / "processed"
RAW_DIR = ANALYSIS_DIR / "results" / "raw"
PNG_DIR = ANALYSIS_DIR / "figures" / "png"
PDF_DIR = ANALYSIS_DIR / "figures" / "pdf"

for directory in [PROCESSED_DIR, RAW_DIR, PNG_DIR, PDF_DIR]:
    directory.mkdir(parents=True, exist_ok=True)


def split_gene_field(value: str) -> list[str]:
    return [item.strip() for item in str(value).split(";") if item.strip()]


genes = (
    pd.read_csv(INPUT_FILE)["gene"]
    .dropna()
    .map(split_gene_field)
    .explode()
    .dropna()
    .drop_duplicates()
    .sort_values()
    .reset_index(drop=True)
)
genes.to_frame("gene").to_csv(PROCESSED_DIR / "turquoise_module_clean_genes_from_python.csv", index=False)

payload = {
    "organism": "hsapiens",
    "query": genes.tolist(),
    "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"],
    "user_threshold": 0.05,
    "all_results": False,
    "ordered": False,
    "no_evidences": True,
    "domain_scope": "annotated",
    "significance_threshold_method": "g_SCS",
}
response = requests.post(
    "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
    json=payload,
    timeout=120,
)
response.raise_for_status()
data = response.json()

results = pd.DataFrame(data["result"])
results["negative_log10_p"] = -results["p_value"].clip(lower=1e-300).map(lambda value: __import__("math").log10(value))
results["intersection_genes"] = results.get("intersections", pd.Series([[]] * len(results))).map(
    lambda value: ";".join(value) if isinstance(value, list) else ""
)
results = results.sort_values("p_value")
results.to_csv(PROCESSED_DIR / "gprofiler_enrichment_all_results_from_python.csv", index=False)
results.query("significant == True and p_value < 0.05").to_csv(
    PROCESSED_DIR / "gprofiler_enrichment_significant_from_python.csv",
    index=False,
)

top_terms = results.query("significant == True and p_value < 0.05").head(25).copy()
top_terms["label"] = "[" + top_terms["source"] + "] " + top_terms["name"]

plt.figure(figsize=(11, 8))
sns.barplot(data=top_terms, x="negative_log10_p", y="label", color="#4C78A8")
plt.title("Turquoise module: top enriched functional terms")
plt.xlabel("-log10(g:SCS adjusted p-value)")
plt.ylabel("")
plt.tight_layout()
plt.savefig(PNG_DIR / "top25_enriched_terms_barplot_from_python.png", dpi=300)
plt.savefig(PDF_DIR / "top25_enriched_terms_barplot_from_python.pdf")
