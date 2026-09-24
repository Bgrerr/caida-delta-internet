# -*- coding: utf-8 -*-
"""
=====================================================================
 ANÁLISIS DE LA TOPOLOGÍA DE INTERNET (NIVEL AS) CON DATOS DE CAIDA
 Cálculo de la "constante universal" Δ (γ): exponente de la
 distribución de grado P(k) ~ k^(-γ) en una red libre de escala.
=====================================================================

Además de imprimir el análisis en consola, el script:
  * guarda figuras individuales en docs/<seccion>/
  * REESCRIBE AUTOMÁTICAMENTE la sección de resultados del README.md
    (entre los marcadores <!-- INICIO:RESULTADOS --> y <!-- FIN:RESULTADOS -->)
    con tablas, figuras e interpretación, para que GitHub muestre todo.

DATASET: CAIDA AS Relationships (serial-1)
  https://publicdata.caida.org/datasets/as-relationships/serial-1/AAAAMMDD.as-rel.txt.bz2
  Formato:  AS1|AS2|-1  (AS1 proveedor de AS2)   |   AS1|AS2|0  (pares)

USO:
  python analisis_caida_delta.py                         # descarga 20240101
  python analisis_caida_delta.py --fecha 20260801
  python analisis_caida_delta.py --archivo data/20240101.as-rel.txt.bz2
  python analisis_caida_delta.py --sin-hiperbolicidad    # más rápido
  python analisis_caida_delta.py --no-readme             # no tocar README.md
"""

import argparse
import bz2
import os
import random
import re
import time
import urllib.request
import warnings
from collections import Counter
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
try:
    import powerlaw
    HAY_POWERLAW = True
except ImportError:
    HAY_POWERLAW = False
    print("[AVISO] 'powerlaw' no instalado: se omitirá el ajuste MLE (pip install powerlaw)")

URL_BASE = "https://publicdata.caida.org/datasets/as-relationships/serial-1/"
DIR_DATOS = "data"
DIR_RESULTADOS = "resultados"
DIR_DOCS = "docs"
README = "README.md"
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 140, "font.size": 11})


# ==================================================================
# 1. CARGA DEL DATASET
# ==================================================================
def descargar_dataset(fecha: str) -> str:
    os.makedirs(DIR_DATOS, exist_ok=True)
    archivo = f"{fecha}.as-rel.txt.bz2"
    ruta = os.path.join(DIR_DATOS, archivo)
    if os.path.exists(ruta):
        print(f"[OK] Usando archivo local existente: {ruta}")
        return ruta
    url = URL_BASE + archivo
    print(f"[..] Descargando {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (analisis academico)"})
    with urllib.request.urlopen(req, timeout=180) as r, open(ruta, "wb") as f:
        f.write(r.read())
    print(f"[OK] Guardado como {ruta}")
    return ruta


def cargar_relaciones(ruta: str):
    """Devuelve (DataFrame de enlaces, lista de comentarios de cabecera)."""
    abrir = bz2.open if ruta.endswith(".bz2") else open
    filas, cabecera = [], []
    with abrir(ruta, "rt", encoding="utf-8", errors="ignore") as f:
        for linea in f:
            if linea.startswith("#"):
                cabecera.append(linea.strip("# \n"))
                continue
            partes = linea.strip().split("|")
            if len(partes) >= 3:
                filas.append((int(partes[0]), int(partes[1]), int(partes[2])))
    df = pd.DataFrame(filas, columns=["as1", "as2", "rel"])
    return df, cabecera


# ==================================================================
# 2. GRAFO Y ANÁLISIS ESTRUCTURAL
# ==================================================================
def construir_grafo(df):
    G = nx.Graph()
    G.add_edges_from(zip(df["as1"], df["as2"]))
    G.remove_edges_from(nx.selfloop_edges(G))
    return G


def analisis_basico(G, df, muestras_camino=300):
    grados = np.array([d for _, d in G.degree()], dtype=float)
    N, E = G.number_of_nodes(), G.number_of_edges()
    k_med, k2_med = 2 * E / N, np.mean(grados ** 2)
    comp = max(nx.connected_components(G), key=len)
    Gc = G.subgraph(comp)

    C = nx.algorithms.approximation.average_clustering(Gc, trials=20000, seed=42)

    rng = random.Random(42)
    muestra = rng.sample(list(Gc.nodes()), min(muestras_camino, len(comp)))
    dist_total, pares, diam, hist_dist = 0, 0, 0, Counter()
    for s in muestra:
        L = nx.single_source_shortest_path_length(Gc, s)
        dist_total += sum(L.values()); pares += len(L) - 1
        diam = max(diam, max(L.values()))
        hist_dist.update(v for v in L.values() if v > 0)

    kappa = k2_med / k_med
    st = {
        "N": N, "E": E,
        "E_p2c": int((df["rel"] == -1).sum()), "E_p2p": int((df["rel"] == 0).sum()),
        "k_med": k_med, "k2_med": k2_med, "k_max": int(grados.max()), "k_min": int(grados.min()),
        "k_mediana": float(np.median(grados)),
        "densidad": nx.density(G), "frac_gigante": len(comp) / N,
        "n_componentes": nx.number_connected_components(G),
        "clustering": C, "l_med": dist_total / pares, "diametro": diam,
        "asortatividad": nx.degree_assortativity_coefficient(G),
        "kappa": kappa, "fc_aleatorio": 1 - 1 / (kappa - 1) if kappa > 1 else float("nan"),
        "hist_dist": hist_dist,
        "ln_N": np.log(N), "lnln_N": np.log(np.log(N)),
    }
    print("\n" + "=" * 62 + "\n ANÁLISIS ESTRUCTURAL\n" + "=" * 62)
    for k, v in st.items():
        if k != "hist_dist":
            print(f"  {k:<16}: {v:,.4f}" if isinstance(v, float) else f"  {k:<16}: {v:,}")
    return st


def clases_de_grado(G):
    grados = np.array([d for _, d in G.degree()])
    bins = [(1, 1, "1 (stub puro)"), (2, 2, "2"), (3, 10, "3 – 10"), (11, 100, "11 – 100"),
            (101, 1000, "101 – 1 000"), (1001, 10**9, "> 1 000 (hubs)")]
    filas = []
    for lo, hi, et in bins:
        m = (grados >= lo) & (grados <= hi)
        filas.append({"clase": et, "n_AS": int(m.sum()), "pct_AS": 100 * m.mean(),
                      "pct_enlaces": 100 * grados[m].sum() / grados.sum()})
    return pd.DataFrame(filas)


def top_as(G, df, n=15):
    clientes = df[df.rel == -1].groupby("as1").size()
    proveedores = df[df.rel == -1].groupby("as2").size()
    pares = pd.concat([df[df.rel == 0].as1, df[df.rel == 0].as2]).value_counts()
    top = sorted(G.degree(), key=lambda x: -x[1])[:n]
    N = G.number_of_nodes()
    return pd.DataFrame([{
        "rango": i + 1, "ASN": a, "grado": d, "pct_nodos": 100 * d / (N - 1),
        "clientes": int(clientes.get(a, 0)), "proveedores": int(proveedores.get(a, 0)),
        "pares": int(pares.get(a, 0))} for i, (a, d) in enumerate(top)])


# ==================================================================
# 3. CONSTANTE UNIVERSAL Δ (γ)
# ==================================================================
def binning_logaritmico(grados, n_bins=30):
    bordes = np.unique(np.logspace(np.log10(grados.min()), np.log10(grados.max() + 1), n_bins).astype(int))
    hist, bordes = np.histogram(grados, bins=bordes)
    pdf = hist / (np.diff(bordes) * len(grados))
    centros = np.sqrt(bordes[:-1] * bordes[1:])
    m = pdf > 0
    return centros[m], pdf[m]


def calcular_delta(G):
    grados = np.array([d for _, d in G.degree() if d > 0])
    R = {}
    k_b, p_b = binning_logaritmico(grados)
    lr = stats.linregress(np.log10(k_b), np.log10(p_b))
    R["A"] = dict(gamma=-lr.slope, b=lr.intercept, R2=lr.rvalue ** 2, x=k_b, y=p_b,
                  err=lr.stderr)

    k_ord = np.sort(np.unique(grados))
    ccdf = np.array([(grados >= k).mean() for k in k_ord])
    m = (k_ord >= 2) & (k_ord <= np.percentile(grados, 99.9))
    lr2 = stats.linregress(np.log10(k_ord[m]), np.log10(ccdf[m]))
    R["B"] = dict(gamma=1 - lr2.slope, pend=lr2.slope, b=lr2.intercept, R2=lr2.rvalue ** 2,
                  x=k_ord, y=ccdf, err=lr2.stderr)

    rk = np.sort(grados)[::-1]; rg = np.arange(1, len(rk) + 1)
    lr3 = stats.linregress(np.log10(rg), np.log10(rk))
    R["C"] = dict(R=lr3.slope, gamma=1 - 1 / lr3.slope, b=lr3.intercept, R2=lr3.rvalue ** 2,
                  x=rg, y=rk)

    if HAY_POWERLAW:
        fit = powerlaw.Fit(grados, discrete=True, verbose=False)
        comps = {}
        for alt in ["exponential", "lognormal", "truncated_power_law", "stretched_exponential"]:
            try:
                comps[alt] = fit.distribution_compare("power_law", alt, normalized_ratio=True)
            except Exception:
                comps[alt] = (float("nan"), float("nan"))
        R["D"] = dict(gamma=fit.power_law.alpha, sigma=fit.power_law.sigma, kmin=fit.power_law.xmin,
                      KS=fit.power_law.D, fit=fit, comps=comps,
                      frac_cola=(grados >= fit.power_law.xmin).mean(),
                      n_cola=int((grados >= fit.power_law.xmin).sum()))
    R["final"] = R["D"]["gamma"] if "D" in R else R["A"]["gamma"]

    print("\n" + "=" * 62 + "\n CONSTANTE UNIVERSAL Δ (γ)\n" + "=" * 62)
    print(f"  [A] PDF log-log      γ = {R['A']['gamma']:.3f}  (R²={R['A']['R2']:.3f})")
    print(f"  [B] CCDF             γ = {R['B']['gamma']:.3f}  (R²={R['B']['R2']:.3f})")
    print(f"  [C] Rango Faloutsos  R = {R['C']['R']:.3f} → γ = {R['C']['gamma']:.3f}")
    if "D" in R:
        d = R["D"]
        print(f"  [D] MLE              γ = {d['gamma']:.3f} ± {d['sigma']:.3f}  (k_min={d['kmin']:.0f})")
        for alt, (r, p) in d["comps"].items():
            print(f"      power-law vs {alt:<22}: R={r:+.2f}  p={p:.3g}")
    print(f"\n  >> Δ = γ ≈ {R['final']:.2f}")
    return R


# ==================================================================
# 4. ROBUSTEZ Y 5. HIPERBOLICIDAD DE GROMOV
# ==================================================================
def robustez(G, fracciones=(0, .001, .002, .005, .01, .02, .03, .05, .1, .2, .3, .5), seed=42):
    N = G.number_of_nodes()
    por_grado = [n for n, _ in sorted(G.degree(), key=lambda x: -x[1])]
    aleatorio = list(G.nodes()); random.Random(seed).shuffle(aleatorio)

    def S(orden, f):
        H = G.copy(); H.remove_nodes_from(orden[:int(f * N)])
        return len(max(nx.connected_components(H), key=len)) / N if H.number_of_nodes() else 0

    df = pd.DataFrame({"f": fracciones,
                       "S_aleatorio": [S(aleatorio, f) for f in fracciones],
                       "S_ataque": [S(por_grado, f) for f in fracciones]})
    print("\n" + "=" * 62 + "\n ROBUSTEZ\n" + "=" * 62 + "\n" + df.to_string(index=False))
    return df


def hiperbolicidad_gromov(G, n_fuentes=120, n_cuadruplas=200000, seed=7):
    Gc = G.subgraph(max(nx.connected_components(G), key=len))
    rng = random.Random(seed)
    fuentes = rng.sample(list(Gc.nodes()), min(n_fuentes, Gc.number_of_nodes()))
    D = {s: nx.single_source_shortest_path_length(Gc, s) for s in fuentes}
    deltas = []
    for _ in range(n_cuadruplas):
        x, y, u, v = rng.sample(fuentes, 4)
        s = sorted([D[x][y] + D[u][v], D[x][u] + D[y][v], D[x][v] + D[y][u]], reverse=True)
        deltas.append((s[0] - s[1]) / 2)
    deltas = np.array(deltas)
    print("\n" + "=" * 62 + f"\n δ DE GROMOV: máx = {deltas.max():.1f}, media = {deltas.mean():.3f}\n" + "=" * 62)
    return deltas


# ==================================================================
# 6. FIGURAS
# ==================================================================
def _guardar(fig, carpeta, nombre, figs):
    ruta = os.path.join(carpeta, nombre)
    fig.tight_layout(); fig.savefig(ruta); plt.close(fig)
    figs[nombre.split(".")[0]] = ruta
    return ruta


def generar_figuras(st, R, rob, clases, deltas, etiqueta, carpeta):
    os.makedirs(carpeta, exist_ok=True)
    figs = {}

    # 1. Composición del dataset
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ax[0].bar(["Proveedor → cliente", "Par ↔ par"], [st["E_p2c"], st["E_p2p"]],
              color=["#1f77b4", "#ff7f0e"])
    for i, v in enumerate([st["E_p2c"], st["E_p2p"]]):
        ax[0].text(i, v, f"{v:,}\n({100*v/st['E']:.1f}%)", ha="center", va="bottom")
    ax[0].set(title="Tipos de relación entre AS", ylabel="Número de enlaces")
    ax[0].set_ylim(0, max(st["E_p2c"], st["E_p2p"]) * 1.2)
    ax[1].bar(clases["clase"], clases["pct_AS"], color="#2ca02c", label="% de AS")
    ax[1].plot(clases["clase"], clases["pct_enlaces"], "ro-", label="% de extremos de enlace")
    ax[1].set(title="Heterogeneidad: pocos AS concentran los enlaces", ylabel="%", xlabel="Grado k")
    ax[1].tick_params(axis="x", rotation=25); ax[1].legend()
    _guardar(fig, carpeta, "01_dataset.png", figs)

    # 2. PDF log-log
    a = R["A"]
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.loglog(a["x"], a["y"], "o", ms=7, label="P(k) empírica (binning log)")
    xx = np.logspace(np.log10(a["x"].min()), np.log10(a["x"].max()), 100)
    ax.loglog(xx, 10 ** a["b"] * xx ** (-a["gamma"]), "r-", lw=2,
              label=f"Regresión: pendiente = −{a['gamma']:.2f} (R² = {a['R2']:.3f})")
    ax.set(xlabel="Grado k", ylabel="P(k)", title="Distribución de grado en escala log-log")
    ax.legend(); ax.grid(True, which="both", alpha=.3)
    _guardar(fig, carpeta, "02_pdf_loglog.png", figs)

    # 3. CCDF + MLE
    b = R["B"]
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.loglog(b["x"], b["y"], ".", color="gray", label="CCDF empírica P(K ≥ k)")
    xx = np.logspace(np.log10(2), np.log10(b["x"].max()), 100)
    ax.loglog(xx, 10 ** b["b"] * xx ** b["pend"], "b--", lw=2,
              label=f"OLS: pendiente = {b['pend']:.2f} → γ = {b['gamma']:.2f}")
    if "D" in R:
        d = R["D"]
        kk = np.logspace(np.log10(d["kmin"]), np.log10(b["x"].max()), 100)
        ax.loglog(kk, d["frac_cola"] * (kk / d["kmin"]) ** (1 - d["gamma"]), "r-", lw=2.5,
                  label=f"MLE: γ = {d['gamma']:.2f} ± {d['sigma']:.2f}")
        ax.axvline(d["kmin"], color="r", ls=":", label=f"k_min = {d['kmin']:.0f}")
    ax.set(xlabel="Grado k", ylabel="P(K ≥ k)", title="Función de distribución complementaria (CCDF)")
    ax.legend(); ax.grid(True, which="both", alpha=.3)
    _guardar(fig, carpeta, "03_ccdf_mle.png", figs)

    # 4. Rango
    c = R["C"]
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.loglog(c["x"], c["y"], ".", ms=3, color="green", label="Grado vs. rango")
    ax.loglog(c["x"], 10 ** c["b"] * c["x"] ** c["R"], "k-",
              label=f"R = {c['R']:.2f} → γ = {c['gamma']:.2f}")
    ax.set(xlabel="Rango r del AS", ylabel="Grado d_r", title="Ley de rango (Faloutsos et al., 1999)")
    ax.legend(); ax.grid(True, which="both", alpha=.3)
    _guardar(fig, carpeta, "04_rango.png", figs)

    # 5. Distancias
    h = st["hist_dist"]; tot = sum(h.values())
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    xs = sorted(h); ax.bar(xs, [100 * h[x] / tot for x in xs], color="#9467bd")
    ax.axvline(st["l_med"], color="k", ls="--", label=f"⟨ℓ⟩ = {st['l_med']:.2f}")
    ax.set(xlabel="Distancia (saltos AS)", ylabel="% de pares", title="Distribución de distancias: mundo pequeño")
    ax.legend()
    _guardar(fig, carpeta, "05_distancias.png", figs)

    # 6. Robustez
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(rob.f, rob.S_aleatorio, "o-", label="Fallos aleatorios")
    ax.plot(rob.f, rob.S_ataque, "s-", color="crimson", label="Ataque dirigido a hubs")
    ax.set_xscale("symlog", linthresh=1e-3)
    ax.set(xlabel="Fracción f de AS eliminados", ylabel="Tamaño relativo de la componente gigante S",
           title="Robusta ante fallos, frágil ante ataques")
    ax.legend(); ax.grid(True, alpha=.3)
    _guardar(fig, carpeta, "06_robustez.png", figs)

    # 7. Gromov
    if deltas is not None:
        cnt = Counter(deltas); xs = sorted(cnt)
        fig, ax = plt.subplots(figsize=(7.5, 4.5))
        ax.bar([str(x) for x in xs], [100 * cnt[x] / len(deltas) for x in xs], color="#8c564b")
        ax.set(xlabel="δ(x, y, u, v)", ylabel="% de cuádruplas", yscale="log",
               title=f"Hiperbolicidad de Gromov (δ máx = {deltas.max():.1f})")
        _guardar(fig, carpeta, "07_gromov.png", figs)

    # 8. Resumen 2x2
    fig, axs = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle(f"Topología AS de Internet ({etiqueta}) — Δ = γ ≈ {R['final']:.2f}",
                 fontsize=15, fontweight="bold")
    axs[0, 0].loglog(a["x"], a["y"], "o"); xx = np.logspace(np.log10(a["x"].min()), np.log10(a["x"].max()), 50)
    axs[0, 0].loglog(xx, 10 ** a["b"] * xx ** (-a["gamma"]), "r-", label=f"pendiente −{a['gamma']:.2f}")
    axs[0, 0].set(title="(a) PDF log-log", xlabel="k", ylabel="P(k)"); axs[0, 0].legend()
    axs[0, 1].loglog(b["x"], b["y"], ".", color="gray")
    if "D" in R:
        axs[0, 1].loglog(kk, d["frac_cola"] * (kk / d["kmin"]) ** (1 - d["gamma"]), "r-",
                         label=f"MLE γ = {d['gamma']:.2f}")
        axs[0, 1].legend()
    axs[0, 1].set(title="(b) CCDF + MLE", xlabel="k", ylabel="P(K ≥ k)")
    axs[1, 0].loglog(c["x"], c["y"], ".", ms=3, color="green")
    axs[1, 0].set(title=f"(c) Ley de rango, R = {c['R']:.2f}", xlabel="rango", ylabel="grado")
    axs[1, 1].plot(rob.f, rob.S_aleatorio, "o-", label="aleatorio")
    axs[1, 1].plot(rob.f, rob.S_ataque, "s-", color="crimson", label="ataque")
    axs[1, 1].set_xscale("symlog", linthresh=1e-3); axs[1, 1].legend()
    axs[1, 1].set(title="(d) Robustez", xlabel="f", ylabel="S")
    for axx in axs.flat: axx.grid(True, which="both", alpha=.3)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _guardar(fig, carpeta, "00_resumen.png", figs)
    print(f"\n[OK] {len(figs)} figuras guardadas en {carpeta}/")
    return figs


# ==================================================================
# 7. GENERACIÓN AUTOMÁTICA DEL README
# ==================================================================
def _md_tabla(df, formatos=None):
    formatos = formatos or {}
    cab = "| " + " | ".join(df.columns) + " |\n|" + "|".join(["---"] * len(df.columns)) + "|\n"
    filas = []
    for _, r in df.iterrows():
        celdas = [formatos[c](r[c]) if c in formatos else str(r[c]) for c in df.columns]
        filas.append("| " + " | ".join(celdas) + " |")
    return cab + "\n".join(filas)


def construir_markdown(st, R, rob, clases, top, deltas, figs, etiqueta, ruta_datos, cabecera, t_ejec, seccion):
    rel = lambda p: p.replace(os.sep, "/")
    g = R["final"]
    fuente = next((c for c in cabecera if "source" in c.lower()), "")
    if 2 < g < 3:
        regimen = ("**2 < γ < 3 → régimen libre de escala.** ⟨k⟩ es finito pero ⟨k²⟩ diverge al crecer N; "
                   "la red es un *mundo ultra-pequeño* y extremadamente robusta ante fallos aleatorios.")
    elif g <= 2:
        regimen = ("**γ ≤ 2 → régimen anómalo.** Incluso ⟨k⟩ diverge con N; los hubs se conectan con una "
                   "fracción finita de todos los nodos (estructura tipo estrella). Valores ≈ 2 son habituales "
                   "en snapshots recientes de CAIDA por el crecimiento del *peering* de grandes redes.")
    else:
        regimen = ("**γ > 3 → heterogeneidad moderada.** ⟨k²⟩ es finito; la red se comporta de forma más "
                   "parecida a una red aleatoria en robustez y distancias.")

    top_md = top.copy()
    top_md["ASN"] = top_md["ASN"].map(lambda a: f"[AS{a}](https://asrank.caida.org/asns/{a})")
    rob_fc_ataque = next((f for f, s in zip(rob.f, rob.S_ataque) if s < 0.05), None)

    L = []
    titulo = "Resultados con datos de CAIDA" if seccion == "resultados" else "Validación con datos sintéticos"
    L.append(f"### {titulo}\n")
    L.append(f"> Sección generada automáticamente por `analisis_caida_delta.py` el "
             f"{datetime.now():%Y-%m-%d %H:%M} (tiempo de ejecución: {t_ejec/60:.1f} min).\n")
    L.append(f"**Archivo analizado:** `{os.path.basename(ruta_datos)}`  ")
    if fuente:
        L.append(f"**Fuente declarada en el archivo:** `{fuente}`  ")
    L.append(f"**Snapshot:** `{etiqueta}`\n")

    L.append(f"#### Resumen\n\n![Resumen]({rel(figs['00_resumen'])})\n")
    L.append("| Resultado principal | Valor |\n|---|---|")
    L.append(f"| **Constante universal Δ = γ (MLE)** | **{g:.3f}**" +
             (f" ± {R['D']['sigma']:.3f}" if "D" in R else "") + " |")
    L.append(f"| Sistemas Autónomos (nodos) | {st['N']:,} |")
    L.append(f"| Enlaces | {st['E']:,} |")
    L.append(f"| Grado promedio ⟨k⟩ | {st['k_med']:.2f} |")
    L.append(f"| Grado máximo | {st['k_max']:,} |\n")
    L.append(regimen + "\n")

    L.append("#### 1. Estadísticas del dataset\n")
    L.append(f"![Composición del dataset]({rel(figs['01_dataset'])})\n")
    L.append("| Métrica | Valor | Comentario |\n|---|---|---|")
    filas = [
        ("Nodos N (AS)", f"{st['N']:,}", "Redes con política de enrutamiento propia"),
        ("Enlaces E", f"{st['E']:,}", "Relaciones BGP inferidas"),
        ("Enlaces proveedor → cliente", f"{st['E_p2c']:,} ({100*st['E_p2c']/max(st['E'],1):.1f} %)", "Tránsito pagado (rel = −1)"),
        ("Enlaces par ↔ par", f"{st['E_p2p']:,} ({100*st['E_p2p']/max(st['E'],1):.1f} %)", "Peering, típicamente en IXPs (rel = 0)"),
        ("Grado promedio ⟨k⟩", f"{st['k_med']:.3f}", "2E / N"),
        ("Grado mediano", f"{st['k_mediana']:.0f}", "Mucho menor que ⟨k⟩: distribución sesgada"),
        ("Segundo momento ⟨k²⟩", f"{st['k2_med']:,.1f}", "Enorme frente a ⟨k⟩² = " + f"{st['k_med']**2:,.1f}"),
        ("Grado mínimo / máximo", f"{st['k_min']} / {st['k_max']:,}", f"El mayor hub toca el {100*st['k_max']/(st['N']-1):.1f} % de la red"),
        ("Densidad", f"{st['densidad']:.2e}", "Red muy dispersa"),
        ("Componentes conexas", f"{st['n_componentes']:,}", ""),
        ("Fracción en componente gigante", f"{100*st['frac_gigante']:.2f} %", ""),
        ("Clustering medio C", f"{st['clustering']:.4f}", "Aproximado (20 000 muestras)"),
        ("Longitud media de camino ⟨ℓ⟩", f"{st['l_med']:.3f}", f"ln N = {st['ln_N']:.2f}; ln ln N = {st['lnln_N']:.2f}"),
        ("Diámetro (cota inferior)", f"{st['diametro']}", "Por muestreo BFS"),
        ("Asortatividad r", f"{st['asortatividad']:.4f}", "r < 0: disasortativa" if st['asortatividad'] < 0 else "r ≥ 0: asortativa"),
        ("κ = ⟨k²⟩ / ⟨k⟩", f"{st['kappa']:.2f}", "Criterio de Molloy-Reed (> 2: hay componente gigante)"),
    ]
    L += [f"| {a} | {b} | {c} |" for a, b, c in filas]
    L.append("\n**Distribución de los AS por clase de grado**\n")
    L.append(_md_tabla(clases.rename(columns={"clase": "Grado k", "n_AS": "N.º de AS",
                                              "pct_AS": "% de AS", "pct_enlaces": "% de extremos de enlace"}),
                       {"N.º de AS": lambda v: f"{int(v):,}", "% de AS": lambda v: f"{v:.2f} %",
                        "% de extremos de enlace": lambda v: f"{v:.2f} %"}))
    L.append(f"\n**Los {len(top)} AS con mayor grado** (enlace a su ficha en CAIDA AS Rank)\n")
    L.append(_md_tabla(top_md.rename(columns={"rango": "#", "grado": "Grado", "pct_nodos": "% de la red",
                                              "clientes": "Clientes", "proveedores": "Proveedores",
                                              "pares": "Pares"}),
                       {"Grado": lambda v: f"{int(v):,}", "% de la red": lambda v: f"{v:.2f} %",
                        "Clientes": lambda v: f"{int(v):,}", "Proveedores": lambda v: f"{int(v):,}",
                        "Pares": lambda v: f"{int(v):,}"}))
    L.append(f"\n![Distancias]({rel(figs['05_distancias'])})\n")

    L.append("#### 2. Cálculo de la constante universal Δ (γ)\n")
    L.append(f"![PDF log-log]({rel(figs['02_pdf_loglog'])})\n")
    L.append(f"La recta de la figura es $\\log P(k) = -\\gamma\\,\\log k + c$. Su pendiente es "
             f"**−{R['A']['gamma']:.3f}**, así que Δ = γ ≈ {R['A']['gamma']:.2f} por este método.\n")
    L.append(f"![CCDF y MLE]({rel(figs['03_ccdf_mle'])})\n")
    L.append(f"![Ley de rango]({rel(figs['04_rango'])})\n")
    L.append("| Método | Pendiente medida | γ estimado | Calidad |\n|---|---|---|---|")
    L.append(f"| A. Regresión PDF log-log | {-R['A']['gamma']:.3f} | {R['A']['gamma']:.3f} | R² = {R['A']['R2']:.4f} |")
    L.append(f"| B. Regresión CCDF | {R['B']['pend']:.3f} | {R['B']['gamma']:.3f} | R² = {R['B']['R2']:.4f} |")
    L.append(f"| C. Ley de rango (Faloutsos) | R = {R['C']['R']:.3f} | {R['C']['gamma']:.3f} | R² = {R['C']['R2']:.4f} |")
    if "D" in R:
        d = R["D"]
        L.append(f"| **D. Máxima verosimilitud (CSN 2009)** | — | **{d['gamma']:.3f} ± {d['sigma']:.3f}** | "
                 f"KS = {d['KS']:.4f}; k_min = {d['kmin']:.0f}; cola = {d['n_cola']:,} AS ({100*d['frac_cola']:.1f} %) |")
        nombres = {"exponential": "Exponencial", "lognormal": "Lognormal",
                   "truncated_power_law": "Ley de potencias truncada", "stretched_exponential": "Exponencial estirada"}
        L.append("\n**Pruebas de razón de verosimilitud** (ley de potencias contra alternativas)\n")
        L.append("| Alternativa | R normalizado | p-valor | Veredicto |\n|---|---|---|---|")
        for alt, (r, p) in d["comps"].items():
            if np.isnan(r):
                ver = "No evaluable"
            elif p >= 0.1:
                ver = "No concluyente (ambas plausibles)"
            else:
                ver = "Favorece ley de potencias" if r > 0 else f"Favorece {nombres[alt].lower()}"
            L.append(f"| {nombres[alt]} | {r:+.3f} | {p:.3g} | {ver} |")
        L.append("\nSi R > 0 con p < 0.1, la ley de potencias ajusta mejor. Es habitual que la ley de potencias "
                 "truncada o la lognormal resulten igual de plausibles en datos reales: la cola de Internet "
                 "decae algo más rápido en los grados más altos por límites físicos y económicos.\n")

    L.append("#### 3. Robustez: consecuencia práctica de Δ\n")
    L.append(f"![Robustez]({rel(figs['06_robustez'])})\n")
    L.append(_md_tabla(rob.rename(columns={"f": "Fracción eliminada f", "S_aleatorio": "S (fallos aleatorios)",
                                           "S_ataque": "S (ataque a hubs)"}),
                       {"Fracción eliminada f": lambda v: f"{100*v:.1f} %",
                        "S (fallos aleatorios)": lambda v: f"{v:.3f}", "S (ataque a hubs)": lambda v: f"{v:.3f}"}))
    L.append(f"\n- Umbral teórico ante fallos aleatorios (Molloy-Reed): $f_c = 1 - 1/(\\kappa-1)$ = "
             f"**{100*st['fc_aleatorio']:.1f} %** de los AS.")
    L.append(f"- Ante ataque dirigido, la componente gigante cae por debajo del 5 % al eliminar "
             f"**{100*rob_fc_ataque:.1f} %** de los AS." if rob_fc_ataque is not None else
             "- Ante ataque dirigido, la componente gigante se reduce drásticamente con pocos AS eliminados.")
    L.append("")

    if deltas is not None:
        cnt = Counter(deltas)
        L.append("#### 4. Hiperbolicidad δ de Gromov\n")
        L.append(f"![Gromov]({rel(figs['07_gromov'])})\n")
        L.append(f"Se muestrearon {len(deltas):,} cuádruplas de nodos de la componente gigante. "
                 f"**δ máximo = {deltas.max():.1f}** y δ medio = {deltas.mean():.3f}. "
                 "Un valor tan pequeño frente al diámetro indica que la red es *similar a un árbol* a gran escala, "
                 "con curvatura negativa.\n")
        L.append("| δ | % de cuádruplas |\n|---|---|")
        L += [f"| {v:.1f} | {100*c/len(deltas):.3f} % |" for v, c in sorted(cnt.items())]
        L.append("")
    return "\n".join(L)


def actualizar_readme(md, seccion):
    ini, fin = f"<!-- INICIO:{seccion.upper()} -->", f"<!-- FIN:{seccion.upper()} -->"
    if not os.path.exists(README):
        print(f"[AVISO] No existe {README}; se guarda la sección en {DIR_RESULTADOS}/seccion_{seccion}.md")
        open(os.path.join(DIR_RESULTADOS, f"seccion_{seccion}.md"), "w", encoding="utf-8").write(md)
        return
    txt = open(README, encoding="utf-8").read()
    # Los marcadores deben ocupar una línea completa (evita coincidir con menciones en el texto)
    patron = re.compile(r"^" + re.escape(ini) + r"$.*?^" + re.escape(fin) + r"$", flags=re.S | re.M)
    if not patron.search(txt):
        txt += f"\n\n{ini}\n{fin}\n"
    nuevo = patron.sub(lambda _: f"{ini}\n\n{md}\n\n{fin}", txt)
    open(README, "w", encoding="utf-8").write(nuevo)
    print(f"[OK] README.md actualizado (sección {seccion})")


# ==================================================================
# MAIN
# ==================================================================
def main():
    p = argparse.ArgumentParser(description="Topología AS de CAIDA y constante universal Δ")
    p.add_argument("--fecha", default="20240101", help="Snapshot AAAAMMDD (día 01 de un mes)")
    p.add_argument("--archivo", default=None, help="Ruta a un .as-rel.txt(.bz2) ya descargado")
    p.add_argument("--sin-hiperbolicidad", action="store_true", help="Omitir δ de Gromov")
    p.add_argument("--no-readme", action="store_true", help="No modificar README.md")
    p.add_argument("--seccion", default="resultados", choices=["resultados", "validacion"],
                   help="Sección del README a reescribir")
    args = p.parse_args()
    t0 = time.time()

    ruta = args.archivo or descargar_dataset(args.fecha)
    etiqueta = os.path.basename(ruta).split(".")[0]
    df, cabecera = cargar_relaciones(ruta)
    print(f"[OK] {len(df):,} relaciones AS leídas de {ruta}")

    G = construir_grafo(df)
    st = analisis_basico(G, df)
    clases = clases_de_grado(G)
    top = top_as(G, df)
    print("\n" + clases.to_string(index=False) + "\n\n" + top.to_string(index=False))
    R = calcular_delta(G)
    rob = robustez(G)
    deltas = None if args.sin_hiperbolicidad else hiperbolicidad_gromov(G)

    figs = generar_figuras(st, R, rob, clases, deltas, etiqueta, os.path.join(DIR_DOCS, args.seccion))

    os.makedirs(DIR_RESULTADOS, exist_ok=True)
    pd.DataFrame(G.degree(), columns=["ASN", "grado"]).sort_values("grado", ascending=False) \
        .to_csv(os.path.join(DIR_RESULTADOS, f"grados_{etiqueta}.csv"), index=False)

    md = construir_markdown(st, R, rob, clases, top, deltas, figs, etiqueta, ruta, cabecera,
                            time.time() - t0, args.seccion)
    open(os.path.join(DIR_RESULTADOS, f"resumen_{etiqueta}.md"), "w", encoding="utf-8").write(md)
    if not args.no_readme:
        actualizar_readme(md, args.seccion)
    print(f"\n[OK] Terminado en {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
