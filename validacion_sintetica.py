# -*- coding: utf-8 -*-
"""
Validación del método: genera un grafo Barabási–Albert (gamma teórico = 3)
en el MISMO formato que CAIDA y ejecuta el análisis completo sobre él.
Si el script estima gamma ≈ 3, el cálculo de Delta es correcto.

Uso:  python validacion_sintetica.py
"""
import bz2
import os
import random
import subprocess
import sys

import networkx as nx

N, M, SEMILLA = 20000, 3, 1
os.makedirs("data", exist_ok=True)
ruta = os.path.join("data", "sintetico_BA.as-rel.txt.bz2")

random.seed(SEMILLA)
G = nx.barabasi_albert_graph(N, M, seed=SEMILLA)
with bz2.open(ruta, "wt") as f:
    f.write("# Grafo sintético Barabási–Albert en formato CAIDA: AS1|AS2|rel\n")
    for u, v in G.edges():
        f.write(f"{u}|{v}|{random.choice([-1, -1, 0])}\n")

print(f"[OK] Grafo BA (N={N}, m={M}) escrito en {ruta}. Gamma teórico = 3.0\n")
subprocess.run([sys.executable, "analisis_caida_delta.py", "--archivo", ruta,
                "--seccion", "validacion"], check=True)
