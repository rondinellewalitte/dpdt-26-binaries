# -*- coding: utf-8 -*-
"""Em que direcao a alavanca enviesa o estimador de `qual_barra_26`?

A revisao pede que a nota diga que "o vies da alavanca joga a favor: 2,73 e
piso, 1,03 e teto". Antes de escrever isso, mede-se. O estimador e
Sigma z^2 / Sigma (1 - h) por conjunto, com z e h calculados com as barras
ASSUMIDAS. Quando as barras verdadeiras sao k_T x (TESS) e k_S x (SuperWASP)
as assumidas, o ajuste ponderado usa pesos errados, a matriz-chapeu muda, e
o que cada conjunto devolve nao e necessariamente k^2.

Simulacao nos 26 desenhos reais (E, barras assumidas, fonte de cada epoca):
residuos gerados como ruido gaussiano com barras verdadeiras k x assumidas
em torno de uma parabola exata (o modelo e certo por construcao - so as
barras erram), ajuste quadratico ponderado com as barras assumidas,
estimador por conjunto, 2000 realizacoes por par (k_T, k_S). Grade: k_T em
{1, 1,34, 2}, k_S em {1, 1,65, 2,15}.

EXPECTATIVA (escrita antes de rodar): em (1, 1) os dois devolvem 1
(identidade). Em (1, 1,65) o SuperWASP devolve ~2,7 e o TESS ~1 se nao
houver vazamento. A pergunta e o par (k_T > 1, k_S): se o excesso do TESS
VAZA para o SuperWASP (o ajuste, sobrepesando o TESS, empurra o desajuste
para as temporadas), o SuperWASP recuperado sera MAIOR que k_S^2 - e 2,73
seria teto, nao piso - e o TESS recuperado MENOR que k_T^2 - e 1,03 seria
piso, nao teto. Nao sei a direcao; a frase da nota vai ser a que a
simulacao autorizar.

RODADA NO ESTADO B (2026-09-13): mesmos desenhos, barras SuperWASP
corrigidas. EXPECTATIVA antes de rodar: nada muda de forma - o estimador
com barras FIXAS e imparcial (dentro de 2%) para fatores 1, 1,34 e 2 no
TESS e 1, 1,65 e 2,15 no SuperWASP, sem vazamento; a senoide de 2 a no
TESS levanta o TESS primeiro. O que este script NAO ve - o vies de
co-estimacao com barras re-estimadas - esta em nulo_barras_reestimadas_26
(nulo 0,77), e a nota passa a citar os dois lado a lado.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
from qual_barra_26 import alavancas  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
N_SIM = 2000


def desenhos():
    out = []
    for p in sorted(BASE.glob("oc__*.json")):
        j = json.loads(p.read_text(encoding="utf-8"))
        if j.get("status") != "ok" or "parabola" not in j:
            continue
        pts = pd.DataFrame(j["pontos"])
        out.append((pts.E.values.astype(float), pts.sig_min.values / 1440.0, pts.fonte.str.startswith("TESS").values))
    assert len(out) == 26
    return out


def estimador(des, kT, kS, rng):
    zT, hT, zS, hS = [], [], [], []
    for E, sig, tess in des:
        verdadeira = sig * np.where(tess, kT, kS)
        r = rng.normal(0, verdadeira)                      # residuos em torno de uma parabola exata
        w = 1 / sig
        A = np.vander(E, 3) * w[:, None]
        coef, *_ = np.linalg.lstsq(A, r * w, rcond=None)   # o ajuste "absorve" parte do ruido
        res = r - np.polyval(coef, E)
        z = res / sig
        h = alavancas(E, sig)
        zT += list(z[tess] ** 2); hT += list(1 - h[tess]); zS += list(z[~tess] ** 2); hS += list(1 - h[~tess])
    return sum(zT) / sum(hT), sum(zS) / sum(hS)


if __name__ == "__main__":
    des = desenhos()
    rng = np.random.default_rng(26)
    print(f"{'k_T':>5s} {'k_S':>5s} | {'verdade T':>10s} {'recup. T':>9s} {'IC 95%':>13s} | {'verdade S':>10s} {'recup. S':>9s} {'IC 95%':>13s}")
    linhas = []
    for kT in (1.0, 1.34, 2.0):
        for kS in (1.0, 1.65, 2.15):
            fT, fS = np.array([estimador(des, kT, kS, rng) for _ in range(N_SIM)]).T
            linhas.append({"kT": kT, "kS": kS, "verdade_T": kT ** 2, "recup_T": fT.mean(), "T_lo": np.percentile(fT, 2.5), "T_hi": np.percentile(fT, 97.5),
                           "verdade_S": kS ** 2, "recup_S": fS.mean(), "S_lo": np.percentile(fS, 2.5), "S_hi": np.percentile(fS, 97.5)})
            L = linhas[-1]
            print(f"{kT:5.2f} {kS:5.2f} | {L['verdade_T']:10.2f} {L['recup_T']:9.2f} {L['T_lo']:6.2f}–{L['T_hi']:<6.2f} | {L['verdade_S']:10.2f} {L['recup_S']:9.2f} {L['S_lo']:6.2f}–{L['S_hi']:<6.2f}", flush=True)
    # e se o excesso do TESS nao for barra, mas ESTRUTURA fora do espaco do
    # modelo? Senoide de amplitude A e periodo 2 anos so nas epocas TESS
    # (barras corretas nos dois conjuntos): o que a parabola nao absorve,
    # para onde vai?
    def estimador_estrutura(des, A_min, rng):
        zT, hT, zS, hS = [], [], [], []
        for E, sig, tess in des:
            r = rng.normal(0, sig)
            t_rel = (E - E[tess].mean()) * 0.5 / 365.25     # P ~ 0,5 d: E para anos (aproximacao suficiente)
            r = r + np.where(tess, (A_min / 1440.0) * np.sin(2 * np.pi * t_rel / 2.0 + rng.uniform(0, 2 * np.pi)), 0.0)
            w = 1 / sig
            Amat = np.vander(E, 3) * w[:, None]
            coef, *_ = np.linalg.lstsq(Amat, r * w, rcond=None)
            res = r - np.polyval(coef, E)
            z = res / sig
            h = alavancas(E, sig)
            zT += list(z[tess] ** 2); hT += list(1 - h[tess]); zS += list(z[~tess] ** 2); hS += list(1 - h[~tess])
        return sum(zT) / sum(hT), sum(zS) / sum(hS)
    print("\n  estrutura coerente so no TESS (senoide, P = 2 a), barras corretas nos dois conjuntos:")
    for A in (1.0, 3.0, 6.0):
        fT, fS = np.array([estimador_estrutura(des, A, rng) for _ in range(N_SIM)]).T
        linhas.append({"kT": np.nan, "kS": np.nan, "estrutura_T_min": A, "recup_T": fT.mean(), "T_lo": np.percentile(fT, 2.5), "T_hi": np.percentile(fT, 97.5),
                       "recup_S": fS.mean(), "S_lo": np.percentile(fS, 2.5), "S_hi": np.percentile(fS, 97.5)})
        print(f"    A = {A:.0f} min: recup. T {fT.mean():.2f} ({np.percentile(fT, 2.5):.2f}–{np.percentile(fT, 97.5):.2f}) | recup. S {fS.mean():.2f} ({np.percentile(fS, 2.5):.2f}–{np.percentile(fS, 97.5):.2f})", flush=True)
    pd.DataFrame(linhas).to_parquet(BASE / "vies_alavanca_26.parquet", index=False)
    print("\n  leitura: recup. > verdade = o estimador SUPERestima aquele conjunto (o valor medido e teto); recup. < verdade = SUBestima (piso).")
    print(f"  -> {BASE / 'vies_alavanca_26.parquet'}")
