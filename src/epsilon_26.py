# -*- coding: utf-8 -*-
"""epsilon por alvo: a fracao de orbitas de terceiro corpo que, amostradas
nas epocas REAIS de cada um dos 26 com as barras reais, viram "curvatura
detectada" pela cadeia desta nota (p_curv < 0,05 E p_adv < 0,05), e a
fracao que alem disso passa na bondade de ajuste (p_gof >= 0,05).

Segunda rodada da revisao, item 1: o epsilon = 0,83 de `v527_ltte.py` foi
medido para UMA amplitude, UM P3 e UM padrao de seis epocas em tres grupos
e estava sendo emprestado a 26 alvos com N = 4-7 e padroes diferentes.
Aqui ele e medido em cada alvo sobre uma populacao de orbitas.

Populacao de orbitas (declarada, nao ajustada):
  P3        log-uniforme em 0,1-20 anos (entre um setor e a alavanca)
  M_bin     uniforme 1,5-2,5 Msol (binarias de P < 5 d, Tmag 9-12)
  M3        uniforme 0,1-1,0 Msol
  cos i     uniforme em [0, 1] (orientacao isotropica)
  fase      uniforme em [0, 1)
  e = 0     (circular; a excentricidade nao entra no que se testa)
  A = a3 * M3/(M_bin+M3) * sin i / c, com a3^3 = (M_bin+M3) P3^2 (UA, Msol, anos)
Conferencia da relacao: M_bin 2,5, M3 1,0, i 70 graus, P3 2,733 a -> A =
6,6 min, o valor publicado do V527 Dra.

Injecao: tau(t) = A sin(2 pi (t - t_ref)/P3 + 2 pi fase) sobre efemeride
exatamente linear nos instantes reais do alvo; ajuste linear e quadratico
com `cadeia_oc.ajustar` e as barras reais; p_curv por delta-chi2;
adversarial de `cadeia_oc.adversarial`; p_gof do chi2 da parabola com N-3
gl. Sem ruido alem do LTTE. N_DRAW orbitas por alvo, semente fixa.

Saida por alvo: N epocas, numero de grupos (epocas separadas por > 0,5 a),
barra mediana, epsilon (deteccao), epsilon3 (deteccao E p_gof >= 0,05),
fracao das orbitas com A > barra mediana, e epsilon condicionado a A >
barra. Resumo: media e soma de epsilon nos 26, e a correlacao (Spearman)
de epsilon com N e com o numero de grupos.

EXPECTATIVA (escrita antes de rodar, commitada antes de rodar):
  (a) epsilon varia bastante entre alvos - espero de < 0,3 a > 0,7;
  (b) epsilon e maior em alvos com epocas concentradas em poucos grupos
      (e o mecanismo do V527 Dra: tres parametros passam por tres grupos)
      e menor em alvos com mais epocas espalhadas - correlacao negativa
      com o numero de grupos, e negativa ou nula com N;
  (c) a media de epsilon nos 26 fica ABAIXO de 0,83, porque a populacao
      inclui orbitas de amplitude menor que a barra (M3 pequeno, i baixo,
      P3 curto), que nao produzem deteccao;
  (d) epsilon3 < epsilon em todos, e a diferenca e maior onde ha mais gl.
Desvio em qualquer direcao vai para a nota como esta.

RODADA NO ESTADO B (2026-09-13): as barras SuperWASP dos 26 cresceram por
1,0-1,73 (mediana 1,40) com o conserto do sqrt n. EXPECTATIVA (antes de
rodar): Sigma epsilon cai de 6,7 para ~6,0 (a rodada x1,65 deu 5,7 e a
x1,34 deu ~6,0); a queda concentra-se nos alvos cuja barra e dominada
pela dispersao entre temporadas; a ordem dos alvos por epsilon e a
correlacao com a barra (rho -0,64) mudam pouco.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cadeia_oc as C  # noqa: E402
import config  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
N_DRAW = 2000
SEMENTE = 20260913
P3_MIN, P3_MAX = 0.1, 20.0
MBIN = (1.5, 2.5)
M3 = (0.1, 1.0)
S_POR_ANO = 365.25 * 86400.0
UA_MIN = 499.004784 / 60.0     # tempo-luz de 1 UA, em minutos
P_DET = 0.05


def amplitude_min(P3_anos, m_bin, m3, cos_i):
    m_tot = m_bin + m3
    a3 = (m_tot * P3_anos ** 2) ** (1.0 / 3.0)        # UA
    return a3 * m3 / m_tot * np.sqrt(1 - cos_i ** 2) * UA_MIN


def grupos(t_dias, gap_anos=0.5):
    t = np.sort(np.asarray(t_dias))
    return int(1 + np.sum(np.diff(t) > gap_anos * 365.25))


def epsilon_alvo(pontos, P, rng):
    E = pontos.E.values.astype(float)
    t_lin = pontos.t0.values[0] + P * (E - E[0])          # efemeride linear exata
    sig_d = pontos.sig_min.values / 1440.0
    n = len(E)
    P3 = 10 ** rng.uniform(np.log10(P3_MIN), np.log10(P3_MAX), N_DRAW)
    mb = rng.uniform(*MBIN, N_DRAW)
    m3 = rng.uniform(*M3, N_DRAW)
    ci = rng.uniform(0, 1, N_DRAW)
    fase = rng.uniform(0, 1, N_DRAW)
    A = amplitude_min(P3, mb, m3, ci)
    det = np.zeros(N_DRAW, bool)
    det3 = np.zeros(N_DRAW, bool)
    for k in range(N_DRAW):
        tau = (A[k] / 1440.0) * np.sin(2 * np.pi * (t_lin - t_lin[0]) / (P3[k] * 365.25) + 2 * np.pi * fase[k])
        t = t_lin + tau
        _, _, chi2l, _ = C.ajustar(E, t, sig_d, 1)
        _, _, chi2p, _ = C.ajustar(E, t, sig_d, 2)
        p_curv = stats.chi2.sf(chi2l - chi2p, 1)
        if p_curv >= P_DET:
            continue
        ag = pd.DataFrame({"E": E, "t0": t, "sig_min": pontos.sig_min.values})
        p_adv = C.adversarial(ag, P, float(t[0]))["p_adversarial"]
        if p_adv >= P_DET:
            continue
        det[k] = True
        det3[k] = stats.chi2.sf(chi2p, n - 3) >= P_DET
    barra = float(np.median(pontos.sig_min.values))
    acima = A > barra
    return {"n_epocas": n, "n_grupos": grupos(pontos.t0.values), "barra_mediana_min": barra,
            "epsilon": float(det.mean()), "epsilon3": float(det3.mean()),
            "frac_A_acima_barra": float(acima.mean()),
            "epsilon_dado_A_acima": float(det[acima].mean()) if acima.any() else np.nan,
            "A_mediana_min": float(np.median(A)), "n_draw": N_DRAW}


if __name__ == "__main__":
    # conferencia da relacao de amplitude no caso publicado
    a_v527 = amplitude_min(2.733, 2.5, 1.0, np.cos(np.radians(70)))
    print(f"conferencia: M_bin 2,5 + M3 1,0, i 70 graus, P3 2,733 a -> A = {a_v527:.2f} min (publicado 6,6)")
    rng = np.random.default_rng(SEMENTE)
    t26 = pd.read_parquet(BASE / "tabela_26.parquet")
    linhas = []
    for tic in sorted(t26.TIC):
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        pontos = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
        r = epsilon_alvo(pontos, float(j["P_escada_d"]), rng)
        r["tic"] = int(tic)
        linhas.append(r)
        print(f"TIC {tic}: N={r['n_epocas']} grupos={r['n_grupos']} barra {r['barra_mediana_min']:.2f} min | "
              f"A mediana {r['A_mediana_min']:.1f} min, {r['frac_A_acima_barra']:.0%} acima da barra | "
              f"epsilon {r['epsilon']:.2f}  epsilon3 {r['epsilon3']:.2f}  epsilon|A>barra {r['epsilon_dado_A_acima']:.2f}", flush=True)
    d = pd.DataFrame(linhas)
    d.to_parquet(BASE / "epsilon_26.parquet", index=False)
    rho_g = stats.spearmanr(d.epsilon, d.n_grupos)
    rho_n = stats.spearmanr(d.epsilon, d.n_epocas)
    print("\n== resumo (expectativa: epsilon de < 0,3 a > 0,7; maior com menos grupos; media < 0,83; epsilon3 < epsilon)")
    print(f"  epsilon: min {d.epsilon.min():.2f}, mediana {d.epsilon.median():.2f}, media {d.epsilon.mean():.2f}, max {d.epsilon.max():.2f}; soma {d.epsilon.sum():.2f}")
    print(f"  epsilon3: mediana {d.epsilon3.median():.2f}, media {d.epsilon3.mean():.2f}; soma {d.epsilon3.sum():.2f}; epsilon3 < epsilon em {(d.epsilon3 < d.epsilon).sum()} de 26")
    print(f"  Spearman epsilon x grupos: rho {rho_g.statistic:+.2f} (p {rho_g.pvalue:.3f}); epsilon x N: rho {rho_n.statistic:+.2f} (p {rho_n.pvalue:.3f})")
    for g, s in d.groupby("n_grupos"):
        print(f"  grupos={g}: n={len(s)} epsilon mediana {s.epsilon.median():.2f}")
    print(f"  fracao das orbitas com A > barra: mediana {d.frac_A_acima_barra.median():.0%}; epsilon dado A > barra: mediana {d.epsilon_dado_A_acima.median():.2f}")
    print(f"  -> {BASE / 'epsilon_26.parquet'}")
