# -*- coding: utf-8 -*-
"""Passo 3 da rodada 'quatro medidas': epsilon por alvo para orbitas de terceiro
corpo LONGAS e EXCENTRICAS, em tres faixas de P3, com tres eficiencias.

Expectativa em `expectativas/epsilon_longo_26.md`, commit c6ac0d5, escrita
ANTES de rodar. Estende `epsilon_26.py` (mesma cadeia: `cadeia_oc.ajustar`,
p_curv por Delta chi2, `cadeia_oc.adversarial`, p_gof com N - 3 gl; efemeride
exatamente linear nos instantes reais com as barras reais).

Populacao por faixa (0,1-20, 20-50, 50-100 a; P3 log-uniforme dentro da
faixa), e uniforme 0-0,8, omega uniforme, anomalia media em t_ref uniforme,
M_bin U(1,5, 2,5), M3 U(0,1, 1,0), cos i uniforme (Apendice C). LTTE
excentrico (Irwin 1952): Delta t = A [ (1 - e^2)/(1 + e cos nu) sin(nu + omega)
+ e sin omega ], A = a12 sin i / c = `epsilon_26.amplitude_min`; Kepler por
Newton. Tres eficiencias: eps (barras fixas, sem ruido), eps_re (sinal +
ruido, barras re-estimadas pela regra do pipeline -
`nulo_barras_reestimadas_26.barras_reestimadas`, como completude_secular_26),
eps3 (fixas, deteccao e p_gof >= 0,05). Controles: (A) epsilon_26 reproduz
6,37 / 4,25; (B) gerador novo com e = 0 em 0,1-20 a a < 0,15 de 6,37; (C) V527
Dra em 360 fases reproduz v527_ltte (93% / 84%). Fracao por faixa contada na
Tabela 5 do Tokovinin 2006 (cache data/catalogs/tokovinin2006_table5.parquet),
definicao estrita (terciario mais proximo = par cujo componente e o pai do
par interno) e frouxa (qualquer par do sistema; a da nota).

RODADA CINCO, Passo C (`--gp`): as mesmas orbitas (tres faixas, e <= 0,8,
mesma populacao e semente) injetadas COM o ruido GP do Passo B - g ~ N(0, K)
nas epocas TESS do desenho com os hiperparametros agrupados de
deriva_gp_26.json (Matern-3/2, o preferido, e exponencial) mais ruido branco
N(0, sigma_i) em todas as epocas - e julgadas pelo MESMO teste do Passo B:
C = diag + K(TESS), GLS linear/quadratica, p_curv por Delta chi2, adversarial
marginal na metrica de C (`gls_vies_26.adversarial_gls`), p_gof. eps_gp por
alvo e faixa; eps3_gp (p_gof >= 0,05); e a taxa de falso alarme do teste com
o mesmo ruido e sem orbita (controle, 500 sorteios). f por faixa = contagem
ESTRITA (8/0/3 de 42); intervalo de N = sum_b f_b Sigma eps_b por Monte Carlo
com f_b ~ Beta de Jeffreys (k + 1/2, n - k + 1/2), 20 000 sorteios. N esperado
contra os 7 sobreviventes de B. Saidas: epsilon_gp_26.parquet, epsilon_gp_26.json.

EXPECTATIVA (Passo C, escrita e commitada ANTES de rodar):
  - Falso alarme (sem orbita, ruido GP + branco, teste GP): 2-6% por alvo,
    mediana ~4%; soma nos 26: 0,5-1,5.
  - Sigma eps_gp (Matern): 0,1-20 a 4,0-5,5 (era 5,79 com barras fixas: o GP
    tira peso da curvatura interna ao TESS); 20-50 a 13-17 (era 16,5: a
    alavanca SuperWASP carrega); 50-100 a 8-11 (era 10,8). Exponencial
    dentro de +-10% do Matern.
  - N central (f estrito): 1,5-2,0 no total (0,1-20 a 0,8-1,1; 20-50 a 0;
    50-100 a 0,6-0,8); intervalo Jeffreys 95%: 0,7-4,5, o limite superior
    dominado pelo 0/42 de 20-50 a (Beta(0,5; 42,5): p97,5 ~ 0,06 -> x16 ~ 1).
  - Entre os 7 sobreviventes de B: N ~ 1,5-2,0 esperados como orbitas
    subamostradas (V527 Dra e um, conhecido) - "at most one or two" segue
    valendo no central e nao no superior do intervalo.
Desvio em qualquer direcao vai ao relatorio.

RESULTADO (Passo C, 2026-09-17, contra a expectativa de 2ba3628): falso alarme
mediana 0,000, soma 0,00, max 0,001 (esperado 2-6%: DESVIO da expectativa,
nao do teste - o adversarial de 1 sigma na metrica de C zera o falso
alarme; completude_secular_26 ja dava 0,2% na cadeia diagonal). Matern:
Sigma eps_gp 3,88 / 13,66 / 7,18 (razao contra as fixas 0,67 / 0,83 /
0,66; esperado 4-5,5 / 13-17 / 8-11: a curta e a longa abaixo); exponencial
a < 1%. N central 1,25 (0,74 + 0 + 0,51; esperado 1,5-2,0: abaixo),
Jeffreys 95% [0,80, 2,46], mediana 1,46; N eps3 0,91. Entre os 7
sobreviventes de B, 1,25 esperados como orbitas (18%). "At most one or
two" vale no central e falha por 0,46 no superior de 97,5%.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cadeia_oc as co  # noqa: E402
import config  # noqa: E402
import epsilon_26 as EPS  # noqa: E402
import nulo_barras_reestimadas_26 as NULO  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
BANDAS = (("0,1-20 a", 0.1, 20.0), ("20-50 a", 20.0, 50.0), ("50-100 a", 50.0, 100.0))
N_DRAW = 2000
SEMENTE = 20260916
E_MAX = 0.8
P_DET = 0.05
EXPECTATIVA = "c6ac0d5"
S_POR_ANO = 365.25 * 86400.0


def kepler(M, e, n_iter=30):
    Ecc = M + e * np.sin(M)
    for _ in range(n_iter):
        Ecc = Ecc - (Ecc - e * np.sin(Ecc) - M) / (1 - e * np.cos(Ecc))
    return Ecc


def ltte_min(t_d, t_ref_d, P3_anos, A_min, e, omega, M0):
    """Irwin 1952. t em dias; devolve Delta t em minutos."""
    M = M0 + 2 * np.pi * (t_d - t_ref_d) / (P3_anos * 365.25)
    Ecc = kepler(M, e)
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(Ecc / 2), np.sqrt(1 - e) * np.cos(Ecc / 2))
    return A_min * ((1 - e ** 2) / (1 + e * np.cos(nu)) * np.sin(nu + omega) + e * np.sin(omega))


def populacao(rng, lo, hi, n, e_max):
    P3 = 10 ** rng.uniform(np.log10(lo), np.log10(hi), n)
    mb = rng.uniform(*EPS.MBIN, n)
    m3 = rng.uniform(*EPS.M3, n)
    ci = rng.uniform(0, 1, n)
    M0 = rng.uniform(0, 2 * np.pi, n)
    omega = rng.uniform(0, 2 * np.pi, n)
    e = rng.uniform(0, e_max, n) if e_max > 0 else np.zeros(n)
    return P3, EPS.amplitude_min(P3, mb, m3, ci), e, omega, M0


def testes(E, t, sig_d, P):
    """p_curv, p_adv (so onde a curvatura passou, como epsilon_26), p_gof, dP/dt."""
    _, _, chi2l, _ = co.ajustar(E, t, sig_d, 1)
    coef, _, chi2p, _ = co.ajustar(E, t, sig_d, 2)
    p_curv = float(stats.chi2.sf(chi2l - chi2p, 1))
    p_gof = float(stats.chi2.sf(chi2p, len(E) - 3))
    if p_curv < P_DET:
        ag = pd.DataFrame({"E": E, "t0": t, "sig_min": sig_d * 1440.0})
        p_adv = float(co.adversarial(ag, P, float(t[0]))["p_adversarial"])
    else:
        p_adv = 1.0
    return p_curv, p_adv, p_gof, 2 * coef[0] / P * S_POR_ANO


def desenho(tic):
    j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
    pts = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
    return {"tic": int(tic), "P": float(j["P_escada_d"]), "E": pts.E.values.astype(float), "t0": pts.t0.values.astype(float),
            "sig": pts.sig_min.values.astype(float), "formal": pts.sig_formal_min.values.astype(float),
            "tess": pts.fonte.str.startswith("TESS").values, "pontos": pts}


def epsilon_alvo_banda(d, rng, lo, hi, e_max=E_MAX, n_draw=N_DRAW):
    E, P, sig = d["E"], d["P"], d["sig"]
    t_lin = d["t0"][0] + P * (E - E[0])
    sig_d = sig / 1440.0
    P3, A, e, omega, M0 = populacao(rng, lo, hi, n_draw, e_max)
    det = np.zeros(n_draw, bool); det3 = np.zeros(n_draw, bool); det_re = np.zeros(n_draw, bool)
    for k in range(n_draw):
        tau = ltte_min(t_lin, t_lin[0], P3[k], A[k], e[k], omega[k], M0[k])
        # barras fixas, sem ruido (epsilon_26)
        pc, pa, pg, _ = testes(E, t_lin + tau / 1440.0, sig_d, P)
        det[k] = pc < P_DET and pa < P_DET
        det3[k] = det[k] and pg >= P_DET
        # barras re-estimadas: sinal + ruido, regra do pipeline (completude_secular_26)
        r_min = tau + rng.normal(0, sig)
        sig_re = NULO.barras_reestimadas(d, r_min, rng)
        pc, pa, _, _ = testes(E, r_min / 1440.0, sig_re / 1440.0, P)
        det_re[k] = pc < P_DET and pa < P_DET
    barra = float(np.median(sig))
    return {"epsilon": float(det.mean()), "epsilon3": float(det3.mean()), "epsilon_re": float(det_re.mean()),
            "A_mediana_min": float(np.median(A)), "frac_A_acima_barra": float((A > barra).mean()), "barra_mediana_min": barra,
            "e_media": float(e.mean()), "n_draw": n_draw}


def epsilon_alvo_banda_gp(d, rng, lo, hi, nome, A, tc, e_max=E_MAX, n_draw=N_DRAW):
    """Passo C: orbita + ruido GP (K nas TESS) + branco, julgada com a C do Passo B (GLS). Tambem o falso alarme."""
    import deriva_gp_26 as GP
    import gls_vies_26 as G
    E, P, sig = d["E"], d["P"], d["sig"]
    t_lin = d["t0"][0] + P * (E - E[0])
    C, it = GP.cov_desenho(d["pontos"], nome, A, tc)
    K_min = GP.kernel(nome, np.abs(np.subtract.outer(d["t0"][it], d["t0"][it])), A, tc)
    Lk = np.linalg.cholesky(K_min + 1e-9 * np.eye(len(it)))
    P3, Aamp, e, omega, M0 = populacao(rng, lo, hi, n_draw, e_max)
    det = np.zeros(n_draw, bool); det3 = np.zeros(n_draw, bool); fa = np.zeros(n_draw, bool)
    dof = len(E) - 3

    def julgar(t):
        _, _, chi2l, _ = G.ajustar_gls(E, t, C, 1)
        _, _, chi2p, _ = G.ajustar_gls(E, t, C, 2)
        pc = float(stats.chi2.sf(chi2l - chi2p, 1))
        pa = G.adversarial_gls(E, t, C) if pc < P_DET else 1.0
        return pc, pa, float(stats.chi2.sf(chi2p, dof))
    for k in range(n_draw):
        ruido = rng.normal(0, sig)
        ruido[it] += Lk @ rng.normal(0, 1, len(it))
        tau = ltte_min(t_lin, t_lin[0], P3[k], Aamp[k], e[k], omega[k], M0[k])
        pc, pa, pg = julgar(t_lin + (tau + ruido) / 1440.0)
        det[k] = pc < P_DET and pa < P_DET
        det3[k] = det[k] and pg >= P_DET
        if k < n_draw // 4:                                   # falso alarme: mesmo ruido, sem orbita (500 sorteios)
            pc0, pa0, _ = julgar(t_lin + ruido / 1440.0)
            fa[k] = pc0 < P_DET and pa0 < P_DET
    return {"epsilon_gp": float(det.mean()), "epsilon3_gp": float(det3.mean()), "falso_alarme_gp": float(fa[: n_draw // 4].mean()),
            "A_mediana_min": float(np.median(Aamp)), "n_draw": n_draw}


def passo_c(tics, des, tok):
    hip = json.loads((BASE / "deriva_gp_26.json").read_text(encoding="utf-8"))["hiper_agrupados"]
    linhas = []
    for nome in ("matern32", "exp"):
        A, tc = hip[nome]["A_min"], hip[nome]["tau_c_d"]
        for rot, lo, hi in BANDAS:
            rng = np.random.default_rng(SEMENTE + int(lo * 10))          # a MESMA semente do Passo 3 por faixa
            for t in tics:
                r = epsilon_alvo_banda_gp(des[t], rng, lo, hi, nome, A, tc)
                linhas.append({"kernel": nome, "faixa": rot, "tic": t, **r})
                print(f"{nome:8s} {rot:9s} TIC {t}: eps_gp {r['epsilon_gp']:.2f} eps3_gp {r['epsilon3_gp']:.2f} FA {r['falso_alarme_gp']:.3f}", flush=True)
    R = pd.DataFrame(linhas); R.to_parquet(BASE / "epsilon_gp_26.parquet", index=False)
    B = pd.read_parquet(BASE / "deriva_gp_26_desenho.parquet")
    Rl = pd.read_parquet(BASE / "epsilon_longo_26.parquet")
    rng = np.random.default_rng(7); NMC = 20000
    out = {"kernels": {}}
    for nome, g in R.groupby("kernel", sort=False):
        sob = B[(B.kernel == nome) & B.em_D11 & B.detecta]
        fa = g.groupby("tic").falso_alarme_gp.mean()
        res = {"falso_alarme_mediana": float(fa.median()), "falso_alarme_soma_26": float(fa.sum()), "falso_alarme_max": float(fa.max()), "faixas": {}}
        N_tot = np.zeros(NMC); N3_tot = np.zeros(NMC); Nc = 0.0; N3c = 0.0
        print(f"\n== {nome}: falso alarme mediana {fa.median():.3f}, soma nos 26 {fa.sum():.2f}, max {fa.max():.3f} ({int(fa.idxmax())})")
        for rot, lo, hi in BANDAS:
            gb = g[g.faixa == rot]; S = float(gb.epsilon_gp.sum()); S3 = float(gb.epsilon3_gp.sum())
            k, n = tok[rot]["k_estrito"], tok[rot]["n"]
            f_mc = rng.beta(k + 0.5, n - k + 0.5, NMC)
            N_tot += f_mc * S; N3_tot += f_mc * S3; Nc += k / n * S; N3c += k / n * S3
            S_fix = float(Rl[Rl.faixa == rot].epsilon.sum())
            res["faixas"][rot] = {"sigma_eps_gp": S, "sigma_eps3_gp": S3, "sigma_eps_fixas_passo3": S_fix, "razao_gp_fixas": S / S_fix, "k": k, "n": n,
                                  "N_central": k / n * S, "N_jeffreys_p2.5": float(np.percentile(f_mc * S, 2.5)), "N_jeffreys_p50": float(np.percentile(f_mc * S, 50)),
                                  "N_jeffreys_p97.5": float(np.percentile(f_mc * S, 97.5))}
            print(f"  {rot:9s}: Sigma eps_gp {S:.2f} (fixas Passo 3 {S_fix:.2f}, razao {S / S_fix:.2f}), eps3_gp {S3:.2f} | f {k}/{n} | N {k / n * S:.2f}, Jeffreys 95% "
                  f"[{np.percentile(f_mc * S, 2.5):.2f}, {np.percentile(f_mc * S, 97.5):.2f}]")
        res["total"] = {"N_central": Nc, "N_jeffreys_p2.5": float(np.percentile(N_tot, 2.5)), "N_jeffreys_p50": float(np.percentile(N_tot, 50)), "N_jeffreys_p97.5": float(np.percentile(N_tot, 97.5)),
                        "N3_central": N3c, "N3_jeffreys_p97.5": float(np.percentile(N3_tot, 97.5)), "sobreviventes_B": int(len(sob)), "quais_B": sorted(int(t) for t in sob.tic)}
        frase = "sobrevive" if Nc <= 2 else "NAO sobrevive"
        frase_sup = "sobrevive" if np.percentile(N_tot, 97.5) <= 2 else "NAO sobrevive"
        print(f"  TOTAL: N {Nc:.2f} (Jeffreys 95% [{np.percentile(N_tot, 2.5):.2f}, {np.percentile(N_tot, 97.5):.2f}], mediana {np.percentile(N_tot, 50):.2f}); N eps3 {N3c:.2f} (ate {np.percentile(N3_tot, 97.5):.2f}) "
              f"| sobreviventes de B: {len(sob)} de 11 -> {Nc:.2f} deles esperados como orbitas ({Nc / len(sob):.0%}); at most one or two: {frase} no central, {frase_sup} no superior")
        out["kernels"][nome] = res
    out["tokovinin_estrito"] = {rot: {"k": tok[rot]["k_estrito"], "n": tok[rot]["n"]} for rot, _, _ in BANDAS}
    (BASE / "epsilon_gp_26.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print("  -> epsilon_gp_26.parquet, epsilon_gp_26.json")


def wilson(k, n, z=1.959964):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    den = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return float(c - h), float(c + h)


def contar_tokovinin():
    t5 = pd.read_parquet(config.CATALOGS / "tokovinin2006_table5.parquet")
    t5["Per_yr"] = np.where(t5.x_Per == "d", t5.Per / 365.25, np.where(t5.x_Per == "y", t5.Per, t5.Per * 1000.0))
    t5[["c1", "c2", "pai"]] = t5.Comp.str.split(",", expand=True)
    inner = t5[(t5.x_Per == "d") & (t5.Per < 3.0)]
    assert len(inner) == 42 and inner.HIP.nunique() == 42, len(inner)
    linhas = []
    for r in inner.itertuples():
        mesmo = t5[t5.HIP == r.HIP]
        if r.pai == "*":
            p3 = np.nan
        else:
            outer = mesmo[(mesmo.c1 == r.pai) | (mesmo.c2 == r.pai)]
            p3 = float(outer.Per_yr.min()) if len(outer) else np.nan
        outros = mesmo[~((mesmo.x_Per == "d") & (mesmo.Per < 3.0))]
        linhas.append({"HIP": int(r.HIP), "Comp": r.Comp, "P_in_d": float(r.Per), "P3_estrito_a": p3,
                       "P3_frouxo_min_a": float(outros.Per_yr.min()) if len(outros) else np.nan,
                       "algum_par_0.1_20": bool(((outros.Per_yr >= 0.1) & (outros.Per_yr < 20)).any())})
    R = pd.DataFrame(linhas)
    out = {"n_sistemas": int(len(R)), "n_com_terciario": int(R.P3_estrito_a.notna().sum()), "n_sem_terciario": int(R.P3_estrito_a.isna().sum())}
    for rot, lo, hi in BANDAS:
        k = int(((R.P3_estrito_a >= lo) & (R.P3_estrito_a < hi)).sum())
        w = wilson(k, len(R))
        out[rot] = {"k_estrito": k, "n": int(len(R)), "f_estrito": k / len(R), "wilson_lo": w[0], "wilson_hi": w[1],
                    "HIP": sorted(int(h) for h in R.HIP[(R.P3_estrito_a >= lo) & (R.P3_estrito_a < hi)])}
    k_frouxo = int(R["algum_par_0.1_20"].sum())
    out["0,1-20 a"]["k_frouxo"] = k_frouxo
    out["0,1-20 a"]["f_frouxo"] = k_frouxo / len(R)
    out["0,1-20 a"]["HIP_so_no_frouxo"] = sorted(int(h) for h in R.HIP[R["algum_par_0.1_20"] & ~((R.P3_estrito_a >= 0.1) & (R.P3_estrito_a < 20))])
    return out, R


def controle_v527():
    """A orbita publicada do V527 Dra pelo gerador novo (e = 0), nas convencoes de v527_ltte.py."""
    j = json.loads((BASE / "oc__424461577.json").read_text(encoding="utf-8"))
    pontos = pd.DataFrame(j["pontos"]); P = float(j["P_escada_d"])
    E = pontos.E.values.astype(float); sig = pontos.sig_min.values / 1440.0
    t_lin = pontos.t0.values[pontos.fonte.values == "TESS s19"][0] + P * E
    linhas = []
    for fase in np.arange(0, 1, 1 / 360):
        tau = ltte_min(t_lin, t_lin[2], 2.733, 6.6, 0.0, 0.0, 2 * np.pi * fase)
        t = t_lin + tau / 1440.0
        _, _, chi2l, _ = co.ajustar(E, t, sig, 1)
        coef, _, chi2p, _ = co.ajustar(E, t, sig, 2)
        ag = pd.DataFrame({"E": E, "t0": t, "sig_min": pontos.sig_min.values})
        linhas.append({"fase": fase, "dPdt": 2 * coef[0] / P * S_POR_ANO, "p_curv": float(stats.chi2.sf(chi2l - chi2p, 1)),
                       "p_adv": float(co.adversarial(ag, P, float(t[2]))["p_adversarial"])})
    r = pd.DataFrame(linhas)
    ref = pd.read_parquet(BASE / "v527_ltte_injetado.parquet")
    return {"frac_p_curv": float((r.p_curv < 0.05).mean()), "frac_p_adv": float((r.p_adv < 0.05).mean()),
            "dPdt_min": float(r.dPdt.min()), "dPdt_max": float(r.dPdt.max()),
            "ref_frac_p_curv": float((ref.p_curv < 0.05).mean()), "ref_frac_p_adv": float((ref.p_adv < 0.05).mean()),
            "max_dif_dPdt_vs_ref": float(np.abs(r.dPdt.values - ref.dPdt.values).max())}


if __name__ == "__main__":
    t26 = pd.read_parquet(BASE / "tabela_26.parquet")
    tics = sorted(int(t) for t in t26.TIC)
    des = {t: desenho(t) for t in tics}
    if "--gp" in sys.argv:
        passo_c(tics, des, contar_tokovinin()[0])
        sys.exit()
    eps_ref = pd.read_parquet(BASE / "epsilon_26.parquet").set_index("tic")
    # controle A: epsilon_26 tal qual (mesma semente) reproduz o parquet gravado
    rng = np.random.default_rng(EPS.SEMENTE)
    ctrlA = {t: EPS.epsilon_alvo(des[t]["pontos"], des[t]["P"], rng) for t in tics}
    sA = sum(v["epsilon"] for v in ctrlA.values()); sA3 = sum(v["epsilon3"] for v in ctrlA.values())
    difA = max(abs(ctrlA[t]["epsilon"] - eps_ref.loc[t, "epsilon"]) for t in tics)
    print(f"controle A (epsilon_26 re-rodado): Sigma eps {sA:.3f} (gravado {eps_ref.epsilon.sum():.3f}), Sigma eps3 {sA3:.3f} (gravado {eps_ref.epsilon3.sum():.3f}); max |dif| por alvo {difA:.4f}")
    assert difA < 1e-9, "epsilon_26 nao reproduz o parquet gravado"
    # controle C: V527
    cC = controle_v527()
    print(f"controle C (V527 Dra, gerador novo e = 0): p_curv < 0,05 em {cC['frac_p_curv']:.0%} (v527_ltte {cC['ref_frac_p_curv']:.0%}), "
          f"p_adv < 0,05 em {cC['frac_p_adv']:.0%} ({cC['ref_frac_p_adv']:.0%}); dP/dt {cC['dPdt_min']:+.4f}..{cC['dPdt_max']:+.4f}; max |dif dP/dt| {cC['max_dif_dPdt_vs_ref']:.2e}")
    # controle B: gerador novo com e = 0 na faixa curta
    rng = np.random.default_rng(SEMENTE)
    ctrlB = {t: epsilon_alvo_banda(des[t], rng, 0.1, 20.0, e_max=0.0) for t in tics}
    sB = sum(v["epsilon"] for v in ctrlB.values()); sB3 = sum(v["epsilon3"] for v in ctrlB.values()); sBre = sum(v["epsilon_re"] for v in ctrlB.values())
    print(f"controle B (gerador novo, e = 0, 0,1-20 a): Sigma eps {sB:.3f} (esperado 6,37 +- 0,15), eps3 {sB3:.3f}, eps_re {sBre:.3f}", flush=True)
    # as tres faixas com excentricidade
    linhas = []
    for rot, lo, hi in BANDAS:
        rng = np.random.default_rng(SEMENTE + int(lo * 10))
        for t in tics:
            r = epsilon_alvo_banda(des[t], rng, lo, hi)
            linhas.append({"faixa": rot, "P3_lo": lo, "P3_hi": hi, "tic": t, "n_epocas": len(des[t]["E"]), **r,
                           "epsilon_curto_gravado": float(eps_ref.loc[t, "epsilon"])})
            print(f"{rot:9s} TIC {t}: A med {r['A_mediana_min']:6.1f} min ({r['frac_A_acima_barra']:.0%} > barra {r['barra_mediana_min']:.1f}) | "
                  f"eps {r['epsilon']:.2f} eps3 {r['epsilon3']:.2f} eps_re {r['epsilon_re']:.2f}", flush=True)
    R = pd.DataFrame(linhas)
    R.to_parquet(BASE / "epsilon_longo_26.parquet", index=False)
    tok, Rt = contar_tokovinin()
    Rt.to_parquet(BASE / "tokovinin_42.parquet", index=False)
    print(f"\n== Tokovinin 2006 Tabela 5: {tok['n_sistemas']} pares P < 3 d, {tok['n_com_terciario']} com terciario, {tok['n_sem_terciario']} sem")
    resumo = {"expectativa_commit": EXPECTATIVA, "controles": {"A_sigma_eps": sA, "A_sigma_eps3": sA3, "A_max_dif": difA, "B_sigma_eps_e0": sB, "B_sigma_eps3_e0": sB3, "B_sigma_eps_re_e0": sBre, "C_v527": cC},
              "tokovinin": tok, "faixas": {}}
    tot = {"fix": 0.0, "re": 0.0, "e3": 0.0, "fix_hi": 0.0, "re_hi": 0.0, "fix_lo": 0.0}
    for rot, lo, hi in BANDAS:
        g = R[R.faixa == rot]; k = tok[rot]
        S, S3, Sre = float(g.epsilon.sum()), float(g.epsilon3.sum()), float(g.epsilon_re.sum())
        f, wlo, whi = k["f_estrito"], k["wilson_lo"], k["wilson_hi"]
        resumo["faixas"][rot] = {"sigma_eps": S, "sigma_eps3": S3, "sigma_eps_re": Sre, "eps_mediana": float(g.epsilon.median()), "eps_re_mediana": float(g.epsilon_re.median()),
                                 "f_estrito": f, "k": k["k_estrito"], "wilson": [wlo, whi], "N_fixas": f * S, "N_fixas_wilson": [wlo * S, whi * S],
                                 "N_re": f * Sre, "N_re_wilson": [wlo * Sre, whi * Sre], "N_eps3": f * S3, "N_eps3_wilson": [wlo * S3, whi * S3]}
        if rot == "0,1-20 a":
            resumo["faixas"][rot].update({"f_frouxo": k["f_frouxo"], "N_fixas_frouxo": k["f_frouxo"] * S, "N_re_frouxo": k["f_frouxo"] * Sre, "N_eps3_frouxo": k["f_frouxo"] * S3})
        tot["fix"] += f * S; tot["re"] += f * Sre; tot["e3"] += f * S3; tot["fix_hi"] += whi * S; tot["re_hi"] += whi * Sre; tot["fix_lo"] += wlo * S
        print(f"  {rot:9s}: f estrito {k['k_estrito']}/42 = {f:.3f} (Wilson {wlo:.3f}-{whi:.3f}) | Sigma eps {S:.2f}, eps3 {S3:.2f}, eps_re {Sre:.2f} | "
              f"N fixas {f * S:.2f} ({wlo * S:.2f}-{whi * S:.2f}), N re {f * Sre:.2f} ({wlo * Sre:.2f}-{whi * Sre:.2f}), N eps3 {f * S3:.2f}"
              + (f" | frouxo {k['k_frouxo']}/42: N fixas {k['f_frouxo'] * S:.2f}, re {k['f_frouxo'] * Sre:.2f}" if rot == "0,1-20 a" else ""))
    resumo["total"] = {"N_fixas": tot["fix"], "N_fixas_wilson_lo": tot["fix_lo"], "N_fixas_wilson_hi": tot["fix_hi"], "N_re": tot["re"], "N_re_wilson_hi": tot["re_hi"], "N_eps3": tot["e3"]}
    print(f"  TOTAL (f estrito): N fixas {tot['fix']:.2f} (Wilson {tot['fix_lo']:.2f}-{tot['fix_hi']:.2f}), N re-estimadas {tot['re']:.2f} (ate {tot['re_hi']:.2f}), N eps3 {tot['e3']:.2f}")
    print(f"  'at most one or two of the eleven': {'sobrevive' if tot['fix'] <= 2.0 else 'NAO sobrevive'} com barras fixas ({tot['fix']:.2f}); "
          f"{'sobrevive' if tot['re'] <= 2.0 else 'NAO sobrevive'} com re-estimadas ({tot['re']:.2f}); Wilson superior fixas {tot['fix_hi']:.2f}")
    (BASE / "epsilon_longo_26.json").write_text(json.dumps(resumo, indent=1, default=float), encoding="utf-8")
    print(f"  -> {BASE / 'epsilon_longo_26.parquet'}, tokovinin_42.parquet, epsilon_longo_26.json")
