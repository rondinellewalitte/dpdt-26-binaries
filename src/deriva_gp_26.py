# -*- coding: utf-8 -*-
"""Rodada cinco, Passo B: ruido ENTRE setores como processo estacionario (GP),
com hiperparametros medidos fora do D11 e aplicados ao desenho.

O Passo 1 da rodada anterior mostrou o variograma dos residuos primarios no
bloco TESS estendido subindo de gama_norm 1,3 (< 30 d) a ~18 (1-4 a) e
saturando - nao e passeio aleatorio (o gate q > 0 a 2 sigma nao fechou). Aqui
o modelo e estacionario: t_i = a + b E_i + c E_i^2 (media quadratica por alvo)
+ g(t_i) + ruido branco sigma_i, com g ~ GP de kernel exponencial
k = A^2 exp(-tau/tau_c) ou Matern-3/2 k = A^2 (1 + sqrt3 tau/tau_c)
exp(-sqrt3 tau/tau_c). Dados: as epocas primarias do bloco TESS estendido de
`jitter_primarios_26.parquet` (799 setores, estado D), com E, fase e
criterios de `jitter_primarios_26.classificar_setores`.

1. Hiperparametros (A, tau_c) por MV AGRUPADA sobre os 15 alvos fora do D11
   (10 do norte informativos + 5 do sul com 3-6 epocas), com a media
   quadratica de cada alvo perfilada (GLS a cada C): -2 ln L = sum_alvos
   [ r^T C^-1 r + ln|C| ]. Grade A 0-6 min x tau_c 5-3000 d (log), refino
   Nelder-Mead; IC 68% por perfil 1D (Delta = 1). Os dois kernels.
2. Checagem: (A, tau_c) por alvo, cada um sozinho com a propria quadratica,
   nos 26 - os do D11 contra os agrupados.
3. Desenho (os `pontos` de cada um dos 26, nada re-extraido): C = diag(sigma^2)
   + K nas epocas TESS (K so entre TESS; SuperWASP diagonal, sem o bloco 0,78
   - a comparacao e com a cadeia diagonal e com os cenarios A/B/C do Passo 1).
   GLS linear e quadratica (`gls_vies_26.ajustar_gls`), p_curv por Delta chi2,
   adversarial de 1 sigma MARGINAL na metrica de C (`gls_vies_26.adversarial_gls`),
   p_gof; vao sP x N com sP da GLS linear so-TESS. Sobrevivencia de D11 e D9.
4. TIC 232634196: refit com 18 epocas (7 + 11 dos setores 76-86) e C com o GP
   nas 15 TESS; predicao para 2027-07-02 (E = 2206, o mesmo de
   previsao_232634196): media quadratica + posterior do GP no instante
   (kernel condicionado aos residuos TESS), sigma = sqrt(var modelo + var
   posterior GP); contra a separacao linear - quadratica.
Saidas: deriva_gp_26.json, deriva_gp_26_hiper.parquet (por alvo e kernel),
deriva_gp_26_desenho.parquet (os 26 x kernel).

EXPECTATIVA (escrita e commitada ANTES de rodar):
  - Agrupados (15 fora do D11), exponencial: A = 0,8-2,0 min, tau_c = 0,2-1,5 a
    (o variograma sem os 4 alvos de sinal real satura em gama_norm ~6 a
    partir de ~1 a; 229476285, sigma_j 2,5, e 359552377, 4,5, puxam A para
    cima). Matern-3/2: A dentro de +-30% e tau_c dentro de x2 do exponencial;
    |Delta(-2 ln L)| entre kernels < 6.
  - Por alvo no D11: A tipico 0,3-1,5 min; 232634196, 424461577 e 198388252
    acima de 2 min (curvatura real/LTTE dentro do bloco - nao e ruido);
    tau_c mal restringido (IC atravessa a grade) na maioria.
  - Desenho, kernel exponencial: D11 sobrevivem 3-6 de 11 (<= o cenario A,
    5: setores adjacentes de 2019-20 ficam correlacionados e valem menos que
    epocas brancas); 232634196 sobrevive; D9 adequados 3-5; nenhuma deteccao
    nova; vao nunca ambiguo (sP x N < 0,1). Matern: mesma contagem +-1.
  - 232634196, 2027-07-02: sigma da janela 2,5-4,0 min (era 2,2 na
    quadratica das 7 e 1,4 na linear); separacao linear - quadratica 16-18
    min, >= 4 sigma. A posterior do GP em 2027,5 (2,5 a apos s86) volta a
    media: |media posterior| < 0,5 min.
Desvio em qualquer direcao vai ao relatorio.

RESULTADO (2026-09-17, contra a expectativa de 3cce8cd): agrupado fora do D11
(392 epocas): exp A 0,82 min [0,8, 0,8], tau_c 91 d = 0,25 a [81, 96];
Matern A 0,84, tau_c 68 d = 0,19 a; -2lnL 414 vs 396 - Matern preferido por
18 (esperado < 6: DESVIO; o exponencial tem cuspide em tau = 0 e o residuo
entre setores adjacentes nao a tem). Por alvo no D11: 232634196 2,4,
424461577 3,8 (5,5 Matern), 198388252 1,8; os demais 0-1,0; tau_c na borda
em 14 de 26. Desenho: D11 sobrevivem 7 de 11 nos dois kernels (144194304,
198388252, 198408416, 232634196, 329246824, 392536812, 424461577), esperado
3-6: DESVIO - o A agrupado (0,82) e MENOR que o sigma_j proprio dos alvos
que caiam na inflacao branca (144194304 2,8, 424461577 3,6) e MAIOR que o
dos de barra minuscula (230386284: A proprio 0, barras 0,1 min), que
passam a cair (p_adv 0,06); caem 229914020, 230386284, 377253090,
390021728; D9 adequados 6 de 9; nenhuma nova; vao < 0,006. 232634196:
7 epocas quad +- 3,06 min (modelo 2,95, GP 0,82), separacao 21,0 min =
5,6 sigma; 18 epocas quad +- 1,94 (esperado 2,5-4,0: abaixo), p_gof 0,025
(o GP de 0,25 a nao absorve a tendencia de 10 meses), separacao 14,3 min
= 6,0 sigma; media posterior do GP em 2027,5 = 0,00 min.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
import gls_vies_26 as G  # noqa: E402
import jitter_primarios_26 as J  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
S_POR_ANO = 365.25 * 86400.0
P_DET = 0.05
GRADE_A = np.arange(0.0, 6.0 + 1e-9, 0.1)                          # min
GRADE_TC = np.exp(np.linspace(np.log(5.0), np.log(3000.0), 40))    # dias
E_PRED_232634196 = 2206
EXPECTATIVA = "3cce8cd"


def kernel(nome, tau_d, A, tc):
    if nome == "exp":
        return A ** 2 * np.exp(-tau_d / tc)
    if nome == "matern32":
        x = np.sqrt(3.0) * tau_d / tc
        return A ** 2 * (1 + x) * np.exp(-x)
    raise ValueError(nome)


def bloco_alvo(tic, S):
    """Epocas primarias ok do bloco TESS estendido: E, t (d), sigma (min)."""
    j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
    P, sP, T0, E_ult = J.escada_do_registro(j)
    d = J.classificar_setores(S[(S.tic == tic) & (S.minimo == "primario")], P, sP, T0, E_ult, 0.0)
    ok = d[d.ok].sort_values("t0_btjd")
    return {"tic": tic, "P": P, "E": ok.E.values.astype(float), "t_d": ok.t0_btjd.values.astype(float), "sig_min": ok.sig_min.values.astype(float), "n": int(len(ok))}


def menos2lnL_bloco(b, nome, A, tc, grau=2):
    """Perfil na media polinomial (GLS) para (A, tc) fixos; tudo em minutos."""
    n = b["n"]
    if n <= grau + 1:
        return 0.0
    t_min = (b["t_d"] - b["t_d"][0]) * 1440.0
    tau = np.abs(np.subtract.outer(b["t_d"], b["t_d"]))
    C = np.diag(b["sig_min"] ** 2) + kernel(nome, tau, A, tc)
    X = np.vander(b["E"], grau + 1)
    L = np.linalg.cholesky(C)
    Ci_X = np.linalg.solve(C, X); Ci_t = np.linalg.solve(C, t_min)
    beta = np.linalg.solve(X.T @ Ci_X, X.T @ Ci_t)
    r = t_min - X @ beta
    return float(r @ np.linalg.solve(C, r) + 2 * np.sum(np.log(np.diag(L))))


def ajustar_hiper(blocos, nome):
    """Grade + refino; IC68 por perfil 1D na grade."""
    grade = np.array([[sum(menos2lnL_bloco(b, nome, A, tc) for b in blocos) for tc in GRADE_TC] for A in GRADE_A])
    iA, iT = np.unravel_index(np.argmin(grade), grade.shape)
    f = lambda p: sum(menos2lnL_bloco(b, nome, abs(p[0]), np.exp(p[1])) for b in blocos)  # noqa: E731
    r = optimize.minimize(f, [GRADE_A[iA] if GRADE_A[iA] > 0 else 0.05, np.log(GRADE_TC[iT])], method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-4, "maxiter": 400})
    A, tc, L0 = abs(r.x[0]), float(np.exp(r.x[1])), float(r.fun)
    pA = grade.min(axis=1) - grade.min(); pT = grade.min(axis=0) - grade.min()
    okA = np.where(pA <= 1.0)[0]; okT = np.where(pT <= 1.0)[0]
    return {"A_min": A, "tau_c_d": tc, "tau_c_anos": tc / 365.25, "menos2lnL": L0, "menos2lnL_A0": float(grade[0].min()),
            "A_ic68": [float(GRADE_A[okA[0]]), float(GRADE_A[okA[-1]])], "tau_c_ic68_d": [float(GRADE_TC[okT[0]]), float(GRADE_TC[okT[-1]])],
            "tau_c_na_borda": bool(okT[0] == 0 or okT[-1] == len(GRADE_TC) - 1), "n_alvos": len(blocos), "n_epocas": int(sum(b["n"] for b in blocos))}


def cov_desenho(pontos, nome, A, tc):
    """C em dias^2: diagonal com as barras do desenho + K do GP so entre as epocas TESS."""
    sig = pontos.sig_min.values / 1440.0
    C = np.diag(sig ** 2)
    it = np.where(pontos.fonte.str.startswith("TESS").values)[0]
    tt = pontos.t0.values[it]
    C[np.ix_(it, it)] += kernel(nome, np.abs(np.subtract.outer(tt, tt)), A, tc) / 1440.0 ** 2
    return C, it


def refazer_gp(pontos, P, nome, A, tc):
    E = pontos.E.values.astype(float); t = pontos.t0.values.astype(float)
    C, it = cov_desenho(pontos, nome, A, tc)
    _, _, chi2l, _ = G.ajustar_gls(E, t, C, 1)
    coef, _, chi2p, cov = G.ajustar_gls(E, t, C, 2)
    dof = len(E) - 3
    # vao: sP da GLS linear so-TESS com a C do bloco
    _, _, _, cov_t = G.ajustar_gls(E[it], t[it], C[np.ix_(it, it)], 1)
    sP = float(np.sqrt(cov_t[0, 0])); N = float(np.max(np.abs(np.delete(E, it) - E[it].mean())))
    return {"dPdt": 2 * coef[0] / P * S_POR_ANO, "sdPdt": 2 * np.sqrt(cov[0, 0]) / P * S_POR_ANO, "chi2": float(chi2p), "dof": dof,
            "chi2r": chi2p / dof, "p_gof": float(stats.chi2.sf(chi2p, dof)), "p_curv": float(stats.chi2.sf(chi2l - chi2p, 1)),
            "p_adv": G.adversarial_gls(E, t, C), "sP_vezes_N": sP * N, "vao_ok": bool(sP * N < 0.25), "delta_chi2": float(chi2l - chi2p)}


def previsao_232634196(hip, nome):
    j = json.loads((BASE / "oc__232634196.json").read_text(encoding="utf-8"))
    pts = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
    nov = pd.read_parquet(BASE / "setores_76_86_232634196.parquet")
    P = float(j["P_escada_d"])
    p18 = pd.concat([pts[["fonte", "t0", "sig_min", "E"]], pd.DataFrame({"fonte": ["TESS s%d" % s for s in nov.setor], "t0": nov.t0_btjd, "sig_min": nov.sig_min, "E": nov.E})], ignore_index=True).sort_values("t0").reset_index(drop=True)
    out = {}
    for rot, pp in (("7", pts), ("18", p18)):
        E = pp.E.values.astype(float); t = pp.t0.values.astype(float)
        C, it = cov_desenho(pp, nome, hip["A_min"], hip["tau_c_d"])
        res = {}
        for grau, nome_g in ((1, "lin"), (2, "quad")):
            coef, r, chi2, cov = G.ajustar_gls(E, t, C, grau)
            x = np.array([E_PRED_232634196 ** k for k in range(grau, -1, -1)], float)
            m = float(x @ coef); var_m = float(x @ cov @ x)
            # posterior do GP no instante da predicao, condicionado aos residuos TESS
            tt = t[it]; Ct = C[np.ix_(it, it)]
            kst = kernel(nome, np.abs(m - tt), hip["A_min"], hip["tau_c_d"]) / 1440.0 ** 2
            mu_gp = float(kst @ np.linalg.solve(Ct, r[it])); var_gp = float(hip["A_min"] ** 2 / 1440.0 ** 2 - kst @ np.linalg.solve(Ct, kst))
            res[nome_g] = {"t_pred_bjd": 2457000.0 + m + mu_gp, "t_modelo_bjd": 2457000.0 + m, "mu_gp_min": mu_gp * 1440.0, "sig_modelo_min": np.sqrt(var_m) * 1440.0,
                           "sig_gp_min": np.sqrt(max(var_gp, 0.0)) * 1440.0, "sig_total_min": np.sqrt(var_m + max(var_gp, 0.0)) * 1440.0,
                           "chi2": float(chi2), "dof": len(E) - grau - 1, "p_gof": float(stats.chi2.sf(chi2, len(E) - grau - 1))}
            if grau == 2:
                res[nome_g]["dPdt"] = 2 * coef[0] / P * S_POR_ANO; res[nome_g]["sdPdt"] = 2 * np.sqrt(cov[0, 0]) / P * S_POR_ANO
        sep = (res["lin"]["t_pred_bjd"] - res["quad"]["t_pred_bjd"]) * 1440.0
        res["separacao_lin_quad_min"] = sep; res["separacao_sobre_sigma"] = sep / np.hypot(res["lin"]["sig_total_min"], res["quad"]["sig_total_min"])
        out[rot] = res
    return out


def por_alvo(tics, t26, D11, D9):
    """POS-HOC (sugestao do revisor interno, sem expectativa previa): o desenho com o (A, tau_c) PROPRIO de cada alvo
    (MV do bloco estendido, deriva_gp_26_hiper.parquet) em vez do agrupado - para separar 'correlacao' de 'agrupamento'
    na contagem de sobreviventes (o agrupado fica abaixo do sigma_j proprio de 144194304/424461577 e acima do de 230386284)."""
    H = pd.read_parquet(BASE / "deriva_gp_26_hiper.parquet")
    linhas = []
    for nome in ("matern32", "exp"):
        for t in tics:
            h = H[(H.kernel == nome) & (H.tic == t)].iloc[0]
            j = json.loads((BASE / f"oc__{t}.json").read_text(encoding="utf-8"))
            pontos = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
            r = refazer_gp(pontos, float(j["P_escada_d"]), nome, float(h.A_min), float(h.tau_c_d))
            det = r["p_curv"] < P_DET and r["p_adv"] < P_DET and r["vao_ok"]
            linhas.append({"kernel": nome, "tic": t, "em_D11": t in D11, "em_D9": t in D9, "A_proprio": float(h.A_min), "tau_c_proprio_d": float(h.tau_c_d),
                           "tau_c_na_borda": bool(h.tau_c_na_borda), **r, "detecta": det, "detecta_adequado": det and r["p_gof"] >= P_DET})
    R = pd.DataFrame(linhas); R.to_parquet(BASE / "deriva_gp_26_desenho_por_alvo.parquet", index=False)
    print("== POS-HOC: desenho com (A, tau_c) PROPRIO de cada alvo (sem expectativa previa)")
    out = {}
    for nome, g in R.groupby("kernel", sort=False):
        d11 = g[g.em_D11]; d9 = g[g.em_D9]
        print(f"  {nome}: D11 sobrevivem {int(d11.detecta.sum())} de 11 ({', '.join(str(t) for t in d11[d11.detecta].tic)}) | caem: {', '.join(str(t) for t in d11[~d11.detecta].tic)} | "
              f"D9 adequados {int(d9.detecta_adequado.sum())} de 9 | novas {int(g[~g.em_D11].detecta.sum())} | vao ambiguo {int((~g.vao_ok).sum())}")
        print(d11[["tic", "A_proprio", "tau_c_proprio_d", "dPdt", "sdPdt", "p_curv", "p_adv", "p_gof", "detecta"]].round(4).to_string(index=False))
        out[nome] = {"D11_sobrevivem": int(d11.detecta.sum()), "quais": sorted(int(t) for t in d11[d11.detecta].tic), "D9_adequados": int(d9.detecta_adequado.sum()), "novas": int(g[~g.em_D11].detecta.sum())}
    jj = json.loads((BASE / "deriva_gp_26.json").read_text(encoding="utf-8"))
    jj["pos_hoc_por_alvo"] = {"nota": "sugestao do revisor interno apos o Passo B; sem expectativa previa", **out}
    (BASE / "deriva_gp_26.json").write_text(json.dumps(jj, indent=1, default=float), encoding="utf-8")


if __name__ == "__main__":
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    tics = sorted(int(t) for t in t26.index)
    D11 = {t for t in tics if t26.loc[t, "curvatura"]}; D9 = {t for t in D11 if not t26.loc[t, "inadequada"]}
    if "--por-alvo" in sys.argv:
        por_alvo(tics, t26, D11, D9)
        sys.exit()
    S = pd.read_parquet(BASE / "jitter_primarios_26.parquet")
    blocos = {t: bloco_alvo(t, S) for t in tics}
    fora = [blocos[t] for t in tics if t not in D11]
    print(f"blocos: {sum(b['n'] for b in blocos.values())} epocas primarias ok em {len(tics)} alvos; fora do D11: {len(fora)} alvos, {sum(b['n'] for b in fora)} epocas")
    hiper, linhas_h = {}, []
    for nome in ("exp", "matern32"):
        h = ajustar_hiper(fora, nome)
        hiper[nome] = h
        print(f"== {nome}: agrupado fora do D11: A {h['A_min']:.3f} min [{h['A_ic68'][0]:.1f}, {h['A_ic68'][1]:.1f}], tau_c {h['tau_c_d']:.0f} d = {h['tau_c_anos']:.2f} a "
              f"[{h['tau_c_ic68_d'][0]:.0f}, {h['tau_c_ic68_d'][1]:.0f}]{' (borda)' if h['tau_c_na_borda'] else ''}; -2lnL {h['menos2lnL']:.1f} (A = 0: {h['menos2lnL_A0']:.1f})", flush=True)
        for t in tics:
            ht = ajustar_hiper([blocos[t]], nome)
            linhas_h.append({"kernel": nome, "tic": t, "em_D11": t in D11, "n": blocos[t]["n"], **{k: v for k, v in ht.items() if k not in ("n_alvos", "n_epocas")}})
            if t in D11:
                print(f"   D11 {t}: A {ht['A_min']:.2f} [{ht['A_ic68'][0]:.1f}, {ht['A_ic68'][1]:.1f}] tau_c {ht['tau_c_d']:.0f} d [{ht['tau_c_ic68_d'][0]:.0f}, {ht['tau_c_ic68_d'][1]:.0f}]"
                      f"{' (borda)' if ht['tau_c_na_borda'] else ''} | Delta(-2lnL) vs A = 0: {ht['menos2lnL_A0'] - ht['menos2lnL']:.1f}", flush=True)
    H = pd.DataFrame(linhas_h); H.to_parquet(BASE / "deriva_gp_26_hiper.parquet", index=False)
    print(f"   |Delta(-2lnL)| entre kernels (agrupado): {abs(hiper['exp']['menos2lnL'] - hiper['matern32']['menos2lnL']):.1f}")
    for nome in ("exp", "matern32"):
        g = H[H.kernel == nome]
        print(f"   {nome}: A por alvo - fora do D11 mediana {g[~g.em_D11].A_min.median():.2f}, D11 mediana {g[g.em_D11].A_min.median():.2f}; tau_c na borda em {int(g.tau_c_na_borda.sum())} de 26")
    # desenho
    linhas = []
    for nome in ("exp", "matern32"):
        h = hiper[nome]
        for t in tics:
            j = json.loads((BASE / f"oc__{t}.json").read_text(encoding="utf-8"))
            pontos = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
            r = refazer_gp(pontos, float(j["P_escada_d"]), nome, h["A_min"], h["tau_c_d"])
            det = r["p_curv"] < P_DET and r["p_adv"] < P_DET and r["vao_ok"]
            linhas.append({"kernel": nome, "tic": t, "em_D11": t in D11, "em_D9": t in D9, **r, "detecta": det, "detecta_adequado": det and r["p_gof"] >= P_DET,
                           "dPdt_original": float(t26.loc[t, "dPdt"]), "sdPdt_original": float(t26.loc[t, "s"])})
    R = pd.DataFrame(linhas); R.to_parquet(BASE / "deriva_gp_26_desenho.parquet", index=False)
    infl = json.loads((BASE / "jitter_primarios_26.json").read_text(encoding="utf-8"))
    print("\n== desenho com C = diag + GP(TESS):")
    for nome, g in R.groupby("kernel", sort=False):
        d11 = g[g.em_D11]; d9 = g[g.em_D9]
        print(f"  {nome} (A {hiper[nome]['A_min']:.2f} min, tau_c {hiper[nome]['tau_c_anos']:.2f} a): D11 sobrevivem {int(d11.detecta.sum())} de 11 ({', '.join(str(t) for t in d11[d11.detecta].tic)}) | "
              f"caem: {', '.join(str(t) for t in d11[~d11.detecta].tic)} | D9 adequados {int(d9.detecta_adequado.sum())} de 9 | novas {int(g[~g.em_D11].detecta.sum())} | vao ambiguo {int((~g.vao_ok).sum())} | sPxN max {g.sP_vezes_N.max():.3f}")
        print(d11[["tic", "dPdt", "sdPdt", "sdPdt_original", "p_curv", "p_adv", "p_gof", "sP_vezes_N", "detecta"]].round(4).to_string(index=False))
    print(f"  cenarios do Passo 1 (inflacao branca): A {infl['inflacao']['A: MV (informativos), IC95 sup (nao informativos)']['D11_sobrevivem']}, "
          f"B {infl['inflacao']['B: IC95 sup em todos']['D11_sobrevivem']}, C {infl['pos_hoc_quadratica']['D11_sobrevivem']} de 11")
    # 232634196
    prev = {nome: previsao_232634196(hiper[nome], nome) for nome in ("exp", "matern32")}
    for nome, pv in prev.items():
        for rot in ("7", "18"):
            q, l = pv[rot]["quad"], pv[rot]["lin"]
            print(f"  232634196 {nome} {rot} epocas, E = {E_PRED_232634196}: quad {q['t_pred_bjd']:.4f} +- {q['sig_total_min']:.2f} min (modelo {q['sig_modelo_min']:.2f}, GP {q['sig_gp_min']:.2f}, mu_gp {q['mu_gp_min']:+.2f}), "
                  f"dP/dt {q['dPdt']:+.4f} +- {q['sdPdt']:.4f}, p_gof {q['p_gof']:.3f}; lin {l['t_pred_bjd']:.4f} +- {l['sig_total_min']:.2f}; separacao {pv[rot]['separacao_lin_quad_min']:+.1f} min = {pv[rot]['separacao_sobre_sigma']:.1f} sigma")
    out = {"expectativa_commit": EXPECTATIVA, "hiper_agrupados": hiper, "desenho": {nome: {"D11_sobrevivem": int(g[g.em_D11].detecta.sum()), "quais": sorted(int(t) for t in g[g.em_D11 & g.detecta].tic),
                                                                                           "D9_adequados": int(g[g.em_D9].detecta_adequado.sum()), "novas": int(g[~g.em_D11].detecta.sum())} for nome, g in R.groupby("kernel", sort=False)},
           "previsao_232634196": prev}
    (BASE / "deriva_gp_26.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(f"  -> deriva_gp_26.json, deriva_gp_26_hiper.parquet, deriva_gp_26_desenho.parquet em {BASE}")
