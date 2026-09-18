# -*- coding: utf-8 -*-
"""Item 2 da revisao externa: o vies de cadeia do SuperWASP (-2,14 +- 0,78 min)
entra hoje como ruido INDEPENDENTE por epoca (hypot(barra, 0,78) na diagonal).
Ele nao e: a incerteza do vies e um deslocamento COMUM de todo o bloco
arquival do alvo contra o bloco TESS - a forma degenerada com a curvatura -
e na diagonal a media de n temporadas o reduz por sqrt(n).

Aqui o ajuste e refeito por minimos quadrados GENERALIZADOS com a matriz
de covariancia certa: C = diag(sigma_interna^2) + 0,78^2 x J no bloco
SuperWASP (J = matriz de uns entre as temporadas do alvo; sigma_interna =
sqrt(barra^2 - 0,78^2), o que sobra da barra sem o termo do vies), e
diagonal no TESS (independente). Linear e quadratico com a mesma C;
Delta chi2 -> p_curv; chi2 do quadratico com N - 3 gl -> p_gof; adversarial
com a MESMA regra de cadeia_oc.adversarial (cada epoca 1 sigma marginal na
direcao que mais enfraquece a curvatura), avaliado com a C completa;
dP/dt e sigma pela covariancia GLS. Tudo o mais igual a reinflar_tess_26.
refazer. O que se compara e a tabela do estado B (diagonal) com a GLS.

EXPECTATIVA (escrita e commitada antes de rodar; a do autor, e a razao
dele): o bloco fora da diagonal so existe DENTRO do arquivo, logo o efeito
recai sobre os alvos cuja curvatura e carregada pelo arquivo, que e
essencialmente um - TIC 232634196: sigma(dP/dt) dele cresce (a media das 3
temporadas passa a ter 0,78 min de incerteza comum em vez de 0,78/sqrt 3)
e p_curv/p_adv continuam < 1e-3 (excursao de ~110 min contra barras de
20). Os outros nove de D10 sao carregados pelo TESS e NAO se movem:
mesmos vereditos (curvatura + adversarial), sigma(dP/dt) dentro de ~10%.
Nenhuma deteccao nova. Se algum dos nove se mexer de veredito, ha algo na
implementacao que nao e o vies correlacionado - parar e olhar antes de
escrever. O revisor, ao contrario, espera perda parcial nos apertados
(StKM 1-1676: 0,78 min comum contra ~1,4 min de excursao em 16 a): a
medida decide. Sanidade: com o bloco fora da diagonal zerado a GLS tem de
reproduzir a tabela do estado B numero a numero.
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

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
S_POR_ANO = 86400.0 * 365.25
P_DET = 0.05


def covariancia(sig_min, tess, correlacionado=True, piso_min=None):
    """C em dias^2: diagonal com as barras atribuidas; no bloco SuperWASP, o PISO da barra
    arquival sai da diagonal e entra como bloco comum.

    Ate o estado E o piso era so a incerteza do vies da cadeia (0,75 min). No estado F ele e
    hypot(0,75; sigma da media do vies de forma daquele alvo) - e as DUAS partes sao de modo
    comum (uma e a calibracao da cadeia, a outra e um numero unico por alvo), entao o piso inteiro
    vai para o bloco. `piso_min=None` mantem o comportamento antigo, para quem chama sem alvo."""
    sig = np.asarray(sig_min, float) / 1440.0
    C = np.diag(sig ** 2)
    if correlacionado:
        v = float(co.VIES_CADEIA_SIG if piso_min is None else piso_min) / 1440.0
        isw = np.where(~tess)[0]
        for i in isw:
            C[i, i] = max(sig[i] ** 2 - v ** 2, 1e-12)      # sigma interna
        C[np.ix_(isw, isw)] += v ** 2                        # bloco comum
    return C


def ajustar_gls(E, t, C, grau):
    X = np.vander(np.asarray(E, float), grau + 1)
    Ci = np.linalg.inv(C)
    cov = np.linalg.inv(X.T @ Ci @ X)
    coef = cov @ X.T @ Ci @ t
    res = t - X @ coef
    chi2 = float(res @ Ci @ res)
    return coef, res, chi2, cov


def adversarial_gls(E, t0, C):
    """A regra de cadeia_oc.adversarial, com a C completa: 1 sigma MARGINAL por epoca, na direcao que mais enfraquece."""
    sig = np.sqrt(np.diag(C))
    dchi = lambda tt: ajustar_gls(E, tt, C, 1)[2] - ajustar_gls(E, tt, C, 2)[2]
    base = dchi(t0)
    pior = t0.copy()
    for i in range(len(t0)):
        melhor_d, melhor = None, base                    # como em cadeia_oc.adversarial: cada passo compara com a BASE
        for sgn in (-1, 1):
            tt = pior.copy(); tt[i] += sgn * sig[i]
            d = dchi(tt)
            if d < melhor:
                melhor, melhor_d = d, sgn
        if melhor_d is not None:
            pior[i] += melhor_d * sig[i]
    return float(stats.chi2.sf(dchi(pior), 1))


def refazer_gls(pontos, P, correlacionado=True, piso_min=None):
    E = pontos.E.values.astype(float)
    t = pontos.t0.values.astype(float)
    tess = pontos.fonte.str.startswith("TESS").values
    C = covariancia(pontos.sig_min.values, tess, correlacionado, piso_min)
    _, _, chi2l, _ = ajustar_gls(E, t, C, 1)
    coef, _, chi2p, cov = ajustar_gls(E, t, C, 2)
    dof = len(E) - 3
    return {"dPdt": 2 * coef[0] / P * S_POR_ANO, "sdPdt": 2 * np.sqrt(cov[0, 0]) / P * S_POR_ANO,
            "chi2r": chi2p / dof, "p_gof": float(stats.chi2.sf(chi2p, dof)), "p_curv": float(stats.chi2.sf(chi2l - chi2p, 1)),
            "p_adv": adversarial_gls(E, t, C)}


if __name__ == "__main__":
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    D10 = set(t26.index[t26.curvatura])
    linhas = []
    for tic in sorted(t26.index):
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        pontos = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
        P = float(j["P_escada_d"])
        diag = refazer_gls(pontos, P, correlacionado=False)
        gls = refazer_gls(pontos, P, correlacionado=True, piso_min=co.vies_forma(tic)[1])   # estado F: o piso do alvo
        det_B = bool(t26.loc[tic, "curvatura"])
        det_diag = diag["p_curv"] < P_DET and diag["p_adv"] < P_DET
        det_gls = gls["p_curv"] < P_DET and gls["p_adv"] < P_DET
        # sanidade: a GLS sem bloco reproduz a cadeia (diagonal)
        assert abs(diag["dPdt"] - t26.loc[tic, "dPdt"]) < 1e-6 and abs(diag["sdPdt"] - t26.loc[tic, "s"]) < 1e-6, (tic, diag, t26.loc[tic, ["dPdt", "s"]])
        assert det_diag == det_B, (tic, "veredito diagonal difere da tabela", diag, det_B)
        assert abs(diag["p_adv"] - t26.loc[tic, "padv"]) < 1e-6 and abs(diag["p_curv"] - t26.loc[tic, "p"]) < 1e-6, (tic, "p da GLS diagonal difere da cadeia")
        linhas.append({"tic": int(tic), "em_D10": tic in D10, "dPdt_B": diag["dPdt"], "s_B": diag["sdPdt"], "p_curv_B": diag["p_curv"], "p_adv_B": diag["p_adv"],
                       "dPdt_gls": gls["dPdt"], "s_gls": gls["sdPdt"], "p_curv_gls": gls["p_curv"], "p_adv_gls": gls["p_adv"], "p_gof_gls": gls["p_gof"], "chi2r_gls": gls["chi2r"],
                       "detecta_B": det_B, "detecta_gls": det_gls, "razao_s": gls["sdPdt"] / diag["sdPdt"]})
        r = linhas[-1]
        print(f"TIC {tic} {'D10' if r['em_D10'] else '   '} | dP/dt {r['dPdt_B']:+.4f}+-{r['s_B']:.4f} -> {r['dPdt_gls']:+.4f}+-{r['s_gls']:.4f} (s x{r['razao_s']:.2f}) "
              f"| p_curv {r['p_curv_B']:.3g} -> {r['p_curv_gls']:.3g} | p_adv {r['p_adv_B']:.3g} -> {r['p_adv_gls']:.3g} | "
              f"{'sobrevive' if det_gls else ('CAI' if det_B else ('NOVA' if det_gls else ''))}", flush=True)
    R = pd.DataFrame(linhas)
    R.to_parquet(BASE / "gls_vies_26.parquet", index=False)
    d10 = R[R.em_D10]
    print(f"\n== resumo (expectativa: os 9 carregados pelo TESS nao mudam de veredito, s dentro de ~10%; 232634196 alarga e sobrevive; nenhuma nova)")
    print(f"  D10 sobrevivem: {int(d10.detecta_gls.sum())} de {len(d10)} | caem: {sorted(d10[~d10.detecta_gls].tic.tolist())} | novas: {sorted(R[~R.em_D10 & R.detecta_gls].tic.tolist())}")
    print(f"  razao s_gls/s_B: mediana {R.razao_s.median():.2f}; nos D10 sem 232634196: max {d10[d10.tic != 232634196].razao_s.max():.2f}; 232634196: {float(R[R.tic == 232634196].razao_s.iloc[0]):.2f}")
    print(f"  chi2_red mediano GLS {R.chi2r_gls.median():.2f}; p_gof < 0,05: {int((R.p_gof_gls < 0.05).sum())}")
    print(f"  -> {BASE / 'gls_vies_26.parquet'}")
