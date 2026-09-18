# -*- coding: utf-8 -*-
"""TIC 232634196, Passo 0 da rodada 'quatro medidas': jitter entre setores
medido nos 11 setores 76-86 e levado as 18 epocas.

Expectativa em `expectativas/jitter_232634196.md`, commit cf0ef4e, escrita
ANTES de rodar. sigma_j por maxima verossimilhanca contra o modelo linear
LOCAL ajustado so aos 11 (nao contra a quadratica das 7): a cada sigma_j da
grade o ajuste ponderado por 1/(sigma_i^2 + sigma_j^2) e refeito e
-2 ln L = sum [ r^2/(s_i^2 + s_j^2) + ln(s_i^2 + s_j^2) ]; IC por perfil
(Delta 1 e 3,84). Depois as 18 epocas (7 do registro + 11) com sigma_j em
quadratura em TODAS as epocas TESS, SuperWASP inalterado: dP/dt, chi2/nu,
p_gof, Delta chi2 lin - quad, p_curv, p_adv; inclinacao interna dos 11
contra a quadratica das 7 e das 18 com as barras novas; e as 7 do desenho com
sigma_j nas 4 TESS, so para registro (o Passo 1 faz isso nos 26).
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
TIC = 232634196
S_POR_ANO = 365.25 * 86400.0
GRADE_MIN = np.arange(0.0, 10.0 + 1e-9, 0.005)


def menos2lnL(E, t_min, sig_min, sj):
    """Perfil em beta (linear em E) para um sigma_j fixo; t e sigma em minutos."""
    v = sig_min ** 2 + sj ** 2
    w = 1.0 / v
    X = np.vstack([np.ones_like(E), E]).T
    beta = np.linalg.solve(X.T @ (X * w[:, None]), X.T @ (w * t_min))
    r = t_min - X @ beta
    return float(np.sum(r ** 2 / v + np.log(v))), r


def perfil_sigma_j(E, t_min, sig_min, grade=GRADE_MIN):
    L = np.array([menos2lnL(E, t_min, sig_min, s)[0] for s in grade])
    i = int(np.argmin(L))
    dL = L - L[i]

    def ic(delta):
        ok = np.where(dL <= delta)[0]
        return float(grade[ok[0]]), float(grade[ok[-1]]), bool(ok[-1] < len(grade) - 1)
    lo68, hi68, _ = ic(1.0)
    lo95, hi95, f95 = ic(3.84)
    return {"sigma_j_min": float(grade[i]), "ic68": [lo68, hi68], "ic95": [lo95, hi95], "ic95_fecha": f95,
            "menos2lnL_min": float(L[i]), "menos2lnL_em_0": float(L[0]), "n": int(len(E)), "dof_linear": int(len(E) - 2)}


def refit(E, t, sig_d, P, t0_ref):
    """Linear e quadratica (cadeia_oc.ajustar), curvatura, adversarial, gof."""
    cl, _, chi2l, covl = co.ajustar(E, t, sig_d, 1)
    cq, _, chi2q, covq = co.ajustar(E, t, sig_d, 2)
    dof = len(E) - 3
    ag = pd.DataFrame({"E": E, "t0": t, "sig_min": sig_d * 1440.0})
    adv = co.adversarial(ag, P, float(t0_ref))
    return {"dPdt_s_por_ano": 2 * cq[0] / P * S_POR_ANO, "sdPdt_s_por_ano": 2 * np.sqrt(covq[0, 0]) / P * S_POR_ANO,
            "chi2_quad": float(chi2q), "dof": dof, "chi2_red": float(chi2q / dof), "p_gof": float(stats.chi2.sf(chi2q, dof)),
            "chi2_lin": float(chi2l), "delta_chi2_lin_quad": float(chi2l - chi2q),
            "p_curv": float(stats.chi2.sf(chi2l - chi2q, 1)), "p_adv": float(adv["p_adversarial"]),
            "coef_quad": cq.tolist(), "cov_quad": covq.tolist(), "coef_lin": cl.tolist(), "cov_lin": covl.tolist()}


def prever(coef, cov, E):
    x = np.array([E ** k for k in range(len(coef) - 1, -1, -1)], float)
    return float(x @ coef), float(np.sqrt(x @ cov @ x))


def inclinacao(ano, oc_min, sig_min):
    A = np.vstack([ano - ano.mean(), np.ones(len(ano))]).T
    w = 1.0 / sig_min ** 2
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    beta = cov @ (A.T @ (w * oc_min))
    return float(beta[0]), float(np.sqrt(cov[0, 0]))


if __name__ == "__main__":
    j = json.loads((BASE / f"oc__{TIC}.json").read_text(encoding="utf-8"))
    assert "/ 2" in j["barra_tess"], j["barra_tess"]
    novos = pd.read_parquet(BASE / "setores_76_86_232634196.parquet").sort_values("E").reset_index(drop=True)
    assert len(novos) == 11 and novos.setor.tolist() == list(range(76, 87)), novos.setor.tolist()
    P = float(j["P_escada_d"])
    # 1. sigma_j por MV contra a linear local aos 11 (em minutos, t relativo ao primeiro)
    E11 = novos.E.values.astype(float)
    t11_min = (novos.t0_btjd.values - novos.t0_btjd.values[0]) * 1440.0
    s11 = novos.sig_min.values
    perfil = perfil_sigma_j(E11, t11_min, s11)
    sj = perfil["sigma_j_min"]
    _, r_local = menos2lnL(E11, t11_min, s11, sj)
    print(f"TIC {TIC}: 11 setores 76-86, barras {s11.min():.2f}-{s11.max():.2f} min (mediana {np.median(s11):.2f})")
    print(f"  sigma_j (MV, linear local aos 11) = {sj:.3f} min; IC68 [{perfil['ic68'][0]:.3f}, {perfil['ic68'][1]:.3f}]; "
          f"IC95 [{perfil['ic95'][0]:.3f}, {perfil['ic95'][1]:.3f}] (fecha: {perfil['ic95_fecha']}); "
          f"-2lnL(0) - min = {perfil['menos2lnL_em_0'] - perfil['menos2lnL_min']:.2f}")
    print(f"  residuos locais (min): {np.round(r_local, 2).tolist()}; rms {np.sqrt(np.mean(r_local ** 2)):.2f}")
    # 2. as 7 do registro e as 18 com sigma_j nas TESS
    pts = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
    tess7 = pts.fonte.str.startswith("TESS").values
    E7, t7, s7 = pts.E.values.astype(float), pts.t0.values, pts.sig_min.values
    E18 = np.concatenate([E7, E11]); t18 = np.concatenate([t7, novos.t0_btjd.values])
    s18 = np.concatenate([s7, s11]); tess18 = np.concatenate([tess7, np.ones(11, bool)])
    s18_j = np.where(tess18, np.hypot(s18, sj), s18)
    s7_j = np.where(tess7, np.hypot(s7, sj), s7)
    t0_ref18 = float(t18[tess18].max()); t0_ref7 = float(t7[tess7].max())
    base7 = refit(E7, t7, s7 / 1440.0, P, t0_ref7)
    assert abs(base7["dPdt_s_por_ano"] - j["parabola"]["dPdt_s_por_ano"]) < 1e-9, "as 7 nao reproduzem o registro"
    sem_j18 = refit(E18, t18, s18 / 1440.0, P, t0_ref18)
    com_j18 = refit(E18, t18, s18_j / 1440.0, P, t0_ref18)
    com_j7 = refit(E7, t7, s7_j / 1440.0, P, t0_ref7)
    for rot, r in (("7 do registro (controle)", base7), ("18 sem jitter (30c1edf)", sem_j18),
                   ("18 COM sigma_j nas TESS", com_j18), ("7 COM sigma_j nas 4 TESS", com_j7)):
        print(f"  {rot:28s}: dP/dt {r['dPdt_s_por_ano']:+.4f} +- {r['sdPdt_s_por_ano']:.4f} s/ano | chi2/nu {r['chi2_quad']:.1f}/{r['dof']} = {r['chi2_red']:.2f}, "
              f"p_gof {r['p_gof']:.3f} | Delta chi2 lin-quad {r['delta_chi2_lin_quad']:.1f} | p_curv {r['p_curv']:.2e} p_adv {r['p_adv']:.2e}")
    # 3. inclinacao interna dos 11 com barras novas: contra a quadratica das 7 (fora da amostra) e das 18
    ano = 2000.0 + (novos.t0_btjd.values + 2457000.0 - 2451544.5) / 365.25
    incl = {}
    for rot, r in (("quad_7", base7), ("quad_18_com_j", com_j18), ("lin_7", None)):
        if r is None:
            coef, cov = np.array(base7["coef_lin"]), np.array(base7["cov_lin"])
        else:
            coef, cov = np.array(r["coef_quad"]), np.array(r["cov_quad"])
        pred = np.array([prever(coef, cov, e) for e in E11])
        oc = (novos.t0_btjd.values - pred[:, 0]) * 1440.0
        sb = np.hypot(np.hypot(s11, sj), pred[:, 1] * 1440.0)
        b, sb_ = inclinacao(ano, oc, sb)
        chi2 = float(np.sum((oc / sb) ** 2))
        incl[rot] = {"inclinacao_min_por_ano": b, "s": sb_, "oc_medio_min": float(oc.mean()), "chi2_11": chi2,
                     "p_gof_11": float(stats.chi2.sf(chi2, 11)), "oc_min": oc.tolist()}
        print(f"  O-C dos 11 contra {rot:14s} (barras hypot(s_i, s_j, pred)): media {oc.mean():+.2f} min, chi2 {chi2:.1f}/11 (p {incl[rot]['p_gof_11']:.3f}); "
              f"inclinacao {b:+.2f} +- {sb_:.2f} min/ano ({abs(b / sb_):.1f} sigma)")
    frase = "(b) quadratica consistente, dado o jitter medido; linear excluida" if com_j18["p_gof"] >= 0.05 and com_j18["p_curv"] < 0.05 \
        else "(a) nao polinomial no nivel do minuto" if com_j18["p_gof"] < 0.05 else "nenhuma das duas (linear nao excluida)"
    print(f"\n  FRASE: {frase}")
    out = {"tic": TIC, "expectativa_commit": "cf0ef4e", "perfil_sigma_j": perfil, "residuos_locais_min": r_local.tolist(),
           "barras_11_min": s11.tolist(), "refit": {"7_registro": base7, "18_sem_jitter": sem_j18, "18_com_jitter": com_j18, "7_com_jitter": com_j7},
           "inclinacao_11": incl, "frase": frase}
    (BASE / "jitter_232634196.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(f"  -> {BASE / 'jitter_232634196.json'}")
