# -*- coding: utf-8 -*-
"""Compara o dP/dt dos 26 com o que os minimos da base de Brno (VarAstro,
herdeira do O-C Gateway) dao por si sos, nos alvos em que ha minimos
suficientes. E a verificacao externa que um referee faria: dados de outras
pessoas, outra cadeia, mesma estrela.

Entrada: data/orquestra/oc_lote/minimos_brno_26.parquet (literatura_26.py).

Regras da comparacao, fixadas antes de olhar os numeros:
  - so minimos primarios, so CCD/fotoeletrico (visual e fotografico sao
    outra classe de precisao e entram so como contagem);
  - E = arredondamento de (JD - m0)/P com a efemeride do proprio VarAstro
    (o O-C ja vem contra ela);
  - parabola por minimos quadrados sem pesos; a barra de cada ponto e
    1,4826 x MAD dos residuos da parabola (barra do proprio dado, como no
    lote), e a incerteza de dP/dt sai dela;
  - dP/dt [s/ano] = (2c / P) x 365,25 x 86400, c em dias/ciclo^2; se o
    P do VarAstro e P/2 do catalogo (V391 Dra), dP/dt e multiplicado por 2;
  - entra quem tem >= 4 primarios CCD e >= 10 anos de alavanca CCD.

EXPECTATIVA (antes de rodar): nos alvos em que a nossa parabola e adequada
e Brno tem >= 10 anos de CCD - V564 Dra (-0,0199 +- 0,0018), CV Dra
(-0,0006 +- 0,0021), BL Dra (+0,0098 +- 0,0049) - os dois dP/dt concordam
dentro de 2 sigma combinados. Nos marcados como parabola inadequada (V527
Dra, que tem LTTE publicado sem termo secular; V391 Dra; V353 Dra) nao ha
expectativa de concordancia: a parabola nao e o modelo. V584 Dra (8 anos)
nao testa. Desvio em qualquer direcao vai para a nota como esta.

Ressalva de independencia: a base de Brno pode conter minimos derivados do
TESS por terceiros (o campo "method" nao distingue); o relato conta quantos
minimos caem depois de 2018,5.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
MIN_PRIM_CCD = 4
MIN_ALAVANCA_ANOS = 10.0
S_POR_ANO = 365.25 * 86400.0


def ajustar_parabola(E, oc):
    A = np.vstack([np.ones_like(E), E, E ** 2]).T
    coef, *_ = np.linalg.lstsq(A, oc, rcond=None)
    res = oc - A @ coef
    sig = 1.4826 * np.median(np.abs(res - np.median(res)))
    if sig <= 0:
        raise RuntimeError("MAD zero: residuos degenerados")
    cov = sig ** 2 * np.linalg.inv(A.T @ A)
    return coef, np.sqrt(np.diag(cov)), res, sig


def _ano(jd):
    return 2000.0 + (np.asarray(jd) - 2451545.0) / 365.25


def _janela_nossa(tic):
    import json
    j = json.loads((config.DATA / "orquestra" / "oc_lote" / f"oc__{tic}.json").read_text(encoding="utf-8"))
    t = np.array([q["t0"] for q in j["pontos"]]) + 2457000.0
    return float(_ano(t.min())), float(_ano(t.max()))


def comparar(tic, m, P_cat, t0_cat_btjd, dpdt_nosso, s_nosso, classe):
    m = m[(m.tipo == "P")].copy()
    n_total = len(m)
    m = m[m.metodo.str.contains("ccd", case=False, na=False)]
    if len(m) < MIN_PRIM_CCD:
        return {"tic": tic, "va_nome": m.va_nome.iloc[0] if len(m) else None, "n_prim_total": n_total,
                "n_prim_ccd": len(m), "motivo": f"< {MIN_PRIM_CCD} primarios CCD"}
    # os "primarios" de Brno sao primarios no periodo do CATALOGO? Quando o
    # VarAstro dobra em P/2, metade deles e secundario, e a comparacao nao
    # e entre as mesmas coisas (V391 Dra: 5 de 7 em fase 0,5)
    fase = ((m.jd.values - (t0_cat_btjd + 2457000.0)) / P_cat) % 1.0
    frac_sec = float(np.mean(np.abs(fase - 0.5) < 0.1))
    if frac_sec > 0.2:
        return {"tic": tic, "va_nome": m.va_nome.iloc[0], "n_prim_total": n_total, "n_prim_ccd": len(m),
                "frac_sec_no_P_catalogo": frac_sec,
                "motivo": f"{frac_sec:.0%} dos 'primarios' de Brno em fase 0,5 no P do catalogo"}
    P_va, m0 = float(m.P_va.iloc[0]), float(m.m0_va.iloc[0])
    E = np.round((m.jd.values - m0) / P_va)
    oc = m.oc_d.values
    anos = (m.jd.max() - m.jd.min()) / 365.25
    if anos < MIN_ALAVANCA_ANOS:
        return {"tic": tic, "va_nome": m.va_nome.iloc[0], "n_prim_total": n_total, "n_prim_ccd": len(m),
                "alavanca_anos": anos, "motivo": f"alavanca CCD {anos:.1f} a < {MIN_ALAVANCA_ANOS:.0f}"}
    coef, err, res, sig = ajustar_parabola(E, oc)
    # a parabola descreve os minimos de Brno? Compara o residuo da parabola
    # com a dispersao entre minimos do MESMO ano (anos com >= 3 pontos):
    # se o residuo e bem maior que a dispersao intra-ano, a forma nao e
    # parabola - e o dP/dt de cada janela e um coeficiente local
    anos_pt = np.floor(2000.0 + (m.jd.values - 2451545.0) / 365.25)
    intra = [np.std(oc[anos_pt == a], ddof=1) for a in np.unique(anos_pt) if (anos_pt == a).sum() >= 3]
    sig_intra = float(np.sqrt(np.mean(np.square(intra)))) if intra else np.nan
    fator = 2.0 if abs(P_va * 2 - P_cat) < 0.01 * P_cat else 1.0
    if fator == 1.0 and abs(P_va - P_cat) > 0.01 * P_cat:
        raise RuntimeError(f"TIC {tic}: P do VarAstro {P_va} nao e P nem P/2 do catalogo {P_cat}")
    dpdt = fator * (2 * coef[2] / P_va) * S_POR_ANO
    s_dpdt = fator * (2 * err[2] / P_va) * S_POR_ANO
    z = (dpdt_nosso - dpdt) / np.hypot(s_nosso, s_dpdt)
    # a mesma parabola so nos minimos de Brno dentro da NOSSA janela (+-0,5 a):
    # se o dP/dt de Brno se move para o nosso, a diferenca e a janela sobre um
    # O-C que nao e parabola - nao a cadeia
    a0, a1 = _janela_nossa(tic)
    mj = m[(_ano(m.jd) >= a0 - 0.5) & (_ano(m.jd) <= a1 + 0.5)]
    jan = {"n_janela": len(mj), "dPdt_brno_janela": np.nan, "s_dPdt_brno_janela": np.nan, "z_janela": np.nan}
    if len(mj) >= MIN_PRIM_CCD:
        cj, ej, *_ = ajustar_parabola(np.round((mj.jd.values - m0) / P_va), mj.oc_d.values)
        dj, sj = fator * (2 * cj[2] / P_va) * S_POR_ANO, fator * (2 * ej[2] / P_va) * S_POR_ANO
        jan.update({"dPdt_brno_janela": dj, "s_dPdt_brno_janela": sj, "z_janela": (dpdt_nosso - dj) / np.hypot(s_nosso, sj)})
    return {"tic": tic, "va_nome": m.va_nome.iloc[0], "n_prim_total": n_total, "n_prim_ccd": len(m),
            "frac_sec_no_P_catalogo": frac_sec, **jan,
            "n_apos_2018_5": int((m.jd > 2458300).sum()), "alavanca_anos": anos, "dof": len(m) - 3,
            "sigma_ponto_min": sig * 1440, "rms_res_min": float(np.sqrt(np.mean(res ** 2))) * 1440,
            "max_res_min": float(np.max(np.abs(res))) * 1440, "sigma_intra_ano_min": sig_intra * 1440,
            "n_anos_com_3": len(intra), "dPdt_brno": dpdt, "s_dPdt_brno": s_dpdt,
            "dPdt_nosso": dpdt_nosso, "s_dPdt_nosso": s_nosso, "z": z, "classe_nossa": classe,
            "fator_P": fator, "motivo": "comparado"}


TT_UTC_D = 69.184 / 86400.0  # TDB(~TT) - UTC desde 2017, em dias; HJD-BJD (< 8 s) ignorado
JANELA_ANOS = 1.0


def sobrepor(tic, m, sig_ponto_d):
    """Poe as NOSSAS epocas (JSON do lote) na efemeride do VarAstro e compara
    cada uma com a media dos minimos CCD primarios de Brno a menos de 1 ano.
    Nao ha modelo aqui: e dado contra dado."""
    import json
    m = m[(m.tipo == "P") & m.metodo.str.contains("ccd", case=False, na=False)]
    P_va, m0 = float(m.P_va.iloc[0]), float(m.m0_va.iloc[0])
    j = json.loads((config.DATA / "orquestra" / "oc_lote" / f"oc__{tic}.json").read_text(encoding="utf-8"))
    p = pd.DataFrame(j["pontos"])
    jd_utc = p.t0.values + 2457000.0 - TT_UTC_D
    E = np.round((jd_utc - m0) / P_va)
    oc_n = (jd_utc - (m0 + P_va * E)) * 1440
    out = []
    for f, jd, o, s in zip(p.fonte, jd_utc, oc_n, p.sig_min):
        viz = m[np.abs(m.jd - jd) < JANELA_ANOS * 365.25]
        if not len(viz):
            out.append({"tic": tic, "fonte": f, "oc_nosso_min": o, "sig_nosso_min": s, "n_brno_1ano": 0})
            continue
        ob = viz.oc_d.mean() * 1440
        sb = sig_ponto_d * 1440 / np.sqrt(len(viz))
        out.append({"tic": tic, "fonte": f, "oc_nosso_min": o, "sig_nosso_min": s, "n_brno_1ano": len(viz),
                    "oc_brno_min": ob, "sig_brno_min": sb, "dif_min": o - ob, "z_ponto": (o - ob) / np.hypot(s, sb)})
    return pd.DataFrame(out)


if __name__ == "__main__":
    base = config.DATA / "orquestra" / "oc_lote"
    mins = pd.read_parquet(base / "minimos_brno_26.parquet")
    t26 = pd.read_parquet(base / "tabela_26.parquet").set_index("TIC")
    alvos = pd.read_parquet(config.DATA / "orquestra" / "alvos_oc_88.parquet").set_index("tic")
    linhas = []
    for tic, m in mins.groupby("tic"):
        r = t26.loc[tic]
        linhas.append(comparar(int(tic), m, float(alvos.loc[tic, "period_d"]), float(alvos.loc[tic, "t0_btjd"]),
                               float(r["dPdt"]), float(r["s"]), str(r["classe"])))
    df = pd.DataFrame(linhas)
    df.to_parquet(base / "comparacao_brno_26.parquet", index=False)
    pd.set_option("display.width", 250)
    cols = ["tic", "va_nome", "n_prim_total", "n_prim_ccd", "n_apos_2018_5", "alavanca_anos", "dof",
            "sigma_ponto_min", "rms_res_min", "max_res_min", "sigma_intra_ano_min", "n_anos_com_3",
            "dPdt_brno", "s_dPdt_brno", "dPdt_nosso", "s_dPdt_nosso", "z", "classe_nossa", "motivo"]
    print(df.reindex(columns=cols).to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    c = df[df.motivo == "comparado"]
    print(f"\ncomparados: {len(c)} de {len(df)} com minimos; |z| <= 2 em {(c.z.abs() <= 2).sum()}")
    print("\n== a mesma parabola de Brno restrita a nossa janela (+-0,5 a)")
    print(c[["tic", "va_nome", "n_prim_ccd", "dPdt_brno", "s_dPdt_brno", "n_janela", "dPdt_brno_janela",
             "s_dPdt_brno_janela", "dPdt_nosso", "s_dPdt_nosso", "z", "z_janela"]]
          .to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(f"  |z| <= 2 no registro inteiro: {(c.z.abs() <= 2).sum()} de {len(c)}; "
          f"na nossa janela: {(c.z_janela.abs() <= 2).sum()} de {c.z_janela.notna().sum()}")
    print("\n== sobreposicao dado contra dado (nossas epocas na efemeride do VarAstro; Brno a < 1 ano)")
    sob = pd.concat([sobrepor(int(r.tic), mins[mins.tic == r.tic], r.sigma_ponto_min / 1440) for r in c.itertuples()],
                    ignore_index=True)
    sob.to_parquet(base / "sobreposicao_brno_26.parquet", index=False)
    print(sob.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    com = sob[sob.n_brno_1ano > 0].copy()
    # separado por era: os minimos de Brno de 2019+ podem ter vindo do TESS
    # (o campo "method" nao distingue); so a era SuperWASP e independente
    com["era"] = np.where(com.fonte.str.startswith("TESS"), "TESS 2019+", "SuperWASP 2004-08")
    print(f"\nepocas nossas com vizinho de Brno a < 1 ano: {len(com)} de {len(sob)}")
    for era, g in com.groupby("era"):
        print(f"  {era}: n={len(g)}  |z| <= 2 em {(g.z_ponto.abs() <= 2).sum()}  mediana da diferenca "
              f"{g.dif_min.median():+.2f} min  rms de z {np.sqrt(np.mean(g.z_ponto ** 2)):.2f}")
    print("  (a diferenca mediana e da ordem da ambiguidade HJD_UTC/BJD_TDB de ~1 min: o que se compara e a dispersao)")
