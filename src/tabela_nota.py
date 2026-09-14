# -*- coding: utf-8 -*-
"""Tabelas 3 e 4 da nota, construidas dos JSON do lote por UM caminho de
codigo - antes disso a classificacao vivia num heredoc de sessao, e cada
mudanca de criterio era arqueologia. Regra da revisao: "cada correcao de
numero roda pela cadeia, nao e editada a mao; a Tabela 4 e recomputada da
Tabela 3 depois de qualquer mudanca".

CRITERIO DE ADEQUACAO DA PARABOLA (ponto 1.1 da revisao), declarado antes
de rodar: o limiar fixo chi2_red >= 2 tem probabilidade de falso alarme
que depende dos graus de liberdade (p = 0,157 com 1 gl; 0,092 com 4 gl).
Passa a valer o p do chi2 de bondade de ajuste com os gl do proprio alvo,
no MESMO nivel usado para a curvatura (p < 0,05) - o nivel e fixado por
coerencia com o outro teste da nota, nao pelo que ele marca. Declaro que
ja sei, do texto anterior, que V527 Dra tem chi2 = 6,54 com 3 gl (p =
0,088): sob p < 0,05 ele deixa de ser marcado, e a "validacao cega" da
versao anterior cai. O criterio nao e escolhido para o preservar. O que
sai da mudanca e o que sai. A contagem sob chi2_red >= 2 e reportada ao
lado, com o p que ela implica por gl, para que o leitor veja as duas.

Saidas: data/orquestra/oc_lote/tabela_26.parquet (regenerado),
reports/tabelas_nota.md (blocos markdown das Tabelas 3 e 4, teste de
sinal, lista dos marginais) - a nota copia dali, nao digita.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
from orquestra.contrato import Taxa  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
P_ADEQUACAO = 0.05      # bondade de ajuste da parabola, mesmo nivel da curvatura
P_CURVATURA = 0.05
CHI2R_ANTIGO = 2.0      # so para reportar ao lado


def nome_curto(va_nome, simbad_id):
    if isinstance(va_nome, str) and va_nome:
        m = re.match(r"^(V)0*(\d+) ([A-Z][a-z]{2})$", va_nome)
        if m:
            return f"V{m.group(2)} {m.group(3)}"
        if re.match(r"^[A-Z]{1,2} [A-Z][a-z]{2}$", va_nome):
            return va_nome
    return re.sub(r"\s+", " ", simbad_id.replace("V* ", "")).strip() if isinstance(simbad_id, str) else ""


def tabela3():
    lit = pd.read_parquet(BASE / "literatura_26.parquet").set_index("tic")
    alvos = pd.read_parquet(config.DATA / "orquestra" / "alvos_oc_88.parquet").set_index("tic")
    from astropy.coordinates import SkyCoord, get_constellation
    import astropy.units as u
    linhas = []
    for p in sorted(BASE.glob("oc__*.json")):
        x = json.loads(p.read_text(encoding="utf-8"))
        if x.get("status") != "ok" or "parabola" not in x:
            continue
        tic = int(x["tic"])
        chi2, dof = x["parabola"]["chi2"], x["parabola"]["dof"]
        p_gof = float(stats.chi2.sf(chi2, dof))
        forte = (x["p_curvatura"] < P_CURVATURA) and (x["adversarial"]["p_adversarial"] < P_CURVATURA)
        inad = p_gof < P_ADEQUACAO
        c = SkyCoord(float(alvos.loc[tic, "ra"]) * u.deg, float(alvos.loc[tic, "dec"]) * u.deg)
        linhas.append({
            "TIC": tic, "nome": nome_curto(lit.loc[tic, "va_nome"], lit.loc[tic, "simbad_id"]),
            "ra": float(alvos.loc[tic, "ra"]), "dec": float(alvos.loc[tic, "dec"]),
            "const": get_constellation(c, short_name=True), "tmag": float(alvos.loc[tic, "tmag"]),
            "P": x["P_ref_d"], "n": x["n_pontos"], "nT": x["n_aglomerados_tess"], "nS": x["n_temporadas_swasp"],
            "dof": dof, "span": (max(q["t0"] for q in x["pontos"]) - min(q["t0"] for q in x["pontos"])) / 365.25,
            "dPdt": x["parabola"]["dPdt_s_por_ano"], "s": x["parabola"]["sdPdt_s_por_ano"],
            "chi2": chi2, "chi2r": chi2 / dof, "p_gof": p_gof, "p": x["p_curvatura"],
            "padv": x["adversarial"]["p_adversarial"], "curvatura": forte, "inadequada": inad,
            "inadequada_chi2r2": chi2 / dof >= CHI2R_ANTIGO,
            "classe": ("curvature (parabola inadequate)" if forte and inad else
                       "curvature (parabola adequate)" if forte else
                       ("no significant curvature (fails adversarial" + ("; parabola inadequate)" if inad else ")"))
                       if x["p_curvatura"] < P_CURVATURA else
                       "no significant curvature (parabola inadequate)" if inad else
                       "no significant curvature"),
        })
    d = pd.DataFrame(linhas).sort_values("TIC").reset_index(drop=True)
    assert len(d) == 26, len(d)
    return d


def fmt_p(p):
    return "< 10⁻³" if p < 1e-3 else f"{p:.3f}"


def reinflado():
    """p_curv e p_adv sob barras infladas, por alvo: SuperWASP x fator medido
    (reinflar_swasp_26.parquet, cenario central) e TESS +1,0 / +2,6 min/ano
    (reinflar_tess_26.parquet)."""
    # a coluna "SuperWASP x fator" saiu em 2026-09-13: o fator medido (1,35) e o que o nulo com barras
    # re-estimadas preve (0,77) para o mesmo estimador nos 23 adequados; inflar por ele nao e correcao
    out = {}
    rt = BASE / "reinflar_tess_26.parquet"
    if rt.exists():
        R = pd.read_parquet(rt)
        out["tess"] = {(int(r.tic), float(r.cenario_min_por_ano)): r for r in R.itertuples()}
    return out or None


def _pp(r):
    return f"{fmt_p(r.p_curv)}, {fmt_p(r.p_adv)}{'*' if r.detecta else ''}"


def md_tabela3(d):
    R = reinflado() or {}
    extra = ""
    if "tess" in R:
        extra += " | TESS bars +1.0 min yr⁻¹: p_curv, p_adv | TESS bars +2.6: p_curv, p_adv"
    out = ["| TIC | Name | RA, Dec (J2000, deg) | Tmag | P (d) | N (TESS+SW) | d.o.f. | span (yr) | quadratic coefficient as dP/dt (s yr⁻¹) | χ²_red | p_gof | p_curv | p_adv | class" + extra + " |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|" + "---|" * (2 * ("tess" in R))]
    for x in d.itertuples():
        col = ""
        if "tess" in R:
            col += " | " + _pp(R["tess"][(x.TIC, 1.0)]) + " | " + _pp(R["tess"][(x.TIC, 2.6)])
        out.append(f"| {x.TIC} | {x.nome} | {x.ra:.4f}, {x.dec:+.4f} | {x.tmag:.1f} | {x.P:.4f} | {x.n} ({x.nT}+{x.nS}) | {x.dof} | "
                   f"{x.span:.1f} | {x.dPdt:+.4f} ± {x.s:.4f} | {x.chi2r:.2f} | {fmt_p(x.p_gof)} | {fmt_p(x.p)} | {fmt_p(x.padv)} | {x.classe}{col} |")
    return "\n".join(out)


def coluna4(d):
    z = d.dPdt.abs() / d.s
    surv = d.curvatura
    tx = Taxa.de(int(surv.sum()), len(d))
    pos, neg = int((d.dPdt > 0).sum()), int((d.dPdt < 0).sum())
    p_sinal = float(stats.binomtest(pos, pos + neg, 0.5).pvalue)
    return {
        "n": len(d), "med": d.dPdt.median(), "abs_med": d.dPdt.abs().median(),
        "q25": d.dPdt.quantile(.25), "q75": d.dPdt.quantile(.75), "min": d.dPdt.min(), "max": d.dPdt.max(),
        "sig_med": d.s.median(), "chi2r_med": d.chi2r.median(), "chi2r_q25": d.chi2r.quantile(.25), "chi2r_q75": d.chi2r.quantile(.75),
        "z2": int((z > 2).sum()), "z3": int((z > 3).sum()), "z5": int((z > 5).sum()),
        "p05": int((d.p < P_CURVATURA).sum()), "surv": int(surv.sum()), "surv_lo": tx.ic95[0], "surv_hi": tx.ic95[1],
        "surv_pos": int((surv & (d.dPdt > 0)).sum()), "surv_neg": int((surv & (d.dPdt < 0)).sum()),
        "pos": pos, "neg": neg, "p_sinal": p_sinal,
        "n_dof_lt1": int((d.dof < 1).sum()), "n_guard": int((d.dPdt.abs() > 10).sum()),
    }


SIGMA_PRECISO = 0.01    # s/ano: subconjunto limitado pela precisao, nao pelo valor


def md_tabela4(d):
    """Tres colunas: todos os 26; os adequados (p_gof); os bem medidos
    (sigma < 0,01 s/ano). A mediana de |dP/dt| muda com a particao - a
    particao por adequacao tira coeficientes pequenos bem medidos e deixa
    grandes mal medidos - e a tabela mostra isso em vez de esconder."""
    cols = [("all", d), (f"quadratic fit adequate, p_gof ≥ {P_ADEQUACAO}", d[~d.inadequada]),
            (f"adequate and σ(dP/dt) < {SIGMA_PRECISO} s yr⁻¹", d[(d.s < SIGMA_PRECISO) & ~d.inadequada])]
    cs = [coluna4(sub) for _, sub in cols]
    eps = BASE / "epsilon_26.parquet"
    if eps.exists():
        # dois alvos tem epsilon = 0 (barras SuperWASP grandes): a ausencia de
        # curvatura neles nao testa nada, e eles nao entram no denominador
        # honesto da fracao com curvatura
        e = pd.read_parquet(eps).set_index("tic").epsilon
        for (_, sub), c in zip(cols, cs):
            sens = sub[sub.TIC.map(e) > 0]
            k, n = int(sens.curvatura.sum()), len(sens)
            tx = Taxa.de(k, n)
            c["surv_sens"] = f"{k} of {n} ({k / n:.0%}, 95% CI {tx.ic95[0]:.0%}–{tx.ic95[1]:.0%})"
            c["n_sem_sens"] = len(sub) - n
    rf = BASE / "reinflar_tess_26.parquet"
    if rf.exists():
        # quantos de cada particao continuam detectados com as barras do
        # TESS reinfladas pela deriva entre setores (reinflar_tess_26.py)
        R = pd.read_parquet(rf)
        for (_, sub), c in zip(cols, cs):
            ks = [int(R[(R.cenario_min_por_ano == dd) & R.tic.isin(sub.TIC)].detecta.sum()) for dd in (1.0, 2.6)]
            c["surv_reinf"] = f"{ks[0]} / {ks[1]} of {len(sub)}"
    # as linhas "SuperWASP x1,65", "TESS x1,34", conjunta e teto conjunto sairam em 2026-09-13 (ver reinflado())

    def linha(q, f):
        return "| " + q + " | " + " | ".join(f(c) for c in cs) + " |"

    rows = [
        linha("dP/dt, median", lambda c: f"{c['med']:+.4f} s yr⁻¹"),
        linha("\\|dP/dt\\|, median", lambda c: f"{c['abs_med']:.4f} s yr⁻¹"),
        linha("dP/dt, interquartile range", lambda c: f"[{c['q25']:+.4f}, {c['q75']:+.4f}] s yr⁻¹"),
        linha("dP/dt, minimum and maximum", lambda c: f"[{c['min']:+.3f}, {c['max']:+.3f}] s yr⁻¹"),
        linha("Sign of dP/dt: positive / negative (two-sided sign test p)", lambda c: f"{c['pos']} / {c['neg']} (p = {c['p_sinal']:.2f})"),
        linha("σ(dP/dt), median", lambda c: f"{c['sig_med']:.4f} s yr⁻¹"),
        linha("χ²_red of quadratic fit, median (25th–75th)", lambda c: f"{c['chi2r_med']:.2f} ({c['chi2r_q25']:.2f}–{c['chi2r_q75']:.2f})"),
        linha("\\|dP/dt\\|/σ > 2 / > 3 / > 5", lambda c: f"{c['z2']} / {c['z3']} / {c['z5']}"),
        linha("Curvature with p < 0.05", lambda c: str(c["p05"])),
        linha("Curvature surviving the adversarial test",
              lambda c: f"{c['surv']} ({c['surv'] / c['n']:.0%}, 95% CI {c['surv_lo']:.0%}–{c['surv_hi']:.0%}) — {c['surv_pos']} > 0, {c['surv_neg']} < 0"),
        linha("Same, over targets with detection efficiency ε > 0 (Section 5.4)", lambda c: c.get("surv_sens", "—") if c.get("n_sem_sens", 0) else "—"),
        linha("Same, with TESS bars inflated by 1.0 / 2.6 min yr⁻¹ (Section 5.3, sensitivity ceiling)", lambda c: c.get("surv_reinf", "—")),
        linha("Measured but not tested (d.o.f. < 1)", lambda c: str(c["n_dof_lt1"])),
        linha("\\|dP/dt\\| > 10 s yr⁻¹ (cycle-count guard)", lambda c: str(c["n_guard"])),
    ]
    out = ["| Quantity | " + " | ".join(f"n = {c['n']} ({rot})" for (rot, _), c in zip(cols, cs)) + " |", "|---|---|---|---|"]
    out += rows
    return "\n".join(out), cs[1], cs[0], cs[2]


if __name__ == "__main__":
    d = tabela3()
    d.to_parquet(BASE / "tabela_26.parquet", index=False)
    t3 = md_tabela3(d)
    t4, a, b, c_prec = md_tabela4(d)
    inad = d[d.inadequada]
    inad_old = d[d.inadequada_chi2r2]
    marg = d[(d.p < P_CURVATURA) & (d.p > 0.01)]
    consts = d.const.value_counts()
    resumo = [
        f"Flag p_gof < {P_ADEQUACAO}: {len(inad)} of 26 flagged: " + ", ".join(f"{x.TIC} (χ²_red {x.chi2r:.2f}, {x.dof} dof, p {x.p_gof:.3f})" for x in inad.itertuples()),
        f"  of which with significant surviving curvature: {int(inad.curvatura.sum())}; without: {int((~inad.curvatura).sum())}",
        f"Old rule χ²_red ≥ 2: {len(inad_old)} of 26 flagged: " + ", ".join(f"{x.TIC} ({x.dof} dof, p_gof {x.p_gof:.3f})" for x in inad_old.itertuples()),
        f"  implied false-flag probability of χ²_red ≥ 2 by dof: " + ", ".join(f"{k} dof: {stats.chi2.sf(2 * k, k):.3f}" for k in sorted(d.dof.unique())),
        f"  chance that a flag of {len(inad_old)}/26 hits one given target by coincidence: {len(inad_old) / 26:.0%}",
        f"Marginal curvature (0.01 < p_curv < 0.05): " + ", ".join(f"{x.TIC} (p {x.p:.3f}, p_adv {x.padv:.3f})" for x in marg.itertuples())
        + f" — surviving adversarial: {int(marg.curvatura.sum())}",
        f"Constellations: " + ", ".join(f"{k} {v}" for k, v in consts.items()),
        "D11 (curvature surviving adversarial, all 26): " + ", ".join(str(t) for t in d[d.curvatura].TIC),
        "D9 (D11 with adequate fit): " + ", ".join(str(t) for t in d[d.curvatura & ~d.inadequada].TIC),
        "D11 minus D9 (surviving but flagged): " + ", ".join(str(t) for t in d[d.curvatura & d.inadequada].TIC),
        "Flagged without significant curvature: " + ", ".join(str(t) for t in d[~d.curvatura & d.inadequada].TIC),
        f"Column 3 of Table 4 (adequate and sigma < {SIGMA_PRECISO}): n={c_prec['n']}, surviving {c_prec['surv']}: "
        + ", ".join(str(t) for t in d[(d.s < SIGMA_PRECISO) & ~d.inadequada & d.curvatura].TIC),
        f"Sign test n={a['n']}: {a['pos']}+/{a['neg']}- p={a['p_sinal']:.2f}; n=26: {b['pos']}+/{b['neg']}- p={b['p_sinal']:.2f}",
        f"Surviving curvature: {a['surv']}/{a['n']} = {a['surv'] / a['n']:.0%} (CI {a['surv_lo']:.0%}–{a['surv_hi']:.0%}); {b['surv']}/26 = {b['surv'] / 26:.0%} (CI {b['surv_lo']:.0%}–{b['surv_hi']:.0%})",
        f"Surviving curvature over targets with epsilon > 0: all {b.get('surv_sens', '-')} | adequate {a.get('surv_sens', '-')} | col3 {c_prec.get('surv_sens', '-')}",
        f"Median |dP/dt| by partition: all {b['abs_med']:.4f} (n=26) | adequate {a['abs_med']:.4f} (n={a['n']}) | sigma<{SIGMA_PRECISO} {c_prec['abs_med']:.4f} (n={c_prec['n']}); "
        f"{b['abs_med']:.4f} s/yr = {b['abs_med'] / 86400:.2e} d/yr",
        f"Sign test sigma<{SIGMA_PRECISO}: {c_prec['pos']}+/{c_prec['neg']}- p={c_prec['p_sinal']:.2f}; surviving {c_prec['surv']}/{c_prec['n']}",
    ]
    md = "# Tabelas da nota (geradas por src/tabela_nota.py)\n\n## Tabela 3\n\n" + t3 + "\n\n## Tabela 4\n\n" + t4 + "\n\n## Resumo\n\n" + "\n".join("- " + s for s in resumo) + "\n"
    (config.ROOT / "reports" / "tabelas_nota.md").write_text(md, encoding="utf-8")
    # a nota copia dali, nao digita: os blocos entre marcadores sao substituidos
    nota = config.ROOT / "reports" / "dpdt_note.md"
    if nota.exists() and "--inserir" in sys.argv:
        txt = nota.read_text(encoding="utf-8")
        for rot, bloco in (("tabela3", t3), ("tabela4", t4)):
            ini, fim = f"<!-- {rot}:inicio -->", f"<!-- {rot}:fim -->"
            if txt.count(ini) != 1 or txt.count(fim) != 1:
                raise RuntimeError(f"marcadores {rot} ausentes ou duplicados na nota")
            txt = txt[:txt.index(ini) + len(ini)] + "\n" + bloco + "\n" + txt[txt.index(fim):]
        nota.write_text(txt, encoding="utf-8")
        print("nota: Tabelas 3 e 4 inseridas entre marcadores")
    print("\n".join(resumo))
    print(f"\n-> {BASE / 'tabela_26.parquet'} e reports/tabelas_nota.md")
