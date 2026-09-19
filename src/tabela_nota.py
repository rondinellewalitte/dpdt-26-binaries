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


# Tabela 2: o PRIMEIRO guard que cada um dos 62 nao medidos falhou, lido do registro por alvo (campo `erro`).
# A ordem das linhas e a da nota; o texto de "Nature" e fixo por categoria. Categoria fora da lista = falha alta.
GUARDS_TABELA2 = [
    ("<500", lambda m: "pontos validos" in m, "SuperWASP curve with < 500 valid points after cleaning", "archive depth; structural"),
    ("cobertura", lambda m: "cobertura por bloco" in m, "Phase coverage of the eclipse not testable in ≥ 1 season block (Section 3.3; a season whose fitted depth is ≤ 0 counts as not testable, Section 3.1)",
     "archive sampling; recoverable with a second archive"),
    ("download", lambda m: "veio vazio" in m, "SuperWASP light-curve download returned an empty file in four attempts (a different operation from the cone search of Table 1, step 4)",
     "archive-side failure; recoverable if the archive serves the file"),
    ("escada", lambda m: "escada nao fechou" in m, "Cycle count between TESS sectors ambiguous (period ladder did not close)",
     "local sector coverage; intermediate sectors were added from MAST and none was recovered (per-target outcome in the deposit)"),
    ("T14", lambda m: "T14 do catalogo inaplicavel" in m, "No applicable eclipse duration in the catalogue (4 NaN; 1 contact-binary width; 1 W UMa)",
     "catalogue gap; T14 was measured on the TESS curve for five and none was recovered (per-target outcome in the deposit)"),
]


def tabela2():
    """As linhas da Tabela 2 (markdown) a partir dos 88 registros; confere 62 e que nenhum motivo fica fora das categorias."""
    contagem = {k: [] for k, *_ in GUARDS_TABELA2}
    n_ok = 0
    for p in sorted(BASE.glob("oc__*.json")):
        x = json.loads(p.read_text(encoding="utf-8"))
        if x.get("status") == "ok":
            n_ok += 1
            continue
        m = str(x.get("erro", x.get("motivo", "")))
        cats = [k for k, f, *_ in GUARDS_TABELA2 if f(m)]
        if len(cats) != 1:
            raise RuntimeError(f"Tabela 2: TIC {x['tic']} com motivo fora das categorias ({cats}): {m[:120]}")
        contagem[cats[0]].append(int(x["tic"]))
    n_falha = sum(len(v) for v in contagem.values())
    assert n_ok == 26 and n_falha == 62, (n_ok, n_falha)
    linhas = ["| Reason | N | Nature |", "|---|---|---|"]
    for k, _, razao, natureza in GUARDS_TABELA2:
        linhas.append(f"| {razao} | {len(contagem[k])} | {natureza} |")
    return "\n".join(linhas)


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


CLASSE_CURTA = {"contraditado": "contradicted", "nao contraditado": "not contradicted",
                "ambiguo": "ambiguous", "LTTE publicado": "published LTTE"}


def janela():
    """O teste de janela (Seccao 5.5) e a sensibilidade ao vies (Seccao 3.4), por alvo.

    so-TESS SEM as epocas do desenho: caso_completo_d9_sem_desenho.parquet;
    faixa de z entre modelos admitidos, classe e DMD: ruido_admitido_d9.parquet (95%);
    Delta(dP/dt) por passo de 2,14 min no vies: sensib_vies_26.parquet (diagonal, sigma_v da cadeia).
    """
    out = {}
    try:
        sd = pd.read_parquet(BASE / "caso_completo_d9_sem_desenho.parquet").set_index("tic")
        ra = pd.read_parquet(BASE / "ruido_admitido_d9.parquet")
        ra = ra[ra.nivel == 0.95].set_index("tic")
        rec_gp = pd.read_parquet(BASE / "reconcilia_hi_d9.parquet").set_index("tic")
        rec_n = {int(t): int(rec_gp.loc[t, "n_tess"]) for t in rec_gp.index}
        sv = pd.read_parquet(BASE / "sensib_vies_26.parquet")
        sv = sv[(sv.metodo == "diagonal") & (sv.sigma_v_min == sv.sigma_v_min.min())]
        piv = sv.pivot(index="tic", columns="delta_min", values="dPdt")
        col0, col2 = 0.0, -2.14
        out["dvies"] = {int(t): abs(piv.loc[t, col0] - piv.loc[t, col2]) for t in piv.index}
        d9 = pd.read_parquet(BASE / "caso_completo_d9.parquet").set_index("tic")   # traz tambem os nao testaveis
        for t in d9.index:
            if not bool(d9.loc[t, "testavel"]):    # V Gru: as 4 epocas do cache SAO as do desenho
                out.setdefault("janela", {})[int(t)] = {"dPdt": None, "s": None, "n": int(d9.loc[t, "n_tess"]),
                                                        "z_lo": None, "z_hi": None, "classe": "not testable", "dmd": None}
        for t in sd.index:
            if not bool(sd.loc[t, "testavel"]):
                continue
            if int(t) == 424461577:
                # fora da CLASSIFICACAO (z e DMD nao se aplicam), mas a linha tem de ser
                # comensuravel com as outras: coeficiente e N do MESMO bloco e do MESMO modelo
                # admitido que as demais, e nao do ajuste sem as epocas do desenho (rodada
                # dezessete: a coluna trazia 36 epocas onde as outras traziam 40)
                iw = int(np.argmax(ra.loc[t, "s_adm"]))
                out.setdefault("janela", {})[int(t)] = {
                    "dPdt": float(ra.loc[t, "dPdt_adm"][iw]), "s": float(ra.loc[t, "s_adm"][iw]),
                    "n": int(rec_n.get(int(t), sd.loc[t, "n_tess_sem"])), "z_lo": None, "z_hi": None,
                    "gp_dPdt": float(rec_gp.loc[t, "GP_agrupado_dPdt"]), "gp_s": float(rec_gp.loc[t, "GP_agrupado_s"]),
                    "gp_chi2r": float(rec_gp.loc[t, "GP_agrupado_chi2r"]),
                    "classe": "published LTTE", "dmd": None}
                continue
            # COMENSURABILIDADE (rodada catorze, M3): a linha inteira vem do MESMO modelo. Ate
            # aqui o coeficiente saia do bloco so-TESS sem as epocas do desenho e o z saia dos
            # modelos admitidos - duas fontes na mesma linha, e uma delas (o GP agrupado, que a
            # 5.3.2 diz nao ser admitido em nenhum dos oito) na coluna principal. Agora
            # coeficiente, barra, z, classe e DMD saem todos dos admitidos, e o GP agrupado tem
            # coluna propria, rotulada.
            adm = list(ra.loc[t, "dPdt_adm"]), list(ra.loc[t, "s_adm"])
            i_pior = int(np.argmax(adm[1]))          # o menos favoravel: maior sigma (e o DMD)
            out.setdefault("janela", {})[int(t)] = {
                "dPdt": float(adm[0][i_pior]), "s": float(adm[1][i_pior]),
                "dPdt_lo": float(min(adm[0])), "dPdt_hi": float(max(adm[0])), "n_adm": len(adm[0]),
                "modelos": list(ra.loc[t, "admitidos"]), "n": int(rec_n.get(int(t), sd.loc[t, "n_tess_sem"])),
                "gp_dPdt": float(rec_gp.loc[t, "GP_agrupado_dPdt"]), "gp_s": float(rec_gp.loc[t, "GP_agrupado_s"]),
                "gp_chi2r": float(rec_gp.loc[t, "GP_agrupado_chi2r"]),
                "z_lo": float(ra.loc[t, "z_min"]), "z_hi": float(ra.loc[t, "z_max"]),
                "classe": CLASSE_CURTA.get(str(ra.loc[t, "classe"]), str(ra.loc[t, "classe"])),
                "dmd": float(ra.loc[t, "DMD_3sigma_adm"])}

    except FileNotFoundError:
        return {}
    return out


def bloco_minimo():
    """C2 da rodada treze: o menor bloco comum arquival (min) que derrota cada deteccao, da
    varredura de 0,75 a 10 min em `vies_forma_propagado_26.py`. Fora da grade -> "> 10"; alvo que
    nao e deteccao nao tem o que derrotar. A varredura roda DEPOIS de tabela_26.parquet existir,
    entao a primeira passada das tabelas pode nao ter o arquivo: a coluna sai como "--" e
    confere_nota acusa. A ordem canonica e tabela_nota -> vies_forma_propagado_26 -> tabela_nota."""
    arq = BASE / "vies_forma_propagado_varredura.parquet"
    if not arq.exists():
        print("  AVISO: varredura do bloco comum ausente; a coluna da Tabela 3 sai vazia "
              "(rode src/vies_forma_propagado_26.py e regere as tabelas)")
        return {}
    v = pd.read_parquet(arq).set_index("tic")
    if "em_D11" not in v.columns or "bloco_min_para_perder_deteccao" not in v.columns:
        # varredura de um estado anterior (esquema antigo, so os 8 do teste de janela): tratar
        # como ausente, senao a Tabela 3 sai com a coluna de outro estado da cadeia
        print("  AVISO: varredura com esquema antigo; a coluna da Tabela 3 sai vazia "
              "(rode src/vies_forma_propagado_26.py e regere as tabelas)")
        return {}
    out = {}
    for tic, r in v.iterrows():
        if not bool(r.em_D11):
            out[int(tic)] = None
            continue
        b = r.bloco_min_para_perder_deteccao
        out[int(tic)] = "> 10" if b is None or (isinstance(b, float) and np.isnan(b)) else f"{float(b):g}"
    return out


def md_tabela3(d):
    bm = bloco_minimo()
    out = ["| TIC | Name | RA, Dec (J2000, deg) | Tmag | P (d) | N (TESS+SW) | d.o.f. | span (yr) | quadratic coefficient as dP/dt (s yr⁻¹) | χ²_red | p_gof | p_curv | p_adv | class | archival offset that defeats it (min) |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for x in d.itertuples():
        b = bm.get(int(x.TIC), "--") if bm else "--"
        out.append(f"| {x.TIC} | {x.nome} | {x.ra:.4f}, {x.dec:+.4f} | {x.tmag:.1f} | {x.P:.4f} | {x.n} ({x.nT}+{x.nS}) | {x.dof} | "
                   f"{x.span:.1f} | {x.dPdt:+.4f} ± {x.s:.4f} | {x.chi2r:.2f} | {fmt_p(x.p_gof)} | {fmt_p(x.p)} | {fmt_p(x.padv)} | {x.classe} | "
                   f"{'—' if b is None else b} |")
    return "\n".join(out)


def md_tabela6(d):
    """Tabela 6: o teste de janela (Seccao 5.5) e a sensibilidade ao vies (Seccao 3.4), 26 linhas."""
    J = janela()
    jan, dv = J.get("janela", {}), J.get("dvies", {})
    out = ["| TIC | Name | block, admitted (s yr⁻¹) | N | z (adm.) | class | MDD | pooled GP (s yr⁻¹, χ²/ν) | Δ per 2.14 min |",
           "|---|---|---|---|---|---|---|---|---|"]
    for x in d.itertuples():
        w = jan.get(int(x.TIC))
        if w is None:
            col = "— | — | — | not in D9 | — | —"
        else:
            bloco = "—" if w["dPdt"] is None else f"{w['dPdt']:+.4f} ± {w['s']:.4f}"
            if w["z_lo"] is None:
                zz, dmd = "—", "—"
            else:
                zz = f"{w['z_lo']:+.1f}" if abs(w["z_hi"] - w["z_lo"]) < 0.05 else f"{w['z_lo']:+.1f} to {w['z_hi']:+.1f}"
                dmd = f"{w['dmd']:.4f}"
            gp = "—" if w.get("gp_dPdt") is None else f"{w['gp_dPdt']:+.4f} ± {w['gp_s']:.4f} ({w['gp_chi2r']:.2f})"
            col = f"{bloco} | {w['n'] if w['n'] else '—'} | {zz} | {w['classe']} | {dmd} | {gp}"
        out.append(f"| {x.TIC} | {x.nome} | {col} | {dv.get(int(x.TIC), float('nan')):.4f} |")
    return "\n".join(out)


GRADE_BLOCO = ["0.75", "1.0", "1.5", "2.0", "2.5", "3.0", "4.0", "5.0", "7.0", "10.0"]


def md_tabela7():
    """Tabela 6 (Apendice E): a varredura do bloco comum arquival, 26 linhas.

    |z| do teste de janela em cada bloco (— onde o alvo nao esta no teste), asterisco onde a
    deteccao do desenho ja nao sobrevive naquele bloco, e o menor bloco da grade que a derrota."""
    arq = BASE / "vies_forma_propagado_varredura.parquet"
    if not arq.exists():
        print("  AVISO: varredura ausente; Tabela 6 nao gerada (rode src/vies_forma_propagado_26.py)")
        return None
    v = pd.read_parquet(arq).set_index("tic")
    out = ["| TIC | class (window test) | " + " | ".join(f"{float(g):g}" for g in GRADE_BLOCO) + " | offset that defeats it |",
           "|---|---|" + "---|" * (len(GRADE_BLOCO) + 1)]
    for tic in sorted(v.index):
        r = v.loc[tic]
        # PODA (rodada quinze): a tabela impressa traz so as 11 deteccoes - as linhas dos 15 que
        # nao sao deteccao nao tem o que ser derrotado e nao entram no artigo; a varredura
        # completa dos 26 esta no deposito (vies_forma_propagado_varredura.parquet)
        if not bool(r["em_D11"]):
            continue
        cel = []
        for g in GRADE_BLOCO:
            z = r["z_por_bloco"].get(g)
            # o asterisco marca uma DETECCAO que nao sobrevive naquele bloco; num alvo que nao e
            # deteccao nao ha o que marcar (senao a linha inteira sai com asterisco sem sentido)
            marca = "" if (not r["em_D11"] or r["detecta_por_bloco"][g]) else "*"
            cel.append(("—" if z is None else f"{abs(float(z)):.1f}") + marca)
        b = r["bloco_min_para_perder_deteccao"]
        alvo = f"{float(b):g}" if not (b is None or pd.isna(b)) else ("> 10" if r["em_D11"] else "—")
        classe = CLASSE_CURTA.get(r["classe_hoje"], r["classe_hoje"]) if r["no_teste_de_janela"] else "—"
        out.append(f"| {tic} | {classe} | "
                   + " | ".join(cel) + f" | {alvo} |")
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
            c["surv_sens"] = f"{k} of {n}"
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
              lambda c: f"{c['surv']} of {c['n']} — {c['surv_pos']} > 0, {c['surv_neg']} < 0"),
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
    t6 = md_tabela6(d)
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
        f"Surviving curvature: {a['surv']} of {a['n']}; {b['surv']} of 26",
        f"Surviving curvature over targets with epsilon > 0: all {b.get('surv_sens', '-')} | adequate {a.get('surv_sens', '-')} | col3 {c_prec.get('surv_sens', '-')}",
        f"Median |dP/dt| by partition: all {b['abs_med']:.4f} (n=26) | adequate {a['abs_med']:.4f} (n={a['n']}) | sigma<{SIGMA_PRECISO} {c_prec['abs_med']:.4f} (n={c_prec['n']}); "
        f"{b['abs_med']:.4f} s/yr = {b['abs_med'] / 86400:.2e} d/yr",
        f"Sign test sigma<{SIGMA_PRECISO}: {c_prec['pos']}+/{c_prec['neg']}- p={c_prec['p_sinal']:.2f}; surviving {c_prec['surv']}/{c_prec['n']}",
    ]
    t2 = tabela2()
    # A Tabela 7 (varredura de offset, no apendice) tambem vai para o arquivo gerado: ate a rodada
    # dezessete ela era montada so na hora de inserir, e por isso a guarda de blocos gerados nao
    # tinha com que compara-la - uma edicao feita dentro dela na nota passaria batida ate a
    # proxima insercao a desfazer em silencio. E o defeito do item 22 dentro da guarda do item 22.
    t7 = md_tabela7()
    md = ("# Tabelas da nota (geradas por src/tabela_nota.py)\n\n## Tabela 2\n\n" + t2 + "\n\n## Tabela 3\n\n" + t3
          + "\n\n## Tabela 4\n\n" + t4 + "\n\n## Tabela 6\n\n" + t6 + "\n\n## Tabela 7\n\n" + (t7 or "(sem varredura)")
          + "\n\n## Resumo\n\n" + "\n".join("- " + s for s in resumo) + "\n")
    (config.ROOT / "reports" / "tabelas_nota.md").write_text(md, encoding="utf-8")
    # a nota copia dali, nao digita: os blocos entre marcadores sao substituidos
    nota = config.ROOT / "reports" / "dpdt_note.md"
    if nota.exists() and "--inserir" in sys.argv:
        txt = nota.read_text(encoding="utf-8")
        blocos = [("tabela2", t2), ("tabela3", t3), ("tabela4", t4), ("tabela6", t6)] + ([("tabela7", t7)] if t7 else [])
        for rot, bloco in blocos:
            ini, fim = f"<!-- {rot}:inicio -->", f"<!-- {rot}:fim -->"
            if txt.count(ini) != 1 or txt.count(fim) != 1:
                raise RuntimeError(f"marcadores {rot} ausentes ou duplicados na nota")
            txt = txt[:txt.index(ini) + len(ini)] + "\n" + bloco + "\n" + txt[txt.index(fim):]
        nota.write_text(txt, encoding="utf-8")
        print(f"nota: {len(blocos)} tabelas inseridas entre marcadores ({', '.join(r for r, _ in blocos)})")
    print("\n".join(resumo))
    print(f"\n-> {BASE / 'tabela_26.parquet'} e reports/tabelas_nota.md")
