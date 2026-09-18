# -*- coding: utf-8 -*-
"""Rodada nove, Passo L: TODO numero citado na nota, com o arquivo de origem.

Um numero no abstract ou numa legenda e o mais lido e o menos conferido. Este
script le os parquets/JSON da cadeia (estado E) e imprime, por secao, cada
grandeza que o texto cita, com a FONTE ao lado; grava `reports/numeros_da_nota.json`.
A conferencia final da nota compara o que esta escrito com o que sai daqui.

Nao mede nada: so le e formata. Se um arquivo faltar, falha alto.
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
import deriva_gp_26 as GP  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
OUT = config.ROOT / "reports" / "numeros_da_nota.json"


def jl(nome):
    return json.loads((BASE / nome).read_text(encoding="utf-8"))


def pq(nome):
    return pd.read_parquet(BASE / nome)


def n_efetivo(tic="230386284"):
    """n efetivo do bloco TESS sob o GP agrupado: sigma_i^2 * (1^T C^-1 1). Com C diagonal da n."""
    import caso_completo_d9 as H
    import jitter_primarios_26 as J
    j = jl(f"oc__{tic}.json")
    P, sP, T0, E_ult = J.escada_do_registro(j)
    S = pq("jitter_primarios_26.parquet")
    d = J.classificar_setores(S[(S.tic == int(tic)) & (S.minimo == "primario")], P, sP, T0, E_ult, 0.0)
    d = d[d.ok]
    pp = pd.DataFrame({"fonte": ["TESS s%d" % s for s in d.setor], "t0": d.t0_btjd.values,
                       "sig_min": d.sig_min.values, "E": d.E.values}).sort_values("t0").reset_index(drop=True)
    hip = jl("deriva_gp_26.json")["hiper_agrupados"]["matern32"]
    s2 = float(np.median(pp.sig_min.values / 1440.0)) ** 2
    out = {}
    for rot, A in (("diagonal", 0.0), ("agrupado", hip["A_min"])):
        C, _ = GP.cov_desenho(pp, "matern32", A, hip["tau_c_d"])
        um = np.ones(len(pp))
        out[rot] = float(s2 * (um @ np.linalg.solve(C, um)))
    out["n"] = int(len(pp))
    return out


def temporadas_tess(tic="230386284", fator=3.0):
    """Quantos blocos de observacao o bloco TESS tem, com 'bloco' = setores separados por mais de 3 tau_c."""
    import caso_completo_d9 as H  # noqa: F401
    import jitter_primarios_26 as J
    j = jl(f"oc__{tic}.json")
    P, sP, T0, E_ult = J.escada_do_registro(j)
    S = pq("jitter_primarios_26.parquet")
    d = J.classificar_setores(S[(S.tic == int(tic)) & (S.minimo == "primario")], P, sP, T0, E_ult, 0.0)
    t = np.sort(d[d.ok].t0_btjd.values)
    tc = jl("deriva_gp_26.json")["hiper_agrupados"]["matern32"]["tau_c_d"]
    cortes = int((np.diff(t) > fator * tc).sum())
    return {"n_setores": int(len(t)), "n_blocos": cortes + 1, "tau_c_d": tc, "criterio": f"gap > {fator} tau_c"}


def numeros():
    t26 = pq("tabela_26.parquet")
    D11 = sorted(int(t) for t in t26.TIC[t26.curvatura])
    D9 = sorted(int(t) for t in t26.TIC[t26.curvatura & ~t26.inadequada])
    qb = pq("qual_barra_26.parquet")
    nulo = pq("nulo_barras_reestimadas_26.parquet")
    nulo_c = pq("nulo_condicionado_26.parquet")
    nulo_g = pq("nulo_global_26.parquet")
    cs = pq("completude_secular_26_resumo.parquet")
    el = jl("epsilon_longo_26.json"); eg = jl("epsilon_gp_26.json")
    gp = jl("deriva_gp_26.json")
    ra = jl("ruido_admitido.json"); rad9 = pq("ruido_admitido_d9.parquet"); ra26 = pq("ruido_admitido_26.parquet")
    rec = pq("reconcilia_hi_d9.parquet"); var = pq("reconcilia_hi_d9_varredura.parquet")
    d9 = jl("caso_completo_d9.json"); sd = pq("caso_completo_d9_sem_desenho.parquet")
    c232 = jl("caso_232634196_completo.json")
    hj = jl("vies_hj_indep.json"); forma = jl("vies_forma_swasp.json"); tempo = jl("checa_tempo_swasp.json")
    sen = jl("sensib_vies_26.json"); ffi = jl("ffi_fora_amostra_d9.json")
    jit = jl("jitter_primarios_26.json"); sec = pq("secundarios_26.parquet")
    brno = pq("comparacao_brno_26.parquet")
    d9tab = rad9[rad9.nivel == 0.95].set_index("tic")
    b95 = ra26[ra26.nivel == 0.95]
    ne = n_efetivo()
    F = {}   # cada entrada: (valor, fonte)

    def put(sec_, chave, valor, fonte):
        F.setdefault(sec_, {})[chave] = {"valor": valor, "fonte": fonte}

    # --- amostra e cadeia
    put("2", "n_medidos", int(len(t26)), "tabela_26.parquet (linhas)")
    put("2", "n_setores_tess_cache", int(jit["n_setores"]) if "n_setores" in jit else int(pq("jitter_primarios_26.parquet").shape[0] / 2),
        "jitter_primarios_26.parquet: epocas primarias + secundarias")
    put("2", "ultimo_setor", int(max(x["setor_max"] for x in ffi["inventario"])), "ffi_fora_amostra_d9.json: inventario MAST")
    put("2", "n_produtos_depois_s86", int(ffi["n_ffi_depois_corte"]), "ffi_fora_amostra_d9.json")
    put("2", "ffi_tess_spoc_dt_mediana_min", float(np.median(np.abs(pq("ffi_fora_amostra_d9_controle.parquet").query("provenance == 'TESS-SPOC'").dt_ffi_menos_2min_min))),
        "ffi_fora_amostra_d9_controle.parquet")
    put("2", "ffi_qlp_dt_max_min", float(np.max(np.abs(pq("ffi_fora_amostra_d9_controle.parquet").query("provenance == 'QLP'").dt_ffi_menos_2min_min))),
        "ffi_fora_amostra_d9_controle.parquet")
    # --- 3.4 vies
    put("3.4", "vies_cadeia_min", co.VIES_CADEIA_MIN, "src/cadeia_oc.py VIES_CADEIA_MIN (estado E)")
    put("3.4", "vies_cadeia_sig", co.VIES_CADEIA_SIG, "src/cadeia_oc.py VIES_CADEIA_SIG")
    for k, v in hj["conjuntos"].items():
        put("3.4", f"hj[{k}]", {"media": v["media_min"], "sigma": v["sigma_min"], "chi2": v["chi2"], "dof": v["dof"],
                                "z_vs_zero": v["z_vs_zero"], "sem_WASP18": v["sem_WASP18_media"], "sem_WASP18_sigma": v["sem_WASP18_sigma"]},
            "vies_hj_indep.json: conjuntos")
    put("3.4", "forma_ponderado_min", forma["hj_ponderado"]["oc_ponderado_min"] if "oc_ponderado_min" in forma.get("hj_ponderado", {}) else forma["hj_ponderado"],
        "vies_forma_swasp.json: hj_ponderado")
    put("3.4", "tempo_dif_s", tempo.get("resumo", tempo), "checa_tempo_swasp.json")
    put("3.4", "sensib_excursao", sen["excursao_0_a_-4.28_diag"], "sensib_vies_26.json")
    # --- 4 distribuicao
    put("4", "dPdt_mediana", float(t26.dPdt.median()), "tabela_26.parquet")
    put("4", "abs_mediana", float(t26.dPdt.abs().median()), "tabela_26.parquet")
    put("4", "iqr", [float(t26.dPdt.quantile(.25)), float(t26.dPdt.quantile(.75))], "tabela_26.parquet")
    put("4", "barra_mediana", float(t26.s.median()), "tabela_26.parquet")
    put("4", "chi2r_mediano", float(t26.chi2r.median()), "tabela_26.parquet")
    put("4", "D11", D11, "tabela_26.parquet: curvatura")
    put("4", "D9", D9, "tabela_26.parquet: curvatura & ~inadequada")
    put("4", "flags_p_gof", sorted(int(t) for t in t26.TIC[t26.inadequada]), "tabela_26.parquet: inadequada")
    put("4", "sinais", {"pos": int((t26.dPdt > 0).sum()), "neg": int((t26.dPdt < 0).sum()),
                        "p": float(stats.binomtest(int((t26.dPdt > 0).sum()), len(t26), 0.5).pvalue)}, "tabela_26.parquet")
    # --- 5.1 / 5.2
    # 5.2: fatores por conjunto, recalculados da mesma forma que qual_barra_26 imprime (razao de somas com 1 - h)
    for f in ("SuperWASP", "TESS"):
        g = qb[qb.fonte == f]
        put("5.2", f"fator_{f}", float((g.z_corr ** 2).sum() / (1.0 - g.h).sum() * (1.0 - g.h).sum() / (1.0 - g.h).sum())
            if False else float((g.z ** 2).sum() / (1.0 - g.h).sum()), "qual_barra_26.parquet: soma z^2 / soma (1 - h)")
        ad = qb[(qb.fonte == f) & (qb.tic.isin(t26.TIC[~t26.inadequada]))]
        put("5.2", f"fator_{f}_23_adequados", float((ad.z ** 2).sum() / (1.0 - ad.h).sum()), "qual_barra_26.parquet, alvos com p_gof >= 0,05")
    put("5.2", "nulo", {f"{r.modo} | {r.fonte}": {"media": r.media, "mediana": r.mediana, "medido": r.medido, "percentil": r.percentil_do_medido}
                        for r in nulo.itertuples()}, "nulo_barras_reestimadas_26.parquet")
    put("5.2", "nulo_global", {"mediana_nulo": nulo_g.mediana_chi2r_nulo.iloc[0] if "mediana_chi2r_nulo" in nulo_g else None,
                               "mediana_medida": nulo_g.mediana_chi2r_medida.iloc[0] if "mediana_chi2r_medida" in nulo_g else None,
                               "agregado_nulo": nulo_g.agregado_nulo.iloc[0] if "agregado_nulo" in nulo_g else None,
                               "agregado_medido": nulo_g.agregado_medido.iloc[0] if "agregado_medido" in nulo_g else None},
        "nulo_global_26.parquet")
    ca = pq("completude_secular_26_alvo.parquet")
    put("5.2", "completude", {"A50_mediana": float(ca.A50.median()), "n_nunca_0.5_em_0.2": int((ca["completude_0.2"] < 0.5).sum()),
                              "supressao_mediana_0.2": float(ca["supressao_0.2"].median()), "supressao_max_0.2": float(ca["supressao_0.2"].max()),
                              "falso_alarme_mediano": float(ca.falso_alarme.median())},
        "completude_secular_26_alvo.parquet")
    put("5.2", "completude_curva", {str(r.dpdt): {"compl_re": r.compl_re_med, "sup": r.sup_med} for r in cs.itertuples()},
        "completude_secular_26_resumo.parquet")
    # --- 5.3.2 ruido entre setores
    put("5.3.2", "GP_agrupado", {"A_min": gp["hiper_agrupados"]["matern32"]["A_min"], "tau_c_d": gp["hiper_agrupados"]["matern32"]["tau_c_d"],
                                 "A_exp": gp["hiper_agrupados"]["exp"]["A_min"], "tau_c_exp_d": gp["hiper_agrupados"]["exp"]["tau_c_d"]},
        "deriva_gp_26.json: hiper_agrupados")
    put("5.3.2", "GP_desenho_agrupado_D11", gp["desenho"]["matern32"]["D11_sobrevivem"] if "desenho" in gp else None, "deriva_gp_26.json: desenho")
    put("5.3.2", "GP_desenho_por_alvo_D11", gp["pos_hoc_por_alvo"]["matern32"]["D11_sobrevivem"], "deriva_gp_26.json: pos_hoc_por_alvo")
    put("5.3.2", "admitido_D11", ra["B_D11_sobrevivem_95"], "ruido_admitido.json")
    put("5.3.2", "admitido_quais", ra["B_quais"], "ruido_admitido.json")
    put("5.3.2", "admitido_caem", sorted(set(D11) - set(ra["B_quais"])), "ruido_admitido.json + tabela_26")
    put("5.3.2", "admitido_modelo", ra["B_modelo_aplicado"], "ruido_admitido.json")
    put("5.3.2", "admitido_marcados", sorted(int(t) for t in b95[b95.marcado_menos_rejeitado].tic), "ruido_admitido_26.parquet")
    put("5.3.2", "chi2r_agrupado_nos_8", {str(int(r.tic)): float(r.GP_agrupado_chi2r) for r in rec.itertuples()}, "reconcilia_hi_d9.parquet")
    put("5.3.2", "n_efetivo", ne, "calculado aqui: sigma^2 (1^T C^-1 1) em oc__230386284")
    put("5.3.2", "temporadas_tess", {"gap_3tau": temporadas_tess(fator=3.0), "gap_1tau": temporadas_tess(fator=1.0),
                                     "gap_1tau_por_alvo": {str(t): temporadas_tess(str(t), 1.0)["n_blocos"] for t in D9 if t != 144194304},
                                     "gap_3tau_por_alvo": {str(t): temporadas_tess(str(t), 3.0)["n_blocos"] for t in D9 if t != 144194304}},
        "calculado aqui de jitter_primarios_26.parquet: blocos de setores separados por mais de 1 ou 3 tau_c")
    put("5.3.2", "sigma_j", {"primario_mediano_min": jit["sigma_j_prim_mediana_min"], "secundario_mediano_min": jit["sigma_j_sec_mediana_min"],
                             "quad_mediana_min": jit["pos_hoc_quadratica"]["sj_quad_mediana_informativos"], "n_setores": jit["n_setores_total"],
                             "n_informativos": jit["n_informativos"], "regra_parada": jit["regra_parada_3min_dispara"]}, "jitter_primarios_26.json")
    put("5.3.2", "variograma", jit["variograma_ajuste"], "jitter_primarios_26.json: variograma_ajuste")
    put("5.3.2", "inflacao_branca", jit["inflacao"], "jitter_primarios_26.json: inflacao")
    put("5.3.2", "secundarios_derivas", int(sec.veredito.str.contains("NAO andam juntos").sum()), "secundarios_26.parquet: veredito")
    put("5.3.2", "razao_sigma_agrupado_diagonal", {str(int(r.tic)): float(r.razao_s_agrupado_sobre_diagonal) for r in rec.itertuples()},
        "reconcilia_hi_d9.parquet")
    # --- 5.4 orbitas
    put("5.4", "faixas", el["faixas"], "epsilon_longo_26.json: faixas")
    put("5.4", "total", el["total"], "epsilon_longo_26.json: total")
    put("5.4", "tokovinin", el["tokovinin"], "epsilon_longo_26.json: tokovinin")
    put("5.4", "epsilon_gp", eg["kernels"]["matern32"], "epsilon_gp_26.json: kernels.matern32")
    # --- 5.5 janela
    put("5.5", "classes", {str(r.tic): r.classe for r in d9tab.reset_index().itertuples()}, "ruido_admitido_d9.parquet (nivel 0,95)")
    put("5.5", "z_faixa", {str(r.tic): [r.z_min, r.z_max] for r in d9tab.reset_index().itertuples()}, "ruido_admitido_d9.parquet")
    put("5.5", "DMD", {str(r.tic): r.DMD_3sigma_adm for r in d9tab.reset_index().itertuples()}, "ruido_admitido_d9.parquet")
    put("5.5", "so_tess_sem_desenho", {str(r.tic): [r.dPdt_sem, r.s_sem, r.n_tess_sem] for r in sd.itertuples() if r.testavel}, "caso_completo_d9_sem_desenho.parquet")
    put("5.5", "excursao_min", {str(a["tic"]): a["oc_excursao_min"] for a in d9["alvos"]}, "caso_completo_d9.json")
    put("5.5", "V527_completo", {m: [float(rec.set_index('tic').loc[424461577, f'tudo_{m}_dPdt']), float(rec.set_index('tic').loc[424461577, f'tudo_{m}_s'])]
                                 for m in ("diagonal", "GP_agrupado", "GP_A_proprio", "branco_proprio")}, "reconcilia_hi_d9.parquet")
    put("5.5", "c232_previsao", c232["4_previsao"], "caso_232634196_completo.json")
    put("5.5", "c232_pos_hoc", c232["pos_hoc"], "caso_232634196_completo.json")
    put("5.5", "c232_local", c232["efemeride_local_2024"], "caso_232634196_completo.json")
    put("5.5", "brno", {"V564_z": float(brno.set_index('tic').loc[329246824, 'z']), "V564_z_janela": float(brno.set_index('tic').loc[329246824, 'z_janela']),
                        "CVDra_z": float(brno.set_index('tic').loc[229687624, 'z']), "CVDra_z_janela": float(brno.set_index('tic').loc[229687624, 'z_janela'])},
        "comparacao_brno_26.parquet")
    put("5.5", "amplitude_varredura", {str(float(r.A_min)): [float(r.dPdt), float(r.s)] for r in var[var.tic == 230386284].itertuples() if not r.proprio},
        "reconcilia_hi_d9_varredura.parquet")
    # a 5.5 nao cita mais hash de commit (rodada doze: o hash so existe no repositorio privado
    # e nao e verificavel por um terceiro). O que ela cita e o nome do script; a verificacao e
    # que o docstring dele continue trazendo a previsao por classe declarada antes de rodar
    doc_j = (config.ROOT / "src" / "ffi_fora_amostra_d9.py").read_text(encoding="utf-8")
    assert "PREVISAO POR CLASSE" in doc_j, "ffi_fora_amostra_d9.py nao traz mais a previsao por classe no docstring"
    put("5.5", "J_previsao_script", "ffi_fora_amostra_d9", "src/ffi_fora_amostra_d9.py: docstring com PREVISAO POR CLASSE")
    return F


if __name__ == "__main__":
    F = numeros()
    for sec_, d in F.items():
        print(f"\n== {sec_}")
        for k, v in d.items():
            val = v["valor"]
            txt = json.dumps(val, default=float) if not isinstance(val, (int, float, str)) else (f"{val:.4f}" if isinstance(val, float) else str(val))
            print(f"   {k:<28s} {txt[:150]:<150s} <- {v['fonte']}")
    OUT.write_text(json.dumps(F, indent=1, default=float), encoding="utf-8")
    print(f"\n  -> {OUT}")
