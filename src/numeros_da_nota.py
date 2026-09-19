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
    # --- 4: as contagens da observacao (iv), que mudaram com o estado F
    put("4", "n_p_curv_abaixo_005", int((t26.p < 0.05).sum()), "tabela_26.parquet")
    put("4", "caem_no_adversarial", sorted(int(x) for x in t26.TIC[(t26.p < 0.05) & ~t26.curvatura]), "tabela_26.parquet")
    # --- 5.2.1: os agregados que a nota cita e que ate a rodada catorze nao eram registrados
    cal = jl("calibracao_chi2r.json")
    ngl = pq("nulo_global_26.parquet").iloc[0]
    put("5.2.1", "mediana_chi2r_26", float(cal["mediana_observada"]), "calibracao_chi2r.json (reinflar_tess_26)")
    put("5.2.1", "P_mediana_esperada", float(cal["P_maior_ou_igual"]), "calibracao_chi2r.json")
    put("5.2.1", "agregado_adequados", float(cal["agregado_adequados"]), "calibracao_chi2r.json")
    put("5.2.1", "nulo_mediana", {"mediana": float(ngl.mediana_chi2r_nulo), "pct": float(ngl.pct_mediana)}, "nulo_global_26.parquet")
    put("5.2.1", "nulo_agregado", {"mediana": float(ngl.agregado_nulo), "medido": float(ngl.agregado_medido),
                                   "pct": float(ngl.pct_agregado)}, "nulo_global_26.parquet")
    # --- 3.4 estado F: a correcao de forma por alvo, e o que ela moveu (rodada treze)
    vforma = json.loads((config.DATA / "orquestra" / "vies_forma_por_alvo.json").read_text(encoding="utf-8"))["alvos"]
    prevF = pd.DataFrame(json.loads((config.DATA / "orquestra" / "previsao_estado_F.json").read_text(encoding="utf-8"))["alvos"]).set_index("tic")
    varre = pq("vies_forma_propagado_varredura.parquet").set_index("tic")
    consF = json.loads((BASE / "vies_forma_propagado_26.json").read_text(encoding="utf-8"))
    d11 = sorted(int(x) for x in t26.TIC[t26.curvatura])   # t26 aqui nao esta indexado por TIC
    vabs = {int(k): abs(float(v["vies_min"])) for k, v in vforma.items()}
    put("3.4", "forma_max_D11", max(vabs[t] for t in d11), "vies_forma_por_alvo.json nos 11 do D11")
    put("3.4", "forma_max_amostra", max(vabs.values()), "vies_forma_por_alvo.json")
    put("3.4", "forma_tic_max", int(max(vabs, key=vabs.get)), "vies_forma_por_alvo.json")
    put("3.4", "forma_desloc_mediano", float(prevF[prevF.tem_forma].desloc_sigma.median()), "previsao_estado_F.json")
    put("3.4", "forma_desloc_max", float(prevF[prevF.tem_forma].desloc_sigma.max()), "previsao_estado_F.json")
    put("3.4", "forma_instaveis", sorted(int(k) for k, v in vforma.items() if v["instavel"]), "vies_forma_por_alvo.json")
    sobc = pq("sobreposicao_corrigida_26.parquet")
    sobc = sobc[(sobc.tic == 229687624) & (sobc.variante == "cadeia")].sort_values("fonte")
    put("3.4", "cvdra_dif_corrigida", [[float(r.dif_min), float(np.hypot(r.sig_nosso_min, r.sig_brno_min))] for r in sobc.itertuples()],
        "sobreposicao_corrigida_26.parquet (variante da cadeia)")
    # --- 5.3.2 e 5.5: o corte conservador e a varredura do bloco comum
    # completude e supressao no proprio coeficiente, por alvo do D11 (rodada treze: nao eram
    # registradas, e por isso a re-simulacao do estado F nao aparecia como divergencia)
    cse = pq("completude_secular_26_alvo.parquet").set_index("tic")
    gD11 = cse.loc[[t for t in d11 if t in cse.index]]
    put("5.2.5", "completude_D11", {"min": float(gD11.completude_no_proprio.min()), "max": float(gD11.completude_no_proprio.max()),
                                    "mediana": float(gD11.completude_no_proprio.median())}, "completude_secular_26_alvo.parquet")
    put("5.2.5", "supressao_D11", {"min": float(gD11.supressao_no_proprio.min()), "max": float(gD11.supressao_no_proprio.max()),
                                   "mediana": float(gD11.supressao_no_proprio.median())}, "completude_secular_26_alvo.parquet")
    put("5.2.5", "completude_0.02_ge_meio", int((cse["completude_0.02"] >= 0.5).sum()), "completude_secular_26_alvo.parquet")
    put("5.3.2", "conservadora_sobrevivem", len(consF["D11_sobrevive"]["conservadora"]), "vies_forma_propagado_26.json")
    put("5.3.2", "conservadora_quais", sorted(consF["D11_sobrevive"]["conservadora"]), "vies_forma_propagado_26.json")
    put("5.5", "bloco_minimo", {str(t): (None if varre.loc[t, "bloco_min_para_perder_deteccao"] is None
                                         else float(varre.loc[t, "bloco_min_para_perder_deteccao"])) for t in d11},
        "vies_forma_propagado_varredura.parquet")
    put("5.5", "z_bloco_3min", {str(t): float(varre.loc[t, "z_por_bloco"]["3.0"]) for t in varre.index
                                if varre.loc[t, "no_teste_de_janela"]}, "vies_forma_propagado_varredura.parquet")
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
    # ---- rodada dezessete: bloco B virado em limite, barras formais nos 26, GP agrupado, deflacao
    md = pq("media_diferenca_oc_d9.parquet").set_index("tic")
    # O conjunto que a 5.3.1 discute. Ate a rodada dezessete ele so existia digitado num patch, e
    # a frase que o apresentava ("the eight the window test classifies") descrevia outro conjunto.
    # Fica aqui, explicito, ate o humano decidir entre mante-lo e adotar a regra do D11 abaixo.
    OITO_DO_TEXTO = [198408416, 229914020, 230386284, 232634196, 329246824, 377253090, 390021728, 392536812]
    t8 = md.loc[OITO_DO_TEXTO]
    # A variante por regra: deteccoes (D11) com bloco completo nos dois minimos, sem a calibracao.
    # Difere do de cima por 198388252, deteccao sinalizada com os mesmos 39 setores de 229914020.
    _d11b = [t for t in sorted(int(x) for x in t26.TIC[t26.curvatura]) if t in md.index and t != 424461577]
    put("5.3.1", "anticorrelada_limite_D11",
        {"n": len(_d11b), "alvos": _d11b,
         "n_abaixo_do_vagar": sum(1 for t in _d11b if md.loc[t, "delta_pico_a_pico_95_min"] < md.loc[t, "amp_P_min"]),
         "n_r_positivo": sum(1 for t in _d11b if md.loc[t, "r_setor"] > 0),
         "n_anticorrelado_p05": sum(1 for t in _d11b if md.loc[t, "r_setor"] < 0 and md.loc[t, "p_r_setor"] < 0.05),
         "acrescentado": {"tic": 198388252, "r_setor": float(md.loc[198388252, "r_setor"]),
                          "limite_min": float(md.loc[198388252, "delta_pico_a_pico_95_min"]),
                          "vagar_min": float(md.loc[198388252, "amp_P_min"])}},
        "media_diferenca_oc_d9.parquet restrito ao D11 (variante por regra, nao adotada no texto)")
    abaixo = sorted(int(t) for t in t8.index if t8.loc[t, "delta_pico_a_pico_95_min"] < t8.loc[t, "amp_P_min"])
    # A fracao do vagar que o limite ocupa e o que diz o quanto ele exclui; a desigualdade sozinha
    # nao diz. E o vagar tem de ser o que um termo suave NAO explica (rodada dezenove): contra a
    # efemeride linear do bloco a excursao carrega a parabola do proprio bloco dentro de si.
    _fr = {int(t): 100.0 * float(t8.loc[t, "delta_pico_a_pico_95_min"]) / float(t8.loc[t, "amp_P_quad_min"])
           for t in t8.index}
    abaixo = sorted(t for t in _fr if _fr[t] < 100.0)
    put("5.3.1", "anticorrelada_limite",
        {"min": float(t8.delta_pico_a_pico_95_min.min()), "max": float(t8.delta_pico_a_pico_95_min.max()),
         "n": int(len(t8)), "n_abaixo_do_vagar": len(abaixo), "abaixo": abaixo,
         "fracao_min_pct": min(_fr.values()), "fracao_max_pct": max(_fr.values()),
         "abaixo_da_metade": sorted(t for t in _fr if _fr[t] < 50.0),
         "fracao_pct": {str(t): _fr[t] for t in _fr},
         "vagar_quad_min": {str(int(t)): float(t8.loc[t, "amp_P_quad_min"]) for t in t8.index},
         "vagar_quad_faixa": [float(t8.amp_P_quad_min.min()), float(t8.amp_P_quad_min.max())],
         "limite_faixa": [float(t8.delta_pico_a_pico_95_min.min()), float(t8.delta_pico_a_pico_95_min.max())],
         "queda_com_deriva_pct": {str(int(t)): 100.0 * (1 - float(t8.loc[t, "sigma_extra_D_95_com_deriva_min"])
                                                        / float(t8.loc[t, "sigma_extra_D_95_min"]))
                                  for t in t8.index if bool(t8.loc[t, "deriva_d_significativa"])},
         "por_alvo": {str(int(t)): [float(t8.loc[t, "delta_pico_a_pico_95_min"]), float(t8.loc[t, "amp_P_min"])]
                      for t in t8.index}},
        "media_diferenca_oc_d9.parquet: excesso branco de D, limite de perfil a 95%")
    put("5.3.1", "anticorrelada_V527",
        {"limite_pico_a_pico_min": float(md.loc[424461577, "delta_pico_a_pico_95_min"]),
         "vagar_min": float(md.loc[424461577, "amp_P_min"]),
         "vagar_quad_min": float(md.loc[424461577, "amp_P_quad_min"]),
         "fracao_pct": 100.0 * float(md.loc[424461577, "delta_pico_a_pico_95_min"]) / float(md.loc[424461577, "amp_P_quad_min"])},
        "media_diferenca_oc_d9.parquet: alvo de calibracao")

    # A passada com barras formais. O resumo cita DOIS numeros dela e o parecer os leu como se
    # fossem de passadas diferentes: os 31 sao quem chega a um ajuste com barras formais, os 26
    # sao os deste artigo - a mesma passada, dois recortes. Lidos do resumo versionado, nao dos
    # JSON por alvo (que ficam fora do repositorio).
    _bf = pq("../oc_lote_barras_formais/resumo_oc_lote.parquet")
    _bf = _bf[_bf.status == "ok"].set_index("tic")
    _n26 = [t for t in t26.TIC if t in _bf.index]

    def _conta(ix):
        s_ = _bf.loc[ix]
        return {"n": len(ix), "significativos": int(((s_.dPdt_s_por_ano.abs() / s_.sdPdt_s_por_ano) > 3).sum()),
                "sobrevivem_adversarial": int(((s_.p_curvatura < 0.05) & (s_.p_adversarial < 0.05)).sum())}

    put("5.1", "barras_formais_nos_26", _conta(_n26),
        "oc_lote_barras_formais/resumo_oc_lote.parquet restrito aos 26 da Tabela 3")
    put("5.1", "barras_formais_todos", _conta(list(_bf.index)),
        "oc_lote_barras_formais/resumo_oc_lote.parquet: todos os que chegam a um ajuste")

    # o GP agrupado como classificador: quantos contradiria (nao e sensibilidade; e comparacao)
    _r = rec.set_index("tic")
    _cl = [t for t in _r.index if t != 424461577]
    put("5.3.2", "GP_agrupado_classes",
        {"n": len(_cl), "contraditados": sum(1 for t in _cl if abs(float(_r.loc[t, "GP_agrupado_z"])) > 3)},
        "reconcilia_hi_d9.parquet: coluna GP_agrupado_z")

    # deflacao das barras do desenho de 230386284 (a barra do BLOCO domina a soma em quadratura)
    _t = 230386284
    _mod = list(rad9[rad9.nivel == 0.95].set_index("tic").loc[_t, "admitidos"])[0]
    _d7, _s7 = float(_r.loc[_t, "dPdt_7"]), float(_r.loc[_t, "s_7"])
    _db, _sb = float(_r.loc[_t, _mod + "_dPdt"]), float(_r.loc[_t, _mod + "_s"])
    _den = {k: float(np.hypot(_sb, _s7 / k)) for k in (1, 2, 3)}
    put("5.5", "deflacao_230386284",
        {"z": {str(k): (_db - _d7) / _den[k] for k in (1, 2, 3)}, "s_bloco": _sb, "s_desenho": _s7,
         "queda_denominador_pct": 100 * (1 - _den[3] / _den[1])},
        "reconcilia_hi_d9.parquet: desenho e bloco do modelo admitido")

    # rodada dezenove: a sobreposicao de CV Dra, que a 3.4 e a 5.5 escreviam com valores
    # diferentes (a 5.5 ficou na variante SEM correcao de forma, do estado E)
    _sob = pq("sobreposicao_corrigida_26.parquet")
    _cv = _sob[(_sob.tic == 229687624) & (_sob.variante == "cadeia")].sort_values("fonte")
    put("3.4", "sobreposicao_CVDra",
        {"dif_min": [float(x) for x in _cv.dif_min], "sig_min": [float(np.hypot(a, b)) for a, b in
                                                                 zip(_cv.sig_nosso_min, _cv.sig_brno_min)]},
        "sobreposicao_corrigida_26.parquet (variante da cadeia)")

    # rodada vinte: o estimador alternativo (largura do perfil livre) contra o da cadeia
    _alt = pq("estimador_alt_d11.parquet")
    _altg = pq("estimador_alt_d11_global.parquet").set_index("tic")
    _altok = _alt[_alt.dt_min.notna() & ~_alt.na_borda]
    put("3.4", "estimador_alt",
        {"n_alvos": int(len(_altg)), "n_temporadas": int(len(_alt)), "n_borda": int(_alt.na_borda.sum()),
         "razao_t14_mediana": float(_altg.razao.median()), "razao_t14_max": float(_altg.razao.max()),
         "n_mais_estreito": int((_altg.razao < 1).sum()),
         "dt_temporada_mediana": float(_altok.dt_min.abs().median()),
         "dt_temporada_p90": float(_altok.dt_min.abs().quantile(0.9)),
         "dt_temporada_max": float(_altok.dt_min.abs().max()),
         "dt_global_mediana": float(_altg.dt_global_min.abs().median()),
         "dt_global_max": float(_altg.dt_global_min.abs().max()),
         "dt_global_232634196": float(_altg.loc[232634196, "dt_global_min"])},
        "estimador_alt_d11.parquet e estimador_alt_d11_global.parquet")
    # o reajuste no espaco de epocas que a cadeia usa (uma por temporada), que substituiu a
    # comparacao por proxy contra a coluna de offset minimo
    _ref = pq("estimador_alt_refit.parquet").set_index("tic")
    put("3.4", "estimador_alt_refit",
        {"n": int(len(_ref)), "sobrevivem": int(_ref.sobrevive_alt.sum()),
         "sobrevivem_cadeia": int(_ref.sobrevive_cadeia.sum()),
         "desloc_sigma_mediana": float(_ref.desloc_sigma.median()),
         "desloc_sigma_max": float(_ref.desloc_sigma.max()),
         "pior_alvo": int(_ref.desloc_sigma.idxmax()),
         "pior_temporada_min": float(_alt[_alt.tic == int(_ref.desloc_sigma.idxmax())].dt_min.abs().max()),
         "cobertura_parcial": sorted(int(t) for t in _ref.index[_ref.cobertura_parcial])},
        "estimador_alt_refit.parquet")

    # rodada vinte: a amostra de calibradores do vies arquival, refeita em codigo e esgotada
    _cx = json.loads((BASE / "calibradores_expandido.json").read_text(encoding="utf-8"))
    put("3.4", "calibradores_expandido",
        {"n_candidatos": _cx["n_candidatos"], "n_elegiveis": _cx["n_elegiveis"],
         "n_rendem_epoca": _cx["n_rendem_epoca"], "n_passam": _cx["n_passam"],
         "IW22": _cx["conjuntos"]["IW22"], "ExoClock3": _cx["conjuntos"]["ExoClock3"]},
        "calibradores_expandido.json (35 candidatos com fonte SuperWASP; IW22 J/ApJS/259/62 e ExoClock III J/ApJS/265/4)")

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
