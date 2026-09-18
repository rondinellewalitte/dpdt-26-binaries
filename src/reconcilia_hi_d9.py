# -*- coding: utf-8 -*-
"""Rodada nove, Passo K: reconciliacao H x I - o mesmo bloco TESS dos 8 testaveis
do D9 sob QUATRO modelos de ruido, para dizer de onde vem cada diferenca.

Modelos (todos com o mesmo estimador, o mesmo E0 e o mesmo GLS de H):
  1. "diagonal"      - C = diag(sigma_i^2), ou seja o GP com A = 0.
  2. "GP agrupado"   - Matern-3/2 com (A, tau_c) agrupados fora do D11
                       (deriva_gp_26.json: A 0,835 min, tau_c 67,7 d) - o de H e I.
  3. "GP A proprio"  - Matern-3/2 com o (A, tau_c) do PROPRIO alvo
                       (deriva_gp_26_hiper.parquet, MV no bloco estendido).
  4. "branco proprio"- C = diag(sigma_i^2 + sigma_j^2) com o sigma_j QUADRATICO do
                       proprio alvo (jitter_sigma_j_quad.parquet; o linear absorve
                       curvatura real e nao serve de ruido para um ajuste quadratico).
Para cada alvo e modelo: A usado, dP/dt so-TESS +- sigma, chi2/nu, z contra as 7
epocas. Mais a "diferenca minima detectavel" DMD = 3 sigma_TESS sob o GP agrupado
(a menor diferenca que o bloco TESS veria a 3 sigma) e o z nessa barra,
(so-TESS - 7 epocas) / DMD.

Duas diferencas nomeadas pelo autor, que este passo tem de explicar:
  a) 230386284: sigma 0,0008 -> 0,0110 (H pos-hoc com A proprio -> H/I com A
     agrupado). Teste do mecanismo: varredura de sigma(dP/dt) contra A com tau_c
     fixo. Se sigma cresce monotonicamente com A e a razao acompanha A/sigma_i,
     o mecanismo E a amplitude.
  b) 424461577: -0,0048 -> +0,0513. Teste: os dois numeros sao do MESMO modelo de
     ruido (GP agrupado) e de CONJUNTOS diferentes - "tudo" (TESS + SuperWASP) e
     "so-TESS". Se o -0,0048 se reproduzir no ajuste "tudo" sob todos os quatro
     modelos e o +0,0513 no "so-TESS" sob todos os quatro, a diferenca e de
     CONJUNTO (a alavanca de 16 anos entra ou nao), nao de amplitude.
REGRA DO AUTOR: se o mecanismo nao for a amplitude, PARAR e reportar (o Passo L
nao roda).

EXPECTATIVA (escrita e commitada ANTES de rodar):
  - Ordenacao das barras: diagonal <= branco proprio <= GP A proprio <= GP
    agrupado em quase todos; a excecao e 424461577 e 232634196, cujo A proprio
    (5,54 e 2,31 min) e MAIOR que o agrupado - neles o GP proprio da a maior barra.
  - Razao sigma(GP agrupado)/sigma(diagonal): esperada grande onde a barra por
    setor e pequena. Previsao quantitativa: a razao cresce com A/sigma_i mediano
    (Spearman > 0,8 nos 8); 230386284 (barra 0,10 min, A/sigma 8,4) x12-15;
    392536812 (barra 0,58) x2,5-3,5; 232634196 (barra 0,78) x2-3.
  - (a) 230386284: a varredura em A deve dar sigma(dP/dt) monotona e ~linear em A
    acima de A ~ sigma_i, reproduzindo 0,0008 em A = 0 e 0,0110 em A = 0,835
    dentro de 2%. Mecanismo = AMPLITUDE (confirmo).
  - (b) 424461577: espero que "tudo" fique em -0,004 a -0,006 nos quatro modelos e
    "so-TESS" em +0,03 a +0,25 nos quatro (o maior com A proprio, que engole o
    LTTE). Ou seja: mecanismo = CONJUNTO (alavanca SuperWASP), NAO amplitude.
    Por isso espero que a regra do autor DISPARE e o Passo L nao rode nesta
    rodada - a diferenca (b) esta explicada, mas nao pela amplitude.
  - DMD (3 sigma sob GP agrupado): 0,033 (230386284) a 0,079 (424461577) s/ano,
    e ~0,14 em 232634196; nenhum dos quatro "estaveis" tem DMD abaixo de 0,03,
    logo "nao contraditado" nao e "confirmado" em nenhum deles.
  - Nenhum veredito das Tabelas 3-5 muda aqui: este passo nao reajusta o desenho.

RESULTADO (2026-09-17, contra a expectativa acima; log logs_E/36_passo_K.txt):
  ORDENACAO das barras: como esperado - diagonal <= branco proprio <= GP A
  proprio <= GP agrupado, com 232634196 e 424461577 invertendo (A proprio 2,31 e
  5,54 min > 0,835). Razao sigma(agrupado)/sigma(diagonal) x2,3 a x13,6, e ela
  segue A/barra: Spearman +0,93 (p 0,001; esperado > 0,8). Valores previstos e
  medidos: 230386284 x13,6 (previ x12-15), 392536812 x2,9 (x2,5-3,5), 232634196
  x2,3 (x2-3).
  (a) 230386284 CONFIRMADO como AMPLITUDE: a varredura em A com tau_c fixo da
  sigma(dP/dt) = 0,00081 (A = 0) -> 0,00279 (0,2) -> 0,00534 (0,4) -> 0,00793
  (0,6) -> 0,01098 (0,835) -> 0,01573 (1,2) -> 0,02615 (2,0) -> 0,06527 (5,0):
  monotona e linear em A acima de A ~ sigma_i (0,10 min), reproduzindo os dois
  numeros do autor (0,0008 e 0,0110) dentro de 0,2%. dP/dt mal se move
  (+0,0091 -> +0,0088). Monotonia em A vale nos 8.
  (b) 424461577 NAO e amplitude - e CONJUNTO, e nem so: "tudo" (TESS + SW) da
  -0,0045 / -0,0048 / +0,0011 / -0,0066 nos quatro modelos (estavel em ~-0,005,
  como previ), mas "so-TESS" da -0,0244 / +0,0513 / +0,2430 / -0,0090: MUDA DE
  SINAL com o modelo de ruido (previ +0,03 a +0,25 nos quatro - errei). Ou seja,
  o par que o autor citou (-0,0048 -> +0,0513) separa-se em duas partes: a
  alavanca do SuperWASP (conjunto) e, dentro do so-TESS, uma sensibilidade
  enorme ao A (o LTTE de 2,73 a e ruido correlacionado para o GP).
  => A REGRA DO AUTOR DISPARA: o mecanismo de (b) nao e a amplitude. PARAR;
  o Passo L (reescrita) nao roda nesta rodada.
  DMD (3 sigma sob GP agrupado): 0,031 a 0,079 s/ano (0,142 em 232634196);
  previ 0,033-0,079 e ~0,14 - dentro. Nenhum "estavel" tem DMD < 0,03: "nao
  contraditado" nao e "confirmado" em nenhum deles.
  POS-HOC (sem expectativa previa, rotulado): escolhendo em cada alvo o modelo
  que os PROPRIOS dados calibram (chi2/nu mais proximo de 1), |z| > 3 em 4 de 8
  (198408416 +6,10; 392536812 -4,05; 390021728 +3,83; 232634196 -3,21) e os
  outros quatro ficam em |z| <= 2,93 (329246824 +2,93 e o de borda). O GP
  agrupado NAO e calibrado em nenhum alvo: chi2/nu 0,05-0,45 nos quietos e 2,16
  e 6,37 em 232634196 e V527 Dra. O z de cada alvo varia entre modelos por 1 a
  14 unidades, e em 377253090 ate o SINAL muda. Consequencia direta para a
  reescrita: a frase "contradiz 3 de 8 e nao contradiz 4" depende do modelo de
  ruido escolhido; com o modelo calibrado seriam 4 contraditos, entre eles um
  alvo do NUCLEO que o Passo I chamava de estavel (392536812).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import caso_completo_d9 as H  # noqa: E402
import cadeia_oc as co  # noqa: E402
import config  # noqa: E402
import jitter_primarios_26 as J  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
S_POR_ANO = H.S_POR_ANO
GRADE_A = (0.0, 0.2, 0.4, 0.6, 0.835, 1.2, 2.0, 3.0, 5.0)
EXPECTATIVA = "ver git log: commit da expectativa (rodada nove, Passo K)"
# classes do Passo I (a495c02), so para rotular a tabela - nao entram em nenhuma conta aqui
CLASSE_I = {329246824: "estavel", 392536812: "estavel", 230386284: "estavel", 377253090: "estavel",
            198408416: "ambiguo", 390021728: "contraditado", 232634196: "contraditado", 424461577: "contraditado"}


def blocos(tic):
    """so-TESS (todos os setores ok do cache), tudo (com SuperWASP), E0 e P - como em H."""
    j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
    P, sP, T0, E_ult = J.escada_do_registro(j)
    S = pd.read_parquet(BASE / "jitter_primarios_26.parquet")
    d = J.classificar_setores(S[(S.tic == tic) & (S.minimo == "primario")], P, sP, T0, E_ult, 0.0)
    d = d[d.ok]
    tess = pd.DataFrame({"fonte": ["TESS s%d" % s for s in d.setor], "t0": d.t0_btjd.values, "sig_min": d.sig_min.values,
                         "E": d.E.values, "setor": d.setor.values, "ano": H.ano(d.t0_btjd.values)}).sort_values("t0").reset_index(drop=True)
    pts = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
    sw = pts[pts.fonte.str.startswith("SuperWASP")][["fonte", "t0", "sig_min", "E"]].copy(); sw["setor"] = -1; sw["ano"] = H.ano(sw.t0.values)
    tudo = pd.concat([sw, tess], ignore_index=True).sort_values("t0").reset_index(drop=True)
    fontes_desenho = set(pts[pts.fonte.str.startswith("TESS")].fonte)
    ref = j["parabola"]
    return P, float(np.mean(tess.E.values)), tess, tudo, fontes_desenho, float(ref["dPdt_s_por_ano"]), float(ref["sdPdt_s_por_ano"])


def com_jitter(pp, sj):
    q = pp.copy(); q["sig_min"] = np.hypot(q.sig_min.values, sj); return q


if __name__ == "__main__":
    assert abs(co.VIES_CADEIA_MIN - (-1.26)) < 1e-9, "K usa o estado E"
    hip = json.loads((BASE / "deriva_gp_26.json").read_text(encoding="utf-8"))["hiper_agrupados"]["matern32"]
    A_G, TC_G = hip["A_min"], hip["tau_c_d"]
    HP = pd.read_parquet(BASE / "deriva_gp_26_hiper.parquet"); HP = HP[HP.kernel == "matern32"].set_index("tic")
    SJ = pd.read_parquet(BASE / "jitter_sigma_j_quad.parquet").set_index("tic")
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    tics = [t for t in H.D9 if t != 144194304]   # V Gru: 4 setores, nao testavel (Passo H)
    print(f"== PASSO K: {len(tics)} testaveis do D9, so-TESS sob quatro modelos de ruido (estado E)")
    print(f"   GP agrupado: A {A_G:.3f} min, tau_c {TC_G:.1f} d\n")
    linhas, varredura = [], []
    for tic in tics:
        P, E0, tess, tudo, f_des, d7, s7 = blocos(tic)
        A_p, tc_p = float(HP.loc[tic, "A_min"]), float(HP.loc[tic, "tau_c_d"])
        sj_q, sj_l = float(SJ.loc[tic, "sj_quad"]), float(SJ.loc[tic, "sj_lin"])
        bar = float(np.median(tess.sig_min.values))
        modelos = (("diagonal", 0.0, TC_G, tess), ("GP agrupado", A_G, TC_G, tess), ("GP A proprio", A_p, tc_p, tess),
                   ("branco proprio", 0.0, TC_G, com_jitter(tess, sj_q)))
        r = {"tic": tic, "nome": str(t26.loc[tic, "nome"]), "classe_I": CLASSE_I[tic], "n_tess": int(len(tess)),
             "barra_mediana_min": bar, "A_agrupado_min": A_G, "A_proprio_min": A_p, "tau_c_proprio_d": tc_p, "sj_quad_min": sj_q, "sj_lin_min": sj_l,
             "dPdt_7": d7, "s_7": s7}
        for nome, A, tc, pp in modelos:
            q = H.ajustar(pp, 2, "matern32", A, tc, P, E0)
            k = nome.replace(" ", "_")
            r[f"{k}_A_min"] = A if nome != "branco proprio" else sj_q
            r[f"{k}_dPdt"] = q["dPdt"]; r[f"{k}_s"] = q["sdPdt"]; r[f"{k}_chi2r"] = q["chi2r"]; r[f"{k}_p_gof"] = q["p_gof"]
            r[f"{k}_z"] = (q["dPdt"] - d7) / np.hypot(q["sdPdt"], s7)
        # "tudo" (com SuperWASP) nos mesmos quatro modelos - o teste (b)
        for nome, A, tc, pp in (("diagonal", 0.0, TC_G, tudo), ("GP agrupado", A_G, TC_G, tudo), ("GP A proprio", A_p, tc_p, tudo),
                                ("branco proprio", 0.0, TC_G, com_jitter(tudo, sj_q))):
            q = H.ajustar(pp, 2, "matern32", A, tc, P, E0)
            r[f"tudo_{nome.replace(' ', '_')}_dPdt"] = q["dPdt"]; r[f"tudo_{nome.replace(' ', '_')}_s"] = q["sdPdt"]
        r["DMD_3sigma"] = 3.0 * r["GP_agrupado_s"]
        r["z_DMD"] = (r["GP_agrupado_dPdt"] - d7) / r["DMD_3sigma"]
        r["razao_s_agrupado_sobre_diagonal"] = r["GP_agrupado_s"] / r["diagonal_s"]
        r["A_sobre_barra"] = A_G / bar
        linhas.append(r)
        # varredura em A (tau_c agrupado): o teste do mecanismo
        for A in GRADE_A + (A_p,):
            q = H.ajustar(tess, 2, "matern32", A, TC_G if A != A_p else tc_p, P, E0)
            varredura.append({"tic": tic, "A_min": A, "tau_c_d": TC_G if A != A_p else tc_p, "proprio": bool(A == A_p),
                              "dPdt": q["dPdt"], "s": q["sdPdt"], "chi2r": q["chi2r"]})
        print(f"  TIC {tic} {r['nome']:<18s} (n {r['n_tess']}, barra mediana {bar:.2f} min, A proprio {A_p:.2f}, sigma_j quad {sj_q:.2f}; 7 epocas {d7:+.4f} +- {s7:.4f})")
        for nome, _, _, _ in modelos:
            k = nome.replace(" ", "_")
            print(f"     {nome:<15s} A {r[k + '_A_min']:.3f} min: dP/dt {r[k + '_dPdt']:+.4f} +- {r[k + '_s']:.4f}  chi2/nu {r[k + '_chi2r']:.2f}  z {r[k + '_z']:+.2f}")
        print(f"     DMD (3 sigma, GP agrupado) {r['DMD_3sigma']:.4f} s/ano; z na barra DMD {r['z_DMD']:+.2f}; sigma(agrupado)/sigma(diagonal) x{r['razao_s_agrupado_sobre_diagonal']:.1f} (A/barra {r['A_sobre_barra']:.1f})")
        print(f"     tudo (TESS + SW): " + "; ".join(f"{n} {r['tudo_' + n.replace(' ', '_') + '_dPdt']:+.4f} +- {r['tudo_' + n.replace(' ', '_') + '_s']:.4f}" for n, _, _, _ in modelos))
    R = pd.DataFrame(linhas); V = pd.DataFrame(varredura)
    R.to_parquet(BASE / "reconcilia_hi_d9.parquet", index=False); V.to_parquet(BASE / "reconcilia_hi_d9_varredura.parquet", index=False)
    rho = stats.spearmanr(R.A_sobre_barra, R.razao_s_agrupado_sobre_diagonal)
    print(f"\n== mecanismo (a): Spearman(A/barra, sigma_agrupado/sigma_diagonal) = {rho.statistic:+.2f} (p {rho.pvalue:.3f}) nos {len(R)} alvos")
    print(R[["tic", "barra_mediana_min", "A_sobre_barra", "diagonal_s", "GP_agrupado_s", "razao_s_agrupado_sobre_diagonal"]].round(4).to_string(index=False))
    v0 = V[V.tic == 230386284].sort_values("A_min")
    print("\n   varredura em A de 230386284 (tau_c agrupado, exceto a linha do A proprio):")
    print(v0[["A_min", "tau_c_d", "proprio", "dPdt", "s", "chi2r"]].round(5).to_string(index=False))
    mono = bool(all((np.diff(V[(V.tic == t) & (~V.proprio)].sort_values("A_min").s.values) >= -1e-9).all() for t in tics))
    print(f"   sigma(dP/dt) monotona em A em todos os {len(tics)}: {mono}")
    rV = R.set_index("tic")
    print(f"\n== mecanismo (b) 424461577: so-TESS nos quatro modelos " +
          "; ".join(f"{n} {rV.loc[424461577, n.replace(' ', '_') + '_dPdt']:+.4f} +- {rV.loc[424461577, n.replace(' ', '_') + '_s']:.4f}"
                    for n in ("diagonal", "GP agrupado", "GP A proprio", "branco proprio")))
    print("   tudo (TESS + SW) nos quatro modelos: " +
          "; ".join(f"{n} {rV.loc[424461577, 'tudo_' + n.replace(' ', '_') + '_dPdt']:+.4f} +- {rV.loc[424461577, 'tudo_' + n.replace(' ', '_') + '_s']:.4f}"
                    for n in ("diagonal", "GP agrupado", "GP A proprio", "branco proprio")))
    # POS-HOC (sem expectativa previa): qual dos quatro modelos os PROPRIOS dados calibram (chi2/nu mais proximo
    # de 1 em log) e o que ele diz. O GP agrupado nao e calibrado em nenhum: cobre demais os quietos (chi2/nu 0,05-0,45)
    # e de menos os barulhentos (2,2 e 6,4).
    MOD = ("diagonal", "GP_agrupado", "GP_A_proprio", "branco_proprio")
    print("\n== POS-HOC (sem expectativa previa): o modelo calibrado pelos proprios dados (chi2/nu mais proximo de 1)")
    cal = []
    for r in linhas:
        m = min(MOD, key=lambda k: abs(np.log(r[k + "_chi2r"])) if r[k + "_chi2r"] > 0 else np.inf)
        cal.append({"tic": r["tic"], "nome": r["nome"], "classe_I": r["classe_I"], "modelo_calibrado": m, "A_min": r[m + "_A_min"],
                    "chi2r": r[m + "_chi2r"], "dPdt": r[m + "_dPdt"], "s": r[m + "_s"], "z": r[m + "_z"],
                    "z_min": min(r[k + "_z"] for k in MOD), "z_max": max(r[k + "_z"] for k in MOD)})
        c = cal[-1]
        print(f"   {c['tic']} {c['nome']:<18s} {c['modelo_calibrado']:<15s} (A {c['A_min']:.2f}, chi2/nu {c['chi2r']:.2f}): dP/dt {c['dPdt']:+.4f} +- {c['s']:.4f}, z {c['z']:+.2f}"
              f"  | z nos quatro modelos: {c['z_min']:+.2f} a {c['z_max']:+.2f} | classe do Passo I: {c['classe_I']}")
    ncont = sum(abs(c["z"]) > 3 for c in cal)
    print(f"   |z| > 3 no modelo calibrado: {ncont} de {len(cal)} ({sorted(int(c['tic']) for c in cal if abs(c['z']) > 3)}); "
          f"alvos em que o sinal de z muda entre modelos: {sorted(int(c['tic']) for c in cal if c['z_min'] * c['z_max'] < 0)}")
    out = {"expectativa": EXPECTATIVA, "pos_hoc_modelo_calibrado": cal, "gp_agrupado": {"A_min": A_G, "tau_c_d": TC_G},
           "spearman_A_sobre_barra_vs_razao_sigma": {"rho": float(rho.statistic), "p": float(rho.pvalue)},
           "sigma_monotona_em_A": bool(mono), "alvos": linhas}
    (BASE / "reconcilia_hi_d9.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(f"\n  -> {BASE / 'reconcilia_hi_d9.parquet'}, _varredura.parquet, .json")
