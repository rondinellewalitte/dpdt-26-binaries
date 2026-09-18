# -*- coding: utf-8 -*-
"""Rodada nove, Passo L (parte medida): o CRITERIO DE MODELO ADMITIDO, aplicado
(A) ao teste de janela so-TESS dos 8 testaveis do D9 e (B) ao desenho dos 26.

A REGRA E POSTERIOR AO PASSO K e usa SOMENTE chi2/nu - nenhum z entra na escolha
do modelo. Foi escrita pelo autor depois de ver que o GP agrupado nao e calibrado
em nenhum dos 8 (chi2/nu 0,05-6,4), e esta registrada assim para que a nota diga
que ela nao e cega.

  Modelo ADMITIDO: chi2/nu dentro do intervalo de 95% de chi2_nu/nu para o nu
  daquele alvo, isto e [chi2_nu(0,025)/nu, chi2_nu(0,975)/nu] (para nu ~ 36-40 o
  intervalo e ~[0,59; 1,51]).
  CLASSE por consenso entre os admitidos: |z| > 3 em TODOS -> "contraditado";
  |z| < 2 em TODOS -> "nao contraditado"; caso contrario -> "ambiguo".
  Se NENHUM for admitido: usar o MENOS rejeitado (menor |log(chi2/nu / centro)|
  em unidades do proprio intervalo) e marcar a linha.
  424461577 (V527 Dra) fica FORA do teste so-TESS: classe "LTTE publicado", e o
  numero reportado e o do ajuste COMPLETO (TESS + SuperWASP).
  SENSIBILIDADE: o mesmo com o intervalo em 90% e 99%; se alguma classe mudar,
  vai ao relatorio.

(A) usa os quatro modelos do Passo K (diagonal; GP agrupado A 0,835 / tau_c 67,7;
GP com A proprio; branco com sigma_j quadratico proprio), lidos de
reconcilia_hi_d9.parquet - nenhum ajuste novo, so a regra.

(B) 5.3.2: para CADA um dos 26, os dois modelos de ruido com parametro proprio
(branco com sigma_j quadratico; GP Matern-3/2 com A proprio) sao julgados pelo
MESMO criterio NO BLOCO TESS ESTENDIDO (todos os setores do cache, quadratica),
e o admitido e aplicado ao DESENHO (as 4-7 epocas das Tabelas 3-5): barras TESS
infladas por sigma_j, ou covariancia GP com o A do alvo. Dai sai a contagem de
sobreviventes do D11. Se os dois forem admitidos, vale o de chi2/nu mais proximo
de 1 (declarado); se nenhum, o menos rejeitado, marcado.

EXPECTATIVA (escrita e commitada ANTES de rodar):
(A) Com os chi2/nu do Passo K e nu = n - 3:
  - 329246824 (nu 39): so o branco proprio admitido (1,06); z +2,93 -> AMBIGUO.
  - 392536812 (nu 39): diagonal, GP proprio e branco admitidos (0,71); z -4,05
    nos tres -> CONTRADITADO.
  - 198408416 (nu 36): GP proprio (0,75) e branco (1,03); z +3,45 e +6,10 ->
    CONTRADITADO.
  - 230386284 (nu 35): NENHUM admitido (0,05-0,28); menos rejeitado = diagonal;
    z +0,97 -> NAO CONTRADITADO, marcado.
  - 377253090 (nu 37): GP proprio (0,98) e branco (0,97); z +0,27 e +0,04 ->
    NAO CONTRADITADO.
  - 390021728 (nu 38): GP proprio (1,06) e branco (1,10); z +3,83 e +5,50 ->
    CONTRADITADO.
  - 232634196 (nu 40): GP proprio (1,00) e branco (1,12); z -3,21 e -9,53 ->
    CONTRADITADO.
  - 424461577: LTTE publicado; ajuste completo -0,0045 a +0,0011 conforme o
    modelo, compativel com zero em todos.
  Contagem esperada: 4 contraditados, 2 nao contraditados (um marcado), 1
  ambiguo, 1 LTTE publicado. SENSIBILIDADE: espero NENHUMA mudanca de classe em
  90% nem em 99% - em 99% o GP agrupado entra em 377253090 (0,56) sem mudar a
  classe (z +0,29), e 230386284 continua sem nenhum admitido.
(B) Espero que o branco proprio seja admitido na maioria (o sigma_j quadratico
  foi ajustado no proprio bloco) e o GP proprio tambem; onde os dois entram, o
  vencedor deve ser o de chi2/nu mais proximo de 1, sem padrao fixo. Sobreviventes
  do D11 esperados: 6 a 8 (o cenario C branco deu 7 e o GP por alvo deu 7 no
  estado E), com o nucleo {198408416, 232634196, 329246824, 392536812} inteiro e
  377253090 dentro; 144194304 (4 setores, nu 1) entra com intervalo largo e
  qualquer modelo admitido - a linha dele fica marcada como nao informativa.
  Se a contagem sair fora de 5-9, ou se o nucleo quebrar, vai ao relatorio.

RESULTADO (2026-09-17, contra a expectativa acima; log logs_E/37_ruido_admitido.txt):
(A) EXATAMENTE o previsto, alvo a alvo: 4 CONTRADITADOS (392536812 z -4,05 com
  tres modelos admitidos; 198408416 +3,45 a +6,10; 390021728 +3,83 a +5,50;
  232634196 -9,53 a -3,21), 2 NAO CONTRADITADOS (230386284, sem nenhum modelo
  admitido - chi2/nu 0,05-0,28 - com o menos rejeitado, diagonal, e z +0,97;
  377253090, dois admitidos, z +0,04 a +0,27), 1 AMBIGUO (329246824: so o branco
  proprio admitido, z +2,93) e 1 LTTE PUBLICADO (424461577; pelo ajuste completo
  -0,0045 a +0,0011, compativel com zero). Sensibilidade: em 90% e em 99%
  NENHUMA classe muda (muda so o numero de admitidos: em 99% entram o GP
  agrupado em 377253090 e o GP proprio em 329246824 e V527, sem efeito).
  DMD pelo modelo admitido (3 sigma): 0,0024 (230386284) a 0,398 (232634196);
  nos nao contraditados 0,0024 e 0,0146 s/ano - o teste de 230386284 e apertado,
  o de 377253090 e duas vezes o proprio |dP/dt|.
(B) D11 sobrevivem 7 de 11 (144194304, 198408416, 230386284, 232634196,
  329246824, 377253090, 392536812) - dentro do 6-8 esperado; caem 198388252,
  229914020, 390021728 e 424461577; NUCLEO INTACTO; nenhuma nova fora do D11.
  Modelo aplicado: GP com A proprio em 15 dos 26, branco com sigma_j proprio em
  11; dois alvos sem nenhum admitido (230386284 e 329248002, ambos com chi2/nu
  0,17-0,28 no bloco estendido: barras por setor conservadoras demais) usam o
  menos rejeitado e ficam marcados. Em 90% e 99% nenhum veredito muda.
  Nota: os dois modelos coincidem quando sigma_j = A = 0 (11 alvos), e ai o
  criterio escolhe qualquer um deles - o resultado e o mesmo.
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
import deriva_gp_26 as GP  # noqa: E402
import jitter_primarios_26 as J  # noqa: E402
import reconcilia_hi_d9 as K  # noqa: E402
import reinflar_tess_26 as RT  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
MODELOS = ("diagonal", "GP_agrupado", "GP_A_proprio", "branco_proprio")
NIVEIS = (0.95, 0.90, 0.99)
P_DET = 0.05
NUCLEO = (198408416, 232634196, 329246824, 392536812)
LTTE_PUBLICADO = 424461577
EXPECTATIVA = "ver git log: commit da expectativa (rodada nove, Passo L: criterio de modelo admitido)"


def intervalo(nu, nivel):
    a = (1.0 - nivel) / 2.0
    return float(stats.chi2.ppf(a, nu) / nu), float(stats.chi2.ppf(1 - a, nu) / nu)


def admitido(chi2r, nu, nivel):
    lo, hi = intervalo(nu, nivel)
    return bool(lo <= chi2r <= hi), lo, hi


def distancia_rejeicao(chi2r, nu, nivel):
    """Quao fora do intervalo esta, em unidades da meia-largura em log - para 'o menos rejeitado'."""
    lo, hi = intervalo(nu, nivel)
    c = np.sqrt(lo * hi)
    return abs(np.log(chi2r / c)) / (0.5 * abs(np.log(hi / lo)))


def classe_por_consenso(zs):
    if not zs:
        return "sem modelo"
    if all(abs(z) > 3.0 for z in zs):
        return "contraditado"
    if all(abs(z) < 2.0 for z in zs):
        return "nao contraditado"
    return "ambiguo"


def parte_A(nivel, verboso=True):
    R = pd.read_parquet(BASE / "reconcilia_hi_d9.parquet").set_index("tic")
    linhas = []
    for tic in R.index:
        nu = int(R.loc[tic, "n_tess"]) - 3
        lo, hi = intervalo(nu, nivel)
        adm = [m for m in MODELOS if admitido(float(R.loc[tic, m + "_chi2r"]), nu, nivel)[0]]
        marcado = False
        if not adm:
            adm = [min(MODELOS, key=lambda m: distancia_rejeicao(float(R.loc[tic, m + "_chi2r"]), nu, nivel))]
            marcado = True
        zs = [float(R.loc[tic, m + "_z"]) for m in adm]
        if tic == LTTE_PUBLICADO:
            classe = "LTTE publicado"
        else:
            classe = classe_por_consenso(zs)
        d = {"tic": int(tic), "nome": R.loc[tic, "nome"], "nivel": nivel, "nu": nu, "int_lo": lo, "int_hi": hi,
             "admitidos": adm, "n_admitidos": len(adm), "marcado_menos_rejeitado": marcado,
             "z_min": float(min(zs)), "z_max": float(max(zs)), "classe": classe,
             "dPdt_adm": [float(R.loc[tic, m + "_dPdt"]) for m in adm], "s_adm": [float(R.loc[tic, m + "_s"]) for m in adm],
             "DMD_3sigma_adm": float(3.0 * max(R.loc[tic, m + "_s"] for m in adm)),
             "chi2r": {m: float(R.loc[tic, m + "_chi2r"]) for m in MODELOS},
             "tudo_dPdt_adm": [float(R.loc[tic, "tudo_" + m + "_dPdt"]) for m in adm],
             "tudo_s_adm": [float(R.loc[tic, "tudo_" + m + "_s"]) for m in adm]}
        linhas.append(d)
        if verboso:
            print(f"  {tic} {d['nome']:<18s} nu {nu} intervalo [{lo:.2f}, {hi:.2f}] | admitidos {', '.join(adm) if not marcado else adm[0] + ' (MENOS REJEITADO)'} "
                  f"| chi2/nu " + ", ".join(f"{m.split('_')[0][:4]} {d['chi2r'][m]:.2f}" for m in MODELOS)
                  + f" | z {d['z_min']:+.2f} a {d['z_max']:+.2f} | DMD {d['DMD_3sigma_adm']:.4f} -> {d['classe'].upper()}")
    return pd.DataFrame(linhas)


def parte_B(nivel, verboso=True):
    HP = pd.read_parquet(BASE / "deriva_gp_26_hiper.parquet"); HP = HP[HP.kernel == "matern32"].set_index("tic")
    SJ = pd.read_parquet(BASE / "jitter_sigma_j_quad.parquet").set_index("tic")
    S = pd.read_parquet(BASE / "jitter_primarios_26.parquet")
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    D11 = {int(t) for t in t26.index if t26.loc[t, "curvatura"]}
    linhas = []
    for tic in sorted(int(t) for t in t26.index):
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        P, sP, T0, E_ult = J.escada_do_registro(j)
        d = J.classificar_setores(S[(S.tic == tic) & (S.minimo == "primario")], P, sP, T0, E_ult, 0.0)
        d = d[d.ok]
        tess = pd.DataFrame({"fonte": ["TESS s%d" % s for s in d.setor], "t0": d.t0_btjd.values, "sig_min": d.sig_min.values,
                             "E": d.E.values, "setor": d.setor.values, "ano": H.ano(d.t0_btjd.values)}).sort_values("t0").reset_index(drop=True)
        A_p, tc_p = float(HP.loc[tic, "A_min"]), float(HP.loc[tic, "tau_c_d"])
        sj = float(SJ.loc[tic, "sj_quad"])
        E0 = float(np.mean(tess.E.values)); nu = len(tess) - 3
        cand = {}
        if nu >= 1:
            cand["branco_proprio"] = H.ajustar(K.com_jitter(tess, sj), 2, "matern32", 0.0, tc_p, P, E0)["chi2r"]
            cand["GP_A_proprio"] = H.ajustar(tess, 2, "matern32", A_p, tc_p, P, E0)["chi2r"]
        adm = [m for m in cand if admitido(cand[m], nu, nivel)[0]]
        marcado = False
        if not adm:
            adm = [min(cand, key=lambda m: distancia_rejeicao(cand[m], nu, nivel))] if cand else ["branco_proprio"]
            marcado = True
        escolhido = min(adm, key=lambda m: abs(np.log(cand[m])) if cand else 0.0) if cand else "branco_proprio"
        pontos = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
        if escolhido == "branco_proprio":
            p = pontos.copy()
            it = p.fonte.str.startswith("TESS").values
            p.loc[it, "sig_min"] = np.hypot(p.sig_min.values[it], sj)
            r = RT.refazer(p, P)
        else:
            r = GP.refazer_gp(pontos, P, "matern32", A_p, tc_p)
        det = r["p_curv"] < P_DET and r["p_adv"] < P_DET and r["vao_ok"]
        linhas.append({"tic": tic, "nome": str(t26.loc[tic, "nome"]), "nivel": nivel, "n_tess_ext": int(len(tess)), "nu": nu,
                       "sj_quad_min": sj, "A_proprio_min": A_p, "tau_c_proprio_d": tc_p,
                       "chi2r_branco": cand.get("branco_proprio", np.nan), "chi2r_gp": cand.get("GP_A_proprio", np.nan),
                       "admitidos": adm, "marcado_menos_rejeitado": marcado, "modelo_aplicado": escolhido,
                       "em_D11": tic in D11, "nucleo": tic in NUCLEO, "nao_informativo": bool(nu <= 2),
                       "dPdt": r["dPdt"], "sdPdt": r["sdPdt"], "p_curv": r["p_curv"], "p_adv": r["p_adv"], "p_gof": r["p_gof"],
                       "vao_ok": r["vao_ok"], "detecta": det})
        if verboso:
            L = linhas[-1]
            print(f"  {tic} {L['nome']:<18s} nu {nu:2d} | chi2/nu branco {L['chi2r_branco']:.2f} (sj {sj:.2f}), GP A proprio {L['chi2r_gp']:.2f} (A {A_p:.2f}) "
                  f"-> aplica {escolhido}{' (MENOS REJEITADO)' if marcado else ''} | desenho: dP/dt {r['dPdt']:+.4f} +- {r['sdPdt']:.4f}, p_curv {r['p_curv']:.4f}, p_adv {r['p_adv']:.4f}, p_gof {r['p_gof']:.3f}"
                  + (" | D11" if L["em_D11"] else "") + (" DETECTA" if det else ""))
    return pd.DataFrame(linhas)


if __name__ == "__main__":
    assert abs(co.VIES_CADEIA_MIN - (-1.26)) < 1e-9, "L usa o estado E"
    print("== PARTE A: classe dos 8 testaveis por consenso entre modelos admitidos (criterio posterior ao Passo K, so chi2/nu)\n")
    A95 = parte_A(0.95)
    cont = A95.classe.value_counts().to_dict()
    print(f"\n  contagem (95%): {cont}")
    print("\n== sensibilidade do limiar (90% e 99%)")
    tabs = {0.95: A95}
    for niv in (0.90, 0.99):
        tabs[niv] = parte_A(niv, verboso=False)
        dif = [(int(t), A95.set_index('tic').loc[t, 'classe'], tabs[niv].set_index('tic').loc[t, 'classe'])
               for t in A95.tic if A95.set_index('tic').loc[t, 'classe'] != tabs[niv].set_index('tic').loc[t, 'classe']]
        nadm = {int(t): int(tabs[niv].set_index('tic').loc[t, 'n_admitidos']) for t in A95.tic}
        print(f"   {int(niv * 100)}%: contagem {tabs[niv].classe.value_counts().to_dict()}; mudam de classe: {dif if dif else 'nenhuma'}; n admitidos por alvo {nadm}")
    A = pd.concat(tabs.values(), ignore_index=True)
    A.to_parquet(BASE / "ruido_admitido_d9.parquet", index=False)

    print("\n== PARTE B: modelo de ruido proprio escolhido no bloco estendido e aplicado ao DESENHO (26 alvos)\n")
    B95 = parte_B(0.95)
    d11 = B95[B95.em_D11]
    print(f"\n  D11 sobrevivem {int(d11.detecta.sum())} de {len(d11)}: {sorted(int(t) for t in d11[d11.detecta].tic)}")
    print(f"  caem: {sorted(int(t) for t in d11[~d11.detecta].tic)} | nucleo intacto: {bool(d11[d11.nucleo].detecta.all())} | "
          f"novas fora do D11: {sorted(int(t) for t in B95[~B95.em_D11 & B95.detecta].tic)}")
    print(f"  modelo aplicado: {B95.modelo_aplicado.value_counts().to_dict()}; marcados (nenhum admitido): {sorted(int(t) for t in B95[B95.marcado_menos_rejeitado].tic)}")
    tabsB = {0.95: B95}
    for niv in (0.90, 0.99):
        tabsB[niv] = parte_B(niv, verboso=False)
        b = tabsB[niv]; d = b[b.em_D11]
        dif = [(int(t), bool(B95.set_index('tic').loc[t, 'detecta']), bool(b.set_index('tic').loc[t, 'detecta'])) for t in B95.tic
               if bool(B95.set_index('tic').loc[t, 'detecta']) != bool(b.set_index('tic').loc[t, 'detecta'])]
        print(f"   {int(niv * 100)}%: D11 {int(d.detecta.sum())} de {len(d)}; modelo {b.modelo_aplicado.value_counts().to_dict()}; mudam de veredito: {dif if dif else 'nenhum'}")
    pd.concat(tabsB.values(), ignore_index=True).to_parquet(BASE / "ruido_admitido_26.parquet", index=False)
    out = {"expectativa": EXPECTATIVA, "regra": "modelo admitido: chi2/nu no intervalo de 95% de chi2_nu/nu; classe por consenso entre admitidos; "
                                                "nenhum admitido -> menos rejeitado, marcado; 424461577 fora do teste so-TESS (LTTE publicado)",
           "posterior_ao_passo_K": True, "usa_apenas_chi2r": True,
           "A_classes_95": {str(r.tic): r.classe for r in A95.itertuples()},
           "A_contagem": {str(k): int(v) for k, v in cont.items()},
           "B_D11_sobrevivem_95": int(d11.detecta.sum()), "B_quais": sorted(int(t) for t in d11[d11.detecta].tic),
           "B_nucleo_intacto": bool(d11[d11.nucleo].detecta.all()),
           "B_modelo_aplicado": {str(k): int(v) for k, v in B95.modelo_aplicado.value_counts().items()}}
    (BASE / "ruido_admitido.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(f"\n  -> {BASE / 'ruido_admitido_d9.parquet'}, ruido_admitido_26.parquet, ruido_admitido.json")
