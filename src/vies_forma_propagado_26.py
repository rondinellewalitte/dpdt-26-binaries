# -*- coding: utf-8 -*-
"""Rodada treze, bloco C: o vies de FORMA POR ALVO (Secao 3.4) entra no orcamento de erro,
das duas maneiras que o parecer pede, e o teste de janela e recalculado sob cada uma.

O PROBLEMA. O teste de forma (`vies_forma_swasp.py`) mede, alvo a alvo, um vies do estimador
SuperWASP de ate +-3 min, com o template empirico e a cobertura de fase REAIS daquele alvo.
Hoje nada disso entra no orcamento: as epocas arquivais carregam so o vies COMUM da cadeia,
-1,26 +- 0,75 min. Um offset comum as estacoes arquivais de um alvo e a forma maximamente
degenerada com curvatura - e por isso que `gls_vies_26.py` existe -, e a variancia de modo
comum escala com o quadrado do termo: 0,75 -> 3 min e fator 16. Nos alvos com barras arquivais
de 1,5 a 5 min o termo passa a ser comparavel ou maior que a barra por epoca.
`sensib_vies_26.py` NAO responde a isto: ele varre o valor central de um vies comum a todos os
alvos, na diagonal, e nao infla sigma nenhum.

DEPOIS DO ESTADO F (rodada treze, C1) a correcao passou a ser da CADEIA: cada epoca arquival ja
sai corrigida pelo vies de forma do proprio alvo, e o piso da barra arquival ja carrega a
incerteza dessa medida. Este script deixa de ser quem corrige e passa a ser quem MEDE A
SENSIBILIDADE em volta, com quatro tratamentos por alvo:

  cadeia        - o estado F como esta, com o piso do alvo movido da diagonal para bloco comum;
  sem_correcao  - desfaz a correcao e volta ao estado E (e a sanidade de que a troca e reversivel);
  conservadora  - bloco comum de 3 min, o maior vies de forma que o teste encontra, como limite;
  forma_alt     - a outra forma simulada do alvo, porque a regra escolhe uma das duas.

A REGRA DA FORMA, declarada antes da primeira corrida: vale a forma cuja DISPERSAO ENTRE
TEMPORADAS SIMULADA mais se aproxima da real (menor |log(disp_sint/disp_real)|); alvo marcado
`instavel` entra sem correcao e com piso de 3 min. O que NAO e propagado, e fica dito: o vies de
forma tem uma parte comum (a media sobre as temporadas, que e o que se corrige) e uma parte que
varia de temporada para temporada; so a comum entra. A variante conservadora e o limite que cobre
as duas. Corrigir pela forma em cima do vies da cadeia nao e contagem dupla porque o mesmo teste,
nos quatro calibradores, da -0,10 min (Secao 3.4).

A EXPECTATIVA da primeira corrida (estado E) esta no commit 17864a2 e o resultado em 0f65c7d: o
portao nao disparou - 4 contraditados em todas as variantes, nenhuma classe mudou -, o que mudou
foi a contagem de DETECCOES sob o limite conservador. A varredura do bloco comum (0,75 a 10 min)
mede, alvo a alvo, qual o menor bloco que derrota a deteccao: e a coluna nova da Tabela 3 e a
tabela do Apendice E.

SANIDADES, todas antes de qualquer numero: a diagonal pura reproduz `tabela_26.parquet`; o bloco
com o piso do proprio alvo reproduz `gls_vies_26.parquet`; e desfazer a correcao reproduz o dP/dt
do estado E registrado em `previsao_estado_F.json`.

Uso:
    python src/vies_forma_propagado_26.py
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
import gls_vies_26 as G  # noqa: E402
import ruido_admitido_26 as RA  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="backslashreplace")
except (AttributeError, ValueError):
    pass

BASE = config.DATA / "orquestra" / "oc_lote"
S_POR_ANO = 86400.0 * 365.25
P_DET = 0.05
BLOCO_CONSERVADOR = 3.0          # min, o maior vies de forma que a Secao 3.4 relata
VARIANTES = ("cadeia", "sem_correcao", "conservadora", "forma_alt")
Z_IC95 = 1.959963984540054       # ic95_media_min -> sigma da media


def covariancia(sig_min, tess, sigma_bloco_min, piso_alvo_min=None):
    """C em dias^2. sigma_bloco_min=None -> diagonal pura (a cadeia). Caso contrario, o PISO da
    barra arquival daquele alvo (0,75 min ate o estado E; hypot(0,75; sigma do vies de forma) no
    estado F) sai da diagonal das epocas SuperWASP e o bloco comum entra no lugar."""
    sig = np.asarray(sig_min, float) / 1440.0
    C = np.diag(sig ** 2)
    if sigma_bloco_min is None:
        return C
    piso = float(co.VIES_CADEIA_SIG if piso_alvo_min is None else piso_alvo_min) / 1440.0
    b = float(sigma_bloco_min) / 1440.0
    isw = np.where(~tess)[0]
    for i in isw:
        C[i, i] = max(sig[i] ** 2 - piso ** 2, 1e-12)
    C[np.ix_(isw, isw)] += b ** 2
    return C


def ajuste(pontos, P, sigma_bloco_min, delta_min=0.0, piso_alvo_min=None):
    """Ajuste do DESENHO (as 4-7 epocas das Tabelas 3-5) com deslocamento das epocas arquivais e
    bloco comum. delta_min > 0 ATRASA as epocas arquivais (subtrai), como a correcao de forma."""
    E = pontos.E.values.astype(float)
    t = pontos.t0.values.astype(float).copy()
    tess = pontos.fonte.str.startswith("TESS").values
    if delta_min:
        t[~tess] -= float(delta_min) / 1440.0          # vies = ajustado - verdadeiro: corrigir e subtrair
    C = covariancia(pontos.sig_min.values, tess, sigma_bloco_min, piso_alvo_min)
    _, _, chi2l, _ = G.ajustar_gls(E, t, C, 1)
    coef, _, chi2p, cov = G.ajustar_gls(E, t, C, 2)
    dof = len(E) - 3
    return {"dPdt": 2 * coef[0] / P * S_POR_ANO, "s": 2 * np.sqrt(cov[0, 0]) / P * S_POR_ANO,
            "chi2r": chi2p / dof, "p_gof": float(stats.chi2.sf(chi2p, dof)),
            "p_curv": float(stats.chi2.sf(chi2l - chi2p, 1)), "p_adv": G.adversarial_gls(E, t, C),
            "n": len(E), "n_sw": int((~tess).sum())}


def forma_por_alvo(resumo):
    """Regra declarada: vale a forma cuja dispersao entre temporadas simulada mais se aproxima da
    real; devolve tambem a outra, como sensibilidade."""
    escolha = {}
    for tic, g in resumo[resumo.tipo == "alvo"].groupby("curva"):
        g = g.copy()
        g["dist"] = np.abs(np.log(g.razao_disp_sint_real.replace(0, np.nan)))
        g = g.sort_values("dist")
        princ, outra = g.iloc[0], g.iloc[1]
        escolha[int(tic)] = {
            "forma": princ.forma, "vies_min": float(princ.vies_medio_min),
            "sig_vies_min": float(princ.ic95_media_min) / Z_IC95,
            "instavel": bool(princ.instavel), "razao_disp": float(princ.razao_disp_sint_real),
            "forma_alt": outra.forma, "vies_alt_min": float(outra.vies_medio_min),
            "sig_vies_alt_min": float(outra.ic95_media_min) / Z_IC95,
            "instavel_alt": bool(outra.instavel), "razao_disp_alt": float(outra.razao_disp_sint_real)}
    return escolha


def variantes(tic, f):
    """As colunas de tratamento, por alvo, DEPOIS do estado F: a cadeia ja aplica a correcao de
    forma do proprio alvo, entao o que resta e (i) a propria cadeia, (ii) desfazer a correcao,
    para ver o estado E ao lado, (iii) o limite conservador de 3 min e (iv) a forma alternativa.
    Cada tupla e (nome, deslocamento adicional das epocas arquivais em min, bloco comum em min)."""
    aplicado, piso = co.vies_forma(tic)                # o que a cadeia ja pos nos registros
    return [("cadeia", 0.0, piso),
            ("sem_correcao", -aplicado, co.VIES_CADEIA_SIG),
            ("conservadora", 0.0, BLOCO_CONSERVADOR),
            # a regra da forma escolhe uma das duas simuladas; em alguns alvos a outra tem vies
            # varias vezes maior, e e por isso que ela e ajustada e nao so tabelada
            ("forma_alt", (0.0 if f["instavel_alt"] else f["vies_alt_min"]) - aplicado,
             BLOCO_CONSERVADOR if f["instavel_alt"] else float(np.hypot(co.VIES_CADEIA_SIG, f["sig_vies_alt_min"])))], piso


def previsao_estado_E():
    """Os dP/dt do estado E, da previsao commitada em comparar_estados_EF --prever, para a
    sanidade de que desfazer a correcao de forma volta exatamente ao estado anterior."""
    p = config.DATA / "orquestra" / "previsao_estado_F.json"
    if not p.exists():
        return None
    return pd.DataFrame(json.loads(p.read_text(encoding="utf-8"))["alvos"]).set_index("tic")


PREV_E = None


def main():
    global PREV_E
    PREV_E = previsao_estado_E()
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    gls = pd.read_parquet(BASE / "gls_vies_26.parquet").set_index("tic")
    resumo = pd.read_parquet(BASE / "vies_forma_swasp_resumo.parquet")
    escolha = forma_por_alvo(resumo)
    D11 = sorted(t26.index[t26.curvatura])

    linhas = []
    print("== C1: o desenho dos 26 sob cada variante (sanidade primeiro)")
    for tic in sorted(t26.index):
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        pontos = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
        P = float(j["P_escada_d"])
        f = escolha[int(tic)]

        # sanidade 1: diagonal pura reproduz a cadeia
        d = ajuste(pontos, P, None)
        # 1e-6 e a tolerancia que gls_vies_26 ja usa para a mesma comparacao: a cadeia ajusta por
        # outro caminho numerico e a diferenca observada e de 1e-12 a 3e-8 s/ano (arredondamento)
        assert abs(d["dPdt"] - t26.loc[tic, "dPdt"]) < 1e-6 and abs(d["s"] - t26.loc[tic, "s"]) < 1e-6, (tic, "diagonal != tabela_26", d)
        assert abs(d["p_curv"] - t26.loc[tic, "p"]) < 1e-6 and abs(d["p_adv"] - t26.loc[tic, "padv"]) < 1e-6, (tic, "p != tabela_26")
        vs, piso = variantes(int(tic), f)
        # sanidade 2: o bloco com o piso do proprio alvo reproduz a GLS da Secao 3.4 (estado F)
        h = ajuste(pontos, P, piso, 0.0, piso)
        assert abs(h["dPdt"] - gls.loc[tic, "dPdt_gls"]) < 1e-9 and abs(h["s"] - gls.loc[tic, "s_gls"]) < 1e-9, (tic, "bloco do piso != gls_vies_26")
        # sanidade 3: desfazer a correcao de forma reproduz o dP/dt do estado E (previsao
        # commitada). Desfazer e as DUAS coisas: a epoca volta e o piso da barra volta a 0,75 -
        # so o deslocamento deixaria os pesos do estado F e o ajuste saia 1e-6 fora. A tolerancia e
        # 1e-8 porque a barra volta por duas raizes (subtrai piso^2, soma 0,75^2): o residuo
        # medido no pior alvo e 1,6e-9 s/ano, arredondamento puro
        if PREV_E is not None and int(tic) in PREV_E.index:
            pe = pontos.copy()
            tess_e = pe.fonte.str.startswith("TESS").values
            delta_ap, piso_ap = co.vies_forma(int(tic))
            pe.loc[~tess_e, "t0"] = pe.t0.values[~tess_e] + delta_ap / 1440.0
            interna_e = np.sqrt(np.maximum(pe.sig_min.values[~tess_e] ** 2 - piso_ap ** 2, 0.0))
            pe.loc[~tess_e, "sig_min"] = np.hypot(interna_e, co.VIES_CADEIA_SIG)
            e = ajuste(pe, P, None)
            alvo_E = float(PREV_E.loc[int(tic), "dPdt_E"])
            # tolerancia RELATIVA: nos tres alvos instaveis o piso do estado F e 3,09 min e
            # reconstruir a barra antiga subtrai dois quadrados quase iguais (9,9 - 9,6), o que
            # custa duas casas de precisao; o pior residuo medido e 2,7e-8 em 0,219 s/ano
            assert abs(e["dPdt"] - alvo_E) <= 1e-6 * max(abs(alvo_E), 1e-3), (tic, "sem_correcao != estado E", e["dPdt"], alvo_E)

        for nome, delta, bloco in vs:
            r = ajuste(pontos, P, bloco, delta, piso)
            linhas.append({"tic": int(tic), "nome": t26.loc[tic, "nome"], "em_D11": tic in D11,
                           "variante": nome, "delta_min": delta, "bloco_min": bloco,
                           "forma": f["forma"], "instavel": f["instavel"], "razao_disp": f["razao_disp"],
                           "vies_forma_min": f["vies_min"], "sig_vies_forma_min": f["sig_vies_min"],
                           "forma_alt": f["forma_alt"], "vies_forma_alt_min": f["vies_alt_min"],
                           "dPdt": r["dPdt"], "s": r["s"], "p_curv": r["p_curv"], "p_adv": r["p_adv"],
                           "p_gof": r["p_gof"], "chi2r": r["chi2r"], "n": r["n"], "n_sw": r["n_sw"],
                           "detecta": bool(r["p_curv"] < P_DET and r["p_adv"] < P_DET),
                           "dPdt_diag": float(t26.loc[tic, "dPdt"]), "s_diag": float(t26.loc[tic, "s"])})
    R = pd.DataFrame(linhas)
    R.to_parquet(BASE / "vies_forma_propagado_26.parquet", index=False)
    print("  sanidade: diagonal == tabela_26, bloco do piso == gls_vies_26, sem_correcao == estado E")

    base = R[R.variante == "cadeia"].set_index("tic")
    print(f"\n== C2: os {len(D11)} do D11, por variante (sigma antes -> depois, veredito)")
    for tic in D11:
        b = base.loc[tic]
        print(f"  {tic} {b['nome']:<18s} forma {b['forma']:<10s}{' INSTAVEL' if b['instavel'] else '':9s} "
              f"vies {b['vies_forma_min']:+.2f} +- {b['sig_vies_forma_min']:.2f} min (alt {b['forma_alt'][:4]} {b['vies_forma_alt_min']:+.2f})")
        for v in VARIANTES:
            r = R[(R.tic == tic) & (R.variante == v)].iloc[0]
            print(f"      {v:<19s} bloco {r['bloco_min']:.2f} | dP/dt {r['dPdt']:+.4f} +- {r['s']:.4f} "
                  f"(s x{r['s'] / b['s']:.2f}, desloc {abs(r['dPdt'] - b['dPdt']) / b['s']:.2f} sigma) | "
                  f"p_curv {r['p_curv']:.3g} p_adv {r['p_adv']:.3g} -> {'DETECTA' if r['detecta'] else 'cai'}")

    print("\n== C3: contagem do D11 por variante")
    for v in VARIANTES:
        sub = R[(R.variante == v) & R.em_D11]
        caem = sorted(sub[~sub.detecta].tic.tolist())
        novas = sorted(R[(R.variante == v) & ~R.em_D11 & R.detecta].tic.tolist())
        razao = sub.s.values / base.loc[sub.tic.values, "s"].values
        print(f"  {v:<19s} sobrevivem {int(sub.detecta.sum())}/{len(sub)} | caem {caem} | novas {novas} | "
              f"s x: mediana {float(np.median(razao)):.2f}, max {float(razao.max()):.2f} ({int(sub.tic.values[razao.argmax()])})")

    # ---------------- teste de janela ----------------
    print("\n== C4: teste de janela - z recalculado com o sigma do desenho de cada variante")
    rec = pd.read_parquet(BASE / "reconcilia_hi_d9.parquet").set_index("tic")
    saida, contagens = [], {}
    for v in VARIANTES:
        sub = R[R.variante == v].set_index("tic")
        print(f"\n  -- {v} --")
        classes = {}
        for tic in rec.index:
            nu = int(rec.loc[tic, "n_tess"]) - 3
            adm = [m for m in RA.MODELOS if RA.admitido(float(rec.loc[tic, m + "_chi2r"]), nu, 0.95)[0]]
            marcado = not adm
            if marcado:
                adm = [min(RA.MODELOS, key=lambda m: RA.distancia_rejeicao(float(rec.loc[tic, m + "_chi2r"]), nu, 0.95))]
            d7, s7 = float(sub.loc[tic, "dPdt"]), float(sub.loc[tic, "s"])
            zs = [(float(rec.loc[tic, m + "_dPdt"]) - d7) / float(np.hypot(rec.loc[tic, m + "_s"], s7)) for m in adm]
            classe = "LTTE publicado" if tic == RA.LTTE_PUBLICADO else RA.classe_por_consenso(zs)
            classes[int(tic)] = classe
            saida.append({"variante": v, "tic": int(tic), "nome": rec.loc[tic, "nome"], "admitidos": adm,
                          "marcado_menos_rejeitado": marcado, "dPdt_desenho": d7, "s_desenho": s7,
                          "s_desenho_publicado": float(rec.loc[tic, "s_7"]), "dPdt_desenho_publicado": float(rec.loc[tic, "dPdt_7"]),
                          "z_min": float(min(zs, key=abs)), "z_max": float(max(zs, key=abs)), "classe": classe,
                          "acima_de_3": bool(all(abs(z) > 3 for z in zs)), "acima_de_4": bool(all(abs(z) > 4 for z in zs))})
            r = saida[-1]
            print(f"    {tic} {r['nome']:<18s} s_desenho {s7:.4f} (publicado {r['s_desenho_publicado']:.4f}, x{s7 / r['s_desenho_publicado']:.2f}) "
                  f"| |z| {abs(r['z_min']):.2f} a {abs(r['z_max']):.2f} -> {classe.upper()}")
        contagens[v] = {"contraditado": sum(c == "contraditado" for c in classes.values()),
                        "ambiguo": sum(c == "ambiguo" for c in classes.values()),
                        "nao contraditado": sum(c == "nao contraditado" for c in classes.values()),
                        "classes": classes}
    J = pd.DataFrame(saida)
    J.to_parquet(BASE / "vies_forma_propagado_janela.parquet", index=False)

    print("\n== C5: PORTAO - contagem de contraditados por variante")
    ref = contagens["cadeia"]["contraditado"]
    mudou = False
    for v, c in contagens.items():
        acima3 = sum(1 for r in saida if r["variante"] == v and r["tic"] != RA.LTTE_PUBLICADO and r["acima_de_3"])
        acima4 = sum(1 for r in saida if r["variante"] == v and r["tic"] != RA.LTTE_PUBLICADO and r["acima_de_4"])
        dif = [t for t in c["classes"] if c["classes"][t] != contagens["cadeia"]["classes"][t]]
        mudou |= bool(dif) and v != "cadeia"
        print(f"  {v:<19s} contraditados {c['contraditado']} | ambiguos {c['ambiguo']} | nao contraditados {c['nao contraditado']} "
              f"| |z|>3 em todos os admitidos: {acima3} | |z|>4: {acima4} | classes que mudam: {dif if dif else 'nenhuma'}")
    print(f"\n  referencia de hoje: {ref} contraditados")
    if mudou:
        print("  PORTAO: a contagem/classe MUDA em alguma variante -> parar e reportar antes de tocar em texto.")
    else:
        print("  PORTAO: nenhuma classe muda em nenhuma variante.")

    print("")
    print("== C6: varredura do bloco comum - onde cada coisa quebra (os 26; z nos 8 do teste)")
    grade = [0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0]
    quebra = []
    for tic in sorted(t26.index):
        j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
        pontos = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
        P = float(j["P_escada_d"])
        no_teste = tic in rec.index
        if no_teste:
            nu = int(rec.loc[tic, "n_tess"]) - 3
            adm = [m for m in RA.MODELOS if RA.admitido(float(rec.loc[tic, m + "_chi2r"]), nu, 0.95)[0]]
            if not adm:
                adm = [min(RA.MODELOS, key=lambda m: RA.distancia_rejeicao(float(rec.loc[tic, m + "_chi2r"]), nu, 0.95))]
        zs_por_bloco, det_por_bloco = {}, {}
        for bloco in grade:
            r = ajuste(pontos, P, bloco, 0.0, co.vies_forma(int(tic))[1])
            if no_teste:
                zs = [(float(rec.loc[tic, m + "_dPdt"]) - r["dPdt"]) / float(np.hypot(rec.loc[tic, m + "_s"], r["s"])) for m in adm]
                zs_por_bloco[bloco] = min(abs(z) for z in zs)
            det_por_bloco[bloco] = bool(r["p_curv"] < P_DET and r["p_adv"] < P_DET)
        z3 = next((b for b in grade if no_teste and zs_por_bloco[b] <= 3.0), None)
        dcai = next((b for b in grade if not det_por_bloco[b]), None)
        detecta = bool(t26.loc[tic, "curvatura"])
        quebra.append({"tic": int(tic), "nome": t26.loc[tic, "nome"], "em_D11": detecta, "no_teste_de_janela": no_teste,
                       "classe_hoje": contagens["cadeia"]["classes"].get(int(tic), "fora do teste"),
                       "bloco_min_para_z_ate_3": z3, "bloco_min_para_perder_deteccao": dcai if detecta else None,
                       "z_por_bloco": {str(b): zs_por_bloco.get(b) for b in grade},
                       "detecta_por_bloco": {str(b): det_por_bloco[b] for b in grade}})
        if no_teste or detecta:
            print(f"    {tic} {quebra[-1]['nome']:<18s} {quebra[-1]['classe_hoje']:<16s} " +
                  ("|z| " + " ".join(f"{b:g}:{zs_por_bloco[b]:.2f}{'' if det_por_bloco[b] else '*'}" for b in grade)
                   if no_teste else "fora do teste de janela") +
                  (f" | cai abaixo de 3 em {z3 if z3 else '> 10'} min" if no_teste else "") +
                  (f"; perde a deteccao em {dcai if dcai else '> 10'} min" if detecta else "; nao e deteccao"))
    V = pd.DataFrame(quebra)
    V.to_parquet(BASE / "vies_forma_propagado_varredura.parquet", index=False)
    d11v = V[V.em_D11]
    print(f"  deteccoes por bloco: " + ", ".join(
        f"{b:g} min: {int(sum(r['detecta_por_bloco'][str(b)] for _, r in d11v.iterrows()))}/{len(d11v)}" for b in grade))

    resumo_json = {"varredura": quebra, "contagens": {v: {k: c[k] for k in ("contraditado", "ambiguo", "nao contraditado")} for v, c in contagens.items()},
                   "classes": {v: {str(t): cl for t, cl in c["classes"].items()} for v, c in contagens.items()},
                   "bloco_conservador_min": BLOCO_CONSERVADOR,
                   "forma_por_alvo": {str(k): v for k, v in escolha.items()},
                   "D11_sobrevive": {v: sorted(R[(R.variante == v) & R.em_D11 & R.detecta].tic.tolist())
                                     for v in R.variante.unique()}}
    (BASE / "vies_forma_propagado_26.json").write_text(json.dumps(resumo_json, indent=1, default=float), encoding="utf-8")
    print(f"\n  -> {BASE / 'vies_forma_propagado_26.parquet'}\n  -> {BASE / 'vies_forma_propagado_janela.parquet'}\n  -> {BASE / 'vies_forma_propagado_26.json'}")


if __name__ == "__main__":
    main()
