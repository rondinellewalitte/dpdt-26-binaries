# -*- coding: utf-8 -*-
"""O teste do revisor (item 1): que Sigma z^2 / Sigma (1 - h) o NULO produz
quando a barra de cada epoca e RE-ESTIMADA dos dados simulados, exatamente
como o pipeline faz - e nao mantida fixa, como em `vies_alavanca_26`.

O pipeline (estado B): a barra de cada temporada SuperWASP e
hypot(max(formal, std dos desvios das n temporadas, ddof=1), 0,78), com
n = 2 (17 alvos) ou 3 (9 alvos); a barra de cada setor TESS e
max(formal, |t0_a - t0_b| / sqrt 2), das duas metades do setor. Em ambos
os casos s e estimado com nu = 1-2 graus de liberdade a partir dos MESMOS
pontos que entram no residuo. Com s ~ sigma sqrt(chi2_nu / nu) independente
de r, z^2 = r^2 / s^2 segue F(1, nu), cuja media e nu / (nu - 2): indefinida
para nu <= 2. O piso (formal, 0,78) trunca a cauda; o quanto sobra e o que
este script mede.

Simulacao nos 26 desenhos reais (estado B): para cada alvo, as epocas nas
posicoes reais; ruido gaussiano com sigma VERDADEIRO = a barra atribuida do
estado B (a melhor estimativa que temos), em todas as epocas; SEM curvatura
e SEM orbitas. As barras SuperWASP sao re-estimadas: std (ddof=1) dos
desvios simulados das n temporadas do alvo em torno da media delas,
max com a formal real, hypot 0,78 - a mesma regra do pipeline. As barras
TESS sao re-estimadas como max(formal real, |a - b| / sqrt 2) com a e b
duas realizacoes de metade de setor (sigma_metade = sigma sqrt 2). Ajuste
quadratico ponderado pelas barras re-estimadas, z e h com elas, e o
estimador Sigma z^2 / Sigma (1 - h) por conjunto. 2000 realizacoes.
Sai: a distribuicao nula do fator de cada conjunto (media, mediana,
2,5-97,5%), e onde o fator MEDIDO (SW 1,35, TESS 1,02) cai nela.

EXPECTATIVA (escrita e commitada antes de rodar):
  - SuperWASP: o nulo NAO e centrado em 1. Media 1,3-1,8 (a cauda de F(1,1)
    truncada pelo piso), mediana 1,1-1,4, 97,5% acima de 2,5. O 1,35 medido
    cai dentro do intervalo central (entre os percentis 25 e 75): o excesso
    SuperWASP que sobrou depois do conserto do sqrt n E o vies de estimar a
    dispersao com 2-3 amostras, e nao pede ruido vermelho nem orbita.
  - TESS: nulo perto de 1 (media 1,0-1,2), porque a formal domina a barra da
    maioria dos setores (o piso corta a cauda); o 1,02 medido cai no meio.
  - Se o nulo SuperWASP der media < 1,15, o 1,35 medido volta a ser excesso
    (fraco) e a leitura muda; se der > 2, o estimador e inutilizavel com
    nu = 1-2 e a nota tem de dizer isso e trocar de teste.
Desvio em qualquer direcao vai para o humano antes de qualquer outra coisa.

RESULTADO (rodado depois): SW 0,77 (0,63-0,93), TESS 0,96 - a expectativa acima
errou de direcao (r e s vem dos mesmos pontos; n = 2 -> z^2 = 0,5 por construcao).

SEGUNDA RODADA, item 2 (revisor: "a particao dos 23 e condicionada e o nulo dela
nao e"): a medida nos 23 adequados (SW 0,73, TESS 0,55) foi feita nos alvos que
PASSARAM o corte p_gof >= 0,05, e o nulo com que se comparou (0,78 nos 23
desenhos, em estrutura_swasp_26) nao aplicou corte nenhum. Aqui o corte entra
DENTRO de cada realizacao: p_gof do quadratico de cada alvo com as barras
re-estimadas (chi2 = Sigma z^2 do alvo, N - 3 gl, exatamente como tabela_nota),
descarta-se quem tem p_gof < 0,05, e o fator e calculado nos sobreviventes.

EXPECTATIVA (escrita e commitada antes de rodar):
  - Quantos o corte marca no nulo: como as barras re-estimadas deflacionam o
    chi2 (fator 0,77), o corte nominal de 5% marca MENOS de 5% dos alvos:
    0,3-1,0 alvo por realizacao (26 x ~1-4%), contra 3 nos dados. P(>= 3
    marcados | nulo) < 0,10. Se sair >= 0,2, os tres marcados nao sao excesso
    e a 5.2.3 tem de dizer isso.
  - Fator SuperWASP nos sobreviventes: 0,72-0,76 (o corte tira as realizacoes
    de chi2 alto, mas tira poucos alvos). Medido 0,73 nos 23: percentil 20-45%.
  - Fator TESS nos sobreviventes: 0,85-0,93. Medido 0,55: percentil 5-15%. Se
    cair abaixo de 5%, "neither archive shows excess" na 5.2.3 vira "SuperWASP
    consistent with its null; TESS below it", e o abstract e a Secao 6 seguem.
  - Sanidade: a mesma realizacao sem corte reproduz a distribuicao acima
    (0,77 / 0,96) dentro do erro de Monte Carlo.
RESULTADO (rodado depois): marcados 0,83 por realizacao, P(>= 3) = 0,057; SW
condicionado 0,77 (0,63-0,92), medido 0,73 no percentil 30%; TESS condicionado
0,70 (0,37-1,08), medido 0,55 no percentil 23%. DESVIO: o nulo TESS cai muito
mais do que o esperado (0,70, nao 0,85-0,93) - quem e marcado no nulo sao as
realizacoes em que uma barra TESS co-estimada de duas metades saiu pequena; o
corte tira a cauda pesada do TESS (diagnostico abaixo: fracao do z^2 TESS nos
marcados). A leitura ("nos 23 nenhum dos dois conjuntos mostra excesso") nao
muda; o percentil do TESS sobe de 13% (nulo errado) para 23% (nulo certo).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cadeia_oc as co  # noqa: E402
import config  # noqa: E402
from qual_barra_26 import alavancas  # noqa: E402
from scipy import stats  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
N_SIM = 2000
SEMENTE = 7


def desenhos():
    out = []
    for p in sorted(BASE.glob("oc__*.json")):
        j = json.loads(p.read_text(encoding="utf-8"))
        if j.get("status") != "ok" or "parabola" not in j:
            continue
        pts = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
        tess = pts.fonte.str.startswith("TESS").values
        formal = pts.sig_formal_min.values.astype(float)     # gravada por oc_lote para TESS e SuperWASP (estado B)
        assert np.isfinite(formal).all(), (j["tic"], "ponto sem barra formal: a cadeia precisa ter rodado com o campo")
        out.append({"tic": int(j["tic"]), "E": pts.E.values.astype(float), "sig": pts.sig_min.values.astype(float),
                    "formal": formal, "tess": tess})
    assert len(out) == 26, len(out)
    return out


def barras_reestimadas(d, r, rng):
    """r: residuos simulados (min) com sigma verdadeiro = d['sig']. Devolve as barras que o pipeline atribuiria."""
    sig = d["sig"].copy()
    tess, formal = d["tess"], d["formal"]
    # SuperWASP: std dos desvios das n temporadas em torno da media, max formal, hypot 0,78
    isw = np.where(~tess)[0]
    if len(isw) >= 2:
        disp = float(np.std(r[isw], ddof=1))
        sig[isw] = np.hypot(np.maximum(formal[isw], disp), co.VIES_CADEIA_SIG)
    # TESS: |a - b| / sqrt 2 de duas metades (cada metade com sigma sqrt 2 x a do setor), max formal
    itess = np.where(tess)[0]
    a = rng.normal(0, d["sig"][itess] * np.sqrt(2))
    b = rng.normal(0, d["sig"][itess] * np.sqrt(2))
    sig[itess] = np.maximum(formal[itess], np.abs(a - b) / np.sqrt(2))
    return sig


def realizacao(des, rng, reestimar=True, vies_comum=False):
    """vies_comum: o 0,78 min NAO e espalhamento por temporada, e um deslocamento COMUM do bloco SuperWASP do
    alvo (item 2 da revisao). A verdade passa a ser sigma_interno por epoca (sqrt(sig^2 - 0,78^2)) mais um
    delta ~ N(0, 0,78) partilhado pelas temporadas; a re-estimacao (std das temporadas, hypot 0,78) e a mesma."""
    zT, hT, zS, hS = [], [], [], []
    for d in des:
        if vies_comum:
            sig_v = np.where(d["tess"], d["sig"], np.sqrt(np.maximum(d["sig"] ** 2 - co.VIES_CADEIA_SIG ** 2, 1e-6)))
            r = rng.normal(0, sig_v) + np.where(d["tess"], 0.0, rng.normal(0, co.VIES_CADEIA_SIG))
        else:
            r = rng.normal(0, d["sig"])
        sig = barras_reestimadas(d, r, rng) if reestimar else d["sig"]
        w = 1 / sig
        X = np.vstack([np.ones_like(d["E"]), d["E"], d["E"] ** 2]).T
        beta, *_ = np.linalg.lstsq(X * w[:, None], r * w, rcond=None)
        res = r - X @ beta
        z = res / sig
        h = alavancas(d["E"], sig)
        zT.extend(z[d["tess"]]); hT.extend(h[d["tess"]]); zS.extend(z[~d["tess"]]); hS.extend(h[~d["tess"]])
    zT, hT, zS, hS = map(np.array, (zT, hT, zS, hS))
    return (zT ** 2).sum() / (1 - hT).sum(), (zS ** 2).sum() / (1 - hS).sum()


def realizacao_condicionada(des, rng, p_corte=0.05):
    """Uma realizacao do nulo com barras re-estimadas; devolve o fator (TESS, SW) SEM corte e COM o corte
    p_gof >= p_corte aplicado alvo a alvo dentro da realizacao, mais o numero de alvos marcados."""
    por_alvo = []
    for d in des:
        r = rng.normal(0, d["sig"])
        sig = barras_reestimadas(d, r, rng)
        w = 1 / sig
        X = np.vstack([np.ones_like(d["E"]), d["E"], d["E"] ** 2]).T
        beta, *_ = np.linalg.lstsq(X * w[:, None], r * w, rcond=None)
        z = (r - X @ beta) / sig
        h = alavancas(d["E"], sig)
        chi2, dof = float((z ** 2).sum()), len(d["E"]) - 3
        p_gof = float(stats.chi2.sf(chi2, dof))              # a mesma regra de tabela_nota (inadequada = p_gof < 0,05)
        por_alvo.append((z, h, d["tess"], p_gof >= p_corte))

    def fator(sel):
        zT = np.concatenate([z[t] for z, h, t, ok in por_alvo if ok or not sel])
        hT = np.concatenate([h[t] for z, h, t, ok in por_alvo if ok or not sel])
        zS = np.concatenate([z[~t] for z, h, t, ok in por_alvo if ok or not sel])
        hS = np.concatenate([h[~t] for z, h, t, ok in por_alvo if ok or not sel])
        return (zT ** 2).sum() / (1 - hT).sum(), (zS ** 2).sum() / (1 - hS).sum()

    n_marc = sum(1 for *_, ok in por_alvo if not ok)
    zT2_marc = sum(float((z[t] ** 2).sum()) for z, h, t, ok in por_alvo if not ok)
    zT2_tot = sum(float((z[t] ** 2).sum()) for z, h, t, ok in por_alvo)
    return fator(False), fator(True), n_marc, (zT2_marc, zT2_tot)


if __name__ == "__main__":
    des = desenhos()
    q = pd.read_parquet(BASE / "qual_barra_26.parquet")
    medido = {}
    for fonte in ("TESS", "SuperWASP"):
        g = q[q.fonte == fonte]
        medido[fonte] = float((g.z ** 2).sum() / (1 - g.h).sum())
    print(f"medido (estado B): SuperWASP {medido['SuperWASP']:.2f}, TESS {medido['TESS']:.2f}")
    rng = np.random.default_rng(SEMENTE)
    linhas = []
    for modo, reest, comum in (("barras fixas (como vies_alavanca)", False, False), ("barras RE-ESTIMADAS (como o pipeline)", True, False),
                               ("re-estimadas, 0,78 como vies COMUM do bloco", True, True)):
        fT, fS = np.array([realizacao(des, rng, reest, comum) for _ in range(N_SIM)]).T
        for fonte, f, m in (("SuperWASP", fS, medido["SuperWASP"]), ("TESS", fT, medido["TESS"])):
            pct = float((f <= m).mean())
            L = {"modo": modo, "fonte": fonte, "media": f.mean(), "mediana": float(np.median(f)), "p2.5": np.percentile(f, 2.5), "p25": np.percentile(f, 25),
                 "p75": np.percentile(f, 75), "p97.5": np.percentile(f, 97.5), "medido": m, "percentil_do_medido": pct}
            linhas.append(L)
            print(f"  {modo:38s} {fonte:9s}: nulo media {L['media']:.2f}, mediana {L['mediana']:.2f}, 2,5-97,5% {L['p2.5']:.2f}-{L['p97.5']:.2f}, "
                  f"25-75% {L['p25']:.2f}-{L['p75']:.2f} | medido {m:.2f} esta no percentil {pct:.0%}", flush=True)
    d = pd.DataFrame(linhas)
    d.to_parquet(BASE / "nulo_barras_reestimadas_26.parquet", index=False)
    print(f"  -> {BASE / 'nulo_barras_reestimadas_26.parquet'}")

    # ---- segunda rodada, item 2: o corte p_gof DENTRO do nulo ----
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    adequados = set(t26.index[~t26.inadequada])
    n_marc_dados = int(t26.inadequada.sum())
    qa = q[q.tic.isin(adequados)]
    med23 = {fonte: float((g.z ** 2).sum() / (1 - g.h).sum()) for fonte, g in qa.groupby("fonte")}
    print("\n== nulo CONDICIONADO ao corte p_gof >= 0,05 dentro de cada realizacao (expectativa: 0,3-1,0 marcados por realizacao; "
          "SW 0,72-0,76; TESS 0,85-0,93)")
    print(f"  dados: {n_marc_dados} marcados de 26; medido nos {len(adequados)} adequados: SW {med23['SuperWASP']:.2f}, TESS {med23['TESS']:.2f}")
    rng = np.random.default_rng(SEMENTE + 1)
    R = [realizacao_condicionada(des, rng) for _ in range(N_SIM)]
    sem = np.array([r[0] for r in R]); com = np.array([r[1] for r in R]); nm = np.array([r[2] for r in R])
    zT2 = np.array([r[3] for r in R]); frac_zT2_marc = float(zT2[:, 0].sum() / zT2[:, 1].sum())
    print(f"  marcados por realizacao: media {nm.mean():.2f}, mediana {np.median(nm):.0f}, max {nm.max()}; "
          f"P(>= {n_marc_dados} marcados | nulo) = {(nm >= n_marc_dados).mean():.3f}; taxa por alvo {nm.mean() / 26:.3%}")
    print(f"  diagnostico: fracao do Sigma z^2 TESS total que esta nos alvos-realizacao marcados: {frac_zT2_marc:.1%} "
          f"(em {nm.mean() / 26:.1%} dos alvos-realizacao) - e o que o corte tira do nulo TESS")
    linhas2 = [{"n_marcados_dados": n_marc_dados, "marcados_media": nm.mean(), "marcados_mediana": float(np.median(nm)),
                "P_ge_marcados_dados": float((nm >= n_marc_dados).mean()), "frac_zT2_TESS_nos_marcados": frac_zT2_marc}]
    for k, rot in ((0, "TESS"), (1, "SuperWASP")):
        f0, f1, m = sem[:, k], com[:, k], med23[rot]
        pct0, pct1 = float((f0 <= medido[rot]).mean()), float((f1 <= m).mean())
        print(f"  {rot:9s}: sem corte (sanidade, medido 26 = {medido[rot]:.2f}) media {f0.mean():.2f}, 2,5-97,5% {np.percentile(f0, 2.5):.2f}-{np.percentile(f0, 97.5):.2f}, percentil {pct0:.0%}"
              f" | COM corte (medido 23 = {m:.2f}) media {f1.mean():.2f}, mediana {np.median(f1):.2f}, 2,5-97,5% {np.percentile(f1, 2.5):.2f}-{np.percentile(f1, 97.5):.2f}, "
              f"25-75% {np.percentile(f1, 25):.2f}-{np.percentile(f1, 75):.2f}, percentil do medido {pct1:.0%}", flush=True)
        linhas2[0].update({f"{rot}_sem_media": f0.mean(), f"{rot}_sem_pct": pct0, f"{rot}_com_media": f1.mean(), f"{rot}_com_mediana": float(np.median(f1)),
                           f"{rot}_com_p2.5": np.percentile(f1, 2.5), f"{rot}_com_p97.5": np.percentile(f1, 97.5), f"{rot}_com_p25": np.percentile(f1, 25),
                           f"{rot}_com_p75": np.percentile(f1, 75), f"{rot}_medido_23": m, f"{rot}_com_pct": pct1})
    pd.DataFrame(linhas2).to_parquet(BASE / "nulo_condicionado_26.parquet", index=False)
    print(f"  -> {BASE / 'nulo_condicionado_26.parquet'}")
