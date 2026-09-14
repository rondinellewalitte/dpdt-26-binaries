# -*- coding: utf-8 -*-
"""O espelho de `vies_alavanca_26`: quanto do fator SuperWASP de 2,73 pode ser
ESTRUTURA astrofisica em vez de barra pequena?

Uma LTTE subamostrada e exatamente o caso: nas 4 epocas TESS (3 parametros)
a parabola a absorve; nas temporadas SuperWASP, a 12-16 anos de distancia,
nao. Se parte do 2,73 for sinal, multiplicar as barras por 1,65 enterra
sinal real em ruido - e o V527 Dra prova que pelo menos um alvo tem essa
estrutura. Com f = 0,24 e amplitude mediana de 2 min (comparavel as barras)
o efeito nao e desprezivel a priori.

Rodada: a populacao declarada de `epsilon_26` (P3 log-uniforme 0,1-20 a,
M_bin 1,5-2,5, M3 0,1-1,0, i isotropico, fase uniforme, e = 0) injetada nos
26 desenhos reais (instantes reais, barras ASSUMIDAS = verdadeiras, sem
erro de barra), em TODAS as epocas do alvo (a orbita e da estrela), com
ocorrencia f = 1 (toda estrela tem orbita) e f = 0,24 (a fracao na janela
de Tokovinin). Ajuste quadratico ponderado, estimador Sigma z^2 / Sigma
(1 - h) por conjunto, 2000 realizacoes por f. Sai o fator SuperWASP e o
fator TESS que a populacao sozinha gera, com intervalos, e a fracao do
excesso medido (2,73 - 1 = 1,73) que ela explicaria. Tambem: que f seria
preciso para a populacao sozinha dar 2,73 no SuperWASP, e que fator TESS
isso implicaria - contra o 1,03 medido.

EXPECTATIVA (escrita e commitada antes de rodar): f = 0,24 -> fator
SuperWASP 1,10-1,30 e TESS 1,03-1,10 (as orbitas com P3 < 4 a nao sao
parabolicas sobre os 4 setores e aparecem no TESS tambem); f = 1 -> SW
1,5-2,2, TESS 1,2-1,5. Se em f = 0,24 o SW ficar <= 1,15, a leitura "barra"
para o 2,73 esta limpa; se ficar >= 1,3, o 1,65 e parte barra e parte
astrofisica e a correcao conjunta passa a teto. O fator TESS e a trava:
uma populacao capaz de gerar 2,73 no SuperWASP tambem levantaria o TESS
acima do 1,03 medido.

RODADA NO ESTADO B (2026-09-13): os fatores medidos passam a ser lidos de
qual_barra_26.parquet (SW 1,35, chi2 0,89-2,29; TESS 1,02, 0,66-1,79) em
vez dos 2,73 / 1,03 / 1,81 gravados no codigo. Os desenhos sao os do
estado B (barras SuperWASP maiores). EXPECTATIVA (antes de rodar): em
f = 0,24 a populacao sozinha da SW ~1,5 (era 1,73: barras maiores,
mesmas amplitudes) e TESS ~1,58 (inalterado: barras TESS iguais);
condicionado a TESS <= 1,02, SW ~1,3 - MAIOR que o excesso medido inteiro
(1,35): a decomposicao aditiva da fator so-barra sqrt(1,35 - 0,3) ~ 1,0,
i.e. nenhuma correcao de barra; a "fracao do excesso" passa de 100% e
deixa de ser um numero util - o que se reporta e que orbitas em f = 0,24
produziriam mais excesso SuperWASP do que o observado. Grade em f: a
travessia de 5% fica onde estava (~0,5), porque o lado TESS nao mudou.

VARIANTE COM BARRAS RE-ESTIMADAS (2026-09-13): a injecao acima usa barras
FIXAS, cujo nulo e 1,00; a medida (1,35 nos 26; 0,73 nos 23 adequados) vem
do estimador do pipeline, cujo nulo e 0,77 (nulo_barras_reestimadas_26).
Comparar os dois e comparar estimadores diferentes: o "f = 0,21" e o "119%"
nao valem. Agora a injecao re-estima as barras como o pipeline (std das
temporadas simuladas + orbita, max formal, hypot 0,78; metades no TESS) e
reporta, por f, o fator SuperWASP e TESS com ESSE estimador, mais
P(SW >= medido | f) nos 26 (contra 1,35) e nos 23 adequados (contra 0,73).
EXPECTATIVA (antes de rodar): sem orbita (f = 0) o nulo reproduz 0,77;
em f = 0,24 o fator SuperWASP re-estimado fica ~0,95-1,10 (a re-estimacao
absorve parte da orbita na barra: excesso menor que o +0,42 do caso fixo)
e o TESS ~1,5 (as metades nao absorvem uma orbita de anos); nos 26,
P(SW >= 1,35 | f) sobe de ~0 em f = 0 para 0,1-0,3 em f = 0,24; nos 23,
P(SW >= 0,73 | f) cai abaixo de 5% so para f >= 0,7 (limite fraco pelo
lado SuperWASP). Desvio em qualquer direcao vai para o humano.

PROPAGACAO DO NULO CONDICIONADO (2026-09-13, segunda rodada, item 2): o bloco
"23 adequados" acima selecionava os 23 DESENHOS e nao aplicava o corte
p_gof dentro da realizacao - a mesma comparacao errada que a 5.2.3 retirou
(nulo_condicionado_26: TESS 0,96 -> 0,70 com o corte). Agora a grade dos
"23" roda nos 26 desenhos com o corte p_gof >= 0,05 aplicado alvo a alvo em
cada realizacao (barras re-estimadas, orbitas injetadas com f), e o fator e
lido nos sobreviventes - contra o medido nos 23 (SW 0,73, TESS 0,55).
EXPECTATIVA (antes de rodar): f = 0 reproduz o nulo condicionado (SW 0,77,
TESS 0,70) dentro do erro de Monte Carlo; P(SW >= 0,73 | f) fica em
0,6-0,9 em todo f (o lado arquival continua sem informacao); P(TESS <=
0,55 | f) parte de ~0,23 em f = 0 e cai com f: 0,10-0,18 em f = 0,24,
0,03-0,10 em f = 0,41, < 0,03 em f = 1 - i.e. o limite pelos 23 fica em
f <~ 0,4-0,7 a 95%, ainda acima do 0,24 de Tokovinin, e a 5.2.4 continua a
chama-lo de nao informativo. Se P(TESS <= 0,55 | f = 0,24) sair < 0,05, os
residuos TESS dos 23 passam a limitar f no valor publicado e a frase da
5.2.4 e o "Not claimed" mudam. O corte, com orbitas, tambem tira alvos cuja
orbita torna a parabola inadequada; e o que o corte real fez nos dados.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
import epsilon_26 as EPS  # noqa: E402
from qual_barra_26 import alavancas  # noqa: E402
import nulo_barras_reestimadas_26 as NULO  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
N_SIM = 2000


def fatores_medidos():
    """Fator de variancia Sigma z^2 / Sigma (1 - h) por conjunto, e o limite superior do IC de chi2 do TESS,
    lidos do que qual_barra_26 gravou - nao gravados no codigo (eram 2,73 / 1,03 / 1,81, do estado A)."""
    from scipy import stats
    d = pd.read_parquet(BASE / "qual_barra_26.parquet")
    out = {}
    for fonte in ("SuperWASP", "TESS"):
        g = d[d.fonte == fonte]
        esperado = float((1 - g.h).sum())
        razao = float((g.z ** 2).sum()) / esperado
        out[fonte] = (razao, esperado * razao / stats.chi2.ppf(0.975, esperado), esperado * razao / stats.chi2.ppf(0.025, esperado))
    return out


_F = fatores_medidos()
F_SW_MEDIDO, F_TESS_MEDIDO, F_TESS_95 = round(_F["SuperWASP"][0], 2), round(_F["TESS"][0], 2), round(_F["TESS"][2], 2)


def desenhos():
    out = []
    for p in sorted(BASE.glob("oc__*.json")):
        j = json.loads(p.read_text(encoding="utf-8"))
        if j.get("status") != "ok" or "parabola" not in j:
            continue
        pts = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
        out.append((pts.E.values.astype(float), pts.t0.values, pts.sig_min.values / 1440.0, pts.fonte.str.startswith("TESS").values,
                    pts.sig_formal_min.values.astype(float), int(j["tic"])))
    assert len(out) == 26
    return out


def realizacao(des, f, rng, reestimar=False):
    """reestimar: as barras sao re-estimadas do residuo simulado (ruido + orbita) como o pipeline faz, com
    `nulo_barras_reestimadas_26.barras_reestimadas` (em minutos) - o estimador da MEDIDA, cujo nulo e 0,77."""
    zT, hT, zS, hS = [], [], [], []
    for E, t, sig, tess, formal, tic in des:
        r = rng.normal(0, sig)
        if rng.uniform() < f:
            P3 = 10 ** rng.uniform(np.log10(EPS.P3_MIN), np.log10(EPS.P3_MAX))
            A = EPS.amplitude_min(P3, rng.uniform(*EPS.MBIN), rng.uniform(*EPS.M3), rng.uniform(0, 1)) / 1440.0
            r = r + A * np.sin(2 * np.pi * (t - t[0]) / (P3 * 365.25) + 2 * np.pi * rng.uniform())
        if reestimar:
            sig = NULO.barras_reestimadas({"sig": sig * 1440.0, "formal": formal, "tess": tess}, r * 1440.0, rng) / 1440.0
        w = 1 / sig
        Amat = np.vander(E, 3) * w[:, None]
        coef, *_ = np.linalg.lstsq(Amat, r * w, rcond=None)
        z = (r - np.polyval(coef, E)) / sig
        h = alavancas(E, sig)
        zT += list(z[tess] ** 2); hT += list(1 - h[tess]); zS += list(z[~tess] ** 2); hS += list(1 - h[~tess])
    return sum(zT) / sum(hT), sum(zS) / sum(hS)


def realizacao_condicionada(des, f, rng, p_corte=0.05):
    """Como `realizacao(..., reestimar=True)`, mas com o corte p_gof >= p_corte aplicado alvo a alvo DENTRO da
    realizacao (o mesmo de nulo_barras_reestimadas_26.realizacao_condicionada); fator nos sobreviventes."""
    from scipy import stats
    zT, hT, zS, hS, n_marc = [], [], [], [], 0
    for E, t, sig, tess, formal, tic in des:
        r = rng.normal(0, sig)
        if rng.uniform() < f:
            P3 = 10 ** rng.uniform(np.log10(EPS.P3_MIN), np.log10(EPS.P3_MAX))
            A = EPS.amplitude_min(P3, rng.uniform(*EPS.MBIN), rng.uniform(*EPS.M3), rng.uniform(0, 1)) / 1440.0
            r = r + A * np.sin(2 * np.pi * (t - t[0]) / (P3 * 365.25) + 2 * np.pi * rng.uniform())
        sig = NULO.barras_reestimadas({"sig": sig * 1440.0, "formal": formal, "tess": tess}, r * 1440.0, rng) / 1440.0
        w = 1 / sig
        Amat = np.vander(E, 3) * w[:, None]
        coef, *_ = np.linalg.lstsq(Amat, r * w, rcond=None)
        z = (r - np.polyval(coef, E)) / sig
        if stats.chi2.sf(float((z ** 2).sum()), len(E) - 3) < p_corte:
            n_marc += 1
            continue
        h = alavancas(E, sig)
        zT += list(z[tess] ** 2); hT += list(1 - h[tess]); zS += list(z[~tess] ** 2); hS += list(1 - h[~tess])
    return sum(zT) / sum(hT), sum(zS) / sum(hS), n_marc


if __name__ == "__main__":
    des = desenhos()
    rng = np.random.default_rng(2024)
    linhas = []
    print(f"populacao de {EPS.__name__}: P3 {EPS.P3_MIN}-{EPS.P3_MAX} a, M_bin {EPS.MBIN}, M3 {EPS.M3}; medido: SW {F_SW_MEDIDO}, TESS {F_TESS_MEDIDO} (limite 95% {F_TESS_95})")
    for f in (0.24, 0.41, 1.0):
        fT, fS = np.array([realizacao(des, f, rng) for _ in range(N_SIM)]).T
        # a trava do TESS: entre as realizacoes compativeis com o fator TESS medido
        # (<= 1,03; e <= 1,81, o limite a 95%), quanto de SuperWASP a populacao gera?
        c1, c2 = fT <= F_TESS_MEDIDO, fT <= F_TESS_95
        L = {"f": f, "fator_T": fT.mean(), "T_med": np.median(fT), "T_lo": np.percentile(fT, 2.5), "T_hi": np.percentile(fT, 97.5),
             "fator_S": fS.mean(), "S_med": np.median(fS), "S_lo": np.percentile(fS, 2.5), "S_hi": np.percentile(fS, 97.5),
             "fracao_do_excesso_SW": (fS.mean() - 1) / (F_SW_MEDIDO - 1),
             "P_T_le_medido": c1.mean(), "S_dado_T_le_medido": fS[c1].mean() if c1.any() else np.nan,
             "S_dado_T_le_medido_p95": np.percentile(fS[c1], 95) if c1.any() else np.nan,
             "P_T_le_181": c2.mean(), "S_dado_T_le_181": fS[c2].mean() if c2.any() else np.nan,
             "S_dado_T_le_181_p95": np.percentile(fS[c2], 95) if c2.any() else np.nan,
             "P_S_ge_medido": (fS >= F_SW_MEDIDO).mean()}
        linhas.append(L)
        print(f"  f = {f:.2f}: fator TESS media {L['fator_T']:.2f}, mediana {L['T_med']:.2f} ({L['T_lo']:.2f}–{L['T_hi']:.2f}) | SuperWASP media {L['fator_S']:.2f}, mediana {L['S_med']:.2f} ({L['S_lo']:.2f}–{L['S_hi']:.2f}) "
              f"| media explica {L['fracao_do_excesso_SW']:.0%} do excesso SW medido | P(SW >= {F_SW_MEDIDO}) = {L['P_S_ge_medido']:.2f}", flush=True)
        print(f"           trava do TESS: P(TESS <= {F_TESS_MEDIDO}) = {L['P_T_le_medido']:.2f}; nessas realizacoes SW medio {L['S_dado_T_le_medido']:.2f} (p95 {L['S_dado_T_le_medido_p95']:.2f}) "
              f"| P(TESS <= {F_TESS_95}) = {L['P_T_le_181']:.2f}; nessas SW medio {L['S_dado_T_le_181']:.2f} (p95 {L['S_dado_T_le_181_p95']:.2f})", flush=True)
    d = pd.DataFrame(linhas)
    d.to_parquet(BASE / "estrutura_swasp_26.parquet", index=False)
    # o fator do TESS como verossimilhanca sobre f (revisao 6, item 2): a
    # populacao declarada preve excesso no TESS que cresce com f; o 1,03
    # medido e uma restricao INTERNA sobre f, independente do Tokovinin.
    # Grade em f, P(fator TESS <= 1,03 | f) e P(<= 1,81 | f), 1000 realizacoes.
    print("\n== o fator TESS medido como restricao sobre f (populacao declarada)")
    grade = []
    for f in (0.05, 0.10, 0.15, 0.20, 0.24, 0.30, 0.41, 0.50, 0.70, 1.00):
        fT, fS = np.array([realizacao(des, f, rng) for _ in range(1000)]).T
        grade.append({"f": f, "P_T_le_103": (fT <= F_TESS_MEDIDO).mean(), "P_T_le_181": (fT <= F_TESS_95).mean(),
                      "T_med": np.median(fT), "S_med": np.median(fS)})
        print(f"  f = {f:.2f}: P(TESS <= {F_TESS_MEDIDO}) = {grade[-1]['P_T_le_103']:.3f}  P(TESS <= {F_TESS_95}) = {grade[-1]['P_T_le_181']:.3f}  (medianas TESS {grade[-1]['T_med']:.2f}, SW {grade[-1]['S_med']:.2f})", flush=True)
    g = pd.DataFrame(grade)
    g.to_parquet(BASE / "estrutura_swasp_26_grade_f.parquet", index=False)
    for nivel, col in ((0.05, "P_T_le_103"), (0.05, "P_T_le_181")):
        acima = g[g[col] < nivel]
        print(f"  limite: P({col}) cai abaixo de {nivel} a partir de f = {acima.f.min():.2f}" if len(acima) else f"  {col} nunca cai abaixo de {nivel} na grade")
    # que f daria 2,73 so por orbitas? (interpolacao linear no excesso, que e ~linear em f)
    exc1 = d[d.f == 1.0].fator_S.iloc[0] - 1
    f_nec = (F_SW_MEDIDO - 1) / exc1 if exc1 > 0 else np.inf
    print(f"\n  f necessario para a populacao sozinha dar {F_SW_MEDIDO} no SuperWASP: {f_nec:.2f}"
          + (f" (impossivel: > 1; em f = 1 o SW so chega a {d[d.f == 1.0].fator_S.iloc[0]:.2f})" if f_nec > 1 else
             f"; nesse f o TESS ficaria em ~{1 + f_nec * (d[d.f == 1.0].fator_T.iloc[0] - 1):.2f} contra {F_TESS_MEDIDO} medido"))
    print(f"  -> {BASE / 'estrutura_swasp_26.parquet'}")

    # ---- o mesmo, com o ESTIMADOR DA MEDIDA (barras re-estimadas; nulo 0,77) ----
    # A comparacao que vale: a injecao e a medida com o mesmo estimador. Nos 26 (medido 1,35) e nos 23
    # com parabola adequada, o criterio ja declarado em 3.5 (medido 0,73) - ver nulo_barras_reestimadas_26.
    print("\n== injecao com barras RE-ESTIMADAS (o estimador do pipeline; nulo esperado 0,77 em f = 0)")
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    adequados = set(t26.index[~t26.inadequada])
    q = pd.read_parquet(BASE / "qual_barra_26.parquet")
    med = {}
    for rot, sub in (("26", None), ("23 adequados", adequados)):
        qq = q if sub is None else q[q.tic.isin(sub)]
        med[rot] = {f: float((g.z ** 2).sum() / (1 - g.h).sum()) for f, g in qq.groupby("fonte")}
    print(f"  medido: 26 -> SW {med['26']['SuperWASP']:.2f}, TESS {med['26']['TESS']:.2f} | 23 adequados -> SW {med['23 adequados']['SuperWASP']:.2f}, TESS {med['23 adequados']['TESS']:.2f}")
    linhas = []
    print("  (o bloco '23 adequados' roda nos 26 desenhos com o corte p_gof >= 0,05 DENTRO de cada realizacao; expectativa: f = 0 -> SW 0,77, TESS 0,70)")
    for rot, sub in (("26", None), ("23 adequados", adequados)):
        for f in (0.0, 0.10, 0.24, 0.41, 0.70, 1.0):
            if sub is None:
                fT, fS = np.array([realizacao(des, f, rng, reestimar=True) for _ in range(1000)]).T
            else:
                R = np.array([realizacao_condicionada(des, f, rng) for _ in range(1000)])
                fT, fS, nmarc = R[:, 0], R[:, 1], R[:, 2]
                print(f"  {rot:12s} f = {f:.2f}: marcados por realizacao {nmarc.mean():.2f}")
            mS, mT = med[rot]["SuperWASP"], med[rot]["TESS"]
            L = {"conjunto": rot, "f": f, "fator_S": fS.mean(), "S_med": float(np.median(fS)), "S_lo": np.percentile(fS, 2.5), "S_hi": np.percentile(fS, 97.5),
                 "fator_T": fT.mean(), "T_med": float(np.median(fT)), "T_lo": np.percentile(fT, 2.5), "T_hi": np.percentile(fT, 97.5),
                 "medido_S": mS, "medido_T": mT, "P_S_ge_medido": float((fS >= mS).mean()), "P_T_le_medido": float((fT <= mT).mean())}
            linhas.append(L)
            print(f"  {rot:12s} f = {f:.2f}: SW media {L['fator_S']:.2f}, mediana {L['S_med']:.2f} ({L['S_lo']:.2f}-{L['S_hi']:.2f}); TESS media {L['fator_T']:.2f}, mediana {L['T_med']:.2f} "
                  f"| P(SW >= {mS:.2f}) = {L['P_S_ge_medido']:.3f}  P(TESS <= {mT:.2f}) = {L['P_T_le_medido']:.3f}", flush=True)
    dr = pd.DataFrame(linhas)
    dr.to_parquet(BASE / "estrutura_swasp_26_reestimada.parquet", index=False)
    for rot in ("26", "23 adequados"):
        g = dr[dr.conjunto == rot]
        baixo = g[g.P_S_ge_medido < 0.05]
        print(f"  {rot}: P(SW >= medido | f) < 0,05 " + (f"a partir de f = {baixo.f.min():.2f}" if len(baixo) else "em nenhum f da grade")
              + f"; P(TESS <= medido | f) < 0,05 " + (f"a partir de f = {g[g.P_T_le_medido < 0.05].f.min():.2f}" if (g.P_T_le_medido < 0.05).any() else "em nenhum f da grade"))
    print(f"  -> {BASE / 'estrutura_swasp_26_reestimada.parquet'}")
