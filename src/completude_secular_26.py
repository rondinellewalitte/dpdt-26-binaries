# -*- coding: utf-8 -*-
"""Segunda rodada de revisao, itens 3 + 4: completude para um termo SECULAR e a
auto-supressao do estimador de barras.

O que se injeta: dP/dt secular puro (sem LTTE), t(E) = Q E^2 + ruido, com
Q = (dP/dt / S_POR_ANO) x P / 2 (d/ciclo^2; a mesma convencao da Tabela 3:
232634196 com -0,214 s/ano e P = 1,2511 d da Q = -4,244e-9, conferido), nas
26 geometrias reais (ciclos E dos registros do estado B', contados do primeiro
setor TESS), ruido gaussiano com sigma VERDADEIRO = a barra atribuida do
estado B, sinal com sinal aleatorio por realizacao. Grade de |dP/dt|:
0 (falso alarme), 0,005, 0,01, 0,02, 0,05, 0,10, 0,20, 0,25 s/ano.

O ponto inteiro: as barras sao RE-ESTIMADAS dentro do laco, das temporadas e
metades simuladas, com a regra do pipeline (`nulo_barras_reestimadas_26.
barras_reestimadas`: SuperWASP = hypot(max(formal, std das temporadas
simuladas - sinal incluido), 0,78); TESS = max(formal, |a - b|/sqrt 2 das
metades)). No pipeline real as temporadas sao medidas contra o periodo da
escada TESS, e um termo secular aparece como espalhamento entre temporadas:
e exatamente isso que a re-estimacao absorve na barra. Com barras fixas
mede-se outra coisa; as duas sao gravadas para a diferenca ser vista.

Duas curvas por alvo e por amplitude (N_SIM realizacoes por celula):
  - completude: fracao que passa curvatura (p_curv < 0,05) E o gate
    adversarial (p_adv < 0,05; a regra de cadeia_oc.adversarial), com as
    barras re-estimadas; e a mesma com as barras verdadeiras fixas;
  - supressao: sigma(dP/dt) da covariancia do quadratico com as barras
    re-estimadas sobre sigma(dP/dt) com as barras verdadeiras (mediana e
    quartis por celula); e o vies do coeficiente (mediana de dP/dt
    recuperado / injetado).
Tambem: fracao de realizacoes com p_gof < 0,05 (parabola "inadequada").

Limite de leitura: as barras atribuidas do estado B sao a melhor estimativa
do ruido, mas nos alvos onde ha estrutura real elas ja a contem; o ruido
simulado e entao MAIOR que o real e a supressao medida e um PISO.

EXPECTATIVA DO REVISOR/AUTOR (escrita para ser conferida): supressao
desprezivel abaixo de ~0,02 s/ano, crescendo para 20-40% em 0,2 s/ano;
completude limitada pelas barras e nao pelo padrao de amostragem, como ε.

EXPECTATIVA MINHA (escrita e commitada antes de rodar; onde difere da acima,
a medida decide):
  - falso alarme em |dP/dt| = 0 com barras re-estimadas: 1-4% por alvo (o
    gate adversarial segura abaixo dos 5% nominais);
  - completude 50% atingida entre 0,01 e 0,03 s/ano na maioria dos 24 com
    ε > 0, e so acima de 0,1 nos dois com ε = 0 (224605072, 229476285); em
    0,2 s/ano completude > 0,9 em todos menos esses dois; a amplitude de
    50% correlaciona com a barra SuperWASP mediana do alvo (Spearman > 0,6)
    e nao com o numero de epocas;
  - supressao: 0,9-1,0 em |dP/dt| <= 0,005 (co-estimacao com nu pequeno
    subestima s); >= 1,2 ja em 0,02 s/ano nos alvos com temporadas a >= 1
    ano de distancia - porque um termo de 0,02 s/ano a E ~ -8000 ciclos
    (P = 0,5 d, 12 anos antes do TESS) move duas temporadas a 730 ciclos
    uma da outra por Q x 2 x 8000 x 730 ~ 2,7 min, o tamanho das barras;
    mediana sobre os 26 em 0,2 s/ano entre 1,5 e 3, e 2-5 nos alvos cuja
    alavanca e o arquivo. Isto e MAIOR que os 20-40% da expectativa acima;
    se a mediana em 0,2 sair < 1,4, a conta do espalhamento entre temporadas
    esta errada em algum lugar e e preciso olhar antes de escrever;
  - vies do coeficiente: nenhum (estimador linear com pesos simetricos):
    mediana dP/dt recuperado / injetado = 1,00 +- 0,05 em toda a grade;
  - completude com barras re-estimadas contra fixas: menor por < 10 pontos
    em toda a grade (a barra cresce, mas o sinal que a faz crescer e grande);
  - p_gof < 0,05 com termo secular real: raro (< 5%) em toda a grade - a
    barra re-estimada absorve o proprio sinal e a parabola nunca fica
    "inadequada" por causa de um secular. E o que torna o corte cego a ele.
Desvio em qualquer direcao vai para o humano antes de qualquer outra coisa.

RESULTADO (rodado depois; 300 realizacoes por celula):
  - falso alarme em 0: mediana 0,2%, <= 6% em 25 de 26, 16% em 229476285
    (barras SW de 359 min: o ajuste fica com 3 epocas TESS);
  - amplitude de 50%: mediana 0,015 s/ano; >= 0,5 em 16 de 26 em 0,02; em
    0,2 > 0,9 em 21, < 0,75 nos 5 de barra SW > 8 min, e 3 nunca chegam a
    0,5 (224605072, 229476285, 359552377) - mais que os 2 esperados;
    Spearman(A50, barra SW) +0,42 (p 0,05), (A50, n epocas) -0,32 (p 0,13);
  - supressao mediana 1,02 / 1,04 / 1,12 / 1,15 / 1,26 / 1,61 / 2,82 / 3,39
    na grade; em 0,2 s/ano 8-30x em 9 alvos de barra SW 1,5-5 min (max 24x
    em 198408416). A expectativa do revisor/autor (20-40% em 0,2) falhou por
    uma ordem de grandeza nos alvos de barra pequena; a minha acertou a
    mediana (2,8 em 1,5-3) e errou o teto (2-5, saiu 8-30);
  - vies do coeficiente 1,00 em toda a grade; p_gof < 0,05 mediana 2-4%,
    max 21% num desenho de 4 epocas (243352373);
  - completude re-estimada contra fixa: -25 pontos em 0,02 (0,65 vs 0,90),
    mais que os < 10 esperados; iguais de 0,05 para cima;
  - 232634196 no proprio coeficiente: completude 0,94, supressao 1,7 (com o
    ruido = as barras de 19,8 min).
TERCEIRA CONSEQUENCIA (revisor, terceira rodada, item 1): se um termo real
infla a barra arquival co-estimada, nos alvos com termo real o fator de
residuos SuperWASP (Sigma z^2 / Sigma (1 - h), o estimador da 5.2.3) cai
ABAIXO do nulo de 0,77 - a barra ficou grande para o ruido que sobra depois
que a parabola absorve o termo. Logo "barras corretas" e "barras infladas
por termos reais" dao a mesma leitura do teste. Modo --fator: o termo
injetado nos 26 (sinal aleatorio) em cada amplitude da grade, barras
re-estimadas, fator SuperWASP e TESS agregados por realizacao, sem corte
(contra 1,35 nos 26) e com o corte p_gof dentro (contra 0,73 nos 23);
1000 realizacoes por amplitude.
EXPECTATIVA (escrita e commitada antes de rodar): SW cai monotonicamente
com a amplitude - ~0,75 em 0,005, ~0,70 em 0,02, ~0,5 em 0,05, ~0,3 em 0,2
(a barra e sinal; o residuo e ruido pequeno); o TESS NAO se move (a parabola
absorve o termo exatamente e as barras de metades nao o veem). Com o corte,
o mesmo. P(SW <= 0,73 | A) sobe de ~0,30 em 0 para > 0,5 em 0,02: o 0,73
medido nos 23 e compativel com nenhum termo E com termos de ~0,02 s/ano
em todos - o teste nao separa. Se o SW NAO cair (ficar em 0,77 na grade
inteira), a leitura do revisor esta errada e e preciso entender por que
antes de escrever.
RESULTADO do --fator (rodado depois): com o corte, SW 0,77 / 0,74 / 0,59 /
0,43 / 0,30 em 0,005 / 0,01 / 0,02 / 0,05 / 0,2 (caiu mais rapido que o
esperado); TESS 0,71-0,72 (parado, esperado). DESVIO: SEM o corte a MEDIA
do fator primeiro sobe (0,81 / 0,85 / 0,79) e so depois cai (0,65 / 0,49 /
0,39), com cauda a 1,5-2,3. A leitura escrita na hora ("o termo infla as
barras arquivais e entrega o ajuste ao bloco TESS, cujo erro de
extrapolacao aparece nos residuos arquivais") era conjectura. Tres
discriminadores, mesmo modo, mesma grade:
  (a) barras FIXAS (verdadeiras) com o termo injetado: o termo esta no
      espaco do modelo, logo residuo = ruido projetado e o fator tem de
      ficar em 1,00 em toda a grade. Se ficar, a subida e da
      re-estimacao; se subir tambem, a leitura esta errada.
  (b) MEDIANA ao lado da media, com e sem corte: se a mediana sem corte
      NAO sobe (fica <= 0,77 e cai), a "subida" e da cauda pesada (poucas
      realizacoes com fator alto) e o texto tem de dizer "a media";
      se a mediana sobe, a subida e do grosso da distribuicao.
  (c) alavanca media das epocas SuperWASP por amplitude (barras
      re-estimadas): se cai de ~0,43 em 0 para < 0,3 em 0,02-0,05, o
      ajuste esta mesmo passando ao bloco TESS; se nao cai, a leitura
      "entrega o ajuste ao TESS" nao se sustenta.
EXPECTATIVA (antes de rodar): (a) 1,00 +- 0,02 na grade inteira; (b) a
mediana sem corte tambem sobe, para ~0,80 em 0,01, e cai depois - i.e. o
efeito e do grosso, nao so da cauda; (c) a alavanca SW cai de 0,43 para
~0,25 em 0,05 e ~0,1 em 0,2. Se (b) der mediana monotona, a frase da 5.2.3
vira "a media" e a cauda; se (c) nao cair, a leitura sai do texto.
Uso: python src/completude_secular_26.py            (simula e grava)
     python src/completude_secular_26.py --resumo   (so le o parquet e deriva por alvo)
     python src/completude_secular_26.py --fator    (fator de residuos conforme o termo injetado)
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
import nulo_barras_reestimadas_26 as NULO  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
S_POR_ANO = 86400.0 * 365.25
GRADE = (0.0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.25)     # |dP/dt| em s/ano
N_SIM = 300
P_DET = 0.05
SEMENTE = 11


def desenhos():
    out = []
    for p in sorted(BASE.glob("oc__*.json")):
        j = json.loads(p.read_text(encoding="utf-8"))
        if j.get("status") != "ok" or "parabola" not in j:
            continue
        pts = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
        d = {"tic": int(j["tic"]), "P": float(j["P_escada_d"]), "E": pts.E.values.astype(float), "sig": pts.sig_min.values.astype(float),
             "formal": pts.sig_formal_min.values.astype(float), "tess": pts.fonte.str.startswith("TESS").values}
        assert np.isfinite(d["formal"]).all(), d["tic"]
        out.append(d)
    assert len(out) == 26, len(out)
    return out


def Q_de_dpdt(dpdt_s_ano, P):
    return 0.5 * (dpdt_s_ano / S_POR_ANO) * P


def testes(E, t, sig_d, P):
    """p_curv, p_adv, p_gof, dP/dt e sigma(dP/dt) (s/ano) com as barras sig_d (dias) - a regra da cadeia."""
    _, _, chi2l, _ = co.ajustar(E, t, sig_d, 1)
    coef, _, chi2p, cov = co.ajustar(E, t, sig_d, 2)
    p_curv = float(stats.chi2.sf(chi2l - chi2p, 1))
    p_gof = float(stats.chi2.sf(chi2p, len(E) - 3))
    dpdt = 2 * coef[0] / P * S_POR_ANO
    sdpdt = 2 * np.sqrt(cov[0, 0]) / P * S_POR_ANO
    if p_curv < P_DET:                                   # o adversarial so decide onde a curvatura passou (como epsilon_26)
        ag = pd.DataFrame({"E": E, "t0": t, "sig_min": sig_d * 1440.0})
        p_adv = co.adversarial(ag, P, float(t[0]))["p_adversarial"]
    else:
        p_adv = 1.0
    return p_curv, p_adv, p_gof, dpdt, sdpdt


def celula(d, A, rng, n_sim=N_SIM):
    E, P = d["E"], d["P"]
    det_re, det_fix, gof_re, sup, rec = [], [], [], [], []
    for _ in range(n_sim):
        sinal = rng.choice((-1.0, 1.0)) if A > 0 else 1.0
        Q = Q_de_dpdt(sinal * A, P)
        r_min = Q * E ** 2 * 1440.0 + rng.normal(0, d["sig"])          # sinal + ruido, em minutos
        sig_re = NULO.barras_reestimadas(d, r_min, rng)                 # a regra do pipeline, sinal incluido
        t = r_min / 1440.0
        pc, pa, pg, dp, sd = testes(E, t, sig_re / 1440.0, P)
        pc0, pa0, _, _, sd0 = testes(E, t, d["sig"] / 1440.0, P)       # barras verdadeiras, fixas
        det_re.append(pc < P_DET and pa < P_DET); det_fix.append(pc0 < P_DET and pa0 < P_DET); gof_re.append(pg < P_DET)
        sup.append(sd / sd0)
        if A > 0:
            rec.append(dp / (sinal * A))
    sup = np.array(sup)
    return {"tic": d["tic"], "dpdt": A, "completude_re": float(np.mean(det_re)), "completude_fixa": float(np.mean(det_fix)),
            "frac_pgof": float(np.mean(gof_re)), "supressao_med": float(np.median(sup)), "supressao_p25": float(np.percentile(sup, 25)),
            "supressao_p75": float(np.percentile(sup, 75)), "vies_med": float(np.median(rec)) if rec else np.nan, "n_sim": n_sim}


def resumo_por_alvo(R, des, t26):
    """Quantidades derivadas por alvo, do parquet: amplitude de 50% (interpolacao linear na grade), falso alarme,
    completude e supressao no |dP/dt| tabulado do proprio alvo (interpolacao), completude em 0,02 e 0,2."""
    barra_sw = {d["tic"]: float(np.median(d["sig"][~d["tess"]])) for d in des}
    n_ep = {d["tic"]: len(d["E"]) for d in des}
    linhas = []
    for tic, g in R.groupby("tic"):
        g = g.sort_values("dpdt")
        c, a, sup = g.completude_re.values, g.dpdt.values, g.supressao_med.values
        a_prop = float(abs(t26.loc[tic, "dPdt"]))
        linhas.append({"tic": int(tic), "em_D10": bool(t26.loc[tic, "curvatura"]), "n_epocas": n_ep[tic], "barra_sw_med_min": barra_sw[tic],
                       "falso_alarme": float(c[a == 0][0]), "A50": float(np.interp(0.5, c, a)) if c.max() >= 0.5 else np.nan,
                       "dpdt_tab": a_prop, "completude_no_proprio": float(np.interp(a_prop, a, c)), "supressao_no_proprio": float(np.interp(a_prop, a, sup)),
                       "completude_0.02": float(c[a == 0.02][0]), "completude_0.2": float(c[a == 0.2][0]), "supressao_0.2": float(sup[a == 0.2][0])})
    T = pd.DataFrame(linhas)
    ok = T.A50.notna()
    rb = stats.spearmanr(T.A50[ok], T.barra_sw_med_min[ok]); rn = stats.spearmanr(T.A50[ok], T.n_epocas[ok])
    print(f"\n== por alvo (derivado do parquet)")
    print(f"  falso alarme: mediana {T.falso_alarme.median():.3f}, <= 0,06 em {int((T.falso_alarme <= 0.06).sum())} de 26, max {T.falso_alarme.max():.2f} em {int(T.loc[T.falso_alarme.idxmax(), 'tic'])}")
    print(f"  falso alarme somado sobre os 26 (deteccoes esperadas de puro ruido): {T.falso_alarme.sum():.2f}; maximo entre os D10: {T[T.em_D10].falso_alarme.max():.3f}")
    print(f"  A50: mediana {T.A50.median():.3f} s/ano; nao atingem 0,5: {sorted(T.tic[~ok].tolist())}; Spearman(A50, barra SW) {rb.statistic:+.2f} (p {rb.pvalue:.2f}), (A50, n) {rn.statistic:+.2f} (p {rn.pvalue:.2f})")
    print(f"  completude >= 0,5 em 0,02: {int((T['completude_0.02'] >= 0.5).sum())} de 26; > 0,9 em 0,2: {int((T['completude_0.2'] > 0.9).sum())}; < 0,75 em 0,2: {sorted((int(t), round(c, 2)) for t, c in zip(T.tic, T['completude_0.2']) if c < 0.75)}")
    print(f"  supressao em 0,2 >= 5: {sorted(((int(t), round(v, 1)) for t, v in zip(T.tic, T['supressao_0.2']) if v >= 5), key=lambda x: -x[1])}")
    d10 = T[T.em_D10].sort_values("tic")
    print(f"  D10 no proprio |dP/dt|: completude {d10.completude_no_proprio.min():.2f}-{d10.completude_no_proprio.max():.2f} (mediana {d10.completude_no_proprio.median():.2f}); "
          f"supressao {d10.supressao_no_proprio.min():.2f}-{d10.supressao_no_proprio.max():.2f} (mediana {d10.supressao_no_proprio.median():.2f})")
    for r in d10.itertuples():
        print(f"    {r.tic} |dP/dt| {r.dpdt_tab:.3f}: completude {r.completude_no_proprio:.2f}, supressao {r.supressao_no_proprio:.2f}")
    T.to_parquet(BASE / "completude_secular_26_alvo.parquet", index=False)
    print(f"  -> {BASE / 'completude_secular_26_alvo.parquet'}")
    return T


def fator_por_amplitude(des, A, rng, n_sim=1000, p_corte=0.05, reestimar=True):
    """Fator de residuos (TESS, SW) por realizacao com o termo A injetado nos 26; barras re-estimadas (o pipeline) ou
    FIXAS nas verdadeiras (discriminador a); sem corte e com o corte p_gof >= p_corte dentro da realizacao (a mesma
    regra de nulo_barras_reestimadas_26). Devolve tambem a alavanca media das epocas SuperWASP (discriminador c)."""
    from qual_barra_26 import alavancas
    sem, com, hsw = [], [], []
    for _ in range(n_sim):
        por_alvo = []
        for d in des:
            E, P = d["E"], d["P"]
            sinal = rng.choice((-1.0, 1.0)) if A > 0 else 1.0
            r = Q_de_dpdt(sinal * A, P) * E ** 2 * 1440.0 + rng.normal(0, d["sig"])
            sig = NULO.barras_reestimadas(d, r, rng) if reestimar else d["sig"].copy()
            w = 1 / sig
            X = np.vstack([np.ones_like(E), E, E ** 2]).T
            beta, *_ = np.linalg.lstsq(X * w[:, None], r * w, rcond=None)
            z = (r - X @ beta) / sig
            h = alavancas(E, sig)
            ok = stats.chi2.sf(float((z ** 2).sum()), len(E) - 3) >= p_corte
            por_alvo.append((z, h, d["tess"], ok))

        def fator(sel):
            zT = np.concatenate([z[t] for z, h, t, ok in por_alvo if ok or not sel]); hT = np.concatenate([h[t] for z, h, t, ok in por_alvo if ok or not sel])
            zS = np.concatenate([z[~t] for z, h, t, ok in por_alvo if ok or not sel]); hS = np.concatenate([h[~t] for z, h, t, ok in por_alvo if ok or not sel])
            return (zT ** 2).sum() / (1 - hT).sum(), (zS ** 2).sum() / (1 - hS).sum()
        sem.append(fator(False)); com.append(fator(True))
        hsw.append(float(np.mean(np.concatenate([h[~t] for z, h, t, ok in por_alvo]))))
    return np.array(sem), np.array(com), np.array(hsw)


if __name__ == "__main__":
    des = desenhos()
    if "--fator" in sys.argv:
        q = pd.read_parquet(BASE / "qual_barra_26.parquet")
        t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
        adequados = set(t26.index[~t26.inadequada])
        m26 = {f: float((g.z ** 2).sum() / (1 - g.h).sum()) for f, g in q.groupby("fonte")}
        m23 = {f: float((g.z ** 2).sum() / (1 - g.h).sum()) for f, g in q[q.tic.isin(adequados)].groupby("fonte")}
        print(f"== fator de residuos conforme o termo secular injetado nos 26 (barras re-estimadas); medido: 26 SW {m26['SuperWASP']:.2f} TESS {m26['TESS']:.2f} | 23 SW {m23['SuperWASP']:.2f} TESS {m23['TESS']:.2f}")
        print("   expectativa: SW ~0,75 em 0,005, ~0,70 em 0,02, ~0,5 em 0,05, ~0,3 em 0,2; TESS parado; P(SW <= 0,73 | A) > 0,5 a partir de 0,02")
        rng = np.random.default_rng(SEMENTE + 5)
        linhas = []
        for A in GRADE:
            sem, com, hsw = fator_por_amplitude(des, A, rng)
            semF, comF, hswF = fator_por_amplitude(des, A, rng, reestimar=False)          # discriminador (a): barras fixas
            L = {"dpdt": A, "SW_sem": sem[:, 1].mean(), "SW_sem_med": float(np.median(sem[:, 1])), "SW_sem_p2.5": np.percentile(sem[:, 1], 2.5), "SW_sem_p97.5": np.percentile(sem[:, 1], 97.5),
                 "TESS_sem": sem[:, 0].mean(), "SW_com": com[:, 1].mean(), "SW_com_med": float(np.median(com[:, 1])), "SW_com_p2.5": np.percentile(com[:, 1], 2.5), "SW_com_p97.5": np.percentile(com[:, 1], 97.5),
                 "TESS_com": com[:, 0].mean(), "P_SW_com_le_med23": float((com[:, 1] <= m23["SuperWASP"]).mean()), "P_SW_sem_ge_med26": float((sem[:, 1] >= m26["SuperWASP"]).mean()),
                 "P_TESS_com_le_med23": float((com[:, 0] <= m23["TESS"]).mean()), "h_SW_medio": float(hsw.mean()),
                 "SW_sem_fixas": semF[:, 1].mean(), "SW_sem_fixas_med": float(np.median(semF[:, 1])), "TESS_sem_fixas": semF[:, 0].mean(), "h_SW_medio_fixas": float(hswF.mean())}
            linhas.append(L)
            print(f"  |dP/dt| = {A:.3f}: sem corte SW media {L['SW_sem']:.2f} MEDIANA {L['SW_sem_med']:.2f} ({L['SW_sem_p2.5']:.2f}-{L['SW_sem_p97.5']:.2f}), TESS {L['TESS_sem']:.2f} | com corte SW media {L['SW_com']:.2f} mediana {L['SW_com_med']:.2f} ({L['SW_com_p2.5']:.2f}-{L['SW_com_p97.5']:.2f}), TESS {L['TESS_com']:.2f} "
                  f"| P(SW_com <= {m23['SuperWASP']:.2f}) = {L['P_SW_com_le_med23']:.2f}, P(TESS_com <= {m23['TESS']:.2f}) = {L['P_TESS_com_le_med23']:.2f}, P(SW_sem >= {m26['SuperWASP']:.2f}) = {L['P_SW_sem_ge_med26']:.3f} "
                  f"| h_SW {L['h_SW_medio']:.2f} | barras FIXAS: SW {L['SW_sem_fixas']:.2f} (mediana {L['SW_sem_fixas_med']:.2f}), TESS {L['TESS_sem_fixas']:.2f}, h_SW {L['h_SW_medio_fixas']:.2f}", flush=True)
        pd.DataFrame(linhas).to_parquet(BASE / "completude_secular_26_fator.parquet", index=False)
        print(f"  -> {BASE / 'completude_secular_26_fator.parquet'}")
        sys.exit(0)
    if "--resumo" in sys.argv:
        R = pd.read_parquet(BASE / "completude_secular_26.parquet")
        t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
        resumo_por_alvo(R, des, t26)
        sys.exit(0)
    # sanidade da convencao de Q contra a Tabela 3
    j = json.loads((BASE / "oc__232634196.json").read_text(encoding="utf-8"))
    assert abs(Q_de_dpdt(j["parabola"]["dPdt_s_por_ano"], j["P_escada_d"]) - j["parabola"]["Q_d_por_ciclo2"]) < 1e-14
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    eps = pd.read_parquet(BASE / "epsilon_26.parquet").set_index("tic") if (BASE / "epsilon_26.parquet").exists() else None
    rng = np.random.default_rng(SEMENTE)
    linhas = []
    for d in des:
        for A in GRADE:
            linhas.append(celula(d, A, rng))
        g = pd.DataFrame([x for x in linhas if x["tic"] == d["tic"]])
        print(f"TIC {d['tic']} {'D10' if t26.loc[d['tic'], 'curvatura'] else '   '} n={len(d['E'])} barra SW mediana {np.median(d['sig'][~d['tess']]):5.1f} min | "
              f"completude(re): {' '.join(f'{c:.2f}' for c in g.completude_re)} | supressao: {' '.join(f'{s:.2f}' for s in g.supressao_med)}", flush=True)
    R = pd.DataFrame(linhas)
    R.to_parquet(BASE / "completude_secular_26.parquet", index=False)

    print(f"\n== por amplitude (grade {GRADE}; expectativa minha: falso alarme 1-4%; 50% entre 0,01-0,03; supressao >= 1,2 em 0,02, mediana 1,5-3 em 0,2)")
    res = []
    for A in GRADE:
        g = R[R.dpdt == A]
        L = {"dpdt": A, "compl_re_med": g.completude_re.median(), "compl_re_min": g.completude_re.min(), "n_compl_ge_0.5": int((g.completude_re >= 0.5).sum()),
             "compl_fixa_med": g.completude_fixa.median(), "sup_med": g.supressao_med.median(), "sup_p25": g.supressao_med.quantile(0.25), "sup_p75": g.supressao_med.quantile(0.75),
             "sup_max": g.supressao_med.max(), "tic_sup_max": int(g.loc[g.supressao_med.idxmax(), "tic"]), "vies_med": g.vies_med.median(), "frac_pgof_med": g.frac_pgof.median(), "frac_pgof_max": g.frac_pgof.max()}
        res.append(L)
        print(f"  |dP/dt| = {A:.3f}: completude(re) mediana {L['compl_re_med']:.2f}, min {L['compl_re_min']:.2f}, >= 0,5 em {L['n_compl_ge_0.5']} de 26 (fixa: mediana {L['compl_fixa_med']:.2f}) | "
              f"supressao mediana {L['sup_med']:.2f} (IQR {L['sup_p25']:.2f}-{L['sup_p75']:.2f}, max {L['sup_max']:.2f} em {L['tic_sup_max']}) | vies {L['vies_med']:.3f} | p_gof<0,05: mediana {L['frac_pgof_med']:.3f}, max {L['frac_pgof_max']:.3f}")
    pd.DataFrame(res).to_parquet(BASE / "completude_secular_26_resumo.parquet", index=False)
    # amplitude de 50% por alvo (interpolacao linear na grade) e o que a explica
    a50 = {}
    for tic, g in R.groupby("tic"):
        g = g.sort_values("dpdt")
        c, a = g.completude_re.values, g.dpdt.values
        a50[tic] = float(np.interp(0.5, c, a)) if c.max() >= 0.5 else np.nan
    a50 = pd.Series(a50)
    barra_sw = pd.Series({d["tic"]: float(np.median(d["sig"][~d["tess"]])) for d in des})
    n_ep = pd.Series({d["tic"]: len(d["E"]) for d in des})
    ok = a50.notna()
    rho_b = stats.spearmanr(a50[ok], barra_sw[ok]); rho_n = stats.spearmanr(a50[ok], n_ep[ok])
    print(f"\n  amplitude de completude 50%: mediana {a50.median():.3f} s/ano, {int(ok.sum())} de 26 atingem 0,5 na grade; nao atingem: {sorted(int(t) for t in a50[~ok].index)}")
    print(f"  Spearman(A50, barra SW mediana) = {rho_b.statistic:.2f} (p = {rho_b.pvalue:.3f}); Spearman(A50, n epocas) = {rho_n.statistic:.2f} (p = {rho_n.pvalue:.3f})")
    D10 = t26.index[t26.curvatura]
    print(f"  D10 em |dP/dt| = 0,02: completude {' '.join(f'{R[(R.tic == t) & (R.dpdt == 0.02)].completude_re.iloc[0]:.2f}' for t in D10)}")
    print(f"  232634196: supressao por amplitude {' '.join(f'{s:.2f}' for s in R[R.tic == 232634196].sort_values('dpdt').supressao_med)}")
    print(f"  -> {BASE / 'completude_secular_26.parquet'}, {BASE / 'completude_secular_26_resumo.parquet'}")
    resumo_por_alvo(R, des, t26)
