# -*- coding: utf-8 -*-
"""Rodada treze/catorze, bloco B: MEDIA e DIFERENCA primario/secundario no bloco TESS completo.

A pergunta. A Secao 5.5 mede um vagar de 2,4 a 8,5 min nos alvos contraditados e diz que nao
atribui origem. Manchas migrantes deslocam o primario e o secundario em sentidos OPOSTOS
(Tran et al. 2013, ~390 binarias Kepler, amplitudes de +-200 a 300 s, anticorrelacao atribuida a
manchas; Balaji et al. 2015 rastreia as longitudes; Kalimeris et al. 2002 da o mecanismo do
deslocamento do centro de luz). Um termo secular, um LTTE ou qualquer coisa que mova o sistema
inteiro deslocam os dois JUNTOS. Entao:

    M = (O-C_primario + O-C_secundario)/2      cancela a componente anticorrelada
    D =  O-C_primario - O-C_secundario         isola a componente anticorrelada

O metodo e o que Borkovits et al. (2016), ja citado, descreve para suprimir manchas.

ASSIMETRIA, declarada antes e para ir ao texto: a media cancela SO a componente anticorrelada.
Manchas polares, ou persistentes sobre o hemisferio nao eclipsado, deslocam os dois minimos no
mesmo sentido e sobrevivem a media. Logo "esta na diferenca -> ha componente de manchas" e
solido; "esta na media -> nao e mancha" NAO e, e nao deve ser escrito.

COMO. Epocas por setor de `jitter_primarios_26.parquet` (799 primarias e 799 secundarias, o
mesmo estimador e a mesma barra do desenho). Efemeride LINEAR ajustada ao proprio bloco de
primarias do alvo; o secundario usa a MESMA efemeride deslocada de P/2. Agrupamento por ANO
CIVIL, como a 5.5 - por setor o ruido intra-ano (0,2 a 2,2 min) domina. Por alvo: amplitude do
vagar em M e em D (maximo menos minimo das medias anuais, e chi2 contra constante), correlacao
de Pearson entre as medias anuais do primario e do secundario, e o coeficiente quadratico
ajustado a M sob o modelo de ruido ADMITIDO (o criterio de chi2/nu da 5.3.2), comparado com o do
primario sozinho.

CALIBRACAO DO PROPRIO TESTE: V527 Dra (TIC 424461577), o unico alvo da amostra em que os DOIS
sinais sao conhecidos - orbita LTTE publicada (que move os dois minimos juntos) e manchas
migrantes documentadas (Ceki et al. 2024). Se o teste separar os dois ali, ele tem credito nos
demais; se nao separar, o que ele disser nos outros vale menos.

EXPECTATIVA (escrita e commitada ANTES de rodar)

(a) GERAL. Se o vagar dos contraditados for de manchas, a amplitude em D deve ser comparavel ou
maior que a excursao do primario (2,4 a 8,5 min) e as medias anuais de primario e secundario
devem anticorrelacionar (r < 0). Espero componente anticorrelada presente em PELO MENOS DOIS dos
oito testaveis - a literatura de Kepler a encontra em boa parte das binarias de contato, e tres
dos nossos ja mostram deriva significativa do secundario no desenho (Secao 5.3.1). Espero
tambem que M seja MAIS QUIETO que o primario nesses alvos (amplitude menor), e que o coeficiente
quadratico de M fique mais proximo de zero que o do primario onde houver anticorrelacao.

(b) V527 Dra. Espero D com vagar de alguns minutos (as manchas de Ceki et al.) e M mantendo o
sinal coerente do LTTE - ou seja, M NAO deve ficar mais quieto que o primario neste alvo, e a
correlacao anual deve ficar entre -0,5 e +0,5 (os dois efeitos somados). Se D sair plano, o
teste nao tem sensibilidade aqui e o credito nao e ganho; digo isso no relatorio.

(c) 126945917 / HD 202042 - o unico com comparacao externa direta. Pawar et al. (2026) reportam
anticorrelacao primario/secundario DENTRO dos setores TESS; a Secao 5.3.1 da, para o desenho,
deriva do secundario de +0,50 +- 0,47 min/ano (compativel com zero); e no estado F o alvo cruzou
p_curv de 0,081 para 0,004 sem virar deteccao (a porta adversarial o rejeita, p_adv 0,67). O
bloco completo dele tem SEIS setores uteis (1, 27, 28, 68, 95, 105), de 2018,6 a 2026,5: um por
ano, com barras de ~0,5 min.
    O que espero ver: como a nossa medida e ENTRE setores (anos) e a de Pawar et al. e DENTRO de
    setores (ciclos, ~27 d), as duas nao medem a mesma escala, e nao ha razao forte para que a
    anticorrelacao intra-setor apareca aqui. Espero D com dispersao da ordem de 1 a 3 min e
    correlacao anual mal determinada (|r| < 0,7 com 6 pontos nao e significativa), e M parecido
    com o primario.
    O que cada desfecho diria: (i) D com amplitude >> barras E r negativo forte -> a componente
    de manchas persiste na escala de anos, o que ESTENDE Pawar et al. para alem do setor e
    sugere que a curvatura recem-aparecida do alvo (p_curv 0,004) e de manchas, nao de periodo;
    (ii) D compativel com ruido -> nao contradiz Pawar et al. (escalas diferentes), e o teste
    simplesmente nao tem alcance nesse alvo com seis epocas; (iii) D grande mas r positivo ->
    algo move os dois minimos juntos e a leitura de manchas nao serve. Nenhum desses desfechos
    autoriza escrever que o alvo tem ou nao tem manchas a partir daqui: com seis epocas o que se
    pode dizer e sobre a escala de anos, e isso e o que a nota vai dizer.

(d) SANIDADE, antes de qualquer leitura: as epocas primarias por setor tem de reproduzir as do
desenho nos setores em que os dois existem (o mesmo controle de 10^-3 min da 5.3.2), e a
diferenca mediana d = t_sec - t_prim - P/2 por alvo tem de bater com `secundarios_26.parquet`
dentro da barra. Se nao bater, parar.

EXPECTATIVA DO DENOMINADOR E DA COBERTURA (rodada dezenove, escrita antes de rodar)

(e) A EXCURSAO EM TORNO DA QUADRATICA, `amp_P_quad_min`. O vagar contra o qual o limite e
comparado nao pode conter o termo secular do proprio bloco: senao a fracao mede quanto a
componente anticorrelada explica de uma parabola, e nao do vagar. Espero, alvo a alvo:
  - 232634196 (bloco -0,65 s/ano, o maior da amostra): 12,36 -> 6 a 9 min. A fracao do limite
    (5,75 min) sobe de 47% para 64-96%, e o alvo SAI do grupo "abaixo da metade".
  - 198408416 (bloco +0,047): 2,76 -> 2,0 a 2,6; fracao 40-52%, no limiar do grupo.
  - 230386284 (+0,0091), 329246824, 392536812 (blocos pequenos): mudanca abaixo de 20%.
  - 424461577 / V527 Dra (LTTE de 6,6 min): 7,34 -> 5 a 7; a fracao do limite (0,49) fica em
    7-10%, e continua sendo o caso mais apertado da amostra.
  - Contagem "abaixo da metade": hoje 2; espero 0 a 2 depois. Se sair 3 ou mais, a mudanca esta
    na direcao que fortalece o resultado e vai ao relatorio como desvio antes de qualquer texto.
Contra a 5.5: `oc_excursao_min` da 5.5 e a excursao anual contra a quadratica ajustada a
DESENHO + ARQUIVO, nao ao bloco. Nos alvos em que as duas janelas concordam (230386284,
377253090) espero as duas quantidades a menos de 30% uma da outra; nos contraditados e em V527
Dra espero desacordo maior, porque a discordancia entre as janelas E o resultado da 5.5. Se o
desacordo for grande tambem nos nao contraditados, as duas quantidades nao sao a mesma coisa em
lugar nenhum e cada uma tem de aparecer com o seu nome.

(f) COBERTURA DO LIMITE. O modelo do limite e D = c + ruido + excesso branco: uma deriva de D
entra no excesso, logo o limite publicado JA a cobre. Ajustando um termo linear em D junto com o
excesso, espero o excesso cair nos quatro alvos com deriva significativa (198408416, 229914020,
230386284, 232634196) - reducao de 20 a 60% - e menos de 10% nos demais. Se cair pouco nos quatro,
a deriva nao domina o excesso; se cair muito, o limite publicado e conservador por essa margem, e
e isso que a 5.3.1 tem de escrever. Em nenhuma hipotese o limite COM deriva ajustada substitui o
publicado: ele responde outra pergunta (quanto sobra de branco depois de tirar a deriva), e a
pergunta do artigo e quanto cabe de componente anticorrelada no total.

Uso:
    python src/media_diferenca_oc_d9.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
import ruido_admitido_26 as RA  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="backslashreplace")
except (AttributeError, ValueError):
    pass

BASE = config.DATA / "orquestra" / "oc_lote"
S_POR_ANO = 86400.0 * 365.25
MIN_SETORES = 20                      # bloco "completo" para este teste
EXTERNO = 126945917                   # HD 202042: comparacao com Pawar et al. (2026)
CALIBRACAO = 424461577                # V527 Dra: LTTE publicado + manchas documentadas


def ano_civil(t_btjd):
    return 2000.0 + (np.asarray(t_btjd, float) + 2457000.0 - 2451544.5) / 365.25


def series(j, tic, P):
    """O-C por setor do primario e do secundario, na MESMA efemeride linear (ajustada as
    primarias), em minutos; devolve um quadro por setor com M e D."""
    d = j[(j.tic == tic) & j.ok]
    pri = d[d.minimo == "primario"].set_index("setor").sort_index()
    sec = d[d.minimo == "secundario"].set_index("setor").sort_index()
    comuns = pri.index.intersection(sec.index)
    pri, sec = pri.loc[comuns], sec.loc[comuns]
    E = np.round((pri.t0_btjd.values - pri.t0_btjd.values[0]) / P)
    w = 1.0 / pri.sig_min.values ** 2
    A = np.vstack([E, np.ones_like(E)]).T
    coef = np.linalg.solve(A.T @ (w[:, None] * A), A.T @ (w * pri.t0_btjd.values))   # P e T0 do bloco
    efem = coef[1] + coef[0] * E
    oc_p = (pri.t0_btjd.values - efem) * 1440.0
    # o secundario cai meio periodo depois do primario do MESMO ciclo (ou meio antes, conforme
    # a ancora): o ciclo e resolvido pelo arredondamento, e o residuo e o O-C
    n = np.round((sec.t0_btjd.values - efem - coef[0] / 2) / coef[0])
    oc_s = (sec.t0_btjd.values - (efem + coef[0] / 2 + n * coef[0])) * 1440.0
    return pd.DataFrame({"tic": tic, "setor": comuns, "t_btjd": pri.t0_btjd.values,
                         "ano": ano_civil(pri.t0_btjd.values), "P_bloco_d": coef[0],
                         "oc_p": oc_p, "sig_p": pri.sig_min.values,
                         "oc_s": oc_s, "sig_s": sec.sig_min.values,
                         "M": (oc_p + oc_s) / 2.0, "sig_M": np.hypot(pri.sig_min.values, sec.sig_min.values) / 2.0,
                         "D": oc_p - oc_s, "sig_D": np.hypot(pri.sig_min.values, sec.sig_min.values)})


def por_ano(s, col, sig):
    """Media ponderada por ano civil, com a barra da media e a dispersao interna."""
    out = []
    for ano, g in s.groupby(s.ano.astype(int)):
        w = 1.0 / g[sig].values ** 2
        m = float((w * g[col].values).sum() / w.sum())
        out.append({"ano": int(ano), "n": len(g), "valor": m, "sigma": float(1.0 / np.sqrt(w.sum())),
                    "disp": float(np.std(g[col].values, ddof=1)) if len(g) > 1 else np.nan})
    return pd.DataFrame(out)


def amplitude(a):
    """Amplitude do vagar anual e chi2 contra constante."""
    if len(a) < 2:
        return dict(amp=np.nan, chi2r=np.nan, p=np.nan)
    w = 1.0 / a.sigma.values ** 2
    m = (w * a.valor.values).sum() / w.sum()
    chi2 = float((w * (a.valor.values - m) ** 2).sum())
    return dict(amp=float(a.valor.max() - a.valor.min()), chi2r=chi2 / (len(a) - 1),
                p=float(stats.chi2.sf(chi2, len(a) - 1)))


def amplitude_quadratica(s, col, sig, P):
    """Excursao das medias anuais DEPOIS de removida uma quadratica ajustada ao proprio bloco.

    O `amplitude()` acima e residual a efemeride LINEAR do bloco, entao inclui a curvatura do
    proprio bloco - e pedir que manchas expliquem um termo parabolico e pedir a coisa errada
    (rodada dezenove). Esta e a quantidade que serve de denominador da fracao da 5.3.1."""
    t = s.t_btjd.values
    E = (t - t[0]) / P
    X = np.vander(E, 3)
    w = 1.0 / s[sig].values ** 2
    coef = np.linalg.solve(X.T @ (w[:, None] * X), X.T @ (w * s[col].values))
    res = s[col].values - X @ coef
    a = por_ano(s.assign(_r=res), "_r", sig)
    return amplitude(a)


def limite_com_deriva(s):
    """O mesmo limite de perfil, com um termo LINEAR em D ajustado junto (rodada dezenove).

    Nao substitui `limite_anticorrelada`: serve para dizer QUANTO do excesso que o limite cobre e
    deriva. O limite publicado nao tem termo de deriva, logo uma deriva real esta dentro dele."""
    D, sg = s.D.values, s.sig_D.values
    x = (s.t_btjd.values - s.t_btjd.values[0]) / 365.25

    def m2lnL(e2):
        v = sg ** 2 + max(e2, 0.0)
        w = 1.0 / v
        A = np.vstack([x, np.ones_like(x)]).T
        c = np.linalg.solve(A.T @ (w[:, None] * A), A.T @ (w * D))
        r = D - A @ c
        return float(((r ** 2) / v + np.log(v)).sum())

    grade = np.concatenate([[0.0], np.geomspace(1e-4, 100.0, 400)])
    val = np.array([m2lnL(e2) for e2 in grade])
    i = int(np.argmin(val))
    acima = np.where((grade > grade[i]) & (val > val[i] + 2.71))[0]
    e95 = float(np.sqrt(grade[acima[0]])) if len(acima) else float(np.sqrt(grade[-1]))
    return {"sigma_extra_D_min_com_deriva": float(np.sqrt(grade[i])), "sigma_extra_D_95_com_deriva_min": e95,
            "delta_pico_a_pico_95_com_deriva_min": 2.83 * e95 / 2.0}


def quadratica(s, col, sig, P):
    """Coeficiente quadratico da serie, em s/ano, sob a diagonal e sob o GP agrupado, com o
    criterio de modelo admitido (chi2/nu no intervalo de 95%) para dizer qual vale."""
    t = s.t_btjd.values
    E = (t - t[0]) / P
    X = np.vander(E, 3)
    saida = {}
    for nome, extra in (("diagonal", 0.0), ("GP_agrupado", 0.835)):
        sg = np.hypot(s[sig].values, extra)
        w = 1.0 / (sg / 1440.0) ** 2
        cov = np.linalg.inv(X.T @ (w[:, None] * X))
        coef = cov @ (X.T @ (w * (s[col].values / 1440.0)))
        res = s[col].values / 1440.0 - X @ coef
        chi2 = float(res @ (w * res))
        nu = len(E) - 3
        ok, lo, hi = RA.admitido(chi2 / nu, nu, 0.95)
        saida[nome] = {"dPdt": 2 * coef[0] / P * S_POR_ANO, "s": 2 * np.sqrt(cov[0, 0]) / P * S_POR_ANO,
                       "chi2r": chi2 / nu, "admitido": bool(ok), "int": [lo, hi]}
    return saida


def limite_anticorrelada(s):
    """Excesso branco da serie D sobre as proprias barras, por maxima verossimilhanca, e o limite
    superior de 95% - convertido em limite sobre a componente ANTICORRELADA (metade de D).

    O modelo e D_i = c + ruido(sig_D_i) + extra(sigma_e): -2 ln L = sum[ (D-c)^2/(s^2+e^2) +
    ln(s^2+e^2) ]. O limite de 95% vem do perfil (Delta(-2 ln L) = 2,71 num so lado)."""
    D, sg = s.D.values, s.sig_D.values

    def m2lnL(e2):
        v = sg ** 2 + max(e2, 0.0)
        c = (D / v).sum() / (1.0 / v).sum()
        return float((((D - c) ** 2) / v + np.log(v)).sum())

    grade = np.concatenate([[0.0], np.geomspace(1e-4, 100.0, 400)])
    val = np.array([m2lnL(e2) for e2 in grade])
    i = int(np.argmin(val))
    e_hat = float(np.sqrt(grade[i]))
    acima = np.where((grade > grade[i]) & (val > val[i] + 2.71))[0]
    e95 = float(np.sqrt(grade[acima[0]])) if len(acima) else float(np.sqrt(grade[-1]))
    return {"sigma_extra_D_min": e_hat, "sigma_extra_D_95_min": e95,
            "delta_rms_95_min": e95 / 2.0, "delta_pico_a_pico_95_min": 2.83 * e95 / 2.0}


def main():
    j = pd.read_parquet(BASE / "jitter_primarios_26.parquet")
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    sec26 = pd.read_parquet(BASE / "secundarios_26.parquet").set_index("tic")

    n_ok = j[j.ok].groupby(["tic", "minimo"]).size().unstack(fill_value=0)
    alvos = sorted(set(n_ok.index[(n_ok.primario >= MIN_SETORES) & (n_ok.secundario >= MIN_SETORES)]) | {EXTERNO})
    print(f"== B: media e diferenca primario/secundario no bloco completo ({len(alvos)} alvos: "
          f"{sum(1 for a in alvos if a != EXTERNO)} com >= {MIN_SETORES} setores, mais {EXTERNO})")

    tudo, resumo = [], []
    for tic in alvos:
        P = float(t26.loc[tic, "P"])
        s = series(j, tic, P)
        if len(s) < 4:
            print(f"  TIC {tic}: {len(s)} setores com os dois minimos - fora")
            continue
        # sanidade (d): a diferenca mediana bate com secundarios_26?
        d_med = float(np.median(s.D.values)) / 2.0 * -1.0        # D = p - s; secundarios_26 usa s - p - P/2
        d_ref = float(sec26.loc[tic, "d_mediano_min"]) if tic in sec26.index else np.nan
        tudo.append(s)

        aM, aD = por_ano(s, "M", "sig_M"), por_ano(s, "D", "sig_D")
        aP = por_ano(s, "oc_p", "sig_p")
        aS = por_ano(s, "oc_s", "sig_s")
        r = float(np.corrcoef(aP.valor.values, aS.valor.values)[0, 1]) if len(aP) > 2 else np.nan
        # a correlacao POR SETOR, com o seu p, e a deriva de D: e ela que pega o caso em que os
        # dois minimos se movem em sentidos opostos setor a setor sem que D, com a barra larga
        # do secundario, chegue a rejeitar uma constante
        rs, prs = stats.pearsonr(s.oc_p.values, s.oc_s.values)
        wD = 1.0 / s.sig_D.values ** 2
        xD = s.t_btjd.values / 365.25
        AD = np.vstack([xD, np.ones_like(xD)]).T
        cD = np.linalg.solve(AD.T @ (wD[:, None] * AD), AD.T @ (wD * s.D.values))
        covD = np.linalg.inv(AD.T @ (wD[:, None] * AD))
        qM, qP = quadratica(s, "M", "sig_M", P), quadratica(s, "oc_p", "sig_p", P)
        ampM, ampD, ampP = amplitude(aM), amplitude(aD), amplitude(aP)
        lim = limite_anticorrelada(s)
        lim_d = limite_com_deriva(s)
        ampPq = amplitude_quadratica(s, "oc_p", "sig_p", P)
        deriva_sig = abs(cD[0]) > 2 * np.sqrt(covD[0, 0])
        linha = {"tic": int(tic), "nome": t26.loc[tic, "nome"], "n_setores": len(s), "n_anos": len(aM),
                 "sinal_r": ("+" if rs > 0 else "−"), "deriva_d_significativa": bool(deriva_sig),
                 **lim,
                 "amp_P_min": ampP["amp"], "amp_M_min": ampM["amp"], "amp_D_min": ampD["amp"],
                 **lim_d,
                 "amp_P_quad_min": ampPq["amp"], "chi2r_P_quad": ampPq["chi2r"], "p_const_P_quad": ampPq["p"],
                 "chi2r_P": ampP["chi2r"], "chi2r_M": ampM["chi2r"], "chi2r_D": ampD["chi2r"],
                 "p_const_D": ampD["p"], "r_anual": r, "r_setor": float(rs), "p_r_setor": float(prs),
                 "deriva_D_min_por_ano": float(cD[0]), "s_deriva_D": float(np.sqrt(covD[0, 0])), "d_mediano_min": -2 * d_med, "d_ref_26": d_ref,
                 "dPdt_P_diag": qP["diagonal"]["dPdt"], "s_P_diag": qP["diagonal"]["s"], "chi2r_P_diag": qP["diagonal"]["chi2r"],
                 "dPdt_M_diag": qM["diagonal"]["dPdt"], "s_M_diag": qM["diagonal"]["s"], "chi2r_M_diag": qM["diagonal"]["chi2r"],
                 "dPdt_P_gp": qP["GP_agrupado"]["dPdt"], "s_P_gp": qP["GP_agrupado"]["s"], "chi2r_P_gp": qP["GP_agrupado"]["chi2r"],
                 "dPdt_M_gp": qM["GP_agrupado"]["dPdt"], "s_M_gp": qM["GP_agrupado"]["s"], "chi2r_M_gp": qM["GP_agrupado"]["chi2r"],
                 "admitido_P": [k for k, v in qP.items() if v["admitido"]], "admitido_M": [k for k, v in qM.items() if v["admitido"]]}
        resumo.append(linha)
        marca = " [CALIBRACAO]" if tic == CALIBRACAO else (" [EXTERNO]" if tic == EXTERNO else "")
        print(f"  TIC {tic} {linha['nome']:<18s}{marca} {len(s)} setores em {len(aM)} anos | "
              f"amplitude anual: P {ampP['amp']:.2f}, M {ampM['amp']:.2f}, D {ampD['amp']:.2f} min | "
              f"D constante? p = {ampD['p']:.3g} | r(P,S) setor {rs:+.2f} (p {prs:.1g}) | deriva de D {cD[0]:+.2f} ± {np.sqrt(covD[0, 0]):.2f} min/ano")
        print(f"      dP/dt do primario {qP['diagonal']['dPdt']:+.4f} ± {qP['diagonal']['s']:.4f} (chi2/nu {qP['diagonal']['chi2r']:.2f})"
              f" | da media {qM['diagonal']['dPdt']:+.4f} ± {qM['diagonal']['s']:.4f} (chi2/nu {qM['diagonal']['chi2r']:.2f})"
              f" | admitidos: P {qP and [k for k, v in qP.items() if v['admitido']]}, M {[k for k, v in qM.items() if v['admitido']]}")

    S = pd.concat(tudo, ignore_index=True)
    S.to_parquet(BASE / "media_diferenca_oc_d9_setores.parquet", index=False)
    R = pd.DataFrame(resumo).set_index("tic", drop=False)
    R.to_parquet(BASE / "media_diferenca_oc_d9.parquet", index=False)
    # DOIS criterios, e os dois vao ao relatorio, porque o segundo foi escrito DEPOIS de ver o
    # resultado do primeiro (que deixava de fora o alvo em que os dois minimos andam em sentidos
    # opostos setor a setor sem que D, com a barra larga do secundario, rejeite uma constante).
    # Trocar um corte olhando para os sobreviventes e o que a regra da casa proibe: entao nenhum
    # dos dois e "o" criterio ate o autor decidir - os dois ficam lado a lado.
    anticorr_ini = R[(R.r_anual < 0) & (R.p_const_D < 0.05)]        # escrito ANTES de rodar
    anticorr = R[(R.r_setor < 0) & (R.p_r_setor < 0.05)]            # escrito DEPOIS de rodar
    print(f"\n== resumo, com os DOIS criterios (o segundo foi escrito depois de rodar):")
    print(f"   (i)  declarado antes: r anual < 0 E D nao constante a 5%  -> {len(anticorr_ini)} de {len(R)}: {sorted(anticorr_ini.tic.tolist())}")
    print(f"   (ii) escrito depois:  r por setor < 0 com p < 0,05        -> {len(anticorr)} de {len(R)}: {sorted(anticorr.tic.tolist())}")
    print(f"   nos dois: {sorted(set(anticorr_ini.tic) & set(anticorr.tic))}")
    print("")
    print("== C2: limite superior (95%) sobre a componente anticorrelada, e o vagar do primario")
    for tic in sorted(R.index):
        r = R.loc[tic]
        marca = "  <- limite ABAIXO do vagar" if r.delta_pico_a_pico_95_min < r.amp_P_min else ""
        print(f"   {tic} {r['nome']:<20s} sinal r {r['sinal_r']}  deriva em d {'SIM' if r['deriva_d_significativa'] else 'nao'} | "
              f"sigma_extra(D) {r['sigma_extra_D_min']:.2f} (95%: {r['sigma_extra_D_95_min']:.2f}) min | "
              f"delta_95 {r['delta_rms_95_min']:.2f} rms, {r['delta_pico_a_pico_95_min']:.2f} pico a pico | "
              f"vagar do primario {r['amp_P_min']:.2f} min{marca}")
    print(f"   M mais quieto que o primario em {int((R.amp_M_min < R.amp_P_min).sum())} de {len(R)}")
    (BASE / "media_diferenca_oc_d9.json").write_text(json.dumps({
        "alvos": resumo, "anticorrelados_criterio_previo": sorted(int(x) for x in anticorr_ini.tic),
        "anticorrelados_criterio_posterior": sorted(int(x) for x in anticorr.tic),
        "criterio_previo": "r anual < 0 e D nao constante a 5% (escrito antes de rodar)",
        "criterio_posterior": "r(P,S) por setor < 0 com p < 0,05 (escrito depois de rodar)",
        "criterio_anticorrelacao": "r(P,S) por setor < 0 com p < 0,05",
        "calibracao": int(CALIBRACAO), "externo": int(EXTERNO)}, indent=1, default=float), encoding="utf-8")
    print(f"\n  -> {BASE / 'media_diferenca_oc_d9.parquet'} e _setores.parquet")


if __name__ == "__main__":
    main()
