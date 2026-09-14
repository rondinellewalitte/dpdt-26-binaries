"""Epocas do SuperWASP (2006-2008) e a calibracao da cadeia que as torna usaveis.

PASSO 0 - VIABILIDADE, ja fechada e positiva para o TIC 142874476:

    1SWASP J032525.83-342938.6, a 0,15 arcsec do alvo
    9.185 pontos, 2006-07-02 a 2008-01-19 (566,8 d), duas temporadas
    dispersao robusta 12,9 mmag por ponto, erro declarado 12,5 mmag
    dobrando na efemeride linear do TESS: eclipse de 18,3 mmag, SNR 18,6

A profundidade de 18,3 mmag contra 34 mmag no TESS e o esperado: a abertura do
SuperWASP e de dezenas de segundos de arco e engole a vizinha, e a binagem em
1/60 de ciclo (35 min) borra um eclipse de 118 min.

O SISTEMA DE TEMPO E O PRIMEIRO CANDIDATO A VIES, E E DA ORDEM DO SINAL. O
arquivo entrega HJD em UTC; o TESS entrega BJD em TDB. A diferenca tem duas
partes:

    TDB - UTC em 2006-2008        +65,18 s = +1,086 min   (constante)
    helio contra baricentro        +-0,008 min            (desprezivel)
    -----------------------------------------------------------------
    total mediano                 +1,104 min

Isso e MAIOR que a barra de erro da epoca do setor 4 (0,65 min). Comparar um
HJD_UTC com um BJD_TDB sem converter poe 1,1 min de erro sistematico no ponto
de maior alavanca do ajuste - e o sinal inteiro que se mede tem 17,7 min.

PASSO 1 - A CALIBRACAO DA CADEIA. Nao ha verdade conhecida interna ao
SuperWASP. A externa parcial e sobre o INSTRUMENTO: binarias eclipsantes da
mesma epoca e do mesmo hemisferio cuja efemeride moderna seja precisa o
bastante para propagar 18 anos. O criterio de precisao nao e negociavel - a
mesma exigencia que impediu o cruzamento com arquivo historico antes do setor
31 - e por isso o calibrador precisa de setores TESS separados por ANOS, nao de
um setor so.

    sigma_P vindo de N ciclos com epocas de sigma_t: sigma_P ~ sigma_t / N
    propagando M ciclos para tras: erro ~ sigma_t * M / N

    Com s4-s105 (N ~ 1900 ciclos para P = 1,5 d) e sigma_t = 0,4 min, propagar
    mais M ~ 1600 ciclos ate 2007 custa ~0,34 min. Serve.
    Com um setor so, N ~ 18 ciclos, o mesmo M custa 35 min. Nao serve.

A mediana dos O-C dos calibradores e o VIES da cadeia; a dispersao e a
incerteza sistematica a somar em quadratura na epoca do alvo.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import config
import ids as ids_mod
import epocas_142874476 as ep

import numpy as np
import pandas as pd

CACHE = Path(config.CACHE) / "superwasp"
OUT_CAL = config.RESULTS / "superwasp_calibracao.parquet"
OUT_ALVO = config.RESULTS / "superwasp_alvo.parquet"

# SuperWASP-Sul fica no SAAO, Sutherland.
SAAO_LAT, SAAO_LON, SAAO_ALT = -32.3783, 20.8105, 1798.0

# Abaixo disto o ponto nao e fotometria, e erro de reducao.
MAG_MAX_VALIDA = 20.0


def baixar_curva(sourceid, force=False):
    """CSV do arquivo publico (CERIT/Masaryk). Falha alto."""
    import urllib.parse
    import urllib.request
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / (sourceid.replace(" ", "_") + ".csv")
    if p.exists() and not force and p.stat().st_size > 1000:
        return p
    url = "https://wasp.cerit-sc.cz/csv?object=" + urllib.parse.quote(sourceid)
    # Resposta vazia (138 bytes) e o arquivo respondendo com erro, nao ausencia
    # de dado - e sem repeticao ela entrava no lote como "nao mediu" (4 dos 62).
    # Mesmo desenho do cone search: tenta de novo com espera, e so depois falha.
    import time as _time
    ult = None
    for tentativa in range(4):
        try:
            urllib.request.urlretrieve(url, p)
            if p.stat().st_size >= 1000:
                return p
            ult = "veio vazio (%d bytes)" % p.stat().st_size
            p.unlink(missing_ok=True)
        except Exception as e:
            ult = str(e)
        _time.sleep(min(2.0 * (2 ** tentativa), 30.0))
    raise RuntimeError("FALHA [superwasp]: %s nao baixou em 4 tentativas: %s. Isso e "
                       "resposta ruim, nao ausencia de dado." % (sourceid, ult))


def hjd_utc_para_bjd_tdb(hjd, ra_deg, dec_deg):
    """Converte a escala E o referencial. Ver o cabecalho: vale +1,1 min.

    Desfaz a correcao heliocentrica, converte UTC->TDB e aplica a baricentrica.
    """
    from astropy import units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    loc = EarthLocation(lat=SAAO_LAT * u.deg, lon=SAAO_LON * u.deg,
                        height=SAAO_ALT * u.m)
    alvo = SkyCoord(ra_deg * u.deg, dec_deg * u.deg)
    t = Time(np.asarray(hjd, float), format="jd", scale="utc", location=loc)
    jd_utc = np.asarray(hjd, float) - t.light_travel_time(alvo, "heliocentric").jd
    t2 = Time(jd_utc, format="jd", scale="utc", location=loc)
    return t2.tdb.jd + t2.light_travel_time(alvo, "barycentric").jd


def carregar(sourceid, ra, dec):
    """Curva limpa, em BJD_TDB. Devolve (t_bjd, fluxo_relativo, sigma_mag)."""
    p = baixar_curva(sourceid)
    d = pd.read_csv(p)
    d.columns = ["hjd", "camera", "mag", "err"]
    n0 = len(d)
    d = d[(d["mag"] < MAG_MAX_VALIDA) & np.isfinite(d["mag"])]
    if len(d) < 500:
        raise RuntimeError("FALHA [superwasp]: %s ficou com %d pontos validos de "
                           "%d" % (sourceid, len(d), n0))
    # CLIPA OS BRILHANTES TAMBEM, e nao so os invalidos. O corte em mag > 20
    # tirava os quatro pontos de 29 mag e deixava passar fotometria em 11,64
    # contra mediana de 12,04 - fluxo 1,45, uma excursao POSITIVA de 45%. Com o
    # modelo de caixa (valores entre -1 e 0), um outlier brilhante dentro da
    # caixa resolve a profundidade para NEGATIVA: foi assim que a temporada de
    # 2007 saiu com -27.729 ppm e a epoca inteira desandou 290 min. Assimetrico
    # de proposito: eclipse e escuro, entao o lado fraco fica com folga maior.
    med0 = float(np.median(d["mag"]))
    s0 = float(np.median(np.abs(d["mag"] - med0)) * 1.4826)
    manter = (d["mag"] > med0 - 4 * s0) & (d["mag"] < med0 + 10 * s0)
    n1 = int((~manter).sum())
    d = d[manter]
    if len(d) < 500:
        raise RuntimeError("FALHA [superwasp]: %s ficou com %d pontos apos o "
                           "clip" % (sourceid, len(d)))
    t = hjd_utc_para_bjd_tdb(d["hjd"].values, ra, dec)
    med = float(np.median(d["mag"]))
    f = 10 ** (-0.4 * (d["mag"].values - med))
    sig = float(np.median(np.abs(d["mag"] - med)) * 1.4826)
    return t, f, sig


def largura_do_minimo(t, f, P, nb=200):
    """T14 medida na curva, nao suposta como fracao fixa do periodo.

    `t14 = min(0,10*P*24, 4 h)` errava nos dois sentidos: numa EA de P = 2,4 d
    dava uma caixa de 4 h onde o eclipse tem 2, e o ajuste devolvia
    profundidade de 73% - a caixa larga demais varre a linha de base junto. Com
    T14 medido isso nao acontece.
    """
    idx = np.clip((((t / P) % 1.0) * nb).astype(int), 0, nb - 1)
    prof = np.array([np.median(f[idx == b]) if (idx == b).sum() >= 3 else np.nan
                     for b in range(nb)])
    v = prof[np.isfinite(prof)]
    if len(v) < nb * 0.5:
        return np.nan
    topo = float(np.nanpercentile(prof, 95))
    d = topo - float(np.nanmin(prof))
    if d <= 0:
        return np.nan
    frac = float(np.mean(prof < topo - 0.5 * d))
    return float(np.clip(frac, 0.005, 0.35)) * P * 24


def medir_epoca(t, f, P, t0_previsto, span_min=None, t14_h=None, passos=3):
    """Epoca por perfil unitario com profundidade livre, em GRADE REFINADA.

    A grade fixa de 0,1 min devolvia sigma = 0,0 em setores com 28.000 pontos -
    a regiao de delta-chi2 = 1 fica MENOR que o passo, e o zero virava peso
    infinito no ajuste de periodo, que saia NaN. Agora a grade refina ate o
    sigma ser resolvido, e o piso e meio passo em vez de zero.

    O erro sai do delta-chi2 = 1 com o chi2 reescalado para chi2/dof = 1 no
    minimo. Isso corrige erro fotometrico otimista mas NAO corrige ruido
    vermelho - para isso existe `epoca_por_subconjunto`.
    """
    t = np.asarray(t, float)
    f = np.asarray(f, float)
    if t14_h is None or not np.isfinite(t14_h):
        t14_h = largura_do_minimo(t, f, P)
    if not np.isfinite(t14_h):
        return np.nan, np.nan, np.nan
    if span_min is None:
        # MEIO PERIODO, NAO 120 min FIXOS. Com janela fixa, um alvo cuja epoca
        # propagada erra mais que 120 min encosta na BORDA e devolve exatamente
        # -120,00 - numero que parece medida e e o limite da grade. Foi o que
        # produziu O-C de -120,00 e +81,79 na primeira calibracao. Meio periodo
        # e a janela sem ambiguidade: alem dela o proprio ciclo muda.
        span_min = P * 1440.0 / 2
    k = np.round((np.mean(t) - t0_previsto) / P)
    centro = t0_previsto + k * P
    y = f - 1.0
    n = 801
    for it in range(passos):
        grade = centro + np.linspace(-span_min, span_min, n) / 1440.0
        chi2 = np.empty(n)
        for i, g in enumerate(grade):
            m = ep.sim.modelo(t, g, P=P, t14_h=t14_h, prof=1.0) - 1.0
            den = float(m @ m)
            a = float(m @ y) / den if den > 1e-12 else 0.0
            chi2[i] = float(np.sum((y - a * m) ** 2))
        i = int(np.argmin(chi2))
        centro = float(grade[i])
        passo_min = 2 * span_min / (n - 1)
        dof = max(len(y) - 2, 1)
        esc = chi2[i] / dof
        c = chi2 / esc if esc > 0 else chi2
        alvo = c[i] + 1.0
        esq = grade[:i + 1][c[:i + 1] <= alvo]
        dir_ = grade[i:][c[i:] <= alvo]
        sig = ((dir_.max() - esq.min()) / 2 * 1440.0
               if len(esq) and len(dir_) else np.nan)
        # so refina enquanto o sigma nao estiver resolvido pela grade
        if not np.isfinite(sig) or sig > 3 * passo_min:
            break
        span_min = max(4 * passo_min, 3 * (sig if np.isfinite(sig) else passo_min))
    if np.isfinite(sig):
        sig = max(sig, passo_min / 2)       # piso: meio passo, nunca zero
    if abs(centro - (t0_previsto + np.round((centro - t0_previsto) / P) * P))             > 0.45 * P:
        return np.nan, np.nan, np.nan      # encostou na borda: nao e medida
    m = ep.sim.modelo(t, centro, P=P, t14_h=t14_h, prof=1.0) - 1.0
    den = float(m @ m)
    prof = float(m @ y) / den if den > 1e-12 else np.nan
    return centro, float(sig), prof


def noites_distintas(t, offset=0.5):
    """Quantas noites cobrem estes pontos. Diagnostico, nao unidade de erro."""
    return int(np.unique(np.floor(np.asarray(t, float) - offset)).size)


def amostras_independentes(t, dt_h):
    """Conta VISITAS: grupos separados por mais que `dt_h`. Diagnostico.

    NAO E A UNIDADE DO ERRO DE EPOCA, e a medida diz isso. A hipotese razoavel
    era que exposicoes agrupadas em minutos amostrassem a mesma fase da rampa e
    so contassem como uma. Testada no SuperWASP do alvo, onde a resposta e
    conhecida (barra medida de 1,88 min):

        N = pontos    8.809  ->  1,63 min   razao 1,15   <- reproduz
        N = visitas   1.005  ->  4,84 min   razao 0,39
        N = noites      186  -> 11,24 min   razao 0,17

    O PONTO E A UNIDADE, e as outras duas erram por 2,6x e 5,8x. Duas coisas
    explicam. Primeiro, dentro de uma visita o ruido e dominado por foton e a
    media dos pontos desce mesmo com sqrt(N) - nao ha sistematica compartilhada
    grande o bastante em poucos minutos. Segundo, e mais importante, a
    informacao de tempo esta na FORMA DA RAMPA amostrada ao longo do eclipse
    inteiro: sistematica comum de noite desloca o NIVEL do fluxo, e nivel nao
    carrega fase.

    O risco de ruido correlacionado nao desaparece - ele so nao se trata
    deflacionando N por uma regra a priori. Trata-se medindo:
    `epoca_por_subconjunto` compara a dispersao entre blocos independentes com a
    barra formal, e a maior vale. E o mesmo principio de sempre neste projeto -
    a barra sai do dado, nao de um modelo de quanto o dado deveria ser ruim.

    Fica como diagnostico porque a cadencia importa para OUTRA coisa: se um
    levantamento so visitar o alvo uma vez por noite, os pontos podem nunca
    cair no ingresso, e ai o problema nao e a barra, e a cobertura de fase.
    """
    t = np.sort(np.asarray(t, float))
    if len(t) == 0:
        return 0
    return int(1 + np.sum(np.diff(t) > dt_h / 24.0))


def epoca_por_subconjunto(t, f, P, t0_previsto, t14_h=None, n_sub=2,
                          modo="temporada", exigir_cobertura=True):
    """A barra de erro que o ruido VERMELHO exige, e nao a formal.

    Reescalar o chi2 para chi2/dof = 1 conserta erro fotometrico otimista, mas
    nao conserta correlacao: ruido vermelho atua entre pontos VIZINHOS dentro do
    proprio eclipse, e infla o erro do tempo de minimo mais do que o
    reescalonamento captura. O teste que fecha e o mesmo de sempre neste
    projeto - dispersao entre subconjuntos independentes.

    Aqui o corte natural sao as TEMPORADAS: o SuperWASP observa em blocos
    separados por meses, e dois blocos nao compartilham nem noite nem
    calibracao de campo. Devolve (t0_global, sigma_formal, sigma_entre_sub,
    lista de epocas por subconjunto).

    A regra e explicita: se a dispersao entre subconjuntos exceder a barra
    formal, e ela que vale.
    """
    t = np.asarray(t, float)
    f = np.asarray(f, float)
    t0g, sigf, _ = medir_epoca(t, f, P, t0_previsto, t14_h=t14_h)
    o = np.argsort(t)
    t, f = t[o], f[o]
    if modo == "ano":
        # Levantamento com cobertura sazonal continua por anos - ASAS-SN, ATLAS -
        # nao tem lacuna de 30 d que separe blocos, e o corte por gap devolveria
        # um grupo so. O ano e a granularidade que existe ali, e ele tambem
        # separa calibracao: refeitura de flat, troca de camera, mudanca de
        # pipeline acontecem entre temporadas anuais.
        ano = np.floor((t - t.min()) / 365.25).astype(int)
        grupos = [np.where(ano == a)[0] for a in np.unique(ano)]
        grupos = [g for g in grupos if len(g) > 100]
    else:
        # blocos separados por mais de 30 d = temporadas
        corte = np.where(np.diff(t) > 30.0)[0]
        grupos = np.split(np.arange(len(t)), corte + 1)
        grupos = [g for g in grupos if len(g) > 200]
    if len(grupos) < n_sub:                 # sem blocos: parte pela metade
        grupos = np.array_split(np.arange(len(t)), n_sub)
        grupos = [g for g in grupos if len(g) > 100]
    # O BLOCO SO ENTRA SE FOR MEDIVEL. O filtro acima e por CONTAGEM TOTAL de
    # pontos, que nao ve cobertura: 657 pontos num ano com 30 no eclipse e 5
    # bins de fase vazios passam por `len(g) > 100`. Medido no ATLAS do alvo,
    # e o resultado foi O-C anuais espalhados por +-950 min - praticamente
    # aleatorio dentro de +-P/2 - com sigma formal de 7 min citado ao lado. No
    # SuperWASP o WASP-42 tinha DOIS blocos com n_dentro = 0 e 12 bins vazios,
    # e `medir_epoca` devolvia epoca finita para eles. A barra de 951 min que a
    # guarda de dispersao pegou era o sintoma; isto e a causa.
    subs, rejeitados = [], []
    for g in grupos:
        cob = cobertura_de_fase(t[g], f[g], P, t0_previsto, t14_h)             if t14_h and np.isfinite(t14_h) else {"aprovado": True,
                                                  "motivo": "sem t14",
                                                  "n_dentro": -1,
                                                  "bins_vazios": -1}
        if exigir_cobertura and not cob["aprovado"]:
            rejeitados.append({"n": len(g), "n_dentro": cob["n_dentro"],
                               "bins_vazios": cob["bins_vazios"],
                               "motivo": cob["motivo"]})
            continue
        t0s, sgs, pfs = medir_epoca(t[g], f[g], P, t0_previsto, t14_h=t14_h)
        # Profundidade livre que sai <= 0 e o modelo dizendo "nao ha eclipse
        # aqui": a cobertura de fase passou, a grade devolveu um minimo de
        # chi2, e o que ele ajustou foi um CLAREAMENTO. Medido no lote de
        # dP/dt: 5 de 31 alvos tinham temporada assim, com dispersao entre
        # temporadas de 82 min contra 1,8 nos outros - e dois deles entravam
        # como "curvatura forte" (199688472 a 13 sigma) porque a parabola se
        # curvava ate dois pontos de barra enorme a 19 anos. E o caso do
        # secundario/primario trocado de `epocas_142874476.medir`, na versao
        # do arquivo: epoca finita de um eclipse que nao esta la.
        if not np.isfinite(t0s) or not (pfs > 0):
            rejeitados.append({"n": len(g), "n_dentro": cob["n_dentro"],
                               "bins_vazios": cob["bins_vazios"],
                               "motivo": "epoca nao medivel" if not np.isfinite(t0s)
                               else "profundidade ajustada <= 0 (%.4f): nao ha eclipse" % pfs})
            continue
        k = np.round((t0s - t0g) / P)
        subs.append({"n": len(g), "n_noites": noites_distintas(t[g]),
                     "t_medio": float(np.mean(t[g])),
                     "n_dentro": cob["n_dentro"],
                     "bins_vazios": cob["bins_vazios"],
                     "t0": t0s, "sigma_formal_min": sgs, "prof": pfs,
                     "desvio_min": float((t0s - (t0g + k * P)) * 1440)})
    if len(subs) >= 2:
        dv = np.array([s["desvio_min"] for s in subs])
        # DUAS quantidades, para dois desenhos: `sig_sub` e o erro da MEDIA
        # dos subconjuntos (std/sqrt n) - a barra da epoca combinada t0g,
        # que e o que `cadeia_oc.medir_alvo_de` (o controle) usa. Quem poe
        # CADA subconjunto como epoca propria precisa da dispersao de UM
        # subconjunto, std(dv), e ela vai em cada bloco. `oc_lote` usava
        # `sig_sub` por temporada: barra sqrt(n) pequena em 49 das 61
        # epocas SuperWASP dos 26, e "fator 2,73" na nota.
        dispersao = float(np.std(dv, ddof=1))
        for s in subs:
            s["dispersao_entre_blocos_min"] = dispersao
        sig_sub = float(dispersao / np.sqrt(len(dv)))
        status = "testado (%d blocos, %d rejeitados)" % (len(subs),
                                                        len(rejeitados))
    else:
        # NAO E A MESMA COISA QUE "testei e o ruido e pequeno". O chamador tem
        # que distinguir: cair na barra formal aqui e usar uma barra cuja
        # premissa - ausencia de correlacao - NAO FOI VERIFICADA. Devolver
        # `nan` sozinho fazia os tres chamadores fazerem exatamente isso, em
        # silencio, e a correcao acima AUMENTA a frequencia deste caso.
        sig_sub = np.nan
        status = ("NAO TESTAVEL: %d de %d blocos com cobertura de fase"
                  % (len(subs), len(subs) + len(rejeitados)))
    return t0g, float(sigf), sig_sub, subs, status

def semente_local(t, f, P, nb=300):
    """O minimo mais profundo da PROPRIA curva, em tempo absoluto.

    SEMEAR PELA EFEMERIDE PROPAGADA NAO SERVE AQUI. O periodo do catalogo tem 5
    casas; propagado 1.160 ciclos do setor 3 ao 105 isso vira dezenas de
    minutos de erro de fase, e o ajuste do setor antigo caia FORA do eclipse -
    sintoma: profundidade ajustada de 49.484 ppm num alvo cuja dobra direta da
    730.343 ppm nos quatro setores. Cada setor cobre poucas dezenas de ciclos,
    entao o minimo dele mesmo e inequivoco e nao depende de propagacao nenhuma.
    """
    t = np.asarray(t, float)
    f = np.asarray(f, float)
    fase = (t / P) % 1.0
    idx = np.clip((fase * nb).astype(int), 0, nb - 1)
    prof = np.array([np.median(f[idx == b]) if (idx == b).sum() >= 3 else np.nan
                     for b in range(nb)])
    if not np.isfinite(prof).any():
        return np.nan
    b = int(np.nanargmin(prof))
    fase_min = (b + 0.5) / nb
    # o ciclo observado mais proximo do centro do conjunto
    k = np.round((np.mean(t) / P) - fase_min)
    return float((k + fase_min) * P)


def escada_de_periodo(epocas, P0, max_res_min=3.0):
    """Resolve a CONTAGEM DE CICLOS antes de ajustar o periodo.

    O periodo de catalogo tem 5 casas decimais. Propagado do setor 3 ao 105 -
    1.160 ciclos - isso erra ate meio periodo, e atribuir o numero do ciclo com
    ele produz residuos de 24 a 203 min num ajuste que deveria fechar em menos
    de 1. O sintoma nao e ruido: e alias de ciclo.

    A saida e uma escada. Comeca pelo par mais PROXIMO no tempo, onde o numero
    de ciclos entre as duas epocas e inequivoco mesmo com o P grosseiro; refina
    P com esse par; usa o P refinado para numerar a proxima epoca mais proxima;
    repete. Cada degrau so estende a alavanca ate onde o P do degrau anterior
    aguenta.

    `epocas`: DataFrame com t0 e sig (minutos). Devolve (P, T0, sigma_P, E,
    residuos_min) ou None se algum degrau nao fechar dentro de `max_res_min`.
    """
    # A escada trabalha em ORDEM DE TEMPO, mas devolve E e residuos na ordem
    # em que o chamador passou as epocas. Antes devolvia na ordem interna,
    # sem dizer: no lote, TESS seguido de SuperWASP entrou desordenado e o E
    # de cada linha foi parar na linha errada - um setor a 1.756 ciclos saiu
    # com E = -5628, e dP/dt de -24.122 s/ano com barra de 0,004 min. Lixo
    # com cara de resultado, e a funcao "fechou" porque, na ordem interna,
    # ela estava certa.
    # (o nome e `ordem_entrada`, nao `ordem`: a escada ja usa `ordem` para a
    # sequencia em que os degraus sao adicionados, e a colisao devolveu E na
    # ordem de construcao dos degraus - vista falhando antes desta linha)
    ordem_entrada = np.argsort(epocas["t0"].values.astype(float), kind="stable")
    e = epocas.iloc[ordem_entrada].reset_index(drop=True)
    if len(e) < 2:
        return None
    t = e["t0"].values.astype(float)
    s = e["sig"].values.astype(float)
    # par mais proximo no tempo
    i = int(np.argmin(np.diff(t)))
    ordem = [i, i + 1]
    P = (t[i + 1] - t[i]) / max(round((t[i + 1] - t[i]) / P0), 1)
    T0 = t[i]
    while len(ordem) < len(t):
        # a proxima epoca mais proxima de alguma ja usada
        rest = [j for j in range(len(t)) if j not in ordem]
        j = min(rest, key=lambda k: min(abs(t[k] - t[m]) for m in ordem))
        cand = ordem + [j]
        E = np.round((t[cand] - T0) / P)
        A = np.vstack([E, np.ones_like(E)]).T
        w = 1 / (s[cand] / 1440) ** 2
        W = np.diag(w)
        try:
            c = np.linalg.solve(A.T @ W @ A, A.T @ W @ t[cand])
        except np.linalg.LinAlgError:
            return None
        r = (t[cand] - A @ c) * 1440
        if np.max(np.abs(r)) > max_res_min:
            return None                      # degrau nao fechou: alias nao resolvido
        P, T0 = float(c[0]), float(c[1])
        ordem = cand
    E = np.round((t - T0) / P)
    A = np.vstack([E, np.ones_like(E)]).T
    w = 1 / (s / 1440) ** 2
    W = np.diag(w)
    c = np.linalg.solve(A.T @ W @ A, A.T @ W @ t)
    cov = np.linalg.inv(A.T @ W @ A)
    r = (t - A @ c) * 1440

    # ALIAS GLOBAL: cada degrau pode fechar em < 3 min e o ENCADEAMENTO estar
    # errado. Um ciclo a mais num degrau intermediario propaga para todos os
    # posteriores e cada degrau individual continua fechando. O que pega isso e
    # perturbar a contagem de CADA epoca em +-1 e exigir que a nominal seja a
    # melhor por margem. Medido no TIC 129764561: nominal chi2 2,15e7 contra
    # 1,05e12 na melhor perturbacao - cinco ordens de grandeza.
    def _chi2(Ev):
        Av = np.vstack([Ev, np.ones_like(Ev)]).T
        try:
            cv = np.linalg.solve(Av.T @ W @ Av, Av.T @ W @ t)
        except np.linalg.LinAlgError:
            return np.inf
        rv = t - Av @ cv
        return float(rv @ W @ rv)
    if len(t) > 2:                       # com 2 epocas nao ha o que testar
        base = _chi2(E.astype(float))
        for i in range(len(t)):
            for dk in (-1.0, 1.0):
                Ev = E.astype(float).copy()
                Ev[i] += dk
                if _chi2(Ev) <= base:
                    return None          # alias global nao resolvido
    E_out, r_out = np.empty_like(E), np.empty_like(r)
    E_out[ordem_entrada], r_out[ordem_entrada] = E, r      # de volta a ordem do chamador
    return (float(c[0]), float(c[1]), float(np.sqrt(cov[0, 0])), E_out, r_out)


# -----------------------------------------------------------------------------
# PASSO 1 FALHOU COMO DESENHADO, E O MOTIVO E FISICO
#
# Dos 10 calibradores: 6 descartados por propagacao imprecisa (2 epocas TESS
# so, ou span < 2,1 anos), 1 por ter menos de 3 epocas. Sobraram 3:
#
#   TIC 129764561  O-C  +22,3 min   barra total 2,17
#   TIC 215336287  O-C  +71,0 min   barra total 1,32
#   TIC 401869914  O-C  +81,8 min   barra total 7,94
#
#   mediana +71,0   dispersao 31,7 min   -> 30x a maior barra
#
# ISSO NAO E RUIDO DE MEDIDA, E MUDANCA DE PERIODO REAL NOS PROPRIOS
# CALIBRADORES. Convertendo cada O-C em dP/dt implicado:
#
#   129764561  +0,051 s/ano      215336287  +0,171 s/ano
#   401869914  +0,202 s/ano      (o alvo:   +0,079 s/ano)
#
# Taxa tipica de Algol com transferencia de massa: 0,01 a 1 s/ano. Os tres
# estao dentro, e cercam o valor do alvo.
#
# O metodo supunha calibrador com relogio estavel. Binaria eclipsante de
# periodo curto NAO tem: e justamente a classe onde o periodo muda. Para
# extrair um vies de ~1 min de baixo de uma dispersao astrofisica de 31,7 min
# seriam precisos mais de 1.000 calibradores. Nao ha.
#
# O QUE A CALIBRACAO AINDA ESTABELECE, que nao e nada:
#   - nao ha erro GROSSO na cadeia. Um deslocamento de meia hora, de meio dia
#     ou de convencao de epoca apareceria como offset comum enorme, e os tres
#     caem na faixa que a fisica ja preve.
#   - os tres sao POSITIVOS, como o alvo. Sinal comum e fraco, mas nao
#     contradiz.
#
# O que ela NAO estabelece: uma barra sistematica de poucos minutos. Ela fica
# em aberto, e qualquer numero que dependa dela precisa dizer isso.
#
# A UNICA SISTEMATICA DE CADEIA IDENTIFICADA E QUANTIFICADA e a de sistema de
# tempo, +1,104 min, calculada de primeiros principios e aplicada em
# `hjd_utc_para_bjd_tdb`. Nao ha candidata conhecida na escala de 20 min.
#
# -----------------------------------------------------------------------------
# PASSOS 2 E 3 — a epoca de 2006-2008 e o que ela decide
#
#   t0 = 2454293,067494 BJD_TDB     8.809 pontos apos clip, 12,3 mmag
#   barra formal 1,88 min           dispersao entre temporadas 0,52 min
#   temporadas: +1,03 e -0,00 min, profundidades +13.651 e +12.346 ppm
#
# Aqui a dispersao entre temporadas ficou ABAIXO da formal, entao a formal e a
# barra honesta - a regra ja estava escrita para o caso contrario e nao
# precisou ser acionada.
#
#   E = -4699, 19,0 anos antes do setor 105
#
#   residuo no modelo LINEAR   dos 4 pontos TESS:  +20,74 min   (11 sigma)
#   residuo no modelo PARABOLA dos 4 pontos TESS:  -18,18 min   +- 1,88
#
# A CURVATURA E CONFIRMADA, A PARABOLA ESPECIFICA NAO. Com periodo constante o
# ponto de 2006 cairia sobre a reta; ele esta 20,7 min fora, a 11 sigma. Mas a
# parabola ajustada nos 4 pontos TESS previa +38,9 min para 2007 e o medido e
# +20,7 - ela SUPERESTIMA por fator 1,9 quando extrapolada de 7,65 para 19
# anos. Isso e o que se espera se o O-C nao for parabola pura: mudanca
# episodica de periodo, ou termo de tempo de luz por terceiro corpo somado a
# uma deriva menor.
#
# A RESSALVA QUE NAO SE FECHA. Um vies de cadeia de +20 min tornaria o residuo
# linear compativel com zero, e a calibracao do Passo 1 nao exclui isso -
# dispersao de 31,7 min nos calibradores. O argumento contra e fisico e nao
# estatistico: a unica sistematica de cadeia identificada vale 1,1 min, e nao
# ha mecanismo conhecido de 20 min. Enquanto o Passo 1 nao tiver calibradores
# de relogio estavel - transitos de Jupiteres quentes, nao binarias -, o
# numero de 2006 vale condicionado a isso, e o resultado dos 4 pontos TESS
# continua sendo o que se sustenta sozinho.


def cobertura_de_fase(t, f, P, t0, t14_h, n_bins=12):
    """A CHECAGEM QUE VEM ANTES DA FORMULA. Cadencia nao decide a barra - decide
    se ha rampa amostrada.

    `sigma_t = (rms/prof)*tau/sqrt(2*N_tau)` supoe que os pontos cobrem o
    ingresso e o egresso. Um levantamento que visite o alvo uma vez por noite
    pode acumular centenas de pontos NO FUNDO do eclipse e nenhum nas bordas, e
    a formula devolveria uma barra bonita para uma epoca que na pratica nao esta
    medida. A barra formal do ajuste tambem nao denuncia: com o fundo bem
    amostrado o chi2 tem minimo definido, so que largo e deslocavel.

    Padrao de referencia, o SuperWASP do alvo, que passa com folga: 564 pontos
    dentro do eclipse, 12 bins de fase, NENHUM vazio -
    [51,43,40,53,47,52,47,47,37,50,45,52].

    Devolve dict com n_dentro, histograma, bins vazios, profundidade MEDIDA (a
    diluicao depende da abertura de cada instrumento) e o veredito.
    """
    t = np.asarray(t, float)
    f = np.asarray(f, float)
    largura = (t14_h / 24.0) / P
    ph = ((t - t0) / P) % 1.0
    ph = np.where(ph > 0.5, ph - 1, ph)
    dentro = np.abs(ph) < largura / 2
    n = int(dentro.sum())
    h = np.histogram(ph[dentro], bins=n_bins,
                     range=(-largura / 2, largura / 2))[0] if n else np.zeros(n_bins, int)
    fora = ~dentro & (np.abs(ph) < 3 * largura)
    base = float(np.median(f[fora])) if fora.sum() > 10 else float(np.median(f))
    prof = base - float(np.median(f[dentro])) if n else np.nan
    vazios = int((h == 0).sum())
    # rms fora de eclipse, que e o que entra na formula
    rms = float(np.median(np.abs(f[fora] - base)) * 1.4826) if fora.sum() > 10 else np.nan
    ok = (n >= 50) and (vazios == 0) and np.isfinite(prof) and prof > 0
    return {"n_dentro": n, "histograma": h.tolist(), "bins_vazios": vazios,
            "prof_medida": prof, "prof_mmag": (2.5 * np.log10(1 / (1 - prof)) * 1000
                                               if np.isfinite(prof) and prof < 1 else np.nan),
            "rms_fora": rms, "rms_mmag": rms * 2.5 / np.log(10) * 1000 if np.isfinite(rms) else np.nan,
            "aprovado": bool(ok),
            "motivo": ("ok" if ok else
                       "poucos pontos no eclipse (%d)" % n if n < 50 else
                       "%d de %d bins de fase vazios: a rampa nao esta amostrada"
                       % (vazios, n_bins) if vazios else "profundidade nao medivel")}
