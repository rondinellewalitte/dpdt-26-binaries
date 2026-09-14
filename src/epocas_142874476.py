"""Epocas do TIC 142874476 em quatro aglomerados, e o que elas distinguem.

O PROBLEMA QUE ISTO RESOLVE. Com tres aglomerados - SuperWASP (2006-2008),
setor 31 (2020), setor 105 (2026) - e tres parametros, a parabola passa pelos
tres pontos por construcao. Uma senoide de tempo de luz por terceiro corpo
tambem passa. Nao ha grau de liberdade e portanto nao ha teste. As FFI dos
setores 4 (2018) e 97 (2025) caem NOS INTERVALOS, e e isso que separa as duas
hipoteses: parabola tem curvatura de sinal unico, senoide inverte.

A SIMULACAO VEIO ANTES (`src/simular_ffi.py`). Cadencia de 30 min poe so ~4
pontos num eclipse de 1,967 h, e o ponto e uma INTEGRAL sobre a janela, nao uma
amostra - o borramento e distorcao deterministica, nao ruido. Medido em 200
realizacoes com fase sorteada, gap de downlink e ruido correlacionado: a epoca
media do aglomerado sai com espalhamento de 0,47 min no cenario tipico e 0,82
no pessimista, contra o criterio de 2 min. O controle negativo de 2 h de
cadencia degrada para 5,2 min, o que mostra que a simulacao E CAPAZ de reprovar.

TRES REDUCOES INDEPENDENTES NO SETOR 4 - QLP, GSFC-eleanor-lite e TGLC. Elas
compartilham os pixels e nao compartilham a fotometria: mascara, fundo e
detrending sao diferentes. Concordancia entre elas nao prova ausencia de
sistematica de FFI, mas discordancia denuncia. O espalhamento entre as tres e
reportado junto e entra como piso da incerteza da epoca.
"""
import glob
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import config
import ids as ids_mod
import simular_ffi as sim

import numpy as np
import pandas as pd

TIC = 142874476
P_REF = 1.4740475            # d, reanalise com o setor 31
FFI = Path(config.CACHE) / "ffi_142874476"
OUT = config.RESULTS / "epocas_142874476.parquet"


def _norm(t, f, q=None):
    ok = np.isfinite(t) & np.isfinite(f)
    if q is not None:
        ok &= (np.asarray(q) == 0)
    t, f = np.asarray(t)[ok], np.asarray(f)[ok]
    base = np.percentile(f, 95)
    return t, f / base


def carregar():
    """Uma entrada por (setor, reducao). Falha alto se um arquivo esperado sumir."""
    from astropy.io import fits
    fontes = []
    for p in sorted(glob.glob(str(FFI / "mastDownload" / "**" / "*.fits"),
                              recursive=True)):
        nome = os.path.basename(p)
        with fits.open(p) as h:
            hdr, d = h[0].header, h[1].data
            setor = int(hdr.get("SECTOR"))
            cols = [c.upper() for c in d.columns.names]
            if "KSPSAP_FLUX" in cols:
                red, fl = "QLP", d["KSPSAP_FLUX"]
            elif "DET_FLUX" in cols:
                red, fl = "QLP", d["DET_FLUX"]
            elif "CORR_FLUX" in cols:
                red, fl = "eleanor", d["CORR_FLUX"]
            elif "CAL_APER_FLUX" in cols:
                red, fl = "TGLC", d["cal_aper_flux"]
            elif "PDCSAP_FLUX" in cols:
                # TESS-SPOC e SPOC compartilham o esquema de colunas; o que
                # separa e a cadencia. Chamar os dois de "SPOC" faria a
                # deduplicacao por (setor, reducao) descartar um deles em
                # silencio - foi assim que o s105 apareceu duplicado antes.
                cad_s = float(h[1].header.get("TIMEDEL", 0)) * 86400
                red = "SPOC" if cad_s < 300 else "TESS-SPOC"
                fl = d["PDCSAP_FLUX"]
            else:
                raise RuntimeError("FALHA [epocas]: %s sem coluna de fluxo "
                                   "reconhecida: %s" % (nome, cols))
            tempo = d["TIME"] if "TIME" in d.columns.names else d["time"]
            q = None
            for k in ("QUALITY", "TESS_flags"):
                if k in d.columns.names:
                    q = d[k]
                    break
            cad = float(h[1].header.get("TIMEDEL", np.nan)) * 1440
            t, f = _norm(tempo, fl, q)
            if len(t) < 200:
                raise RuntimeError("FALHA [epocas]: %s com so %d pontos validos"
                                   % (nome, len(t)))
            fontes.append({"setor": setor, "reducao": red, "arquivo": nome,
                           "cadencia_min": cad, "n": len(t), "t": t, "f": f})
    return fontes


def medir(t, f, P=P_REF, span_min=120.0, n=2401, semente=None, t14_h=None):
    """Epoca media com PROFUNDIDADE LIVRE em cada t0 tentado.

    Reducao diferente dilui diferente - a QLP e a TGLC nao entregam a mesma
    profundidade para a mesma estrela. Fixar a profundidade do modelo faria o
    chi2 medir desajuste de amplitude junto com desajuste de fase. Com a
    profundidade resolvida analiticamente em cada t0, o chi2 mede so fase.

    A SEMENTE NAO PODE SER O MINIMO GLOBAL DA CURVA. A TGLC dilui mais que as
    outras: o primario dela sai com 9.047 ppm ajustados contra 22.495 da QLP, e
    o minimo global caiu no SECUNDARIO - 11,019 d de diferenca, que sao 7,475
    ciclos, meio ciclo exato. A epoca resultante ficava 1.008 min fora e parecia
    discordancia entre reducoes quando era eclipse errado. A semente passa a ser
    o minimo previsto pela efemeride de referencia, e a profundidade ajustada
    fica gravada para que trocar primario por secundario seja visivel.
    """
    # `t14_h` entrou para a cadeia em lote: o trapezio de `sim.modelo` tinha a
    # duracao DESTE alvo (1,967 h) como default, e outra binaria com eclipse de
    # 1 h ou 4 h seria medida com a forma errada. Sem argumento, nada muda.
    t14_h = sim.T14_H if t14_h is None else float(t14_h)
    if semente is not None:
        k = np.round((np.mean(t) - semente) / P)
        t0c = float(semente + k * P)
    else:
        t0c = float(t[int(np.argmin(f))])
    grade = t0c + np.linspace(-span_min, span_min, n) / 1440.0
    y = f - 1.0
    chi2 = np.empty(n)
    for i, g in enumerate(grade):
        m = sim.modelo(t, g, P=P, t14_h=t14_h, prof=1.0) - 1.0   # perfil unitario
        den = float(m @ m)
        a = float(m @ y) / den if den > 0 else 0.0
        chi2[i] = float(np.sum((y - a * m) ** 2))
    var = float(np.var(y))
    c = chi2 / var
    i = int(np.argmin(c))
    t0 = float(grade[i])
    alvo = c[i] + 1.0
    esq = grade[:i + 1][c[:i + 1] <= alvo]
    dir_ = grade[i:][c[i:] <= alvo]
    sig = ((dir_.max() - esq.min()) / 2 * 1440.0
           if len(esq) and len(dir_) else np.nan)
    m = sim.modelo(t, t0, P=P, t14_h=t14_h, prof=1.0) - 1.0
    prof = float(m @ y) / float(m @ m)
    return t0, float(sig), prof


def rodar():
    """Duas passadas: a primeira ancora no s105 (2 min, sem ambiguidade), a
    segunda usa essa ancora como semente para todas as outras."""
    fontes = carregar()
    vistos, unicas = set(), []
    for s in fontes:
        chave = (s["setor"], s["reducao"])
        if chave in vistos:          # o mesmo FITS pode aparecer por dois caminhos
            continue
        vistos.add(chave)
        unicas.append(s)
    fontes = unicas
    ancora = [s for s in fontes if s["setor"] == 105]
    if not ancora:
        raise RuntimeError("FALHA [epocas]: sem s105 para ancorar a semente")
    t0_anc, _, _ = medir(ancora[0]["t"], ancora[0]["f"])
    linhas = []
    for s in fontes:
        t0, sig, prof = medir(s["t"], s["f"], semente=t0_anc)
        # epoca no ciclo mais proximo do centro do conjunto de dados
        linhas.append({k: s[k] for k in
                       ("setor", "reducao", "arquivo", "cadencia_min", "n")}
                      | {"t0_btjd": t0, "sigma_min": sig,
                         "prof_ajustada_ppm": prof * 1e6,
                         "t_medio": float(np.mean(s["t"])),
                         "baseline_d": float(s["t"].max() - s["t"].min())})
    d = pd.DataFrame(linhas).sort_values(["setor", "reducao"])
    d.to_parquet(OUT, index=False)
    return d


def oc(d, P=P_REF, t0_ref=None):
    """O-C contra a efemeride linear, com a epoca de referencia no s105."""
    if t0_ref is None:
        s105 = d[d["setor"] == 105]
        if not len(s105):
            raise RuntimeError("FALHA [epocas]: sem s105 para ancorar a efemeride")
        t0_ref = float(s105["t0_btjd"].iloc[0])
    d = d.copy()
    d["ciclo"] = np.round((d["t0_btjd"] - t0_ref) / P).astype(int)
    d["oc_min"] = (d["t0_btjd"] - (t0_ref + d["ciclo"] * P)) * 1440.0
    return d, t0_ref


# -----------------------------------------------------------------------------
# RESULTADO (4 aglomerados, E de -1896 a 0, 7,65 anos)
#
#   setor  reducoes  epoca BTJD     sigma
#      4        3    1424.80937   0,51 min  (espalhamento entre QLP/eleanor/TGLC)
#     31        1    2157.40971   0,30 min
#     97        1    3961.64447   0,30 min
#    105        1    4219.60278   0,30 min
#
# O espalhamento de 0,51 min entre as tres reducoes independentes do setor 4 bate
# com os 0,47 min que a simulacao previu ANTES da extracao. E a confirmacao de
# que a simulacao mediu a coisa certa.
#
#   LINEAR      residuos +1,30  -0,61  +0,10  +0,07 min   chi2 10,73 / 2 dof
#   PARABOLA    residuos +0,20  -0,13  +0,31  -0,24 min   chi2  2,04 / 1 dof
#   delta chi2 = 8,68 com 1 parametro a mais -> p = 0,003
#
#   Q = 1,368e-09 d/ciclo^2  ->  dP/dt = 6,78e-07 d/ano = 5,9e-02 s/ano
#
# O QUE ISTO ESTABELECE E O QUE NAO ESTABELECE.
#
# ESTABELECE: ha curvatura no O-C dentro do TESS sozinho, sem depender do
# SuperWASP. Antes das FFI havia 3 aglomerados e 3 parametros - a parabola
# passava por construcao. Agora sao 4 pontos e 1 grau de liberdade, e a parabola
# melhora de forma afirmavel.
#
# NAO ESTABELECE que a causa seja transferencia de massa. Uma senoide de tempo
# de luz por terceiro corpo precisa de 4+ parametros e com 4 pontos fica com 0
# dof - NAO TESTAVEL. A epoca do SuperWASP e o quinto ponto que falta, e ela e
# que decide entre parabola e senoide.
#
# O SETOR 4 CARREGA O SINAL SOZINHO. Sensibilidade medida deslocando so a epoca
# do s4:
#
#   -1,0 min -> delta chi2 2,36  p = 0,124   (o resultado some)
#   -0,5 min -> delta chi2 5,03  p = 0,025
#    0,0 min -> delta chi2 8,68  p = 0,003
#   +0,5 min -> delta chi2 13,3  p < 0,001
#
# Ou seja, o resultado exige que a epoca do s4 esteja certa dentro de ~1 min. O
# espalhamento entre reducoes e 0,51 min, o que e tranquilizador mas NAO fecha:
# as tres compartilham os pixels, e uma sistematica comum de FFI nao apareceria
# nesse espalhamento. Cadencia de 30 min, eclipse de 118 min, ~17 eclipses no
# setor - e o unico ponto do ajuste com essa fragilidade.


# -----------------------------------------------------------------------------
# ROBUSTEZ DO delta chi2 = 12,23
#
# OS QUATRO PONTOS SAO TODOS DO TESS. O SuperWASP NAO esta neste ajuste - ele
# continua sendo o quinto ponto que falta, e e ele que separaria parabola de
# senoide. O residuo de +2,16 min no ajuste linear e do SETOR 4, nao de 2006.
#
#   setor    ano     E      sigma
#       4  2018,84  -1896   0,65 min   (espalhamento entre 3 reducoes + correcao)
#      31  2020,84  -1399   0,30 min
#      97  2025,78   -175   0,30 min
#     105  2026,49      0   0,30 min
#
# A estrutura de alavanca e a que se espera: s97 e s105 estao a 0,71 ano um do
# outro e fixam o vertice juntos; a curvatura vem do trio s4/s31/(s97+s105).
#
# DESLOCANDO CADA PONTO EM 1 SIGMA, UM DE CADA VEZ:
#
#     s4   +1sig -> 18,18      s31  +1sig ->  9,32
#     s97  +1sig -> 10,36      s105 +1sig -> 14,95
#
#   Nenhum ponto sozinho derruba: o pior caso e 9,32, p = 0,0023.
#
# TESTE ADVERSARIAL - os quatro deslocados 1 sigma na combinacao que mais
# enfraquece (s4-, s31+, s97+, s105-):
#
#     1 sigma -> delta chi2 2,68   p = 0,10   o resultado dissolve
#     2 sigma -> delta chi2 0,05   p = 0,83
#
# E o resultado NAO depende da correcao de vies para existir: removendo os
# +0,9 min por completo, delta chi2 = 5,93, p = 0,015.
#
# LEITURA HONESTA: significativo, robusto a erro em qualquer ponto isolado, e
# nao robusto a uma conspiracao coerente de 1 sigma nos quatro. Isso e o estado
# normal de um ajuste com 4 pontos e 1 grau de liberdade, e e por isso que o
# quinto ponto - a epoca do SuperWASP - vale mais que refinar os quatro.
#
# O QUE FALTARIA PARA INCLUIR O SuperWASP. Ele e fotometria de solo de 2006-2008
# com abertura de dezenas de segundos de arco e banda larga sem filtro, e este
# alvo tem vizinha que entra inteira nessa abertura. A epoca dele precisaria do
# mesmo tratamento que o s4 recebeu - vies medido contra verdade conhecida, nao
# suposto zero -, e nao ha aqui dado com verdade conhecida naquela cadeia. Ate
# la, qualquer ajuste que o inclua e condicionado a uma acuracia nao verificada.
