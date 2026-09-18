# -*- coding: utf-8 -*-
"""Rodada sete, Passo H: o teste com TODOS os setores do cache para os 9 do D9
(generalizacao de caso_232634196_completo.py), no estado E da cadeia.

Para cada alvo do D9, com as epocas primarias TESS de jitter_primarios_26.parquet
(estimador do estado D, TESS-only: nao dependem do vies), E e criterios por
`jitter_primarios_26.classificar_setores` (P, T0 da escada do registro), as
temporadas SuperWASP dos `pontos` do registro (estado E) e o GP agrupado
Matern-3/2 (deriva_gp_26.json, C = diag + K so entre TESS, GLS):
  1. dP/dt so-TESS (quadratica no bloco 2019-2025) e dP/dt em tudo (TESS + SW);
     p_gof de cada; Delta chi2 cubica x quadratica (so-TESS e tudo) com p pelo
     teste F (F = Delta chi2 / (chi2_cub / nu_cub) ~ F(1, nu_cub): o Delta chi2
     cru superestima quando chi2/nu > 1, como em 232634196).
  2. O-C das epocas TESS agrupado por ano civil contra a quadratica-tudo
     (media, sd, n por ano).
  3. P_local(2024) - P_escada: linear (GP) aos setores 76-86 presentes (>= 5),
     contra o que a quadratica-tudo exige no centro do bloco 2024:
     P_quad(E) = c1 + 2 c0 E -> P_quad(E_2024) - P_escada, com sigma dos dois.
  V564 Dra (329246824): tambem contra o VarAstro (registro inteiro, 26 anos:
     dPdt_brno de comparacao_brno_26.parquet) e a janela.

CLASSIFICACAO (declarada antes de rodar):
  "estavel": |dP/dt_so-TESS - dP/dt_7| / hypot(s, s_7) <= 2 E |dP/dt_tudo -
     dP/dt_7| / hypot(s, s_7) <= 2 E cubica nao exigida (p_F >= 0,05 no so-TESS);
  "dependente da janela": caso contrario;
  "nao testavel": < 8 setores primarios ok no cache (a cubica precisaria de >= 4
     graus de liberdade e a efemeride local de >= 5 setores em 2024): V Gru
     (144194304, 4 setores: 1, 28, 68, 105 - os do sul) entra assim; os itens
     computaveis (so-TESS com 1 gl, tudo com 7 pontos) sao impressos para registro.
  Prioridade: 329246824, 392536812, 198408416; depois os demais.

EXPECTATIVA (escrita e commitada ANTES de rodar; informada pelo sigma_j linear /
quadratico por alvo do Passo 1 e pelo (A, tau_c) proprio do Passo B, que sao
medidas TESS-only ja feitas - e NAO cega para 232634196, que e o Passo F):
  232634196: dependente da janela (F: so-TESS -0,70 vs -0,215; cubica exigida
     Delta chi2 39; P_local - P_escada +0,46 +- 0,44 s contra -2,3 s). No estado
     E os numeros movem < 1% (dP/dt_7 -0,2167).
  424461577 (V527 Dra, LTTE 2,73 a, 6,6 min; sigma_j 3,6, A 5,5 min, tau_c 1 a):
     dependente da janela com certeza - o bloco TESS cobre 2 ciclos: cubica
     exigida (p_F < 1e-3), so-TESS longe de -0,036 (|z| > 3).
  392536812 (NSVS 2953494, +0,054 +- 0,003, P 0,59 d): a curvatura do desenho
     daria 7,8 min de ponta a ponta em 2019-2025 (sagita 1,9 min) - mas sigma_j
     LINEAR = 0 (IC95 < 0,4 min) nos 42 setores com barras 0,58: o bloco TESS e
     linear dentro das barras. Espero so-TESS |dP/dt| < 0,02, incompativel com
     +0,054 a > 2 sigma (70%); cubica nao exigida (70%) -> "dependente da janela"
     pela clausula de compatibilidade: o +0,054 e carregado pela alavanca de
     2004-08, nao pelo TESS.
  230386284 (StKM 1-1676, +0,0087 +- 0,0006, P 0,34 d, barras 0,1 min, sigma_j
     0,18 -> 0,00 na quadratica, A = 0): a sagita esperada de ~0,5 min e vista
     (a quadratica esgota o sigma_j). Espero so-TESS +0,006 a +0,011 +- 0,002-
     0,003, compativel; cubica nao exigida -> "estavel" (70%). E o alvo com o
     teste mais forte da rodada.
  329246824 (V564 Dra, -0,020 +- 0,0026, P 0,59 d, sigma_j 0,41 = 0,41, A 0,6,
     tau_c 0,5 a): a quadratica nao absorve nada do sigma_j -> a curvatura
     interna e menor que a do desenho; so-TESS entre -0,020 e 0,000 +- ~0,010,
     compativel a 2 sigma (70%); cubica exigida 50% (tau_c de meio ano) ->
     "dependente" 50%. Tudo (45): -0,015 a -0,020. VarAstro: -0,0066 +- 0,0011
     em 26 anos fica a 3-5 sigma do tudo e a < 2 sigma do so-TESS.
  198408416 (NSVS 2910034, +0,0103 +- 0,0011, P 0,36 d, sigma_j 0,93 -> 0,51,
     A 0,6, tau_c 0,47 a): a quadratica absorve metade do sigma_j - ha curvatura
     interna (sagita esperada 0,6 min). so-TESS compativel a 2 sigma (60%);
     cubica exigida 50% -> "estavel" 35%.
  377253090 (NSVS 3056525, -0,0057 +- 0,0011, sagita 0,3 min, sigma_j 0,34 =
     0,32, A 0,32, tau_c 0,09 a): teste fraco (so-TESS +- 0,008); compativel
     quase certo; cubica nao exigida 65% -> "estavel" 60%, com a ressalva de que
     "estavel" aqui e "nao contradito".
  390021728 ([GGM2006] 3034154, -0,0048 +- 0,0010, sagita 0,35 min, sigma_j 1,15
     -> 0,89, A 0,9, tau_c 0,09 a): teste fraco; cubica nao exigida 60% ->
     "estavel" 55%.
  144194304 (V Gru): nao testavel; so-TESS com 1 gl compativel (90%).
  Contagem esperada entre os 8 testaveis: estavel 2-4, dependente 4-6; dos tres
  prioritarios, 1 estavel (230386284 nao e prioritario; dos tres, 198408416 e o
  unico com chance real).
Desvio em qualquer direcao vai ao relatorio.

RESULTADO (2026-09-17, contra a expectativa acima; log em logs_E/31_caso_completo_d9.txt):
  8 testaveis (38-43 setores cada, 2019,58-2024,93; V Gru nao testavel com 4) -
  3 ESTAVEIS (392536812, 230386284, 377253090) e 5 DEPENDENTES DA JANELA
  (329246824, 198408416, 390021728, 232634196, 424461577); esperava 2-4 / 4-6.
  Nenhum setor > 2025,0 em nenhum dos 9 (o cache termina em s86).
  Por alvo (so-TESS | tudo | excursao do O-C por ano | P_local - P_escada
  contra o que a quadratica exige):
  - 230386284 +0,0088 +- 0,0110 (z +0,07) | +0,0081 +- 0,0013 (z -0,02) |
    excursao 0,11 min em 6 anos | +0,030 +- 0,090 s contra +0,020 +- 0,011:
    ESTAVEL, e o unico com confirmacao POSITIVA nos tres testes (como previsto).
  - 377253090 -0,0022 +- 0,0139 | -0,0075 +- 0,0016 | 0,66 min | -0,005 +- 0,115
    contra -0,029: ESTAVEL por teste fraco (nao contradito), como previsto.
  - 392536812 +0,0184 +- 0,0216 (z -1,60) | +0,0491 +- 0,0031 | 0,82 min |
    +0,087 +- 0,191 contra +0,095 +- 0,022: ESTAVEL - DESVIO: eu previa
    "dependente" por incompatibilidade, mas a barra so-TESS (0,022) e grande
    demais para excluir +0,053. O teste nao contradiz; tambem nao confirma.
  - 329246824 (V564 Dra) -0,0036 +- 0,0196 (z +0,86) | -0,0172 +- 0,0029 |
    1,49 min | +0,291 +- 0,167 contra -0,014 (+1,8 sigma): DEPENDENTE - mas
    SO pelo teste F da cubica (Delta chi2 2,7, p_chi2 0,099, F 11,4 porque
    chi2/nu = 0,30). Com o criterio de Delta chi2 seria ESTAVEL. VarAstro
    (26,3 a): -0,0066 +- 0,0011 fica a +0,15 sigma do so-TESS e a -3,4 sigma
    do tudo; na janela (-0,0130 +- 0,0042), +0,47 e -0,83 sigma.
  - 198408416 +0,0452 +- 0,0120 (z +2,94) | +0,0122 +- 0,0014 | 2,40 min, com
    minimo em 2021 (-1,31) e maximo em 2019 (+1,08) | +0,203 +- 0,099 contra
    +0,046: DEPENDENTE por incompatibilidade; a cubica e exigida no ajuste com
    tudo (Delta chi2 8,0, p_F 5e-5). E um alvo do NUCLEO com a mesma assinatura
    de 232634196, tres vezes menor.
  - 390021728 +0,0272 +- 0,0103 (z +3,15) | -0,0033 +- 0,0012 | 2,61 min |
    -0,029 +- 0,059 contra +0,005: DEPENDENTE por incompatibilidade; so-TESS
    com chi2/nu 1,98 e p_gof 0,000 (o unico alem de 232634196 e V527 com o
    bloco TESS mal descrito pelo GP agrupado). O sinal do dP/dt so-TESS e
    OPOSTO ao do desenho.
  - 232634196 -0,7022 +- 0,0473 (z -9,41) | -0,2940 +- 0,0185, chi2/nu 4,10 |
    8,47 min | +0,464 +- 0,441 contra -0,754 +- 0,061 (+2,7 sigma): DEPENDENTE,
    como no Passo F (os numeros do estado E movem < 0,5%).
  - 424461577 (V527 Dra) +0,0513 +- 0,0264 (z +3,31) | -0,0048 +- 0,0035
    (z +6,73 do -0,0370 do desenho) | 8,22 min | +1,708 +- 0,251 contra +0,227
    (+5,9 sigma): DEPENDENTE. DESVIO na causa: eu previa cubica exigida
    (p_F < 1e-3) e deu p_F 0,47 - com chi2/nu 6,4 o denominador do F absorve o
    ciclo; o LTTE de 2,73 a aparece como dispersao, nao como termo cubico. O
    resultado importante: com os 40 setores a curvatura de V527 Dra CAI de
    -0,0370 +- 0,0033 para -0,0048 +- 0,0035, compativel com zero - o falso
    positivo conhecido (O-C publicado sem termo secular) e corrigido pelo bloco
    TESS completo.
  - 144194304 (V Gru): NAO TESTAVEL (4 setores, 2018,6-2026,5); para registro,
    so-TESS +0,0283 +- 0,0127 (z +0,84) e tudo +0,0157 +- 0,0020 (z -0,76).
  Leitura: dos 4 do nucleo, 1 estavel (392536812, por teste fraco), 3
  dependentes; dos 8 testaveis, so 230386284 tem confirmacao positiva. A
  dependencia da janela nao e propriedade de 232634196: e da maioria do D9.

POS-HOC (--por-alvo; sugestao do revisor interno, sem expectativa previa): o mesmo
teste com o (A, tau_c) PROPRIO de cada alvo em vez do agrupado. Duas classes mudam:
  - 392536812 ESTAVEL -> DEPENDENTE: com A proprio = 0 a barra so-TESS cai de 0,0216
    para 0,0074 e o bloco TESS passa a CONTRADIZER o desenho (+0,0209 contra +0,0532:
    -4,05 sigma). O "estavel" do agrupado era barra grande, como suspeitado.
  - 424461577 DEPENDENTE -> ESTAVEL: com A proprio = 5,54 min (ajustado no proprio
    LTTE) a barra vai a 0,185 e nada e contradito. E circular - a amplitude do GP foi
    estimada do sinal que se quer testar - e NAO deve ser lido como estabilidade.
  Os outros seis nao mudam de classe: 230386284 fica estavel e MUITO mais apertado
  (+0,0091 +- 0,0008 contra +0,0081 +- 0,0006 do desenho: +0,97 sigma - a confirmacao
  positiva mais forte do conjunto); 198408416 +0,0468 +- 0,0107 (+3,45 sigma);
  390021728 +0,0286 +- 0,0088 (+3,83 sigma); 232634196 -0,6485 +- 0,1328 (-3,21) com
  cubica ainda exigida (p_F 5e-4); 329246824 +0,0001 +- 0,0184 com cubica exigida
  (p_F 4e-3); 377253090 -0,0050 +- 0,0049 (+0,27). V Gru, nao testavel pela contagem,
  daria z +2,28 com A proprio = 0.
  Leitura: a contagem 3 estaveis / 5 dependentes e a mesma nos dois, mas QUEM esta em
  cada lado muda; o unico alvo estavel nos dois com barra apertada e 230386284.

PASSO I (rodada oito, --sem-desenho): o mesmo so-TESS com as 4 epocas do DESENHO
(as que entram nas Tabelas 3-5) REMOVIDAS - o que sobra do bloco TESS e um
conjunto independente daquele com que se compara. Criterio unico: o ajuste
so-TESS sob o GP AGRUPADO. Regra de classe declarada antes de rodar
(`classe_I`): "estavel" = |z| < 2 E cubica nao exigida pelo Delta chi2
(p_chi2 > 0,05; NAO pelo teste F); "contraditado" = |z| > 3; "ambiguo" = o
resto. As classes saem lado a lado, com e sem as epocas do desenho.

EXPECTATIVA DO PASSO I (escrita e commitada ANTES de rodar):
  Com desenho, a regra nova aplicada aos numeros do Passo H da: estavel 4
  (329246824 z +0,86 p_cub 0,099; 392536812 -1,60 / 0,896; 230386284 +0,07 /
  0,910; 377253090 +0,29 / 0,740), ambiguo 1 (198408416 +2,94), contraditado 3
  (390021728 +3,15; 232634196 -9,41; 424461577 +3,31). Note que 329246824 muda
  de lado em relacao ao Passo H so pela troca do F pelo Delta chi2.
  Sem desenho: 4 de ~40 epocas saem, e sao as das PONTAS do bloco (a escada
  escolhe aglomerados separados), entao espero barra so-TESS x1,05-1,35
  (mediana ~1,15) e |dz| tipico < 0,6, maximo ~1,5 em 232634196 (onde as 4 do
  desenho ancoram a parabola). Contagem esperada: estavel 4 +- 1, ambiguo 1-2,
  contraditado 2-3; os que podem trocar sao os de borda - 198408416 (2,94:
  ambiguo -> contraditado em ~40%), 424461577 (3,31: contraditado -> ambiguo em
  ~50%, porque a barra cresce), 390021728 (3,15: -> ambiguo em ~35%).
  230386284, 377253090 e 232634196 NAO devem mudar de classe (margem grande).
  Se algum "estavel" virar "contraditado" (ou vice-versa) sem ser de borda, e
  sinal de que as epocas do desenho carregavam o resultado - vai ao relatorio.

RESULTADO DO PASSO I (2026-09-17, contra a expectativa acima; log
logs_E/33_passo_I.txt): as classes sao IDENTICAS com e sem as epocas do desenho
nos 8 testaveis - estavel 4 (329246824, 392536812, 230386284, 377253090),
ambiguo 1 (198408416), contraditado 3 (390021728, 232634196, 424461577);
nenhum muda. A barra so-TESS quase nao se mexe (x1,00 mediana, max x1,01, contra
os x1,05-1,35 que eu esperava) e |dz| mediana 0,02, max 0,09 (esperava < 0,6 e
ate ~1,5 em 232634196). DESVIO, e na direcao boa: as 4 epocas do desenho pesam
~1/10 do bloco e nao estao nas pontas em termos de peso - o resultado so-TESS e
independente delas. 144194304 (V Gru) nao entra: as 4 epocas do cache SAO as do
desenho. A regra nova (Delta chi2 em vez do teste F) move 329246824 de
"dependente da janela" (Passo H) para "estavel": e o alvo cuja classe depende do
criterio de cubica, como ja registrado.
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
import gls_vies_26 as G  # noqa: E402
import jitter_primarios_26 as J  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
BASE = config.DATA / "orquestra" / "oc_lote"
S_POR_ANO = 365.25 * 86400.0
D9 = [329246824, 392536812, 198408416, 230386284, 377253090, 390021728, 232634196, 424461577, 144194304]
LOCAIS = list(range(76, 87))
N_MIN_TESTAVEL = 8
N_MIN_LOCAL = 5
P_F = 0.05
P_CUB = 0.05   # Passo I: limiar do Delta chi2 da cubica
Z_COMPAT = 2.0
V564 = 329246824
EXPECTATIVA = "ver git log: commit da expectativa (rodada sete, Passo H)"
EXPECTATIVA_I = "ver git log: commit da expectativa (rodada oito, Passo I)"


def ano(bt):
    return 2000.0 + (np.asarray(bt, float) + 2457000.0 - 2451544.5) / 365.25


def ajustar(pp, grau, nome, A, tc, P, E0):
    """GLS de t(E - E0) com C = diag + K entre TESS; coef na ordem do vander (grau mais alto primeiro)."""
    E = pp.E.values.astype(float) - E0; t = pp.t0.values.astype(float)
    C, it = GP.cov_desenho(pp, nome, A, tc)
    coef, r, chi2, cov = G.ajustar_gls(E, t, C, grau)
    dof = len(E) - grau - 1
    out = {"n": int(len(E)), "grau": grau, "E0": float(E0), "chi2": float(chi2), "dof": int(dof), "chi2r": float(chi2 / dof) if dof > 0 else np.nan,
           "p_gof": float(stats.chi2.sf(chi2, dof)) if dof > 0 else np.nan, "coef": coef.tolist(), "cov": cov.tolist(), "res_min": (r * 1440.0).tolist()}
    if grau == 2:
        out["dPdt"] = 2 * coef[0] / P * S_POR_ANO; out["sdPdt"] = 2 * np.sqrt(cov[0, 0]) / P * S_POR_ANO
    return out


def teste_F(q, c):
    """Cubica x quadratica: Delta chi2 cru (p qui2 com 1 gl) e teste F (Delta chi2 / (chi2_cub / nu_cub) ~ F(1, nu_cub))."""
    d = q["chi2"] - c["chi2"]
    if c["dof"] < 1:
        return {"delta_chi2": float(d), "p_chi2": np.nan, "F": np.nan, "p_F": np.nan, "aplicavel": False}
    F = d / (c["chi2"] / c["dof"])
    return {"delta_chi2": float(d), "p_chi2": float(stats.chi2.sf(d, 1)), "F": float(F), "p_F": float(stats.f.sf(F, 1, c["dof"])), "aplicavel": True}


def periodo_da_quadratica(q, E, P):
    """P_quad(E) - P_escada em segundos, com sigma: P(E) = c1 + 2 c0 (E - E0)."""
    c = np.asarray(q["coef"]); cv = np.asarray(q["cov"]); x = E - q["E0"]
    g = np.array([2 * x, 1.0, 0.0])
    return float((g @ c - P) * 86400.0), float(np.sqrt(g @ cv @ g) * 86400.0)


def classe_I(z, p_cub_chi2, testavel):
    """Regra do Passo I (rodada oito), declarada antes de rodar: o unico criterio e o ajuste
    so-TESS sob o GP AGRUPADO, e a cubica e julgada pelo Delta chi2 (p_chi2), nao pelo teste F
    (cujo denominador e ajustado nos mesmos residuos e dispara quando as barras cobrem demais)."""
    if not testavel or not np.isfinite(z):
        return "nao testavel"
    if abs(z) > 3.0:
        return "contraditado"
    if abs(z) < 2.0 and (not np.isfinite(p_cub_chi2) or p_cub_chi2 > P_CUB):
        return "estavel"
    return "ambiguo"


def um_alvo(tic, S, A, tc, brno, sem_desenho=False):
    j = json.loads((BASE / f"oc__{tic}.json").read_text(encoding="utf-8"))
    P, sP, T0, E_ult = J.escada_do_registro(j)
    d = J.classificar_setores(S[(S.tic == tic) & (S.minimo == "primario")], P, sP, T0, E_ult, 0.0)
    fora = d[~d.ok]
    d = d[d.ok]
    tess = pd.DataFrame({"fonte": ["TESS s%d" % s for s in d.setor], "t0": d.t0_btjd.values, "sig_min": d.sig_min.values, "E": d.E.values,
                         "setor": d.setor.values, "ano": ano(d.t0_btjd.values)}).sort_values("t0").reset_index(drop=True)
    pts = pd.DataFrame(j["pontos"]).sort_values("t0").reset_index(drop=True)
    sw = pts[pts.fonte.str.startswith("SuperWASP")][["fonte", "t0", "sig_min", "E"]].copy(); sw["setor"] = -1; sw["ano"] = ano(sw.t0.values)
    # controle: as epocas TESS do desenho reproduzem os pontos do registro (mesmo estimador, mesmo E)
    n_des = 0
    for r in pts[pts.fonte.str.startswith("TESS")].itertuples():
        m = tess[tess.fonte == r.fonte]
        assert len(m) == 1 and abs(float(m.t0.iloc[0]) - r.t0) * 1440 < 0.01 and int(m.E.iloc[0]) == int(r.E), (tic, r.fonte)
        n_des += 1
    fontes_desenho = set(pts[pts.fonte.str.startswith("TESS")].fonte)
    if sem_desenho:   # Passo I: fora as epocas que entram nas Tabelas 3-5 - o resto do bloco TESS e independente delas
        tess = tess[~tess.fonte.isin(fontes_desenho)].reset_index(drop=True)
    ref = j["parabola"]; d7, s7 = float(ref["dPdt_s_por_ano"]), float(ref["sdPdt_s_por_ano"])
    E0 = float(np.mean(tess.E.values))
    tudo = pd.concat([sw, tess], ignore_index=True).sort_values("t0").reset_index(drop=True)
    res = {"tic": tic, "nome": j.get("nome", ""), "P_escada_d": P, "n_tess": int(len(tess)), "n_tess_desenho": n_des, "n_fora_criterio": int(len(fora)),
           "motivos_fora": fora.motivo.value_counts().to_dict() if len(fora) else {}, "setores": tess.setor.tolist(),
           "ano_min": float(tess.ano.min()), "ano_max": float(tess.ano.max()), "n_sw": int(len(sw)), "dPdt_7": d7, "s_7": s7,
           "testavel": bool(len(tess) >= N_MIN_TESTAVEL)}
    # 1. so-TESS e tudo: linear, quadratica, cubica
    fits = {}
    for rot, pp in (("so_tess", tess), ("tudo", tudo)):
        f = {g: ajustar(pp, g, "matern32", A, tc, P, E0) for g in (1, 2, 3) if len(pp) - g - 1 >= 0}
        fits[rot] = f
        q = f[2]
        res[f"{rot}_dPdt"] = q["dPdt"]; res[f"{rot}_s"] = q["sdPdt"]; res[f"{rot}_chi2"] = q["chi2"]; res[f"{rot}_dof"] = q["dof"]
        res[f"{rot}_chi2r"] = q["chi2r"]; res[f"{rot}_p_gof"] = q["p_gof"]
        res[f"{rot}_delta_chi2_lin_quad"] = f[1]["chi2"] - q["chi2"] if 1 in f else np.nan
        res[f"{rot}_z_vs_7"] = (q["dPdt"] - d7) / np.hypot(q["sdPdt"], s7)
        tf = teste_F(q, f[3]) if 3 in f else {"delta_chi2": np.nan, "p_chi2": np.nan, "F": np.nan, "p_F": np.nan, "aplicavel": False}
        res[f"{rot}_cub_delta_chi2"] = tf["delta_chi2"]; res[f"{rot}_cub_p_chi2"] = tf["p_chi2"]; res[f"{rot}_cub_F"] = tf["F"]; res[f"{rot}_cub_p_F"] = tf["p_F"]; res[f"{rot}_cub_aplicavel"] = tf["aplicavel"]
        if 3 in f and f[3]["dof"] >= 1:
            # curvatura local da cubica no centro do bloco (E = E0): d2t/dE2 = 2 c2 (c3 anula em E0)
            c3 = np.asarray(f[3]["coef"]); cv3 = np.asarray(f[3]["cov"]); g = np.array([0.0, 2.0, 0.0, 0.0])
            res[f"{rot}_cub_dPdt_centro"] = float(g @ c3) / P * S_POR_ANO; res[f"{rot}_cub_s_centro"] = float(np.sqrt(g @ cv3 @ g)) / P * S_POR_ANO
    # 2. O-C das epocas TESS por ano civil contra a quadratica-tudo
    qa = fits["tudo"][2]; ca = np.asarray(qa["coef"])
    oc = (tess.t0.values - np.polyval(ca, tess.E.values - E0)) * 1440.0
    anos = np.floor(tess.ano.values).astype(int)
    por_ano = {int(a): {"n": int((anos == a).sum()), "media_min": float(oc[anos == a].mean()), "sd_min": float(oc[anos == a].std(ddof=0)),
                        "barra_mediana_min": float(np.median(tess.sig_min.values[anos == a]))} for a in sorted(set(anos))}
    res["oc_por_ano"] = por_ano
    res["oc_excursao_min"] = float(max(v["media_min"] for v in por_ano.values()) - min(v["media_min"] for v in por_ano.values()))
    # 3. efemeride local 2024 x o que a quadratica-tudo exige
    local = tess[tess.setor.isin(LOCAIS)].reset_index(drop=True)
    res["n_local"] = int(len(local)); res["setores_local"] = local.setor.tolist()
    if len(local) >= N_MIN_LOCAL:
        ll = ajustar(local, 1, "matern32", A, tc, P, E0)
        P_loc = ll["coef"][0]; sP_loc = np.sqrt(ll["cov"][0][0])
        E24 = float(np.mean(local.E.values))
        dq, sq = periodo_da_quadratica(qa, E24, P)
        res.update({"P_local_menos_escada_s": float((P_loc - P) * 86400.0), "s_P_local_s": float(sP_loc * 86400.0), "P_quad2024_menos_escada_s": dq, "s_P_quad2024_s": sq,
                    "z_local_vs_quad": float(((P_loc - P) * 86400.0 - dq) / np.hypot(sP_loc * 86400.0, sq)), "local_chi2": ll["chi2"], "local_dof": ll["dof"], "local_p_gof": ll["p_gof"]})
    else:
        res.update({"P_local_menos_escada_s": np.nan, "s_P_local_s": np.nan, "P_quad2024_menos_escada_s": np.nan, "s_P_quad2024_s": np.nan, "z_local_vs_quad": np.nan})
    # classificacao declarada
    if not res["testavel"]:
        res["classe"] = "nao testavel"
    else:
        compat = abs(res["so_tess_z_vs_7"]) <= Z_COMPAT and abs(res["tudo_z_vs_7"]) <= Z_COMPAT
        cub = res["so_tess_cub_aplicavel"] and res["so_tess_cub_p_F"] < P_F
        res["classe"] = "estavel" if (compat and not cub) else "dependente da janela"
        res["compat_2sigma"] = bool(compat); res["cubica_exigida"] = bool(cub)
    res["sem_desenho"] = bool(sem_desenho); res["n_desenho_removidas"] = int(n_des if sem_desenho else 0)
    res["classe_I"] = classe_I(res["so_tess_z_vs_7"], res["so_tess_cub_p_chi2"], res["testavel"])
    if tic == V564 and brno is not None:
        b = brno.set_index("tic").loc[V564]
        for rot in ("so_tess", "tudo"):
            res[f"{rot}_z_vs_varastro_26a"] = float((res[f"{rot}_dPdt"] - b.dPdt_brno) / np.hypot(res[f"{rot}_s"], b.s_dPdt_brno))
            res[f"{rot}_z_vs_varastro_janela"] = float((res[f"{rot}_dPdt"] - b.dPdt_brno_janela) / np.hypot(res[f"{rot}_s"], b.s_dPdt_brno_janela))
        res["varastro_26a"] = float(b.dPdt_brno); res["varastro_26a_s"] = float(b.s_dPdt_brno); res["varastro_janela"] = float(b.dPdt_brno_janela); res["varastro_janela_s"] = float(b.s_dPdt_brno_janela)
        res["varastro_alavanca_anos"] = float(b.alavanca_anos)
    res["fits"] = {k: {str(g): v for g, v in f.items()} for k, f in fits.items()}
    return res, tudo


if __name__ == "__main__":
    assert abs(co.VIES_CADEIA_MIN - (-1.26)) < 1e-9, "vies comum da cadeia trocado"   # o estado F manteve o vies comum e acrescentou a correcao de forma por alvo (Apendice D)
    S = pd.read_parquet(BASE / "jitter_primarios_26.parquet")
    hip = json.loads((BASE / "deriva_gp_26.json").read_text(encoding="utf-8"))["hiper_agrupados"]["matern32"]
    A, tc = hip["A_min"], hip["tau_c_d"]
    t26 = pd.read_parquet(BASE / "tabela_26.parquet").set_index("TIC")
    d9_tab = sorted(int(t) for t in t26.index if t26.loc[t, "curvatura"] and not t26.loc[t, "inadequada"])
    assert set(d9_tab) == set(D9), (d9_tab, "o D9 do estado E nao e o esperado")
    brno = pd.read_parquet(BASE / "comparacao_brno_26.parquet")
    print(f"GP agrupado Matern-3/2: A {A:.3f} min, tau_c {tc:.1f} d; cadeia estado E (vies {co.VIES_CADEIA_MIN:+.2f} +- {co.VIES_CADEIA_SIG:.2f})\n")
    linhas, pontos = [], []
    for tic in D9:
        r, tudo = um_alvo(tic, S, A, tc, brno)
        r["nome"] = str(t26.loc[tic, "nome"])
        linhas.append(r); pontos.append(tudo.assign(tic=tic))
        print(f"== TIC {tic} {r['nome']} (P {r['P_escada_d']:.4f} d): {r['n_tess']} setores ok no cache ({r['n_tess_desenho']} no desenho; {r['n_fora_criterio']} fora do criterio {r['motivos_fora']}), "
              f"{r['ano_min']:.2f}-{r['ano_max']:.2f}; 7 epocas: {r['dPdt_7']:+.4f} +- {r['s_7']:.4f} -> {'TESTAVEL' if r['testavel'] else 'NAO TESTAVEL'}")
        for rot in ("so_tess", "tudo"):
            print(f"   {rot:8s}: dP/dt {r[rot + '_dPdt']:+.4f} +- {r[rot + '_s']:.4f} (z vs 7 ep {r[rot + '_z_vs_7']:+.2f}); chi2 {r[rot + '_chi2']:.1f}/{r[rot + '_dof']} = {r[rot + '_chi2r']:.2f}, p_gof {r[rot + '_p_gof']:.3f}; "
                  f"Delta chi2 lin-quad {r[rot + '_delta_chi2_lin_quad']:.1f}; cubica: Delta chi2 {r[rot + '_cub_delta_chi2']:.1f}, p_chi2 {r[rot + '_cub_p_chi2']:.2e}, F {r[rot + '_cub_F']:.1f}, p_F {r[rot + '_cub_p_F']:.2e}"
                  + (f"; curvatura local no centro {r[rot + '_cub_dPdt_centro']:+.4f} +- {r[rot + '_cub_s_centro']:.4f}" if rot + "_cub_dPdt_centro" in r else ""))
        print("   O-C TESS por ano vs quadratica-tudo (min): " + "; ".join(f"{a}: {v['media_min']:+.2f} (sd {v['sd_min']:.2f}, n {v['n']}, barra {v['barra_mediana_min']:.2f})" for a, v in r["oc_por_ano"].items())
              + f" | excursao {r['oc_excursao_min']:.2f} min")
        if np.isfinite(r["P_local_menos_escada_s"]):
            print(f"   local 2024 ({r['n_local']} setores {r['setores_local'][0]}-{r['setores_local'][-1]}): P_local - P_escada {r['P_local_menos_escada_s']:+.3f} +- {r['s_P_local_s']:.3f} s; "
                  f"quadratica-tudo exige {r['P_quad2024_menos_escada_s']:+.3f} +- {r['s_P_quad2024_s']:.3f} s; diferenca {r['z_local_vs_quad']:+.1f} sigma")
        else:
            print(f"   local 2024: nao aplicavel ({r['n_local']} setores em 76-86)")
        if tic == V564:
            print(f"   VarAstro: registro inteiro ({r['varastro_alavanca_anos']:.1f} a) {r['varastro_26a']:+.4f} +- {r['varastro_26a_s']:.4f}: z so-TESS {r['so_tess_z_vs_varastro_26a']:+.2f}, tudo {r['tudo_z_vs_varastro_26a']:+.2f}; "
                  f"janela {r['varastro_janela']:+.4f} +- {r['varastro_janela_s']:.4f}: z so-TESS {r['so_tess_z_vs_varastro_janela']:+.2f}, tudo {r['tudo_z_vs_varastro_janela']:+.2f}")
        print(f"   -> {r['classe'].upper()}" + (f" (compativel a 2 sigma: {r['compat_2sigma']}; cubica exigida: {r['cubica_exigida']})" if "compat_2sigma" in r else "") + "\n")
    if "--sem-desenho" in sys.argv:
        print("== PASSO I: so-TESS SEM as epocas do desenho (Tabelas 3-5), GP agrupado como unico criterio\n")
        lado = []
        for tic in D9:
            com = next(x for x in linhas if x["tic"] == tic)
            if not com["testavel"]:   # V Gru: as 4 epocas do cache SAO as do desenho; sem elas nao sobra nada
                print(f"  {tic} {com['nome']:<18s} NAO TESTAVEL: {com['n_tess']} epocas no cache, todas do desenho - sem ajuste possivel")
                continue
            sem, _ = um_alvo(tic, S, A, tc, brno, sem_desenho=True)
            sem["nome"] = com["nome"]
            lado.append({"tic": tic, "nome": com["nome"], "testavel": com["testavel"], "n_tess_com": com["n_tess"], "n_tess_sem": sem["n_tess"], "n_removidas": com["n_tess_desenho"],
                         "dPdt_7": com["dPdt_7"], "s_7": com["s_7"],
                         "dPdt_com": com["so_tess_dPdt"], "s_com": com["so_tess_s"], "z_com": com["so_tess_z_vs_7"], "cub_p_chi2_com": com["so_tess_cub_p_chi2"], "chi2r_com": com["so_tess_chi2r"], "classe_I_com": com["classe_I"],
                         "dPdt_sem": sem["so_tess_dPdt"], "s_sem": sem["so_tess_s"], "z_sem": sem["so_tess_z_vs_7"], "cub_p_chi2_sem": sem["so_tess_cub_p_chi2"], "chi2r_sem": sem["so_tess_chi2r"], "classe_I_sem": sem["classe_I"],
                         "muda": com["classe_I"] != sem["classe_I"], "razao_sigma": sem["so_tess_s"] / com["so_tess_s"], "dz": sem["so_tess_z_vs_7"] - com["so_tess_z_vs_7"]})
            L = lado[-1]
            print(f"  {tic} {com['nome']:<18s} 7 epocas {L['dPdt_7']:+.4f} +- {L['s_7']:.4f} | COM desenho ({L['n_tess_com']:2d} ep) {L['dPdt_com']:+.4f} +- {L['s_com']:.4f} z {L['z_com']:+.2f} p_cub {L['cub_p_chi2_com']:.3f} -> {L['classe_I_com']:<12s}"
                  f" | SEM ({L['n_tess_sem']:2d} ep) {L['dPdt_sem']:+.4f} +- {L['s_sem']:.4f} z {L['z_sem']:+.2f} p_cub {L['cub_p_chi2_sem']:.3f} -> {L['classe_I_sem']:<12s}" + ("  <- MUDA" if L["muda"] else ""))
        LD = pd.DataFrame(lado); LD.to_parquet(BASE / "caso_completo_d9_sem_desenho.parquet", index=False)
        t = LD[LD.testavel]
        cont = lambda col: {c: int((t[col] == c).sum()) for c in ("estavel", "ambiguo", "contraditado")}  # noqa: E731
        print(f"\n  contagem nos {len(t)} testaveis: COM desenho {cont('classe_I_com')} | SEM desenho {cont('classe_I_sem')}; mudam {sorted(int(x) for x in t[t.muda].tic)}")
        print(f"  barra so-TESS sem/com: mediana x{t.razao_sigma.median():.2f} (min {t.razao_sigma.min():.2f}, max {t.razao_sigma.max():.2f}); |dz| mediana {t.dz.abs().median():.2f}, max {t.dz.abs().max():.2f}")
        jj = json.loads((BASE / "caso_completo_d9.json").read_text(encoding="utf-8"))
        jj["passo_I_sem_desenho"] = {"expectativa": EXPECTATIVA_I, "regra": "estavel: |z| < 2 e p_chi2 da cubica > 0,05; contraditado: |z| > 3; ambiguo: o resto", "lado_a_lado": lado}
        (BASE / "caso_completo_d9.json").write_text(json.dumps(jj, indent=1, default=float), encoding="utf-8")
        print(f"  -> {BASE / 'caso_completo_d9_sem_desenho.parquet'}, caso_completo_d9.json")
        sys.exit()
    # POS-HOC (sugestao do revisor interno, sem expectativa previa): o GP AGRUPADO (A 0,835 min) cobre demais os alvos
    # quietos - chi2/nu de 0,05 (230386284) a 0,58 - e de menos os barulhentos (chi2/nu 6,4 em V527). Isto infla a barra
    # so-TESS (que decide a compatibilidade) e deflaciona o F da cubica. Aqui o mesmo teste com o (A, tau_c) PROPRIO de
    # cada alvo (deriva_gp_26_hiper.parquet, MV no bloco estendido) - a versao principiada de "reescalar por sqrt(chi2/nu)".
    if "--por-alvo" in sys.argv:
        H = pd.read_parquet(BASE / "deriva_gp_26_hiper.parquet")
        H = H[H.kernel == "matern32"].set_index("tic")
        ph = []
        for tic in D9:
            h = H.loc[tic]
            r2, _ = um_alvo(tic, S, float(h.A_min), float(h.tau_c_d), brno)
            r2["A_proprio"] = float(h.A_min); r2["tau_c_proprio_d"] = float(h.tau_c_d); r2["tau_c_na_borda"] = bool(h.tau_c_na_borda)
            base = next(x for x in linhas if x["tic"] == tic)
            r2["classe_agrupado"] = base["classe"]; r2["muda"] = r2["classe"] != base["classe"]
            ph.append(r2)
            print(f"   {tic}: A proprio {h.A_min:.2f} min tau_c {h.tau_c_d:.0f} d | so-TESS {r2['so_tess_dPdt']:+.4f} +- {r2['so_tess_s']:.4f} (z {r2['so_tess_z_vs_7']:+.2f}; agrupado {base['so_tess_dPdt']:+.4f} +- {base['so_tess_s']:.4f}, z {base['so_tess_z_vs_7']:+.2f}) | "
                  f"chi2/nu {r2['so_tess_chi2r']:.2f} | cubica p_F {r2['so_tess_cub_p_F']:.2e} | {r2['classe']}" + (" <- MUDA" if r2["muda"] else ""))
        PH = pd.DataFrame([{k: v for k, v in r.items() if not isinstance(v, (dict, list))} for r in ph])
        PH.to_parquet(BASE / "caso_completo_d9_por_alvo.parquet", index=False)
        jj = json.loads((BASE / "caso_completo_d9.json").read_text(encoding="utf-8"))
        jj["pos_hoc_por_alvo"] = {"nota": "sugestao do revisor interno apos o Passo H; sem expectativa previa", "alvos": ph}
        (BASE / "caso_completo_d9.json").write_text(json.dumps(jj, indent=1, default=float), encoding="utf-8")
        print(f"   POS-HOC: mudam de classe {sorted(int(r['tic']) for r in ph if r['muda'])}; estaveis {sum(r['classe'] == 'estavel' for r in ph)} de {sum(r['testavel'] for r in ph)} testaveis")
        sys.exit()
    R = pd.DataFrame([{k: v for k, v in r.items() if not isinstance(v, (dict, list))} for r in linhas])
    R.to_parquet(BASE / "caso_completo_d9.parquet", index=False)
    pd.concat(pontos, ignore_index=True).to_parquet(BASE / "caso_completo_d9_pontos.parquet", index=False)
    (BASE / "caso_completo_d9.json").write_text(json.dumps({"expectativa": EXPECTATIVA, "gp": {"A_min": A, "tau_c_d": tc}, "vies_cadeia": [co.VIES_CADEIA_MIN, co.VIES_CADEIA_SIG], "alvos": linhas}, indent=1, default=float), encoding="utf-8")
    print("== resumo: " + "; ".join(f"{r.tic} {r.classe}" for r in R.itertuples()))
    print(f"   testaveis {int(R.testavel.sum())}: estavel {int((R.classe == 'estavel').sum())}, dependente {int((R.classe == 'dependente da janela').sum())}; nao testavel {int((~R.testavel).sum())}")
    print(R[["tic", "nome", "n_tess", "dPdt_7", "so_tess_dPdt", "so_tess_s", "so_tess_z_vs_7", "so_tess_cub_p_F", "tudo_dPdt", "tudo_s", "tudo_z_vs_7", "tudo_chi2r", "oc_excursao_min", "P_local_menos_escada_s", "P_quad2024_menos_escada_s", "z_local_vs_quad", "classe"]].round(4).to_string(index=False))
    print(f"  -> {BASE / 'caso_completo_d9.parquet'}, _pontos.parquet, .json")
