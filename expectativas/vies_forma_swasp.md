# Expectativa: `src/vies_forma_swasp.py` (Passo 4)

Escrita e commitada ANTES da primeira execucao. Nada rodado. `superwasp.py`
NAO muda: o script chama `superwasp.carregar` (clip assimetrico, HJD_UTC ->
BJD_TDB, normalizacao) e `superwasp.epoca_por_subconjunto` (trapezio com
profundidade livre, grade refinada, temporadas por lacuna > 30 d, gate de
cobertura) exatamente como `run_hj_calib.py` e `oc_lote.medir_alvo_lote`.

## O que o script faz

1. Curvas reais: os 4 calibradores (hj_calibradores.parquet: WASP-18 b, -4 b,
   -19 b, -6 b) e os 26 alvos (alvos_oc_88.parquet: sourceid, ra, dec, P, t14_h,
   depth_frac, morph). Cada CSV do cache passa por `carregar` -> t (BJD_TDB), f.
2. Controle embutido (a curva REAL pelo mesmo caminho): os 4 HJ reproduzem o
   O-C de hj_calibracao.parquet (-2,85, -2,55, +0,85, -0,64 min) a < 0,05 min;
   os 26 reproduzem as epocas por temporada dos `pontos` a < 0,05 min.
3. Residuos reais: r_i = f_i - curva media em fase (mediana em 200 bins de fase
   do proprio alvo, P da escada/ efemeride do HJ) - tira a forma real, deixa o
   ruido com a correlacao real. Bootstrap em blocos de NOITE (lacunas > 0,5 d):
   para cada noite real, o bloco de residuos de uma noite sorteada (com
   reposicao; encadeada se mais curta). Instantes, cadencia e temporadas reais.
4. Formas injetadas, todas SIMETRICAS em torno do minimo (o vies procedural que
   se testa e o da forma contra o trapezio + a amostragem real):
   (i) transito Mandel-Agol (`batman` 2.5.1, LD quadratico u1 0,45 u2 0,20
       declarados) com rp/rs, a/R*, i publicados (hj_transito_params.parquet:
       Cortes-Zuleta 2020 para WASP-18/19, Basturk 2025 para WASP-4, Kokori 2023
       para WASP-6) nas curvas dos 4 HJ;
   (ii) eclipse destacado profundo/curto: ocultacao com escurecimento de bordo
        (batman, i = 90, rp/rs = sqrt(depth_frac), a/R* = (1 + k)/sin(pi T14/P))
        com a profundidade e a duracao do catalogo, nas curvas dos 26;
   (iii) contato/W UMa: elipsoidal 1 - A_ell (1 - cos 4 pi phi)/2 com A_ell =
        0,5 x profundidade, mais eclipses primario (depth_frac) e secundario
        (0,8 x) de perfil cosseno de duracao T14, nas curvas dos 26.
   Verdade: efemeride linear (P da escada, T0 = epoca global real medida);
   t0_previsto dado ao pipeline = a verdade (como a propagacao TESS na cadeia,
   cujo erro de minutos e irrelevante na janela de P/2).
5. 500 realizacoes por curva e forma (semente fixa); vies = epoca ajustada -
   verdadeira (global e por temporada), barra formal, dispersao entre
   temporadas, profundidade, status. Paralelo (10 processos, ~1-1,5 h).
6. Resumo: media do vies por curva e forma com IC 95% (2,5-97,5 das
   realizacoes); vies medio ponderado dos 4 HJ como na cadeia (pesos
   1/sig_total^2 de hj_calibracao) e por realizacao -> IC; diferenca (ii) - (i)
   e (iii) - (i) com IC (bootstrap sobre curvas e realizacoes); razao dispersao
   sintetica entre temporadas / real.
7. Saidas: vies_forma_swasp.parquet (curva x forma x realizacao),
   vies_forma_swasp_resumo.parquet (curva x forma), vies_forma_swasp.json.

## Numeros esperados

- Controle: 4 de 4 HJ a < 0,05 min do O-C gravado; 61 de 61 temporadas dos 26
  a < 0,05 min (mesmo codigo, mesmo cache; desvio aqui = o cache ou o codigo
  mudou desde a cadeia, e para tudo).
- (i) vies do transito sintetico por HJ: **|media| < 0,3 min** (forma simetrica;
  a amostragem real quebra a simetria pouco); dispersao entre realizacoes
  **0,5-2 min** (da ordem da barra formal + temporadas reais). Media ponderada
  dos 4: **0,0 +- 0,4 min** -> NAO reproduz -2,14 +- 0,78: o -2,14 nao e
  procedural de forma. Vai a 3.4 com destaque. Hipoteses que sobram (a
  listar, nao a testar aqui): sistema de tempo (TDB - UTC = 65 s em 2006-2008;
  -2,14 min ~ 2x), propagacao TESS dos HJ (WASP-4b tem decaimento publicado,
  B8), astrofisica dos 4.
- (ii) eclipse destacado: **|media| < 0,3 min** por alvo; (iii) contato:
  **|media| < 0,5 min** (a base elipsoidal curva contra o trapezio plano, mas
  simetrica); diferencas (ii) - (i) e (iii) - (i): **compativeis com 0 a 2
  sigma, |dif| < 0,5 min**. Dispersao entre realizacoes maior em (iii)
  (profundidade efetiva menor sobre base curva): x1,2-2 de (ii).
- Dispersao sintetica entre temporadas / real: **0,2-0,8** (o real tem
  astrofisica e sistematica de temporada que o bootstrap de noite nao gera; em
  232634196 o real e 19,9 min contra formal 0,5).
- Se (i) der vies de -1 a -3 min: o -2,14 E procedural (trapezio sobre
  transito com amostragem real) e 3.4 fica como esta com o mecanismo nomeado;
  se (ii)/(iii) diferirem de (i) em > 1 min, o vies dos HJ nao transfere as
  binarias e a correcao de 3.4 muda por forma. Desvio em qualquer direcao vai ao
  relatorio.
