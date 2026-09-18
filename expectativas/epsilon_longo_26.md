# Expectativa: `src/epsilon_longo_26.py` (Passo 3)

Escrita e commitada ANTES da primeira execucao. Nada rodado. Estende
`epsilon_26.py` (mesma cadeia: `cadeia_oc.ajustar` linear e quadratica,
p_curv por Delta chi2, `cadeia_oc.adversarial`, p_gof com N - 3 gl; efemeride
exatamente linear nos instantes reais com as barras reais).

## O que o script faz

1. Populacao por FAIXA de P3 (log-uniforme dentro de cada faixa, N_DRAW = 2 000
   por alvo e faixa, semente fixa): 0,1-20, 20-50, 50-100 anos. e uniforme em
   0-0,8; omega uniforme em 0-2 pi; fase (anomalia media em t_ref) uniforme;
   M_bin U(1,5, 2,5), M3 U(0,1, 1,0), cos i uniforme - Apendice C.
2. LTTE excentrico (Irwin 1952): Delta t = A [ (1 - e^2)/(1 + e cos nu) sin(nu + omega)
   + e sin omega ], A = a12 sin i / c = `epsilon_26.amplitude_min` (a12 = a3 M3/(M_bin+M3)),
   Kepler resolvido por Newton, nu da anomalia excentrica.
3. Tres eficiencias por alvo e faixa: eps (barras fixas, sem ruido, deteccao =
   p_curv < 0,05 e p_adv < 0,05 - como epsilon_26); eps_re (sinal + ruido
   N(0, barra), barras RE-ESTIMADAS com `nulo_barras_reestimadas_26.barras_reestimadas`
   - a regra de completude_secular_26 - e os testes com elas); eps3 (barras
   fixas, deteccao e p_gof >= 0,05).
4. Controles: (A) `epsilon_26.epsilon_alvo` com a semente 20260913 reproduz
   Sigma eps = 6,37 e Sigma eps3 = 4,25 exatamente (o codigo e o mesmo); (B) o
   gerador novo com e = 0 forcado e faixa 0,1-20 da Sigma eps a < 0,15 de 6,37
   (ruido de Monte Carlo: sigma por alvo ~0,01, na soma ~0,05); (C) a orbita
   do V527 Dra (A 6,6 min, P3 2,733 a, e = 0) injetada pelo gerador novo nas 6
   epocas do TIC 424461577 em 360 fases reproduz `v527_ltte.py`: p_curv < 0,05
   em 93% +- 2, p_adv < 0,05 em 84% +- 2, dP/dt de -0,058 a +0,058.
5. Fracao por faixa na Tabela 5 do Tokovinin 2006 (VizieR J/A+A/450/681, cacheada
   em data/catalogs/tokovinin2006_table5.parquet), contada pelo script: para
   cada um dos 42 pares com P_in < 3 d, o terciario MAIS PROXIMO = o par cujo
   componente e o pai do par interno (campo Comp "primario,secundario,pai"),
   periodo minimo entre eles; "sem terciario" quando o pai e a raiz. Wilson 95%.
   Registra tambem a contagem FROUXA (qualquer par do mesmo sistema com P na
   faixa, mesmo de outro subsistema) porque e a que a nota usou ate aqui.
6. Numero esperado por faixa e total: f x Sigma eps, f x Sigma eps_re, f x Sigma eps3.
7. Saidas: epsilon_longo_26.parquet (alvo x faixa), epsilon_longo_26.json
   (controles, contagens, numeros esperados).

## Numeros esperados

- Tokovinin, contagem estrita (feita hoje na tabela cacheada, antes de escrever
  isto - e uma contagem de literatura, nao uma medida nossa): 42 pares com
  P < 3 d, 33 com terciario listado, 9 sem; terciario mais proximo em
  0,1-20 a: **8 de 42** (f 0,19, Wilson 0,10-0,33); 20-50 a: **0 de 42**
  (Wilson 0-0,084); 50-100 a: **3 de 42** (0,071, Wilson 0,025-0,19; HIP 47053
  64 a, 86187 90 a, 35487 93 a). A contagem frouxa da **10 de 42** em 0,1-20 a
  (HIP 16846 e 61910 entram por um par de OUTRO subsistema - 3,2 a e 0,12 a -
  que nao move o par eclipsante): e a que esta na nota (f = 0,24). DESVIO JA
  IDENTIFICADO: a nota conta 10 onde a definicao fisica (LTTE sobre o par
  eclipsante) da 8; vai ao relatorio, decisao do humano.
- Sigma eps por faixa (barras fixas): 0,1-20 a **6,0-6,8** (e ate 0,8 muda
  pouco: a fase e a amplitude dominam); 20-50 a **10-18** (A ~ 10-25 min, o
  arco de 16-20 anos e meia orbita e a parabola por tres grupos absorve o
  ciclo); 50-100 a **7-15** (A maior, mas a parte quadratica cai como
  P3^(-4/3): de 50 para 100 a, x0,4). Os dois alvos com eps = 0 na faixa
  curta (224605072, 229476285) sobem de 0 nas faixas longas (A >> barra).
- Sigma eps3 / Sigma eps: 0,1-20 a **0,60-0,72** (4,25/6,37 = 0,67 hoje); nas
  faixas longas **0,7-0,95** (o residuo da parabola e menor quando o arco e
  suave).
- Sigma eps_re / Sigma eps: 0,1-20 a **0,6-0,9** (5.2.5: o lado arquival
  encolhe); 20-100 a **0,3-0,7** (a deriva dentro do bloco SuperWASP de 1-4
  anos infla a barra re-estimada: 5-10 min para A 15 min e P3 30 a).
- Numeros esperados (f estrito x Sigma): 0,1-20 a **1,1-1,3** (fixas), 20-50 a
  **0** (limite superior Wilson x Sigma: **< 1,5**), 50-100 a **0,5-1,1**
  (Wilson 0,2-2,9); total fixas **1,6-2,4**, re-estimadas **1,0-1,6**; com
  eps3: total **1,2-1,9** fixas. Com f frouxo (0,24) a faixa curta da 1,5.
- Resposta a pergunta "at most one or two of the eleven": **sobrevive com
  barras re-estimadas (1-2) e fica no limite com barras fixas (2, ate 3 no
  Wilson superior de 50-100 a)** - a frase passa a "one to three" ou ganha a
  qualificacao "with re-estimated bars". Desvio em qualquer direcao vai ao
  relatorio.
