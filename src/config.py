"""Parametros globais, dimensionados a partir da sondagem da maquina.

Sondagem (2026-09-06): i5-13420H, 8 nucleos fisicos / 12 logicos, 31.7 GB RAM,
76.7 GB livres em C:, Python 3.12.14 no .venv.
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# SANDBOX: isola cache, resultados e fila num diretorio proprio, mantendo os
# artefatos caros (catalogos e alvo-set) compartilhados e somente-leitura.
#
# Existe porque testar a retomada exige matar o processo, e rodar isso contra o
# cache e o banco de uma execucao real produziria DOIS escritores - violando a
# invariante de escritor unico que `state.reconcile` assume - e faria o
# reconcile do teste apagar como orfaos os FITS em voo da execucao real.
_SANDBOX = os.environ.get("NASA_SANDBOX")
DATA = ROOT / "data"
CACHE = DATA / "cache"            # FITS de light curve, apagados apos processar
TPF_CACHE = DATA / "cache" / "tpf"  # target pixel files, fila tardia (30-50 MB cada)
# Separacao deliberada por CUSTO DE REPRODUCAO, nao por tipo de conteudo.
#
# CATALOGS e TARGETSET guardam o que e caro de refazer: o dump do VSX sozinho
# leva 26 min. RESULTS guarda o que se regenera rodando o pipeline de novo.
# Estavam juntos, e um `rm data/results/*.parquet` levou o targetset junto -
# barato naquele momento porque os catalogos estavam em cache, caro na rodada
# do setor inteiro. Limpeza de resultado nunca pode alcancar artefato caro.
CATALOGS = DATA / "catalogs"      # catalogos de exclusao, cacheados (caro)
TARGETSET = DATA / "targetset"    # alvo-set e crossmatch do TIC (caro)
RESULTS = DATA / "results"        # metricas, curvas, lotes (descartavel)
# CURVAS DE CONJUNTO ROTULADO, FORA DO ALCANCE DO `reconcile`.
#
# `state.reconcile` apaga todo FITS do CACHE que o banco nao reclama - e correto,
# e a invariante de escritor unico. Mas TODO conjunto de validacao deste projeto
# vem de alvos EXCLUIDOS do escopo (VSX, Gaia variavel, Prsa), e alvo excluido
# por definicao nunca entra na fila. Entao curva de validacao e sempre orfa, e
# morre no proximo arranque: aconteceu duas vezes em duas rodadas, a ultima
# levando 6.152 arquivos.
#
# Depender de ordem de execucao para nao perder trabalho e o mesmo erro que
# juntar RESULTS com TARGETSET. Aqui a separacao e por QUEM PODE APAGAR.
VALIDACAO = DATA / "validacao"    # rotulados: caros, e ninguem os reclama

if _SANDBOX:
    # so o que e descartavel e reescrito muda de lugar; catalogos e alvo-set
    # continuam apontando para os originais, lidos e nunca escritos
    _SBOX = Path(_SANDBOX)
    CACHE = _SBOX / "cache"
    TPF_CACHE = _SBOX / "cache" / "tpf"
    RESULTS = _SBOX / "results"
REPORTS = ROOT / "reports"
STATE_DB = (Path(_SANDBOX) / "state.db") if _SANDBOX else (ROOT / "state.db")

for _d in (CACHE, TPF_CACHE, CATALOGS, TARGETSET, RESULTS, REPORTS, VALIDACAO):
    _d.mkdir(parents=True, exist_ok=True)

# --- paralelismo -----------------------------------------------------------
# Downloads: I/O-bound, teto do MAST.
N_DOWNLOAD_THREADS = 4
DOWNLOAD_BATCH = 200
# Busca: CPU-bound. 6 processos deixa 2 nucleos fisicos de folga.
# NUMBA_NUM_THREADS=1 e obrigatorio: sem isso cada processo abre 8 threads
# numba e a maquina thrasha.
N_SEARCH_WORKERS = 6
NUMBA_THREADS_PER_WORKER = 1

# --- disco -----------------------------------------------------------------
DISK_CAP_GB = float(os.environ.get("NASA_DISK_CAP_GB", 20.0))
TPF_CAP_GB = 8.0   # subteto dentro do teto geral

# --- MAST ------------------------------------------------------------------
MAST_MAX_RETRIES = 5
MAST_BACKOFF_BASE = 2.0
MAST_BACKOFF_CAP = 120.0

# --- selecao de alvos ------------------------------------------------------
TMAG_MAX = 13.0
# CROWDSAP e FLAG, nao exclusao. Diluicao por vizinha nao invalida a deteccao,
# invalida a ATRIBUICAO do sinal aquela estrela - e isso e pergunta da Fase 6,
# que decide com a profundidade observada na mao. No lote de 50, 30% dos alvos
# ficaram abaixo de 0.8: descartar um terco do setor por uma medida que nem
# sequer diz que o sinal e falso seria jogar fora candidato bom.
CROWDSAP_MIN = 0.8          # so verificavel apos o download (header do FITS)
EXCLUSION_RADIUS_ARCSEC = 20.0   # cruzamento posicional -> exclui
BLEND_RADIUS_ARCSEC = 42.0       # 2 pixels TESS -> apenas marca blend

# --- pre-processamento -----------------------------------------------------
# Limiar de planicidade do teste de forma, CALIBRADO no gate de 80 objetos
# (37 binarias do Prsa, 31 planetas CP/KP), nao escolhido a olho.
#
# O valor inicial era 0.25 e nunca disparava: a planicidade medida fica entre
# 0.5 e 0.9 nas duas classes, e so 1 de 68 objetos caia abaixo de 0.25. O teste
# parecia conservador e estava inerte.
#
# Varredura do poder discriminante (rej_binarias - rej_planetas):
#   0.50 -> +8%   IC [-4%, +18%]   (cruza zero)
#   0.55 -> +16%  IC [+2%, +28%]
#   0.60 -> +57%  IC [+37%, +70%]  <- escolhido: 0% de planeta perdido
#   0.65 -> +58%  IC [+37%, +73%]  (mesmo poder, mas ja perde 6% dos planetas)
#   0.70 -> +56%  IC [+33%, +72%]  (perde 23%)
# 0.60 e 0.65 empatam dentro do ruido; fica o que nao derruba planeta.
SHAPE_FLATNESS_MIN = 0.60
# AJUSTADO IN-SAMPLE. O limiar foi varrido no MESMO conjunto em que o poder de
# +57% foi reportado, entao esse numero e otimista por selecao. O valor honesto
# so aparece quando este limiar encontrar as binarias de OUTRO setor - trate
# +57% como teto, nao como desempenho esperado.

# Razao maxima entre a largura a meia profundidade e a duracao do TLS antes de
# o resultado ser considerado inconsistente. Era 1.5 e descartava 9 objetos do
# gate como DESCONHECIDO: a duracao do TLS subestima sistematicamente a largura
# a meia profundidade da curva dobrada, e o corte apertado transformava isso em
# perda de dado em vez de medida.
SHAPE_W50_MAX_RATIO = 3.0

WOTAN_METHOD = "biweight"
WOTAN_WINDOW_D = 0.5
SIGMA_CLIP_UPPER = 5.0      # so para cima; cortar para baixo apaga o transito
BIN_MINUTES_FOR_BLS = 10.0  # transitos duram horas; principal ganho de wall-clock

# --- busca -----------------------------------------------------------------
LS_PERIOD_MIN_D, LS_PERIOD_MAX_D = 0.05, 20.0
LS_FAP_MAX = 1e-4
LS_AMPLITUDE_SNR_MIN = 3.0
# SNR local do pico: potencia do pico sobre o continuo local em volta dele.
# A FAP do Lomb-Scargle assume ruido branco, e residuo de PDCSAP e ruido
# vermelho - contra ele, FAP < 1e-4 passa quase toda curva. Este e o corte que
# discrimina de verdade: sistematica levanta uma banda larga, sinal real faz um
# pico estreito acima do proprio continuo.
# Limiar do SNR local, CALIBRADO por ROC contra conjunto rotulado (VSX), nao
# escolhido a olho. O valor fixo anterior era 8.0 - quatorze vezes menor que o
# calibrado, e por isso nao filtrava nada. Ver src/validate_variability.py.
#
# E funcao do RUIDO, nao um numero unico: o SNR local de fundo cresce com o
# ruido do alvo, e um limiar fixo aprovava 4.9% dos negativos no regime de
# ruido baixo contra 26.8% no regime alto - fator 5.5 de desequilibrio, com o
# excesso todo caindo nos alvos fracos, que e onde o ruido passa por sinal.
# Com o limiar adaptativo a aprovacao fica em ~10% em todas as faixas.
#
#   limiar(ruido_ppm) = LS_LOCAL_SNR_A * ruido_ppm ** LS_LOCAL_SNR_ALPHA
LS_LOCAL_SNR_A = 0.0378
LS_LOCAL_SNR_ALPHA = 1.063
# Minimo de ciclos dentro do baseline para o periodo estar restrito pelos
# dados. Num setor de ~27 d, 3 ciclos limitam o ramo A a P <~ 9 d: acima
# disso a variabilidade pode ser real, mas o periodo nao e mensuravel com um
# setor so e a deteccao precisa de outro setor para se sustentar.
MIN_CYCLES_FOR_PERIOD = 3.0
# Fracao maxima do lote que pode compartilhar o mesmo periodo antes de ele ser
# tratado como sistematica de campo: variavel real nao e compartilhada entre
# alvos, sistematica e.
# Excesso minimo sobre a densidade tipica de periodos para um bin contar como
# sistematica de campo. Substitui um limiar de FRACAO ABSOLUTA (2% do total),
# que funcionava com 200 alvos e virava inalcancavel com 6.017: os periodos se
# espalham por ~290 bins e nenhum chega a 2% do total nem tendo 3,6x o esperado.
# O excesso relativo nao depende do tamanho da amostra.
SHARED_PERIOD_MIN_EXCESS = 2.5
SHARED_PERIOD_MAX_FRACTION = 0.02   # mantido para compatibilidade
SHARED_PERIOD_TOLERANCE = 0.01

BLS_PERIOD_MIN_D, BLS_PERIOD_MAX_D = 0.5, 15.0
TLS_TOP_FRACTION = 0.01     # TLS so nos ~1% melhores por SDE
SDE_MIN = 9.0          # entra no vetting e no parquet: nada e descartado aqui

# CORTE DE APRESENTACAO do relatorio - quem vira PAGINA, nao quem vira dado.
# O gargalo real do projeto nao e CPU, e revisao humana: com SDE 9 no setor
# inteiro o relatorio teria da ordem de mil paginas, e ninguem olha mil paginas.
#
# Calibrado contra a distribuicao medida de SDE (39 CP/KP do setor 105 contra
# 32 alvos aleatorios do alvo-set):
#   confirmados   p10 =  9.8   mediana = 17.7   p90 = 27.3
#   aleatorios    p10 =  4.3   mediana =  5.6   p90 =  8.2   max = 12.8
# As distribuicoes quase nao se sobrepoem: o fundo morre onde os confirmados
# comecam.
#
#   SDE>= 9  retem 90% dos confirmados [IC 76-97%]   <- escolhido
#   SDE>=12  retem 82%                 [IC 66-92%]
#   SDE>=15  retem 74%                 [IC 58-87%]
#   SDE>=18  retem 41%                 [IC 26-58%]   perda alta demais
#
# O CORTE FOI DE 12 PARA 9 DEPOIS DE MEDIR EM VOLUME. Com 200 alvos ORDINARIOS
# do tier 1 - nao mais o conjunto rotulado e enriquecido - a distribuicao de
# fundo ficou:  p50 = 5.5   p90 = 7.5   p99 = 9.1   max = 10.0
# e passam do corte:
#   SDE>= 9  ->  4/200 = 2.0%  IC [0.5%, 5.0%]  ->  ~120 candidatos no tier 1
#   SDE>=12  ->  0/200 = 0.0%  IC [0.0%, 1.8%]  ->  0 a 110
# Com 12 NENHUM alvo ordinario passa, e 120 paginas cabem folgadamente no teto
# de 200. Manter 12 custaria 8 pontos de retencao de planeta confirmado (90%
# contra 82%) para resolver um problema de volume que a medicao mostrou nao
# existir - e a perda cairia justamente sobre transito raso, que e onde estao
# os planetas pequenos que o projeto procura. Os CP/KP que calibraram o corte
# original tem mediana de SDE 17.7 porque sao transitos faceis por construcao;
# eles nao representam um planeta pequeno em ana M.
#
# RESSALVA que precisa ser dita: com apenas 32 alvos aleatorios, "0 de 32" acima
# de SDE 15 NAO significa zero. O IC de 95% ainda admite ~11% do alvo-set, isto
# e, ate ~1.200 candidatos. A projecao pontual de "0 paginas" seria falsa
# precisao, e por isso o corte nao foi empurrado para 15. O numero de paginas
# real da Fase 7 so se conhece rodando o setor.
#
# Nada e perdido por este corte: quem fica entre 9 e 12 continua em
# results_full.parquet com todos os vereditos, apenas nao gera pagina.
REPORT_SDE_MIN = 9.0

# Teto de paginas do relatorio. Se o corte de SDE produzir mais que isto, o
# relatorio SOBE o corte sozinho ate caber e REGISTRA que subiu - gerar 1.200
# paginas que ninguem vai abrir e o mesmo que nao gerar relatorio nenhum, e a
# projecao de quantas paginas o setor produz tem IC largo demais (com 32 alvos
# aleatorios, "0 de 32" ainda admite ~11% do alvo-set) para confiar num numero
# fixo escolhido de antemao.
REPORT_MAX_PAGES = 200
MIN_TRANSITS = 3

# --- sistematicas do TESS a rejeitar --------------------------------------
TESS_ORBIT_PERIOD_D = 13.7
# Tolerancia POR periodo, nao uma so. 13.7 d (orbita do TESS) e 6.85 d (seu
# harmonico) sao especificos do instrumento e merecem janela larga. Ja 1 d nao
# tem no espaco a forca que tem em survey terrestre, e uma janela de 5% ali
# cobre 0.95-1.05 d, faixa povoada por hot Jupiters reais: no gate de 80
# objetos esse corte rejeitou um planeta CONFIRMADO (P = 1.0077 d) e nenhuma
# binaria. Poder negativo. Janela estreita para ele, e 0.5 d sai da lista.
# O harmonico de 6.85 d saiu de 5% para 1%: com 5%, a janela cobria 6,51-7,19 d
# e rejeitou 53 dos 317 candidatos do setor - 17% deles - numa faixa densamente
# povoada por planetas reais. Mesmo defeito que o corte de 1 d ja tinha tido.
# 13.7 d fica com 5% porque a orbita e um efeito forte e a faixa e menos
# povoada; 1 d e 6.85 d ficam estreitos.
SYSTEMATIC_PERIODS_TOL = ((13.7, 0.05), (6.85, 0.01), (1.0, 0.005))
SYSTEMATIC_PERIODS_D = tuple(p for p, _ in SYSTEMATIC_PERIODS_TOL)
SYSTEMATIC_TOLERANCE = 0.05  # usado para a duracao do setor


def apply_worker_env():
    """Chamar no inicio de cada processo worker, antes de importar numba."""
    n = str(NUMBA_THREADS_PER_WORKER)
    for var in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ.setdefault(var, n)
