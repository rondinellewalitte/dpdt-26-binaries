# -*- coding: utf-8 -*-
"""Item 3 da revisao externa (2026-09-13): o arquivo TESS completo dos 26.

O censo de `vies_cache_288` mostrou que o arquivo tem mediana de 39 setores
SPOC 2-min para os 26 alvos medidos, contra 4 no cache dos 22 manifestos:
824 setores no total (6-43 por alvo), ~1,6 GB. Com 25-40 epocas primarias
por alvo, a dispersao entre setores dos PRIMARIOS e medida em vez de
extrapolada dos secundarios de 5 alvos, e a "estrutura interna ao TESS" das
dez deteccoes vira testavel em 6 anos. E download, nao metodo.

O que este script faz: para cada um dos 26, pede ao MAST a lista de produtos
LC SPOC 2-min de TODOS os setores e escreve
`data/catalogs/manifest_suplementar_completo26_2min.parquet` - que
`oc_lote.manifestos()` ja apanha pelo glob `manifest_suplementar*`. O
download acontece dentro da cadeia (`oc_lote.baixar_alvo`), por URL, com
falha alta. Escopo: os 26 (o veiculo), nao os 88 - os 62 nao medidos
poderiam recuperar com mais setores, mas isso muda a definicao da amostra
(o passo 3 do funil) e e outra decisao.

Depois: `oc_lote.py --alvos data/orquestra/alvos_oc_88.parquet` (os 62 seguem
parando nas guardas; os 26 ganham os setores novos), e toda a cadeia de
scripts de medida.

EXPECTATIVA (escrita e commitada ANTES de rodar), sobre a cadeia com os
manifestos completos e as barras SuperWASP ja corrigidas:
  - MAST responde para 26 de 26; o manifesto tem ~800 curvas (o censo
    contou 824 setores; alguns podem nao ter produto LC).
  - a escada de periodo E TERRITORIO NOVO: desenhada incremental e testada
    em 4-7 epocas; com 30-40 e setores de dispersao alta (o caso dos 128 min
    do Apendice A) a tolerancia pode nao fechar. Expectativa: fecha em >= 20
    dos 26; os que nao fecharem sao relatados com o setor que quebrou, nao
    forcados. Os aglomerados por setor sobem de 4-7 para 25-40 por alvo.
  - alavanca media do TESS cai de 0,62 para ~0,1: a decomposicao passa a ter
    poder no TESS. Expectativa para o fator de variancia TESS: 1,5-3 (a
    deriva entre setores de 1-2,6 min/ano vista nos secundarios de 5 alvos
    nao esta na barra intra-setor de nenhum). Fator SuperWASP: o mesmo da
    rodada anterior dentro do intervalo (as barras SuperWASP nao mudam).
  - mediana do chi2_red da parabola com ~30 graus de liberdade: 1,5-3 - e
    ESSA e a medida direta da dispersao entre setores dos primarios,
    equivalente a 1-3 min por setor em torno da parabola.
  - deteccoes: com 4 epocas TESS e 3 parametros a estrutura interna vira
    coeficiente; com 30 vira residuo. Das 10 deteccoes internas ao TESS,
    expectativa: 3-6 sobrevivem a curvatura + adversarial; TIC 232634196
    (alavanca de 20 anos) sobrevive. Deteccao que some e resultado, nao
    perda. Nenhuma deteccao NOVA acima de 2 de 26 (o que 30 epocas
    permitem ver e curvatura interna aos 6 anos de TESS, e a ordem de
    grandeza 1e-7 d/ano curva ~0,5 min em 6 anos: abaixo da barra).
  - CROWDSAP (cabecalho SPOC) por setor: mediana > 0,9 nos 26 (Tmag 8,8-13,8);
    2-4 alvos com algum setor < 0,8; relatado por alvo, ao lado da epoca.
  - controle TIC 142874476 (A e B) nao se move: esta so num manifesto (s105)
    e o caminho do controle nao passa por este manifesto.
Desvio em qualquer direcao vai para o humano antes de qualquer outra coisa.
"""
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
# DECISAO (2026-09-13): o item 3 fica para DEPOIS do arXiv (a submissao nao depende dele; um referee que o
# pedir recebe-o com o pedido na mao). O manifesto foi gerado (799 curvas) e esta em data/catalogs/pendente/,
# FORA do glob `manifest_suplementar*` de oc_lote.manifestos(): mover para data/catalogs/ com o nome abaixo
# e o que liga o arquivo completo na cadeia. Com ele ligado, a cadeia baixa ~1,6 GB dentro de baixar_alvo.
MAN = config.CATALOGS / "pendente" / "completo26_2min.parquet"        # ativo: config.CATALOGS / "manifest_suplementar_completo26_2min.parquet"


def produtos_lc(tic):
    """Todos os produtos LC SPOC 2-min do alvo no MAST, com URL de download."""
    from astroquery.mast import Observations
    ult = None
    for tentativa in range(3):
        try:
            o = Observations.query_criteria(target_name=str(tic), obs_collection="TESS",
                                            dataproduct_type="timeseries").to_pandas()
            o = o[(o.t_exptime == 120) & (o.provenance_name == "SPOC")]
            if not len(o):
                return []
            prods = Observations.get_product_list(Observations.query_criteria(
                obsid=list(o.obsid.astype(str)))).to_pandas()
            prods = prods[prods.productSubGroupDescription == "LC"]
            out = []
            for r in prods.itertuples():
                nome = Path(r.dataURI).name
                setor = int(nome.split("-s")[1][:4])
                out.append({"tic": int(tic), "filename": nome, "sector": setor,
                            "url": "https://mast.stsci.edu/api/v0.1/Download/file?uri=" + r.dataURI})
            return out
        except Exception as e:  # noqa: BLE001 - contado e repetido, nao escondido
            ult = e
            time.sleep(3.0 * (tentativa + 1))
    raise RuntimeError(f"MAST falhou 3x para TIC {tic}: {str(ult)[:80]}")


if __name__ == "__main__":
    t26 = pd.read_parquet(config.DATA / "orquestra" / "oc_lote" / "tabela_26.parquet")
    censo = pd.read_parquet(config.DATA / "results" / "vies_cache_288.parquet").set_index("tic")
    linhas, falhas = [], []
    t0 = time.time()
    for i, tic in enumerate(t26.TIC, 1):
        try:
            p = produtos_lc(int(tic))
        except Exception as e:  # noqa: BLE001
            falhas.append((int(tic), str(e)[:80]))
            print(f"[{i}/26] TIC {tic}: FALHA {str(e)[:80]}", flush=True)
            continue
        setores = sorted({x["sector"] for x in p})
        n_censo = int(censo.n_real.get(int(tic), -1))
        print(f"[{i}/26] TIC {tic}: {len(p)} curvas LC em {len(setores)} setores (censo: {n_censo})  {(time.time() - t0) / 60:.1f} min", flush=True)
        linhas.extend(p)
    man = pd.DataFrame(linhas).drop_duplicates(["tic", "filename"]).reset_index(drop=True)
    man.to_parquet(MAN, index=False)
    print(f"\nmanifesto completo: {len(man)} curvas, {man.tic.nunique()} alvos, {man.sector.nunique()} setores distintos -> {MAN.name}")
    if falhas:
        print(f"FALHAS de MAST: {falhas}")
        sys.exit(1)
