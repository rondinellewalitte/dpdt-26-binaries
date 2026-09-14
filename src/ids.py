"""Identificador nunca passa por ponto flutuante - como asserção, não comentário.

O defeito que motivou este módulo custou uma fase inteira sem produzir sintoma:
`pd.to_numeric(col, errors="coerce")` sobre uma coluna `object` que contenha
QUALQUER ausente devolve float64. float64 tem 15-16 dígitos significativos, e
todo `source_id` do Gaia DR3 tem 19. O TIC devolvia 6662513802550562048 e o
parquet guardava 6662513802550561792 - 256 a menos, um identificador que não
existe. 89% dos IDs gravados não existiam no `gaiadr3.gaia_source`, e o único
que existia apontava para uma estrela a 32" do alvo.

O que torna esse defeito pior que os outros nove do projeto: ele produziu um
resultado PLAUSÍVEL. "11% dos alvos têm contrapartida Gaia" parece um número
astrofísico - completude de catálogo, corte de magnitude, qualquer coisa. Foi
lido como tal. Todos os outros defeitos produziram algo estranho o bastante
para investigar.

O GATILHO EXATO, medido (não suposto):

    object + np.nan   -> float64   CORROMPE
    object + None     -> float64   CORROMPE
    object + ""       -> float64   CORROMPE
    object, só ints   -> float64   CORROMPE   (já-inteiro-Python também!)
    object, sem NA    -> int64     ok
    string dtype      -> Int64     ok, mesmo com pd.NA

A `astroquery` devolve tabelas mascaradas, que viram `object` com ausentes.
Todo caminho de catálogo deste projeto passa por ali.

POR QUE MAGNITUDE NÃO É A DEFESA. TIC (11 dígitos) e OID do VSX (8) atravessam
float64 sem perda, porque estão abaixo de 2^53 = 9.007.199.254.740.992. Isso é
verdade hoje e é uma propriedade do catálogo, não do código: nada impede um TIC
de 17 dígitos amanhã, e nada avisaria. A defesa é estrutural - coluna de ID com
dtype de ponto flutuante é erro, independente dos valores que contenha.
"""
import re

import numpy as np
import pandas as pd

FLOAT64_EXACT_MAX = 2 ** 53   # 9.007.199.254.740.992

# Nomes que SAO identificador. A deteccao automatica cobre tambem qualquer
# coluna terminada em `_id`, `_oid` ou `_source_id` - a convencao do projeto.
ID_NAMES = frozenset({"tic", "id", "oid", "source_id", "gaia_source_id",
                      "vsx_oid", "eb_source_id", "blend_de", "tic_vizinho"})


def id_columns(df):
    """Colunas que sao identificador, por nome."""
    return [c for c in df.columns
            if str(c).lower() in ID_NAMES
            or re.search(r"(^|_)(id|oid|source_id)$", str(c).lower())]


def to_int_id(series, name="id", strip_prefix=None):
    """object/str/Int -> Int64 sem tocar em float, com round-trip verificado.

    Aceita 'TIC 12345', '12345', '12345.0' e 12345. O que não casar com um
    inteiro vira <NA> - ausência explícita, nunca um valor arredondado.
    """
    if series is None:
        return pd.array([], dtype="Int64")
    s = pd.Series(series)
    if pd.api.types.is_float_dtype(s.dtype):
        # Chegou aqui já como float: o estrago, se houve, é anterior a este
        # módulo. Só é seguro se todo valor for exato em float64.
        big = s.abs() >= FLOAT64_EXACT_MAX
        if big.any():
            raise ValueError(
                f"IDENTIFICADOR CORROMPIDO [{name}]: a coluna chegou como "
                f"{s.dtype} com {int(big.sum())} valor(es) >= 2^53. "
                "A conversão para float já perdeu dígitos - conserte na origem, "
                "convertendo de string direto para Int64.")
    s = s.astype("string").str.strip()
    if strip_prefix:
        s = s.str.replace(strip_prefix, "", regex=False).str.strip()
    s = s.str.replace(r"\.0*$", "", regex=True)          # '12345.0' -> '12345'
    ok = s.str.fullmatch(r"[+-]?\d+").fillna(False)
    out = s.where(ok).astype("Int64")

    # Round-trip: o valor gravado tem que reproduzir a string de origem.
    back = out.astype("string")
    bad = ok & (back != s)
    if bad.any():
        i = int(np.flatnonzero(bad.values)[0])
        raise ValueError(
            f"IDENTIFICADOR CORROMPIDO [{name}]: round-trip falhou em "
            f"{int(bad.sum())} linha(s). Origem {s.iloc[i]!r}, gravaria "
            f"{back.iloc[i]!r}.")
    return out


def assert_int_ids(df, cols=None, where="", verificar_fonte=False):
    """Antes de gravar: nenhuma coluna de ID pode ter dtype de ponto flutuante.

    Chame imediatamente antes de todo `to_parquet` que contenha identificador.

    O dtype sozinho NÃO teria pego o defeito do Gaia: o código era
    `to_numeric(...).astype("Int64")`, e o float existia só entre as duas
    chamadas - o que chegava ao parquet já era Int64, com os valores errados
    dentro. Por isso vem junto a checagem de rastro de arredondamento, que
    olha os valores e não o tipo.
    """
    ctx = f" ({where})" if where else ""
    for c in (id_columns(df) if cols is None else cols):
        if c not in df.columns:
            continue
        dt = df[c].dtype
        if pd.api.types.is_float_dtype(dt):
            v = pd.to_numeric(df[c], errors="coerce")
            big = int((v.abs() >= FLOAT64_EXACT_MAX).sum())
            raise TypeError(
                f"IDENTIFICADOR EM PONTO FLUTUANTE [{c}]{ctx}: dtype={dt}. "
                f"{big} valor(es) acima de 2^53 já perderam dígitos. "
                "Use ids.to_int_id() a partir da string original.")
        if not (pd.api.types.is_integer_dtype(dt) or dt == "string"
                or pd.api.types.is_object_dtype(dt)):
            raise TypeError(f"IDENTIFICADOR COM DTYPE INESPERADO [{c}]{ctx}: {dt}")
        _assert_no_rounding_trace(df[c], c, ctx,
                                  verificar_fonte=verificar_fonte)
    return df


def read_parquet(path, source_cols=None, **kw):
    """`pd.read_parquet` com as mesmas checagens da gravacao.

    A asseracao na gravacao nao basta: parquet escrito por versao anterior do
    codigo - ou por outra ferramenta - so e detectavel ao ser LIDO. Foi o caso
    dos tres arquivos que a varredura acusou, todos gravados antes de ids.py
    existir. Um deles, `targetset.parquet`, era lido por praticamente todo o
    resto do pipeline sem nenhuma checagem no caminho.

    `source_cols` mapeia coluna de ID -> coluna com a string bruta preservada
    (ex.: {"gaia_source_id": "gaia_raw"}). Quando a origem existe, a
    verificacao deixa de ser estrutural e passa a ser contra a fonte, que e a
    unica evidencia forte: 12.943 de 12.943 exatos no tic_xmatch.
    """
    df = pd.read_parquet(path, **kw)
    assert_int_ids(df, where=f"leitura de {getattr(path, 'name', path)}")
    for id_col, raw_col in (source_cols or {}).items():
        if id_col in df.columns and raw_col in df.columns:
            assert_matches_source(df, id_col, raw_col,
                                  ctx=f" (leitura de {getattr(path, 'name', path)})")
    return df


def assert_matches_source(df, id_col, raw_col, ctx=""):
    """O ID gravado reproduz a string de origem preservada ao lado dele?

    Verificacao contra a FONTE, nao contra estrutura nem contra taxa de
    casamento - e taxa de casamento foi exatamente o que leu um defeito como
    resultado astrofisico.
    """
    ok = df[id_col].notna() & df[raw_col].notna()
    if not ok.any():
        return df
    raw = df.loc[ok, raw_col].astype("string").str.strip()
    back = df.loc[ok, id_col].astype("Int64").astype("string")
    bad = raw != back
    if bad.any():
        i = int(np.flatnonzero(bad.values)[0])
        raise ValueError(
            f"IDENTIFICADOR NAO REPRODUZ A ORIGEM [{id_col}]{ctx}: "
            f"{int(bad.sum()):,} de {int(ok.sum()):,} divergem. "
            f"Origem {raw.iloc[i]!r}, gravado {back.iloc[i]!r}.")
    return df


def verificar_no_gaia(ids_amostra, timeout=180):
    """Os source_id existem mesmo no Gaia DR3? Devolve (n_existem, n_testados).

    E a unica checagem que separa ID genuino de ID arredondado, porque o padrao
    de bits NAO separa - ver `_assert_no_rounding_trace`. Um valor arredondado
    por float64 cai num inteiro que quase nunca corresponde a uma fonte real.
    """
    import pyvo
    ids_amostra = [int(x) for x in ids_amostra]
    if not ids_amostra:
        return 0, 0
    svc = pyvo.dal.TAPService("https://gea.esac.esa.int/tap-server/tap")
    q = ("SELECT source_id FROM gaiadr3.gaia_source WHERE source_id IN (%s)"
         % ",".join(str(x) for x in ids_amostra))
    d = svc.search(q).to_table().to_pandas()
    return len(d), len(ids_amostra)


def _assert_no_rounding_trace(col, name, ctx="", verificar_fonte=False,
                              n_amostra=40, frac_min=0.90):
    """Rastro de float: bits baixos zerados num valor acima de 2^53.

    O TESTE DE BITS SOZINHO NAO SERVE, E ISSO FOI MEDIDO. A premissa era que os
    bits baixos de um identificador genuino sao aleatorios, entao ser multiplo
    de 2^(b-53) teria probabilidade 2^(53-b). Falso para o Gaia: consultando
    2.094 fontes reais do DR3 no campo do setor 20, direto no arquivo do Gaia e
    sem MAST no caminho, 100% sao multiplas de 2^(b-53), com minimo de SETE
    zeros a direita e nenhuma abaixo.

    E o efeito depende da REGIAO DO CEU. No campo sul do setor 105 os ids tem 62
    bits, o limiar vira 2^9 e so 12,2% disparam; no campo norte do setor 20 eles
    tem 60 bits, o limiar vira 2^7 e 100% disparam. A guarda passava no sul por
    acidente de magnitude, nao por estar certa.

    Entao o teste de bits fica como TRIAGEM BARATA e quem decide e a fonte: com
    `verificar_fonte=True` uma amostra vai ao Gaia DR3 e o veredito e se os ids
    existem. ID arredondado quase nunca corresponde a fonte real - foi assim que
    o defeito original apareceu, com 6662513802550561792 no lugar de
    6662513802550562048.

    Sem `verificar_fonte`, avisa e nao levanta: levantar seria trocar um falso
    negativo por um falso positivo garantido no hemisferio norte.
    """
    v = pd.Series(col)
    if pd.api.types.is_integer_dtype(v.dtype):
        mx = v.abs().max()
        if pd.isna(mx) or mx < FLOAT64_EXACT_MAX:
            return
    v = pd.to_numeric(v, errors="coerce").dropna()
    if v.empty:
        return
    try:
        v = v.astype("int64").abs()
    except (TypeError, ValueError, OverflowError):
        return
    big = v[v >= FLOAT64_EXACT_MAX]
    if len(big) < 20:
        return
    bits = big.map(int).map(lambda x: x.bit_length())
    need = (bits - 53).clip(lower=1)
    trailing = big.map(int).map(lambda x: (x & -x).bit_length() - 1)
    frac = float((trailing >= need).mean())
    if frac <= 0.5:
        return
    if not verificar_fonte:
        print("  [ids] AVISO %s%s: %.0f%% dos %d valores acima de 2^53 tem "
              "rastro de bits. Isso e estrutural no Gaia e NAO decide sozinho; "
              "rode com verificar_fonte=True para checar contra o arquivo."
              % (name, ctx, frac * 100, len(big)), flush=True)
        return
    amostra = big.sample(min(n_amostra, len(big)), random_state=0).tolist()
    n_ok, n_tot = verificar_no_gaia(amostra)
    if n_tot and n_ok / n_tot < frac_min:
        raise ValueError(
            "IDENTIFICADOR CORROMPIDO [%s]%s: so %d de %d source_id sorteados "
            "EXISTEM no Gaia DR3. Combinado com o rastro de bits em %.0f%% da "
            "coluna, isso e arredondamento por float64. Exemplo: %d. Converta "
            "de string direto para Int64 (ids.to_int_id)."
            % (name, ctx, n_ok, n_tot, frac * 100, int(amostra[0])))
    print("  [ids] %s%s: rastro de bits em %.0f%%, mas %d de %d source_id "
          "existem no Gaia DR3 - estrutural, nao corrupcao."
          % (name, ctx, frac * 100, n_ok, n_tot), flush=True)




def read_results(sector, path=None, **kw):
    """Le `results.parquet` EXIGINDO o setor. Sem default.

    OITAVA OCORRENCIA DO CAMINHO FIXO, e a mais perigosa das oito. Nas outras
    sete o caminho compartilhado SOBRESCREVIA e o sintoma era perda visivel -
    um arquivo com menos linhas do que devia. Aqui o pipeline ACUMULA: rodar o
    setor 20 deixou `results.parquet` com 26.335 linhas, 15.356 do s20 e 10.979
    do s105, corretamente separadas pela coluna `sector` e sem duplicata.

    Nada se perde, e por isso e pior. Quem ler sem filtrar recebe a UNIAO de
    dois campos cujas distribuicoes nao se parecem:

        Tmag mediana        9,53 (s20)  contra  11,89 (s105)
        amplitude mediana   0,66 mmag   contra   7,82 mmag
        SDE >= 9            26,4%       contra    5,3%

    E como o s20 e maior, qualquer taxa "do s105" calculada sobre o arquivo
    misturado vira media ponderada com o s20 dominando. A analise que mais sofre
    e justamente a comparacao ENTRE setores.

    Por isso `sector` e posicional e obrigatorio: quem esquecer recebe
    TypeError, nao numero plausivel. `sector=None` e permitido so explicitamente,
    para o caso raro de querer os dois - e ai a escolha esta escrita na chamada.
    """
    import config
    d = read_parquet(path or (config.RESULTS / "results.parquet"), **kw)
    if "sector" not in d.columns:
        raise RuntimeError("FALHA [read_results]: arquivo sem coluna `sector` - "
                           "nao da para separar os campos")
    if sector is None:
        return d
    out = d[d["sector"].astype("int64") == int(sector)]
    if not len(out):
        raise RuntimeError(
            "FALHA [read_results]: setor %s nao esta no arquivo. Presentes: %s"
            % (sector, sorted(d["sector"].dropna().astype(int).unique())))
    return out
