# -*- coding: utf-8 -*-
"""Contrato dos subagentes: as cinco regras onde dao para virar assercao.

REGRA 1  nenhum subagente ajusta limiar. `GuardaDeLimiares` fotografa toda
         constante MAIUSCULA numerica dos modulos de decisao antes da chamada e
         confere depois; diferenca e `LimiarAlterado`. Quando o limiar nao se
         aplica ao caso, o subagente levanta `Inaplicavel` - reporta e para.
REGRA 2  todo subagente devolve DISTRIBUICAO, nao conclusao. `Mediana` sem `n`
         nao constroi; `Taxa` sem intervalo nao constroi; `Comparacao` sem as
         duas amostras descritas (nome, n, selecao) nao constroi. Observacoes
         com vocabulario de conclusao sao rejeitadas pelo validador.
REGRA 3  nenhum artefato sem identidade do conjunto. `Identidade.artefato()` e
         o unico jeito de nomear saida, e o validador exige que cada caminho
         relatado carregue o rotulo E a tag do setor.
REGRA 4  o orquestrador nunca conclui descoberta: `Classe` tem CANDIDATO como
         valor maximo, e nao existe outro.
REGRA 5  contradicao com a expectativa escala para o humano, em qualquer
         direcao: `comparar_com_expectativa` nao sabe o que e "favoravel".

Este modulo e a parte que precisa ter sido VISTA FALHANDO: `test_contrato.py`
constroi cada violacao e exige a excecao.
"""
import copy
import enum
import json
import math
import re
import types
from dataclasses import dataclass, field, asdict
from pathlib import Path

import config


class ContratoViolado(Exception):
    """O subagente devolveu algo fora do contrato. Nao se corrige: para."""


class LimiarAlterado(ContratoViolado):
    """Regra 1: uma constante de decisao mudou durante a chamada."""


class Inaplicavel(Exception):
    """Regra 1, segunda metade: o limiar/criterio nao se aplica a este caso.

    O subagente REPORTA (motivo) e PARA. Nao adapta, nao escolhe outro."""


PALAVRAS_DE_CONCLUSAO = ("descoberta", "descobrimos", "descoberto", "confirmad",
                         "novo planeta", "planeta novo", "detectamos", "e um planeta")


# ----------------------------------------------------------------- identidade
@dataclass(frozen=True)
class Identidade:
    setor: int
    fonte: str      # "spoc" (2 min, tier 1) ou "qlp" (FFI, 200 s)
    rotulo: str     # rotulo da rodada; entra em TODO caminho de saida

    def __post_init__(self):
        if not isinstance(self.setor, int) or isinstance(self.setor, bool) or self.setor <= 0:
            raise ContratoViolado(f"setor invalido: {self.setor!r}")
        if self.fonte not in ("spoc", "qlp"):
            raise ContratoViolado(f"fonte invalida: {self.fonte!r} (spoc|qlp)")
        if not isinstance(self.rotulo, str) or not self.rotulo.strip():
            raise ContratoViolado("rotulo vazio: caminho sem identidade do conjunto")
        if self.tag not in self.rotulo:
            raise ContratoViolado(f"rotulo {self.rotulo!r} nao carrega a tag do setor {self.tag!r}")
        # so [a-z0-9_]: e o que `cruzar_variaveis._sufixo` e `cruzar_solo._sufixo`
        # preservam - qualquer outra coisa viraria "_" e o nome perderia o rotulo
        if not re.fullmatch(r"[a-z0-9_]+", self.rotulo):
            raise ContratoViolado(f"rotulo so pode ter [a-z0-9_]: {self.rotulo!r}")

    @property
    def tag(self):
        return "s%04d" % self.setor

    @property
    def pasta(self):
        return config.DATA / "orquestra" / self.rotulo

    def artefato(self, nome, ext="parquet"):
        """O UNICO jeito de nomear uma saida: <pasta da rodada>/<nome>__<rotulo>.<ext>."""
        self.pasta.mkdir(parents=True, exist_ok=True)
        return self.pasta / f"{nome}__{self.rotulo}.{ext}"

    def confere_caminho(self, p):
        p = Path(p)
        if self.rotulo not in p.name or self.tag not in p.name:
            raise ContratoViolado(f"artefato sem identidade do conjunto: {p} "
                                  f"(precisa de {self.rotulo!r} e {self.tag!r} no nome)")
        return p


# --------------------------------------------------------------- distribuicoes
def _e_inteiro(x):
    return isinstance(x, int) and not isinstance(x, bool)


@dataclass(frozen=True)
class Mediana:
    valor: float
    n: int
    unidade: str = ""

    def __post_init__(self):
        if not _e_inteiro(self.n) or self.n < 0:
            raise ContratoViolado(f"Mediana sem n inteiro: n={self.n!r}")
        v = self.valor
        if self.n == 0:
            if not (isinstance(v, float) and math.isnan(v)):
                raise ContratoViolado(f"Mediana com n=0 precisa de valor NaN, veio {v!r}")
        else:
            if not isinstance(v, (int, float)) or isinstance(v, bool) or not math.isfinite(float(v)):
                raise ContratoViolado(f"Mediana com n={self.n} e valor nao finito: {v!r}")

    def texto(self):
        if self.n == 0:
            return "NaN (n=0)"
        u = (" " + self.unidade) if self.unidade else ""
        return f"{self.valor:g}{u} (n={self.n})"


def wilson(k, n, z=1.96):
    """IC 95% de Wilson para k de n. n=0 -> (nan, nan)."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z * z / n
    centro = (p + z * z / (2 * n)) / den
    meia = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    lo, hi = max(0.0, centro - meia), min(1.0, centro + meia)
    # nos extremos o intervalo contem k/n por definicao; o arredondamento de
    # (centro - meia) dava 1e-19 > 0 para k = 0 e o validador rejeitava a taxa
    if k == 0:
        lo = 0.0
    if k == n:
        hi = 1.0
    return (lo, hi)


@dataclass(frozen=True)
class Taxa:
    k: int
    n: int
    ic95: tuple

    @classmethod
    def de(cls, k, n):
        k, n = int(k), int(n)
        if n < 0 or not 0 <= k <= n:
            raise ContratoViolado(f"Taxa com k/n invalidos: k={k!r} n={n!r}")
        return cls(k, n, wilson(k, n))

    def __post_init__(self):
        if not _e_inteiro(self.k) or not _e_inteiro(self.n) or self.n < 0 or not 0 <= self.k <= self.n:
            raise ContratoViolado(f"Taxa com k/n invalidos: k={self.k!r} n={self.n!r}")
        ic = self.ic95
        if not (isinstance(ic, tuple) and len(ic) == 2):
            raise ContratoViolado(f"Taxa sem intervalo: ic95={ic!r}")
        lo, hi = ic
        if self.n == 0:
            if not (isinstance(lo, float) and isinstance(hi, float) and math.isnan(lo) and math.isnan(hi)):
                raise ContratoViolado("Taxa com n=0 precisa de IC NaN")
            return
        p = self.k / self.n
        if not (math.isfinite(lo) and math.isfinite(hi) and 0 <= lo <= p <= hi <= 1):
            raise ContratoViolado(f"Taxa com IC que nao contem k/n: {self.k}/{self.n} ic={ic}")

    @property
    def p(self):
        return float("nan") if self.n == 0 else self.k / self.n

    def texto(self):
        if self.n == 0:
            return "0/0 (sem base)"
        lo, hi = self.ic95
        return f"{self.k}/{self.n} = {100*self.p:.1f}% (IC95 {100*lo:.1f}-{100*hi:.1f}%)"


@dataclass(frozen=True)
class Amostra:
    nome: str
    n: int
    selecao: str    # COMO a amostra foi selecionada; obrigatorio

    def __post_init__(self):
        for campo in ("nome", "selecao"):
            v = getattr(self, campo)
            if not isinstance(v, str) or not v.strip():
                raise ContratoViolado(f"Amostra sem {campo}")
        if not _e_inteiro(self.n) or self.n < 0:
            raise ContratoViolado(f"Amostra {self.nome!r} sem n inteiro: {self.n!r}")


@dataclass(frozen=True)
class Comparacao:
    grandeza: str
    a: Amostra
    valor_a: object     # Mediana ou Taxa
    b: Amostra
    valor_b: object

    def __post_init__(self):
        if not isinstance(self.a, Amostra) or not isinstance(self.b, Amostra):
            raise ContratoViolado("Comparacao precisa das DUAS amostras descritas (Amostra)")
        for lado, v, am in (("a", self.valor_a, self.a), ("b", self.valor_b, self.b)):
            if not isinstance(v, (Mediana, Taxa)):
                raise ContratoViolado(f"Comparacao: valor_{lado} precisa ser Mediana ou Taxa")
            if v.n != am.n:
                raise ContratoViolado(f"Comparacao: n da amostra {am.nome!r} ({am.n}) "
                                      f"difere do n do valor ({v.n})")
        if not isinstance(self.grandeza, str) or not self.grandeza.strip():
            raise ContratoViolado("Comparacao sem grandeza")

    def texto(self):
        return (f"{self.grandeza}: {self.a.nome} [{self.a.selecao}] {self.valor_a.texto()}  vs  "
                f"{self.b.nome} [{self.b.selecao}] {self.valor_b.texto()}")


# --------------------------------------------------------------------- classe
class Classe(enum.Enum):
    """Regra 4: nao existe valor acima de CANDIDATO."""
    VETADO = "vetado"
    INDECISO = "indeciso"
    CANDIDATO = "candidato"


# ------------------------------------------------------------------ relatorio
@dataclass
class Relatorio:
    subagente: str
    identidade: Identidade
    contagens: dict = field(default_factory=dict)       # nome -> int
    distribuicoes: dict = field(default_factory=dict)   # nome -> Mediana|Taxa|Comparacao
    artefatos: dict = field(default_factory=dict)       # nome -> Path (com identidade)
    observacoes: list = field(default_factory=list)     # status, NUNCA conclusao
    status: str = "ok"                                  # ok | inaplicavel
    segundos: float = float("nan")

    def validar(self):
        if not isinstance(self.identidade, Identidade):
            raise ContratoViolado("Relatorio sem Identidade")
        if self.status not in ("ok", "inaplicavel"):
            raise ContratoViolado(f"status invalido: {self.status!r}")
        for k, v in self.contagens.items():
            if not _e_inteiro(v):
                raise ContratoViolado(f"contagem {k!r} nao e inteiro: {v!r}")
        for k, v in self.distribuicoes.items():
            if not isinstance(v, (Mediana, Taxa, Comparacao)):
                raise ContratoViolado(f"distribuicao {k!r} de tipo {type(v).__name__}: "
                                      "so Mediana, Taxa ou Comparacao (regra 2)")
        for k, p in self.artefatos.items():
            self.identidade.confere_caminho(p)
            if not Path(p).exists():
                raise ContratoViolado(f"artefato {k!r} relatado mas ausente em disco: {p}")
        for o in self.observacoes:
            if not isinstance(o, str):
                raise ContratoViolado(f"observacao nao textual: {o!r}")
            baixo = o.lower()
            for w in PALAVRAS_DE_CONCLUSAO:
                if w in baixo:
                    raise ContratoViolado(f"observacao com vocabulario de conclusao ({w!r}): {o!r}")
        return self

    def como_dict(self):
        def conv(v):
            if isinstance(v, (Mediana, Taxa, Amostra)):
                return asdict(v)
            if isinstance(v, Comparacao):
                return {"grandeza": v.grandeza, "a": asdict(v.a), "valor_a": conv(v.valor_a),
                        "b": asdict(v.b), "valor_b": conv(v.valor_b)}
            if isinstance(v, Path):
                return str(v)
            return v
        return {"subagente": self.subagente, "identidade": asdict(self.identidade),
                "status": self.status, "segundos": self.segundos,
                "contagens": dict(self.contagens),
                "distribuicoes": {k: conv(v) for k, v in self.distribuicoes.items()},
                "artefatos": {k: str(v) for k, v in self.artefatos.items()},
                "observacoes": list(self.observacoes)}

    def texto(self):
        seg = "" if math.isnan(self.segundos) else f"  ({self.segundos:.0f} s)"
        L = [f"[{self.subagente}] {self.identidade.rotulo}  status={self.status}{seg}"]
        for k, v in self.contagens.items():
            L.append(f"  {k:<58} {v}")
        for k, v in self.distribuicoes.items():
            L.append(f"  {k:<58} {v.texto()}")
        for k, p in self.artefatos.items():
            L.append(f"  -> {k}: {Path(p).name}")
        for o in self.observacoes:
            L.append(f"  * {o}")
        return "\n".join(L)


# ------------------------------------------------------------- regra 1: guarda
class GuardaDeLimiares:
    """Fotografa toda constante MAIUSCULA dos modulos e confere depois.

    Nao tenta adivinhar quais constantes sao limiares: TODAS as maiusculas
    entram, de qualquer tipo que nao seja funcao, classe ou modulo. A primeira
    versao so fotografava (int, float, str, tuple, bool) - e deixava de fora
    `run_vetting.FREE_TESTS` (lista: QUAIS testes rodam) e os dois `CATALOGOS`
    (dict: QUAIS catalogos sao consultados). Um subagente podia tirar um teste
    da bateria e a guarda passava limpa.

    A foto e uma COPIA PROFUNDA, nao a referencia: guardar `m.FREE_TESTS` e
    comparar com `m.FREE_TESTS` depois de um `.remove()` no lugar compara o
    objeto com ele mesmo. Onde deepcopy nao da (objeto exotico), entra o repr.
    """
    IGNORAR = (types.FunctionType, types.BuiltinFunctionType, types.ModuleType, type)

    def __init__(self, modulos):
        self.modulos = list(modulos)
        self.antes = None

    @staticmethod
    def _congelar(v):
        try:
            return copy.deepcopy(v)
        except Exception:  # noqa: BLE001 - objeto que nao copia: compara pelo texto
            return ("<repr>", repr(v))

    def fotografar(self):
        foto = {}
        for m in self.modulos:
            for nome in dir(m):
                if nome.isupper() and not nome.startswith("_"):
                    v = getattr(m, nome)
                    if not isinstance(v, self.IGNORAR):
                        foto[(m.__name__, nome)] = self._congelar(v)
        return foto

    def __enter__(self):
        self.antes = self.fotografar()
        return self

    def __exit__(self, *exc):
        depois = self.fotografar()
        difs = [(k, self.antes.get(k), depois.get(k)) for k in set(self.antes) | set(depois)
                if self.antes.get(k) != depois.get(k)]
        if difs:
            desc = "; ".join(f"{m}.{n}: {a!r} -> {d!r}" for (m, n), a, d in sorted(difs))
            raise LimiarAlterado(f"limiar alterado durante a chamada (regra 1): {desc}")
        return False


# ------------------------------------------------------- regra 5: expectativa
def comparar_com_expectativa(observado, expectativa):
    """Compara campo a campo. Devolve lista de divergencias (campo, esperado, observado).

    Nao distingue direcao: o numero maior que o esperado escala igual ao menor.
    Chaves da expectativa que comecam com "_" sao documentacao. Campo esperado
    e ausente no observado tambem e divergencia - ausencia nao e concordancia.
    """
    difs = []

    def _cmp(prefixo, esp, obs):
        if isinstance(esp, dict):
            for k, v in esp.items():
                if k.startswith("_"):
                    continue
                if not isinstance(obs, dict) or k not in obs:
                    difs.append((prefixo + k, v, "<ausente>"))
                else:
                    _cmp(prefixo + k + ".", v, obs[k])
        elif isinstance(esp, list):
            obs_l = obs if isinstance(obs, list) else []
            if sorted(map(str, esp)) != sorted(map(str, obs_l)):
                difs.append((prefixo.rstrip("."), esp, obs))
        else:
            if esp != obs:
                difs.append((prefixo.rstrip("."), esp, obs))

    _cmp("", expectativa, observado)
    return difs


def contar_campos(expectativa):
    """Quantos campos-folha a expectativa compara (chaves "_" nao contam).

    Vai para o relatorio junto do status: "REPRODUZIU" com 1 campo comparado
    e com 22 sao resultados diferentes, e o status sozinho nao distingue.
    """
    if isinstance(expectativa, dict):
        return sum(contar_campos(v) for k, v in expectativa.items() if not k.startswith("_"))
    return 1


def carregar_expectativa(caminho):
    with open(caminho, encoding="utf-8") as f:
        return json.load(f)
