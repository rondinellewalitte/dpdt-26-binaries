"""Orquestra: um orquestrador e cinco subagentes para busca em escala.

Os "subagentes" sao MODULOS PYTHON com contrato tipado, nao agentes de modelo de
linguagem. Quatro das cinco regras sao garantias ("nunca ajusta limiar", "sempre
devolve distribuicao", "caminho sempre carrega identidade", "nunca conclui
descoberta") e garantia que so se espera nao e garantia: em codigo, cada uma e
uma assercao que dispara. O orquestrador e um driver Python que encadeia os
cinco e compara o resultado com uma expectativa escrita ANTES da rodada.

    contrato.py       o que um subagente PODE devolver; o validador que rejeita
    auditor.py        roda os testes estruturais antes e depois; para tudo se falha
    dimensionador.py  mede alvo-set, disco, tempo e sensibilidade; devolve numeros
    buscador.py       roda (ou le) a busca; nao decide nada
    vetador.py        bateria completa; veredito POR TESTE, nunca so o agregado
    cruzador.py       TESS + solo + TOI/confirmados, por ID e por posicao
    orquestrador.py   escolhe setor por criterio explicito, encadeia, consolida
"""
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1]
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
