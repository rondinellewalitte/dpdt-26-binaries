"""Download resiliente das light curves, por HTTP direto.

As URLs vem do manifesto do setor (bulk_manifest.py), entao nao ha uma consulta
de produto por alvo - so o GET do arquivo. No maximo 4 conexoes simultaneas, com
backoff exponencial, respeitando o rate limit do MAST.

HTTP 200 nao prova FITS valido: um proxy ou uma pagina de erro do servidor
devolve 200 com HTML e o arquivo corrompido viraria um alvo silenciosamente
perdido. Por isso todo download e verificado pelos primeiros bytes (`SIMPLE  =`,
assinatura obrigatoria de FITS) e pelo Content-Length quando anunciado.
"""
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import config

UA = {"User-Agent": "nasa_discover/0.1 (astronomia amadora; contato via github)"}
FITS_MAGIC = b"SIMPLE  ="


class DiskCapExceeded(RuntimeError):
    pass


def check_disk_cap(cache_dir=None, cap_gb=None):
    cache = Path(cache_dir or config.CACHE)
    cap = (cap_gb or config.DISK_CAP_GB) * 1024 ** 3
    used = sum(f.stat().st_size for f in cache.rglob("*") if f.is_file())
    if used > cap:
        raise DiskCapExceeded(
            f"cache em {used/1024**3:.2f} GB, acima do teto de {cap/1024**3:.2f} GB")
    return used


def download_one(url, dest, timeout=180, max_retries=None):
    """Baixa e valida um FITS. Devolve (ok, nbytes, erro)."""
    max_retries = max_retries or config.MAST_MAX_RETRIES
    dest = Path(dest)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last = None

    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                announced = r.headers.get("Content-Length")
                data = r.read()

            if len(data) < len(FITS_MAGIC) or not data.startswith(FITS_MAGIC):
                head = data[:80].decode("utf-8", "replace").strip()
                raise ValueError(f"resposta nao e FITS (comeca com {head!r})")
            if announced and int(announced) != len(data):
                raise ValueError(f"tamanho {len(data)} != Content-Length {announced}")

            tmp.write_bytes(data)
            tmp.replace(dest)
            return True, len(data), None

        except Exception as e:
            last = f"{type(e).__name__}: {e}"
            if attempt < max_retries - 1:
                time.sleep(min(config.MAST_BACKOFF_BASE ** attempt, config.MAST_BACKOFF_CAP))
        finally:
            if tmp.exists():
                try:
                    tmp.unlink()
                except OSError:
                    pass
    return False, 0, last


def download_batch(items, cache_dir=None, n_threads=None, on_done=None):
    """items: [(tic, filename, url)]. Devolve {tic: (ok, nbytes, erro, segundos)}.

    O ThreadPool fica no processo principal de proposito: download e I/O, o GIL
    nao atrapalha, e assim o processo principal segue sendo o unico escritor do
    state.db.
    """
    cache = Path(cache_dir or config.CACHE)
    cache.mkdir(parents=True, exist_ok=True)
    n_threads = n_threads or config.N_DOWNLOAD_THREADS

    def work(item):
        tic, fname, url = item
        t0 = time.perf_counter()
        dest = cache / fname
        # FITS ja em disco e valido nao se baixa de novo: numa retomada isso
        # evita re-transferir o lote inteiro que ja tinha vindo antes da queda.
        if dest.exists():
            try:
                with open(dest, "rb") as fh:
                    if fh.read(len(FITS_MAGIC)) == FITS_MAGIC:
                        return tic, (True, dest.stat().st_size, None, 0.0)
            except OSError:
                pass
        ok, nb, err = download_one(url, dest)
        return tic, (ok, nb, err, time.perf_counter() - t0)

    out = {}
    with ThreadPoolExecutor(max_workers=n_threads) as ex:
        futs = {ex.submit(work, it): it for it in items}
        for fut in as_completed(futs):
            tic, res = fut.result()
            out[tic] = res
            if on_done:
                on_done(tic, res)
    return out
