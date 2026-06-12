"""
Core functions for the SIMPApy package.

This module contains the main functions for running SOPA on ranked gene data.
"""

# Cap native math/threadpool sizes BEFORE importing numpy/gseapy so that, when SOPA is
# parallelized across samples (one worker process per CPU), each worker's gseapy uses a
# single Rust/rayon thread instead of trying to grab every core. Without this, N workers
# each spawning a full thread pool oversubscribes the CPU and runs >20x slower. The
# explicit `threads=` argument to gp.prerank still overrides this, so the sequential
# path keeps full multithreading.
import os as _os
for _v in ("RAYON_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    _os.environ.setdefault(_v, "1")

import gseapy as gp
import pandas as pd
import numpy as np
import os
import sys
import glob
import multiprocessing as _mp
import contextlib
from typing import Dict, List, Union, Optional, Tuple
import logging # Import the logging library


@contextlib.contextmanager
def _quiet_gseapy():
    """Suppress gseapy's per-call WARNING spam (e.g. "Duplicated values found in
    preranked stats ... order of those genes will be arbitrary") around a prerank call.

    gseapy builds a *fresh* logger per call (name = ``module + id(self)``) with its own
    handler and ``propagate=False``, and resets its level on every construction, so a
    module-level ``setLevel`` cannot catch it and ``verbose=False`` only toggles INFO (the
    warning is a hardcoded ``logger.warning``). ``logging.disable`` gates emission via
    ``isEnabledFor`` regardless of per-instance handlers, so it reliably mutes the noise;
    we save/restore the prior level to leave global logging state untouched. At 20k
    samples this warning is pure noise and measurable I/O overhead.
    """
    prev = logging.root.manager.disable
    logging.disable(logging.WARNING)
    try:
        yield
    finally:
        logging.disable(prev)

# Recycle each worker process after this many samples. Resets the gradual per-sample
# slowdown (~0.85s -> ~5s over a long run) caused by state/memory accumulation in the
# long-lived gseapy/Rust runtime, and bounds worker memory growth.
_MAXTASKS_PER_CHILD = 25

# Per-worker state, populated by `_worker_init` in each spawned worker process. Holds the
# attached shared-memory view of the ranking matrix plus the static inputs.
_WORKER: dict = {}

def _sopa(
    ranking: pd.Series,
    gene_set: Union[Dict, str],
    minisz: int = 3,
    seeder: int = 7,
    threads: int = 8,
    permutation_num: int = 1000,
    **kwargs
) -> pd.DataFrame:
    """
    Run SOPA on a ranked gene list.
    
    Args:
        ranking: A pandas Series with gene names as index and ranking values.
        gene_set: Gene set database in GMT format or a dictionary.
        minisz: Minimum size of gene sets to consider. Default is 3.
        seeder: Random seed for reproducibility. Default is 7.
        threads: Number of threads to use. Default is 8.
        permutation_num: Number of permutations for calculating FDR. Default is 1000.
        **kwargs: Additional arguments passed to gp.prerank().
        
    Returns:
        A pandas DataFrame with SOPA GSEA results sorted by FDR.
    """
    # Pass all arguments to prerank, including any additional arguments
    with _quiet_gseapy():
        pre_res = gp.prerank(
            rnk=ranking,
            gene_sets=gene_set,
            min_size=minisz,
            seed=seeder,
            threads=threads,
            **kwargs
        )
    
    out = []
    for term in list(pre_res.results):
        out.append([
            term,
            pre_res.results[term]['fdr'],
            pre_res.results[term]['es'],
            pre_res.results[term]['nes'],
            pre_res.results[term]['pval'],
            pre_res.results[term]['matched_genes'],
            pre_res.results[term]['gene %'],
            pre_res.results[term]['lead_genes'],
            pre_res.results[term]['tag %']
        ])
    
    out_df = pd.DataFrame(
        out,
        columns=['Term', 'fdr', 'es', 'nes', 'pval', 'matched_genes', 'gene %', 'lead_genes', 'tag %']
    ).sort_values('fdr').reset_index(drop=True)
    
    return out_df


def _output_path(output_dir: str, col: str) -> str:
    """Result file path for a single sample (used by both the sequential and parallel paths)."""
    return os.path.join(output_dir, f"{col}_gsea_results.csv")


def _parallel_supported() -> bool:
    """Whether the parallel fast path may be used on this platform.

    Uses the 'spawn' start method (fresh interpreter per worker) so workers never inherit
    a warm/locked gseapy Rust runtime from the parent — that inheritance deadlocks under
    'fork' if the parent already ran gp.prerank. We honor the project decision to leave
    Windows on the unchanged sequential path.
    """
    return sys.platform != "win32" and "spawn" in _mp.get_all_start_methods()


# Heuristics for the RAM-aware auto worker count. Each spawned worker is a fresh
# Python+gseapy process whose steady-state RSS on this codebase's Hallmark workload
# (9k genes, 50 sets, 1000 perms) measured ~0.5 GB; 0.6 leaves margin. We also keep a
# fixed slice of RAM free so the parent, the shared-memory matrix, and per-call
# allocation spikes never push the box into swap (which froze it once). These are
# calibrated for Hallmark — a much larger gene-set library inflates per-worker RSS and
# would need a bigger _PER_WORKER_GB.
_PER_WORKER_GB = 0.6
_RAM_RESERVE_GB = 2.0


def _usable_cores() -> int:
    """CPU cores this process may actually use. ``sched_getaffinity`` (Linux) respects
    cgroup/taskset/SLURM CPU pinning, so on a constrained HPC allocation we don't count
    host cores we aren't allowed to run on. Falls back to ``os.cpu_count``."""
    try:
        return max(1, len(os.sched_getaffinity(0)))  # Linux only
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 2)


def _available_ram_gb() -> Optional[float]:
    """Best-effort available RAM in GB, accounting for a cgroup memory limit if one is set
    (common under SLURM/containers, where /proc/meminfo reports the whole host, not the
    slice we're allowed). Returns the min of host MemAvailable and (cgroup limit − usage).
    ``None`` when it can't be determined (e.g. non-Linux) so the caller falls back to the
    old fixed cap."""
    host = None
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    host = int(line.split()[1]) / (1024 ** 2)  # kB -> GB
                    break
    except OSError:
        return None

    def _read_int(path):
        try:
            with open(path) as f:
                return int(f.read().strip())
        except (OSError, ValueError):
            return None

    cgroup = None
    # cgroup v2 then v1; "max"/absurdly-large values mean "no limit" -> ignore.
    v2_max, v2_use = _read_int("/sys/fs/cgroup/memory.max"), _read_int("/sys/fs/cgroup/memory.current")
    v1_max = _read_int("/sys/fs/cgroup/memory/memory.limit_in_bytes")
    v1_use = _read_int("/sys/fs/cgroup/memory/memory.usage_in_bytes")
    SANE = 2 ** 53  # ~9 PB; cgroup "unlimited" is typically a near-int64 sentinel
    if v2_max is not None and v2_max < SANE and v2_use is not None:
        cgroup = (v2_max - v2_use) / (1024 ** 3)
    elif v1_max is not None and v1_max < SANE and v1_use is not None:
        cgroup = (v1_max - v1_use) / (1024 ** 3)

    if host is None:
        return cgroup
    return host if cgroup is None else min(host, cgroup)


def _resolve_processes(processes: int, worker_threads: int = 1, matrix_gb: float = 0.0) -> int:
    """Resolve the worker count. ``processes <= 0`` auto-selects, bounded by BOTH cores and
    free RAM: leave 2 cores free (divided by ``worker_threads`` so total OS threads don't
    oversubscribe), and fit within available RAM after reserving headroom + the shared
    matrix. When RAM can't be measured (non-Linux), falls back to the old fixed cap of 8.
    A positive ``processes`` is honored as-is (the user's explicit choice)."""
    if processes is None or processes <= 0:
        core_budget = max(1, (_usable_cores() - 2) // max(1, worker_threads))
        avail = _available_ram_gb()
        if avail is None:
            n = min(core_budget, 8)
            logging.info("SOPA auto: RAM unknown; using %d worker(s) (cores=%d, capped at 8).",
                         n, core_budget)
            return max(1, n)
        ram_budget = int(max(0.0, avail - _RAM_RESERVE_GB - matrix_gb) / _PER_WORKER_GB)
        n = max(1, min(core_budget, ram_budget))
        logging.info(
            "SOPA auto: %d worker(s) [core budget %d, RAM budget %d @ %.1f GB avail, "
            "%.1f GB reserved + %.2f GB shared matrix].",
            n, core_budget, ram_budget, avail, _RAM_RESERVE_GB, matrix_gb,
        )
        return n
    return processes


def _worker_init(shm_name, shape, dtype_str, genes, columns,
                 gene_set, minisz, seeder, output_dir, resume, worker_threads, kwargs):
    """Initializer run once in each spawned worker (and on each recycle).

    Attaches the parent's shared-memory block as a read-only NumPy view — zero-copy, so
    the ranking matrix is shared rather than re-loaded per worker. Stores static inputs.
    """
    from multiprocessing import shared_memory
    shm = shared_memory.SharedMemory(name=shm_name)
    arr = np.ndarray(shape, dtype=np.dtype(dtype_str), buffer=shm.buf)  # rows = samples
    _WORKER.update(
        shm=shm, arr=arr, genes=pd.Index(genes), columns=columns,
        gene_set=gene_set, minisz=minisz, seeder=seeder,
        output_dir=output_dir, resume=resume, worker_threads=worker_threads, kwargs=kwargs,
    )


def _run_one(idx: int):
    """Worker entry point: run SOPA for one sample (by column index) and write its CSV."""
    col = _WORKER["columns"][idx]
    output_file = _output_path(_WORKER["output_dir"], col)
    if _WORKER["resume"] and os.path.exists(output_file):
        return col, "skipped"

    # Build this sample's ranking from the shared-memory row (copy out before mutating).
    ranking = pd.Series(np.array(_WORKER["arr"][idx]), index=_WORKER["genes"])
    ranking = ranking.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)

    call_kw = dict(_WORKER["kwargs"])
    call_kw.pop("threads", None)          # thread count is governed by worker_threads (below)
    call_kw.setdefault("no_plot", True)   # suppress plot files; does not affect result values
    call_kw.setdefault("outdir", None)

    gsea_result = _sopa(
        ranking=ranking,
        gene_set=_WORKER["gene_set"],
        minisz=_WORKER["minisz"],
        seeder=_WORKER["seeder"],
        threads=_WORKER["worker_threads"],
        **call_kw,
    )
    gsea_result.to_csv(output_file, sep=',')
    return col, "done"


def _sopa_parallel(
    ranks: pd.DataFrame,
    gene_set: Union[Dict, str],
    output_dir: str,
    minisz: int,
    seeder: int,
    processes: int,
    resume: bool,
    worker_threads: int,
    kwargs: dict,
) -> None:
    """Parallel SOPA across samples using a 'spawn' pool + one shared-memory ranking matrix.

    The parent copies the ranking matrix (transposed to samples-as-rows) into a single
    shared-memory block once; each spawned worker attaches a zero-copy view and reads only
    its sample's row. This avoids both the fork-after-threads deadlock and the per-worker
    data reload that would exhaust RAM. Workers are recycled every ``_MAXTASKS_PER_CHILD``
    samples to keep per-sample time flat and memory bounded.

    Each worker runs gseapy with ``worker_threads`` threads. On this codebase's profiling
    the per-sample cost is bound by memory bandwidth, so ``worker_threads=1`` with more
    worker processes scales best; ``worker_threads>1`` is exposed for machines/regimes
    where that trade-off differs.
    """
    from multiprocessing import shared_memory

    columns = list(ranks.columns)
    todo = [i for i, c in enumerate(columns)
            if not (resume and os.path.exists(_output_path(output_dir, c)))]
    if not todo:
        logging.info("SOPA: all samples already have results; nothing to do (resume=True).")
        return

    values = ranks.values                       # genes x samples (no copy)
    shape = (len(columns), len(ranks.index))     # samples x genes (rows = samples)
    shm = shared_memory.SharedMemory(create=True, size=values.dtype.itemsize * shape[0] * shape[1])
    try:
        buf = np.ndarray(shape, dtype=values.dtype, buffer=shm.buf)
        np.copyto(buf, values.T)                 # single transpose-copy into shared memory

        ctx = _mp.get_context("spawn")
        n_workers = min(processes, len(todo))
        initargs = (shm.name, shape, values.dtype.str, list(ranks.index), columns,
                    gene_set, minisz, seeder, output_dir, resume, worker_threads, kwargs)
        logging.info(f"SOPA: running {len(todo)} samples across {n_workers} worker process(es).")
        with ctx.Pool(processes=n_workers, initializer=_worker_init, initargs=initargs,
                      maxtasksperchild=_MAXTASKS_PER_CHILD) as pool:
            for _col, _status in pool.imap_unordered(_run_one, todo, chunksize=1):
                pass
    finally:
        shm.close()
        shm.unlink()


def sopa(
    ranks: pd.DataFrame,
    gene_set: Union[Dict, str],
    output_dir: str,
    minisz: int = 3,
    seeder: int = 7,
    processes: int = 1,
    worker_threads: int = 1,
    resume: bool = True,
    **kwargs
) -> None:
    """
    Run SOPA on all samples in the ranking dataframe and save results as CSV files.

    Args:
        ranks: DataFrame with genes as index and samples as columns.
        gene_set: Gene set database in GMT format or a dictionary.
        output_dir: Directory where results will be saved.
        minisz: Minimum size of gene sets to consider. Default is 3.
        seeder: Random seed for reproducibility. Default is 7.
        processes: Number of worker processes for the parallel fast path.
            ``1`` (default) runs the original sequential loop (works everywhere,
            including Windows). ``>1`` parallelizes across samples using a 'spawn'
            pool (Linux/macOS); on unsupported platforms it logs a warning and falls
            back to sequential. ``<=0`` auto-selects a worker count bounded by BOTH
            usable cores (leaving 2 free) and available RAM (~0.6 GB/worker after a 2 GB
            reserve + the shared matrix); it reads cgroup limits and CPU affinity so it
            adapts on a constrained HPC/SLURM allocation, and falls back to
            ``min(n_cores-2, 8)`` where RAM can't be measured (non-Linux). The chosen
            count is logged at INFO. Per-sample numerical output is identical to the
            sequential path (same gseapy, seed, and parameters).
            NOTE: each worker holds a fresh gseapy process, so peak memory grows with
            ``processes``. On this 16-core/15GB dev box ~8-12 workers is the sweet spot
            (W=12 ~3.0h/20k, ~4.8GB free); a 24c/64GB HPC VM auto-scales to ~22. Watch
            ``free -h`` when overriding this manually.
        worker_threads: gseapy ``threads`` per worker in the parallel path. Default ``1``
            (best on memory-bandwidth-bound hardware: prefer more processes over more
            threads each). Raise it only when worker count is RAM-limited and there are
            spare cores. Total OS threads ≈ ``processes * (worker_threads + 1)``; avoid
            oversubscribing the cores. Ignored on the sequential path (use the ``threads``
            kwarg there).
        resume: If True (default), skip any sample whose ``{col}_gsea_results.csv``
            already exists, so an interrupted run can be restarted without redoing work.
        **kwargs: Additional arguments passed through to ``gp.prerank`` (via ``_sopa``).

    Returns:
        None. Results are saved to files in the output directory.
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    use_parallel = processes is not None and processes != 1
    if use_parallel:
        # Estimate the shared ranking matrix size (float64) so the RAM-aware auto count can
        # reserve room for it; cheap shape arithmetic, no array materialization.
        matrix_gb = (ranks.shape[0] * ranks.shape[1] * 8) / (1024 ** 3)
        n_workers = _resolve_processes(processes, worker_threads, matrix_gb)
        if n_workers > 1 and not _parallel_supported():
            logging.warning(
                "Parallel SOPA is not enabled on this platform (Windows stays on the "
                "sequential path). Falling back to processes=1."
            )
            n_workers = 1
        if n_workers > 1:
            _sopa_parallel(ranks, gene_set, output_dir, minisz, seeder, n_workers,
                           resume, worker_threads, kwargs)
            return

    # Sequential path (original behavior; now also honors `resume`).
    for col in ranks.columns:
        output_file = _output_path(output_dir, col)
        if resume and os.path.exists(output_file):
            continue

        # Sort rankings and handle infinities
        ranking = ranks[col].replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)

        # Run sopa with all the parameters
        gsea_result = _sopa(
            ranking=ranking,
            gene_set=gene_set,
            minisz=minisz,
            seeder=seeder,
            **kwargs
        )

        # Save the GSEA results to a CSV file
        gsea_result.to_csv(output_file, sep=',')

        # Clean up
        del gsea_result, ranking

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def load_sopa(directory: str) -> pd.DataFrame:
    """
    Loads and processes sopa results from a directory of CSV files.

    Args:
        directory: The path to the directory containing the SOPA results files.

    Returns:
        A pandas DataFrame with columns: sample_name, term, fdr, pval.
        Returns an empty DataFrame if no valid files are found or processed.
    """

    all_results = []
    required_columns = ['Term', 'fdr', 'es', 'nes', 'matched_genes', 'gene %', 'tag %', 'lead_genes']

    # Use glob to find all CSV files matching the patterns
    file_pattern_tm = os.path.join(directory, "tm*_gsea_results.csv")
    file_pattern_tw = os.path.join(directory, "tw*_gsea_results.csv")
    file_paths = glob.glob(file_pattern_tm) + glob.glob(file_pattern_tw) # Combine lists directly

    if not file_paths:
        logging.warning(f"No files matching 'tm*_gsea_results.csv' or 'tw*_gsea_results.csv' found in directory: {directory}")
        return pd.DataFrame(columns=['sample_name'] + required_columns) # Return empty DataFrame with expected columns

    for file_path in file_paths:
        # Extract sample name from filename
        file_name = os.path.basename(file_path)
        sample_name = file_name.split("_gsea_results")[0]  # Extract tm(n) or tw(n)

        try:
            # Load the CSV file into a DataFrame
            df = pd.read_csv(file_path)

            # Check for and select required columns
            try:
                df_selected = df[required_columns].copy() # Use .copy() to avoid SettingWithCopyWarning
                df_selected['sample_name'] = sample_name
                all_results.append(df_selected)
            except KeyError as e:
                # Handle missing columns after successful parsing
                missing_cols = list(set(required_columns) - set(df.columns))
                logging.error(f"Error: File {file_path} is missing required columns: {missing_cols}. Skipping.")
                continue # Skip this file

        except pd.errors.ParserError:
            # Handle files that cannot be parsed as CSV
            logging.error(f"Error: Could not parse {file_path} as a CSV file. Skipping.")
            continue
        except FileNotFoundError:
            # Handle case where file disappears between glob and read_csv (unlikely but possible)
            logging.error(f"Error: File {file_path} not found. Skipping.")
            continue
        except Exception as e:
             # Catch any other unexpected errors during file processing
             logging.error(f"An unexpected error occurred processing {file_path}: {e}. Skipping.")
             continue
        # handle when columns are not correctly formatted (e.g. missing 'Term' or 'fdr' columns) after parsing
        except KeyError as e:
            missing_cols = list(set(required_columns) - set(df.columns))
            logging.error(f"Error: File {file_path} is missing required columns: {missing_cols}. Skipping.")
            continue # Skip this file

    # Concatenate all results into a single DataFrame
    if not all_results:
        logging.warning(f"No valid SOPA result files were processed successfully in directory: {directory}")
        return pd.DataFrame(columns=['sample_name'] + required_columns) # Return empty DataFrame if no files were valid

    final_df = pd.concat(all_results, ignore_index=True)

    return final_df
