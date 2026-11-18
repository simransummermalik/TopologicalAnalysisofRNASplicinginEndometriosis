"""Disk budgeting and supervised subprocess groups for local preprocessing."""

from __future__ import annotations

import math
import os
import shlex
import shutil
import signal
import subprocess
from contextlib import ExitStack, contextmanager
from datetime import datetime, timezone
from pathlib import Path

GIB = 1024**3


class RunInterrupted(KeyboardInterrupt):
    def __init__(self, signum):
        self.signum = signum
        super().__init__(f"received signal {signum}")


@contextmanager
def interruption_handlers():
    previous = signal.getsignal(signal.SIGTERM)

    def stop(signum, _frame):
        raise RunInterrupted(signum)

    signal.signal(signal.SIGTERM, stop)
    try:
        yield
    finally:
        signal.signal(signal.SIGTERM, previous)


def add_storage_arguments(parser):
    parser.add_argument(
        "--storage-mode",
        choices=["retain-sam", "stream-bam"],
        default="retain-sam",
        help="retain-sam keeps the existing workflow; stream-bam never writes a SAM file",
    )
    parser.add_argument(
        "--sort-memory-mb",
        type=int,
        default=256,
        help="stream-bam samtools sort memory per thread in MiB (default: 256)",
    )
    parser.add_argument(
        "--min-free-gb",
        type=float,
        default=2.0,
        help="free-space reserve in GiB, added to estimated outputs and rechecked between stages (default: 2)",
    )


def storage_settings(args):
    mode = getattr(args, "storage_mode", "retain-sam")
    memory = getattr(args, "sort_memory_mb", 256)
    reserve = getattr(args, "min_free_gb", 2.0)
    if mode not in {"retain-sam", "stream-bam"}:
        raise ValueError("storage-mode must be retain-sam or stream-bam")
    if not isinstance(memory, int) or memory < 1:
        raise ValueError("sort-memory-mb must be a positive integer")
    if not math.isfinite(reserve) or reserve < 0:
        raise ValueError("min-free-gb must be finite and nonnegative")
    return mode, memory, math.ceil(reserve * GIB)


def disk_budget(args, output, jobs=None):
    mode, _, reserve = storage_settings(args)
    jobs = jobs or [args]
    input_bytes = sum(
        path.stat().st_size
        for job in jobs
        for path in [job.read1, job.read2]
        if path is not None
    )
    # Deliberately conservative file-size heuristics, not guarantees. Original
    # inputs already occupy disk; only newly created outputs/scratch are added.
    multiplier = (
        (3 if args.skip_cleaning else 4)
        if mode == "stream-bam"
        else (12 if args.skip_cleaning else 13)
    )
    estimated = input_bytes * multiplier
    path = Path(output).resolve()
    while not path.exists():
        path = path.parent
    if not path.is_dir():
        raise ValueError(f"output parent is not a directory: {path}")
    free = shutil.disk_usage(path).free
    return {
        "storage_mode": mode,
        "filesystem_path": str(path),
        "input_read_bytes": input_bytes,
        "output_multiplier": multiplier,
        "estimated_additional_bytes": estimated,
        "reserve_bytes": reserve,
        "required_free_bytes": estimated + reserve,
        "available_free_bytes": free,
        "fits": free >= estimated + reserve,
        "basis": "input file bytes times a mode-specific conservative heuristic plus reserve; includes outputs/sort spill, not a guarantee",
    }


def require_budget(budget):
    if not budget["fits"]:
        raise ValueError(
            f"not enough disk space on {budget['filesystem_path']}: "
            f"{budget['available_free_bytes'] / GIB:.2f} GiB free; "
            f"estimated outputs/scratch plus reserve need {budget['required_free_bytes'] / GIB:.2f} GiB "
            f"({budget['storage_mode']}). No processing started. Use a smaller reviewed input or a larger output filesystem."
        )


def print_budget(budget):
    print(
        f"# Disk preflight: {budget['available_free_bytes'] / GIB:.2f} GiB available; "
        f"{budget['required_free_bytes'] / GIB:.2f} GiB estimated including reserve; "
        f"{'fits' if budget['fits'] else 'INSUFFICIENT'} ({budget['storage_mode']})",
        flush=True,
    )


def stop_processes(processes):
    for process in processes:
        try:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGTERM)
            elif process.poll() is None:
                process.terminate()
        except ProcessLookupError:
            pass
    for process in processes:
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            try:
                if os.name == "posix":
                    os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
            except ProcessLookupError:
                pass
            process.wait()


def run_stage(name, commands, output, reserve_bytes, stages, persist):
    """Pipe one or more argv commands; check every exit, terminate peers on error.

    No shell is used. Separate process groups let an interrupted parent stop
    wrappers and their children. Pipe endpoints are closed in the parent so a
    failed consumer cannot leave its producer blocked on a still-open pipe.
    """
    free = shutil.disk_usage(output).free
    stage = {
        "name": name,
        "commands": commands,
        "status": "running",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "free_bytes_before": free,
        "return_codes": [],
        "logs": [],
        "pids": [],
    }
    stages.append(stage)
    persist()
    processes = []
    try:
        if free < reserve_bytes:
            raise ValueError(
                f"disk reserve reached before {name}: {free / GIB:.2f} GiB free; reserve {reserve_bytes / GIB:.2f} GiB"
            )
        print(
            f"[{name}] " + " | ".join(shlex.join(command) for command in commands),
            flush=True,
        )
        with ExitStack() as handles:
            try:
                for number, command in enumerate(commands, 1):
                    log_path = output / f"{len(stages):02d}_{name}_{number}.log"
                    log = handles.enter_context(log_path.open("wb"))
                    previous = processes[-1].stdout if processes else None
                    try:
                        process = subprocess.Popen(
                            command,
                            cwd=output,
                            stdin=previous,
                            stdout=subprocess.PIPE if number < len(commands) else log,
                            stderr=log,
                            start_new_session=os.name == "posix",
                        )
                    finally:
                        if previous is not None:
                            previous.close()
                    processes.append(process)
                    stage["pids"].append(process.pid)
                    stage["logs"].append(log_path.name)
                    persist()
                while any(process.poll() is None for process in processes):
                    if any(process.poll() not in (None, 0) for process in processes):
                        stop_processes(processes)
                        break
                    for process in processes:
                        try:
                            process.wait(timeout=0.1)
                        except subprocess.TimeoutExpired:
                            pass
                codes = [process.wait() for process in processes]
                if any(codes):
                    raise ValueError(
                        f"stage {name} failed (exit codes {codes}); see "
                        + ", ".join(stage["logs"])
                    )
            except BaseException:
                stop_processes(processes)
                raise
            finally:
                for process in processes:
                    if process.stdout is not None:
                        process.stdout.close()
        stage["status"] = "complete"
    except BaseException as error:
        stage.update(
            status="interrupted" if isinstance(error, KeyboardInterrupt) else "failed",
            error=str(error) or type(error).__name__,
        )
        raise
    finally:
        stage["return_codes"] = [process.poll() for process in processes]
        stage["finished_utc"] = datetime.now(timezone.utc).isoformat()
        persist()
