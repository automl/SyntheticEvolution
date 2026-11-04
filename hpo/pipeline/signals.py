"""Turn Slurm's stop signals into a clean stop of the running trial.

Slurm sends SIGUSR1 shortly before the time limit (``--signal=B:USR1@600``)
and SIGTERM on ``scancel``. While a trial runs, the handler raises
ControllerStop, which also kills a running AF3/DSSR child process
(``subprocess.run`` kills its child when interrupted). Between trials it only
records the request, so it never interrupts NePS while NePS updates its state.
"""

import logging
import signal
from contextlib import contextmanager
from typing import Iterator, Optional

logger = logging.getLogger(__name__)

STOP_SIGNALS = (signal.SIGUSR1, signal.SIGTERM)
EXIT_CODE = 4

_requested: Optional[str] = None
_in_trial = False


class ControllerStop(BaseException):
    """Raised inside a running trial when Slurm asks the controller to stop.

    A BaseException, so neither the pipeline nor NePS handle it as an
    ordinary trial failure and continue with the next trial.
    """


def _handle(signum, frame) -> None:
    global _requested
    name = signal.Signals(signum).name
    first = _requested is None
    if first:
        _requested = name
    logger.warning('Received %s (Slurm time limit or cancellation); stopping', name)
    # Raise only once: a second signal must not interrupt the clean stop.
    if first and _in_trial:
        raise ControllerStop(f'received {name} (Slurm time limit or cancellation)')


def install() -> None:
    """Handle the stop signals; call once at controller start."""
    for signum in STOP_SIGNALS:
        signal.signal(signum, _handle)


def requested() -> Optional[str]:
    """Name of the stop signal received so far, if any."""
    return _requested


@contextmanager
def trial_running() -> Iterator[None]:
    """Mark the code in which a stop signal may interrupt the trial."""
    global _in_trial
    _in_trial = True
    try:
        yield
    finally:
        _in_trial = False
