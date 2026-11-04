"""Stop signals from Slurm (time limit, scancel) during and between trials."""
import os
import signal
import subprocess
import sys

import pytest

from hpo.pipeline import signals


@pytest.fixture
def handlers(monkeypatch):
    previous = {signum: signal.getsignal(signum) for signum in signals.STOP_SIGNALS}
    monkeypatch.setattr(signals, '_requested', None)
    monkeypatch.setattr(signals, '_in_trial', False)
    signals.install()
    yield
    for signum, handler in previous.items():
        signal.signal(signum, handler)


def test_signal_between_trials_is_only_recorded(handlers):
    os.kill(os.getpid(), signal.SIGUSR1)
    assert signals.requested() == 'SIGUSR1'


def test_signal_during_trial_kills_the_running_tool(handlers, tmp_path):
    pid_file = tmp_path / 'pid'
    child = f'import os, time; open({str(pid_file)!r}, "w").write(str(os.getpid())); time.sleep(60)'
    signal.signal(signal.SIGALRM, lambda *args: os.kill(os.getpid(), signal.SIGUSR1))
    try:
        signal.setitimer(signal.ITIMER_REAL, 0.5)
        with pytest.raises(signals.ControllerStop, match='SIGUSR1'):
            with signals.trial_running():
                subprocess.run([sys.executable, '-c', child], timeout=30)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, signal.SIG_DFL)
    # subprocess.run kills and reaps its child when interrupted.
    with pytest.raises(ProcessLookupError):
        os.kill(int(pid_file.read_text()), 0)


def test_second_signal_does_not_interrupt_the_stop(handlers):
    with pytest.raises(signals.ControllerStop):
        with signals.trial_running():
            os.kill(os.getpid(), signal.SIGUSR1)
    with signals.trial_running():
        os.kill(os.getpid(), signal.SIGTERM)  # must not raise again
    assert signals.requested() == 'SIGUSR1'
