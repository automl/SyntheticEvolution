"""Failure limit, stop file and stop signals; no NePS or Slurm needed."""
import json
import os
import signal
import subprocess
import sys

import pytest

from hpo.pipeline import failures, signals

SETTINGS = dict(min_trials=5, max_rate=0.2, window=50)


def outcomes(pattern):
    """'s' = success, 'f' = failure, in order of finishing."""
    return [(float(index), char == 'f') for index, char in enumerate(pattern)]


def test_limit_waits_for_min_trials():
    assert failures.failure_limit_exceeded(outcomes('ffff'), SETTINGS) is None
    stop = failures.failure_limit_exceeded(outcomes('fffff'), SETTINGS)
    assert (stop['failed'], stop['counted'], stop['rate']) == (5, 5, 1.0)


def test_limit_uses_rolling_window_in_finishing_order():
    # 100 successes, then 11 failures in the last 50: 22% > 20%.
    late_cluster = outcomes('s' * 100 + 's' * 39 + 'f' * 11)
    stop = failures.failure_limit_exceeded(late_cluster, SETTINGS)
    assert (stop['failed'], stop['counted'], stop['finished']) == (11, 50, 150)
    assert stop['reason'] == '22.0% of the last 50 finished trials failed (max_rate 20%)'
    # Exactly 20% does not stop; the order comes from time_end, not list order.
    assert failures.failure_limit_exceeded(outcomes('s' * 40 + 'f' * 10), SETTINGS) is None
    assert failures.failure_limit_exceeded(list(reversed(late_cluster)), SETTINGS) is not None


def test_cumulative_rate_hides_a_late_cluster():
    late_cluster = outcomes('s' * 100 + 's' * 39 + 'f' * 11)
    assert failures.failure_limit_exceeded(late_cluster, dict(SETTINGS, window=None)) is None


def test_max_rate_one_disables_the_limit():
    assert failures.failure_limit_exceeded(outcomes('f' * 20), dict(SETTINGS, max_rate=1.0)) is None


def test_stop_file_is_exclusive(tmp_path):
    assert failures.write_stop(tmp_path, dict(reason='first')) is True
    assert failures.write_stop(tmp_path, dict(reason='second')) is False
    assert json.loads((tmp_path / 'STOPPED.json').read_text())['reason'] == 'first'
    assert [path.name for path in tmp_path.iterdir()] == ['STOPPED.json']


def test_stopped_run_refuses_to_start_and_names_the_rename(tmp_path):
    failures.check_not_stopped(tmp_path)
    failures.write_stop(tmp_path, dict(reason='30.0% of the last 50 finished trials failed', time=100.0))
    (tmp_path / 'STOPPED-1.json').write_text(json.dumps(dict(time=50.0)))
    with pytest.raises(RuntimeError) as error:
        failures.check_not_stopped(tmp_path)
    message = str(error.value)
    assert '30.0% of the last 50' in message
    assert f'mv {tmp_path / "STOPPED.json"} {tmp_path / "STOPPED-2.json"}' in message


def test_restart_time_is_the_newest_renamed_stop(tmp_path):
    assert failures.stop_time(tmp_path) is None
    (tmp_path / 'STOPPED-1.json').write_text(json.dumps(dict(time=50.0)))
    (tmp_path / 'STOPPED-2.json').write_text(json.dumps(dict(time=80.0)))
    # The active stop file does not count until it is renamed.
    (tmp_path / 'STOPPED.json').write_text(json.dumps(dict(time=99.0)))
    assert failures.stop_time(tmp_path) == 80.0
    (tmp_path / 'STOPPED-3.json').write_text('{}')
    with pytest.raises(ValueError, match='STOPPED-3.json'):
        failures.stop_time(tmp_path)


def test_stop_marker():
    assert failures.is_stop('[stopped] failure limit: ...')
    assert not failures.is_stop('[timeout] received SIGUSR1')
    assert not failures.is_stop(None)


@pytest.fixture
def handlers():
    previous = {signum: signal.getsignal(signum) for signum in signals.STOP_SIGNALS}
    signals.reset()
    signals.install()
    yield
    for signum, handler in previous.items():
        signal.signal(signum, handler)
    signals.reset()


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
