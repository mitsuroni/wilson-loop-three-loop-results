"""Per-solve CAS allocation, charged only around owned process execution.

The measured wall interval includes process launch and mandatory tree cleanup,
but excludes script generation, parsing, cached moments and symbolic assembly.
An exception still charges this interval before propagating to the solver.
"""
import math
import time


class CASBudgetExhausted(TimeoutError):
    def __init__(self, code):
        self.code = code
        super().__init__('CAS allocation exhausted' if code == 'cas_budget'
                         else 'calculation wall-time allocation exhausted')


class CASBudget:
    def __init__(self, seconds, cas_seconds, *, provenance=None):
        if any(isinstance(x, bool) or not math.isfinite(x) for x in (seconds, cas_seconds)):
            raise ValueError('finite resource allocations required')
        if seconds <= 0 or cas_seconds < 0:
            raise ValueError('positive seconds and nonnegative cas_seconds required')
        self.deadline = time.monotonic() + seconds
        self.cas_seconds = cas_seconds
        self.elapsed = 0.0
        self.runs = []
        self.provenance = provenance
        self._publish()

    @property
    def seconds_remaining(self):
        return max(0.0, self.deadline - time.monotonic())

    @property
    def cas_remaining(self):
        return max(0.0, self.cas_seconds - self.elapsed)

    def _publish(self):
        if self.provenance is not None:
            self.provenance.update(cas_elapsed_seconds=self.elapsed,
                cas_process_runs=list(self.runs),
                cas_accounting='owned-process wall interval including launch/cleanup')

    def run(self, runner, *args, timeout_seconds, **kwargs):
        remaining = self.seconds_remaining
        if remaining <= 0:
            raise CASBudgetExhausted('calculation_time_budget')
        if self.cas_remaining <= 0:
            raise CASBudgetExhausted('cas_budget')
        allowed = min(timeout_seconds, remaining, self.cas_remaining)
        if allowed <= 0:
            raise CASBudgetExhausted('calculation_time_budget')
        started = time.monotonic()
        record = {'timeout_seconds': allowed}
        try:
            answer = runner(*args, timeout_seconds=allowed, **kwargs)
            record.update({key: answer[key] for key in ('returncode', 'stop_reason') if key in answer})
            return answer
        except BaseException as exc:
            record['exception_type'] = type(exc).__name__
            raise
        finally:
            elapsed = max(0.0, time.monotonic() - started)
            self.elapsed += elapsed
            record['elapsed_seconds'] = elapsed
            self.runs.append(record)
            self._publish()
