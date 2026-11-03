"""Immutable byte objects, local single-writer lock, transactional receipt index."""
import contextlib
import json
import os
import re
import sqlite3
import tempfile
import uuid
from pathlib import Path

from .contracts import ROOT, canonical, sha, require, validate, RECEIPT, InterfaceError


class Store:
    def __init__(self, path=None):
        self.path = Path(path or ROOT / '.wilson/proof-store').resolve()
        require(self.path.is_relative_to(ROOT.resolve()), 'store must be in the workspace')
        self.path.mkdir(parents=True, exist_ok=True)
        (self.path / 'objects').mkdir(exist_ok=True)
        (self.path / 'runs').mkdir(exist_ok=True)
        self.db = sqlite3.connect(self.path / 'index.sqlite3', isolation_level=None, timeout=5)
        self.db.row_factory = sqlite3.Row
        self.db.execute('PRAGMA foreign_keys=ON')
        self.db.execute('PRAGMA journal_mode=WAL')
        self.db.executescript('''
        CREATE TABLE IF NOT EXISTS runs(
          id TEXT PRIMARY KEY, request TEXT NOT NULL, status TEXT NOT NULL
          CHECK(status IN ('running','succeeded','failed','aborted')), report TEXT);
        CREATE TABLE IF NOT EXISTS receipts(
          id TEXT PRIMARY KEY, goal TEXT NOT NULL, request TEXT NOT NULL, profile TEXT NOT NULL,
          run TEXT NOT NULL REFERENCES runs(id));
        CREATE TABLE IF NOT EXISTS edges(
          child TEXT NOT NULL REFERENCES receipts(id), parent TEXT NOT NULL REFERENCES receipts(id),
          kind TEXT NOT NULL CHECK(kind IN ('proof','planning')), PRIMARY KEY(child,parent,kind));
        CREATE INDEX IF NOT EXISTS request_profile ON receipts(request,profile);
        ''')

    def close(self):
        self.db.close()

    @contextlib.contextmanager
    def writer(self):
        path = self.path / 'writer.lock'
        with path.open('a+b') as lock:
            lock.seek(0, 2)
            if lock.tell() == 0:
                lock.write(b'0')
                lock.flush()
            lock.seek(0)
            try:
                if os.name == 'nt':
                    import msvcrt
                    msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
                else:
                    import fcntl
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError as exc:
                raise InterfaceError('another local writer is active') from exc
            try:
                yield
            finally:
                lock.seek(0)
                if os.name == 'nt':
                    msvcrt.locking(lock.fileno(), msvcrt.LK_UNLCK, 1)
                else:
                    fcntl.flock(lock, fcntl.LOCK_UN)

    def object_path(self, key):
        require(isinstance(key, str) and re.fullmatch('[0-9a-f]{64}', key), 'invalid object digest')
        return self.path / 'objects' / key[:2] / key[2:]

    def put(self, data):
        key = sha(data)
        path = self.object_path(key)
        path.parent.mkdir(exist_ok=True)
        if path.exists():
            require(self.read(key) == data, 'object collision')
            return key
        fd, tmp = tempfile.mkstemp(dir=path.parent, prefix='.tmp-')
        try:
            with os.fdopen(fd, 'wb') as f:
                f.write(data)
                f.flush()
                os.fsync(f.fileno())
            try:
                os.link(tmp, path)  # Atomic create, never overwrite an existing object.
            except FileExistsError:
                require(self.read(key) == data, 'object collision')
        finally:
            Path(tmp).unlink(missing_ok=True)
        return key

    def put_json(self, value):
        return self.put(canonical(value))

    def read(self, key):
        try:
            data = self.object_path(key).read_bytes()
        except OSError as exc:
            raise InterfaceError('missing object: ' + key) from exc
        require(sha(data) == key, 'corrupt object: ' + key)
        return data

    def json(self, key):
        return json.loads(self.read(key))

    def start(self, request):
        run = uuid.uuid4().hex
        self.db.execute('INSERT INTO runs VALUES(?,?,?,NULL)', (run, request, 'running'))
        folder = self.path / 'runs' / run
        folder.mkdir()
        return run, folder

    def fail(self, run, report):
        self.db.execute('UPDATE runs SET status=?, report=? WHERE id=? AND status=?', ('failed', report, run, 'running'))

    def recover(self):
        with self.writer():
            # The OS lock is released on process death, so there is no live writer here.
            return self.db.execute("UPDATE runs SET status='aborted' WHERE status='running'").rowcount

    def commit_checked(self, run, report, receipts, fault=None):
        """Called only by the replay service; object/import APIs cannot grant acceptance."""
        row = self.db.execute('SELECT status FROM runs WHERE id=?', (run,)).fetchone()
        require(row and row['status'] == 'running', 'invalid run transition')
        require(self.json(report)['passed'], 'cannot commit failed replay')
        keys = []
        for r in receipts:
            validate(r, RECEIPT)
            require(r['run'] == run and r['report'] == report, 'receipt run binding')
            for field in ('goal', 'request', 'profile', 'source', 'report'):
                self.read(r[field])
            keys.append(self.put_json(r))
        if fault:
            fault('objects_written')
        self.db.execute('BEGIN IMMEDIATE')
        try:
            for key, r in zip(keys, receipts):
                self.db.execute('INSERT INTO receipts VALUES(?,?,?,?,?)', (key, r['goal'], r['request'], r['profile'], run))
                for kind, field in [('proof', 'proof_dependencies'), ('planning', 'planning_dependencies')]:
                    for parent in r[field]:
                        self.db.execute('INSERT INTO edges VALUES(?,?,?)', (key, parent, kind))
            if fault:
                fault('index_written')
            self.db.execute('UPDATE runs SET status=?,report=? WHERE id=?', ('succeeded', report, run))
            self.db.execute('COMMIT')
        except BaseException:
            self.db.execute('ROLLBACK')
            raise
        return keys

    def receipt(self, key):
        row = self.db.execute('SELECT * FROM receipts WHERE id=?', (key,)).fetchone()
        require(row is not None, 'receipt not admitted in this index: ' + key)
        r = self.json(key)
        validate(r, RECEIPT)
        require(all(row[f] == r[f] for f in ('goal', 'request', 'profile', 'run')), 'index binding mismatch')
        run = self.db.execute('SELECT * FROM runs WHERE id=?', (r['run'],)).fetchone()
        require(run and run['status'] == 'succeeded' and run['report'] == r['report'], 'receipt has no successful run')
        return r

    def state(self, key, profile, stack=None):
        stack = set() if stack is None else stack
        if key in stack:
            return dict(state='damaged', reason='dependency cycle')
        try:
            r = self.receipt(key)
            for field in ('goal', 'request', 'profile', 'source', 'report'):
                self.read(r[field])
            snap = self.json(r['profile'])
            for blob in snap['source_objects'].values():
                self.read(blob)
            if r['profile'] != profile:
                return dict(state='stale', reason='profile changed')
            for parent in r['proof_dependencies'] + r['planning_dependencies']:
                result = self.state(parent, profile, stack | {key})
                if result['state'] != 'current':
                    return dict(state=result['state'], reason='dependency ' + parent + ': ' + result['reason'])
            return dict(state='current', reason='checked under current profile')
        except (InterfaceError, ValueError, KeyError) as exc:
            return dict(state='damaged', reason=str(exc))

    def export(self):
        return dict(format='wilson-store-export/0.1', admission='replay_required',
                    receipts=[r['id'] for r in self.db.execute('SELECT id FROM receipts ORDER BY rowid')],
                    runs=[dict(r) for r in self.db.execute('SELECT * FROM runs ORDER BY rowid')])
