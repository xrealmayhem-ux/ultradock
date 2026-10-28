"""Utilidades para el almacén raw inmutable: checksums y manifest append-only."""
import csv
import datetime
import hashlib
from pathlib import Path

MANIFEST_COLS = [
    "source", "relative_path", "sha256", "size_bytes",
    "retrieved_at", "origin", "registered_at",
]


def now_utc() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


class RawStore:
    """Manifest append-only de data/raw. Registra archivos nuevos y verifica
    que el contenido de lo ya registrado no cambia nunca."""

    def __init__(self, raw_dir: Path):
        self.raw = raw_dir
        self.manifest = raw_dir / "manifest_sha256.csv"
        self.rows = {}
        if self.manifest.exists():
            with open(self.manifest, newline="", encoding="utf-8") as f:
                for row in csv.DictReader(f):
                    self.rows[row["relative_path"]] = row

    def verify_registered(self):
        bad = []
        for rp, row in self.rows.items():
            p = self.raw / rp
            if not p.exists():
                bad.append((rp, "missing"))
            elif sha256_file(p) != row["sha256"]:
                bad.append((rp, "hash_mismatch"))
        return bad

    def registered(self, path: Path) -> bool:
        return path.relative_to(self.raw).as_posix() in self.rows

    def register(self, path: Path, source: str, retrieved_at: str, origin: str) -> bool:
        """Devuelve True si se registró; False si ya estaba (verificando hash)."""
        rel = path.relative_to(self.raw).as_posix()
        if rel in self.rows:
            if sha256_file(path) != self.rows[rel]["sha256"]:
                raise RuntimeError(f"raw modificado respecto al manifest: {rel}")
            return False
        if not path.exists():
            raise RuntimeError(f"no existe para registrar: {path}")
        self.rows[rel] = {
            "source": source,
            "relative_path": rel,
            "sha256": sha256_file(path),
            "size_bytes": str(path.stat().st_size),
            "retrieved_at": retrieved_at,
            "origin": origin,
            "registered_at": now_utc(),
        }
        return True

    def unregistered_files(self):
        files = []
        for p in sorted(self.raw.rglob("*")):
            if not p.is_file():
                continue
            if p == self.manifest or p.name == "README.md":
                continue
            rel = p.relative_to(self.raw).as_posix()
            if rel not in self.rows:
                files.append(rel)
        return files

    def save(self):
        self.manifest.parent.mkdir(parents=True, exist_ok=True)
        with open(self.manifest, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=MANIFEST_COLS)
            w.writeheader()
            for row in self.rows.values():
                w.writerow(row)
