r"""Password-protected dataset links for the static viewer.

A `<id>.link.txt` normally holds a plain URL (like the Nerve page's CSV links).
With this tool it can instead hold that URL encrypted with a shared password:

    enc:v1:<salt b64>:<iv b64>:<ciphertext b64>

The page asks for the password once, derives the key in the browser
(PBKDF2-SHA256, 200k rounds) and decrypts the URL with AES-GCM. Combined with
an unguessable folder name on the blob (container listing is off), the data is
unreachable without the password, and there is no server to run.

    python links.py secret                                   -> prints a new random folder name
    python links.py encrypt --password <pw> --url https://.../SpatialViewer/<secret>/slide1/ --out slide1.link.txt
    python links.py decrypt --password <pw> --file slide1.link.txt

Needs the `cryptography` package (in requirements.txt).
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import secrets
import sys

ROUNDS = 200_000


def encrypt(password: str, url: str) -> str:
    from cryptography.hazmat.primitives.ciphers.aead import AESGCM
    salt = secrets.token_bytes(16)
    iv = secrets.token_bytes(12)
    key = hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), salt, ROUNDS, dklen=32)
    ct = AESGCM(key).encrypt(iv, url.encode("utf-8"), None)
    b = lambda x: base64.b64encode(x).decode()
    return f"enc:v1:{b(salt)}:{b(iv)}:{b(ct)}"


def decrypt(password: str, payload: str) -> str:
    from cryptography.hazmat.primitives.ciphers.aead import AESGCM
    tag, ver, salt, iv, ct = payload.strip().split(":")
    assert tag == "enc" and ver == "v1", "not an enc:v1 link"
    d = base64.b64decode
    key = hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), d(salt), ROUNDS, dklen=32)
    return AESGCM(key).decrypt(d(iv), d(ct), None).decode("utf-8")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("secret", help="print a new random folder name for the blob")
    e = sub.add_parser("encrypt")
    e.add_argument("--password", required=True)
    e.add_argument("--url", required=True)
    e.add_argument("--out", help="write to this file (default: print)")
    dd = sub.add_parser("decrypt")
    dd.add_argument("--password", required=True)
    dd.add_argument("--file", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "secret":
        print(secrets.token_urlsafe(18).replace("-", "x").replace("_", "y"))
    elif a.cmd == "encrypt":
        s = encrypt(a.password, a.url)
        if a.out:
            open(a.out, "w", encoding="utf-8", newline="\n").write(s + "\n")
            print(f"wrote {a.out}")
        else:
            print(s)
    else:
        print(decrypt(a.password, open(a.file, encoding="utf-8").read()))


if __name__ == "__main__":
    sys.exit(main())
