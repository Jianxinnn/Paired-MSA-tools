from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Optional

import requests
from requests.auth import HTTPBasicAuth


@dataclass
class RemoteClient:
    host_url: str
    user_agent: str = "paired-msa-tools/0.1"
    email: str = ""
    auth_user: Optional[str] = None
    auth_pass: Optional[str] = None

    def _auth(self) -> Optional[HTTPBasicAuth]:
        if self.auth_user is None or self.auth_pass is None:
            return None
        return HTTPBasicAuth(self.auth_user, self.auth_pass)

    def submit(self, query: str, endpoint: str, mode: str) -> dict:
        headers = {"User-Agent": self.user_agent}
        payload = {"q": query, "mode": mode, "email": self.email}
        while True:
            try:
                res = requests.post(
                    f"{self.host_url}/{endpoint}",
                    data=payload,
                    timeout=6.02,
                    headers=headers,
                    auth=self._auth(),
                )
            except requests.exceptions.Timeout:
                continue
            except Exception:
                time.sleep(5)
                continue
            break
        try:
            out = res.json()
        except Exception:
            out = {"status": "ERROR", "text": res.text}
        return out

    def status(self, job_id: str) -> dict:
        headers = {"User-Agent": self.user_agent}
        while True:
            try:
                res = requests.get(
                    f"{self.host_url}/ticket/{job_id}",
                    timeout=6.02,
                    headers=headers,
                    auth=self._auth(),
                )
            except requests.exceptions.Timeout:
                continue
            except Exception:
                time.sleep(5)
                continue
            break
        try:
            out = res.json()
        except Exception:
            out = {"status": "ERROR"}
        return out

    def download(self, job_id: str) -> bytes:
        headers = {"User-Agent": self.user_agent}
        while True:
            try:
                res = requests.get(
                    f"{self.host_url}/result/download/{job_id}",
                    timeout=6.02,
                    headers=headers,
                    auth=self._auth(),
                )
            except requests.exceptions.Timeout:
                continue
            except Exception:
                time.sleep(5)
                continue
            break
        res.raise_for_status()
        return res.content

