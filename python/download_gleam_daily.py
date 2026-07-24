"""
Download GLEAM v4.3a daily AET NetCDF files via SFTP from hydras.ugent.be.

Usage:
    python scripts/download_gleam_daily.py

Before running:
    pip install paramiko

The script connects to the GLEAM SFTP server, lists available yearly files
in the daily AET directory, and downloads any that are not already present
locally. Files are saved to data/GLEAM/AET_daily/.

Credentials:
    Set GLEAM_USER and GLEAM_PASS environment variables, or enter them
    interactively at the prompt.
"""

import os
import sys
import socket
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import paramiko
except ImportError:
    sys.exit("paramiko is required. Install with: pip install paramiko")


def load_config(filepath):
    '''Load configuration file into dict'''
    df = pd.read_csv(filepath, header=None, index_col=False)
    config = {}
    for ii in np.arange(len(df)):
        string = df.iloc[ii, 0].replace(" ", "")
        varname = string.rpartition('=')[0]
        varcontents = string.rpartition('=')[2]
        try:
            varcontents = float(varcontents)
        except:
            pass
        config[varname] = varcontents
    return config


# key text file containing (1) proxy key, (2) CDS API key,
# (3) working directory of the data
config = load_config("proxy_config.cfg")

# Proxy settings
proxy = config['proxyKey']
os.environ['http_proxy'] = proxy
os.environ['HTTP_PROXY'] = proxy
os.environ['https_proxy'] = proxy
os.environ['HTTPS_PROXY'] = proxy

# --- Configuration ------------------------------------------------
SFTP_HOST = "hydras.ugent.be"
SFTP_PORT = 2225  # GLEAM SFTP uses a non-standard port

# Proxy configuration ----------------------------------------------
# Corporate proxies block direct SSH/SFTP, causing WinError 10060.
# The proxy is read from proxy_config.cfg (proxyKey) above and used to
# tunnel the SSH connection via HTTP CONNECT.
PROXY_URL = proxy

# Remote directory structure on the GLEAM server:
#   /data/v4.3a/daily/<YYYY>/E_<YYYY>_GLEAM_v4.3a.nc
# We only want the E (AET) variable file from each year folder.
REMOTE_BASE = "/data/v4.3a/daily"

# Local destination
BASE_DIR = Path("D:/tilloal/Documents/01_Projects/RegimeShifts")
LOCAL_DIR = BASE_DIR / "data" / "GLEAM" / "AET_daily"

# --- Credentials (hardcoded) --------------------------------------
user = "gleamuser"
password = "GLEAM4#h-cel_111"

# --- Connect and download -----------------------------------------
LOCAL_DIR.mkdir(parents=True, exist_ok=True)


def open_socket():
    """Open a socket to the SFTP host, tunneling through PROXY_URL if set."""
    if PROXY_URL:
        try:
            from python_socks.sync import Proxy
        except ImportError:
            sys.exit(
                "A proxy is configured but python-socks is not installed.\n"
                "Install with: pip install python-socks"
            )
        print(f"Tunneling through proxy: {PROXY_URL}")
        proxy = Proxy.from_url(PROXY_URL)
        # HTTP CONNECT tunnel to the SFTP host:port
        return proxy.connect(dest_host=SFTP_HOST, dest_port=SFTP_PORT, timeout=30)
    # Direct connection
    return socket.create_connection((SFTP_HOST, SFTP_PORT), timeout=30)


sock = open_socket()
transport = paramiko.Transport(sock)
transport.connect(username=user, password=password)
sftp = paramiko.SFTPClient.from_transport(transport)

print(f"Connected to {SFTP_HOST}")
print(f"Remote base: {REMOTE_BASE}")

# List year folders under /data/v4.3a/daily/
try:
    entries = sftp.listdir(REMOTE_BASE)
except FileNotFoundError:
    sys.exit(f"Remote directory not found: {REMOTE_BASE}\n"
             "Check your account permissions or adjust REMOTE_BASE.")

year_dirs = sorted([e for e in entries if e.isdigit() and len(e) == 4])
if not year_dirs:
    sys.exit(f"No year folders found in {REMOTE_BASE}")

print(f"Found {len(year_dirs)} year folder(s): {year_dirs[0]} .. {year_dirs[-1]}")

# In each year folder, grab only the E (AET) file:
#   E_<YYYY>_GLEAM_v4.3a.nc
remote_files = []
for yr in year_dirs:
    yr_path = f"{REMOTE_BASE}/{yr}"
    for fname in sftp.listdir(yr_path):
        if fname.startswith("E_") and fname.endswith(".nc"):
            remote_files.append((f"{yr_path}/{fname}", fname))

print(f"Found {len(remote_files)} E (AET) file(s) to download")

# Download files not already present locally
downloaded = 0
skipped = 0
for remote_path, fname in sorted(remote_files):
    local_path = LOCAL_DIR / fname
    if local_path.exists():
        # Optionally check size to detect partial downloads
        remote_size = sftp.stat(remote_path).st_size
        if local_path.stat().st_size == remote_size:
            skipped += 1
            continue
    print(f"  Downloading: {fname} ...")
    sftp.get(remote_path, str(local_path))
    downloaded += 1

sftp.close()
transport.close()

print(f"\nDone. Downloaded: {downloaded} | Already present: {skipped}")
print(f"Local directory: {LOCAL_DIR}")
