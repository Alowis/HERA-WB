"""
Download GlobSnow v3.0 monthly SWE NetCDF4 files.
Source: https://www.globsnow.info/swe/archive_v3.0/L3B_monthly_SWE/NetCDF4/

Requirements:
    pip install requests beautifulsoup4
"""
import os, sys, glob, time, pdb
import requests
from bs4 import BeautifulSoup
from pathlib import Path
import pandas as pd
import numpy as np
import re

BASE_URL = "https://www.globsnow.info/swe/archive_v3.0/L3A_daily_SWE/NetCDF4/"
OUT_DIR = Path("./data/globsnow_swe_daily")
OUT_DIR.mkdir(exist_ok=True)
YEARS     = 42    

def load_config(filepath):
    '''Load configuration file into dict'''
    
    df = pd.read_csv(filepath,header=None,index_col=False)
    config = {}
    for ii in np.arange(len(df)): 
        string = df.iloc[ii,0].replace(" ","")
        varname = string.rpartition('=')[0]
        varcontents = string.rpartition('=')[2]
        try:
            varcontents = float(varcontents)
        except:
            pass
        config[varname] = varcontents
    return config
    

config = load_config("proxy_config.cfg")
# key text file containing (1) proxy key, (2) CDS API key, (3) working directory of the data

# Proxy settings
proxy = config ['proxyKey']
os.environ['http_proxy'] = proxy 
os.environ['HTTP_PROXY'] = proxy
os.environ['https_proxy'] = proxy
os.environ['HTTPS_PROXY'] = proxy
print("Fetching file list...")
resp = requests.get(BASE_URL, timeout=30)
resp.raise_for_status()
soup = BeautifulSoup(resp.text, "html.parser")
all_files = sorted(a["href"] for a in soup.find_all("a", href=True) if a["href"].endswith(".nc"))

if not all_files:
    print("No .nc files found — the server may have changed its structure.")
    raise SystemExit(1)
 
# ── Extract year from filename and filter ──────────────────────────────────────
def extract_year(fname):
    m = re.search(r"(\d{4})", fname)
    return int(m.group(1)) if m else None
 
years_available = sorted(set(extract_year(f) for f in all_files if extract_year(f)))
 
if YEARS is not None:
    selected_years = set(years_available[:YEARS])
    files = [f for f in all_files if extract_year(f) in selected_years]
    #print(f"Years in archive : {years_available[0]}–{years_available[-1]}")
    print(f"Downloading first {YEARS} years: {sorted(selected_years)}")
else:
    files = all_files
 
# ── Show size of all files (HEAD requests) ────────────────────────────────────
# print(f"\n{'File':<55} {'Size':>10}")
# print("-" * 67)
# total_bytes = 0
# file_sizes = {}
# for fname in all_files:
#     try:
#         h = requests.head(BASE_URL + fname, timeout=15, allow_redirects=True)
#         size = int(h.headers.get("content-length", 0))
#     except Exception:
#         size = 0
#     file_sizes[fname] = size
#     total_bytes += size
#     marker = " ◀ will download" if fname in files else ""
#     print(f"{fname:<55} {size/1e6:>8.1f} MB{marker}")
 


h = requests.head(BASE_URL + all_files[1], timeout=15, allow_redirects=True)
size = int(h.headers.get("content-length", 0))
print(f"first file {size/1e6:>8.1f} MB")
print("-" * 67)
total_bytes=len(all_files)*size
print(f"{'Total (all files)':<55} {total_bytes/1e6:>8.1f} MB")

 
# ── Download selected files ───────────────────────────────────────────────────
print(f"\nSaving to: {OUT_DIR.resolve()}\n")
downloaded = skipped = errors = 0
 
for i, fname in enumerate(files, 1):
    out_path = OUT_DIR / fname
    if out_path.exists():
        print(f"[{i}/{len(files)}] Skipping (exists): {fname}")
        skipped += 1
        continue
 
    url = BASE_URL + fname
    print(f"[{i}/{len(files)}] Downloading: {fname}", end=" ... ", flush=True)
    try:
        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(out_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 256):
                    f.write(chunk)
        size_mb = out_path.stat().st_size / 1e6
        print(f"OK ({size_mb:.1f} MB)")
        downloaded += 1
    except Exception as e:
        print(f"FAILED ({e})")
        out_path.unlink(missing_ok=True)
        errors += 1
 
print(f"\nDone.  Downloaded: {downloaded}  |  Skipped: {skipped}  |  Errors: {errors}")