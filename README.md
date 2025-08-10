# hydrus-duplicate-finder-matlab

MATLAB tool for finding and marking **potential duplicate videos** in **Hydrus** using:
- **Duration** and optional **aspect ratio** prefilters
- Optional **perceptual hash (aHash)** prefilter
- **SSIM / MS-SSIM** thumbnail comparison

It also prints a **Top Tags** summary in the console and writes a **full CSV** tag report for files flagged as potential duplicates.

> ⚠️ The script **only posts “potential duplicate”** relationships (`relationship = 0`).  
> Already-marked pairs may return `400 Bad Request`; these are counted as **skipped**, not failures.

---

## Features

- Fast candidate generation by **duration** and optional **aspect ratio**
- Optional **aHash prefilter** to massively reduce SSIM workload
- **SSIM / MS-SSIM** thumbnail comparison (MS-SSIM if MATLAB R2021a+)
- Robust, batched posting to Hydrus with retries/backoff
- **Top 10 tags** (console) + **full CSV** of tags across duplicate files
- All helpers embedded—**single `.m` script** for easy sharing

---

## Requirements

- **Hydrus Client** with **Client API** enabled  
- **MATLAB** (R2020b+) with **Image Processing Toolbox**  
  - MS-SSIM requires **R2021a+** (`multissim`); otherwise script falls back to `ssim`

---

## Quick Start

1. **Enable Hydrus Client API**  
   `Services → Manage Services → Client API`

2. **Create an API key**  
   `Services → Manage Services → Client API → Permissions → API Key`

3. **Clone the repo**
   ```bash
   git clone https://github.com/<your-user>/hydrus-duplicate-finder-matlab.git
   cd hydrus-duplicate-finder-matlab
