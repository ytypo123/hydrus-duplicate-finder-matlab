# Hydrus Duplicate Finder & Potential Duplicate Setter (MATLAB)

A MATLAB script for detecting and marking **potential duplicate** files in [Hydrus Network](https://hydrusnetwork.github.io/hydrus/) using the Hydrus Client API.  
It performs duration, aspect ratio, and perceptual hash filtering, then uses SSIM/MS-SSIM thumbnail comparison to identify likely duplicates and post them to Hydrus.

---

## 🚀 Features
- **Automated duplicate detection** for Hydrus-managed files.
- **Multi-step candidate filtering**:
  - Duration filter (seconds tolerance)
  - Aspect ratio filter (fractional tolerance)
  - Optional aHash prefilter (perceptual hash)
  - SSIM / MS-SSIM thumbnail comparison
- **Hydrus Client API integration** to set "potential duplicate" relationships.
- **Tag frequency reporting** among detected duplicates.
- Saves **full CSV reports** with duplicate-tag statistics.

---

## 📦 Requirements
- **MATLAB R2021a+** (R2021a or newer recommended for `multissim` / MS-SSIM)
- **Image Processing Toolbox**
- **Hydrus Client API** enabled and running

---

## ⚡ Quick Start

1. **Enable Hydrus Client API**  
   In Hydrus Client:  
   `Services -> Manage Services -> Client API` → Enable it.

2. **Get an API Key**  
   In Hydrus Client:  
   `Services -> Manage Services -> Client API -> Permissions -> API Key`  
   Copy the key.

3. **Configure the Script**  
   Edit the **CONFIGURATION** section at the top:
   - `api_key` → Your Hydrus API key
   - `tags` → Tags to search (e.g. `{'system:filetype is video', 'system:filesize > 700MB'}`)
   - Tolerances → `tolSeconds` (duration), `arTol` (aspect ratio)
   - Similarity threshold → `similarity_threshold`
   - Toggle filters → `USE_AR_FILTER`, `USE_AHASH_PREFILTER`

4. **Run in MATLAB**  
   The script will:
   - Fetch file hashes & metadata from Hydrus
   - Filter by duration/aspect ratio
   - (Optional) Prefilter with perceptual hash
   - Compare thumbnails with SSIM / MS-SSIM
   - Post **"potential duplicate"** relationships to Hydrus
   - Print top duplicate tags and save CSV report

---

## 🛠 Configuration Parameters

| Parameter                | Description | Example |
|--------------------------|-------------|---------|
| `tags`                   | Tags to search in Hydrus | `{'system:filetype is video'}` |
| `api_key`                | Hydrus API key | `'xxxxxxxx...'` |
| `api_URL`                | Hydrus API base URL | `'http://127.0.0.1:42001'` |
| `USE_AR_FILTER`          | Enable aspect ratio filtering | `true` |
| `USE_AHASH_PREFILTER`    | Enable perceptual hash prefilter | `true` |
| `tolSeconds`             | Duration tolerance in seconds | `1.0` |
| `arTol`                  | Aspect ratio fractional tolerance | `0.01` |
| `similarity_threshold`   | SSIM similarity threshold (0–1) | `0.50` |
| `P` struct               | SSIM preprocessing & parameters | — |

---

## 📊 Output
- **Console**:  
  - Step-by-step log of filtering & comparisons  
  - Summary of duplicates found and posted to Hydrus  
  - Top tags among duplicates (configurable count & minimum frequency)
- **CSV file**:  
  - `hydrus_duplicate_tags_YYYY-MM-DD_HHMMSS.csv`  
  - Full list of tags, counts, and percentages for duplicate files

---

## ⚠ Notes
- Only posts **"potential duplicate"** (relationship = `0`) to Hydrus.
- Already-marked pairs may result in `400 Bad Request` — counted as skipped.
- Aggressive filtering (low `hammingMax`) may miss matches; loosen for recall.
- Weak filtering (high `hammingMax`) may increase SSIM workload.
- Running without aHash prefilter is slower but may find more matches.

---

## 📜 License
This project is provided **as-is** without warranty.  
Refer to `LICENSE` file for details (if applicable).

---

## 💡 Example Usage
```matlab
tags = {'system:filetype is video', 'system:filesize > 700MB'};
api_key = 'your-hydrus-api-key';
api_URL = 'http://127.0.0.1:42001';
USE_AR_FILTER = true;
USE_AHASH_PREFILTER = true;
tolSeconds = 1.0;
arTol = 0.01;
similarity_threshold = 0.50;

% Run script
hydrus_duplicate_finder;
