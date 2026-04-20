"""
scrape_bioinformatiker_salary.py
================================
Scrapes bioinformatician salary data for Germany by Bundesland from
gehaltsvergleich.com (live + Wayback Machine 2019-present), then
automatically generates three publication-quality figure variants via R.

Pipeline
--------
  1. Live scrape of gehaltsvergleich.com (today's data)
  2. Wayback Machine CDX API -> enumerate all archived snapshots
  3. Fetch + parse each unique historical snapshot (cached between runs)
  4. Merge into data/gehaltsvergleich_all.csv
  5. Write embedded R script to disk and run it via subprocess

Output (data/)
--------------
  gehaltsvergleich_live_YYYY-MM-DD.csv   -- today's scrape
  gehaltsvergleich_wayback_history.csv   -- historical snapshots (cache)
  gehaltsvergleich_all.csv               -- merged dataset, input to R

Output (figures, 3 variants x PNG + PDF)
-----------------------------------------
  fig_v1_map_dumbbell_spaghetti.png/.pdf
  fig_v2_map_dumbbell_ribbon.png/.pdf
  fig_v3_mapmultiple_dumbbell.png/.pdf   (A4 portrait)

Usage
-----
  python scrape_bioinformatiker_salary.py

  Re-run once a year to add a new live snapshot and regenerate figures.
  Pass --skip-scrape to regenerate figures only (no network requests).
  Pass --skip-plot  to scrape only (no R execution).

Environment
-----------
  micromamba install -c conda-forge \\
    python requests beautifulsoup4 lxml pandas \\
    r-ggplot2 r-dplyr r-tidyr r-sf r-rnaturalearth r-rnaturalearthdata \\
    r-patchwork r-scales r-ggtext r-glue r-ggrepel r-readr -y
"""

import argparse
import re
import subprocess
import sys
import textwrap
import time
from collections import Counter
from datetime import date
from pathlib import Path

import pandas as pd
import requests
from bs4 import BeautifulSoup

# ── Constants ─────────────────────────────────────────────────────────────────

GV_NATIONAL_URL = (
    "https://www.gehaltsvergleich.com/gehalt/Bioinformatiker-Bioinformatikerin"
)

WAYBACK_CDX_URL = "http://web.archive.org/cdx/search/cdx"
WAYBACK_BASE    = "https://web.archive.org/web"

DATA_DIR = Path("data")

CRAWL_DELAY = 1.5   # seconds between requests — be polite to both servers


# ── Helpers ───────────────────────────────────────────────────────────────────

def make_session() -> requests.Session:
    s = requests.Session()
    s.headers.update({
        "User-Agent": (
            "Mozilla/5.0 (X11; Linux x86_64; rv:124.0) "
            "Gecko/20100101 Firefox/124.0"
        ),
        "Accept-Language": "de,en-US;q=0.7,en;q=0.3",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    })
    return s


def safe_get(session, url, params=None, timeout=30, retries=3, backoff=3.0):
    """GET with exponential backoff retry. Returns Response or raises."""
    for attempt in range(retries):
        try:
            r = session.get(url, params=params, timeout=timeout)
            r.raise_for_status()
            return r
        except requests.RequestException as e:
            if attempt == retries - 1:
                raise
            wait = backoff * (2 ** attempt)
            print(f"    Retry {attempt+1}/{retries} ({e}) -- waiting {wait:.0f}s")
            time.sleep(wait)


def parse_euro(text: str):
    """Parse '4.385 EUR' or '3.352EUR' -> 4385. Returns None if unparseable."""
    if not text:
        return None
    cleaned = re.sub(r"[^\d]", "", text)
    return int(cleaned) if cleaned else None


# ── HTML parser ───────────────────────────────────────────────────────────────

def parse_salary_table(html: str, scrape_date: str, source_url: str) -> list:
    """
    Parse the national gehaltsvergleich.com salary table.

    The national page has one table with columns:
        Region | Q1 | Oe | Q3 | Offene Jobs

    Returns a list of row dicts, one per Bundesland + Deutschland.
    Returns empty list if no matching table found.
    """
    soup = BeautifulSoup(html, "lxml")
    rows = []

    for table in soup.find_all("table"):
        headers = [th.get_text(strip=True) for th in table.find_all("th")]
        # Only process the table that has salary quartile columns
        if not any(h in ("Q1", "Oe", "1. Quartil", "Mittelwert") for h in headers):
            continue

        for tr in table.find_all("tr")[1:]:   # skip header row
            cells = [td.get_text(strip=True) for td in tr.find_all("td")]
            if len(cells) < 3:
                continue

            region = cells[0].rstrip(":").strip()
            q1     = parse_euro(cells[1]) if len(cells) > 1 else None
            avg    = parse_euro(cells[2]) if len(cells) > 2 else None
            q3     = parse_euro(cells[3]) if len(cells) > 3 else None

            if avg is None or not region:
                continue

            rows.append({
                "region":      region,
                "q1":          q1,
                "avg":         avg,
                "q3":          q3,
                "scrape_date": scrape_date,
                "source_url":  source_url,
            })

    return rows


# ── Live scrape ───────────────────────────────────────────────────────────────

def scrape_live(session: requests.Session) -> pd.DataFrame:
    """Scrape today's live gehaltsvergleich.com national page."""
    today = str(date.today())
    print(f"\n{'='*60}")
    print("SOURCE 1: Live scrape")
    print(f"  URL: {GV_NATIONAL_URL}")
    print(f"{'='*60}")

    r = safe_get(session, GV_NATIONAL_URL)
    rows = parse_salary_table(r.text, today, GV_NATIONAL_URL)

    if not rows:
        print("  FAILED: No rows parsed -- check HTML structure")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    print(f"  OK: {len(df)} rows parsed for {today}")
    print(df[["region", "q1", "avg", "q3"]].to_string(index=False))
    return df


# ── Wayback CDX: enumerate snapshots ─────────────────────────────────────────

def get_wayback_snapshots(session: requests.Session) -> list:
    """
    Query Wayback CDX API for all unique-content snapshots of the
    gehaltsvergleich.com Bioinformatiker national page.

    Returns list of dicts with keys: timestamp, year, wayback_url
    Sorted chronologically oldest-first.
    """
    print(f"\n{'='*60}")
    print("SOURCE 2: Wayback Machine CDX -- enumerating snapshots")
    print(f"{'='*60}")

    r = safe_get(
        session,
        WAYBACK_CDX_URL,
        params={
            "url":      "www.gehaltsvergleich.com/gehalt/Bioinformatiker-Bioinformatikerin",
            "output":   "json",
            "fl":       "timestamp,statuscode,digest",
            "filter":   "statuscode:200",
            "collapse": "digest",   # deduplicate by content hash
        },
        timeout=30,
    )

    data = r.json()
    if len(data) <= 1:
        print("  No snapshots found")
        return []

    headers_row = data[0]
    records     = data[1:]
    ts_idx      = headers_row.index("timestamp")

    snapshots = []
    for rec in records:
        ts  = rec[ts_idx]
        yr  = ts[:4]
        snapshots.append({
            "timestamp":   ts,
            "year":        yr,
            "wayback_url": f"{WAYBACK_BASE}/{ts}/{GV_NATIONAL_URL}",
        })

    snapshots.sort(key=lambda x: x["timestamp"])

    print(f"  Found {len(snapshots)} unique-content snapshots:")
    for yr, cnt in sorted(Counter(s["year"] for s in snapshots).items()):
        print(f"    {yr}: {cnt} snapshot(s)")

    return snapshots


# ── Wayback: fetch and parse each snapshot ────────────────────────────────────

def scrape_wayback_history(
    session:    requests.Session,
    snapshots:  list,
    cache_path: Path,
) -> pd.DataFrame:
    """
    Fetch and parse each Wayback snapshot.

    Results are cached to cache_path so re-runs skip already-fetched
    timestamps. Saves incrementally after each successful fetch.
    """
    print(f"\n{'='*60}")
    print(f"Fetching {len(snapshots)} Wayback snapshots")
    print(f"Cache: {cache_path}")
    print(f"{'='*60}")

    # Load existing cache
    if cache_path.exists():
        existing       = pd.read_csv(cache_path)
        done_stamps    = set(existing["timestamp"].astype(str))
        all_rows       = existing.to_dict("records")
        print(f"  Cache loaded: {len(existing)} rows, {len(done_stamps)} timestamps done")
    else:
        done_stamps = set()
        all_rows    = []

    new_rows = 0
    for snap in snapshots:
        ts  = snap["timestamp"]
        url = snap["wayback_url"]
        yr  = snap["year"]

        if ts in done_stamps:
            print(f"  SKIP  {ts[:8]}  ({yr})  [cached]")
            continue

        print(f"  FETCH {ts[:8]}  ({yr})  ", end="", flush=True)
        time.sleep(CRAWL_DELAY)

        try:
            r = safe_get(session, url, timeout=45)
        except requests.RequestException as e:
            print(f"FAILED -- {e}")
            continue

        # Derive scrape_date from timestamp (YYYYMMDDHHmmss -> YYYY-MM-DD)
        scrape_date = f"{ts[:4]}-{ts[4:6]}-{ts[6:8]}"
        rows = parse_salary_table(r.text, scrape_date, url)

        if not rows:
            print("no table found in archive page")
            continue

        for row in rows:
            row["timestamp"] = ts

        all_rows.extend(rows)
        new_rows += len(rows)
        print(f"OK  ({len(rows)} rows)")

        # Save incrementally -- safe against interruptions
        pd.DataFrame(all_rows).to_csv(cache_path, index=False)

    total = len(all_rows)
    print(f"\n  Done: {total} total rows ({new_rows} new, {total - new_rows} cached)")

    return pd.DataFrame(all_rows) if all_rows else pd.DataFrame()


# ── Merge live + history ──────────────────────────────────────────────────────

def merge_all(df_live: pd.DataFrame, df_history: pd.DataFrame) -> pd.DataFrame:
    """
    Combine live and historical DataFrames.
    Deduplicates by (region, scrape_date), live data wins over wayback.
    Adds a 'year' int column for easy grouping in R.
    """
    frames = []

    if not df_live.empty:
        df_l = df_live.copy()
        df_l["timestamp"] = df_l["scrape_date"].str.replace("-", "") + "000000"
        df_l["data_type"] = "live"
        frames.append(df_l)

    if not df_history.empty:
        df_h = df_history.copy()
        df_h["data_type"] = "wayback"
        frames.append(df_h)

    if not frames:
        return pd.DataFrame()

    df = pd.concat(frames, ignore_index=True)

    # live wins over wayback when dates overlap
    df = df.sort_values(
        ["region", "scrape_date", "data_type"],
        ascending=[True, True, False],
    )
    df = df.drop_duplicates(subset=["region", "scrape_date"], keep="first")
    df = df.sort_values(["region", "scrape_date"]).reset_index(drop=True)
    df["year"] = df["scrape_date"].str[:4].astype(int)

    return df


# ── Embedded R script ─────────────────────────────────────────────────────────

R_SCRIPT = r"""
# =============================================================================
# Bioinformatiker Salary in Germany — Multi-variant figures
# Auto-generated by scrape_bioinformatiker_salary.py — do not edit by hand.
# Source: data/gehaltsvergleich_all.csv
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(scales)
library(ggtext)
library(glue)
library(ggrepel)
library(readr)

options(rnaturalearth.use_cache = FALSE)

OUT_DPI  <- 300
DATA_FILE <- "data/gehaltsvergleich_all.csv"

# ── 1. Load & clean ──────────────────────────────────────────────────────────

raw <- read_csv(DATA_FILE, show_col_types = FALSE)

east_states <- c(
  "Brandenburg", "Mecklenburg-Vorpommern", "Sachsen",
  "Sachsen-Anhalt", "Thüringen", "Berlin"
)

salary <- raw |>
  mutate(year = as.integer(year)) |>
  group_by(region, year) |>
  slice_max(scrape_date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(group = ifelse(region %in% east_states, "East", "West"))

salary_states <- salary |> filter(region != "Deutschland")
salary_de     <- salary |> filter(region == "Deutschland")
snap_2026     <- salary_states |> filter(year == 2026)

west_order <- snap_2026 |> filter(group == "West") |>
  arrange(desc(avg)) |> pull(region)
east_order <- snap_2026 |> filter(group == "East") |>
  arrange(desc(avg)) |> pull(region)
state_order <- c(west_order, east_order)

snap_2026     <- snap_2026     |> mutate(region_f = factor(region, levels = rev(state_order)))
salary_states <- salary_states |> mutate(region_f = factor(region, levels = rev(state_order)))

# ── 2. Shapefile ─────────────────────────────────────────────────────────────

de_sf_base <- ne_states(country = "Germany", returnclass = "sf") |>
  select(name, geometry) |>
  mutate(name = recode(name,
    "North Rhine-Westphalia" = "Nordrhein-Westfalen",
    "Bavaria"                = "Bayern",
    "Hesse"                  = "Hessen",
    "Lower Saxony"           = "Niedersachsen",
    "Rhineland-Palatinate"   = "Rheinland-Pfalz",
    "Saxony"                 = "Sachsen",
    "Saxony-Anhalt"          = "Sachsen-Anhalt",
    "Thuringia"              = "Thüringen"
  ))

make_map_sf <- function(yr_data) {
  de_sf_base |> left_join(yr_data, by = c("name" = "region"))
}

unmatched <- make_map_sf(snap_2026) |> filter(is.na(avg)) |> pull(name)
if (length(unmatched) > 0)
  warning("Unmatched states (grey on map): ", paste(unmatched, collapse = ", "))

# ── 3. Palette & theme ────────────────────────────────────────────────────────

CLR_WEST <- "#2166AC"
CLR_EAST <- "#D6604D"
CLR_DE   <- "#1B7837"
CLR_BG   <- "#F7F7F7"

y_colours <- ifelse(
  levels(snap_2026$region_f) %in% east_states, CLR_EAST, CLR_WEST
)

euro_labels <- label_comma(big.mark = ".", decimal.mark = ",", prefix = "\u20ac")

base_theme <- theme_minimal(base_size = 9) +
  theme(
    plot.background   = element_rect(fill = CLR_BG, colour = NA),
    panel.background  = element_rect(fill = CLR_BG, colour = NA),
    plot.title        = element_markdown(size = 10, face = "bold", margin = margin(b = 2)),
    plot.subtitle     = element_markdown(size = 7.5, colour = "grey45", margin = margin(b = 4)),
    axis.title        = element_text(size = 7.5, colour = "grey35"),
    axis.text         = element_text(size = 7),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey90", linewidth = 0.3),
    legend.position   = "bottom",
    legend.key.width  = unit(1.4, "cm"),
    legend.key.height = unit(0.28, "cm"),
    legend.text       = element_text(size = 7),
    legend.title      = element_text(size = 7.5)
  )

# ── 4. Panel builders ─────────────────────────────────────────────────────────

build_map_single <- function(yr = 2026) {
  yr_data <- salary_states |> filter(year == yr)
  sf_data <- make_map_sf(yr_data)
  de_avg  <- salary_de |> filter(year == yr) |> pull(avg)
  ggplot(sf_data) +
    geom_sf(aes(fill = avg), colour = "white", linewidth = 0.55) +
    scale_fill_distiller(
      palette = "Blues", direction = 1,
      name    = "\u00d8 Bruttogehalt (\u20ac/Mo.) | Avg. gross salary",
      labels  = euro_labels,
      breaks  = c(3800, 4100, 4400, 4700),
      na.value = "grey80"
    ) +
    labs(
      title    = glue("**A** Average gross salary — Bioinformatiker/in, {yr}"),
      subtitle = glue(
        "National avg. (**DE**): **\u20ac{comma(de_avg, big.mark='.', decimal.mark=',')}**/mo. \u00b7 ",
        "Source: gehaltsvergleich.com"
      )
    ) +
    coord_sf(expand = FALSE) +
    base_theme +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      panel.grid.major = element_blank(), legend.position = "bottom"
    )
}

build_map_facet <- function() {
  years_to_show <- sort(unique(salary_states$year))
  sf_all <- lapply(years_to_show, function(yr) {
    make_map_sf(salary_states |> filter(year == yr)) |>
      mutate(year_label = as.character(yr))
  }) |> bind_rows()
  ggplot(sf_all) +
    geom_sf(aes(fill = avg), colour = "white", linewidth = 0.3) +
    scale_fill_distiller(
      palette = "Blues", direction = 1,
      name    = "\u00d8 Bruttogehalt (\u20ac/Mo.)",
      labels  = euro_labels,
      breaks  = c(3500, 4000, 4500, 5000),
      na.value = "grey80"
    ) +
    facet_wrap(~year_label, nrow = 2) +
    labs(
      title    = "**A** Avg. gross salary by year — Bioinformatiker/in",
      subtitle = "Source: gehaltsvergleich.com + Wayback Machine"
    ) +
    coord_sf(expand = FALSE) +
    base_theme +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      strip.text = element_text(size = 7, face = "bold"),
      legend.position = "bottom"
    )
}

build_dumbbell <- function() {
  de_avg_2026 <- salary_de |> filter(year == 2026) |> pull(avg)
  ggplot(snap_2026) +
    geom_segment(
      aes(x = q1, xend = q3, y = region_f, yend = region_f, colour = group),
      linewidth = 1.1, alpha = 0.5
    ) +
    geom_point(aes(x = q1,  y = region_f, colour = group),
               shape = 21, fill = "white", size = 2.2, stroke = 1.0) +
    geom_point(aes(x = q3,  y = region_f, colour = group),
               shape = 21, fill = "white", size = 2.2, stroke = 1.0) +
    geom_point(aes(x = avg, y = region_f, colour = group),
               shape = 18, size = 3.0) +
    geom_vline(xintercept = de_avg_2026,
               linetype = "dashed", colour = CLR_DE, linewidth = 0.5) +
    annotate("text", x = de_avg_2026 + 60, y = 0.6,
             label = "DE avg.", colour = CLR_DE, size = 2.3, hjust = 0) +
    scale_colour_manual(values = c("West" = CLR_WEST, "East" = CLR_EAST), name = NULL) +
    scale_x_continuous(labels = euro_labels,
                       breaks = c(3500, 4000, 4500, 5000, 5500, 6000, 9000)) +
    labs(
      title    = "**B** Salary distribution by state, 2026",
      subtitle = "\u25cb Q1 &nbsp;\u25c6 Mean &nbsp;\u25cb Q3 \u00b7 West = blue, East = red",
      x = "Gross monthly salary (\u20ac)", y = NULL
    ) +
    base_theme +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position    = "none",
      axis.text.y        = element_text(size = 7, colour = y_colours)
    )
}

build_spaghetti <- function() {
  de_line <- salary_de |> filter(!is.na(avg))
  ggplot(salary_states |> filter(!is.na(avg)),
         aes(x = year, y = avg, group = region, colour = group)) +
    geom_line(linewidth = 0.55, alpha = 0.7) +
    geom_point(size = 1.4, alpha = 0.85) +
    geom_line(data = de_line, aes(x = year, y = avg, group = 1),
              colour = CLR_DE, linewidth = 0.9, linetype = "dashed",
              inherit.aes = FALSE) +
    geom_label_repel(
      data = salary_states |> filter(year == max(year), !is.na(avg)),
      aes(label = region, colour = group),
      size = 2.1, label.padding = 0.1, label.size = 0,
      fill = alpha(CLR_BG, 0.7), show.legend = FALSE,
      nudge_x = 0.2, direction = "y", hjust = 0, segment.size = 0.2
    ) +
    annotate("text", x = max(de_line$year) + 0.15,
             y = tail(de_line$avg[order(de_line$year)], 1) + 80,
             label = "DE avg.", colour = CLR_DE, size = 2.2, hjust = 0) +
    scale_colour_manual(values = c("West" = CLR_WEST, "East" = CLR_EAST), name = NULL) +
    scale_x_continuous(
      breaks = sort(unique(salary_states$year)),
      expand = expansion(mult = c(0.02, 0.22))
    ) +
    scale_y_continuous(labels = euro_labels) +
    labs(
      title    = "**C** Avg. salary trend by state, 2019\u20132026",
      subtitle = "Dashed line = national average",
      x = NULL, y = "Gross monthly salary (\u20ac)"
    ) +
    base_theme +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 6.5, angle = 30, hjust = 1)
    )
}

build_ribbon <- function() {
  df <- salary_states |> filter(!is.na(avg), !is.na(q1), !is.na(q3))
  ggplot(df, aes(x = year)) +
    geom_ribbon(aes(ymin = q1, ymax = q3, fill = group), alpha = 0.25) +
    geom_line(aes(y = avg, colour = group), linewidth = 0.6) +
    geom_point(aes(y = avg, colour = group), size = 1.2) +
    facet_wrap(~region_f, ncol = 4) +
    scale_fill_manual(values   = c("West" = CLR_WEST, "East" = CLR_EAST), name = NULL) +
    scale_colour_manual(values = c("West" = CLR_WEST, "East" = CLR_EAST), name = NULL) +
    scale_x_continuous(breaks = c(2019, 2022, 2026), labels = c("'19", "'22", "'26")) +
    scale_y_continuous(
      labels = label_number(scale = 1e-3, suffix = "k", prefix = "\u20ac"),
      limits = c(3000, 6000),
      oob    = scales::squish
    ) +
    labs(
      title    = "**C** Q1\u2013avg\u2013Q3 salary band per state, 2019\u20132026",
      subtitle = "Shaded = Q1\u2013Q3 range \u00b7 Line = mean \u00b7 Shared y-axis \u20ac3k\u2013\u20ac6k (BaWü Q3 outlier squished)",
      x = NULL, y = "Gross monthly salary (\u20ac)"
    ) +
    base_theme +
    theme(
      strip.text      = element_text(size = 6, face = "bold"),
      axis.text       = element_text(size = 6),
      legend.position = "bottom",
      panel.spacing   = unit(0.3, "lines")
    )
}

# ── 5. Shared annotation ──────────────────────────────────────────────────────

anno <- plot_annotation(
  title    = "Bioinformatiker/in salaries in Germany by federal state",
  subtitle = paste0(
    "Gross monthly salary \u00b7 ",
    "Source: gehaltsvergleich.com (scraped 2019\u20132026 via Wayback Machine)"
  ),
  caption  = paste0(
    "\u25cb = Q1 \u00b7 \u25c6 = Mean \u00b7 \u25cb = Q3 \u00b7 ",
    "Blue = West Germany \u00b7 Red = East Germany"
  ),
  theme = theme(
    plot.background = element_rect(fill = CLR_BG, colour = NA),
    plot.title      = element_text(size = 13, face = "bold",    margin = margin(b = 3)),
    plot.subtitle   = element_text(size = 9,  colour = "grey40", margin = margin(b = 6)),
    plot.caption    = element_text(size = 7,  colour = "grey55", margin = margin(t = 6))
  )
)

save_fig <- function(fig, stem, w = 297, h = 210) {
  for (ext in c("png", "pdf")) {
    fn <- paste0(stem, ".", ext)
    ggsave(fn, plot = fig, width = w, height = h, units = "mm",
           dpi = OUT_DPI, bg = CLR_BG,
           device = if (ext == "pdf") cairo_pdf else NULL)
    message("  Saved: ", fn)
  }
}

# ── 6. Build and save ─────────────────────────────────────────────────────────

message("\n--- Building panels ---")
p_map_single <- build_map_single(2026)
p_map_facet  <- build_map_facet()
p_dumbbell   <- build_dumbbell()
p_spaghetti  <- build_spaghetti()
p_ribbon     <- build_ribbon()
message("Panels OK\n")

message("--- Variant 1: single map | dumbbell | spaghetti ---")
fig_v1 <- (p_map_single | (p_dumbbell / p_spaghetti)) +
  plot_layout(widths = c(1, 1.1)) + anno
save_fig(fig_v1, "fig_v1_map_dumbbell_spaghetti")

message("--- Variant 2: single map | dumbbell | ribbon ---")
fig_v2 <- (p_map_single | (p_dumbbell / p_ribbon)) +
  plot_layout(widths = c(0.9, 1.2)) + anno
save_fig(fig_v2, "fig_v2_map_dumbbell_ribbon")

message("--- Variant 3: faceted maps | dumbbell (A4 portrait) ---")
fig_v3 <- (p_map_facet / p_dumbbell) +
  plot_layout(heights = c(1.4, 1)) + anno
save_fig(fig_v3, "fig_v3_mapmultiple_dumbbell", w = 210, h = 297)

message("\nAll figures saved. Review fig_v1, fig_v2, fig_v3 and pick your favourite.")
"""


# ── R runner ──────────────────────────────────────────────────────────────────

R_SCRIPT_PATH = Path("bioinformatiker_plot.R")


def run_r_plots():
    """Write the embedded R script to disk and execute it via Rscript."""
    print(f"\n{'='*60}")
    print("PLOTTING: Writing R script and running via Rscript")
    print(f"{'='*60}")

    R_SCRIPT_PATH.write_text(R_SCRIPT.strip())
    print(f"  R script written to: {R_SCRIPT_PATH}")

    # Find Rscript — try PATH first, then common conda locations
    import shutil
    rscript = shutil.which("Rscript")
    if rscript is None:
        # Try micromamba env bin
        import os
        conda_prefix = os.environ.get("CONDA_PREFIX", "")
        candidate = Path(conda_prefix) / "bin" / "Rscript"
        if candidate.exists():
            rscript = str(candidate)

    if rscript is None:
        print("  ERROR: Rscript not found in PATH or CONDA_PREFIX/bin.")
        print("  Install R packages and ensure Rscript is on your PATH, then run:")
        print(f"    Rscript {R_SCRIPT_PATH}")
        return False

    print(f"  Using Rscript: {rscript}")
    print("  Running... (this may take 1-2 minutes for the first run)\n")

    result = subprocess.run(
        [rscript, "--vanilla", str(R_SCRIPT_PATH)],
        capture_output=False,   # stream output directly to terminal
    )

    if result.returncode != 0:
        print(f"\n  ERROR: Rscript exited with code {result.returncode}")
        print(f"  Check output above for R errors.")
        return False

    print("\n  Plotting complete.")
    return True


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Scrape Bioinformatiker salary data and generate figures."
    )
    parser.add_argument(
        "--skip-scrape", action="store_true",
        help="Skip scraping — use existing data/gehaltsvergleich_all.csv"
    )
    parser.add_argument(
        "--skip-plot", action="store_true",
        help="Skip R plotting — scrape only"
    )
    args = parser.parse_args()

    DATA_DIR.mkdir(exist_ok=True)
    today   = str(date.today())
    session = make_session()
    errors  = []

    # ── Scraping ──────────────────────────────────────────────────────────────
    if not args.skip_scrape:

        # 1. Live scrape
        try:
            df_live = scrape_live(session)
            if not df_live.empty:
                out = DATA_DIR / f"gehaltsvergleich_live_{today}.csv"
                df_live.to_csv(out, index=False)
                print(f"\n  Saved: {out}")
            else:
                errors.append("Live scrape returned no data")
        except Exception as e:
            print(f"\n  Live scrape error: {e}")
            errors.append(f"Live scrape: {e}")
            df_live = pd.DataFrame()

        # 2. Wayback history
        df_history = pd.DataFrame()
        cache_path = DATA_DIR / "gehaltsvergleich_wayback_history.csv"

        try:
            time.sleep(CRAWL_DELAY)
            snapshots = get_wayback_snapshots(session)
            if snapshots:
                time.sleep(CRAWL_DELAY)
                df_history = scrape_wayback_history(session, snapshots, cache_path)
            else:
                errors.append("Wayback: no snapshots found")
        except Exception as e:
            print(f"\n  Wayback error: {e}")
            errors.append(f"Wayback: {e}")
            if cache_path.exists():
                print(f"  Loading cached history from {cache_path}")
                df_history = pd.read_csv(cache_path)

        # 3. Merge and save
        df_all = merge_all(df_live, df_history)
        if not df_all.empty:
            out_all = DATA_DIR / "gehaltsvergleich_all.csv"
            df_all.to_csv(out_all, index=False)
            print(f"\n{'='*60}")
            print("FINAL DATASET")
            print(f"  Rows:    {len(df_all)}")
            print(f"  Years:   {sorted(df_all['year'].unique())}")
            print(f"  Regions: {df_all['region'].nunique()}")
            print(f"  Saved:   {out_all}")

    else:
        all_csv = DATA_DIR / "gehaltsvergleich_all.csv"
        if not all_csv.exists():
            print(f"ERROR: --skip-scrape set but {all_csv} does not exist.")
            sys.exit(1)
        print(f"\n  Scraping skipped. Using existing {all_csv}")

    # ── Plotting ──────────────────────────────────────────────────────────────
    if not args.skip_plot:
        plot_ok = run_r_plots()
        if not plot_ok:
            errors.append("R plotting failed — see output above")

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    if errors:
        print(f"Completed with {len(errors)} issue(s):")
        for e in errors:
            print(f"  WARNING: {e}")
    else:
        print("All steps completed successfully.")


if __name__ == "__main__":
    main()
