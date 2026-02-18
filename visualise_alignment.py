#!/usr/bin/env python3
from pathlib import Path
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# --- Simple ClustalX-like residue class colors (approx, journal-acceptable) ---
# Keep it conservative: chemical classes, not rainbow chaos.
CLUSTAL_LIKE = {
    # Hydrophobic
    "A": "#8dd3c7", "V": "#8dd3c7", "I": "#8dd3c7", "L": "#8dd3c7", "M": "#8dd3c7",
    "F": "#8dd3c7", "W": "#8dd3c7", "Y": "#8dd3c7",
    # Polar (uncharged)
    "S": "#ffffb3", "T": "#ffffb3", "N": "#ffffb3", "Q": "#ffffb3", "C": "#ffffb3",
    # Positive
    "K": "#bebada", "R": "#bebada", "H": "#bebada",
    # Negative
    "D": "#fb8072", "E": "#fb8072",
    # Special cases
    "G": "#80b1d3", "P": "#fdb462",
}
DEFAULT_COLOR = "#ffffff"   # background for unknowns/gaps
GAP_COLOR = "#f5f5f5"
EPI_BG = "#e8f0ff"          # subtle epitope column shading

def read_fasta(path: str):
    records = []
    header = None
    seq_parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq_parts)))
    return records

def build_ref_position_map(ref_aln: str):
    ref_pos = 0
    col_to_refpos = []
    for ch in ref_aln:
        if ch != "-":
            ref_pos += 1
            col_to_refpos.append(ref_pos)
        else:
            col_to_refpos.append(None)
    return col_to_refpos

def build_refpos_to_col_map(ref_aln: str):
    ref_pos = 0
    refpos_to_col = {}
    for i, ch in enumerate(ref_aln):
        if ch != "-":
            ref_pos += 1
            refpos_to_col[ref_pos] = i
    return refpos_to_col

def nearest_refpos(col_to_refpos, col_idx: int):
    if col_to_refpos[col_idx] is not None:
        return col_to_refpos[col_idx]
    for off in range(1, len(col_to_refpos)):
        l = col_idx - off
        r = col_idx + off
        if l >= 0 and col_to_refpos[l] is not None:
            return col_to_refpos[l]
        if r < len(col_to_refpos) and col_to_refpos[r] is not None:
            return col_to_refpos[r]
    return None

def aln_span_in_ref(seq_aln: str, col_to_refpos):
    cols = [i for i, ch in enumerate(seq_aln) if ch != "-"]
    if not cols:
        return (None, None)
    start = nearest_refpos(col_to_refpos, cols[0])
    end = nearest_refpos(col_to_refpos, cols[-1])
    return (start, end)

def shorten_id(s: str, maxlen: int = 46):
    if len(s) <= maxlen:
        return s
    return s[:maxlen-10] + "…" + s[-9:]

def find_epitope_cols(ep_aln: str):
    cols = [i for i, ch in enumerate(ep_aln) if ch != "-"]
    if not cols:
        return None
    return (cols[0], cols[-1])

def refpos_window_to_cols(refpos_to_col: dict, start_refpos: int, end_refpos: int):
    start_refpos = max(1, start_refpos)
    end_refpos = max(start_refpos, end_refpos)
    while start_refpos not in refpos_to_col and start_refpos < end_refpos:
        start_refpos += 1
    while end_refpos not in refpos_to_col and end_refpos > start_refpos:
        end_refpos -= 1
    if start_refpos not in refpos_to_col or end_refpos not in refpos_to_col:
        return None
    return (refpos_to_col[start_refpos], refpos_to_col[end_refpos])

def plot_coverage_page(ax, spans, ref_len, reference_id, ep_span, epitope_id, shorten_labels=True):
    palette = plt.get_cmap("tab10")
    bar_lw = 10
    n = len(spans)
    y_positions = list(range(n))[::-1]

    # Wider right margin so end-labels never crop
    ax.figure.subplots_adjust(left=0.42, right=0.93, top=0.88, bottom=0.22)

    pad = max(15, int(ref_len * 0.08))
    ax.set_xlim(1, ref_len + pad)

    for idx, (y, (hid, start, end)) in enumerate(zip(y_positions, spans)):
        ax.hlines(y=y, xmin=start, xmax=end, linewidth=bar_lw,
                  color=palette(idx % 10), alpha=0.95)

        label = hid if not shorten_labels else shorten_id(hid)
        ax.text(-0.02, y, label, transform=ax.get_yaxis_transform(),
                va="center", ha="right", fontsize=9, clip_on=False)

        ax.text(end + max(10, int(ref_len * 0.015)), y, f"{start}-{end}",
                va="center", ha="left", fontsize=8, clip_on=False)

    if ep_span is not None:
        ep_start, ep_end = ep_span
        y_ep = (max(y_positions) + 1) if y_positions else 1
        ax.hlines(y=y_ep, xmin=ep_start, xmax=ep_end, linewidth=6, color="black", alpha=0.9)
        ax.text(-0.02, y_ep, f"[epitope] {epitope_id}",
                transform=ax.get_yaxis_transform(),
                va="center", ha="right", fontsize=9, clip_on=False)
        ax.text((ep_start + ep_end) / 2, y_ep + 0.18, f"{ep_start}-{ep_end}",
                ha="center", va="bottom", fontsize=8, clip_on=False)

    ax.set_ylim(-1, (max(y_positions) + 2 if y_positions else 2))
    ax.set_xlabel(f"Reference coordinates (aa): {reference_id}")
    ax.set_yticks([])
    ax.set_title("KNL1 protein coverage map (alignment-derived)")
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)

def draw_alignment_zoom(ax, records_ordered, window_cols, col_to_refpos, ep_cols, title,
                        block=90, max_seqs=10, font_size=7.8):
    """
    Render a readable alignment window with residue-class coloring and epitope column shading.
    """
    start_col, end_col = window_cols
    window_len = end_col - start_col + 1

    ids = [hid for hid, _ in records_ordered][:max_seqs]
    seqs = [seq for _, seq in records_ordered][:max_seqs]
    window_seqs = [s[start_col:end_col+1] for s in seqs]
    window_ids = [shorten_id(h, 55) for h in ids]

    start_ref = nearest_refpos(col_to_refpos, start_col)
    end_ref = nearest_refpos(col_to_refpos, end_col)

    ax.axis("off")
    ax.set_title(f"{title}\nref {start_ref}-{end_ref} aa (zoom window)", fontsize=11)

    # epitope overlap within window (0-based within window)
    ep_marker = None
    if ep_cols is not None:
        ep_s, ep_e = ep_cols
        ov_s = max(ep_s, start_col)
        ov_e = min(ep_e, end_col)
        if ov_s <= ov_e:
            ep_marker = (ov_s - start_col, ov_e - start_col)

    # layout
    left_pad_chars = 60  # label width
    char_w = 0.0082      # tweak if needed
    line_h = 0.055       # tweak if needed

    y = 0.92
    for bs in range(0, window_len, block):
        be = min(window_len, bs + block)

        # block header
        ax.text(0.01, y, f"cols {start_col+bs}–{start_col+be-1}",
                fontsize=font_size, family="monospace", transform=ax.transAxes)
        y -= line_h * 0.9

        # draw epitope background band for this block (behind sequences)
        if ep_marker is not None:
            ep_s_w, ep_e_w = ep_marker
            ov_s = max(ep_s_w, bs)
            ov_e = min(ep_e_w + 1, be)  # exclusive
            if ov_s < ov_e:
                # convert character positions to axes coords
                x0 = 0.01 + char_w * (left_pad_chars + (ov_s - bs))
                x1 = 0.01 + char_w * (left_pad_chars + (ov_e - bs))
                # band spans the sequence lines for this block
                band_height = line_h * (len(window_ids)) * 1.02
                ax.add_patch(plt.Rectangle((x0, y - band_height + line_h*0.2),
                                           x1 - x0, band_height,
                                           transform=ax.transAxes,
                                           facecolor=EPI_BG, edgecolor="none", zorder=0.1))

        # sequence lines with per-residue background coloring
        for hid, wseq in zip(window_ids, window_seqs):
            seg = wseq[bs:be]
            # write label
            ax.text(0.01, y, f"{hid:<58} ", fontsize=font_size, family="monospace",
                    transform=ax.transAxes, zorder=1.0)

            # draw each residue as a colored cell (fast enough for ~100 cols × few seqs)
            x_start = 0.01 + char_w * left_pad_chars
            for j, ch in enumerate(seg):
                if ch == "-":
                    fc = GAP_COLOR
                else:
                    fc = CLUSTAL_LIKE.get(ch.upper(), DEFAULT_COLOR)
                ax.add_patch(plt.Rectangle((x_start + char_w*j, y - line_h*0.65),
                                           char_w*0.98, line_h*0.78,
                                           transform=ax.transAxes,
                                           facecolor=fc, edgecolor="none", zorder=0.5))
                ax.text(x_start + char_w*j, y, ch,
                        fontsize=font_size, family="monospace",
                        transform=ax.transAxes, va="top", ha="left", zorder=1.0)

            y -= line_h

        y -= line_h * 0.6
        if y < 0.08:
            break

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("aligned_fasta")
    ap.add_argument("--reference_id", required=True)
    ap.add_argument("--epitope_id", default="knl1_epitope")
    ap.add_argument("--include", nargs="*", default=None)
    ap.add_argument("--out_prefix", default="KNL1_covmap")
    ap.add_argument("--no_shorten", action="store_true")
    ap.add_argument("--zoom_pad", type=int, default=60, help="Padding around epitope (aa in reference coords)")
    ap.add_argument("--max_zoom_seqs", type=int, default=10)
    ap.add_argument("--block", type=int, default=90, help="Columns per block on zoom page")
    args = ap.parse_args()

    records = read_fasta(args.aligned_fasta)
    if not records:
        raise SystemExit("No FASTA records found.")
    aln_lens = {len(seq) for _, seq in records}
    if len(aln_lens) != 1:
        raise SystemExit(f"All sequences must be same aligned length; got {sorted(aln_lens)}")

    rec = {hid: seq for hid, seq in records}
    if args.reference_id not in rec:
        raise SystemExit(f"reference_id not found: {args.reference_id}")

    ref_aln = rec[args.reference_id]
    col_to_refpos = build_ref_position_map(ref_aln)
    refpos_to_col = build_refpos_to_col_map(ref_aln)
    ref_len = max([p for p in col_to_refpos if p is not None], default=0)

    ep_span_ref = None
    ep_cols = None
    if args.epitope_id in rec:
        ep_span_ref = aln_span_in_ref(rec[args.epitope_id], col_to_refpos)
        ep_cols = find_epitope_cols(rec[args.epitope_id])
        if ep_span_ref == (None, None):
            ep_span_ref = None
            ep_cols = None

    # coverage spans
    spans = []
    for hid, seq in records:
        if hid == args.epitope_id:
            continue
        if args.include is not None and hid not in args.include:
            continue
        s, e = aln_span_in_ref(seq, col_to_refpos)
        if s is None or e is None:
            continue
        spans.append((hid, s, e))
    spans.sort(key=lambda x: (0 if x[0] == args.reference_id else 1, x[1], x[2]))

    # zoom page ordering
    ordered_ids = [args.reference_id] + [hid for hid, _, _ in spans if hid != args.reference_id]
    seen = set()
    ordered_ids = [x for x in ordered_ids if not (x in seen or seen.add(x))]
    if args.epitope_id in rec:
        ordered_ids = [args.epitope_id] + ordered_ids
    records_ordered = [(hid, rec[hid]) for hid in ordered_ids if hid in rec]

    out_pdf = Path(args.out_prefix + ".pdf")
    out_svg = Path(args.out_prefix + ".svg")

    with PdfPages(out_pdf) as pdf:
        # Page 1
        fig1, ax1 = plt.subplots(figsize=(12, max(2.6, 0.60 * (len(spans) + (1 if ep_span_ref else 0)) + 0.8)))
        plot_coverage_page(
            ax=ax1,
            spans=spans,
            ref_len=ref_len,
            reference_id=args.reference_id,
            ep_span=ep_span_ref,
            epitope_id=args.epitope_id,
            shorten_labels=not args.no_shorten
        )
        pdf.savefig(fig1)
        fig1.savefig(out_svg)
        plt.close(fig1)

        # Page 2: zoomed alignment with epitope highlight + residue-class coloring
        if ep_span_ref is not None:
            ep_start_ref, ep_end_ref = ep_span_ref
            win_start = max(1, ep_start_ref - args.zoom_pad)
            win_end = min(ref_len, ep_end_ref + args.zoom_pad)
            window_cols = refpos_window_to_cols(refpos_to_col, win_start, win_end)
            if window_cols is not None:
                fig2 = plt.figure(figsize=(12, 7.5))
                ax2 = fig2.add_subplot(1, 1, 1)
                draw_alignment_zoom(
                    ax=ax2,
                    records_ordered=records_ordered,
                    window_cols=window_cols,
                    col_to_refpos=col_to_refpos,
                    ep_cols=ep_cols,
                    title="KNL1 epitope region (zoom; Clustal-like coloring)",
                    block=args.block,
                    max_seqs=args.max_zoom_seqs
                )
                pdf.savefig(fig2)
                plt.close(fig2)

    print(f"Wrote multi-page PDF: {out_pdf}")
    print(f"Wrote SVG (page 1): {out_svg}")

if __name__ == "__main__":
    main()

