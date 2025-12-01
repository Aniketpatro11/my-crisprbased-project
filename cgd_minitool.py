import streamlit as st
import pandas as pd
import numpy as np
import textwrap

# -------------------------------
# Basic page config
# -------------------------------
st.set_page_config(
    page_title="CRISPR Guide Design Mini-Tool",
    layout="wide"
)

# -------------------------------
# Helper functions
# -------------------------------

def clean_sequence(seq: str) -> str:
    """Keep only A/T/G/C and uppercase."""
    seq = seq.upper()
    allowed = {"A", "T", "G", "C"}
    return "".join([b for b in seq if b in allowed])


def parse_fasta(text: str) -> str:
    """Very simple FASTA parser: returns concatenated sequence lines."""
    lines = text.strip().splitlines()
    seq_lines = [ln.strip() for ln in lines if not ln.startswith(">")]
    return "".join(seq_lines)


def gc_content(seq: str) -> float:
    if len(seq) == 0:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return 100.0 * gc / len(seq)


def complement(base: str) -> str:
    comp_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return comp_map.get(base, "N")


def rev_complement(seq: str) -> str:
    return "".join(complement(b) for b in seq[::-1])


def self_complementarity_score(seq: str, window: int = 4) -> int:
    """
    Very simplified self-complementarity measure:
    count how many times a small window has its reverse complement
    elsewhere in the sequence.
    """
    n = len(seq)
    score = 0
    for i in range(n - window + 1):
        sub = seq[i:i+window]
        sub_rc = rev_complement(sub)
        # search in remaining sequence
        rest = seq[:i] + seq[i+window:]
        score += rest.count(sub_rc)
    return score


def sliding_windows(seq: str, k: int):
    """Generate all k-length windows and their start indices."""
    for i in range(len(seq) - k + 1):
        yield i, seq[i:i+k]


def off_target_score(seq: str, guide: str, start_idx: int, max_mismatches: int = 5) -> int:
    """
    Simplified off-target estimation:
    count how many other windows in the SAME sequence
    match with <= max_mismatches.
    """
    k = len(guide)
    score = 0
    for i, window in sliding_windows(seq, k):
        # skip the exact on-target region
        if i == start_idx:
            continue
        mismatches = sum(1 for a, b in zip(window, guide) if a != b)
        if mismatches <= max_mismatches:
            score += 1
    return score


def find_guides_forward(seq: str, guide_len: int = 20, pam: str = "NGG"):
    """
    Find guides on the FORWARD strand for SpCas9-like NGG PAM.
    Guide is 20 nt upstream of NGG.
    """
    guides = []

    if pam != "NGG":
        # for now only support NGG pattern
        return guides

    n = len(seq)
    pam_len = 3

    for guide_start in range(0, n - guide_len - pam_len + 1):
        pam_start = guide_start + guide_len
        pam_seq = seq[pam_start:pam_start + pam_len]

        # match N GG (second and third are G)
        if len(pam_seq) == 3 and pam_seq[1:] == "GG":
            guide_seq = seq[guide_start:guide_start + guide_len]
            guides.append({
                "strand": "+",
                "guide_start": guide_start,
                "guide_end": guide_start + guide_len,  # exclusive
                "pam_start": pam_start,
                "pam_end": pam_start + pam_len,       # exclusive
                "guide_seq": guide_seq,
                "pam_seq": pam_seq
            })

    return guides


def compute_guide_table(seq: str, guide_len: int = 20, pam: str = "NGG",
                        advanced_scores: bool = True):
    guides = find_guides_forward(seq, guide_len=guide_len, pam=pam)

    records = []
    for g in guides:
        guide_seq = g["guide_seq"]
        gc = gc_content(guide_seq)

        if advanced_scores:
            self_comp = self_complementarity_score(guide_seq, window=4)
            off_target = off_target_score(seq, guide_seq, g["guide_start"], max_mismatches=5)
        else:
            self_comp = np.nan
            off_target = np.nan

        # Simple heuristic quality scoring
        gc_penalty = abs(gc - 50) / 5.0  # 0 if 50%, grows as we move away
        self_penalty = self_comp
        off_penalty = off_target * 2 if not np.isnan(off_target) else 0
        total_score = gc_penalty + self_penalty + off_penalty

        records.append({
            "Strand": g["strand"],
            "Guide Start (0-based)": g["guide_start"],
            "Guide End (0-based, excl)": g["guide_end"],
            "PAM Start": g["pam_start"],
            "PAM End": g["pam_end"],
            "Guide Sequence (5'→3')": guide_seq,
            "PAM": g["pam_seq"],
            "GC %": round(gc, 2),
            "Self-Complementarity (simplified)": self_comp,
            "Off-target Matches (≤5 mismatches, same seq)": off_target,
            "Total Score (lower is better)": round(total_score, 2)
        })

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame.from_records(records)
    df = df.sort_values(by="Total Score (lower is better)", ascending=True).reset_index(drop=True)
    return df


def highlight_sequence_html(seq: str, guides_meta, line_width: int = 60) -> str:
    """
    Render sequence with <span> highlighting for guides (green) and PAMs (red).
    """
    n = len(seq)
    # For each position, mark what it is
    pos_type = ["normal"] * n

    for g in guides_meta:
        for i in range(g["guide_start"], g["guide_end"]):
            pos_type[i] = "guide"
        for i in range(g["pam_start"], g["pam_end"]):
            pos_type[i] = "pam"

    # Build HTML line by line
    html_lines = []
    for line_start in range(0, n, line_width):
        line_end = min(line_start + line_width, n)
        line_chars = []
        for i in range(line_start, line_end):
            base = seq[i]
            t = pos_type[i]
            if t == "guide":
                span = f"<span style='color:#0b8d0b;font-weight:bold'>{base}</span>"
            elif t == "pam":
                span = f"<span style='color:#d00000;font-weight:bold'>{base}</span>"
            else:
                span = base
            line_chars.append(span)
        line_html = "".join(line_chars)
        # add index label on left
        label = f"<span style='color:gray;font-size:12px'>[{line_start:04d}] </span>"
        html_lines.append(label + line_html)

    return "<br>".join(html_lines)


# -------------------------------
# Sidebar – settings
# -------------------------------
st.sidebar.title("⚙️ Settings")

st.sidebar.markdown(
    """
This is an **educational** CRISPR guide design demo.  
It finds simple SpCas9-like guides (NGG PAM) on the **forward strand only**.
"""
)

pam_choice = st.sidebar.selectbox(
    "PAM pattern",
    options=["NGG"],
    index=0,
    help="Currently only NGG is supported for SpCas9 educational demo."
)

guide_len = st.sidebar.number_input(
    "Guide length (nt)",
    min_value=18,
    max_value=24,
    value=20,
    step=1,
    help="Typical SpCas9 guide RNA length is 20 nt."
)

advanced_scores = st.sidebar.checkbox(
    "Compute advanced scores (self-complementarity & off-target)",
    value=True
)

gc_min, gc_max = st.sidebar.slider(
    "Filter by GC% range",
    min_value=0,
    max_value=100,
    value=(30, 80),
    step=1
)

st.sidebar.markdown("---")
st.sidebar.markdown("Made for learning – not for clinical / lab use.")

# -------------------------------
# Main layout
# -------------------------------
st.title("🧬 CRISPR Guide Design Mini-Tool (Educational)")

st.markdown(
    """
This app teaches the **basics of CRISPR guide RNA selection**:

- Scans a DNA sequence for **SpCas9-like guides** (20 nt upstream of **NGG** PAM).
- Calculates **GC%** and simple **quality scores**.
- Shows **potential off-target-like matches within the same sequence** (very simplified).
- Highlights **guides (green)** and **PAMs (red)** on the sequence.

> ⚠️ **Important:** This is for education only, not for real experimental or clinical use.
"""
)

# -------------------------------
# Input area
# -------------------------------
st.subheader("1️⃣ Input DNA sequence")

col1, col2 = st.columns(2)

with col1:
    seq_textarea = st.text_area(
        "Paste DNA sequence or FASTA here (A/T/G/C only will be kept):",
        height=200,
        help="FASTA headers (lines starting with '>') will be ignored."
    )

with col2:
    fasta_file = st.file_uploader(
        "Or upload a FASTA file",
        type=["fa", "fasta", "txt"]
    )

raw_seq = ""

if fasta_file is not None:
    file_text = fasta_file.read().decode("utf-8")
    raw_seq = parse_fasta(file_text)
elif seq_textarea.strip():
    # Accept both raw DNA or FASTA-like
    if seq_textarea.lstrip().startswith(">"):
        raw_seq = parse_fasta(seq_textarea)
    else:
        raw_seq = seq_textarea

clean_seq = clean_sequence(raw_seq)

if not clean_seq:
    st.info("Provide a DNA sequence (via text area or FASTA upload) to begin.")
    st.stop()

st.success(f"Sequence loaded. Length after cleaning: **{len(clean_seq)} nt**")

# -------------------------------
# Compute guides
# -------------------------------
st.subheader("2️⃣ Detected guides (forward strand, NGG PAM)")

with st.spinner("Scanning sequence for guides and computing scores..."):
    df_guides = compute_guide_table(
        clean_seq,
        guide_len=guide_len,
        pam=pam_choice,
        advanced_scores=advanced_scores
    )

if df_guides.empty:
    st.warning(
        "No guides found with the current parameters. "
        "Try a longer sequence or adjust the settings."
    )
    st.stop()

# Filter by GC%
df_filtered = df_guides[
    (df_guides["GC %"] >= gc_min) & (df_guides["GC %"] <= gc_max)
].reset_index(drop=True)

st.markdown(
    f"Found **{len(df_guides)}** guides total, "
    f"**{len(df_filtered)}** within GC% range {gc_min}–{gc_max}."
)

# Show table
st.dataframe(df_filtered, use_container_width=True)

# -------------------------------
# Download results
# -------------------------------
csv_data = df_filtered.to_csv(index=False).encode("utf-8")
st.download_button(
    "⬇️ Download guide table as CSV",
    csv_data,
    file_name="crispr_guides_educational.csv",
    mime="text/csv"
)

# -------------------------------
# Sequence visualization
# -------------------------------
st.subheader("3️⃣ Sequence map with highlighted guides & PAMs")

# Build guides_meta for visualization from df_guides, not filtered
guides_meta = []
for _, row in df_guides.iterrows():
    guides_meta.append({
        "guide_start": int(row["Guide Start (0-based)"]),
        "guide_end": int(row["Guide End (0-based, excl)"]),
        "pam_start": int(row["PAM Start"]),
        "pam_end": int(row["PAM End"])
    })

seq_html = highlight_sequence_html(clean_seq, guides_meta, line_width=60)
st.markdown(
    """
    <div style="font-family:monospace; white-space:pre; line-height:1.4">
    {}
    </div>
    """.format(seq_html),
    unsafe_allow_html=True
)

st.caption("Guides are shown in **green**, PAMs (NGG) in **red**.")

# -------------------------------
# Simple stats / plots
# -------------------------------
st.subheader("4️⃣ Guide GC% and score distribution")

col_gc, col_score = st.columns(2)

with col_gc:
    st.markdown("**GC% distribution**")
    st.bar_chart(df_guides["GC %"])

with col_score:
    st.markdown("**Total score distribution** (lower is better, heuristic)")
    st.bar_chart(df_guides["Total Score (lower is better)"])

st.markdown(
    """
**How to interpret the scores (simplified heuristic):**

- ✅ **GC%** around 40–60% is usually preferred.
- ⚠ **Self-complementarity** suggests potential hairpins if high.
- ⚠ **Off-target matches** counted here are only within the same sequence and are highly simplified.
- **Total Score** combines these (lower ≈ nicer guide).

Again: this is a **teaching tool**, not a lab-ready design engine.
"""
)
