
# gc_gut_analyzer_app.py

import streamlit as st
from Bio import SeqIO
import matplotlib.pyplot as plt
import random
import io

st.set_page_config(page_title="GC-Gut Analyzer", layout="centered")
st.title("🧬 GC-Gut Analyzer - Microbiome Insight Tool [Developed by: Dr. Keshab Nath]")

# Upload FASTA file
uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta", "fa"])

if uploaded_file is not None:
    # Convert binary file to text stream for Biopython
    text_io = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
    record = next(SeqIO.parse(text_io, "fasta"))
    seq = str(record.seq).upper()
    st.success(f"Sequence loaded: {record.id} ({len(seq)} bases)")

    # GC Content
    gc_content = 100 * (seq.count("G") + seq.count("C")) / len(seq)
    st.metric("GC Content (%)", f"{gc_content:.2f}%")

    # Base frequency bar chart
    st.subheader("Nucleotide Frequency")
    bases = ['A', 'T', 'G', 'C']
    counts = [seq.count(b) for b in bases]
    fig, ax = plt.subplots()
    ax.bar(bases, counts, color=['skyblue', 'orange', 'green', 'red'])
    st.pyplot(fig)

    # Motif search
    st.subheader("Motif Search (e.g., TATA)")
    motif = st.text_input("Enter motif", "TATA").upper()
    positions = [i for i in range(len(seq) - len(motif) + 1) if seq[i:i+len(motif)] == motif]
    st.write(f"Motif '{motif}' found at positions: {positions}")

    # Mutation simulation
    st.subheader("Simulate Mutation")
    mutation_pos = st.slider("Choose mutation position", 0, len(seq)-1, 50)
    original_base = seq[mutation_pos]
    mutated_base = random.choice([b for b in "ATGC" if b != original_base])
    mutated_seq = seq[:mutation_pos] + mutated_base + seq[mutation_pos+1:]
    st.code(f"Original: ...{seq[mutation_pos-5:mutation_pos+6]}...")
    st.code(f"Mutated : ...{mutated_seq[mutation_pos-5:mutation_pos+6]}...")
    st.warning(f"Mutation at position {mutation_pos}: {original_base} → {mutated_base}")

    # Classification rule (simple based on GC content)
    if gc_content > 55:
        label = "Likely Probiotic (High GC%)"
        st.success(label)
    elif gc_content < 40:
        label = "Possibly Pathogenic (Low GC%)"
        st.error(label)
    else:
        label = "Intermediate Classification"
        st.info(label)
