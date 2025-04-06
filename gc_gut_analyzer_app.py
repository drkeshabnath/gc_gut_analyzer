
# gc_gut_analyzer_app.py

import streamlit as st
from Bio import SeqIO, Seq
import matplotlib.pyplot as plt
import random
import io
import pandas as pd
from collections import Counter

st.set_page_config(page_title="GC-Gut Analyzer", layout="centered")
st.title("🧬 GC-Gut Analyzer - Microbiome Insight Tool [Dr. Keshab Nath]")

uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta", "fa"])

def sliding_gc(seq, window=50):
    return [100 * (s.count("G") + s.count("C")) / len(s)
            for s in [seq[i:i+window] for i in range(len(seq)-window + 1)]]

if uploaded_file is not None:
    text_io = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
    record = next(SeqIO.parse(text_io, "fasta"))
    seq = str(record.seq).upper()
    st.success(f"Sequence loaded: {record.id} ({len(seq)} bases)")

    # Sequence Stats
    gc_content = 100 * (seq.count("G") + seq.count("C")) / len(seq)
    at_content = 100 - gc_content
    st.metric("GC Content (%)", f"{gc_content:.2f}%")
    st.write(f"🧬 Sequence Length: {len(seq)} bp")
    st.write(f"🧪 AT Content: {at_content:.2f}%")

    # GC Sliding Window Plot
    st.subheader("📈 GC Content Sliding Window")
    window = st.slider("Window Size", 10, 200, 50)
    gc_vals = sliding_gc(seq, window)
    fig, ax = plt.subplots()
    ax.plot(range(len(gc_vals)), gc_vals, color="green")
    ax.set_title("GC% Sliding Window")
    ax.set_xlabel("Position")
    ax.set_ylabel("GC%")
    st.pyplot(fig)

    # Base Frequency
    st.subheader("🔬 Nucleotide Frequency")
    bases = ['A', 'T', 'G', 'C']
    counts = [seq.count(b) for b in bases]
    fig2, ax2 = plt.subplots()
    ax2.bar(bases, counts, color=['skyblue', 'orange', 'green', 'red'])
    st.pyplot(fig2)

    # Motif search
    st.subheader("🧪 Motif Search")
    motif = st.text_input("Enter motif (e.g., TATA, ATG)", "TATA").upper()
    positions = [i for i in range(len(seq) - len(motif) + 1) if seq[i:i+len(motif)] == motif]
    st.write(f"Found {len(positions)} match(es) for '{motif}': {positions}")

    # Mutation simulation
    st.subheader("🧬 Simulate Mutation")
    mutation_pos = st.slider("Choose mutation position", 0, len(seq)-1, 50)
    original_base = seq[mutation_pos]
    mutated_base = random.choice([b for b in "ATGC" if b != original_base])
    mutated_seq = seq[:mutation_pos] + mutated_base + seq[mutation_pos+1:]
    st.code(f"Original : ...{seq[mutation_pos-5:mutation_pos+6]}...")
    st.code(f"Mutated  : ...{mutated_seq[mutation_pos-5:mutation_pos+6]}...")
    st.warning(f"Mutation: {original_base} → {mutated_base} at position {mutation_pos}")

    # Classification
    st.subheader("🧠 Classification Based on GC Content")
    if gc_content > 55:
        label = "Likely Probiotic (High GC%)"
        st.success(label)
    elif gc_content < 40:
        label = "Possibly Pathogenic (Low GC%)"
        st.error(label)
    else:
        label = "Intermediate Classification"
        st.info(label)

    # Reverse Complement
    st.subheader("🔁 Reverse Complement")
    rev_comp = str(Seq.Seq(seq).reverse_complement())
    st.code(f"{rev_comp[:100]}...")

    # Codon usage (first reading frame)
    st.subheader("🧬 Codon Usage (First Frame)")
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    codon_counts = Counter(codons)
    df_codons = pd.DataFrame(codon_counts.items(), columns=["Codon", "Count"]).sort_values(by="Count", ascending=False)
    st.dataframe(df_codons)

    # Amino Acid Translation
    st.subheader("🧫 Protein Translation (First Frame)")
    protein = str(Seq.Seq(seq).translate(to_stop=True))
    st.code(protein[:100] + ("..." if len(protein) > 100 else ""))

    # Summary Panel
    st.subheader("📊 Summary of Results")
    st.markdown(f'''
- **Sequence Length:** {len(seq)} bp  
- **GC Content:** {gc_content:.2f}%  
- **AT Content:** {at_content:.2f}%  
- **Motif '{motif}' Matches:** {len(positions)}  
- **Codons Counted:** {len(df_codons)}  
- **Translated Protein Length:** {len(protein)} amino acids  
''')
