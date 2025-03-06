import streamlit as st
import pandas as pd
from Bio import SeqIO
from io import StringIO
import re
import plotly.express as px
from concurrent.futures import ProcessPoolExecutor
from reportlab.pdfgen import canvas
from Bio.Seq import Seq

execution_metadata = {
    "execution_count": None,
    "status": "running",
    "error": None
}

st.title('Advanced DNA Promoter Prediction and Non-B DNA Motif Analysis')
st.write('Upload multiple FASTA files to analyze DNA motifs, predict promoter regions, and visualize results.')
st.image('https://github.com/chandrika180898/nnbdna/blob/main/utr_image.jpg', caption='DNA Structure')


uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)

# Updated motifs dictionary to include H-DNA and R-Loop motifs
motifs = {
    "Slipped DNA": re.compile(r'([ATGC]{2,6})\1{1,}'),
    "Z-DNA": re.compile(r'(CG){6,}'),
    "Short Tandem Repeat": re.compile(r'([ATGC]{2,6})\1{2,}'),
    "I-Motif": re.compile(r'((C[A,T]C){3,})'),
    "R-Loop": re.compile(r'(A{4,}[CG]{2,}A{4,})'),
    "Cruciform": re.compile(r'([ATGC]{4,})\1{2,}'),
    "G-Quadruplex": re.compile(r'(G{3,}[ATGC]{1,5}G{3,}[ATGC]{1,5}G{3,}[ATGC]{1,5}G{3,})'),
    "Hairpin": re.compile(r'([ATGC]{4,})\1{1,}'),
    "Triplex": re.compile(r'(A{3,}[ATGC]{1,}A{3,})'),
    "H-DNA": re.compile(r'([AG]{4,}[CT]{4,}[AG]{4,})'),
    "Triplex-forming oligonucleotide (TFO)": re.compile(r'([GATC]{6,}[AG]{4,}[CT]{4,})')
}

# Function to find inverted repeats
def find_inverted_repeats(sequence):
    inverted_repeat_results = []
    pattern = r'([ATGC]{3,})[ATGC]{0,10}([ATGC]{3,})'
    for match in re.finditer(pattern, str(sequence)):
        part1 = match.group(1)
        part2 = match.group(2)[::-1]
        if part1 == part2:
            inverted_repeat_results.append({
                "Motif": "Inverted Repeat",
                "Start": match.start() + 1,
                "End": match.end(),
                "Matched Sequence": sequence[match.start():match.end()]
            })
    return inverted_repeat_results

# Function to find motifs
def find_motifs(sequence):
    results = []
    for motif_name, motif_pattern in motifs.items():
        for match in motif_pattern.finditer(str(sequence)):
            results.append({
                "Motif": motif_name,
                "Start": match.start() + 1,
                "End": match.end(),
                "Matched Sequence": sequence[match.start():match.end()]
            })
    results.extend(find_inverted_repeats(sequence))
    return results

# Parallel sequence analysis
def analyze_sequences_parallel(sequences):
    data = []
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(find_motifs, [record.seq for record in sequences]))
        for record, motif_results in zip(sequences, results):
            for motif in motif_results:
                data.append({
                    "Sequence ID": record.id,
                    **motif,
                    "Length": len(record.seq)
                })
    return pd.DataFrame(data)

# Visualization
def visualize_motifs(df):
    fig = px.scatter(df, x='Start', y='Sequence ID', color='Motif',
                     hover_data=['Matched Sequence'],
                     title="Motif Distribution Across Sequences")
    st.plotly_chart(fig)

# PDF generation
def generate_pdf(df):
    c = canvas.Canvas("motif_report.pdf")
    c.drawString(100, 800, "DNA Motif Analysis Report")
    y = 780
    for i, row in df.iterrows():
        c.drawString(100, y, f"{row['Sequence ID']} | {row['Motif']} | Start: {row['Start']} | End: {row['End']}")
        y -= 20
    c.save()

# Process uploaded files
def process_uploaded_files(uploaded_files):
    all_results = pd.DataFrame()
    for uploaded_file in uploaded_files:
        fasta_sequences = list(SeqIO.parse(StringIO(uploaded_file.getvalue().decode('utf-8')), 'fasta'))
        results_df = analyze_sequences_parallel(fasta_sequences)
        all_results = pd.concat([all_results, results_df], ignore_index=True)
    return all_results

if uploaded_files:
    try:
        results_df = process_uploaded_files(uploaded_files)
        
        if 'Matched Sequence' in results_df.columns:
            results_df['Matched Sequence'] = results_df['Matched Sequence'].apply(lambda x: str(x) if isinstance(x, Seq) else x)
        else:
            st.error("No motifs found or the 'Matched Sequence' column is missing!")

        st.write("### Motif Analysis Results")
        st.dataframe(results_df)
        visualize_motifs(results_df)
        
        if st.button("Generate PDF Report"):
            generate_pdf(results_df)
            with open("motif_report.pdf", "rb") as pdf:
                st.download_button(
                    "Download PDF Report", pdf, file_name="motif_analysis_report.pdf")

        csv = results_df.to_csv(index=False)
        st.download_button(
            label="Download CSV",
            data=csv,
            file_name="motif_analysis_results.csv",
            mime="text/csv"
        )
        execution_metadata["status"] = "completed"
    except Exception as e:
        execution_metadata["status"] = "error"
        execution_metadata["error"] = str(e)
        st.error(f"An error occurred: {e}")
