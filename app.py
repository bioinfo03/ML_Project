import os
import re
from io import BytesIO

import requests
from flask import Flask, render_template, request, send_file
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


app = Flask(__name__)

# --------------------------------------------------
# 1. Kyte-Doolittle hydrophobicity scale
# --------------------------------------------------

KD = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8,
    "C": 2.5, "M": 1.9, "A": 1.8, "G": -0.4,
    "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3,
    "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5,
    "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5
}

# Three wildcard positions can contain only these residues.
ALLOWED_AA = "AVLIMPFYWSTNCHQG"
MOTIF_PATTERN = rf"AL[{ALLOWED_AA}]{{3}}LW"


# --------------------------------------------------
# 2. Calculate sliding-window KD scores
# --------------------------------------------------

def calculate_kd_scores(sequence, window_size=19):
    if window_size <= 0:
        raise ValueError("Window size must be positive.")

    if len(sequence) < window_size:
        return []

    scores = []

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]

        if any(aa not in KD for aa in window):
            scores.append(None)
        else:
            score = sum(KD[aa] for aa in window) / window_size
            scores.append(round(score, 3))

    return scores


# --------------------------------------------------
# 3. Find all motif matches
# --------------------------------------------------

def find_motifs(sequence):
    matches = []

    for match in re.finditer(MOTIF_PATTERN, sequence):
        matches.append({
            "motif": match.group(),
            "start": match.start() + 1,
            "end": match.end()
        })

    return matches


# --------------------------------------------------
# 4. Keep motifs overlapping high-scoring windows
# --------------------------------------------------

def find_qualifying_motifs(
    sequence, window_size=19, threshold=1.6
):
    scores = calculate_kd_scores(sequence, window_size)
    high_positions = set()

    for i, score in enumerate(scores):
        if score is not None and score >= threshold:
            high_positions.update(range(i, i + window_size))

    qualifying = []

    for hit in find_motifs(sequence):
        start = hit["start"] - 1
        end = hit["end"]

        if any(pos in high_positions for pos in range(start, end)):
            qualifying.append(hit)

    return qualifying, scores, high_positions


# --------------------------------------------------
# 5. Retrieve a protein sequence from NCBI
# --------------------------------------------------

def fetch_protein_from_ncbi(accession):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/"
        "entrez/eutils/efetch.fcgi"
    )

    params = {
        "db": "protein",
        "id": accession,
        "rettype": "gp",
        "retmode": "text"
    }

    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()

    record = response.text

    if "ORIGIN" not in record:
        raise ValueError(
            "No protein record found. Check the accession number."
        )

    origin = record.split("ORIGIN", 1)[1].split("//", 1)[0]
    sequence = "".join(re.findall(r"[A-Za-z]+", origin)).upper()

    if not sequence:
        raise ValueError("NCBI returned an empty sequence.")

    if any(aa not in KD for aa in sequence):
        raise ValueError(
            "The protein contains non-standard amino acids "
            "that this version does not support."
        )

    return sequence, record


# --------------------------------------------------
# 6. Create the hydrophobicity plot
# --------------------------------------------------

def create_plot(scores, window_size, threshold, qualifying):
    os.makedirs("static", exist_ok=True)

    positions = [
        i + (window_size + 1) / 2
        for i, score in enumerate(scores)
        if score is not None
    ]

    values = [
        score for score in scores
        if score is not None
    ]

    plt.figure(figsize=(12, 4))

    if values:
        plt.plot(positions, values, linewidth=1.3, label="KD score")

    plt.axhline(
        y=threshold,
        linestyle="--",
        label=f"Threshold = {threshold}"
    )

    for hit in qualifying:
        plt.axvspan(
            hit["start"], hit["end"],
            alpha=0.2
        )

    plt.xlabel("Amino acid position (window center)")
    plt.ylabel("Average Kyte-Doolittle score")
    plt.title("Protein Hydrophobicity Profile")
    plt.legend()
    plt.tight_layout()

    plt.savefig("static/kd_plot.png", dpi=150)
    plt.close()


# --------------------------------------------------
# 7. Home page
# --------------------------------------------------

@app.route("/")
def index():
    return render_template("index.html")


# --------------------------------------------------
# 8. Analyze protein
# --------------------------------------------------

@app.route("/analyze", methods=["POST"])
def analyze():
    accession = request.form.get("accession", "").strip()
    pasted_sequence = request.form.get("sequence", "")

    try:
        window_size = int(request.form.get("window_size", 19))
        threshold = float(request.form.get("threshold", 1.6))

        if window_size <= 0:
            raise ValueError("Window size must be positive.")

        if accession:
            sequence, _ = fetch_protein_from_ncbi(accession)
            source = "NCBI"
        else:
            sequence = "".join(pasted_sequence.split()).upper()
            source = "Pasted sequence"

        if not sequence:
            raise ValueError(
                "Enter an accession number or paste a sequence."
            )

        if any(aa not in KD for aa in sequence):
            raise ValueError(
                "Sequence contains invalid or non-standard amino acids."
            )

        qualifying, scores, high_positions = find_qualifying_motifs(
            sequence, window_size, threshold
        )

        create_plot(scores, window_size, threshold, qualifying)

        # Colour qualifying motif residues blue.
        # Colour other high-window residues red.
        motif_positions = set()

        for hit in qualifying:
            motif_positions.update(
                range(hit["start"] - 1, hit["end"])
            )

        residues = []

        for i, aa in enumerate(sequence):
            if i in motif_positions:
                class_name = "motif"
            elif i in high_positions:
                class_name = "high"
            else:
                class_name = ""

            residues.append({
                "aa": aa,
                "position": i + 1,
                "class_name": class_name
            })

        ncbi_link = None

        if accession:
            ncbi_link = (
                "https://www.ncbi.nlm.nih.gov/protein/?term="
                + requests.utils.quote(accession)
            )

        return render_template(
            "result.html",
            source=source,
            accession=accession,
            sequence=sequence,
            window_size=window_size,
            threshold=threshold,
            qualifying=qualifying,
            residues=residues,
            ncbi_link=ncbi_link
        )

    except (
        ValueError,
        requests.RequestException
    ) as error:
        return render_template(
            "index.html",
            error=str(error)
        ), 400


# --------------------------------------------------
# 9. Export the NCBI GenBank record to Excel
# --------------------------------------------------

@app.route("/export/<accession>")
def export_protein_excel(accession):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/"
        "entrez/eutils/efetch.fcgi"
    )

    params = {
        "db": "protein",
        "id": accession,
        "rettype": "gp",
        "retmode": "text"
    }

    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        record = response.text

        if "ORIGIN" not in record:
            return "No valid protein record found.", 404

        # Parse the record into sections.
        sections = []
        current_section = None
        current_lines = []

        def save_section():
            if current_section is not None:
                content = "\n".join(current_lines).strip()

                if content:
                    sections.append((current_section, content))

        for line in record.splitlines():
            stripped = line.strip()

            if line.startswith("LOCUS"):
                save_section()
                current_section = "LOCUS"
                current_lines = [line]

            elif line.startswith("DEFINITION"):
                save_section()
                current_section = "DEFINITION"
                current_lines = [line]

            elif line.startswith("ACCESSION"):
                save_section()
                current_section = "ACCESSION"
                current_lines = [line]

            elif line.startswith("VERSION"):
                save_section()
                current_section = "VERSION"
                current_lines = [line]

            elif line.startswith("REFERENCE"):
                save_section()
                current_section = "REFERENCE"
                current_lines = [line]

            elif line.startswith("FEATURES"):
                save_section()
                current_section = "FEATURES"
                current_lines = [line]

            elif line.startswith("ORIGIN"):
                save_section()
                current_section = "ORIGIN"
                current_lines = [line]

            elif stripped == "//":
                save_section()
                current_section = None
                current_lines = []

            elif current_section is not None:
                current_lines.append(line)

        save_section()

        workbook = Workbook()
        worksheet = workbook.active
        worksheet.title = "Protein Details"

        worksheet.append(["Section", "Details"])

        for section, content in sections:
            worksheet.append([section, content])

        for cell in worksheet[1]:
            cell.font = Font(bold=True)

        worksheet.column_dimensions["A"].width = 22
        worksheet.column_dimensions["B"].width = 100

        for row in worksheet.iter_rows(min_row=2):
            row[1].alignment = Alignment(
                wrap_text=True,
                vertical="top"
            )

        output = BytesIO()
        workbook.save(output)
        output.seek(0)

        return send_file(
            output,
            as_attachment=True,
            download_name=f"{accession}_protein_details.xlsx",
            mimetype=(
                "application/vnd.openxmlformats-officedocument."
                "spreadsheetml.sheet"
            )
        )

    except requests.RequestException as error:
        return f"NCBI request failed: {error}", 502


# --------------------------------------------------
# 10. Start the Flask application
# --------------------------------------------------

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=False)        

