import os
import requests
import re

# 1. Install Grobid if not installed yet:
# git clone https://github.com/kermitt2/grobid.git
# cd grobid
# ./gradlew clean install
# ./gradlew run
# 2. run GROBID in terminal:
# cd /Users/Iris/grobid
# ./gradlew run or :grobid-service:ru
# 3. Check before running the script: http://localhost:8070/api/isalive -> returns true
# 4. CTRL+C to end connection

# GROBID server URL
GROBID_URL = "http://localhost:8070/api/processHeaderDocument"

# Paths
input_folder = "/Users/Iris/Desktop/BachelorProject/RNA-Review"  # Update this
output_bib_file = "/Users/Iris/Desktop/BachelorProject/RNA-Review/output.bib"

# Function to generate custom BibTeX key
def generate_bibtex_key(bibtex_entry):
    author_match = re.search(r'author\s*=\s*{([^}]+)}', bibtex_entry, re.IGNORECASE)
    year_match = re.search(r'year\s*=\s*{(\d{4})}', bibtex_entry, re.IGNORECASE)
    title_match = re.search(r'title\s*=\s*{([^}]+)}', bibtex_entry, re.IGNORECASE)

    if not (author_match and year_match and title_match):
        return None

    authors = author_match.group(1).split(" and ")
    first_author = authors[0].strip().split(",")[0].lower()

    title_words = re.findall(r'\b[A-Za-z]+\b', title_match.group(1))
    first_word = title_words[0].capitalize() if title_words else "Title"

    year = year_match.group(1)
    return f"{first_author}{year}{first_word}"

# 🔍 Gather all PDFs
pdf_files = []
for root, dirs, files in os.walk(input_folder):
    for file in files:
        if file.lower().endswith(".pdf"):
            pdf_files.append(os.path.join(root, file))

print(f"📄 Found {len(pdf_files)} PDF files.")

# 🧠 Process PDFs
with open(output_bib_file, "w", encoding="utf-8") as bib_out:
    for pdf_path in pdf_files:
        try:
            with open(pdf_path, 'rb') as pdf_file:
                response = requests.post(
                    GROBID_URL,
                    files={'input': pdf_file},
                    data={'consolidateHeader': '1'},
                    headers={"Accept": "application/x-bibtex"}
                )

            if response.status_code == 200:
                bibtex_entry = response.text.strip()

                if bibtex_entry:
                    # Replace BibTeX key
                    bib_key = generate_bibtex_key(bibtex_entry)
                    if bib_key:
                        bibtex_entry = re.sub(r"@(\w+)\{[^,]+,", f"@article{{{bib_key},", bibtex_entry, count=1)

                    bib_out.write(bibtex_entry + "\n\n")
                    print(f"✅ Processed: {os.path.basename(pdf_path)}")
                else:
                    print(f"⚠️ No BibTeX extracted: {os.path.basename(pdf_path)}")
            else:
                print(f"❌ Failed: {os.path.basename(pdf_path)} — Status {response.status_code}")

        except Exception as e:
            print(f"❌ Error with {pdf_path}: {e}")

print(f"\n📚 Done! BibTeX entries saved to: {output_bib_file}")
