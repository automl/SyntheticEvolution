import requests
import json
import time
import os
import urllib

# https://colab.research.google.com/drive/1UKJ-5fY254ALs_nNBvggpwc3raMkNtzo?usp=sharing#scrollTo=UnU26tQ2bBFt
def extract_interpro_entries(uniprot_id):
    UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb"
    url =  '{}/{}.json'.format(UNIPROT_API_URL, uniprot_id)
    uniprot_results = json.load(urllib.request.urlopen(url))
    # sequence = uniprot_results['sequence']['value']

    interpro_entries = []
    if 'uniProtKBCrossReferences' in uniprot_results and isinstance(uniprot_results['uniProtKBCrossReferences'], list):
        for reference in uniprot_results['uniProtKBCrossReferences']:
            if reference.get('database') == 'InterPro':
                # Iterate over properties to find the relevant 'EntryName'
                for prop in reference.get('properties', []):
                    if prop.get('key') == 'EntryName':
                        entry_name = prop.get('value')
                        if entry_name:
                            interpro_entries.append({
                                # "database": reference.get('database'),
                                # "id": reference.get('id'),
                                "FamilyName": entry_name
                            })

    return interpro_entries

# Example:
uniprot_id = 'Q15650'
print(extract_interpro_entries(uniprot_id))

# https://www.ebi.ac.uk/interpro/api/static_files/swagger/#/Entry/get_entry
# https://gustavo-salazar.github.io/ProteinFamiliesTalks/InterProAPI.html#/25
# https://www.youtube.com/watch?v=8i3yszapfRo&t=1559s
def get_protein_family_from_sequence(sequence, pdb_id):
    # if not os.path.exists(output_folder):
    #     os.makedirs(output_folder)
    url_submit = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    data = {'email': 'iris.kramberger@gmail.com', 'sequence': sequence}

    response_submit = requests.post(url_submit, headers=headers, data=data)
    if response_submit.status_code != 200:
        raise Exception(f"Failed to submit sequence: {response_submit.status_code}, {response_submit.text}")
    job_id = response_submit.text
    print(f"Job submitted with ID: {job_id}")

    url_status = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}"

    while True:
        status = requests.get(url_status).text
        print(f"Job status: {status}")
        if status == 'FINISHED':
            break
        elif status in ['RUNNING', 'PENDING']:
            time.sleep(30)

    url_results = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/json"
    response_results = requests.get(url_results)
    if response_results.status_code != 200:
        raise Exception(f"Failed to retrieve results: {response_results.status_code}, {response_results.text}")

    payload = response_results.json()

    protein_families = []

    for result in payload.get('results', []):
        for match in result.get('matches', []):
            signature = match.get('signature', {})
            entry = signature.get('entry', {})
            if entry and entry.get('type') == 'FAMILY':
                family_name = entry.get('name')
                print(family_name)
                protein_families.append(family_name)

    if not protein_families:
        output_folder = '/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data'
        output_filename = os.path.join(output_folder, f"{pdb_id}.json")
        with open(output_filename, 'w') as json_file:
            json.dump(payload, json_file, indent=4)
        return None

    return protein_families


    # To save raw JSON file if no detected family_name:
    # output_folder = '/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data'
    # output_filename = os.path.join(output_folder, f"{job_id}.json")
    # with open(output_filename, 'w') as json_file:
    #     json.dump(payload, json_file, indent=4)

# Example
# sequence = "FKPPPRPDFGTSGRTIKLQANFFEMDIPKIDIYHYELDIKPEKCPRRVNREIVEHMVQHFKTQIFGDRKPVFDGRKNLYTAMPLPIGRDKVELEVTLPGEKDRIFKVSIKWVSCVSLQALHDALSGRLPSVPFETIQALDVVMRHLPSMRYTPVGRSFFTASEGCSNPLGGGREVWFGFHQSVRPSLWKMMLNIDVSATAFYKAQPVIEFVCEVLDFKSIEEQQKPLTDSQRVKFTKEIKGLKVEITHCGKRKYRVCNVTRRPASHQTFPLQQESGQTVECTVAQYFKDRHKLVLRYPHLPCLQVGQEQKHTYLPLEVCNIVAGQRCIKKLTDNQTSTMIRATARSAPDRQEEISKLMRSADFNTDPYVREFGIMVKDEMTDVTGRVLQPPSILYGGRNKAIATPVQGVWDMRNKQFHTGIEIKVWAIACFAPQRQCTEVHLKSFTEQLRKISRDAGMPIQGQPCFCKYAQGADSVEPMFRHLKNTYAGLQLVVVILPGKTPVYAEVKRVGDTVLGMATQCVQMKNVQRTTPQTLSNLCLKINVKLGGVNNILLPQGRPPVFQQPVIFLGADVTHPPGKKPSIAAVVGSMDAHPNRYCATVRVQQHRQEIIQDLAAMVRELLIQFYKSTRFKPTRIIFYRDGVSEGQFQQVLHHELLAIREACIKLEKDYQPGITFIVVQKRHHTRLFCTDKNERVGKSGNIPAGTTVDTKITHPTEFDFYLCSHAGIQGTSRPSHYHVLWDDNRFSSDELQILTYQLCHTYVRCTRSVSIPAPAYYAHLVAFRARYHLHQALAKAVQVHQDTLRTMYFA"
# protein_families = get_protein_family_from_sequence(sequence, "c150")
# print(protein_families)

# def extract_protein_family_from_json(json_filename):
#     with open(json_filename, 'r') as json_file:
#         payload = json.load(json_file)
#
#     protein_families = []
#
#     for result in payload.get('results', []):
#         for match in result.get('matches', []):
#             signature = match.get('signature', {})
#             entry = signature.get('entry', {})
#             if entry and entry.get('type') == 'FAMILY':
#                 family_name = entry.get('name')
#                 print(family_name)
#                 protein_families.append(family_name)
#
#     return protein_families
#
# def process_json_folder(folder_path):
#     for filename in os.listdir(folder_path):
#         if filename.endswith('.json'):
#             json_filepath = os.path.join(folder_path, filename)
#             families = extract_protein_family_from_json(json_filepath)
#             print(f"File: {filename}, Protein Families: {families}")
#
# folder_path = 'results_folder'
# process_json_folder('/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/data')
