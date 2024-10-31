#!/vol/patric3/cli/runtime/bin/python3
from flask import Flask, render_template, request, redirect, send_from_directory
from werkzeug.utils import secure_filename
from waitress import serve
import time
import uuid
import os

print("Starting Hammer Server. Starting timer.")
start = time.time()

application = Flask(__name__)

#Dict of genome id to (genome name, taxonomy)
genome_names = {}
#Dict of hammer to repgen genome id
hammers = {}
#upload Folder
upload_folder = "/disks/tmp/"
num_of_results = 5

@application.errorhandler(404)
def page_not_found(e):
    '''Error message for live production debugging'''
    print(request.url)
    print(e)
    return render_template('404.html'), 404

def hammer_initialization(hammers, repdict):
    '''Intakes the empty dictionaries of hammers and the list of repgen genome names and parses the respective files.
    returns a dictionary of genomeid keys to hammer values, then returns a dictionary of repgen ids to genome names.
    '''
    with open("data/hammers200.worthless20.tbl") as hammer_file:
        lines = hammer_file.readlines()
        for line in lines:
            if "hammer" not in line:
                line = line.strip("\n")
                row = line.split("\t")
                fig_id = row[1].lstrip("fig|")
                fig_id = fig_id.split(".")
                genome_id = fig_id[0] + "." + fig_id[1]
                hammers[row[0]] = genome_id
        print('hammers have been loaded')
    with open("data/rep200.tax.tbl") as repgen_file:
        lines = repgen_file.readlines()
        for line in lines:
            line = line.strip("\n")
            row = line.split("\t")
            repdict[row[0]] = (row[1], row[2])
    return hammers, repdict

hammers, genome_names = hammer_initialization(hammers=hammers, repdict=genome_names)
print(f'This has been running for {round(time.time()-start, 2)} seconds')
@application.route('/services/hammers')
def index():
    print(request.base_url)
    '''Renders the homepage of the website'''
    return render_template('index.html', hammers=hammers, genome_names=genome_names)

def findhammer(fasta_file):
    '''This takes in a fasta file and returns all of the hammers that have hit those contigs.'''
    start = time.time()
    repgen_hits = {}
    whole_contig = ""
    contig_count = 0
    with open(fasta_file) as data_file:
        lines = data_file.readlines()
        for line in lines:
            contig_piece = line.strip("\n")
            if ">" not in contig_piece:
                whole_contig += contig_piece
            #Dict of repgenid to count
            elif ">" in contig_piece and whole_contig != '':
                contig_count += 1
                reverse_contig = whole_contig[::-1]
                repgen_hits.update(hammer_strike(reverse_contig, repgen_hits, fasta_file))
                repgen_hits.update(hammer_strike(whole_contig, repgen_hits, fasta_file))
                whole_contig = ''
        repgen_hits.update(hammer_strike(whole_contig, repgen_hits, fasta_file))
        results = []
        for repgen in repgen_hits.keys():
            try:
                taxonomy = [genome_names[repgen][1]]
                results.append([repgen, genome_names[repgen][0], taxonomy, repgen_hits[repgen]])
            except KeyError:
                print(f"Could not find name in rep100.list.tbl for {repgen}. Please add to the list.")
        print(f'This has been running for {round(time.time() - start, 2)} seconds')
    return results

def hammer_strike(contig, repgen_hits_dict, filename):
    '''inputs a contig and the repgen hits dictionary and adds hammer hits to the repgen dictionary.'''
    hammers_hit = []
    for i in range(len(contig) - 19):
        twenty_mer = contig[i:i + 20].lower()
        if twenty_mer in hammers.keys():
            repgen_hit = hammers[twenty_mer]
            hammers_hit.append(twenty_mer)
            if repgen_hit not in repgen_hits_dict.keys():
                repgen_hits_dict[repgen_hit] = 1
            else:
                repgen_hits_dict[repgen_hit] += 1
    #Used to print the actual 20mers that hit the genome
    print_hammers_hit(hammers_hit, filename)
    return repgen_hits_dict

def print_hammers_hit(hammers_hit_list, filename):
    '''Takes the results of the hammers hit and outputs a file that shows every hammer that was hit.'''
    if ".fna" in filename:
        filename = filename.rstrip(".fna")
    elif ".fasta" in filename:
        filename = filename.rstrip(".fasta")
    with open(f"hammers_hit.{filename}.tbl", "a+") as hammers_output:
        for hammer in hammers_hit_list:
            hammers_output.write(f"{hammer}\n")

def create_download_table(results_list, file_name):
    path = upload_folder + file_name
    with open(path, "w") as download_table:
        download_table.write("Genome ID\tGenome Name\tTaxonomy\tNumber of Hammer Hits\n")
        for list in results_list:
            download_table.write(f"{list[0]}\t{list[1]}\t")
            for i in list[2]:
                download_table.write(i + " ")
            download_table.write(f"\t{list[3]}\n")
    return send_from_directory(directory=path, filename=file_name)

@application.route("/results", methods=['GET', 'POST'])
def upload_file():
    '''Saves the uploaded fasta file uploaded to the website to create a clean path for input later.
    Also renders the results template for the webpage.
    Sorts the hits to have the most hammers hit on top
    Downloads the results as a txt file.
    '''
    if request.method == 'POST':
        if 'fasta_file' not in request.files:
            return redirect(request.url)
        file = request.files['fasta_file']
        num_of_results = int(request.form.get("num_results"))
        if num_of_results == '':
            num_of_results = 5
        if file.filename == '':
            return redirect(request.url)
        if file:
            file_name = secure_filename(file.filename)
            file_path = upload_folder + file_name + str(uuid.uuid4())
            file.save(file_path)
            print("Saved the file " + file_path)
            results = findhammer(file_path)
            #Let's sort!
            top100_counter = [0, 0, 0, 0, 0]
            final_results_list = []
            for result in results:
                for i in range(len(top100_counter)):
                    if result[3] > top100_counter[i]:
                        top100_counter.append(result[3])
                        top100_counter.sort(reverse=True)
                        top100_counter = top100_counter[0:num_of_results]
                        break
            for i in range(len(top100_counter)):
                for result in results:
                    if top100_counter[i] == result[3]:
                        if result not in final_results_list:
                            final_results_list.append(result)
                        break
            create_download_table(final_results_list, file_name)
            os.remove(file_path)
            return render_template("results.html", results=final_results_list, filename=file_path)
    else:
        return render_template("index.html", hammers=hammers, genome_names=genome_names)


if __name__ == "__main__":
    '''Happy Flask'''
    serve(application, host='0.0.0.0', port='5000')