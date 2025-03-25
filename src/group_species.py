from collections import defaultdict

input_file = "../results/combined_summary.txt"
output_file = "../results/combined_summary_grouped.txt"

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    current_report = None
    species_data = defaultdict(lambda: {'reads': 0, 'percent': 0.0})
    total_reads = 0

    for line in infile:
        line = line.strip()

        if line.startswith("=== Summary for"):
            # if we have previous data, write it
            if current_report:
                outfile.write(f"=== Combined Summary for {current_report} ===\n")
                for species, data in species_data.items():
                    percent = round((data['reads'] / total_reads) * 100, 2) if total_reads else 0
                    outfile.write(f"{species}\t{data['reads']}\t{percent}%\n")
                outfile.write("\n")

            current_report = line.split("Summary for ")[-1].replace("===", "").strip()
            species_data = defaultdict(lambda: {'reads': 0, 'percent': 0.0})
            total_reads = 0

        elif line:  # non-empty line with species data
            try:
                species, count, percent = line.split('\t')
                count = int(count)
                species_data[species]['reads'] += count
                total_reads += count
            except ValueError:
                continue  # skip malformed lines

    if current_report:
        outfile.write(f"=== Combined Summary for {current_report} ===\n")
        for species, data in species_data.items():
            percent = round((data['reads'] / total_reads) * 100, 2) if total_reads else 0
            outfile.write(f"{species}\t{data['reads']}\t{percent}%\n")
        outfile.write("\n")


input_file = "../results/combined_summary_grouped.txt"

with open(input_file, 'r') as f:
    current_report = None
    species_list = []
    top_species = ""
    top_percent = 0.0

    for line in f:
        line = line.strip()

        if line.startswith("=== Combined Summary for"):
            # if we've finished a block, print the result
            if current_report:
                print(f"{current_report}")
                print(f"  Total species: {len(species_list)}")
                print(f"  Top species: {top_species} ({top_percent}%)\n")

            # start new report
            current_report = line.replace("=== Combined Summary for ", "").replace("===", "").strip()
            species_list = []
            top_species = ""
            top_percent = 0.0

        elif line:
            try:
                species, count, percent_str = line.split('\t')
                percent = float(percent_str.strip('%'))
                species_list.append(species)

                if percent > top_percent:
                    top_species = species
                    top_percent = percent
            except ValueError:
                continue

    # final report
    if current_report:
        print(f"{current_report}")
        print(f"  Total species: {len(species_list)}")
        print(f"  Top species: {top_species} ({top_percent}%)\n")
