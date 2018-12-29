import os
from pathlib import Path
from subprocess import call
import yaml
configfile: "directories.yaml"

"""
This script trims and QCs fastq files with special modifications for Hiseq.

Basically, the only difference is that this version does not subsample the read files
and does not blast anything.

In addition, this version does not do custom adapter trimming in order to halve the total file output size
"""

# The yaml file should look like this below. Just an object called "directories"
#  and a list of directories in which the pipeline should look for things.
#
#  directories:
#    - "/data/raw/genomic_reads/180315_M00160_0073_000000000-BNBV5/"
#    - "/data/raw/genomic_reads/180123_M00160_0067_000000000-BKDV3/"
#    - "/data/raw/genomic_reads/171018_M00160_0050_000000000-BFNHV/"

#   vvvv  LOOK HERE vvv
#
#  your path to the directory containing trimmomatic-0.35.jar and adapters/TruSeq3-PE-2.fa
#   goes here
trimmomatic = config["trimmomatic"]
trimmomatic = "/usr/local/bin/Trimmomatic-0.35"
bbmergepath = config["bbmergepath"]
bbmergepath = "/usr/local/bin/bbmap/bbmerge.sh"
maxthreads = config["maxthreads"]
maxthreads = 90
#
#
#   ^^^^  LOOK HERE ^^^

def is_fastq(filepath):
    suffix_set = set(Path(filepath).suffixes)
    fastqs = set([".fastq", ".fq"])
    fastq_intersection = suffix_set & fastqs
    if len(fastq_intersection) > 0:
        return True
    else:
        return False

#now we must populate the sample list
config["samples"] = {}
config["adapter_pairs"] = {}
for this_d in config["directories"]:
    these_files = os.listdir(this_d)
    fastqs = [filepath for filepath in these_files if is_fastq(filepath)]
    samplelist = set()
    for filepath in fastqs:
        split = filepath.split("_")
        stop = 1
        if "R1" in split:
            stop = split.index("R1")
        elif "R2" in split:
            stop = split.index("R2")
        temp_samplename = "_".join(split[0:stop])
        add_sample = True
        if "ignore" not in config:
            config["ignore"] = []
        for check in ["Undetermined"] + config["ignore"]:
            if check in temp_samplename:
                add_sample = False
        if add_sample:
            samplelist.add(temp_samplename)
    for locsample in sorted(samplelist):
        #print(locsample.split('_')[0])
        sample_fastqs = sorted([filepath for filepath in fastqs if locsample in filepath])
        f = os.path.join(this_d, sample_fastqs[0])
        r = os.path.join(this_d, sample_fastqs[1])
        if locsample in config["samples"]:
            raise IOError("""Library {} was found in more than two directories.
          Please only include one set of readpairs per library:
           - {}
           - {}""".format(locsample, this_d, config["samples"][locsample]["f"]))
        else:
            config["samples"][locsample] = {}
            config["samples"][locsample]["f"] = f
            config["samples"][locsample]["r"] = r
            readkey = "reads/{}_R1_001.fastq.gz".format(locsample)
            if locsample.split('_')[0] in config["nextera"]:
                config["adapter_pairs"][readkey] = os.path.join(trimmomatic, "adapters/NexteraPE-PE.fa")
            else:
                config["adapter_pairs"][readkey] = os.path.join(trimmomatic, "adapters/TruSeq3-PE-2.fa")


#print(config["adapter_pairs"])
#for key in sorted(config["samples"].keys()):
#    print(key)
#    print("  - f: {}".format(config["samples"][key]["f"]))
#    print("  - r: {}".format(config["samples"][key]["r"]))

#print(config["adapter_pairs"])
print("The libraries identified for Nextera-specific trimming are:")
for key in sorted(config["samples"].keys()):
    if key.split('_')[0] in config["nextera"]:
        print("  - {}".format(key))
        #accession = "reads/{}_R1_001.fastq.gz".format(key)
        #print("  - {}".format(config["adapter_pairs"][accession]))


rule all:
    input:
        # compile fq-jt
        "bin/fqjt",
        "bin/minlen_pair",
        #get the reads
        expand("reads/{sample}_R2_001.fastq.gz", sample=config["samples"]),
        expand("reads/{sample}_R2_001.fastq.gz", sample=config["samples"]),
        #find the adapter sequences present in the files
        #expand("adapters/{sample}_adapters.fa", sample = config["samples"]),
        #expand("adapters/{sample}_merge.log", sample = config["samples"]),
        #fastqc of the raw files
        expand("fastqc/raw/{sample}_R1_001_fastqc.html", sample = config["samples"]),
        expand("fastqc/raw/{sample}_R1_001_fastqc.zip", sample = config["samples"]),
        expand("fastqc/raw/{sample}_R2_001_fastqc.html", sample = config["samples"]),
        expand("fastqc/raw/{sample}_R2_001_fastqc.zip", sample = config["samples"]),
        #trim the files using truseq or nextera adapters
        expand("trimmed/TruSeq3-PE-2/{sample}_R1_001.trim.fastq.gz", sample = config["samples"]),
        expand("trimmed/TruSeq3-PE-2/{sample}_R2_001.trim.fastq.gz", sample = config["samples"]),
        # now do the prox trimming
        #expand("trimmed/temp_final/{sample}_{dir}_001.final_temp.fastq.gz",
        #       sample = config["samples"], dir= ["R1", "R2"]),
        #expand("trimmed/final/{sample}_R1_001.fin.fastq.gz", sample = config["samples"]),
        expand("trimmed/final/{sample}_{dir}_001.fin.fastq.gz", sample = config["samples"], dir = ["R1", "R2"]),
        ##fastqc of the trimmed files
        expand("fastqc/trimmed/TruSeq3-PE-2/{sample}_R1_001.trim_fastqc.html", sample = config["samples"]),
        expand("fastqc/trimmed/TruSeq3-PE-2/{sample}_R1_001.trim_fastqc.zip",  sample = config["samples"]),
        expand("fastqc/trimmed/TruSeq3-PE-2/{sample}_R2_001.trim_fastqc.html", sample = config["samples"]),
        expand("fastqc/trimmed/TruSeq3-PE-2/{sample}_R2_001.trim_fastqc.zip",  sample = config["samples"]),
        # now make a report of the final files
        "report/final_report.txt"

rule compile_fqjt:
    """this executable is used to trim prox ligation data"""
    input:
        f1 = "bin/fq-jt.c"
    output:
        f2 = "bin/fqjt"
    shell:
        """
        gcc -gdwarf-2 -g {input.f1} -lz -o {output.f2}
        """

rule compile_minlen:
    """this executable is used to filter the pairs by minlen"""
    input:
        f1 = "bin/minlen_pair.c"
    output:
        f2 = "bin/minlen_pair"
    shell:
        """
        gcc {input.f1} -lz -o {output.f2}
        """


rule make_fake_reads:
    input:
        [config["samples"][sample][readdir] \
           for readdir in ['f', 'r'] \
           for sample in config["samples"]]
    output:
        forward = ancient(sorted(expand("reads/{sample}_R1_001.fastq.gz", sample=config["samples"]))),
        reverse = ancient(sorted(expand("reads/{sample}_R2_001.fastq.gz", sample=config["samples"])))
    message:
        "making soft links"
    run:
        dir_to_R = {"f": "R1", "r": "R2"}
        for this_sample in sorted(config["samples"]):
            for this_direction in ["f", "r"]:
                new_path = "reads/{}_{}_001.fastq.gz".format(this_sample,
                                                             dir_to_R[this_direction])
                if not os.path.exists(new_path):
                    os.symlink(config["samples"][this_sample][this_direction],
                           new_path)
                else:
                    print("{} already exists. skipping.".format(new_path))

rule find_adapters:
    input:
        f1 = "reads/{sample}_R1_001.fastq.gz",
        f2 = "reads/{sample}_R2_001.fastq.gz",
        bbmerge = bbmergepath
    output:
        adapters = "adapters/{sample}_adapters.fa",
        log = "adapters/{sample}_merge.log"
    threads:
        maxthreads
    shell:
        """{input.bbmerge} t={threads} in1={input.f1} in2={input.f2} \
        outa={output.adapters} > /dev/null 2> {output.log}; \
        sed -i -e 's/Read1_adapter/PrefixPE\/1/g' {output.adapters}; \
        sed -i -e 's/Read2_adapter/PrefixPE\/2/g' {output.adapters}"""

rule trim_pairs:
    input:
        f1 = "reads/{sample}_R1_001.fastq.gz",
        f2 = "reads/{sample}_R2_001.fastq.gz",
        trim_jar = os.path.join(trimmomatic, "trimmomatic-0.35.jar")
    output:
        f_paired =   "trimmed/TruSeq3-PE-2/{sample}_R1_001.trim.fastq.gz",
        r_paired =   "trimmed/TruSeq3-PE-2/{sample}_R2_001.trim.fastq.gz",
        f_unpaired = temp("trimmed/TruSeq3-PE-2/{sample}_1.trim.unpaired.fastq.gz"),
        r_unpaired = temp("trimmed/TruSeq3-PE-2/{sample}_2.trim.unpaired.fastq.gz")
    threads:
        15
    message:
        "adapter and quality trimming files for sample "
    run:
        #"""java -jar {input.trim_jar} PE \
        #-phred33 -threads {threads} \
        #{input.f1} {input.f2} \
        #{output.f_paired} \
        #{output.f_unpaired} \
        #{output.r_paired} \
        #{output.r_unpaired} \
        #ILLUMINACLIP:{input.adapter_file}:2:30:10 \
        #LEADING:3 TRAILING:3 \
        #SLIDINGWINDOW:4:15 MINLEN:36"""
        adapter_path = config["adapter_pairs"][input.f1]
        #print("lol here is adapter path")
        #print(adapter_path)
        shell("""java -jar {} PE \
              -phred33 -threads {} \
              {} {} \
              {} \
              {} \
              {} \
              {} \
              ILLUMINACLIP:{}:2:30:10:1:TRUE \
              LEADING:3 TRAILING:3 \
              SLIDINGWINDOW:4:15 MINLEN:36""".format(
                  input.trim_jar,
                  threads,
                  input.f1, input.f2,
                  output.f_paired,
                  output.f_unpaired,
                  output.r_paired,
                  output.r_unpaired,
                  adapter_path))

rule num_reads:
    input:
        forward_trim =   "trimmed/TruSeq3-PE-2/{sample}_R1_001.trim.fastq.gz",
        f1 = "reads/{sample}_R1_001.fastq.gz"
    output:
        numreads = "report/numreads/{sample}.numreads.txt",
        trimmed = "report/numreads/{sample}.trimmed_TruSeq3-PE-2.numreads.txt",

    shell:
        """echo $(bioawk -cfastx 'END{{print NR}}' {input.forward_trim}) > {output.trimmed}; \
        echo $(bioawk -cfastx 'END{{print NR}}' {input.f1}) > {output.numreads}
        """


rule library_efficiency:
    input:
        numreads = "report/numreads/{sample}.numreads.txt",
        trimmed = "report/numreads/{sample}.trimmed_TruSeq3-PE-2.numreads.txt",
    output:
        truseq = "report/efficiency/{sample}.TruSeq3-PE-2.efficiency.txt",
    message:
        "calculating the library efficiency"
    shell:
        """cat {input.trimmed} {input.numreads} | \
        paste - - | \
        awk '{{print $1/$2}}' > {output.truseq}; \
        """

rule fastqc_raw:
    input:
        f1 = "reads/{sample}_R1_001.fastq.gz",
        f2 = "reads/{sample}_R2_001.fastq.gz",
        adapter_path = os.path.join(trimmomatic, "adapters/TruSeq3-PE-2-fastqc.fa")
    output:
        "fastqc/raw/{sample}_R1_001_fastqc.html",
        "fastqc/raw/{sample}_R1_001_fastqc.zip",
        "fastqc/raw/{sample}_R2_001_fastqc.html",
        "fastqc/raw/{sample}_R2_001_fastqc.zip"
    message:
        "making a fastqc report of the raw reads"
    params:
        raw_dir = "fastqc/raw/"
    shell:
        """fastqc -a {input.adapter_path} \
        -o {params.raw_dir} \
        {input.f1} {input.f2}"""

rule fastqc_trimmed_TruSeqPE3:
    input:
        f1 =   "trimmed/TruSeq3-PE-2/{sample}_R1_001.trim.fastq.gz",
        f2 =   "trimmed/TruSeq3-PE-2/{sample}_R2_001.trim.fastq.gz",
        adapter_path = os.path.join(trimmomatic, "adapters/TruSeq3-PE-2-fastqc.fa")
    output:
        "fastqc/trimmed/TruSeq3-PE-2/{sample}_R1_001.trim_fastqc.html",
        "fastqc/trimmed/TruSeq3-PE-2/{sample}_R1_001.trim_fastqc.zip",
        "fastqc/trimmed/TruSeq3-PE-2/{sample}_R2_001.trim_fastqc.html",
        "fastqc/trimmed/TruSeq3-PE-2/{sample}_R2_001.trim_fastqc.zip"
    message:
        "making a fastqc report of the trimmed reads"
    params:
        raw_dir = "fastqc/trimmed/TruSeq3-PE-2/"
    shell:
        """fastqc -a {input.adapter_path} \
        -o {params.raw_dir} \
        {input.f1} {input.f2}"""

rule GATCGATC_f:
    input:
        f1 = "reads/{sample}_R1_001.fastq.gz",
        f2 = "reads/{sample}_R2_001.fastq.gz",
        f1t =   "trimmed/TruSeq3-PE-2/{sample}_R1_001.trim.fastq.gz",
        f2t =   "trimmed/TruSeq3-PE-2/{sample}_R2_001.trim.fastq.gz",

    output:
        f1o = "report/linkercontent/GATCGATC/{sample}_R1_001.GATCGATC.count",
        f2o = "report/linkercontent/GATCGATC/{sample}_R2_001.GATCGATC.count",
        f1to = "report/linkercontent/GATCGATC/{sample}_R1_001.trim_TruSeq3-PE-2.GATCGATC.count",
        f2to = "report/linkercontent/GATCGATC/{sample}_R2_001.trim_TruSeq3-PE-2.GATCGATC.count",

    shell:
        "set +o pipefail; "
        "bioawk -cfastx '{{print $seq}}' {input.f1} | grep 'GATCGATC'  - | wc -l > {output.f1o}; "
        "bioawk -cfastx '{{print $seq}}' {input.f2} | grep 'GATCGATC'  - | wc -l > {output.f2o}; " 
        "bioawk -cfastx '{{print $seq}}' {input.f1t} | grep 'GATCGATC' - | wc -l > {output.f1to}; "
        "bioawk -cfastx '{{print $seq}}' {input.f2t} | grep 'GATCGATC' - | wc -l > {output.f2to}; "

rule CATGCATG_num:
    input:
        f1 = "reads/{sample}_R1_001.fastq.gz",
        f2 = "reads/{sample}_R2_001.fastq.gz",
        f1t =   "trimmed/TruSeq3-PE-2/{sample}_R1_001.trim.fastq.gz",
        f2t =   "trimmed/TruSeq3-PE-2/{sample}_R2_001.trim.fastq.gz",

    output:
        f1o =       "report/linkercontent/CATGCATG/{sample}_R1_001.CATGCATG.count",
        f2o =       "report/linkercontent/CATGCATG/{sample}_R2_001.CATGCATG.count",
        f1to = "report/linkercontent/CATGCATG/{sample}_R1_001.trim_TruSeq3-PE-2.CATGCATG.count",
        f2to = "report/linkercontent/CATGCATG/{sample}_R2_001.trim_TruSeq3-PE-2.CATGCATG.count",

    shell:
        "set +o pipefail; "
        "bioawk -cfastx '{{print $seq}}' {input.f1} | grep 'CATGCATG'  - | wc -l > {output.f1o}; "
        "bioawk -cfastx '{{print $seq}}' {input.f2} | grep 'CATGCATG'  - | wc -l > {output.f2o}; "
        "bioawk -cfastx '{{print $seq}}' {input.f1t} | grep 'CATGCATG' - | wc -l > {output.f1to}; "
        "bioawk -cfastx '{{print $seq}}' {input.f2t} | grep 'CATGCATG' - | wc -l > {output.f2to};"


rule AATTAATT_num:
    input:
        f1 = "reads/{sample}_R1_001.fastq.gz",
        f2 = "reads/{sample}_R2_001.fastq.gz",
        f1t =   "trimmed/TruSeq3-PE-2/{sample}_R1_001.trim.fastq.gz",
        f2t =   "trimmed/TruSeq3-PE-2/{sample}_R2_001.trim.fastq.gz",

    output:
        f1o =       "report/linkercontent/AATTAATT/{sample}_R1_001.AATTAATT.count",
        f2o =       "report/linkercontent/AATTAATT/{sample}_R2_001.AATTAATT.count",
        f1to = "report/linkercontent/AATTAATT/{sample}_R1_001.trim_TruSeq3-PE-2.AATTAATT.count",
        f2to = "report/linkercontent/AATTAATT/{sample}_R2_001.trim_TruSeq3-PE-2.AATTAATT.count",

    shell:
        "set +o pipefail; "
        "bioawk -cfastx '{{print $seq}}' {input.f1} | grep 'AATTAATT'  - | wc -l > {output.f1o}; "
        "bioawk -cfastx '{{print $seq}}' {input.f2} | grep 'AATTAATT'  - | wc -l > {output.f2o}; "
        "bioawk -cfastx '{{print $seq}}' {input.f1t} | grep 'AATTAATT' - | wc -l > {output.f1to}; "
        "bioawk -cfastx '{{print $seq}}' {input.f2t} | grep 'AATTAATT' - | wc -l > {output.f2to};"


def read_number_from_file(filename):
    with open(filename, "r") as f:
       for line in f:
           if line.strip():
               return line.strip()

rule collate_report:
    """ this method prints out the final qc report of all the libraries."""
    input:
        raw_num =       expand("report/numreads/{sample}.numreads.txt", sample=config["samples"]),
        trim_num =      expand("report/numreads/{sample}.trimmed_TruSeq3-PE-2.numreads.txt", sample=config["samples"]),
        efficiency =    expand("report/efficiency/{sample}.TruSeq3-PE-2.efficiency.txt", sample=config["samples"]),
        f1dpn =         expand("report/linkercontent/GATCGATC/{sample}_R1_001.GATCGATC.count", sample=config["samples"]),
        f2dpn =         expand("report/linkercontent/GATCGATC/{sample}_R2_001.GATCGATC.count", sample=config["samples"]),
        f1mluc =        expand("report/linkercontent/AATTAATT/{sample}_R1_001.AATTAATT.count", sample=config["samples"]),
        f2mluc =        expand("report/linkercontent/AATTAATT/{sample}_R2_001.AATTAATT.count", sample=config["samples"]),
        f1fat =         expand("report/linkercontent/CATGCATG/{sample}_R1_001.CATGCATG.count", sample=config["samples"]),
        f2fat =         expand("report/linkercontent/CATGCATG/{sample}_R2_001.CATGCATG.count", sample=config["samples"]),
        f1tdpn =        expand("report/linkercontent/GATCGATC/{sample}_R1_001.trim_TruSeq3-PE-2.GATCGATC.count", sample=config["samples"]),
        f2tdpn =        expand("report/linkercontent/GATCGATC/{sample}_R2_001.trim_TruSeq3-PE-2.GATCGATC.count", sample=config["samples"]),
        f1tmluc =       expand("report/linkercontent/AATTAATT/{sample}_R1_001.trim_TruSeq3-PE-2.AATTAATT.count", sample=config["samples"]),
        f2tmluc =       expand("report/linkercontent/AATTAATT/{sample}_R2_001.trim_TruSeq3-PE-2.AATTAATT.count", sample=config["samples"]),
        f1tmfat =       expand("report/linkercontent/CATGCATG/{sample}_R1_001.trim_TruSeq3-PE-2.CATGCATG.count", sample=config["samples"]),
        f2tmfat =       expand("report/linkercontent/CATGCATG/{sample}_R2_001.trim_TruSeq3-PE-2.CATGCATG.count", sample=config["samples"]),

    output:
        "report/final_report.txt"
    run:
        finalout_handle = open(output[0], "w")
        print("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}".format(
              "libname",
              "num_reads",
              "num_reads_trimmed_TruSeq3-PE-2",
              "trim_efficiency_TruSeq3-PE-2",
              "f_GATCGATC",
              "r_GATCGATC",
              "f_AATTAATT",
              "r_AATTAATT",
              "f_CATGCATG",
              "r_CATGCATG",
              "f_trim_TruSeq3-PE-2_GATCGATC",
              "r_trim_TruSeq3-PE-2_GATCGATC",
              "f_trim_TruSeq3-PE-2_AATTAATT",
              "r_trim_TruSeq3-PE-2_AATTAATT",
              "f_trim_TruSeq3-PE-2_CATGCATG",
              "r_trim_TruSeq3-PE-2_CATGCATG"),
              file = finalout_handle)
        for thiss in sorted(config["samples"]):
            num_reads = read_number_from_file("report/numreads/{}.numreads.txt".format(thiss))
            nreads_tr = read_number_from_file("report/numreads/{}.trimmed_TruSeq3-PE-2.numreads.txt".format(thiss))
            trim_effi = read_number_from_file("report/efficiency/{}.TruSeq3-PE-2.efficiency.txt".format(thiss))
            f_GATCGAT = read_number_from_file("report/linkercontent/GATCGATC/{}_R1_001.GATCGATC.count".format(thiss))
            r_GATCGAT = read_number_from_file("report/linkercontent/GATCGATC/{}_R2_001.GATCGATC.count".format(thiss))
            f_AATTAAT = read_number_from_file("report/linkercontent/AATTAATT/{}_R2_001.AATTAATT.count".format(thiss))
            r_AATTAAT = read_number_from_file("report/linkercontent/AATTAATT/{}_R2_001.AATTAATT.count".format(thiss))
            f_CATGCAT = read_number_from_file("report/linkercontent/CATGCATG/{}_R1_001.CATGCATG.count".format(thiss))
            r_CATGCAT = read_number_from_file("report/linkercontent/CATGCATG/{}_R2_001.CATGCATG.count".format(thiss))
            f_trim_GA = read_number_from_file("report/linkercontent/GATCGATC/{}_R1_001.trim_TruSeq3-PE-2.GATCGATC.count".format(thiss))
            r_trim_GA = read_number_from_file("report/linkercontent/GATCGATC/{}_R2_001.trim_TruSeq3-PE-2.GATCGATC.count".format(thiss))
            f_trim_AA = read_number_from_file("report/linkercontent/AATTAATT/{}_R1_001.trim_TruSeq3-PE-2.AATTAATT.count".format(thiss))
            r_trim_AA = read_number_from_file("report/linkercontent/AATTAATT/{}_R2_001.trim_TruSeq3-PE-2.AATTAATT.count".format(thiss))
            f_trim_CA = read_number_from_file("report/linkercontent/CATGCATG/{}_R1_001.trim_TruSeq3-PE-2.CATGCATG.count".format(thiss))
            r_trim_CA = read_number_from_file("report/linkercontent/CATGCATG/{}_R2_001.trim_TruSeq3-PE-2.CATGCATG.count".format(thiss))


            print("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}".format(
                    thiss,
                    num_reads,
                    nreads_tr,
                    trim_effi,
                    int(f_GATCGAT)/int(num_reads),
                    int(r_GATCGAT)/int(num_reads),
                    int(f_AATTAAT)/int(num_reads),
                    int(r_AATTAAT)/int(num_reads),
                    int(f_CATGCAT)/int(num_reads),
                    int(r_CATGCAT)/int(num_reads),
                    int(f_trim_GA)/int(nreads_tr),
                    int(r_trim_GA)/int(nreads_tr),
                    int(f_trim_AA)/int(nreads_tr),
                    int(r_trim_AA)/int(nreads_tr),
                    int(f_trim_CA)/int(nreads_tr),
                    int(r_trim_CA)/int(nreads_tr)),
                  file = finalout_handle)
        finalout_handle.close()

rule trim_prox:
    input:
        f_file =   "trimmed/TruSeq3-PE-2/{sample}_{dir}_001.trim.fastq.gz",
        fqjt = "bin/fqjt"
    output:
        final = temp("trimmed/temp_final/{sample}_{dir}_001.final_temp.fastq.gz")
    threads:
        1
    run:
        if not os.path.isdir("trimmed/temp_final"):
            os.mkdir("trimmed/temp_final")
        seq_lookup = wildcards.sample.split('_')[0]
        prox = ""
        if seq_lookup in config["prox"]["AATTAATT"]:
            prox = "AATTAATT"
        elif seq_lookup in config["prox"]["GATCGATC"]:
            prox = "GATCGATC"
        elif seq_lookup in config["prox"]["CATGCATG"]:
            prox = "CATGCATG"

        assert (len(prox) % 2) == 0

        if prox == "":
            os.symlink(os.path.realpath(input.f_file), output.final)
        else:
            callstring = "./{} -f {} -t {} -l {} | gzip > {}".format(
                input.fqjt, input.f_file, prox, int(len(prox)/2), output.final)
            #print(callstring)
            call(callstring, shell=True)

rule trimmomatic_prox_remove_short_seqs:
    input:
        f_file =   "trimmed/temp_final/{sample}_R1_001.final_temp.fastq.gz",
        r_file =   "trimmed/temp_final/{sample}_R2_001.final_temp.fastq.gz",
        minlen_pair = "bin/minlen_pair"
    output:
        f_paired = "trimmed/final/{sample}_R1_001.fin.fastq.gz",
        r_paired = "trimmed/final/{sample}_R2_001.fin.fastq.gz",
    threads:
        1
    params:
        minlen = 36
    run:
        if not os.path.isdir("trimmed/final"):
            os.mkdir("trimmed/final")
        seq_lookup = wildcards.sample.split('_')[0]
        prox = ""
        if seq_lookup in config["prox"]["AATTAATT"]:
            prox = "AATTAATT"
        elif seq_lookup in config["prox"]["GATCGATC"]:
            prox = "GATCGATC"
        elif seq_lookup in config["prox"]["CATGCATG"]:
            prox = "CATGCATG"

        assert (len(prox) % 2) == 0
        if prox == "":
            f_orig = "trimmed/TruSeq3-PE-2/{}_R1_001.trim.fastq.gz".format(wildcards.sample)
            r_orig = "trimmed/TruSeq3-PE-2/{}_R2_001.trim.fastq.gz".format(wildcards.sample)
            os.symlink(os.path.realpath(f_orig), output.f_paired)
            os.symlink(os.path.realpath(r_orig), output.r_paired)
        else:
            callstring = """./{} {} {} {} {} {}""".format(
                input.minlen_pair, params.minlen,
                input.f_file, input.r_file,
                output.f_paired, output.r_paired)
            call(callstring, shell=True)
