from Bio import Entrez
import xml.etree.ElementTree as ET
import sys
from pysradb.sraweb import SRAweb
import os
import re
import nltk
from bs4 import BeautifulSoup
from nltk.tokenize import RegexpTokenizer
import multiprocessing
import getopt, sys
from Bio import SeqIO
import re
import glob


# esearch -db assembly -query "hg19" | esummary
# wget

# This function have been moved up

home_path = "/home/saad18409"
fasterq_path = "/home/saad18409/sratoolkit.2.10.5-ubuntu64/bin"
hisat2_path = "/home/saad18409/hisat2-2.1.0"
samtools_path = "/home/saad18409/bin"
counts_path = "/home/saad18409/subread-2.0.1-Linux-x86_64/bin"
fastqc_path = "/home/saad/FastQC"


def idselector(tree):
    """ Returns top 1 id after searching"""
    x = 0
    req_id = ""
    for node in tree.iter("IdList"):
        for elem in node.iter():
            if x == 1:
                print(elem.tag, elem.text)  # 1st Result found
                id = elem.text
                break
            x += 1
    return id


def check_with_assembly(accession):
    """
        For all the potential genome builds, check one by one which ine actually exists on ncbi assembly
        accession = potential genome Build
    """
    handle = Entrez.esearch(db="assembly", term=accession)
    xml_result = handle.read()
    print(xml_result)
    tree = ET.fromstring(xml_result)
    for item in tree.findall("Count"):
        if int(item.text) > 0:
            return True
        else:
            return False


def findGenomeBuild(htmlFile):
    """
        For Given HTML File this function scrapes that HTML file
        for genome_build and if it exists returns genome_build
        else None
    """
    tk = RegexpTokenizer(
        "((?:[a-zA-Z_]+[0-9_]|[0-9_]+[a-zA-Z_])[a-zA-Z0-9_]*)", gaps=False
    )
    selected_text_for_build = ""
    with open("{}/experimentSample.html".format(home_path)) as htmlFile:
        iter_count = 0
        target_para = 0
        prev = ""
        for line in htmlFile:
            if "data" in line.lower() and "processing" in line.lower():
                if "justify" in line.lower():
                    selected_text_for_build = line
                    break
                target_para = iter_count
            if iter_count == target_para + 1:
                selected_text_for_build = line
            iter_count += 1
        htmlFile.close()
    print(selected_text_for_build)
    selected_text_for_build = selected_text_for_build.split("<br>")
    no_instances_found = False
    index_of_build = len(selected_text_for_build) + 2
    for i, string in enumerate(selected_text_for_build):
        if re.findall("build", string):
            if re.findall("genome", string):
                no_instances_found = True
                index_of_build = i
                break
            else:
                no_instances_found = True
                index_of_build = i
    final_candidates_build = []
    if no_instances_found == True:
        soup = BeautifulSoup(selected_text_for_build[index_of_build], "html.parser")
        stripped_text = soup.get_text()
        geek = tk.tokenize(stripped_text)
        if len(geek) > 0:
            print(geek)
            for j in geek:
                if bool(re.search(r"^([^0-9]*)$", j)) == False:
                    if check_with_assembly(j) == True:
                        final_candidates_build.append(j)
    else:
        for i in selected_text_for_build:
            soup = BeautifulSoup(i, "html.parser")
            stripped_text = soup.get_text()
            geek = tk.tokenize(stripped_text)
            if len(geek) > 0:
                for j in geek:
                    if bool(re.search(r"^([^0-9]*)$", j)) == False:
                        if check_with_assembly(j) == True:
                            final_candidates_build.append(j)
    print(final_candidates_build)
    if len(final_candidates_build) > 0:
        return final_candidates_build[0]
    else:
        return None


def writeMessage(message):
    pythonoutputfile = open("python_script_realtime_log.txt", "a")
    pythonoutputfile.write(message + "\n")
    pythonoutputfile.close()
    return 0


"""
    given GEOID (GSE)
    This Function esearch it's corresponding GSM ID
    and then downloads it's html file and returns that
"""


def getHTML(geo_id):
    paramEutils = {"usehistory": "Y"}
    handle = Entrez.esearch(db="gds", term=geo_id, **paramEutils)
    res = Entrez.read(handle)
    paramEutils["WebEnv"] = res["WebEnv"]
    paramEutils["query_key"] = res["QueryKey"]
    paramEutils["rettype"] = "xml"
    result = Entrez.esummary(db="gds", **paramEutils)
    xml_res = Entrez.read(result)
    SampleID = None
    for id in xml_res:
        SampleID = id["Accession"]
        if SampleID[0:3] == "GSM":
            break
    command = "wget https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={} -O {}/experimentSample.html".format(
        SampleID, home_path
    )
    os.system(command)
    htmlFile = open("{}/experimentSample.html".format(home_path), "r")
    return htmlFile

    """
    This Function takes GSM_ID as input paramater
    and returns whether the result is strand-specific 
    or not 
    """


def isStranded(GSM_ID):
    command = "wget https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={} -O strandedness.html".format(
        GSM_ID
    )
    os.system(command)
    strandFile = open("strandedness.html", "r")
    for line in strandFile:
        if "strand-specific" in line.lower():
            strandFile.close()
            return True
    strandFile.close()
    return False


"""
    Still Experimental (Under Build)
    For Given HTML File this function scrapes that HTML file
    for genome_build and if it exists returns genome_build
    else None
"""


def get_ens_annotations(genomeBuild):
    """
        downloads ensembl annotations from goldenpath
        In-case ensembl annotations exist we will check for pre built indexes
        downloaded annotation file name = hisat2_annotation.gtf.gz
    """
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/{}/bigZips/genes/{}.ensGene.gtf.gz".format(
        genomeBuild.lower(), genomeBuild.lower()
    )
    ret_value = os.system("wget {} -O hisat2_annotation.gtf.gz".format(url))
    if ret_value != 0:
        os.system("rm hisat2_annotation.gtf.gz")
        return False, "hisat2_annotation.gtf.gz"
    else:
        return True, "hisat2_annotation.gtf.gz"


def check_pre_built(genomeBuild, ht2_idx_url, org_name, to_build=True):
    """
        Main function which calls all the other download and build functions
        first check if genome build exists
        if it does the check for pre built indexes
        in all other cases download refseq genome and refseq annotations and then build indexes
    """
    annotationFileName = ""
    refGenomeFileName = ""
    print("in check {}".format(genomeBuild))
    if genomeBuild != None:
        ret, annotationFileName, refGenomeFileName = download_pre_built_idx(
            ht2_idx_url, genomeBuild
        )
        if ret == -1:
            print("building Index")
            (annotationFileName, refGenomeFileName) = downloadRefGenome(
                genomeBuild, org_name
            )
            print(annotationFileName, refGenomeFileName)
            if to_build == True:
                build_status = build_index(refGenomeFileName)
    else:
        print("building Index else")
        annotationFileName, refGenomeFileName = downloadRefGenome(genomeBuild, org_name)
        if to_build == True:
            build_status = build_index(refGenomeFileName)
    return annotationFileName, refGenomeFileName


def download_pre_built_idx(ht2_idx_url, genomeBuild):
    """
        First call annotation function to check if ensembl annotation exists
        If it does then check if pre built indexes exist
        If they do then download and untar them and save with genome build name
        In all other cases return -1
    """
    annotation_exists, annotationFileName = get_ens_annotations(genomeBuild)
    print(annotation_exists)
    if annotation_exists == True:
        return_val = os.system(
            "wget {}/{}.tar.gz".format(ht2_idx_url, genomeBuild.lower())
        )
        if return_val == 0:
            os.system("tar -xf {}.tar.gz".format(genomeBuild))
            os.chdir(genomeBuild)
            os.system("mv * ../")
            to_ret_path = os.getcwd()
            os.chdir("..")
            for i in range(1,9):
                os.system("mv genome.{}.ht2 ht2idxes.{}.ht2".format(i, i))
            return to_ret_path, annotationFileName, genomeBuild.lower()
        else:
            return -1, annotationFileName, genomeBuild.lower()
    else:
        return -1, annotationFileName, genomeBuild.lower()


"""
    This Function Intakes orgainsm name and then returns
    refseq_id for it's latest genome build 
"""


def getLatestRefGenomeName(org_name):
    """
        This Function Intakes orgainsm name and then returns
        refseq_id for it's latest genome build 
    """
    print("ORG ANAME {}".format(org_name))
    command = "{}/datasets assembly-descriptors tax-name '{}' --refseq > orgRefSeqInfo.txt".format(
        home_path, org_name
    )
    os.system(command)
    refSeqAccessionID = None
    file = open("orgRefSeqInfo.txt", "r")
    textOfAccessionFile = []
    for line in file:
        textOfAccessionFile.extend(line.split('"'))
    for words in range(len(textOfAccessionFile)):
        if "assembly_accession" in textOfAccessionFile[words]:
            refSeqAccessionID = textOfAccessionFile[words + 2]
            break
    return refSeqAccessionID


"""
    For Given Refseq id this function downloads reference Genome
    and it's annotations and returns the name of the downloaded
    return annotation File name, reference Genome Filename 
"""


def downloadGenomeUsingWget(refSeqId):
    """
        For Given Refseq id this function downloads reference Genome
        and it's annotations and returns the name of the downloaded
        return annotation File name, reference Genome Filename 
    """
    print("\n REFSEQID {} \n".format(refSeqId))
    paramEutils = {"usehistory": "Y"}
    handle = Entrez.esearch(db="assembly", term=refSeqId, **paramEutils)
    res = Entrez.read(handle)
    paramEutils["WebEnv"] = res["WebEnv"]
    paramEutils["query_key"] = res["QueryKey"]
    paramEutils["rettype"] = "xml"
    result = Entrez.esummary(db="assembly", **paramEutils)
    xml_res = result.read()
    ftpPathToGenome = None
    root = ET.fromstring(xml_res)

    for i in root.iter("FtpPath_RefSeq"):
        val = i.text
        if val != None:
            ftpPathToGenome = val
    name = (ftpPathToGenome).split("/")
    print(str(name[-1]))
    command = "wget -crv -t 45 '{}/{}_genomic.gtf.gz' -O {}_annotations.gtf.gz".format(
        ftpPathToGenome, name[-1], name[-1]
    )
    os.system(command)
    annotationFileName = "{}_annotations.gtf.gz".format(name[-1])
    command = "wget -crv -t 45 '{}/{}_genomic.fna.gz' -O {}_genome.fna.gz".format(
        ftpPathToGenome, name[-1], name[-1]
    )
    os.system(command)
    refgenomeFileName = "{}_genome.fna.gz".format(name[-1])
    os.system("gunzip {}".format(refgenomeFileName))
    refgenomeFileName = str(name[-1]) + "_genome.fna"
    return annotationFileName, refgenomeFileName


"""
    Put GEOID->eSearch->  get GSM_ID from Here
    Put GSM_ID -> wget -> get HTML Page
    Scrape HTML_PAGE -> Genome_Build
    If Genome_Build not Exist
    do 
        ./dataset with organism name -> Refseq Accession ID  (GCF_XXXX )
        use wget to download Reference genome using Accession ID
    else if Genome Build Exists
    do
        esearch -db assembly --query build_name -> Refseq Accession ID 
        use wget to download Reference genome using Accession ID
"""


def downloadRefGenome(genomeBuild, org_name):
    """
        Put .geo_id->eSearch->  get GSM_ID from Here
        Put GSM_ID -> wget -> get HTML Page
        Scrape HTML_PAGE -> Genome_Build
        If Genome_Build not Exist
        do 
            ./dataset with organism name -> Refseq Accession ID  (GCF_XXXX )
            use wget to download Reference genome using Accession ID
        else if Genome Build Exists
        do
            esearch -db assembly --query build_name -> Refseq Accession ID 
            use wget to download Reference genome using Accession ID
    """
    print(genomeBuild)
    annotationFileName, refGenomeFileName = None, None
    if genomeBuild == None:
        refGenomeName = getLatestRefGenomeName(org_name)
        annotationFileName, refGenomeFileName = downloadGenomeUsingWget(refGenomeName)
        print("annotatin done")
        return annotationFileName, refGenomeFileName
    else:
        print(genomeBuild, "<- Genome Build")
        # Search for Refseq id in esearch assmbly Results

        print("EXECUTING COMMANFG LONG COMMAND")
        paramEutils = {"usehistory": "Y"}
        handle = Entrez.esearch(db="assembly", term=genomeBuild, **paramEutils)
        res = Entrez.read(handle)
        paramEutils["WebEnv"] = res["WebEnv"]
        paramEutils["query_key"] = res["QueryKey"]
        paramEutils["rettype"] = "xml"
        result = Entrez.esummary(db="assembly", **paramEutils)
        xml_res = result.read()
        # erret = os.system(command)
        refGenomeName = None
        root = ET.fromstring(xml_res)
        refGenomeName = None
        for i in root.iter("FtpPath_RefSeq"):
            val = i.text
            print(val)
            if val != None:
                refGenomeName = val
        # If RefSeq id doesn't Exists in esearch, then download latest genome
        if refGenomeName == None:
            refGenomeName = getLatestRefGenomeName(org_name)
            print(refGenomeName)
            (annotationFileName, refGenomeFileName,) = downloadGenomeUsingWget(
                refGenomeName
            )
            return annotationFileName, refGenomeFileName
        else:  # If Refseq Id Exists in Esearch scrape then download that build
            print("IN ELSE BLOCK")
            name = (refGenomeName).split("/")
            name_minus_one = name[-1]
            annotationFileName = str(name_minus_one) + "_annotations.gtf.gz"
            refGenomeFileName = str(name_minus_one) + "_genome.fna.gz"

            print(refGenomeFileName)
            command = "wget -crv -t 45 '{}/{}_genomic.gtf.gz' -O {}".format(
                refGenomeName, name[-1], annotationFileName
            )
            os.system(command)
            command = "wget -crv -t 45 '{}/{}_genomic.fna.gz' -O {}".format(
                refGenomeName, name[-1], refGenomeFileName
            )
            os.system(command)
            os.system("gunzip {}".format(refGenomeFileName))
            refGenomeFileName = str(name_minus_one) + "_genome.fna"
            print(name[-1])
            print("IN DOWNLOAD GENOME FUNCTION", refGenomeFileName)
            return (annotationFileName, refGenomeFileName)


def copyFilesFromSubdirs(sraFileNames):
    listOfDirs = os.listdir()
    for i in range(len(listOfDirs)):
        os.system("mv {}/* {}".format(listOfDirs[i], os.getcwd()))
    for dirs in listOfDirs:
        os.system("rm -rf {}".format(dirs))
    sraFileNames = os.listdir()
    return sraFileNames


def download_fq_file(db, sra_id):
    print(os.getcwd())
    # os.system('mkdir {}'.format(sra_id))
    metadata = db.sra_metadata(sra_id, detailed=True)
    print(metadata)
    # os.chdir(sra_id)
    for run_acc in metadata.loc[:, "run_accession"]:
        print(run_acc)
        return_value = os.system(
            "{}/fasterq-dump {} -p -t /home/saad18409/temp_files".format(fasterq_path, str(run_acc))
        )
        print(return_value)


def createFastqFiles(allSRRIds):
    oneNames = []
    secondNames = []
    print(os.getcwd())
    check_paired_list = glob.glob("*_1.fastq")
    check_paired = False
    if len(check_paired_list) > 0:
        check_paired = True
    for i in range(len(allSRRIds)):
        if check_paired == True:
            oneNames.append(allSRRIds[i] + "_1.fastq")
            secondNames.append(allSRRIds[i] + "_2.fastq")
        else:
            oneNames.append(allSRRIds[i] + ".fastq")
    return oneNames, secondNames


def preprocess(
    firstList,
    secondList,
    corr_rRNA,
    corr_trim,
    sortmernaDbDir,
    min_qual1,
    min_read_len1,
    len_window1,
):
    """
    corr_rRNA trims out rRNA reads
    corr_trim trims out bases with QS less than threshold
    """
    if corr_rRNA == True:
        firstList, secondList = removerRnaContamination(
            firstList, secondList, sortmernaDbDir, [], []
        )
    if corr_trim == True:
        firstList, secondList = trimmingBadReads(
            firstName=firstList,
            secondName=secondList,
            min_qual=min_qual1,
            min_read_len=min_read_len1,
            len_window=len_window1,
        )
    return firstList, secondList


def removerRnaContamination(
    firstName: list,
    secondName: list,
    sortmernaDbDir,
    outputFirstName: list,
    outputSecondName: list,
):
    """
    Uses sortmerna library to filter out rRNA contamination
    """
    os.system("cp {}/* {}".format(sortmernaDbDir, os.getcwd()))
    referenceStr = "--ref "  # Complete this string
    if secondName == []:
        for i, read in enumerate(firstName):
            os.system(
                "sortmerna {} --reads {} --fastx --other {}/rRNAfiltered{}".format(
                    referenceStr, read, os.getcwd(), i
                )
            )
            outputFirstName.append("rRNAfiltered{}.fastq".format(i))
    else:
        for i in range(len(firstName)):
            os.system(
                "sortmerna {} --reads {} --reads {} --fastx --other {}/rRNAfiltered{} --out2".format(
                    referenceStr, firstName[i], secondName[i], os.getcwd(), i
                )
            )
            outputFirstName.append("rRNAfiltered{}_fwd.fastq".format(i))
            outputSecondName.append("rRNAfiltered{}_rev.fastq".format(i))
    os.system("rm -f {}").format(referenceStr)  # Complete this
    return (outputFirstName, outputSecondName)


def downloadOnlySRA(db, sra_id):
    df = db.sra_metadata(sra_id, detailed=True)
    print(os.getcwd())
    # os.system('mkdir {}'.format(sra_id))
    metadata = df
    print(metadata)
    # os.chdir(sra_id)
    for run_acc in metadata.loc[:, "run_accession"]:
        print(run_acc)
        return_value = os.system("{}/prefetch -p -O . {}".format(fasterq_path, str(run_acc)))
        print(return_value)


def trimmingBadReads(
    phred=33, firstName=[], secondName=[], min_qual=20, min_read_len=30, len_window=5
):
    """
    uses trimmomatic to trim out bases with low scores,
    uses default parameters i.e. removes bases with quality score < 20 in windows of 5
    """
    trim_home = os.environ["TRIMHOME"]
    illuminaclip_adapters = "ILLUMINACLIP:{}/adapters/TruSeq3-SE.fa:2:30:10".format(
        trim_home
    )
    illuminaclip_Attribute = "SLIDINGWINDOW:{}:{} MINLEN:{}".format(
        len_window, min_qual, min_read_len
    )
    attribute = ""
    if secondName == []:
        attribute = "SE -threads {} -phred{}".format(12, phred)
    else:
        attribute = "PE -threads {} -phred{}".format(12, phred)
    if secondName == []:
        firstNameTrimmed = []
        secondNameTrimmed = []
        for i in range(len(firstName)):
            os.system(
                "java -jar {}/trimmomatic-0.39.jar {} {}/{} {}/trimmedRead{}.fastq {} {}".format(
                    trim_home,
                    attribute,
                    os.getcwd(),
                    firstName[i],
                    os.getcwd(),
                    firstName[i],
                    illuminaclip_adapters,
                    illuminaclip_Attribute,
                )
            )
            firstNameTrimmed.append("trimmedRead{}".format(firstName[i]))
        return firstNameTrimmed, secondNameTrimmed
    else:
        firstNameTrimmed = []
        secondNameTrimmed = []
        for i in range(len(firstName)):
            fwd_input = os.getcwd() + "/" + firstName[i]
            rev_input = os.getcwd() + "/" + secondName[i]
            output_forward_paired = os.getcwd() + "/trimmedRead" + firstName[i]
            output_forward_unpaired = os.getcwd() + "/" + "removethis_1"
            output_reverse_paired = os.getcwd() + "/trimmedRead" + secondName[i]
            output_reverse_unpaired = os.getcwd() + "/" + "removethis_2"
            os.system(
                "java -jar {}/trimmomatic-0.39.jar {} {} {} {}.fastq {}.fastq {}.fastq {}.fastq {} {}".format(
                    trim_home,
                    attribute,
                    fwd_input,
                    rev_input,
                    output_forward_paired,
                    output_forward_unpaired,
                    output_reverse_paired,
                    output_reverse_unpaired,
                    illuminaclip_adapters,
                    illuminaclip_Attribute,
                )
            )
            os.system("rm -f *removethis_1.fastq *removethis_2.fastq")
            firstNameTrimmed.append("trimmedRead{}".format(firstName[i]))
            secondNameTrimmed.append("trimmedRead{}".format(secondName[i]))
        return firstNameTrimmed, secondNameTrimmed


def build_index(refGenome):
    """
        Builds indexes using hisat2-build executable
        refGenome = takes the path to reference genome file
    """
    command = "{}/hisat2-build -p {} {} {}".format(hisat2_path, 12, refGenome, "ht2idxes")
    ret_val = os.system(command)
    return ret_val


def startAlignment(
    indexLocation: str, firstFileNamesList: list, secondFilesNameList=[]
):
    outputFilenames = []
    os.system("touch alignmentSummary")
    if len(secondFilesNameList) == 0:
        for i in range(len(firstFileNamesList)):
            command = "{}/hisat2 --summary-file alignmentSummary -x {} -p {} -U {} -S output_{}.sam".format(
                hisat2_path, indexLocation, 12, firstFileNamesList[i], firstFileNamesList[i][:-6],
            )
            os.system(command)
            outputName = "output_{}.sam".format(firstFileNamesList[i][:-6])
            outputFilenames.append(outputFilenames)
    else:
        for i in range(len(firstFileNamesList)):
            command = "{}/hisat2 --summary-file alignmentSummary -x {} -p {} -1 {} -2 {} -S output_{}.sam".format(
                hisat2_path, indexLocation, 12, firstFileNamesList[i], secondFilesNameList[i], firstFileNamesList[i][:-6],
            )
            os.system(command)
            outputName = "output_{}.sam".format(firstFileNamesList[i][:-6])
            outputFilenames.append(outputFilenames)
    return outputFilenames


def convertSamToBam(samFilenames):
    bamFileNames = []
    for i in range(len(samFilenames)):
        command = "{}/samtools view -S -b {} > {}".format(
            samtools_path, samFilenames[i], samFilenames[i][0:-3] + "bam"
        )
        os.system(command)
        bamFileNames.append(samFilenames[i][0:-3] + "bam")
    return bamFileNames


def sortBamFiles(BamfileNames: list):
    sortedBamFileNames = []
    for i in range(len(BamfileNames)):
        command = "{}/samtools sort {} -o {}".format(
            samtools_path ,BamfileNames[i], BamfileNames[i][0:-3] + "sorted.bam"
        )
        os.system(command)
        sortedBamFileNames.append(BamfileNames[i][0:-3] + "sorted.bam")
    return sortedBamFileNames


def getAllSRRIds(db, sra_id):
    metadata = db.sra_metadata(sra_id, detailed=True)
    print(metadata)
    # os.chdir(sra_id)
    to_ret = []
    for run_acc in metadata.loc[:, "run_accession"]:
        print(run_acc)
        to_ret.append(str(run_acc))
    return to_ret


def qualityControl():
    cwd = os.getcwd()
    os.system("mkdir qcReports")
    os.system("{}/fastqc -t {} -o qcReports/ *.fastq".format(fastqc_path,12))
    os.chdir("qcReports")
    os.system("multiqc .")
    os.chdir("..")
    # send(multiqc_report.html)     IMPLEMENT THIS


def createCountMatrix(pathToGTF, bamInputFile):
    """
    Exon specific expression and alternative splicing not implemented
    Uses featureCounts

    Reference - http://genomespot.blogspot.com/2015/01/generate-rna-seq-count-matrix-with.html
    """
    bam_list = ""
    for i in bamInputFile:
        bam_list += i
        bam_list += " "

    os.system(
        "{}/featureCounts -Q 10 -T {} -a {} -o countMatrix {}".format(
            counts_path, 12, pathToGTF, bam_list
        )
    )


def initial_step(geo_id):
    Entrez.email = "bhavay18384@iiitd.ac.in"
    writeMessage("Process ID: {}".format(os.getpid()))
    refgenome = ""
    search_id = geo_id
    print("Input GEOID is", search_id)
    writeMessage("Input GEOID is {}".format(search_id))
    handle = Entrez.esearch(db="gds", term=search_id)
    pp = handle.read()
    tree = ET.fromstring(pp)
    writeMessage("Selecting the following ID")
    id = idselector(tree)
    writeMessage(str(id))
    handle = Entrez.esummary(db="gds", id=id)
    pp = handle.read()
    tree = ET.fromstring(pp)
    targetsra = ""
    for item in tree.findall("DocSum"):
        for item2 in item.findall("Item"):
            if item2.attrib["Name"] == "taxon":
                refgenome = item2.text
            if item2.attrib["Name"] == "ExtRelations":
                for item3 in item2.findall("Item"):
                    if item3.attrib["Name"] == "ExtRelation":
                        for item4 in item3.findall("Item"):
                            if item4.attrib["Name"] == "TargetObject":
                                targetsra = item4.text
    print("SRA in relation is ", targetsra)
    writeMessage("SRA in relation is {}".format(targetsra))
    print("Selecting the following ID")
    writeMessage("Selecting the following ID")
    writeMessage(str(id))
    print("Organism is", refgenome)
    # os.system("find . -name " + targetsra + "> downloadPath.txt")
    # downloadFileLocation = open("downloadPath.txt", "r").readline()
    return (targetsra, refgenome)


if __name__ == "__main__":
    bkp = None
    toTrim = None
    torRNA = None
    geo_id = None
    min_qual = None
    len_window = None
    min_read_len = None
    argumentList = sys.argv[1:]
    options = "b:t:r:g:q:w:l"
    try:
        arguments, values = getopt.getopt(argumentList, options)
        for currentArgument, currentValue in arguments:
            if currentArgument in ("-b"):
                print(currentArgument, currentValue)
                bkp = int(currentValue)
            elif currentArgument in ("-t"):
                print(currentArgument, currentValue)
                toTrim = currentValue
            elif currentArgument in ("-r"):
                print(currentArgument, currentValue)
                torRNA = currentValue
            elif currentArgument in ("-g"):
                print(currentArgument, currentValue)
                geo_id = str(currentValue)
            elif currentArgument in ("-q"):
                print(currentArgument, currentValue)
                min_qual = int(currentValue)
            elif currentArgument in ("-w"):
                print(currentArgument, currentValue)
                len_window = int(currentValue)
            elif currentArgument in ("-l"):
                print(currentArgument, currentValue)
                min_read_len = int(currentValue)
    except getopt.error as err:
        print(str(err))

    os.system("mkdir {}".format(geo_id))
    os.chdir(geo_id)
    sra_id, org_name = initial_step(geo_id)
    forward_reads = []
    reverse_reads = []
    db = SRAweb()
    ht2_idx_url = "ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data"
    htmlFile = getHTML(geo_id)
    genomeBuild = findGenomeBuild(htmlFile)
    print("SRP, ORG Name, build")
    print(sra_id)
    print(org_name)
    print(genomeBuild)

    if bkp == 1:
        downloadOnlySRA(db, sra_id)
    elif bkp == 2:
        download_fq_file(db, sra_id)
        srr_ids = getAllSRRIds(db, sra_id)
        forward_reads, reverse_reads = createFastqFiles(srr_ids)
        if toTrim != None or torRNA != None:
            forward_reads, reverse_reads = preprocess(
                forward_reads, reverse_reads, torRNA, toTrim, None, min_qual, min_read_len, len_window
            )  # Add sortmerna db
        qualityControl()
    elif bkp == 3:
        download_fq_file(db, sra_id)
        srr_ids = getAllSRRIds(db, sra_id)
        forward_reads, reverse_reads = createFastqFiles(srr_ids)
        if toTrim != None or torRNA != None:
            forward_reads, reverse_reads = preprocess(
                forward_reads, reverse_reads, torRNA, toTrim, None, min_qual, min_read_len, len_window
            )
        gtfName, fnaName = check_pre_built(
            genomeBuild, ht2_idx_url, org_name, False
        )  # downloadRefGenome(genomeBuild, org_name)
        qualityControl()
    elif bkp == 4:
        download_fq_file(db, sra_id)
        srr_ids = getAllSRRIds(db, sra_id)
        forward_reads, reverse_reads = createFastqFiles(srr_ids)
        if toTrim != None or torRNA != None:
            forward_reads, reverse_reads = preprocess(
                forward_reads, reverse_reads, torRNA, toTrim, None, min_qual, min_read_len, len_window
            )
        check_pre_built(genomeBuild, ht2_idx_url, org_name, True)
        qualityControl()
    elif bkp == 5:
        download_fq_file(db, sra_id)
        srr_ids = getAllSRRIds(db, sra_id)
        forward_reads, reverse_reads = createFastqFiles(srr_ids)
        if toTrim != None or torRNA != None:
            forward_reads, reverse_reads = preprocess(
                forward_reads, reverse_reads, torRNA, toTrim, None, min_qual, min_read_len, len_window
            )
        check_pre_built(genomeBuild, ht2_idx_url, org_name, True)
        startAlignment("ht2idxes", forward_reads, reverse_reads)
        qualityControl()
    elif bkp == 6:
        download_fq_file(db, sra_id)
        srr_ids = getAllSRRIds(db, sra_id)
        forward_reads, reverse_reads = createFastqFiles(srr_ids)
        if toTrim != None or torRNA != None:
            forward_reads, reverse_reads = preprocess(
                forward_reads, reverse_reads, torRNA, toTrim, None, min_qual, min_read_len, len_window
            )
        check_pre_built(genomeBuild, ht2_idx_url, org_name, True)
        output_sams = startAlignment("ht2idxes", forward_reads, reverse_reads)
        output_bams = convertSamToBam(output_sams)
        sorted_bams = sortBamFiles()
        qualityControl()
    elif bkp == 7:
        download_fq_file(db, sra_id)
        writeMessage("fq done")
        srr_ids = getAllSRRIds(db, sra_id)
        forward_reads, reverse_reads = createFastqFiles(srr_ids)
        writeMessage(forward_reads, reverse_reads)
        if toTrim != None or torRNA != None:
           forward_reads, reverse_reads = preprocess(
               forward_reads, reverse_reads, torRNA, toTrim, None, min_qual, min_read_len, len_window
           )
        writeMessage("Index Start")
        ref, annot = check_pre_built(genomeBuild, ht2_idx_url, org_name, True)
        writeMessage("Index done")
        output_sams = startAlignment("ht2idxes", forward_reads, reverse_reads)
        writeMessage("Alignment done")
        output_bams = convertSamToBam(output_sams)
        writeMessage("bams done")
        sorted_bams = sortBamFiles(output_bams)
        writeMessage("bams sorted")
        createCountMatrix(annot, sorted_bams)
        writeMessage("count matrix done")
