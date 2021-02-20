import logging
from Utility.parallel_generator import *
from Utility.generators_utilities import class_generator
from Utility.Pileup_class import Pileup_line
import subprocess
from collections import Counter
import traceback

def create_combined_negative(negative_pileups, filtered_negatives):    
            #combine all negatives into single filtered negative, to save runtime in snp_removal
            #parameters for parallel generator
            neg_obj_list = [open(p) for p in negative_pileups]
            neg_gen_list = []
            get_pos_and_id = lambda x: (x.reference_id, x.gene_pos)
            
            #filter out bad lines from negatives
            for file_obj in neg_obj_list:
                neg_gen_list.append(class_generator(Pileup_line, file = file_obj))
            with open(filtered_negatives, 'w') as out:
                parallel_gen = parallel_generator(neg_gen_list, [get_pos_and_id for i in range(len(neg_gen_list))])
                for parallel_line_list in parallel_gen: 
                    snp = False
                    for pileup_list in parallel_line_list:
                        if pileup_list is not None:
                           if pileup_list[0].is_with_any_change():
                                snp = True
                                
                                break
                    if not snp:
                        line_to_write = next(item for item in parallel_line_list if item is not None)
                        out.write(str(line_to_write[0]) + "\n")
            for obj in neg_obj_list:
                obj.close()
        
        
def subtract_pileup_line(pileup_to_filter, pileup_to_subtract):
    reads_to_remove = set() #create empty set for indices to remove
    base_count = int(pileup_to_filter._base_count)
    for i in range(base_count):
        current_read = pileup_to_filter.reads_string[i]
        if (current_read in (pileup_to_subtract.reads_string)) and (current_read in "agctAGCT"):
            reads_to_remove.add(i)
    if(len(reads_to_remove) == base_count):
        return None
    new_read_string = "".join([char for idx, char in enumerate(pileup_to_filter.reads_string) if idx not in reads_to_remove])
    new_qual_string = "".join([char for idx, char in enumerate(pileup_to_filter.quality_string) if idx not in reads_to_remove])
    pileup_to_filter.reads_string = new_read_string
    pileup_to_filter.quality_string = new_qual_string
    pileup_to_filter._base_count = str(base_count - len(reads_to_remove))
    return str(pileup_to_filter)

def create_consensus_pileup_line(pileup_lines, k, number_of_pileups):
    non_empty_pileup_list_list = [item for item in pileup_lines if item is not None]
    number_of_non_empty = len(non_empty_pileup_list_list)
    #if the read isn't in enough pileups, we discard it. (times_read_is_in_line/number_of_pileups >= k)
    seen_in_line_dict = Counter()
    nucleotide_line_occurence_list = {nucleotide:[0 for _ in range(number_of_non_empty)] for nucleotide in "AGCTagct"}
    #seen_in_line_dict makes sure we don't add the same read from the same line more than
    #once, so we can find out how many lines the read is in
    for i in range(number_of_non_empty):
        seen_in_line_dict.clear()
        current_line = non_empty_pileup_list_list[i]
        if len(current_line) > 1:
            print("too many pileups in a list: ", line_list)
        else:
            nuc_read_list = [read for read in current_line[0].reads_string if read in "AGCTagct"]
            for read in nuc_read_list:
                nucleotide_line_occurence_list[read][i] += 1
                if seen_in_line_dict[read] is 0:
                    seen_in_line_dict[read] += 1
    #create combined read:
    big_line = non_empty_pileup_list_list[0][0]
    big_line._base_count = '0'
    big_line.reads_string = ""
    big_line.quality_string = ""
    #create string of number of shared nucleotide reads:
    for nucleotide,lines in nucleotide_line_occurence_list.items(): #['G':[3,1,4,2]
        max_num_reads_above_consensus = 0
        for i in range(max(lines),0,-1):
            number_of_occurences_in_line = sum(x >= i for x in lines)
            if number_of_occurences_in_line/number_of_pileups >= k:
                big_line._base_count = str(int(big_line._base_count) + number_of_occurences_in_line)
                big_line.reads_string += nucleotide*number_of_occurences_in_line
                big_line.quality_string += 'F'*number_of_occurences_in_line
                break #exits to looping over other nucleotides
    for pileup_list in non_empty_pileup_list_list:
        pileup = pileup_list[0]
        for i in range(int(pileup._base_count)):
            read = pileup.reads_string[i]
            if read in "AGCTagct":
                continue
            else:
                big_line.reads_string += read
                big_line.quality_string += pileup.quality_string[i]
                big_line._base_count = str(int(big_line._base_count)+1)
    return str(big_line)    
    

    
def create_consensus_file_by_reads(pileup_filename_list, k, consensus_filename, sorted_input=True):
    '''
    :param pileup_filename_list: list of pileup files to create a consensus list of pileup lines out of
    :param k: the percentage of pileup files a line must be included in to be within the consensus. 0<k<1
    :param consensus_filename: name of output file
    :param sorted_input: binary flag, set to False if input pileups are not sorted
    :return: consensus file
    '''
    sorted_pileups = []
    if not sorted_input:
        for pileup in pileup_filename_list:
            output_pileup = pileup + "_sorted.pileup"
            command = f"sed -r '/^[\t]*$/d' <{pileup} | LC_ALL=C sort -k1,1 -k2,2n -o {output_pileup} "
            try:
                outout=subprocess.check_output(command,shell=True)
            except subprocess.CalledProcessError as e:
                logging.error(f"creating consensus file for {pileup_filename_list} failed")
                raise e
            sorted_pileups.append(output_pileup)
    else:
        sorted_pileups = pileup_filename_list

    pileup_files = [open(pileup, 'r') for pileup in sorted_pileups]

    try:
        generators = [class_generator(Pileup_line, file=file) for file in pileup_files]
        get_pos_and_id = lambda x: (x.reference_id, x.gene_pos)
        parallel_gen = parallel_generator(generators, [get_pos_and_id for i in range(len(generators))])
        number_of_pileups = len(pileup_files)
        with open(consensus_filename, 'w') as out:
            for pileuplist_list in parallel_gen:
                hits = sum(x is not None for x in pileuplist_list) #number of non-None lists, the number of pileups which have the entry
                if (hits/number_of_pileups) >= k: #if the position isn't in enough files, no need to continue
                    try:
                        output_line = create_consensus_pileup_line(pileuplist_list,k,number_of_pileups)
                    except Exception as e:
                        traceback.print_exc()
                        raise e                    
                    out.write(output_line + "\n")
    finally:
        for file in pileup_files:
            file.close()
        return consensus_filename



    
def subtract_sites(pileup_filename, negative_pileup_list, output_file=None, sorted_input=False):
    """
    Reads that are in the negative_pileup_list are removed from pileup_filename.
    :param pileup_filename: pileup to be filtered
    :param negative_pileup_list: list of negative pileups to filter positive pileup by
    :param output_file: file to write output to
    :param sorted_input: binary flag, set to False if input pileups are not sorted
    :return
    """

    # sorting pileups if necessary and preparing filenames
    sorted_negative_pileups = []
    if not sorted_input:
        output_pileup = pileup_filename + "_sorted.pileup"
        command = f"sed -r '/^[\t]*$/d' <{pileup_filename} | LC_ALL=C  sort -k1,1 -k2,2n -o {output_pileup} "
        try:
            outout=subprocess.check_output(command,shell=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"pileup conversions of {sam_name} failed")
            raise e
        sorted_positive_pileup = pileup_filename + "_sorted.pileup"
        for pileup in negative_pileup_list:
            output_pileup = pileup + "_sorted.pileup"
            command = f"sed -r '/^[\t]*$/d' <{pileup} | LC_ALL=C sort -k1,1 -k2,2n -o {output_pileup} "
            try:
                outout=subprocess.check_output(command,shell=True)
            except subprocess.CalledProcessError as e:
                logging.error(f"pileup conversions of {sam_name} failed")
                raise e
            sorted_negative_pileups.append(pileup + "_sorted.pileup")
    else:
        sorted_positive_pileup = pileup_filename
        sorted_negative_pileups = negative_pileup_list
    if output_file is None:
        splited = sorted_positive_pileup.split(".")
        output_file = ".".join(splited[:-1]) + "_no-snp." + splited[-1]

    # parameters for parallel generator
    neg_obj_list = [open(pile) for pile in sorted_negative_pileups]
    neg_gen_list = []
    get_pos_and_id = lambda x: (x.reference_id, x.gene_pos)

    # logic of the function
    try:
        for file_obj in neg_obj_list:
            neg_gen_list.append(class_generator(Pileup_line, file=file_obj))
        with open(sorted_positive_pileup, 'r') as pileup_obj, open(output_file, 'w') as out:
            pos_gen = class_generator(Pileup_line, file=pileup_obj)
            parallel_gen = parallel_generator([pos_gen, *neg_gen_list],
                                              [get_pos_and_id for i in range(len([pos_gen, *neg_gen_list]))])
            for listlist in parallel_gen:  # listlist: list of lists, each entry in listlist is a list of items from one generator
                if listlist[0] is not None:  # listlist[0]: list of pileup_lines(should be 1 line) from the input pileup
                    output_line = None
                    for pileup_list in listlist[1:]:
                        if pileup_list is not None:
                            if (len(listlist[0]) is not 1) or (len(pileup_list) is not 1):
                                #for l in listlist:
                                #    for ll in l:
                                #        print(get_pos_and_id(ll))
                                print("pos list: ", listlist[0], "\nneg lists: ", listlist[1:])
                                raise Exception
                            if output_line is None:
                                pileup_to_filter = listlist[0][0]
                            elif output_line == "":
                                
                                    continue
                            else:
                                pileup_to_filter = Pileup_line(output_line)
                            output_line = subtract_pileup_line(pileup_to_filter, pileup_list[0])
                            if output_line is None:
                                output_line = ""
                            #print(output_line)
                    if output_line is not None and output_line != "":
                        out.write(str(output_line) + "\n")

    except Exception as e:
        print(e)

    finally:
        for fd in neg_obj_list:
            fd.close()
        return output_file

def remove_non_change(pileup_to_filter):
    reads_to_remove = set() #create empty set for indices to remove
    base_count = int(pileup_to_filter._base_count)
    for i in range(base_count):
        current_read = pileup_to_filter.reads_string[i]
        if current_read not in "agctAGCT":
            reads_to_remove.add(i)
    if(len(reads_to_remove) == base_count):
        return None
    new_read_string = "".join([char for idx, char in enumerate(pileup_to_filter.reads_string) if idx not in reads_to_remove])
    new_qual_string = "".join([char for idx, char in enumerate(pileup_to_filter.quality_string) if idx not in reads_to_remove])
    pileup_to_filter.reads_string = new_read_string
    pileup_to_filter.quality_string = new_qual_string
    pileup_to_filter._base_count = str(base_count - len(reads_to_remove))
    return str(pileup_to_filter)

    
def filter_non_change_reads(pileup_file, output_file):
    with open(pileup_file, 'r') as pileup_obj, open(output_file, 'w') as out:
        line_gen = class_generator(Pileup_line, file=pileup_obj)
        for line in line_gen:
            new_line = remove_non_change(line)
            if new_line is not None:
                out.write(new_line+'\n')

def return_nuc_change_count(pileup_line, nucleotide):
    base_count = int(pileup_line._base_count)
    return pileup_line.reads_string.count(nucleotide)
    

def calculate_edit_percent(edited_pileup, full_pileup, output_file, edit_type=None):
    '''
    
    :param edit_type: if None, all edit percentage, if "G", G editing (G in pileup read)
    '''
    with open(edited_pileup, 'r') as edited, open(full_pileup, 'r') as full, open(output_file, 'w') as out:
        edit_gen = class_generator(Pileup_line, file=edited)
        full_gen = class_generator(Pileup_line, file=full)
        get_pos_and_id = lambda x: (x.reference_id, x.gene_pos)
        parallel_gen = parallel_generator([edit_gen, full_gen], [get_pos_and_id, get_pos_and_id])
        for listlist in parallel_gen:
            if listlist[0] is not None and listlist[1] is not None:
                if len(listlist[0]) is not 1 or len(listlist[1]) is not 1:
                    print("edited: ", listlist[0], "\nfull: ", listlist[1])
                else:
                    output_line = listlist[0][0]
                    if edit_type is None:
                        out.write(str(output_line)+"\t"+str((int(listlist[0][0]._base_count))/(int(listlist[1][0]._base_count)))+'\t'+listlist[1][0]._base_count+"\n")
                    else:
                        edit_amount = return_nuc_change_count(output_line, edit_type)
                        out.write(str(output_line)+"\t"+str(edit_amount/int(listlist[1][0]._base_count))+'\t'+listlist[1][0]._base_count+"\n")

def calculate_edit_count(edited_pileup, full_pileup, output_file, edit_type=None):
    '''
    
    :param edit_type: if None, all edit percentage, if "G", G editing (G in pileup read)
    '''
    with open(edited_pileup, 'r') as edited, open(full_pileup, 'r') as full, open(output_file, 'w') as out:
        edit_gen = class_generator(Pileup_line, file=edited)
        full_gen = class_generator(Pileup_line, file=full)
        get_pos_and_id = lambda x: (x.reference_id, x.gene_pos)
        parallel_gen = parallel_generator([edit_gen, full_gen], [get_pos_and_id, get_pos_and_id])
        for listlist in parallel_gen:
            if listlist[0] is not None and listlist[1] is not None:
                if len(listlist[0]) is not 1 or len(listlist[1]) is not 1:
                    print("edited: ", listlist[0], "\nfull: ", listlist[1])
                else:
                    output_line = listlist[0][0]
                    if edit_type is None:
                        out.write(str(output_line)+"\t"+str((int(listlist[0][0]._base_count))/(int(listlist[1][0]._base_count)))+'\t'+listlist[1][0]._base_count+"\n")
                    else:
                        edit_amount = return_nuc_change_count(output_line, edit_type)
                        out.write(str(output_line)+"\t"+str(edit_amount)+'\t'+listlist[1][0]._base_count+"\n")
                        
def calculate_edit_type_count(edited_pileup, output_file, reference_nucleotide, read_nucleotide):
    '''
    
    '''
    with open(edited_pileup, 'r') as edited, open(output_file, 'w') as out:
        edit_gen = class_generator(Pileup_line, file=edited)
        for pileup_line in edit_gen:
            if pileup_line.reference is from_nc:
                edit_type = reference_nucleotide+read_nucleotide
                edit_amount = return_nuc_change_count(pileup_line, to_nc.upper())
                edit_amount += return_nuc_change_count(pileup_line, to_nc.lower())
                if edit_amount > 0:
                    out.write(str(pileup_line)+"\t"+str(edit_type)+'\t'+str(edit_amount)+"\n")
             
   
def get_pileup_positions_and_edit_percent_from_file_list(pileup_filename, pileup_list,from_nc, to_nc, output_file=None, sorted_input=False):
    """
    Reads that are in the negative_pileup_list are removed from pileup_filename.
    :param pileup_filename: pileup to be filtered
    :param negative_pileup_list: list of negative pileups to filter positive pileup by
    :param output_file: file to write output to
    :param sorted_input: binary flag, set to False if input pileups are not sorted
    :return
    """

    # sorting pileups if necessary and preparing filenames
    sorted_negative_pileups = []
    if not sorted_input:
        output_pileup = pileup_filename + "_sorted.pileup"
        command = f"sed -r '/^[\t]*$/d' <{pileup_filename} | LC_ALL=C  sort -k1,1 -k2,2n -o {output_pileup} "
        try:
            outout=subprocess.check_output(command,shell=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"pileup conversions of {sam_name} failed")
            raise e
        sorted_positive_pileup = pileup_filename + "_sorted.pileup"
        for pileup in pileup_list:
            output_pileup = pileup + "_sorted.pileup"
            command = f"sed -r '/^[\t]*$/d' <{pileup} | LC_ALL=C sort -k1,1 -k2,2n -o {output_pileup} "
            try:
                outout=subprocess.check_output(command,shell=True)
            except subprocess.CalledProcessError as e:
                logging.error(f"pileup conversions of {sam_name} failed")
                raise e
            sorted_negative_pileups.append(pileup + "_sorted.pileup")
    else:
        sorted_positive_pileup = pileup_filename
        sorted_negative_pileups = pileup_list
    if output_file is None:
        splited = sorted_positive_pileup.split(".")
        output_file = ".".join(splited[:-1]) + "_no-snp." + splited[-1]

    # parameters for parallel generator
    neg_obj_list = [open(pile) for pile in sorted_negative_pileups]
    neg_gen_list = []
    get_pos_and_id = lambda x: (x.reference_id, x.gene_pos)
    get_pos_and_id_from_str = lambda x: (x.split("\t")[0],int(x.split("\t")[1]))
    # logic of the function
    try:
        for file_obj in neg_obj_list:
            neg_gen_list.append(class_generator(Pileup_line, file=file_obj))
        with open(sorted_positive_pileup, 'r') as pileup_obj, open(output_file, 'w') as out:
            pos_gen = class_generator(str, file=pileup_obj)
            parallel_gen = parallel_generator([pos_gen, *neg_gen_list],
                                              [get_pos_and_id_from_str] + [get_pos_and_id for i in range(len(neg_gen_list))])
            for listlist in parallel_gen:  # listlist: list of lists, each entry in listlist is a list of items from one generator
                if listlist[0] is not None:  # listlist[0]: list of pileup_lines(should be 1 line) from the input pileup
                    output_line = None
                    for pileup_list in listlist[1:]:
                        if pileup_list is not None:
                            if (len(listlist[0]) is not 1) or (len(pileup_list) is not 1):
                                print("pos list: ", listlist[0], "\nneg lists: ", listlist[1:])
                                raise Exception
                            
                            if output_line is None:
                                output_line = pileup_list[0]
                            else:
                                output_line._base_count = str(int(pileup_list[0]._base_count) + int(output_line._base_count))
                                output_line.reads_string += pileup_list[0].reads_string
                                output_line.quality_string += pileup_list[0].quality_string
                            
                    if output_line is not None and output_line != "":
                        total_to_nc = 0
                        for read in output_line.reads_string:
                            if read.upper() == to_nc.upper():
                                total_to_nc += 1
                        if output_line.reference == from_nc:
                            out.write(str(output_line) + "\t" + str(total_to_nc/int(output_line._base_count)) + "\n")

    except Exception as e:
        print(e)
        traceback.print_exc()

    finally:
        for fd in neg_obj_list:
            fd.close()
        return output_file


             
if __name__=="__main__":
    root_dir = "/Data/user/mahmud/resik/mirna_sams/"
    pos_pileups = ["11321_N2_Emb_Primary_siRNA_biorep2_1B.filtered_norep_pileup_generation_sorted_pileup.pileup"
        ,"11331_N2_Emb_Primary_siRNA_biorep3_1C.filtered_norep_pileup_generation_sorted_pileup.pileup"
        ,"11341_N2_Emb_Primary_siRNA_biorep4_1D.filtered_norep_pileup_generation_sorted_pileup.pileup"]
    neg_pileups = ["21311_BB21_Emb_Primary_siRNA_biorep1_5A.filtered_norep_pileup_generation_sorted_pileup.pileup"
        ,"21321_BB21_Emb_Primary_siRNA_biorep2_5B.filtered_norep_pileup_generation_sorted_pileup.pileup"
        ,"21331_BB21_Emb_Primary_siRNA_biorep3_5C.filtered_norep_pileup_generation_sorted_pileup.pileup"]
    pos_pileup_paths = [root_dir+"pileups/"+pos for pos in pos_pileups]
    neg_pileup_paths = [root_dir+"pileups/"+neg for neg in neg_pileups]
    consensus_filename = root_dir+"positive_consensus_sites.pileup"
    subtracted_filename =  consensus_filename+".subtracted"

    #create positive consensus:
    create_consensus_file_by_reads(pos_pileup_paths, 0.5, consensus_filename, sorted_input=True)
    #subtract sites
    subtract_sites(consensus_filename, neg_pileup_paths, subtracted_filename,True)
    
    #filter out non changed reads:
    
    filter_non_change_reads(subtracted_filename, subtracted_filename+".edit_sites")
    
    #create edit count files 
    for from_nc in "AGCT":
        for to_nc in "AGCT":
            if from_nc == to_nc:
                continue
            else:
                output_filename = root_dir + "edit_counts/positive_edit_counts" + "." + from_nc + "_" + to_nc + "_edit_count"
                calculate_edit_type_count(subtracted_filename+".edit_sites", output_filename, from_nc, to_nc)
    #create edit percent files from edit count files      
    for from_nc in "AGCT":
        for to_nc in "AGCT":
            if from_nc == to_nc:
                continue
            else:
                input_filename = root_dir + "edit_counts/positive_edit_counts" + "." + from_nc + "_" + to_nc + "_edit_count"
                output_filename = root_dir + "edit_percents/positive_edit_pcts" + "." + from_nc + "_" + to_nc + "_edit_pct"
                get_pileup_positions_and_edit_percent_from_file_list(input_filename, pos_pileup_paths, from_nc, to_nc, output_filename, True)

        
        
