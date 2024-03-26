import numpy as np
import pandas as pd
import os
import csv
from rpy2 import robjects
import rpy2.robjects.packages as rpackages

def common_member(a, b):
    """
    Returns a list of common elements between two lists.

    Args:
    - a: First list.
    - b: Second list.

    Returns:
    - result: List of common elements between the two input lists.
    """
    result = [i for i in a if i in b]
    return result

def cluster_metrics(local_seg_predictions, true_bp):
    """
    Computes cluster metrics such as recall and precision.

    Args:
    - local_seg_predictions: Dictionary of local segment predictions.
    - true_bp: Dictionary of true breakpoints.

    Returns:
    - recall: List of recall values.
    - precision: List of precision values.
    """
    recall = []
    precision = []
    for i in range(1, 1+len(local_seg_predictions)):
        true_positive = len(common_member(local_seg_predictions['cell'+str(i)],true_bp['cell'+str(i)]))
        current_recall = true_positive / len(true_bp['cell'+str(i)])
        current_precision = true_positive / len(local_seg_predictions['cell'+str(i)])
        recall.append(current_recall)
        precision.append(current_precision)

    print("Mean Recall: "+ str(np.mean(recall)))
    print("Mean Precision: " + str(np.mean(precision)))
    return recall,precision

def save_clusting_dict(clustering_dict, file_name):
    """
    Saves the clustering dictionary to a file.

    Args:
    - clustering_dict: Dictionary of clustering results.
    - file_name: Name of the file to save the dictionary.
    """
    with open(file_name,"w+") as file:
        for _, cell in enumerate(clustering_dict):
            file.write('\t'.join([cell] + [str(i) for i in clustering_dict[cell]]) + '\n')
        file.close()


def local_segmentation_generator(filename, readcount_file_path, num_chrom, num_bin_per_chrom=[], num_permutation=1000, alpha=0.05, min_width=3, k_max=10):
    """
    Generates local segmentation using the CBS algorithm.

    Args:
    - filename: Path to the output local segmentation file.
    - readcount_file_path: Path to the readcount file.
    - num_bin_per_chrom: List of number of bins per chromosome. 
    - num_chrom: Number of chromosomes.
    - num_permutation: Number of permutations.
    - alpha: Alpha value for CBS algorithm.
    - min_width: Minimum width for CBS algorithm.
    - k_max: Maximum number of change-points for CBS algorithm.
    """
    ## check if the output file already exists
    if os.path.exists(filename):
        os.remove(filename)
    
    robjects.r['library']('utils')
    readcount_data = robjects.r['read.csv'](readcount_file_path, sep='\t', header=True)

    ## convert the input dataset to CNA object in R
    input_dimension = robjects.r('dim(%s)'%readcount_data.r_repr())
    ## if num_bin_per_chrom is empty,
    ## default number of bin per chromosome is the total number of bins divided by the number of chromosomes
    if len(num_bin_per_chrom) == 0:
        num_bin_per_chrom = (input_dimension[1] - 1) // num_chrom
        rcode_for_cna_object = 'CNA(t(%s[,2:%d]), rep(1,each=%d),1:%d)'%(readcount_data.r_repr(), input_dimension[1], input_dimension[1]-1, input_dimension[1]-1)
    else:
        try: 
            assert len(num_bin_per_chrom) == num_chrom, "number of bins per chromosome is not equal to the number of chromosomes"
        except AssertionError as error:
            print(error)
            return
        
        bin_index_string = ','.join(map(str, range(1,num_chrom+1)))
        num_bin_per_chrom_string = ','.join(map(str,num_bin_per_chrom))
        rcode_for_cna_object = 'CNA(t(%s[,2:%d]), rep(c(%s),times=c(%s)),1:%d)'%(readcount_data.r_repr(), input_dimension[1], bin_index_string, num_bin_per_chrom_string, (input_dimension[1]-1))
    robjects.r['library']('DNAcopy')
    cna_object = robjects.r(rcode_for_cna_object)

    ## runing the cbs segmentation algorithm on CNA object
    rcode_for_seg = 'segment(%s, nperm=%d, alpha=%f, min.width=%d, kmax=%d)'%(cna_object.r_repr(),num_permutation, alpha, min_width, k_max)
    local_segmentation = robjects.r(rcode_for_seg)
    local_segmentation_output = robjects.r('%s$output'%(local_segmentation.r_repr()))

    rcode_for_convert_seg_output_to_list = '''local_segmentation_result <- list()
                                            for (i in 1:%d){
                                                current_output <- %s[%s$ID == paste("Sample.",i, sep = ""),]
                                                local_segmentation_result[[i]] <- current_output$loc.end
                                                }'''%(input_dimension[0], local_segmentation_output.r_repr(), local_segmentation_output.r_repr())
    robjects.r(rcode_for_convert_seg_output_to_list)

    rcode_save_the_final_result = '''for (i in 1:length(local_segmentation_result)){
            write(paste(c(paste("cell",i,sep=""),local_segmentation_result[i][[1]]), collapse = "\t"), file="%s",append=TRUE)
        }'''%(filename)
    robjects.r(rcode_save_the_final_result)
    return



def local_segmentation_generator_old(filename, readcount_file_path, num_chrom, num_permutation=1000, alpha=0.05, min_width=5, k_max=10):
    """
    Generates local segmentation using the CBS algorithm.

    Args:
    - filename: Path to the output local segmentation file.
    - readcount_file_path: Path to the readcount file.
    - num_chrom: Number of chromosomes.
    - num_permutation: Number of permutations.
    - alpha: Alpha value for CBS algorithm.
    - min_width: Minimum width for CBS algorithm.
    - k_max: Maximum number of change-points for CBS algorithm.
    """
    robjects.r['library']('utils')
    readcount_data = robjects.r['read.csv'](readcount_file_path, sep='\t', header=True)

    ## convert the input dataset to CNA object in R
    input_dimension = robjects.r('dim(%s)'%readcount_data.r_repr())
    #if num_chrom == 1:
    #    num_bin_per_chrom = (input_dimension[1] - 1) // num_chrom
    #else:
    #    num_bin_per_chrom = input_dimension[1] // num_chrom
    num_bin_per_chrom = (input_dimension[1] - 1) // num_chrom
    rcode_for_cna_object = 'CNA(t(%s[,2:%d]), rep(1:%d,each=%d),rep(1:%d,%d))'%(readcount_data.r_repr(), input_dimension[1], num_chrom, num_bin_per_chrom, num_bin_per_chrom, num_chrom)
    robjects.r['library']('DNAcopy')
    cna_object = robjects.r(rcode_for_cna_object)

    ## runing the cbs segmentation algorithm on CNA object
    rcode_for_seg = 'segment(%s, nperm=%d, alpha=%f, min.width=%d, kmax=%d)'%(cna_object.r_repr(),num_permutation, alpha, min_width, k_max)
    local_segmentation = robjects.r(rcode_for_seg)
    local_segmentation_output = robjects.r('%s$output'%(local_segmentation.r_repr()))

    rcode_for_convert_seg_output_to_list = '''local_segmentation_result <- list()
                                            for (i in 1:%d){
                                                current_output <- %s[%s$ID == paste("Sample.",i, sep = ""),]
                                                current_seg <- current_output$loc.end + %d * (current_output$chrom-1)
                                                local_segmentation_result[[i]] <- current_seg
                                                }'''%(input_dimension[0], local_segmentation_output.r_repr(),
                                                    local_segmentation_output.r_repr(), num_bin_per_chrom)
    robjects.r(rcode_for_convert_seg_output_to_list)

    rcode_save_the_final_result = '''for (i in 1:length(local_segmentation_result)){
            write(paste(c(paste("cell",i,sep=""),local_segmentation_result[i][[1]]), collapse = "\t"), file="%s",append=TRUE)
        }'''%(filename)
    robjects.r(rcode_save_the_final_result)
    return