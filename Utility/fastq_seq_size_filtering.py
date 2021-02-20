"""Fastq Sequence Size Filtering

    Usage:
        fastq_trimming param
        fastq_trimming example
        fastq_trimming <fastq_filename> <start_range> <end_range> <matchedfilename> <nonmatchedfilename>
        fastq_trimming -h | --help

    Options:
        -h --help   Show this screen
"""

from docopt import docopt
import Fastq_class 
from generators_utilities import *

# TODO - when returning to akal PUT THOSE 3 LINES IN COMMENT AND ENABLE
# THE FIRST LINE IN THIS FILE
"""import sys
sys.path.append("/home/lammlab4/git/data_type_classes")
import Fastq_class"""


class IllegalArgumentFilterFastqFileBySize(Exception):

    def __init__(self):
        self.message = "wrong parameters entered - should be name of file and 2 numbers to indicate the size we want \"from\" and \"to\"\n the from number should be smaller then the to number"

    def __str__(self):
        return self.message

    def __repr__(self):
        return self.message


class IllegalArgumentFilterBySize(Exception):

    def __init__(self):
        self.message = "wrong parameters entered - should be list of objects ,the objects functors in a dict , 2 numbers to indicate the size we want \"from\" and \"to\"\n the from number should be smaller then the to number"

    def __str__(self):
        return self.message

    def __repr__(self):
        return self.message
#rangeFunctors={ "init":range_init,"id":range_id_getter,"start":range_start_getter,"end":range_end_getter }
# functions for functor


def initseq(seq):
    return seq


def getId(seq):
    return seq.split("\n")[0]


def start(seq):
    return 1


def end(seq):
    list = seq.split("\n")
    return len(list[1])

#=======================================

# main functions :


def filterBySize(listOfObjects, rangeFunctors, startSize, endSize):
    """
    :param listOfObjects: list of objects\ranges that we can get their range with the functors
    :param rangeFunctors: list of functions to work with the objects in the list
    :param startSize: the start range we want to filter by
    :param endSize:  the end range we want to filter by
    :return: new list of objects with only the objects that their range intersect with the range (startSize,endSize),
            and a list of objects that didnt intersect with the refRange
    """
    if (type(listOfObjects) != list) or (type(rangeFunctors) != dict) or (type(startSize) != int) or (type(endSize) != int) \
            or (startSize > endSize) or (startSize < 0):
        raise IllegalArgumentFilterBySize
    filteredList = []
    notPassList = []

    for object in listOfObjects:
        if startSize <= rangeFunctors["end"](object) <= endSize:
            filteredList.append(object)
        else:
            notPassList.append(object)

    if len(filteredList) != 0:
        lastFilterd = filteredList[len(filteredList) - 1]
        if lastFilterd[-1:] == "\n":
            filteredList.pop(len(filteredList) - 1)
            filteredList.append(lastFilterd[:-1])
    if len(notPassList) != 0:
        lastnot = notPassList[len(notPassList) - 1]
        if lastnot[-1:] == "\n":
            notPassList.pop(len(notPassList) - 1)
            notPassList.append(lastnot[:-1])

    return filteredList, notPassList


def filterFastqFileBySize(fastqFileName, smallestSize, largestSize, matchedfilename, nonmatchedfilename):
    """
    :param fastqFileName: the file name that we want to filter by seq sizes
    :param smallestSize: the smallest seq size we want to include in the result
    :param largestSize:  the largest seq size we want to include in the result
    :return: fastq file with only the seqs which sizes were in the range between smallestSize and largesrSize
            , and a file with that seqs which sizes werent in the range
    """
    if (type(fastqFileName) != str) or (type(smallestSize) != int ) or (type(largestSize) != int) or \
            (smallestSize > largestSize) or (smallestSize < 0):
        raise IllegalArgumentFilterFastqFileBySize

    filteredFile = matchedfilename
    noMatchFile = nonmatchedfilename
    open(filteredFile, 'w').close()
    open(noMatchFile, 'w').close()
    lastFlag = False
    with open(fastqFileName, "r") as fastqName:
        generator = generates4Lines(fastqName)
        try:
            lines = next(generator)
        except Exception :
            lines = None
        while lines:
            listOfItems = []
            for index in range(10000):
                if lines[0]:
                    # we put one string that represents the curr seq
                    item = "".join(lines)
                    if lines[3][-1:] != "\n":
                        lines[3] = lines[3] + "\n"
                    Fastq_class.Fastq(lines)
                    listOfItems.append(item)
                try:
                    lines = next(generator)
                except Exception :
                    lines = None
                if lines is None:
                    lastFlag = True
                    break

            filteredFileList, noMatchFileList = filterBySize(listOfItems, {
                                                             "init": initseq, "id": getId, "start": start, "end": end}, smallestSize, largestSize)

            # writing the filtered data to files in chunks
            with open(filteredFile, "a+") as file:
                if len(filteredFileList) > 0:
                    data = "".join(filteredFileList) + "\n"
                    if lastFlag:
                        if data[-1:] == "\n":
                            data = data[:-1]
                    file.write(data)

            with open(noMatchFile, "a+") as file:
                if len(noMatchFileList) > 0:
                    data = "".join(noMatchFileList) + "\n"
                    if lastFlag:
                        if data[-1:] == "\n":
                            data = data[:-1]
                    file.write(data)
    return filteredFile, noMatchFile



def param_description(): 
    """
    Description for the script
    """
    print("The parameters are\n" +
          "fastq_filename: the name of the input fastq file\n"
          "start_range: start of filtering \n"
          "end_range: end of filtering \n"
          "matchedfilename: the filename for the objects in the range \n"
          "nonmatchedfilename: the filename for the objects not in the range \n"
          "Output: filters the fastq into 2 files by intersection with the gien range")


def example_description():
    """
    Usage example for the script
    """
    print("Fastq file input:\n" +
          "@ABC\nGATCATC\n+\n!#!!$!@\n" +
          "@BCD\nGATC\n+\n\"\"\"\"\n" +
          "@GHJ\nTCGAGCAG\n+\n####!!!!\n\n" +
          "Calling the script:\n" +
          "fastq_seq_size_filtering fastq_input 5 7 fastq_matched fastq_nonmatched\n\n" +
          "Output fastq_matched:\n" +
          "@ABC\nGATCATC\n+\n!#!!$!@\n" +
          "Output fastq_nonmatched:\n" +
          "@BCD\nGATC\n+\n\"\"\"\"\n" +
          "@GHJ\nTCGAGCAG\n+\n####!!!!")


#=======================================

if __name__ == "__main__":
    """
    usage : this script receives a fastq file , a start range (in size of sequance) ,
    a end range , filename for matched line of the requested size according to the given range
    and filename for the lines in the fastq that didnt match the given range.

    the output will be given in the 2 files that their name was supplied according to the explanetion above. 
    """
    import sys
    arguments = docopt(__doc__)

    if arguments["param"]:
        param_description()
        sys.exit()

    if arguments["example"]:
        example_description()
        sys.exit()

    fastq_filename = arguments["fastq_filename"]
    start_range = arguments["start_range"]
    end_range = arguments["end_range"]
    matchedfilename = arguments["matchedfilename"]
    nonmatchedfilename = arguments["nonmatchedfilename"]

    filterFastqFileBySize(fastq_filename, int(start_range), int(
        end_range), matchedfilename, nonmatchedfilename)
Ã«† ‘¡Î∆<´»•ÛŒ‘¡Î P