import sys
import argparse
from Bio import Entrez

def query_sequences(database, query):
    
    """Fetches sequences from a chosen Entrez database, using the given search query, and prints them to the Terminal
        Takes database and query as input, returns a list of sequences
    """
    
    # Set the Entrez email parameter 
    email = 'ncbipesquisas@gmail.com'
    Entrez.email = email
    # Esearch searches and retrieves the IDs of the results, which we store for later use by enabling the history feature
    search = Entrez.esearch(db = database, term = query, usehistory = "y", idtype = "acc")
    # Read Parses the XML results returned by any of the above functions. Returns a Python dictionary (in this case) 
    query_results = Entrez.read(search)
    # If a file is not properly closed, it can lead to memory leaks, which can slow down your program or even cause it to crash.
    search.close()
    # Storing the value of the keys "WebEnv" and "QueryKey" from the "query_results" dictionary 
    webenv = query_results["WebEnv"]
    query_key = query_results["QueryKey"]
    # efetch Retrieves records in the requested format from a list of one or more primary IDs or from the user’s environment
    search = Entrez.efetch(db = database, rettype = "fasta", retmode = "json", query_key = query_key, webenv = webenv)
    # reads the result of the fetch, (the desired sequence), and appends it to the empty sequences list
    sequences = search.read()
    # closes the file, for the above explained reasons
    search.close()
    #return the desired sequences, in FASTA format
    return sequences

def arguments_parser():
    
    """Parses the command line arguments, after the user specifies the database and search query"""
    
    # The support for command-line interfaces is built around an instance of argparse.ArgumentParser. It is a container for argument specifications 
    parser = argparse.ArgumentParser(description="Retrieve the sequences from a search query in the Entrez database. Please input the database and the search query")
    # The ArgumentParser.add_argument() method attaches individual argument specifications to the parser
    # I chose two argumentsm one for the database to query, and another for the search query itself
    parser.add_argument("database", help="\nInput the Entrez database to query (nucleotide, protein, genome or gene)")
    parser.add_argument("search_query", help="\nInput your search query")
    
    # The ArgumentParser.parse_args() method runs the parser and places the extracted data in a argparse.Namespace object
    return parser.parse_args()
    
""" Run the functions """

if __name__ == "__main__":
    # Function that scans the users input, storing it for later use
    args = arguments_parser()
    # Main function of the program, will fetch sequences from the NCBI database inputted, using the search term provided
    # Makes use of argparse.Namespace ('.database' and '.search_query'), which stored the user's input in the previous function
    sequences = query_sequences(args.database, args.search_query)

    # Will write the fetched sequences in the 'sequences' list to the standard output, the Terminal
    for sequence in sequences:
        sys.stdout.write(sequence)
