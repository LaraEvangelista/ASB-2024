import sys
from Bio import Entrez

def query_sequences(database, search_query):
    """
    Fetches sequences from a chosen Entrez database using the given search query.
    """
    email = 'ncbipesquisas@gmail.com'
    Entrez.email = email

    # Esearch searches and retrieves the IDs of the results
    search_handle = Entrez.esearch(db=database, term=search_query, usehistory="y", idtype="acc")
    query_results = Entrez.read(search_handle)
    search_handle.close()

    # Storing the value of the keys "WebEnv" and "QueryKey" from the "query_results" dictionary 
    webenv = query_results["WebEnv"]
    query_key = query_results["QueryKey"]

    # Efetch retrieves records in the requested format from a list of one or more primary IDs or from the user’s environment
    fetch_handle = Entrez.efetch(db=database, rettype="fasta", retmode="json", query_key=query_key, webenv=webenv)
    sequences = fetch_handle.read()
    fetch_handle.close()

    # Return the desired sequences, in FASTA format
    return sequences

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("ERROR!!! try: python3 asb_exemplo2.py <database> <search_query>")
        sys.exit(1)

    # Extract the database and search query from command line arguments
    database = sys.argv[1]
    search_query = sys.argv[2]
    sequences = query_sequences(database, search_query)

    # Write the fetched sequences
    sys.stdout.write(sequences)


    
