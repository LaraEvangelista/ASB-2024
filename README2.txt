This is a simple program for fetching biological sequences from the NCBI (National Center for Biotechnology Information) database using the Entrez API. 
It is designed to be used in the discipline of Biological Sequence Analysis.
It provides a convenient way to retrieve sequences from the NCBI database based on user-specified search queries.

How to use:

Run it from the Terminal. 
It should have two arguments, "database", and "term", since the email has been set to a default value.
Database is the NCBI database you want to use, and the term is your search query.

-<database>: The name of the Entrez database to search ("nucleotide", "protein", "genome", "gene").
-<search_query>: The search term to use in the query.

### Example

To fetch nucleotide sequences for the CytB gene from the Passer domesticus organism, you can use the following command:

python3 asb_exemplo2.py nucleotide "Passer domesticus[organism] CytB[gene]"

This will retrieve the sequences and print them to the terminal!

-------------------------------------------------------------------------------------------------------
                                                Thank you!
-------------------------------------------------------------------------------------------------------




