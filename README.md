# ASB-2024
A repository for projects developed in 'Análise de Sequências Biológicas' (ASB), in context of pursuing a Bioinformatics Degree.
These are two simple programs for fetching biological sequences from the NCBI (National Center for Biotechnology Information) database using the Entrez API, and they provide a convenient way to retrieve sequences from the NCBI database based on user-specified search queries. If one doesn't work as intended, please try the other.

*******************************************************************
Feel free to add, contribute to, reproduce, or distribute the program
*******************************************************************

Thank you for taking an interest in this project!

To use it, run it from the Terminal. It should have two arguments, "database", and "term", since the email has been set to a default value. Database is the NCBI database you want to use, and the term is your search query for said database, normally organism name and desired gene.

For example, if you want to fetch the nucleotide sequences for the CytB gene, from the species Psammodromus algirus, your terminal should look like this:

C:\Users\YourName\Program_Folder> python .\ASB_Homework_1_Lara_Evangelista.py nucleotide "Psammodromus algirus[organism], cytb[gene]"

This will retrieve the sequences and print them to the terminal. For Linux users, if the above command doesn't work, substitute "python" for "python3", in the beginning of the command.

-----------------------------------------------------------------------------------

If there is any problem executing the program as intended, please send me a message, or open a pull request.
