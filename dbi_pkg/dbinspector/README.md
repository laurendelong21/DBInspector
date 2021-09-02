## package modules - Cross Verification of Protein Information

### startup

The `startup.py` file:  
initiates the cache directories

### parse
This module calls `startup.py`, which initiates the cache directories.

Calling the parse wrapper function `parse_all` initiates the data downloads from RefSeq's and 
UniProt's API - if necessary. The data is then processed into json files where a single entry has the metadata info assigned to it's accession ID.

+ initiate downloads from RefSeq and UniProt
+ parses data from both databases
+ write parsed results to JSON files to cache with the following structures:

`refseq.json`
```python
{refseq accession.version:  
    {symbol: list(corresponding gene symbols),
     seq: str(amino acid sequence),
     uniprot: str(corresponding uniprot ID)}
```
`uniprot.json`:
```python
{uniprot ID:
    {symbol: list(corresponding gene symbols),
     seq: str(amino acid sequence),
     refseq: list(corresponding refseq accession.versions)}
```

---
### map
Functions in the module `map` can be accessed via the wrapper function `find_entries`. This was designed for 
interaction with `compare`, which in turn interacts with the GUI and CLI.  
`find_entries` takes a UniProt or RefSeq accession ID, or gene symbol as input and fetches the corresponding entries from **both** databases.  
  
The return of a successful query will have the structure that is shown below, with nested dictionaries and lists to reliably distinguish between several entries from seperate databases. This example shows the structure of a result with one UniProt and two RefSeq entries. 
```python
results = {'UniProt': [{'symbol': list(str),
                        'UniProt ID': str,
                        'RefSeq ID': list(str),
                        'sequence': str
                        }],
           'RefSeq': [{'symbol': list(str),
                       'UniProt ID': str,
                       'RefSeq ID': str,
                       'sequence': str
                       }
                      {'symbol': list(str),
                       'UniProt ID': str,
                       'RefSeq ID': str,
                       'sequence': str
                       }
           }
```


     
