# CBB script from DATOS

from Bio import Entrez, SeqIO
import sys
Entrez.email = 'hola@example.com'
taxon = sys.argv[1]
# Retrieve taxonomy ID by taxon name
handle = Entrez.esearch (db="Taxonomy", term=f'{taxon}[All Names]', retmode="xml")
record = Entrez.read(handle)  # retrieve taxon ID
taxon_id = record['IdList']
print(taxon_id[0])



