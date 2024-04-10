import argparse
import csv
import gzip
import json


# The complex tsv file for human was downloaded from EBI complex portal:http://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv
# An example line with header:
# Complex ac	Recommended name	Aliases for complex	Taxonomy identifier	Identifiers (and stoichiometry) of molecules in complex	Evidence Code	Experimental evidence	Go Annotations	Cross references	Description	Complex properties	Complex assembly	Ligand	Disease	Agonist	Antagonist	Comment	Source	Expanded participant list
# CPX-1	SMAD2-SMAD3-SMAD4 complex	SMAD2/SMAD3/SMAD4 transcription factor complex	9606	P84022(1)|Q13485(1)|Q15796(1)	ECO:0005547(biological system reconstruction evidence based on inference from background scientific knowledge used in manual assertion)	-	\
# GO:0071144(heteromeric SMAD protein complex)|GO:0003690(double-stranded DNA binding)|GO:0003700(DNA-binding transcription factor activity)|GO:0006355(regulation of DNA-templated transcription)|GO:0032924(activin receptor signaling pathway)|GO:0007179(transforming growth factor beta receptor signaling pathway)	\
# reactome:R-HSA-9736938(identity)|reactome:R-HSA-9736929(identity)|pubmed:35359452(see-also)|pubmed:16322555(see-also)|complex portal:CPX-1(complex-primary)|wwpdb:1U7V(subset)	\
# A transcription factor complex which binds to the promoters of target genes and recruits co-activators and histone acetyltransferases, such as p300, CBP and P300/CBP-associated factor, \
# facilitating transcription. In response to TGF-beta/activin-family protein binding, TGF-beta type II receptors phosphorylate TGF-beta type I receptors (ALK4, 5 and 7) which in turn phosphorylates SMAD2 on two Ser-465 and Ser-467, and SMAD3 on Ser-423 and Ser-425. \
# This enables binding to SMAD4 to form heteromeric SMAD complexes that enter the nucleus to initiate gene transcription. Because of their relatively low DNA-binding affinity, SMAD complexes interact with a wide variety of DNA-binding proteins. Crosstalk with other signalling pathways and interaction with other DNA-binding cofactors define the specific binding patterns of SMADs; \
# in addition, interaction with coactivators/corepressors modulates their transcriptional activity.	Preferential formation of the regulatory R-Smad/SMAD4 heterotrimer over the R-Smad homotrimer is largely enthalpy driven, contributed by the unique presence of strong electrostatic interactions within the heterotrimeric interfaces. \
# These electrostatic interactions exist only in the heterotrimer due to specific charged residues in the SMAD4 subunit, Asp-493 and Arg-378, mediating complementary electrostatic interactions with the neighbouring R-Smad subunits.	\
# Heterotrimer	-	-	-	-	-	psi-mi:"MI:0469"(IntAct)	P84022(1)|Q13485(1)|Q15796(1)

SOURCE = 'EBI'
SOURCE_URL = 'https://www.ebi.ac.uk/complexportal/'


def process_file(input_file, output_file, source=SOURCE, source_url=SOURCE_URL):
        with gzip.open(input_file, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
            complex_tsv = csv.reader(infile, delimiter='\t')
            next(complex_tsv)
            for complex_row in complex_tsv:
                complex_ac = complex_row[0]
                molecules = complex_row[4].split('|')
                # skip lines containing chemicals from CHEBI or RNAs from RNACentral
                # i.e. only load complexes where all participants have uniprot protein ids
                if any([molecule.startswith('CHEBI:') or molecule.startswith('URS') for molecule in molecules]):
                    continue

                go_terms = complex_row[7].split('|')
                xrefs = complex_row[8].split('|')
                # each molecule/participant in column 5th is mostly a single protein_uniprot_id (stoichiometry), e.g. O60506(0)
                # but sometimes more complicated:
                # 1) a molecule set for paralogs grouped in an array, e.g. [O15342,Q8NHE4](1)
                # 2) have extra str after the uniprot ids, e.g. O60506-1(0) ('-1' for specifying the isoform);
                # O60895-PRO_0000030172 ('-PRO_0000030172' is the chain id)
                # 3) another complex accession id, e.g. CPX-973(1).

                # just load the original str for each molecule in complexes collection

                alias = [] if complex_row[2] == '-' else complex_row[2].split('|')
                experimental_evidence = None if complex_row[6] == '-' else complex_row[6]
                complex_assembly = [] if complex_row[11] == '-' else complex_row[11]
                reactome_xref = []
                for xref in xrefs:
                    if xref.startswith('reactome'):
                        reactome_xref.append(xref.replace('reactome:', ''))

                props = {
                    '_key': complex_ac,
                    'name': complex_row[1],
                    'alias': alias,
                    'molecules': molecules,
                    'evidence_code': complex_row[5],
                    'experimental_evidence': experimental_evidence,
                    'description': complex_row[9],
                    'complex_assembly': complex_assembly,
                    'complex_source': complex_row[17],
                    'reactome_xref': reactome_xref,
                    'source': source,
                    'source_url': source_url
                }
                outfile.write(json.dumps(props) + '\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file', type=str, required=True, help='In tsv.gz format')
    parser.add_argument('--output-file', type=str, required=True, help='Output file, in jsonl.gz')
    args = parser.parse_args()
    process_file(args.input_file, args.output_file)

if __name__ == '__main__':
    main()
