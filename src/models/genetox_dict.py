
import numpy as np
import pandas as pd


genetox = {'Ames': 'bacterial reverse mutation test', 
           'Ames study' :'bacterial reverse mutation test',
           'Ames II' : 'bacterial reverse mutation test',
           'bacterial reverse mutation assay (e.g. Ames test)' : 'bacterial reverse mutation test',
           'Bacterial Mutagenesis' : 'bacterial reverse mutation test (single strain)',
           'Histidine reverse gene mutation, Ames assay' : 'bacterial reverse mutation test',
           'according to B Ames et al. (1975) Mutat Res 31:347-364': 'bacterial reverse mutation test',
           'Cell Transformation' : 'cell transformation',
           'Cell transformation' : 'in vitro cell transformation assay',
           'Chromosome aberrations, in vivo' : 'in vivo chromosome aberrations',
           'Chromosome aberrations in vivo' : 'in vivo chromosome aberrations',
           'Chromosome aberrations in vitro' : 'in vitro mammalian chromosome aberration test',
           'chromosomal aberration test in Chinese hamster lung cells (CHL/IU)' : 'in vitro mammalian chromosome aberration test',
           'In Vitro Chromosome Aberration' : 'in vitro mammalian chromosome aberration test',
           'Chromosomal aberation' : 'in vitro mammalian chromosome aberration test',
           'Chromosomal aberration assay' : 'in vitro mammalian chromosome aberration test',
           'Micronucleus test in vitro, chromosome aberrations' : 'in vitro mammalian chromosome aberration test',
           'Chromosome aberrations' : 'chromosome aberrations (plant)',
           'in vitro mammalian cytogenicity (B10)' :'in vitro mammalian chromosome aberration test',
           'in vitro mammalian cell transformation assay' : 'in vitro cell transformation assay',
           'Chinese hamster ovary cell/hypoxanthine-guanine-phosphoribosyl transferase (CHO/HGPRT) forward gene mutation assay': 'in vitro mammalian cell gene mutation test using the Hprt and xprt genes',
           'in vitro mammalian cell gene mutation': 'in vitro mammalian cell gene mutation test using the Hprt and xprt genes',
           'mammalian cell gene mutation assay' : 'in vitro mammalian cell gene mutation test',
          'in vitro mammalian cell gene mutation (B.17)' :'in vitro mammalian cell gene mutation test using the Hprt and xprt genes', 'in vitro mammalian cell gene mutation (B17)' : 'in vitro mammalian cell gene mutation test using the Hprt and xprt genes',
           'in vitro mammalian cytogenecity (B10)' : 'in vitro mammalian chromosome aberration test', 
           'in vitro mammalian cytogenicity' : 'in vitro mammalian chromosome aberration test', 'in vitro mammalian cytogenicity (B10)' : 'in vitro mammalian chromosome aberration test',
           'in vitro mammallian cytogenicity (B10)' : 'in vitro mammalian chromosome aberration test',
           'in vitro mammalian chromosome aberration test' : 'in vitro mammalian chromosome aberration test',
           'in vitro mammalian chromosome aberration test, human lymphocytes from healthy, non smoking donors' : 'in vitro mammalian chromosome aberration test',
           'Mouse lymphoma assay' : 'in vitro mammalian cell gene mutation test using the thymidine kinase gene',
           'In vitro L5178Y TK+/- Mouse Lymphoma Cell Assay': 'in vitro mammalian cell gene mutation test using the thymidine kinase gene',
           'DNA Damage/Repair' : 'DNA damage/repair',
           'DNA repair' : 'DNA damage/repair',
           'DNA damage' : 'DNA damage/repair',
           'unscheduled DNA synthesis' : 'unscheduled DNA synthesis',
           'DNA damage/gene conversion' : 'DNA damage/repair',
           'Transgenic' : 'transgenic rodent somatic and germ cell gene mutation assays (TGR)',
           'GENE MUTATION ASSAY IN CHINESE HAMSTER V79 CELLS IN VITRO' :'in vitro mammalian cell gene mutation test',
           'bacterial gene mutation assay' : 'bacterial reverse mutation test',
           'Forward gene mutation at the HPRT or ouabain locus' : 'in vitro mammalian cell gene mutation test using the Hprt and xprt genes',
           'Forward gene mutation at the HPRT locus': 'in vitro mammalian cell gene mutation test using the Hprt and xprt genes',
           'Forward gene mutation at the thymidine kinase (TK) locus; chromosome aberrations': 'in vitro mammalian cell gene mutation test using the thymidine kinase gene',
           'Gene mutation' : 'gene mutation (plant)',
           'In Vivo Micronucleus' : 'in vivo micronucleus test',
           'In Vitro Micronucleus' : 'in vitro mammalian cell micronucleus test',
           'in vitro micronucleus assay' : 'in vitro mammalian chromosome aberration test',
           'micronucleus assay' : 'in vitro mammalian chromosome aberration test',
           'in vitro micronucleus assay' : 'in vitro mammalian cell micronucleus test',
           'In vitro micronucleus test in mouse lymphoma L5178Y cells' : 'in vitro mammalian cell gene mutation test using the thymidine kinase gene',
           'In-vitro micronucleus test in cultured human lymphocytes' : 'in vitro mammalian chromosome aberration test',
           'In Vitro Mammalian Cell Micronucleus Test (MNvit)' : 'in vitro mammalian chromosome aberration test',
           'Unscheduled DNA synthesis (UDS) in vitro, DNA effects' : 'unscheduled DNA synthesis (UDS) in vitro',
           'unscheduled DNA synthesis in mammalian cells in vitro' : 'unscheduled DNA synthesis (UDS) in vitro',
           'Unscheduled DNA synthesis (UDS) in vivo' : 'unscheduled DNA synthesis (UDS) in vivo',
           'Unscheduled DNA synthesis' : 'unscheduled DNA synthesis (UDS) in vivo',
           'DNA damage and repair assay, unscheduled DNA synthesis in mammalian cells in vitro' :'unscheduled DNA synthesis (UDS) in vitro',
           'In Vivo Chromosome Aberration' : 'in vivo chromosome aberrations',
'Micronucleus' : 'in vivo micronucleus test',
'Micronucleus test, chromosome aberrations': 'mammalian erythrocyte micronucleus test',
           'Unscheduled DNA synthesis (UDS) in vivo, DNA effects' :'unscheduled DNA synthesis (UDS) in vivo',
           'Unscheduled DNA synthesis (UDS) in vivo; DNA effects' :'unscheduled DNA synthesis (UDS) in vivo',
           'in vitro mammalian cell gene mutation assay' : 'in vitro mammalian cell gene mutation test',
           'In vitro mammalian cell gene mutation test' : 'in vitro mammalian cell gene mutation test',
           }


assay_result_std = {'equivocal' : 'inconclusive', 'inconsistent (cancer in vivo)' : 'inconclusive', 'ambiguous' :'inconclusive', 'negative (cancer in vivo)' : 'negative', 'positive (cancer in vivo)' : 'positive', 'positive (weak)' : 'positive', 'no data' : 'not determined', 'technically compromised' : 'inadequate', 'uninterpretable' : 'inadequate'}

outcome = {'inadequate' : 2, 'inconclusive' : 2, 'not determined' : 2, 'positive' : 1, 'negative' : 0}

def clean_up(df):
    df['standard_assay_type'] = df['assay_type']
    df['standard_assay_type'].replace(genetox,  inplace = True)
    df['assay_result_std'] = df['assay_result']
    df['assay_result_std'].replace(assay_result_std, inplace = True)
    df['assay_outcome'] = df['assay_result_std']
    df['assay_outcome'].replace(outcome, inplace = True)
    return df
