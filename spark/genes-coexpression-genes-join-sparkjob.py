import argparse
import pyspark.sql.functions as F

from pyspark.sql import SparkSession
from pyspark.sql.function import col


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genes', type=str, required=True, help='Genes jsonl.gz')
    parser.add_argument('--coexpressions', type=str, required=True, help='Coexpressions jsonl.gz')
    parser.add_argument('--output-parquet', type=str, required=True, help='Location to write the output parquet')
    parser.add_argument('--output-partitions', type=int, default=1, help='Number of output partitions. Default 32')
    parser.add_argument('--logit-score-cutoff', type=float, default=3.0, help='Include coexpressions with abs >=  cutoff')
    args = parser.parse_args()
    spark = SparkSession.builder.appName('Join genes and coexpressions').getOrCreate()
    genes = spark.read.json(args.genes)
    coexpressions = spark.read.json(args.coexpressions)
    filtered_coexpressions = coexpressions.filter((col('logit_score') <= -args.logit_score_cutoff) | (col('logit_score_cutoff') >= args.logit_score_cutoff))
    jointable = (
        coexpressions.alias('coex')
        .join(genes.alias('g1'), col('coex._source') == col('g1._key'), how='inner')
        .join(genes.alias('g2'), col('coex._target') == col('g2._key'), how='inner')
        .select(
            col('coex.label').alias('coexpression_label'),
            col('coex.logit_score'),
            col('coex.source').alias('coexpression_source'),
            col('coex.source_url').alias('coexpression_source_url'),
            col('g1._key').alias('source_key'),
            col('g1.alias').alias('source_alias'),
            col('g1.chr').alias('source_chr'),
            col('g1.end:long').alias('source_end'),
            col('g1.entrez').alias('source_entrez'),
            col('g1.gene_id').alias('source_gene_id'),
            col('g1.gene_type').alias('source_gene_type'),
            col('g1.hgnc').alias('source_hgnc'),
            col('g1.name').alias('source_name'),
            col('g1.source').alias('source_source'),
            col('g1.source_url').alias('source_source_url'),
            col('g1.start:long').alias('source_start'),
            col('g1.version').alias('source_version'),
            col('g2._key').alias('target_key'),
            col('g2.alias').alias('target_alias'),
            col('g2.chr').alias('target_chr'),
            col('g2.end:long').alias('target_end'),
            col('g2.entrez').alias('target_entrez'),
            col('g2.gene_id').alias('target_gene_id'),
            col('g2.gene_type').alias('target_gene_type'),
            col('g2.hgnc').alias('target_hgnc'),
            col('g2.name').alias('target_name'),
            col('g2.source').alias('target_source'),
            col('g2.source_url').alias('target_source_url'),
            col('g2.start:long').alias('target_start'),
            col('g2.version').alias('target_version')
        )
    )
    jointable.repartition(args.output_partitions).write.mode('overwrite').parquet(args.output_parquet)
