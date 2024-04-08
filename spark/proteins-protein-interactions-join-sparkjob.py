import argparse
import pyspark.sql.functions as F
from pyspark.sql import SparkSession
from pyspark.sql.functions import col


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--proteins', type=str, required=True, help='Proteins as jsonl.gz')
    parser.add_argument('--protein-interactions', type=str, required=True, help='Protein interactions jsonl.gz')
    parser.add_argument('--output-parquet', type=str, required=True, help='Location to write the output parquet')
    parser.add_argument('--output-partitions', type=int, default=1, help='Number of output partitions. Default: 1')
    args = parser.parse_args()
    spark = SparkSession.builder.appName('Join proteins and protein interactions').getOrCreate()
    proteins = spark.read.json(args.proteins)
    protein_interactions = spark.read.json(args.protein_interactions)
    # make dbxrefs an array here, so that it can be ingested as multi value dimension in druid
    reformatted_proteins = proteins.withColumn('dbxrefs', F.transform(col('dbxrefs'), lambda x: x['id']))
    jointable = (
        protein_interactions.alias('inter')
        .join(reformatted_proteins.alias('p1'), col('inter._from') == col('p1._key'), how='inner')
        .join(reformatted_proteins.alias('p2'), col('inter._to') == col('p2._key'), how='inner')
        .select(
            col('inter._from').alias('source_protein'),
            col('inter._to').alias('target_protein'),
            col('inter.confidence_value_biogrid:long').alias('confidence_value_biogrid'),
            col('inter.confidence_value_intact:long').alias('confidence_value_intact'),
            col('inter.detection_method'),
            col('inter.detection_method_code'),
            col('inter.interaction_type'),
            col('inter.interaction_type_code'),
            col('inter.organism'),
            col('inter.pmids'),
            col('inter.source').alias('interaction_source'),
            col('p1.dbxrefs').alias('source_dbxrefs'),
            col('p1.full_name').alias('source_full_name'),
            col('p1.name').alias('source_name'),
            col('p1.source').alias('source_source'),
            col('p1.source_url').alias('source_source_url'),
            col('p2.dbxrefs').alias('target_dbxrefs'),
            col('p2.full_name').alias('target_full_name'),
            col('p2.name').alias('target_name'),
            col('p2.source').alias('target_source'),
            col('p2.source_url').alias('target_source_url')
        )
    )
    jointable.repartition(args.output_partitions).write.parquet(args.output_parquet)

if __name__ == '__main__':
    main()
