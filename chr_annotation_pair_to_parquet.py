import sys
from hashlib import sha256

from awsglue.transforms import *
from awsglue.utils import getResolvedOptions
from pyspark.context import SparkContext
from awsglue.context import GlueContext
from awsglue.job import Job
from pyspark.sql.functions import col, lit, udf
from pyspark.sql.types import StringType

def parse_base_pairs_from_uniq_id(uniq_id):
    return ':'.join(uniq_id.split(':')[1:3])

def generate_unique_key(*args):
    concatenated_value = ''.join(args) + 'GRCh38'
    return sha256(concatenated_value.encode()).hexdigest()

generate_unique_key_udf = udf(generate_unique_key, StringType())

# Get command-line arguments
args = getResolvedOptions(sys.argv,
            [
                'JOB_NAME',
                'annotations_path',
                'lds_path',
                'label',
                'source',
                'source_url',
                'chr',
                'ancestry',
                'output_path',
                'write_mode',
            ]
       )

# Initialize GlueContext and SparkSession
sc = SparkContext()
glueContext = GlueContext(sc)
job = Job(glueContext)
job.init(args['JOB_NAME'], args)

spark = glueContext.spark_session

process_uniq_id_udf = udf(parse_base_pairs_from_uniq_id, StringType())

# Read input files from S3
annotations = spark.read.option('header', 'true').csv(args['annotations_path'])
lds = spark.read.option('header', 'true').csv(args['lds_path'])

joined_lds_1 = lds.join(annotations, col('SNP1') == col('Position'), 'left').select(col('rsID').alias('variant_1_rsid'), col('SNP2'), col('Uniq_ID_1'), col('Uniq_ID_2'), col('R2'), col('Dprime'), col('+/-corr'))
joined_lds_2 = joined_lds_1.join(annotations, col('SNP2') == col('Position'), 'left').select(col('variant_1_rsid'), col('rsID').alias('variant_2_rsid'), col('Uniq_ID_1'), col('Uniq_ID_2'), col('R2'), col('Dprime'), col('+/-corr'))

parsed_lds = joined_lds_2.withColumn('inverted', col('+/-corr') == '+').withColumn('variant_1_base_pair', process_uniq_id_udf(col('Uniq_ID_1'))).withColumn('variant_2_base_pair', process_uniq_id_udf(col('Uniq_ID_2')))

parsed_lds = parsed_lds.drop('Uniq_ID_1', 'Uniq_ID_2', '+/-corr').withColumn('R2', col('R2').cast('float')).withColumn('Dprime', col('Dprime').cast('float'))

parsed_lds = parsed_lds.withColumn('label', lit(args['label'])).withColumn('source', lit(args['source'])).withColumn('source_url', lit(args['source_url'])).withColumn('chr', lit(args['chr'])).withColumn('ancestry', lit(args['ancestry']))

parsed_lds = parsed_lds.withColumn('key', generate_unique_id_udf(col('ancestry'), col('chr'), col('variant_1_base_pair'), col('variant_2_base_pair'))) 

# Write output as parquet to S3
parsed_lds.write.mode(args['write_mode']).parquet(args['output_path'])

job.commit()
