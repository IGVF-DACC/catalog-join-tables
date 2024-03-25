import sys
from awsglue.transforms import *
from awsglue.utils import getResolvedOptions
from pyspark.context import SparkContext
from awsglue.context import GlueContext
from awsglue.job import Job

from pyspark.sql.functions import col

#get command line args
args = getResolvedOptions(sys.argv,
            [
                'JOB_NAME',
                'output_prefix',
                'write_mode',
                'ancestry',
                'variants_input_path',
                'lds_input_prefix',
            ]
        )

sc = SparkContext()
glueContext = GlueContext(sc)
job = Job(glueContext)
job.init(args['JOB_NAME'], args)
spark = glueContext.spark_session

variants_variants = spark.read.parquet(args['lds_input_prefix']+ '/' + args['ancestry'])
variants = spark.read.parquet(args['variants_input_path'])
variants = variants.withColumnRenamed('source', 'variants_source').withColumnRenamed('source_url', 'variants_source_url').withColumnRenamed('chr', 'variants_chr')

joined_df1 = variants_variants.join(
    variants.select(col('rsid').alias('variant_1_rsid'), col('*')),
    on='variant_1_rsid',
    how='inner'
).select(
    col('variant_1_rsid'),
    col('variant_2_rsid'),
    col('R2'),
    col('Dprime'),
    col('inverted'),
    col('variant_1_base_pair'),
    col('variant_2_base_pair'),
    col('label'),
    col('source'),
    col('source_url'),
    col('chr'),
    col('ancestry'),
    *[col(c).alias(f'variant_1_{c}') for c in variants.columns if c != 'rsid']
)

joined_df2 = joined_df1.join(
    variants.select(col('rsid').alias('variant_2_rsid'), col('*')),
    on='variant_2_rsid',
    how='inner'
).select(
    col('variant_1_rsid'),
    col('variant_2_rsid'),
    col('R2'),
    col('Dprime'),
    col('inverted'),
    col('variant_1_base_pair'),
    col('variant_2_base_pair'),
    col('label'),
    col('source'),
    col('source_url'),
    col('chr'),
    col('ancestry'),
    *[col(c) for c in joined_df1.columns if c.startswith('variant_1_') and c not in ['variant_1_base_pair', 'variant_1_rsid']],
    *[col(c).alias(f'variant_2_{c}') for c in variants.columns if c != 'rsid']
)

joined_df2.write.mode('overwrite').parquet(args['output_prefix'] + '/' + 'variants-variants_variants-' + args['ancestry'] + '.parquet')
job.commit()
