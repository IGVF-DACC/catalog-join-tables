import sys
from awsglue.transforms import *
from awsglue.utils import getResolvedOptions
from pyspark.context import SparkContext
from awsglue.context import GlueContext
from awsglue.job import Job
  
sc = SparkContext.getOrCreate()
glueContext = GlueContext(sc)
spark = glueContext.spark_session
job = Job(glueContext)
genes_genes_df = spark.read.json('s3://catalog-columnar/jsonl/genes_genes/*.jsonl.gz')
genes_genes_df.repartition(8).write.mode('overwrite').parquet('s3://catalog-columnar/parquet/genes_genes/')
job.commit()
