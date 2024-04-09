import argparse
import pyspark.sql.functions as F

from pyspark.sql import SparkSession
from pyspark.sql.function import col


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genes', type=str, required=True, help='Genes jsonl.gz')
    parser.add_argument('--coexpressions', type=str, required=True, help='Coexpressions jsonl.gz')
    parser.add_argument('--output-parquet', type=str, required=True, help='Location to write the output parquet')
    parser.add_argument('--output-partitions', type=int, default=32, help='Number of output partitions. Default 32')
    args = parser.parse_args()

