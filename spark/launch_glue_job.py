import argparse
import boto3
import json

def get_glue_client(profile_name, region_name):
    session = boto3.session.Session(profile_name=profile_name, region_name=region_name)
    client = session.client('glue')
    return client

def get_glue_args_from_json(json_filename):
    with open(json_filename) as fp:
        glue_job_args = json.load(fp)
    return glue_job_args

def main(args):
    client = get_glue_client(args.profile_name, args.region_name)
    glue_args = get_glue_args_from_json(args.input_json_path)
    for argument_item in glue_args:
        client.start_job_run(JobName=args.glue_job_name, Arguments=argument_item)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile-name', type=str, required=True)
    parser.add_argument('--region-name', type=str, required=True)
    parser.add_argument('--input-json-path', type=str, help='json file with list of json objects with args', required=True)
    parser.add_argument('--glue-job-name', type=str, required=True)
    args = parser.parse_args()
    main(args)
