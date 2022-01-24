import sys
from pyspark import SparkContext
from operator import add
import json
import helper_function
import boto3

if __name__=="__main__":

    input_file = sys.argv[1]
    output_s3 = sys.argv[2]
    output_filename = sys.argv[3]

    sc = SparkContext().getOrCreate()

    lines=sc.textFile(input_file).filter(lambda x: x.startswith("@") == False)
    objects = lines.map(helper_function.get_sam_object)
    k_v = objects.flatMap(helper_function.mapper)
    result = k_v.reduceByKey(add).collect()

    data = helper_function.format_result(result)

    ###output option 1: rdd saveAsTextFile
    #rdd = sc.parallelize([data])

    #rdd.coalesce(1, shuffle = True).saveAsTextFile(output_s3)

    ###output option 2: boto3
    s3 = boto3.resource(
        's3',
        region_name= 'us-east-1'
    )

    s3.Object(output_s3, output_filename).put(Body=json.dumps(data, indent=2))
