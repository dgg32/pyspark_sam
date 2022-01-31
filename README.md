# pyspark_sam

This repo hosts my code for the article "[Analyze Big Sequence Alignments with PySpark in AWS EMR](https://medium.com/p/a044acaa60af)".


# Prerequisite

1. Spark

2. AWS CLI

3. AWS Account

# Run

Follow the instruction in the article. Once you have uploaded the files into your S3 bucket, run

```console
aws emr create-cluster --name "Spark_step_pip" \
    --release-label emr-6.5.0 \
    --applications Name=Spark \
    --log-uri s3://[your_S3_bucket]/logs/ \
    --instance-type m5.xlarge \
    --instance-count 3 \
    --bootstrap-actions Path=s3://[your_S3_bucket]/emr_bootstrap.sh \
    --use-default-roles --auto-terminate \
    --steps "Type=Spark,Name=SparkProgram,ActionOnFailure=CONTINUE,Args=[--deploy-mode,cluster,--master,yarn,--py-files,s3://[your_S3_bucket]/helper_function.py,s3://[your_S3_bucket]/spark_3mer.py,s3://[your_S3_bucket]/test.sam,[your_S3_bucket],sankey.json]" 
```

When the job finishes, download the sankey.json. And run this command to visualize:

```console
python sankey.py sankey.json
```

## Authors

  

*  **Sixing Huang** - *Concept and Coding*

  

## License

  

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
