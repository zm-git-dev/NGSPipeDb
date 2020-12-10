# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
rule sampling_data_by_links:
    input:
        aread1 = join(config["samplesDir"], "{sample}"+config["read1Suffix"]),
        aread2 = join(config["samplesDir"], "{sample}"+config["read2Suffix"])
    output:
        bread1 = join(sampling_data_outdir, "{sample}"+config["read1Suffix"]),
        bread2 = join(sampling_data_outdir, "{sample}"+config["read2Suffix"]),
    threads:
        5
    log:
        join(sampling_data_outdir, "logfile", "{sample}.log")
    run:
        import os
        read1_src = os.path.abspath(input.aread1)
        read2_src = os.path.abspath(input.aread2)
        os.symlink(read1_src, output.bread1)
        os.symlink(read2_src, output.bread2)