[](http://www.network-science.de/ascii/)
<pre>
 **    **     **        *******
//**  **     ****      /**////**
 //****     **//**     /**   /**
  //**     **  //**    /*******
   /**    **********   /**////
   /**   /**//////**   /**
   /**   /**     /**   /**
   //    //      //    //
</pre>

# Install
To install this latest version:
```shell
pip install git+https://github.com/DingWB/cemba_data
```

# Documentation
## 1. Make sure create the right environment
git clone https://github.com/DingWB/cemba_data.git
mamba env create -f cemba_data/env.yaml
## 2. Generate config.ini
```shell
yap default-mapping-config --mode m3c --barcode_version V2 --bismark_ref ${HOME}/Ref/hg38/hg38_ucsc_with_chrL.bismark1 --genome ${HOME}/Ref/hg38/hg38_ucsc_with_chrL.fa --chrom_size_path ${HOME}/Ref/hg38/hg38_ucsc.main.chrom.sizes > config.ini
```
## 3. Demultiplex
```shell
yap demultiplex --fastq_pattern "test_fastq/*.gz" -o mapping -j 4 --aligner bismark --config_path config.ini

```
## 4. Run mapping
```shell
sh mapping/snakemake/qsub/snakemake_cmd.txt
```

# YAP (Yet Another Pipeline)
Pipeline(s) for mapping and cluster-level aggregation of single nucleus methylome and multi-omic datasets.
Technologies supported:
- snmC-seq(1/2/3)
- snmCT-seq (mC + RNA)
- snmC2T-seq (mC + RNA + Chromatin Accessibility)
- snm3C-seq (mC + Chromatin Conformation)
- any NOMe treated version of the above

[See Documentation](https://hq-1.gitbook.io/mc/)
