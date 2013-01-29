snpEff with HGVS
=================

Maven/Ant/Ivy have a hard time getting Picard and Samtools jars. Fetch these manually from the Picard distribution, alter the pom.xml version to match the version, then:

```
mvn install:install-file -DgroupId=net.sf.samtools -DartifactId=Sam -Dversion=1.79 -Dpackaging=jar -Dfile=/path_to_picard/1.79/sam-1.79.jar
mvn install:install-file -DgroupId=net.sf.picard -DartifactId=Picard -Dversion=1.79 -Dpackaging=jar -Dfile=/path_to_picard/1.79/picard-1.79.jar
```

The snpEff make script assumes a specific Eclipsey directory structure. Let's not muck with that.

To build snpEff-with-dependencies.jar in the ~/workspace/target directory:


```
mkdir -p ~/workspace/SnpEff  ~/workspace/SnpSift  ~/workspace/SnpSql
mkdir -p ~/snpEff

cd ~/workspace/SnpEff
git clone git@github.com:CBMi-BiC/snpEff.git .
./scripts/make.sh
```
