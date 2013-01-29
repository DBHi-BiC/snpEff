snpEff with HGVS
=================

The snpEff make script assumes a specific Eclipsey directory structure. Let's not muck with that.

```
mkdir ~/workspace ~/snpEff
cd ~/workspace
git clone git@github.com:CBMi-BiG/snpEff.git SnpEff
cd SnpEff
```

Maven/Ant/Ivy have a hard time getting Picard and Samtools jars, but they are in the Picard distribution. Alter the pom.xml versions as necessary.

```
wget http://downloads.sourceforge.net/project/picard/picard-tools/1.84/picard-tools-1.84.zip
unzip picard-tools-1.84.zip
mvn install:install-file -DgroupId=net.sf.samtools -DartifactId=Sam -Dversion=1.84 \
-Dpackaging=jar -Dfile=$PWD/picard-tools-1.84/sam-1.84.jar
mvn install:install-file -DgroupId=net.sf.picard -DartifactId=Picard -Dversion=1.84 \
-Dpackaging=jar -Dfile=$PWD/picard-tools-1.84/picard-1.84.jar
```



To build snpEff-with-dependencies.jar in the ~/workspace/SnpEff/target directory:


```
./scripts/make.sh
```
